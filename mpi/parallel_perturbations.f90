! One collective event; retain the historical intrinsic RNG and k/j/i draw order.
module parallel_perturbations
  use variables
  use parallel_runtime
  use parallel_timing, only: time_decisions,coordinate_time_actions
  use parallel_diagnostics, only: state_diagnostics,collect_global_state, &
    validate_observable_sample,accumulate_perturbation_budget,reduce_observable_sum
  use domain_decomposition, only: domain_layout,build_domain,validate_domain
  use conditions, only: inject_pressure_noise,inject_density_noise,inject_density_mode,set_reproducible_random_seed
  use equations, only: prim_to_cons
  use halo_exchange, only: fill_primitive_halos
  use output, only: write_parallel_perturbation_event
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: apply_parallel_perturbation
contains

  subroutine agree_controls(layout)
    type(domain_layout), intent(in) :: layout
    character(len=160) :: message
    real*8 :: values(22),minimum(22),maximum(22)
    logical :: valid
#ifdef USE_MPI
    integer :: status
#endif
    call parallel_require_active()
    call validate_domain(layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    valid=layout%rank==parallel_rank.and.layout%ranks==parallel_size.and.layout%halo==nghost
    valid=valid.and.all(layout%local_shape==[nx,ny,nz]).and. &
      all(layout%global_shape==[nx_global,ny_global,nz_global]).and. &
      all(layout%offset(1:2)==[x_global_offset,y_global_offset])
    call parallel_require_collective(valid,'perturbation layout does not match the allocated grid')
    values=[perturbation_amplitude,perturbation_mode,perturbation_time,integration_time,final_time,adb_idx,g1, &
      dble(perturbation_seed),dble(perturbation_type),merge(1.0d0,0.0d0,apply_perturbation), &
      merge(1.0d0,0.0d0,perturbation_applied),dble(n_steps), &
      dble(layout%global_shape),dble(layout%process_shape),rho_floor,p_floor,bh_mass,a_spin, &
      merge(1.0d0,0.0d0,accelerator_enabled)]
    valid=all(ieee_is_finite(values))
    if (valid)valid=perturbation_amplitude>=0.0d0.and.perturbation_amplitude<1.0d0.and. &
      perturbation_time>=0.0d0.and.adb_idx>1.0d0.and.g1>1.0d0
    valid=valid.and.perturbation_type>=PERT_NONE.and.perturbation_type<=PERT_DENSITY_MODE
    if (perturbation_type==PERT_DENSITY_MODE)valid=valid.and.perturbation_mode>=1.0d0
    valid=valid.and.(.not.apply_perturbation.or.perturbation_type/=PERT_NONE)
    call parallel_require_collective(valid,'invalid perturbation controls')
    minimum=values;maximum=values
#ifdef USE_MPI
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(values,minimum,size(values),MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_perturbations:agree_controls:1',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'perturbation metadata minimum')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(values,maximum,size(values),MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_perturbations:agree_controls:2',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'perturbation metadata maximum')
#endif
    call parallel_require_collective(all(minimum==maximum),'perturbation controls differ between ranks')
  end subroutine agree_controls

  subroutine global_noise(layout,rho_peak,noise)
    type(domain_layout), intent(in) :: layout
    real*8, intent(in) :: rho_peak
    real*8, intent(out) :: noise(nx,ny,nz)
    logical, allocatable :: local_mask(:),packed_mask(:),mask(:,:,:)
    real*8, allocatable :: global_values(:,:,:),packed_values(:)
    integer, allocatable :: counts(:),displacements(:)
    type(domain_layout) :: block
    integer :: local_cells,total,rank_id,offset,i,j,k,d,ios
    integer(int64) :: total64
    logical :: valid
    character(len=160) :: message
#ifdef USE_MPI
    integer :: status
#endif
    total64=1_int64;valid=.true.
    do d=1,3
      if (total64>int(huge(1),int64)/layout%global_shape(d)) then
        valid=.false.;exit
      end if
      total64=total64*layout%global_shape(d)
    end do
    call parallel_require_collective(valid,'perturbation mask exceeds classic MPI count range')
    total=int(total64);local_cells=nx*ny*nz
    noise=0.0d0
    if (parallel_size==1) then
      call set_reproducible_random_seed(perturbation_seed)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            if (p(i,j,k,eq_de)>0.05d0*rho_peak)call random_number(noise(i,j,k))
          end do
        end do
      end do
      return
    end if
    allocate(local_mask(local_cells),counts(parallel_size),displacements(parallel_size),stat=ios)
    call parallel_require_collective(ios==0,'cannot allocate local perturbation mask')
    local_mask=reshape(p(1:nx,1:ny,1:nz,eq_de)>0.05d0*rho_peak,[local_cells])
    counts=0;displacements=0
    if (parallel_rank==0) then
      allocate(packed_mask(total),stat=ios)
      if (ios==0) then
        offset=0
        do rank_id=0,parallel_size-1
          call build_domain(layout%global_shape,layout%halo,rank_id,parallel_size,layout%periodic_x, &
            block,message,layout%process_shape,layout%periodic_y)
          counts(rank_id+1)=product(block%local_shape)
          displacements(rank_id+1)=offset;offset=offset+counts(rank_id+1)
        end do
      end if
    else
      allocate(packed_mask(1),stat=ios)
    end if
    call parallel_require_collective(ios==0,'cannot allocate gathered perturbation mask')
#ifdef USE_MPI
    call MPI_Gatherv(local_mask,local_cells,MPI_LOGICAL,packed_mask,counts,displacements, &
      MPI_LOGICAL,0,MPI_COMM_WORLD,status)
    call parallel_check_mpi(status,'gather perturbation body mask')
#endif
    deallocate(local_mask)
    if (parallel_rank==0) then
      allocate(mask(nx_global,ny_global,nz_global),stat=ios)
      if (ios==0) then
        do rank_id=0,parallel_size-1
          call build_domain(layout%global_shape,layout%halo,rank_id,parallel_size,layout%periodic_x, &
            block,message,layout%process_shape,layout%periodic_y)
          mask(block%offset(1)+1:block%offset(1)+block%local_shape(1), &
            block%offset(2)+1:block%offset(2)+block%local_shape(2),:)= &
            reshape(packed_mask(displacements(rank_id+1)+1:displacements(rank_id+1)+counts(rank_id+1)), &
              block%local_shape)
        end do
      end if
    else
      ios=0
    end if
    deallocate(packed_mask)
    call parallel_require_collective(ios==0,'cannot allocate global perturbation mask')
    if (parallel_rank==0) then
      allocate(global_values(nx_global,ny_global,nz_global),stat=ios)
      if (ios==0) then
        global_values=0.0d0
        call set_reproducible_random_seed(perturbation_seed)
        do k=1,nz_global
          do j=1,ny_global
            do i=1,nx_global
              if (mask(i,j,k))call random_number(global_values(i,j,k))
            end do
          end do
        end do
      end if
      deallocate(mask)
    end if
    call parallel_require_collective(ios==0,'cannot allocate global perturbation values')
    allocate(packed_values(merge(total,1,parallel_rank==0)),stat=ios)
    call parallel_require_collective(ios==0,'cannot allocate packed perturbation values')
    if (parallel_rank==0) then
      do rank_id=0,parallel_size-1
        call build_domain(layout%global_shape,layout%halo,rank_id,parallel_size,layout%periodic_x, &
          block,message,layout%process_shape,layout%periodic_y)
        packed_values(displacements(rank_id+1)+1:displacements(rank_id+1)+counts(rank_id+1))= &
          reshape(global_values(block%offset(1)+1:block%offset(1)+block%local_shape(1), &
            block%offset(2)+1:block%offset(2)+block%local_shape(2),:),[counts(rank_id+1)])
      end do
      deallocate(global_values)
    end if
#ifdef USE_MPI
    call MPI_Scatterv(packed_values,counts,displacements,MPI_DOUBLE_PRECISION,noise,local_cells, &
      MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,status)
    call parallel_check_mpi(status,'scatter canonical perturbation values')
#endif
  end subroutine global_noise

  subroutine apply_parallel_perturbation(layout,injected)
    type(domain_layout), intent(in) :: layout
    logical, intent(out), optional :: injected
    type(time_decisions) :: decisions
    type(state_diagnostics) :: state
    real*8, allocatable :: noise(:,:,:)
    real*8 :: change(neq),previous(neq),packet(neq+1)
    integer :: i,j,k,ios
    integer(int64) :: selected
    call agree_controls(layout)
    if (present(injected))injected=.false.
    call coordinate_time_actions(integration_time,n_steps,final_time,[perturbation_time,final_time,final_time], &
      [apply_perturbation.and..not.perturbation_applied,.false.,.false.],[.false.,.false.],decisions)
    if (.not.decisions%perturb.or.decisions%complete) return
    call parallel_require_collective(.not.accelerator_enabled,'parallel perturbation supports CPU only')
    call validate_observable_sample(integration_time)
    call fill_primitive_halos(p)
    call collect_global_state(state)
    call parallel_require_collective(state%invalid_cells==0,'unphysical state before perturbation')
    selected=count(p(1:nx,1:ny,1:nz,eq_de)>0.05d0*state%rho_max,kind=int64)
    if (perturbation_type==PERT_DENSITY_MODE) then
      call inject_density_mode(p,perturbation_mode,perturbation_amplitude,state%rho_max)
    else
      allocate(noise(nx,ny,nz),stat=ios)
      call parallel_require_collective(ios==0,'cannot allocate perturbation noise')
      call global_noise(layout,state%rho_max,noise)
      if (perturbation_type==PERT_PRESSURE_NOISE) then
        call inject_pressure_noise(p,perturbation_amplitude,perturbation_seed,state%rho_max,noise)
      else
        call inject_density_noise(p,perturbation_amplitude,perturbation_seed,state%rho_max,noise)
      end if
    end if
    change=0.0d0
    !$OMP PARALLEL DO PRIVATE(i,j,k,previous) REDUCTION(+:change)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          previous=up(i,j,k,:)
          call prim_to_cons(p(i,j,k,:),u(i,j,k,:),i,j,k)
          up(i,j,k,:)=u(i,j,k,:)
          change=change+(up(i,j,k,:)-previous)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call accumulate_perturbation_budget(change)
    call fill_primitive_halos(p)
    perturbation_applied=.true.
    packet=[change*dx*dy*dz,dble(selected)]
    call reduce_observable_sum(packet)
    call write_parallel_perturbation_event(integration_time,state%rho_max,int(packet(neq+1),int64),packet(1:neq))
    if (present(injected))injected=.true.
  end subroutine apply_parallel_perturbation
end module parallel_perturbations
