program test_domain_decomposition
  use parallel_runtime
  use domain_decomposition
  use variables
  use initialization, only: allocate_and_grid, precalculate_metric_cache
  use metrics
  use conditions, only: set_initial_conditions
  use equations, only: prim_to_cons
  implicit none

  integer :: fixture,depth,report_unit,checks=0
  real*8 :: primitive_error=0.0d0,conserved_error=0.0d0,primitive_ulps=0.0d0,conserved_ulps=0.0d0
  character(len=40) :: mode
  character(len=80) :: report_name
  character(len=160) :: message
  type(domain_layout) :: block

  call parallel_initialize()
  call get_command_argument(1,mode)
  if (trim(mode)=='bad-small') then
    call build_domain([11,1,1],3,parallel_rank,parallel_size,.false.,block,message)
    call parallel_require_collective(len_trim(message)==0,message)
    call parallel_abort('bad-small expected a rejected grid')
  else if (trim(mode)=='bad-one-rank') then
    call build_domain([37,1,1],3,parallel_rank,parallel_size,.false.,block,message)
    if (parallel_rank==parallel_size-1) block%offset(1)=block%offset(1)+1
    call validate_domain(block,message)
    call parallel_require_collective(len_trim(message)==0,message)
    call parallel_abort('bad-one-rank expected a rejected layout')
  end if
  call check_partition()
  do depth=1,3
    do fixture=1,9
      call check_grid(fixture,depth)
    end do
  end do
  if (trim(mode)=='report') then
    write(report_name,'(A,I0,A)') 'domain_rank_',parallel_rank,'.txt'
    open(newunit=report_unit,file=trim(report_name),status='new',action='write')
    write(report_unit,'(3(I0,1X),4(ES24.16E3,1X))') parallel_rank,parallel_size,checks, &
      primitive_error,conserved_error,primitive_ulps,conserved_ulps
    close(report_unit)
  end if
  write(*,'(A,I0,A,I0,A,I0)') 'DOMAIN_TEST rank=',parallel_rank,' size=',parallel_size,' checks=',checks
  write(*,'(A,I0,4(A,ES22.14E3))') 'DOMAIN_ROUNDOFF rank=',parallel_rank, &
    ' primitive_error=',primitive_error,' conserved_error=',conserved_error, &
    ' primitive_ulps=',primitive_ulps,' conserved_ulps=',conserved_ulps
  if (parallel_rank==0) write(*,'(A)') 'DOMAIN DECOMPOSITION TESTS PASSED'
  call parallel_finalize()

contains

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if (.not.ok) call parallel_abort('domain test: '//label)
  end subroutine require

  subroutine check_partition()
    integer :: parts,r,first,last,i,coverage(37),h
    type(domain_layout) :: b,previous
    logical :: periodic
    do parts=1,5
      do h=1,3
        coverage=0
        do r=0,parts-1
          call build_domain([37,1,7],h,r,parts,.false.,b,message)
          call require(len_trim(message)==0,'valid unequal block')
          call validate_domain(b,message)
          call require(len_trim(message)==0,'canonical block')
          first=b%offset(1)+1
          last=first+b%local_shape(1)-1
          coverage(first:last)=coverage(first:last)+1
          do i=first,last
            call require(domain_owner(b,i)==r,'cell owner')
          end do
          if (r>0) call require(b%offset(1)==previous%offset(1)+previous%local_shape(1),'contiguous blocks')
          call require(all(b%active_axes.eqv.[.true.,.false.,.true.]),'global active axes')
          call require(all(b%owns_lower.eqv.[r==0,.true.,.true.]),'lower owner')
          call require(all(b%owns_upper.eqv.[r==parts-1,.true.,.true.]),'upper owner')
          call require(b%lower_rank==r-1,'lower neighbour')
          if (r<parts-1) then
            call require(b%upper_rank==r+1,'upper neighbour')
          else
            call require(b%upper_rank==-1,'outer boundary sentinel')
          end if
          previous=b
        end do
        call require(all(coverage==1),'domain covered exactly once')
      end do
      do r=0,parts-1
        periodic=.true.
        call build_domain([37,1,1],3,r,parts,periodic,b,message)
        call require(len_trim(message)==0,'periodic layout')
        call require(b%lower_rank==modulo(r-1,parts).and.b%upper_rank==modulo(r+1,parts),'wrapped neighbours')
        call require(.not.b%physical_lower(1).and..not.b%physical_upper(1),'periodic has no physical x boundary')
      end do
    end do
    call build_domain([12,1,1],3,0,4,.false.,b,message)
    call require(len_trim(message)==0.and.b%local_shape(1)==3,'minimum admissible block')
    call build_domain([11,1,1],3,0,4,.false.,b,message)
    call require(len_trim(message)>0,'insufficient halo rejected')
    call build_domain([4,1,1],1,0,4,.false.,b,message)
    call require(len_trim(message)>0,'active dimension cannot collapse locally')
    call build_domain([1,3,1],3,0,1,.false.,b,message)
    call require(len_trim(message)==0.and..not.b%active_axes(1),'serial degenerate x')
    call build_domain([1,3,1],3,0,2,.false.,b,message)
    call require(len_trim(message)>0,'cannot split degenerate x')
    call build_domain([0,1,1],3,0,1,.false.,b,message)
    call require(len_trim(message)>0,'empty global grid rejected')
    call build_domain([37,1,1],0,0,1,.false.,b,message)
    call require(len_trim(message)>0,'invalid halo rejected')
    call build_domain([37,1,1],3,0,0,.false.,b,message)
    call require(len_trim(message)>0,'zero ranks rejected')
    call build_domain([37,1,1],3,-1,4,.false.,b,message)
    call require(len_trim(message)>0,'negative rank rejected')
    call build_domain([37,1,1],3,4,4,.false.,b,message)
    call require(len_trim(message)>0,'out-of-range rank rejected')
    call build_domain([huge(1),1,1],3,0,1,.false.,b,message)
    call require(len_trim(message)>0,'integer extent overflow rejected')
    call build_domain([37,1,1],3,0,1,.false.,b,message)
    call require(domain_owner(b,0)==-1.and.domain_owner(b,38)==-1,'nonphysical index has no owner')
    b%local_shape(1)=36
    call validate_domain(b,message)
    call require(len_trim(message)>0,'tampered block rejected')
  end subroutine check_partition

  subroutine check_grid(kind,depth)
    integer, intent(in) :: kind,depth
    integer :: shape_global(3),i,j,k,g,offset,horizon_global,probe_global,expected_horizon,expected_probe
    real*8 :: grid_spacing(3),position_horizon,min_distance,distance,p_ulps,u_ulps
    real*8, allocatable :: xc(:),xf(:),yc(:),yf(:),zc(:),zf(:)
    real*8, allocatable :: centers(:,:,:,:),faces_x(:,:,:,:),faces_y(:,:,:,:),faces_z(:,:,:,:)
    real*8, allocatable :: primitives(:,:,:,:),conserved(:,:,:,:)
    logical :: periodic
    type(domain_layout) :: local

    shape_global=[37,1,1]
    geom_type='Cartesian'
    metric_type='Minkowski'
    case_name='Sod'
    use_log_r=.false.
    r_min=0.0d0; r_max=1.0d0
    y_min=0.2d0; y_max=1.1d0
    z_min=-0.5d0; z_max=0.5d0
    bh_mass=1.0d0; a_spin=0.0d0
    do_mdot_extraction=.false.
    do_numerical_diagnostics=.false.
    adb_idx=1.4d0
    g1=adb_idx/(adb_idx-1.0d0)
    nghost=depth
    periodic=.false.
    select case(kind)
    case(2) ! SRHD2D in the x-z layout used by KHI.
      shape_global=[37,1,7]; case_name='KHI'; periodic=.true.
      r_min=-0.5d0; r_max=0.5d0
    case(3) ! SRHD2D with the other transverse axis active.
      shape_global=[37,5,1]
    case(4) ! Curvilinear flat metric (jet geometry).
      shape_global=[37,1,7]; geom_type='Cylindrical'; r_min=0.1d0; r_max=4.0d0
    case(5) ! GRHD radial, linear EF.
      geom_type='Spherical'; metric_type='Eddington-Finkelstein'; r_min=1.2d0; r_max=12.0d0
    case(6) ! GRHD radial, logarithmic Kerr-Schild a=0.
      geom_type='Spheroidal'; metric_type='Kerr-Schild'; r_min=1.2d0; r_max=40.0d0; use_log_r=.true.
    case(7) ! Three-dimensional static Kerr cache and shared face geometry.
      shape_global=[37,5,6]; geom_type='Spheroidal'; metric_type='Kerr-Schild'
      r_min=1.2d0; r_max=40.0d0; use_log_r=.true.; a_spin=0.9d0
      y_min=0.0d0; y_max=pi; z_min=0.0d0; z_max=2.0d0*pi
    case(8) ! Excised horizon: only the global first cell owns the probe.
      geom_type='Spheroidal'; metric_type='Kerr-Schild'; r_min=3.0d0; r_max=12.0d0
    case(9) ! Horizon beyond domain, no valid probe.
      geom_type='Spherical'; metric_type='Eddington-Finkelstein'; r_min=0.2d0; r_max=1.0d0
    end select
    do_mdot_extraction=kind>=5.and.kind<=8
    call set_metric_type()
    nx=shape_global(1); ny=shape_global(2); nz=shape_global(3)
    call allocate_and_grid()
    call precalculate_metric_cache()
    call fill_state(kind)
    grid_spacing=[dx,dy,dz]
    xc=x(:); xf=x_face(:); yc=y(:); yf=y_face(:); zc=z(:); zf=z_face(:)
    allocate(centers(155,nx,ny,nz),faces_x(14,nx+1,ny,nz))
    do k=1,nz
      do j=1,ny
        do i=1,nx
          centers(:,i,j,k)=cell_cache(i,j,k)
        end do
        do i=1,nx+1
          faces_x(:,i,j,k)=[alpha_f_x(i,j,k),beta_f_x(:,i,j,k), &
            reshape(gamma_f_x(:,:,i,j,k),[9]),sqrt_gamma_f_x(i,j,k)]
        end do
      end do
    end do
    if (ny>1) then
      allocate(faces_y(14,nx,ny+1,nz))
      do k=1,nz
        do j=1,ny+1
          do i=1,nx
            faces_y(:,i,j,k)=[alpha_f_y(i,j,k),beta_f_y(:,i,j,k), &
              reshape(gamma_f_y(:,:,i,j,k),[9]),sqrt_gamma_f_y(i,j,k)]
          end do
        end do
      end do
    end if
    if (nz>1) then
      allocate(faces_z(14,nx,ny,nz+1))
      do k=1,nz+1
        do j=1,ny
          do i=1,nx
            faces_z(:,i,j,k)=[alpha_f_z(i,j,k),beta_f_z(:,i,j,k), &
              reshape(gamma_f_z(:,:,i,j,k),[9]),sqrt_gamma_f_z(i,j,k)]
          end do
        end do
      end do
    end if
    primitives=p(1:nx,1:ny,1:nz,:); conserved=u
    horizon_global=idx_horizon; probe_global=idx_probe
    ! Independent copy of the original serial locator, using the saved coordinates.
    expected_horizon=0; expected_probe=0
    if (trim(metric_type)/='Minkowski') then
      position_horizon=event_horizon_radius(bh_mass,a_spin)
      if (use_log_r) position_horizon=log(position_horizon)
      if (position_horizon<xf(1)) then
        expected_probe=1
      else if (position_horizon<=xf(size(xf))) then
        min_distance=huge(1.0d0)
        do i=1,nx
          distance=abs(xc(i+nghost+1)-position_horizon)
          if (distance<min_distance) then
            min_distance=distance; expected_horizon=i
          end if
          if (expected_probe==0.and.xc(i+nghost+1)>=position_horizon) expected_probe=i
        end do
      end if
    end if
    call require(horizon_global==expected_horizon.and.probe_global==expected_probe,'legacy serial horizon/probe')
    call release_grid()

    call build_domain(shape_global,depth,parallel_rank,parallel_size,periodic,local,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=shape_global(1); ny=shape_global(2); nz=shape_global(3)
    call allocate_and_grid(local)
    call precalculate_metric_cache()
    call fill_state(kind)
    offset=local%offset(1)
    call require(all([nx,ny,nz]==local%local_shape),'local physical allocation')
    call require(all([nx_global,ny_global,nz_global]==shape_global),'global sizes')
    call require(all(global_active_axes.eqv.(shape_global>1)),'active directions preserved')
    call require(all(owns_lower_boundary.eqv.local%physical_lower).and. &
      all(owns_upper_boundary.eqv.local%physical_upper),'physical boundary metadata')
    call require(all([dx,dy,dz]==grid_spacing),'same global spacing')
    call require(all(x==xc(offset+1:offset+size(x))),'all centers and ghost coordinates')
    call require(all(x_face==xf(offset+1:offset+size(x_face))),'all x faces')
    call require(all(y==yc).and.all(y_face==yf).and.all(z==zc).and.all(z_face==zf),'transverse coordinates')
    call require(all(shape(u)==[nx,ny,nz,neq]).and.all(shape(p)==[nx+2*nghost+1,ny+2*nghost+1, &
      nz+2*nghost+1,neq]),'FP64 state allocation layout')
    if (kind==2) then
      ! The unchanged initializer uses vector/scalar libm at different strip lengths.
      ! A pre-M02 manual allocation reproduces this roundoff without domain code.
      primitive_error=max(primitive_error,maxval(abs(p(1:nx,1:ny,1:nz,:)-primitives(offset+1:offset+nx,:,:,:))))
      conserved_error=max(conserved_error,maxval(abs(u-conserved(offset+1:offset+nx,:,:,:))))
      p_ulps=maxval(abs(p(1:nx,1:ny,1:nz,:)-primitives(offset+1:offset+nx,:,:,:))/ &
        spacing(max(abs(p(1:nx,1:ny,1:nz,:)),abs(primitives(offset+1:offset+nx,:,:,:)))))
      u_ulps=maxval(abs(u-conserved(offset+1:offset+nx,:,:,:))/ &
        spacing(max(abs(u),abs(conserved(offset+1:offset+nx,:,:,:)))))
      primitive_ulps=max(primitive_ulps,p_ulps)
      conserved_ulps=max(conserved_ulps,u_ulps)
      call require(p_ulps<=4.0d0,'KHI initial primitive roundoff <=4 ULP')
      call require(u_ulps<=16.0d0,'KHI initial conserved roundoff <=16 ULP')
      call require(all(p(1:nx,1:ny,1:nz,1:4)==primitives(offset+1:offset+nx,:,:,1:4)), &
        'KHI density, pressure and longitudinal velocities exact')
    else
      call require(all(p(1:nx,1:ny,1:nz,:)==primitives(offset+1:offset+nx,:,:,:)),'initial primitive field')
      call require(all(u==conserved(offset+1:offset+nx,:,:,:)),'initial conserved field')
    end if
    call require(idx_horizon_global==horizon_global.and.idx_probe_global==probe_global,'global horizon/probe stable')
    call require(horizon_owner_rank==domain_owner(local,horizon_global).and. &
      probe_owner_rank==domain_owner(local,probe_global),'probe ownership')
    if (horizon_owner_rank==parallel_rank) then
      call require(idx_horizon+offset==horizon_global,'local horizon')
    else
      call require(idx_horizon==0,'nonowner has no local horizon')
    end if
    if (probe_owner_rank==parallel_rank) then
      call require(idx_probe+offset==probe_global,'local probe')
    else
      call require(idx_probe==0,'nonowner has no local probe')
    end if
    do k=1,nz
      do j=1,ny
        do i=1,nx
          g=offset+i
          call require(all(cell_cache(i,j,k)==centers(:,g,j,k)),'full center cache')
        end do
        do i=1,nx+1
          g=offset+i
          call require(all([alpha_f_x(i,j,k),beta_f_x(:,i,j,k), &
            reshape(gamma_f_x(:,:,i,j,k),[9]),sqrt_gamma_f_x(i,j,k)]==faces_x(:,g,j,k)),'x face cache')
        end do
      end do
    end do
    if (ny>1) then
      do k=1,nz
        do j=1,ny+1
          do i=1,nx
            call require(all([alpha_f_y(i,j,k),beta_f_y(:,i,j,k), &
              reshape(gamma_f_y(:,:,i,j,k),[9]),sqrt_gamma_f_y(i,j,k)]==faces_y(:,offset+i,j,k)),'y face cache')
          end do
        end do
      end do
    end if
    if (nz>1) then
      do k=1,nz+1
        do j=1,ny
          do i=1,nx
            call require(all([alpha_f_z(i,j,k),beta_f_z(:,i,j,k), &
              reshape(gamma_f_z(:,:,i,j,k),[9]),sqrt_gamma_f_z(i,j,k)]==faces_z(:,offset+i,j,k)),'z face cache')
          end do
        end do
      end do
    end if
    call release_grid()
  end subroutine check_grid

  function cell_cache(i,j,k) result(values)
    integer, intent(in) :: i,j,k
    real*8 :: values(155)
    values=[alpha_c(i,j,k),beta_c(:,i,j,k),reshape(gamma_c(:,:,i,j,k),[9]), &
      reshape(gamma_inv_c(:,:,i,j,k),[9]),reshape(gmunu_c(:,:,i,j,k),[16]),sqrt_gamma_c(i,j,k), &
      reshape(chris_c(:,:,:,i,j,k),[64]),reshape(dg_c(:,:,:,i,j,k),[48]),dlna_c(:,i,j,k)]
  end function cell_cache

  subroutine fill_state(kind)
    integer, intent(in) :: kind
    integer :: i,j,k
    if (kind<=3) then
      call set_initial_conditions()
    else
      p=0.0d0
      do k=1,nz
        do j=1,ny
          do i=1,nx
            p(i,j,k,eq_de)=1.0d0+0.001d0*x(i)
            p(i,j,k,eq_pr)=0.1d0
          end do
        end do
      end do
    end if
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),u(i,j,k,:),i,j,k)
        end do
      end do
    end do
  end subroutine fill_state

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,s,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y)) deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z)) deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
  end subroutine release_grid
end program test_domain_decomposition
