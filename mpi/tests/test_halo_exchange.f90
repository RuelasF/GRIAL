program test_halo_exchange
  use parallel_runtime
  use domain_decomposition
  use halo_exchange
  use variables
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use metrics, only: set_metric_type
  use conditions, only: set_boundary_conditions
  implicit none

  integer :: kind,depth,alias,phase=0,checks=0,fixtures=0,unit
  character(len=48) :: mode,report_name
  real*8, parameter :: poison=-987654321.0d0
  real*8 :: maximum_error=0.0d0

  call parallel_initialize()
  call get_command_argument(1,mode)
  if (len_trim(mode)>0.and.trim(mode)/='report') call check_failure_mode()
#ifdef GRIAL_OVERLAP
  if (parallel_size==4) call check_polar_overlap()
#endif
  do depth=1,3
    do kind=1,15
      do alias=0,1
        call check_case(kind,depth,alias==1)
      end do
    end do
    if (parallel_size==1) call check_case(16,depth,.false.)
  end do
  if (trim(mode)=='report') then
    write(report_name,'(A,I0,A)') 'halo_rank_',parallel_rank,'.txt'
    open(newunit=unit,file=trim(report_name),status='new',action='write')
    write(unit,'(4(I0,1X),ES24.16E3)') parallel_rank,parallel_size,checks,fixtures,maximum_error
    close(unit)
  end if
  write(*,'(A,I0,A,I0,A,I0)') 'HALO_TEST rank=',parallel_rank,' size=',parallel_size,' checks=',checks
  if (parallel_rank==0) write(*,'(A)') 'HALO EXCHANGE TESTS PASSED'
  call parallel_finalize()

contains

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if (.not.ok) then
      write(*,*) 'HALO fixture/depth/alias/phase/rank:',kind,depth,alias,phase,parallel_rank
      call parallel_abort('halo test: '//label)
    end if
  end subroutine require

  real*8 function coded_value(global_i,j,k,m,epoch) result(value)
    integer, intent(in) :: global_i,j,k,m,epoch
    value=dble(100000*m+1000*global_i+10*j+k+epoch*1000000)/16777216.0d0
    if (m>=eq_vx.and.mod(global_i+j+k,2)==0) value=-value
    if (m==eq_de.and.mod(global_i,7)==0) value=0.1d0*rho_floor
    if (m==eq_pr.and.mod(global_i,11)==0) value=0.2d0*p_floor
  end function coded_value

  subroutine seed_physical(epoch)
    integer, intent(in) :: epoch
    integer :: i,j,k,m
    p=poison
    do m=1,neq
      do k=1,nz
        do j=1,ny
          do i=1,nx
            p(i,j,k,m)=coded_value(x_global_offset+i,j,k,m,epoch)
          end do
        end do
      end do
    end do
  end subroutine seed_physical

  subroutine configure(which,h,numeric,shape_global,periodic)
    integer, intent(in) :: which,h
    logical, intent(in) :: numeric
    integer, intent(out) :: shape_global(3)
    logical, intent(out) :: periodic
    character(len=20), parameter :: names(15)=[character(len=20) :: 'Sod','Strong','Shu','KHI','Jet','Conv', &
      'MichelA','Dust','OffAxis','FontDaigneEqu','FontDaigneSag','FontDaigne3D','FishMoncEqu','FishMoncSag','FishMonc3D']
    case_id=min(which,15)
    if (which==16) case_id=6
    case_name=names(case_id)
    if (numeric) write(case_name,'(I0)') case_id
    shape_global=[37,1,1]
    if (case_id==4.or.case_id==5.or.case_id==9) shape_global(3)=7
    if (case_id==10.or.case_id==13) shape_global(3)=8
    if (case_id==11.or.case_id==14) shape_global(2)=6
    if (case_id==12.or.case_id==15) shape_global=[37,6,8]
    if (which==16) shape_global=[1,1,1]
    periodic=case_id==4.or.case_id==6
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=h
    geom_type='Cartesian';metric_type='Minkowski';use_log_r=.false.
    r_min=0.0d0;r_max=1.0d0
    y_min=0.0d0;y_max=pi
    z_min=0.0d0;z_max=2.0d0*pi
    if (case_id==5) geom_type='Cylindrical'
    if (case_id>=7) then
      geom_type='Spherical';metric_type='Eddington-Finkelstein'
      r_min=1.2d0;r_max=12.0d0
    end if
    if (case_id>=10) then
      geom_type='Spheroidal';metric_type='Kerr-Schild';use_log_r=.true.
    end if
    bh_mass=1.0d0;a_spin=0.9d0
    do_mdot_extraction=.false.;do_numerical_diagnostics=.false.
    rho_floor=1.0d-10;p_floor=1.0d-13
    call set_metric_type()
  end subroutine configure

  subroutine prepare_injector(epoch,needed)
    integer, intent(in) :: epoch
    logical, intent(in) :: needed
    integer :: g,j,k,m
    if (allocated(michel_injector)) deallocate(michel_injector)
    if (case_id/=7.or..not.needed) return
    allocate(michel_injector(nghost,ny,nz,neq))
    do m=1,neq
      do k=1,nz
        do j=1,ny
          do g=1,nghost
            michel_injector(g,j,k,m)=coded_value(nx_global+g,j,k,m,epoch)
          end do
        end do
      end do
    end do
  end subroutine prepare_injector

  subroutine check_case(which,h,numeric)
    integer, intent(in) :: which,h
    logical, intent(in) :: numeric
    type(domain_layout) :: layout
    character(len=160) :: message
    integer :: shape_global(3),i,j,k,m,global_i,epoch
    logical :: periodic,has_data,started
    real*8 :: expected,error
    real*8, allocatable :: reference(:,:,:,:,:),baseline(:,:,:,:)

    call configure(which,h,numeric,shape_global,periodic)
    call allocate_and_grid()
    call precalculate_metric_cache()
    allocate(reference(1-h:nx+h,1-h:ny+h,1-h:nz+h,neq,2))
    do epoch=1,2
      call seed_physical(epoch)
      call prepare_injector(epoch,.true.)
      call set_boundary_conditions(p)
      reference(:,:,:,:,epoch)=p(1-h:nx+h,1-h:ny+h,1-h:nz+h,:)
    end do
    call release_grid()

    call build_domain(shape_global,h,parallel_rank,parallel_size,periodic,layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
    call allocate_and_grid(layout)
    call precalculate_metric_cache()
    call prepare_injector(1,layout%physical_upper(1))
    call initialize_halo_exchange(layout)
    do epoch=1,2
      phase=epoch
      call seed_physical(epoch)
      call exchange_primitive_halos(p)
      ! Check every allocated entry: only physical cells and received x slabs may change.
      do m=1,neq
        do k=-h,nz+h
          do j=-h,ny+h
            do i=-h,nx+h
              expected=poison
              if (j>=1.and.j<=ny.and.k>=1.and.k<=nz) then
                has_data=i>=1.and.i<=nx
                if (i>=1-h.and.i<=0.and.layout%lower_rank>=0) has_data=.true.
                if (i>nx.and.layout%upper_rank>=0) has_data=.true.
                if (has_data) then
                  global_i=x_global_offset+i
                  if (periodic) global_i=modulo(global_i-1,nx_global)+1
                  expected=coded_value(global_i,j,k,m,epoch)
                end if
              end if
              call require(p(i,j,k,m)==expected,'raw transport or untouched cell')
            end do
          end do
        end do
      end do
      call prepare_injector(epoch,layout%physical_upper(1))
      call fill_primitive_halos(p)
      error=maxval(abs(p(1-h:nx+h,1-h:ny+h,1-h:nz+h,:)- &
        reference(x_global_offset+1-h:x_global_offset+nx+h,:,:,:,epoch)))
      maximum_error=max(maximum_error,error)
      if (error/=0.0d0) then
        write(*,*) 'HALO difference:',error,' at ',maxloc(abs(p(1-h:nx+h,1-h:ny+h,1-h:nz+h,:)- &
          reference(x_global_offset+1-h:x_global_offset+nx+h,:,:,:,epoch)))
      end if
      call require(error==0.0d0,'all filled halos, edges and corners match serial boundaries')
      if (which==5.and.parallel_size>1) then
        allocate(baseline,source=p)
        call seed_physical(epoch)
        call begin_radial_overlap(p,started)
        call require(started,'jet radial overlap was not started')
        call end_radial_overlap(p)
        call require(all(p==baseline),'nonblocking jet halos differ from blocking halos')
        deallocate(baseline)
      end if
    end do
    call finalize_halo_exchange()
    call release_grid()
    fixtures=fixtures+1
  end subroutine check_case

#ifdef GRIAL_OVERLAP
  subroutine check_polar_overlap()
    type(domain_layout) :: layout
    character(len=160) :: message
    integer :: shape_global(3),h,model,epoch,grid_index
    logical :: periodic,started,polar_started
    real*8, allocatable :: baseline(:,:,:,:)
    do model=12,15,3
      do h=1,3
        do grid_index=1,2
          call configure(model,h,.false.,shape_global,periodic)
          shape_global=[37,16,8]
          if (grid_index==1) then
            call build_domain(shape_global,h,parallel_rank,parallel_size,.false.,layout,message,[2,2],.false.)
          else
            call build_domain(shape_global,h,parallel_rank,parallel_size,.false.,layout,message,[1,4],.false.)
          end if
          call parallel_require_collective(len_trim(message)==0,message)
          nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
          call allocate_and_grid(layout)
          call precalculate_metric_cache()
          call initialize_halo_exchange(layout)
          do epoch=1,2
            call seed_physical(epoch)
            call fill_primitive_halos(p)
            allocate(baseline,source=p)
            call seed_physical(epoch)
            call begin_radial_overlap(p,started)
            call require(started,'torus overlap not started')
            call end_radial_overlap(p,polar_started)
            call require(polar_started,'torus polar overlap not started')
            call finish_polar_overlap(p)
            call require(all(p==baseline),'nonblocking polar halos, edges or corners differ')
            deallocate(baseline)
          end do
          call finalize_halo_exchange()
          call release_grid()
        end do
      end do
    end do
  end subroutine check_polar_overlap
#endif

  subroutine check_failure_mode()
    type(domain_layout) :: layout
    character(len=160) :: message
    integer :: n_global
    logical :: periodic
    n_global=37
    if (trim(mode)=='bad-global-shape'.and.parallel_rank==parallel_size-1) n_global=38
    periodic=trim(mode)=='bad-periodic'
    nghost=3;case_name='Sod';ny=1;nz=1
    select case(trim(mode))
    case('bad-count')
      ny=huge(1)/2;nz=2
    case('bad-cache')
      case_name='FontDaigneEqu'
    case('bad-injector')
      case_name='MichelA'
    case('bad-polar')
      case_name='FishMonc3D';nz=7
    end select
    call build_domain([n_global,ny,nz],nghost,parallel_rank,parallel_size,periodic,layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=layout%local_shape(1)
    nx_global=n_global;ny_global=ny;nz_global=nz;x_global_offset=layout%offset(1)
    if (trim(mode)=='bad-polar') then
      allocate(alpha_c(nx,ny),beta_c(3,nx,ny))
      alpha_c=1.0d0;beta_c=0.0d0
    end if
    if (trim(mode)=='bad-layout'.and.parallel_rank==parallel_size-1) layout%offset(1)=layout%offset(1)+1
    call initialize_halo_exchange(layout)
    allocate(p(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq))
    p=poison
    select case(trim(mode))
    case('bad-array')
      if (parallel_rank==parallel_size-1) then
        call exchange_primitive_halos(p(-nghost:nx+nghost-1,:,:,:))
      else
        call exchange_primitive_halos(p)
      end if
    case('bad-thread')
      !$OMP PARALLEL
      call exchange_primitive_halos(p)
      !$OMP END PARALLEL
    case('bad-lifecycle')
      call finalize_halo_exchange()
      call exchange_primitive_halos(p)
    end select
    call parallel_abort('failure mode did not reject invalid input')
  end subroutine check_failure_mode

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y)) deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z)) deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    if (allocated(michel_injector)) deallocate(michel_injector)
  end subroutine release_grid
end program test_halo_exchange
