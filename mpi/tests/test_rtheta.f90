program test_rtheta
  use parallel_runtime
  use parallel_timing, only: clip_step_to_target
  use domain_decomposition
  use halo_exchange
  use variables
  use metrics
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use conditions, only: set_initial_conditions,set_boundary_conditions,prepare_case_boundary_data
  use equations, only: prim_to_cons
  use fluxes, only: init_wavespeed_solver
  use evolution, only: update_adaptive_dt,calc_rhs_cpu,calc_rhs_accelerator, &
    initialize_rhs_accelerator,finalize_rhs_accelerator
  use time_integrator, only: advance_rk_step
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  integer, parameter :: fixture_count=8,config_count=10
  integer :: fixture,config,checks=0,cases=0,total_steps=0,unit,trajectory_unit,grid(2)=[1,1],depth
  character(len=32) :: argument
  character(len=64) :: mode,report_name
  real*8 :: field_error=0.0d0,dt_error=0.0d0,balance_error=0.0d0,halo_error=0.0d0,host_error=0.0d0

  call parallel_initialize()
  call get_command_argument(1,argument)
  if (len_trim(argument)>0) read(argument,*) grid(1)
  call get_command_argument(2,argument)
  if (len_trim(argument)>0) read(argument,*) grid(2)
  call get_command_argument(3,mode)
  if (trim(mode)/=''.and.trim(mode)/='report') call check_invalid_call()
  call check_partitions()
  if (trim(mode)=='report') then
    write(report_name,'(A,I0,A)') 'rtheta_cases_rank_',parallel_rank,'.txt'
    open(newunit=trajectory_unit,file=trim(report_name),status='new',action='write')
    write(trajectory_unit,'(A)') '# fixture config reconstruction solver nx ny nz steps events(16) '// &
      'final_time dt_error field_error min_rho min_p repair_delta_u(5)'
  end if
  do fixture=1,fixture_count
    do depth=1,3
      call check_transport()
    end do
    do config=1,config_count
      call check_evolution()
    end do
  end do
  if (trim(mode)=='report') then
    close(trajectory_unit)
    write(report_name,'(A,I0,A)') 'rtheta_rank_',parallel_rank,'.txt'
    open(newunit=unit,file=trim(report_name),status='new',action='write')
    write(unit,'(5(I0,1X),5(ES24.16E3,1X))') parallel_rank,parallel_size,checks,cases,total_steps, &
      field_error,dt_error,balance_error,halo_error,host_error
    close(unit)
  end if
  write(*,'(A,I0,A,I0,A,I0)') 'RTHETA_TEST rank=',parallel_rank,' cases=',cases,' checks=',checks
  if (parallel_rank==0) write(*,'(A)') 'RTHETA TESTS PASSED'
  call parallel_finalize()

contains

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if (.not.ok) then
      write(*,*) 'Rtheta fixture/config/rank:',fixture,config,parallel_rank
      call parallel_abort('Rtheta test: '//label)
    end if
  end subroutine require

  subroutine configure(shape_global,periodic)
    integer, intent(out) :: shape_global(3)
    logical, intent(out) :: periodic
    integer, parameter :: methods(5)=[REC_GODUNOV,REC_TVD,REC_WENO3,REC_MP5,REC_WENO5]
    shape_global=[37,19,8]
    geom_type='Cartesian';metric_type='Minkowski';case_name='Conv'
    r_min=0.0d0;r_max=1.0d0;y_min=0.0d0;y_max=1.0d0;z_min=0.0d0;z_max=1.0d0
    use_log_r=.false.;a_spin=0.0d0;bh_mass=1.0d0;final_time=0.03d0;periodic=.true.
    if (fixture==2) case_name='KHI'
    if (fixture>=3) then
      case_name='FontDaigne3D'
      if (fixture>=7) case_name='FishMonc3D'
      if (fixture==3) then
        case_name='FontDaigneSag';shape_global(3)=1
      end if
      if (fixture==4) then
        case_name='FishMoncSag';shape_global(3)=1
      end if
      if (fixture==4.or.fixture==6.or.fixture==8) a_spin=0.9d0
      geom_type='Spheroidal';metric_type='Kerr-Schild'
      r_min=1.2d0;r_max=40.0d0;use_log_r=.true.
      y_max=pi;z_max=2.0d0*pi;final_time=0.1d0;periodic=.false.
    end if
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=3
    adb_idx=4.0d0/3.0d0
    g1=adb_idx/(adb_idx-1.0d0);CFL=0.3d0
    rho_floor=1.0d-10;p_floor=1.0d-13
    D_floor=1.0d-10;tau_floor=3.0d-13
    do_mdot_extraction=.false.;do_gw_extraction=.false.;do_ppi_diagnostics=.false.
    do_radial_mass_budget=.false.;do_numerical_diagnostics=.true.
    shock_sensor_thresholds=SENSOR_DEFAULTS;use_shock_sensor=.true.
    rec_method_id=methods(mod(config-1,5)+1);tvd_limiter_id=LIM_MC
    riemann_solver_id=RS_HLLE
    if (config>5) riemann_solver_id=RS_HLLC
    call set_metric_type()
    call init_wavespeed_solver()
  end subroutine configure

  subroutine seed_full_problem()
    integer :: i,j,k
    call set_initial_conditions()
    if (fixture>=3) call require(count(p(1:nx,1:ny,1:nz,eq_de)>100.0d0*rho_floor)>0, &
      'torus fixture contains fluid above atmosphere')
    if (fixture==1) then
      do k=1,nz
        do j=1,ny
          do i=1,nx
            p(i,j,k,eq_de)=1.0d0+0.1d0*sin(2.0d0*pi*(x(i)+y(j)+z(k)))
            p(i,j,k,eq_pr)=1.0d0
            p(i,j,k,eq_vx:eq_vz)=[0.2d0,0.1d0,-0.05d0]
          end do
        end do
      end do
    end if
    call set_boundary_conditions(p)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
        end do
      end do
    end do
  end subroutine seed_full_problem

  subroutine evolve(distributed,history,steps)
    logical, intent(in) :: distributed
    real*8, intent(out) :: history(:)
    integer, intent(out) :: steps
    integration_time=0.0d0;n_steps=0
    numerical_events=0_int64;numerical_delta_u=0.0d0
    do while(integration_time<final_time)
      call update_adaptive_dt()
      dt=clip_step_to_target(integration_time,dt,final_time)
      n_steps=n_steps+1
      call require(n_steps<=size(history),'bounded test trajectory')
      history(n_steps)=dt
      integration_time=integration_time+dt
      call advance_rk_step(distributed)
    end do
    steps=n_steps
    call require(integration_time==final_time,'final physical time')
    call require(numerical_events(NUM_C2P_CALLS)==3_int64*steps*nx*ny*nz,'physical recovery count')
    call require(all(ieee_is_finite(up)).and.all(ieee_is_finite(p(1:nx,1:ny,1:nz,:))), &
      'finite final conserved and primitive state')
    call require(minval(p(1:nx,1:ny,1:nz,eq_de))>0.0d0.and. &
      minval(p(1:nx,1:ny,1:nz,eq_pr))>0.0d0,'positive final density and pressure')
  end subroutine evolve

  subroutine check_evolution()
    type(domain_layout) :: layout
    integer :: shape_global(3),reference_steps,local_steps,component
    logical :: periodic
    character(len=160) :: message
    real*8, allocatable :: initial_p(:,:,:,:),initial_u(:,:,:,:),final_p(:,:,:,:),final_u(:,:,:,:)
    real*8 :: reference_dt(4096),local_dt(4096),initial_sum(neq),local_sum(neq),global_sum(neq),error,scale,case_dt_error
#ifdef USE_MPI
    integer :: status
#endif
    call configure(shape_global,periodic)
    call allocate_and_grid();call precalculate_metric_cache();call seed_full_problem()
    initial_p=p(1:nx,1:ny,1:nz,:);initial_u=up
    do component=1,neq
      initial_sum(component)=sum(up(:,:,:,component))*dx*dy*dz
    end do
    call evolve(.false.,reference_dt,reference_steps)
    final_p=p(1:nx,1:ny,1:nz,:);final_u=up
    call release_grid()

    call build_domain(shape_global,3,parallel_rank,parallel_size,periodic,layout,message,grid,periodic)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
    call allocate_and_grid(layout);call precalculate_metric_cache()
    ! Identical physical seed isolates RK/halo equivalence from vectorized libm initialization.
    call prepare_case_boundary_data()
    p=-12345.0d0
    p(1:nx,1:ny,1:nz,:)=initial_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    up=initial_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    call initialize_halo_exchange(layout)
    call fill_primitive_halos(p)
    call initialize_rhs_accelerator(.true.)
    call calc_rhs_cpu(p,rhs);call calc_rhs_accelerator(p,s)
    error=maxval(abs(rhs-s))/max(1.0d0,maxval(abs(rhs)))
    host_error=max(host_error,error)
    call require(error<=2.0d-13,'host face/strip RHS with global polar indices')
    call finalize_rhs_accelerator(.true.)
    call evolve(.true.,local_dt,local_steps)
    call require(local_steps==reference_steps,'common number of full RK steps')
    error=maxval(abs(local_dt(1:local_steps)-reference_dt(1:reference_steps)))
    case_dt_error=error
    dt_error=max(dt_error,error)
    call require(error==0.0d0,'common step sequence')
    error=max(maxval(abs(up-final_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:))), &
      maxval(abs(p(1:nx,1:ny,1:nz,:)-final_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:))))
    field_error=max(field_error,error)
    call require(error==0.0d0,'distributed/full evolved conserved and primitive cells')
    if (trim(mode)=='report') then
      write(trajectory_unit,'(24(I0,1X),10(ES24.16E3,1X))') fixture,config,rec_method_id,riemann_solver_id, &
        nx,ny,nz,local_steps,numerical_events,integration_time,case_dt_error,error, &
        minval(p(1:nx,1:ny,1:nz,eq_de)),minval(p(1:nx,1:ny,1:nz,eq_pr)),numerical_delta_u
    end if
    if (fixture==1) then
      do component=1,neq
        local_sum(component)=sum(up(:,:,:,component))*dx*dy*dz
      end do
      global_sum=local_sum
#ifdef USE_MPI
      call MPI_Allreduce(local_sum,global_sum,neq,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call parallel_check_mpi(status,'test periodic five-variable integrals')
#endif
      scale=max(1.0d0,maxval(abs(initial_sum)))
      error=maxval(abs(global_sum-initial_sum))/scale
      balance_error=max(balance_error,error)
      call require(error<1.0d-11,'periodic five-variable conservation')
      call require(sum(numerical_events(NUM_C2P_LOW_D:NUM_P2C_VELOCITY_CAP))==0,'smooth wave without repairs')
    end if
    call finalize_halo_exchange();call release_grid()
    cases=cases+1;total_steps=total_steps+local_steps
  end subroutine check_evolution

  subroutine check_partitions()
    type(domain_layout) :: block,peer
    integer :: pr,pt,r,d,i,j,period,coverage(37,19)
    logical :: periodic
    character(len=160) :: message
    do depth=1,3
      do period=0,1
        periodic=period==1
        do pt=1,4
          do pr=1,4
            coverage=0
            do r=0,pr*pt-1
              call build_domain([37,19,8],depth,r,pr*pt,periodic,block,message,[pr,pt],periodic)
              call require(len_trim(message)==0,'valid unequal two-axis partition')
              call validate_domain(block,message)
              call require(len_trim(message)==0,'valid two-axis metadata')
              call require(all(block%offset(1:2)>=0).and.block%offset(3)==0,'nonnegative offsets and complete phi')
              call require(block%local_shape(3)==8,'complete azimuth')
              do j=1,block%local_shape(2)
                do i=1,block%local_shape(1)
                  coverage(block%offset(1)+i,block%offset(2)+j)=coverage(block%offset(1)+i,block%offset(2)+j)+1
                  call require(domain_owner(block,block%offset(1)+i,block%offset(2)+j)==r,'unique cell owner')
                end do
              end do
              call require(domain_owner(block,block%offset(1)+1)==r,'radial probe owner in this polar column')
              call require(domain_owner(block,0,1)==-1.and.domain_owner(block,1,20)==-1,'invalid global cell')
              do d=1,2
                if (block%lower_neighbor(d)>=0) then
                  call build_domain([37,19,8],depth,block%lower_neighbor(d),pr*pt,periodic,peer,message,[pr,pt],periodic)
                  call require(peer%upper_neighbor(d)==r,'reciprocal lower neighbor')
                end if
                if (block%upper_neighbor(d)>=0) then
                  call build_domain([37,19,8],depth,block%upper_neighbor(d),pr*pt,periodic,peer,message,[pr,pt],periodic)
                  call require(peer%lower_neighbor(d)==r,'reciprocal upper neighbor')
                end if
              end do
            end do
            call require(all(coverage==1),'complete nonoverlapping partition')
          end do
        end do
      end do
    end do
    call build_domain([37,5,8],3,0,2,.false.,block,message,[1,2])
    call require(index(message,'polar block too small')>0,'reject undersized polar blocks')
    call build_domain([37,19,8],3,0,4,.false.,block,message,[3,2])
    call require(index(message,'process grid')>0,'reject incompatible process grid')
  end subroutine check_partitions

  subroutine seed_transport(epoch)
    integer, intent(in) :: epoch
    integer :: i,j,k
    p=-12345.0d0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          p(i,j,k,eq_de)=1.0d0+epoch+0.01d0*(i+x_global_offset)+0.0001d0*(j+y_global_offset)+0.000001d0*k
          p(i,j,k,eq_pr)=0.1d0*p(i,j,k,eq_de)
          p(i,j,k,eq_vx)=-0.02d0+0.0001d0*(i+x_global_offset)
          p(i,j,k,eq_vy)=0.00001d0*(j+y_global_offset)
          p(i,j,k,eq_vz)=0.00001d0*k
        end do
      end do
    end do
  end subroutine seed_transport

  subroutine check_transport()
    type(domain_layout) :: block
    integer :: shape_global(3),i,j,k,gi,gj,epoch
    logical :: periodic
    character(len=160) :: message
    real*8, allocatable :: states(:,:,:,:,:),centers(:,:,:,:),faces_x(:,:,:,:),faces_y(:,:,:,:),faces_z(:,:,:,:)
    real*8, allocatable :: full_x(:),full_y(:),full_z(:),full_xf(:),full_yf(:),full_zf(:)
    real*8 :: error
    config=1
    call configure(shape_global,periodic)
    nghost=depth
    call allocate_and_grid();call precalculate_metric_cache()
    allocate(states(1-depth:nx+depth,1-depth:ny+depth,1-depth:nz+depth,neq,2))
    do epoch=1,2
      call seed_transport(epoch);call set_boundary_conditions(p)
      states(:,:,:,:,epoch)=p(1-depth:nx+depth,1-depth:ny+depth,1-depth:nz+depth,:)
    end do
    allocate(centers(155,nx,ny,nz),faces_x(14,nx+1,ny,nz),faces_y(14,nx,ny+1,nz))
    if (nz>1) allocate(faces_z(14,nx,ny,nz+1))
    do k=1,nz
      do j=1,ny
        do i=1,nx
          centers(:,i,j,k)=cell_cache(i,j,k)
        end do
        do i=1,nx+1
          faces_x(:,i,j,k)=[alpha_f_x(i,j,k),beta_f_x(:,i,j,k),reshape(gamma_f_x(:,:,i,j,k),[9]),sqrt_gamma_f_x(i,j,k)]
        end do
      end do
      do j=1,ny+1
        do i=1,nx
          faces_y(:,i,j,k)=[alpha_f_y(i,j,k),beta_f_y(:,i,j,k),reshape(gamma_f_y(:,:,i,j,k),[9]),sqrt_gamma_f_y(i,j,k)]
        end do
      end do
    end do
    if (nz>1) then
      do k=1,nz+1
        do j=1,ny
          do i=1,nx
            faces_z(:,i,j,k)=[alpha_f_z(i,j,k),beta_f_z(:,i,j,k),reshape(gamma_f_z(:,:,i,j,k),[9]),sqrt_gamma_f_z(i,j,k)]
          end do
        end do
      end do
    end if
    allocate(full_x(-depth:nx+depth),full_y(-depth:ny+depth),full_z(-depth:nz+depth))
    allocate(full_xf(0:nx),full_yf(0:ny),full_zf(0:nz))
    full_x=x;full_y=y;full_z=z;full_xf=x_face;full_yf=y_face;full_zf=z_face
    call release_grid()
    call build_domain(shape_global,depth,parallel_rank,parallel_size,periodic,block,message,grid,periodic)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
    call allocate_and_grid(block);call precalculate_metric_cache()
    call require(all(x==full_x(x_global_offset-depth:x_global_offset+nx+depth)),'global radial centers')
    call require(all(y==full_y(y_global_offset-depth:y_global_offset+ny+depth)),'global polar centers')
    call require(all(z==full_z),'complete azimuth centers')
    call require(all(x_face==full_xf(x_global_offset:x_global_offset+nx)),'global radial faces')
    call require(all(y_face==full_yf(y_global_offset:y_global_offset+ny)),'global polar faces')
    call require(all(z_face==full_zf),'complete azimuth faces')
    do k=1,nz
      do j=1,ny
        gj=j+y_global_offset
        do i=1,nx
          gi=i+x_global_offset
          call require(all(cell_cache(i,j,k)==centers(:,gi,gj,k)),'global cell metric cache')
        end do
        do i=1,nx+1
          gi=i+x_global_offset
          call require(all([alpha_f_x(i,j,k),beta_f_x(:,i,j,k),reshape(gamma_f_x(:,:,i,j,k),[9]), &
            sqrt_gamma_f_x(i,j,k)]==faces_x(:,gi,gj,k)),'radial face metric cache')
        end do
      end do
      do j=1,ny+1
        gj=j+y_global_offset
        do i=1,nx
          gi=i+x_global_offset
          call require(all([alpha_f_y(i,j,k),beta_f_y(:,i,j,k),reshape(gamma_f_y(:,:,i,j,k),[9]), &
            sqrt_gamma_f_y(i,j,k)]==faces_y(:,gi,gj,k)),'polar face metric cache')
        end do
      end do
    end do
    if (nz>1) then
      do k=1,nz+1
        do j=1,ny
          gj=j+y_global_offset
          do i=1,nx
            gi=i+x_global_offset
            call require(all([alpha_f_z(i,j,k),beta_f_z(:,i,j,k),reshape(gamma_f_z(:,:,i,j,k),[9]), &
              sqrt_gamma_f_z(i,j,k)]==faces_z(:,gi,gj,k)),'azimuth face metric cache')
          end do
        end do
      end do
    end if
    call initialize_halo_exchange(block)
    do epoch=1,2
      call seed_transport(epoch);call fill_primitive_halos(p)
      error=maxval(abs(p(1-depth:nx+depth,1-depth:ny+depth,1-depth:nz+depth,:)- &
        states(x_global_offset+1-depth:x_global_offset+nx+depth, &
        y_global_offset+1-depth:y_global_offset+ny+depth,:,:,epoch)))
      halo_error=max(halo_error,error)
      call require(error==0.0d0,'two-axis halos, corners, pole crossing and parities')
    end do
    call finalize_halo_exchange();call release_grid()
  end subroutine check_transport

  function cell_cache(i,j,k) result(values)
    integer, intent(in) :: i,j,k
    real*8 :: values(155)
    values=[alpha_c(i,j,k),beta_c(:,i,j,k),reshape(gamma_c(:,:,i,j,k),[9]), &
      reshape(gamma_inv_c(:,:,i,j,k),[9]),reshape(gmunu_c(:,:,i,j,k),[16]),sqrt_gamma_c(i,j,k), &
      reshape(chris_c(:,:,:,i,j,k),[64]),reshape(dg_c(:,:,:,i,j,k),[48]),dlna_c(:,i,j,k)]
  end function cell_cache

  subroutine check_invalid_call()
    type(domain_layout) :: layout
    integer :: shape_global(3),local_grid(2)
    logical :: periodic,periodic_y
    character(len=160) :: message
    fixture=1;config=1
    call configure(shape_global,periodic)
    periodic_y=periodic;local_grid=grid
    if (parallel_rank==parallel_size-1) then
      select case(trim(mode))
      case('bad-small');shape_global(2)=5
      case('bad-grid');local_grid=[1,1]
      case('bad-periodic');periodic_y=.false.
      case('bad-config');shape_global(2)=20
      case('bad-case');case_name='Jet';periodic=.false.;periodic_y=.false.
      end select
    end if
    if (trim(mode)=='bad-odd') then
      fixture=5;call configure(shape_global,periodic)
      shape_global(3)=7;periodic_y=.false.
    end if
    call build_domain(shape_global,3,parallel_rank,parallel_size,periodic,layout,message,local_grid,periodic_y)
    call parallel_require_collective(len_trim(message)==0,message)
    if (trim(mode)=='bad-layout'.and.parallel_rank==parallel_size-1) layout%offset(2)=layout%offset(2)+1
    call validate_domain(layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
    call allocate_and_grid(layout);call precalculate_metric_cache()
    call initialize_halo_exchange(layout)
    if (trim(mode)=='bad-metadata'.and.parallel_rank==parallel_size-1) y_global_offset=y_global_offset+1
    call fill_primitive_halos(p)
    call parallel_abort('invalid rtheta call was not rejected')
  end subroutine check_invalid_call

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,s,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y)) deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z)) deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    if (allocated(michel_injector)) deallocate(michel_injector)
  end subroutine release_grid
end program test_rtheta
