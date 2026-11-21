program test_distributed_rk
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
  use evolution, only: update_adaptive_dt
  use time_integrator, only: advance_rk_step
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  integer, parameter :: fixture_count=10,config_count=10
  integer :: fixture,config,checks=0,cases=0,total_steps=0,unit,trajectory_unit
  character(len=64) :: mode,report_name
  real*8 :: field_error=0.0d0,dt_error=0.0d0,balance_error=0.0d0

  call parallel_initialize()
  call get_command_argument(1,mode)
  if (trim(mode)/=''.and.trim(mode)/='report') call check_invalid_call()
  if (trim(mode)=='report') then
    write(report_name,'(A,I0,A)') 'rk_cases_rank_',parallel_rank,'.txt'
    open(newunit=trajectory_unit,file=trim(report_name),status='new',action='write')
    write(trajectory_unit,'(A)') '# fixture config reconstruction solver nx ny nz steps events(16) '// &
      'final_time dt_error field_error min_rho min_p repair_delta_u(5)'
  end if
  do fixture=1,fixture_count
    do config=1,config_count
      call check_evolution()
    end do
  end do
  if (trim(mode)=='report') then
    close(trajectory_unit)
    write(report_name,'(A,I0,A)') 'rk_rank_',parallel_rank,'.txt'
    open(newunit=unit,file=trim(report_name),status='new',action='write')
    write(unit,'(5(I0,1X),3(ES24.16E3,1X))') parallel_rank,parallel_size,checks,cases,total_steps, &
      field_error,dt_error,balance_error
    close(unit)
  end if
  write(*,'(A,I0,A,I0,A,I0)') 'RK_TEST rank=',parallel_rank,' cases=',cases,' checks=',checks
  if (parallel_rank==0) write(*,'(A)') 'DISTRIBUTED RK TESTS PASSED'
  call parallel_finalize()

contains

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if (.not.ok) then
      write(*,*) 'RK fixture/config/rank:',fixture,config,parallel_rank
      call parallel_abort('RK test: '//label)
    end if
  end subroutine require

  subroutine configure(shape_global,periodic)
    integer, intent(out) :: shape_global(3)
    logical, intent(out) :: periodic
    integer, parameter :: methods(5)=[REC_GODUNOV,REC_TVD,REC_WENO3,REC_MP5,REC_WENO5]
    shape_global=[73,1,1]
    geom_type='Cartesian';metric_type='Minkowski';case_name='Conv'
    r_min=0.0d0;r_max=1.0d0;y_min=0.0d0;y_max=1.0d0;z_min=0.0d0;z_max=1.0d0
    use_log_r=.false.;a_spin=0.0d0;bh_mass=1.0d0;final_time=0.8d0;periodic=.true.
    select case(fixture)
    case(2,3)
      case_name='Sod';final_time=0.15d0;periodic=.false.
      if (fixture==3) then
        case_name='Strong';final_time=0.03d0
      end if
    case(4)
      case_name='KHI';shape_global(3)=12;r_min=-0.5d0;r_max=0.5d0
      z_min=-0.5d0;z_max=0.5d0;final_time=0.03d0
    case(5)
      case_name='Jet';geom_type='Cylindrical';shape_global(3)=12
      periodic=.false.;final_time=0.01d0
    case(6)
      case_name='MichelA';geom_type='Spherical';metric_type='Eddington-Finkelstein'
      r_min=1.5d0;r_max=40.0d0;use_log_r=.true.;final_time=0.2d0;periodic=.false.
    case(7:10)
      case_name='FontDaigne3D'
      if (fixture>=9) case_name='FishMonc3D'
      if (mod(fixture,2)==0) a_spin=0.9d0
      geom_type='Spheroidal';metric_type='Kerr-Schild'
      shape_global=[73,8,8];r_min=1.2d0;r_max=40.0d0;use_log_r=.true.
      if (fixture>=9) shape_global(2)=16
      y_max=pi;z_max=2.0d0*pi;final_time=0.1d0;periodic=.false.
    end select
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=3
    adb_idx=4.0d0/3.0d0
    if (fixture<=3) adb_idx=5.0d0/3.0d0
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
    if (fixture>=7) call require(count(p(1:nx,1:ny,1:nz,eq_de)>100.0d0*rho_floor)>0, &
      'torus fixture contains fluid above atmosphere')
    if (fixture==1) then
      do i=1,nx
        p(i,1,1,eq_de)=1.0d0+0.1d0*sin(2.0d0*pi*x(i))
        p(i,1,1,eq_pr)=1.0d0
        p(i,1,1,eq_vx:eq_vz)=[0.35d0,0.2d0,-0.15d0]
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

    call build_domain(shape_global,3,parallel_rank,parallel_size,periodic,layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
    call allocate_and_grid(layout);call precalculate_metric_cache()
    ! Identical physical seed isolates RK/halo equivalence from vectorized libm initialization.
    call prepare_case_boundary_data()
    p=-12345.0d0
    p(1:nx,1:ny,1:nz,:)=initial_p(x_global_offset+1:x_global_offset+nx,:,:,:)
    up=initial_u(x_global_offset+1:x_global_offset+nx,:,:,:)
    call initialize_halo_exchange(layout)
    call evolve(.true.,local_dt,local_steps)
    call require(local_steps==reference_steps,'common number of full RK steps')
    error=maxval(abs(local_dt(1:local_steps)-reference_dt(1:reference_steps)))
    case_dt_error=error
    dt_error=max(dt_error,error)
    call require(error==0.0d0,'common step sequence')
    error=max(maxval(abs(up-final_u(x_global_offset+1:x_global_offset+nx,:,:,:))), &
      maxval(abs(p(1:nx,1:ny,1:nz,:)-final_p(x_global_offset+1:x_global_offset+nx,:,:,:))))
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

  subroutine check_invalid_call()
    type(domain_layout) :: layout
    integer :: shape_global(3)
    logical :: periodic
    character(len=160) :: message
    fixture=1;config=1
    call configure(shape_global,periodic)
    call build_domain(shape_global,3,parallel_rank,parallel_size,periodic,layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    call allocate_and_grid(layout);call precalculate_metric_cache()
    p=1.0d0;up=1.0d0;dt=0.001d0
    call initialize_halo_exchange(layout)
    if (parallel_rank==parallel_size-1) then
      select case(trim(mode))
      case('rk-zero');dt=0.0d0
      case('rk-nan');dt=ieee_value(0.0d0,ieee_quiet_nan)
      case('rk-step');dt=0.002d0
      case('rk-state');deallocate(up)
      case('rk-halo');call finalize_halo_exchange()
      case('rk-shape');deallocate(up);allocate(up(nx-1,ny,nz,neq))
      case('rk-budget');do_radial_mass_budget=.true.
      end select
    end if
    if (trim(mode)=='rk-thread') then
      !$OMP PARALLEL
      call advance_rk_step(.true.)
      !$OMP END PARALLEL
    else
      call advance_rk_step(.true.)
    end if
    call parallel_abort('invalid RK call was not rejected')
  end subroutine check_invalid_call

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,s,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y)) deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z)) deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    if (allocated(michel_injector)) deallocate(michel_injector)
  end subroutine release_grid
end program test_distributed_rk
