! Local validation driver: per-rank profiles, not production MPI output/restart.
program mpi_shock_profiles
  use variables
  use parallel_runtime
  use parallel_timing, only: time_decisions,coordinate_time_actions,clip_step_to_target
  use domain_decomposition
  use halo_exchange
  use metrics, only: set_metric_type
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use equations, only: prim_to_cons
  use conditions, only: set_initial_conditions
  use fluxes, only: init_wavespeed_solver
  use evolution, only: update_adaptive_dt
  use time_integrator, only: advance_rk_step
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(domain_layout) :: layout
  type(time_decisions) :: decisions
  character(len=64) :: problem,method,solver,n_argument,route
  character(len=512) :: prefix,path
  character(len=160) :: message
  integer :: global_n,i,global_i,status,trace_unit,metadata_unit
  real*8 :: left_rho,right_rho,left_p,right_p,fraction
  logical :: distributed,physical

  call parallel_initialize()
  call parallel_require_collective(command_argument_count()==6, &
    'usage: mpi_shock_profiles sod|strong|shu reconstruction solver nx prefix serial|mpi')
  call get_command_argument(1,problem);call get_command_argument(2,method)
  call get_command_argument(3,solver);call get_command_argument(4,n_argument)
  call get_command_argument(5,prefix);call get_command_argument(6,route)
  read(n_argument,*,iostat=status) global_n
  call parallel_require_collective(status==0,'invalid global shock resolution')
  call parallel_require_collective(global_n>=16,'shock resolution must be at least 16')
  call parallel_require_collective(trim(route)=='serial'.or.trim(route)=='mpi','invalid shock route')
  distributed=trim(route)=='mpi'
  call parallel_require_collective(distributed.or.parallel_size==1,'serial shock route requires one process')
  call parallel_require_collective(.not.accelerator_enabled,'shock comparison supports CPU only')

  case_name='Sod';adb_idx=1.4d0
  left_rho=1.0d0;right_rho=0.125d0;left_p=1.0d0;right_p=0.1d0
  select case(trim(problem))
  case('sod')
  case('strong')
    case_name='Strong';adb_idx=5.0d0/3.0d0
    right_rho=1.0d0;left_p=1000.0d0;right_p=0.01d0
  case('shu')
    case_name='Shu';adb_idx=5.0d0/3.0d0
    shu_interface=0.5d0;shu_rho_left=5.0d0;shu_pressure_left=50.0d0
    shu_rho_right=2.0d0;shu_rho_amplitude=0.3d0;shu_wave_number=50.0d0;shu_pressure_right=5.0d0
  case default
    call parallel_abort('unsupported shock case')
  end select
  select case(trim(method))
  case('godunov');rec_method_id=REC_GODUNOV
  case('tvd');rec_method_id=REC_TVD
  case('weno3');rec_method_id=REC_WENO3
  case('mp5');rec_method_id=REC_MP5
  case('weno5');rec_method_id=REC_WENO5
  case default;call parallel_abort('unsupported shock reconstruction')
  end select
  select case(trim(solver))
  case('hlle');riemann_solver_id=RS_HLLE
  case('hllc');riemann_solver_id=RS_HLLC
  case default;call parallel_abort('unsupported shock Riemann solver')
  end select
  nx=global_n;ny=1;nz=1;nghost=3
  r_min=0.0d0;r_max=1.0d0;y_min=0.0d0;y_max=1.0d0;z_min=0.0d0;z_max=1.0d0
  geom_type='Cartesian';metric_type='Minkowski';use_log_r=.false.;bh_mass=0.0d0;a_spin=0.0d0
  g1=adb_idx/(adb_idx-1.0d0);CFL=0.3d0;final_time=0.4d0
  rho_floor=1.0d-10;p_floor=1.0d-13;D_floor=1.0d-10;tau_floor=3.0d-13
  use_shock_sensor=.true.;shock_sensor_thresholds=SENSOR_DEFAULTS;tvd_limiter_id=LIM_MC
  do_numerical_diagnostics=.true.;do_radial_mass_budget=.false.
  do_mdot_extraction=.false.;do_gw_extraction=.false.;do_ppi_diagnostics=.false.
  call set_metric_type();call init_wavespeed_solver()
  call build_domain([global_n,1,1],nghost,parallel_rank,parallel_size,.false.,layout,message)
  call parallel_require_collective(len_trim(message)==0,message)
  call allocate_and_grid(layout);call precalculate_metric_cache()
  if (distributed) call initialize_halo_exchange(layout)
  p=0.0d0
  if (trim(problem)=='shu') call set_initial_conditions()
  do i=1,nx
    global_i=x_global_offset+i
    ! Exact cell fraction: for odd N the discontinuity cuts one cell in half.
    if (trim(problem)/='shu') then
      fraction=max(0.0d0,min(1.0d0,0.5d0*dble(global_n)-dble(global_i-1)))
      p(i,1,1,eq_de)=fraction*left_rho+(1.0d0-fraction)*right_rho
      p(i,1,1,eq_pr)=fraction*left_p+(1.0d0-fraction)*right_p
    end if
    call prim_to_cons(p(i,1,1,:),up(i,1,1,:),i,1,1)
  end do
  if (trim(problem)=='shu') call write_profile_file('_initial_rank_')
  integration_time=0.0d0;n_steps=0
  numerical_events=0_int64;numerical_delta_u=0.0d0
  write(path,'(A,A,I0,A)') trim(prefix),'_trace_rank_',parallel_rank,'.dat'
  open(newunit=trace_unit,file=trim(path),status='new',action='write',iostat=status)
  if (status/=0) call parallel_abort('cannot create shock time trace')
  write(trace_unit,'(A)') '# step actual_time dt'
  do
    call coordinate_time_actions(integration_time,n_steps,final_time,[0.0d0,0.0d0,0.0d0], &
      [.false.,.false.,.false.],[.false.,.false.],decisions)
    if (decisions%complete) exit
    call update_adaptive_dt()
    dt=clip_step_to_target(integration_time,dt,final_time)
    n_steps=n_steps+1;integration_time=integration_time+dt
    call advance_rk_step(distributed)
    write(trace_unit,'(I0,1X,2(ES25.17E3,1X))',iostat=status) n_steps,integration_time,dt
    if (status/=0) call parallel_abort('cannot write shock time trace')
  end do
  close(trace_unit,iostat=status)
  if (status/=0) call parallel_abort('cannot close shock time trace')
  call parallel_require_collective(all(ieee_is_finite(up)).and. &
    all(ieee_is_finite(p(1:nx,1:ny,1:nz,:))),'nonfinite evolved shock state')
  physical=all(p(1:nx,1,1,eq_de)>0.0d0).and.all(p(1:nx,1,1,eq_pr)>0.0d0)
  physical=physical.and.all(sum(p(1:nx,1,1,eq_vx:eq_vz)**2,dim=2)<1.0d0)
  call parallel_require_collective(physical,'nonphysical evolved shock state')
  call parallel_require_collective(numerical_events(NUM_C2P_CALLS)==3_int64*n_steps*nx, &
    'unexpected physical shock recovery count')

  call write_profile_file('_rank_')
  write(path,'(A,A,I0,A)') trim(prefix),'_metadata_rank_',parallel_rank,'.dat'
  open(newunit=metadata_unit,file=trim(path),status='new',action='write',iostat=status)
  if (status/=0) call parallel_abort('cannot create shock metadata')
  write(metadata_unit,'(8(I0,1X))',iostat=status) parallel_rank,parallel_size,global_n,x_global_offset, &
    nx,n_steps,rec_method_id,riemann_solver_id
  if (status==0) write(metadata_unit,'(8(ES25.17E3,1X))',iostat=status) &
    integration_time,dt,CFL,adb_idx,rho_floor,p_floor,D_floor,tau_floor
  if (status==0) write(metadata_unit,'(16(I0,1X))',iostat=status) numerical_events
  if (status==0) write(metadata_unit,'(5(ES25.17E3,1X))',iostat=status) numerical_delta_u
  if (status/=0) call parallel_abort('cannot write shock metadata')
  close(metadata_unit,iostat=status)
  if (status/=0) call parallel_abort('cannot close shock metadata')
  if (distributed) call finalize_halo_exchange()
  if (parallel_rank==0) write(*,'(A,A,A,I0,A,I0,A,ES25.17E3)') &
    'SHOCK_PROFILE PASSED case=',trim(problem),' ranks=',parallel_size,' steps=',n_steps,' time=',integration_time
  call parallel_finalize()

contains

  subroutine write_profile_file(label)
    character(len=*), intent(in) :: label
    integer :: cell,unit_file,io_status
    character(len=512) :: filename
    write(filename,'(A,A,I0,A)') trim(prefix),label,parallel_rank,'.dat'
    open(newunit=unit_file,file=trim(filename),status='new',action='write',iostat=io_status)
    if (io_status/=0) call parallel_abort('cannot create shock profile')
    write(unit_file,'(A)') '# global_i x rho pressure vx vy vz D tau Sx Sy Sz'
    do cell=1,nx
      write(unit_file,'(I0,1X,11(ES25.17E3,1X))',iostat=io_status) &
        x_global_offset+cell,x(cell),p(cell,1,1,:),up(cell,1,1,:)
      if (io_status/=0) call parallel_abort('cannot write shock profile')
    end do
    close(unit_file,iostat=io_status)
    if (io_status/=0) call parallel_abort('cannot close shock profile')
  end subroutine write_profile_file
end program mpi_shock_profiles
