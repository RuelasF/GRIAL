program test_parallel_observer
  use variables
  use parallel_runtime
  use domain_decomposition
  use parallel_diagnostics
  use parallel_checkpoint
  use parallel_observer
  use parallel_timing,only:clip_step_to_target
  use initialization,only:allocate_and_grid,precalculate_metric_cache
  use conditions,only:set_initial_conditions,prepare_case_boundary_data
  use equations,only:prim_to_cons
  use fluxes,only:init_wavespeed_solver
  use metrics,only:set_metric_type
  use halo_exchange
  use evolution,only:update_adaptive_dt
  use time_integrator,only:advance_rk_step
  use output,only:distributed_output_cursor,export_distributed_output_cursor,restore_distributed_output_cursor
  use,intrinsic::iso_c_binding,only:c_int,c_long
  implicit none
  type,bind(C)::resource_limit
    integer(c_long)::soft,hard
  end type resource_limit
  interface
    integer(c_int) function get_file_limit(resource,limits)bind(C,name='getrlimit')
      import c_int,resource_limit
      integer(c_int),value::resource
      type(resource_limit),intent(out)::limits
    end function get_file_limit
    integer(c_int) function set_file_limit(resource,limits)bind(C,name='setrlimit')
      import c_int,resource_limit
      integer(c_int),value::resource
      type(resource_limit),intent(in)::limits
    end function set_file_limit
  end interface
  type(resource_limit)::previous_limit,injected_limit
  type(domain_layout)::layout
  type(observer_calendar)::calendar,saved_calendar
  type(conservation_balance)::before,after
  type(state_diagnostics)::a,b
  type(distributed_output_cursor)::cursor
  integer::grid(2)=[1,1],ios,i,j,k,unit,step
  real*8::time,step_size
  real*8,allocatable::saved_u(:,:,:,:),saved_p(:,:,:,:)
  character(len=512)::path,arg
  character(len=32)::mode
  character(len=160)::message
  call parallel_initialize()
  call get_command_argument(1,arg);read(arg,*)grid(1)
  call get_command_argument(2,arg);read(arg,*)grid(2)
  call get_command_argument(3,mode);call get_command_argument(4,path)
  if(trim(mode)=='reject'.or.trim(mode)=='accept')then
    call read_parallel_checkpoint_configuration(trim(path),layout,ios,grid)
    if(ios/=0)then
      call report('configuration_rejected');call parallel_finalize();stop
    end if
    call allocate_and_grid(layout);call precalculate_metric_cache();call prepare_case_boundary_data()
    call init_wavespeed_solver();call initialize_halo_exchange(layout)
    up=-12;u=-14;p=-13;numerical_events=43_int64;numerical_delta_u=7.0d0
    time=integration_time;step=n_steps;step_size=dt
    calendar=observer_calendar(next_save=-12.0d0,next_checkpoint=-13.0d0,save_index=42_int64)
    cursor=distributed_output_cursor(diagnostic_root='sentinel',diagnostic_step=456)
    call restore_distributed_output_cursor(cursor)
    call load_observed_checkpoint(trim(path),layout,calendar,ios)
    if(trim(mode)=='accept')then
      call require(ios==0.and.diagnostics_active,'previous observed generation remains recoverable')
      call report('recovered');call finalize_parallel_diagnostics()
      call finalize_halo_exchange();call parallel_finalize();stop
    end if
    call require(ios/=0,'corruption must be rejected')
    call require(all(up==-12.0d0).and.all(u==-14.0d0).and.all(p==-13.0d0),'rejection preserves every field and halo')
    call require(all(numerical_events==43_int64).and.all(numerical_delta_u==7.0d0),'rejection preserves raw histories')
    call require(integration_time==time.and.n_steps==step.and.dt==step_size,'rejection preserves clock')
    call require(.not.diagnostics_active.and.calendar%next_save==-12.0d0.and.calendar%save_index==42_int64, &
      'rejection preserves inactive observer and calendar')
    call export_distributed_output_cursor(cursor)
    if(parallel_rank==0)call require(cursor%diagnostic_root=='sentinel'.and.cursor%diagnostic_step==456, &
      'rejection preserves output cursor')
    call report('rejected_without_mutation')
  else
    call setup()
    if(trim(mode)=='partial')call begin_budget_step(0.001d0)
    if(trim(mode)=='active')call load_observed_checkpoint('missing',layout,calendar,ios)
    if(trim(mode)=='calendar-disagree'.and.parallel_rank==parallel_size-1)calendar%next_save=2
    if(trim(mode)=='path-disagree')then
      path='one'
      if(parallel_rank==parallel_size-1)path='two'
    end if
    if(trim(mode)=='status-disagree')then
      if(parallel_rank==0)then
        call save_observed_checkpoint('bad',layout,calendar,ios)
      else
        call save_observed_checkpoint('bad',layout,calendar)
      end if
      call parallel_abort('invalid status call accepted')
    end if
    if(len_trim(path)==0)path='observed'
    call initialize_observer_calendar(calendar)
    if(trim(mode)=='calendar-disagree'.and.parallel_rank==parallel_size-1)calendar%next_save=2
    if(trim(mode)/='io'.and.trim(mode)/='empty-produce')then
      call sample_parallel_observer(calendar,.true.)
      do i=1,5
        call update_adaptive_dt();dt=min(dt,0.001d0)
        n_steps=n_steps+1;integration_time=integration_time+dt
        call advance_rk_step(.true.);call sample_parallel_observer(calendar)
      end do
      call sample_parallel_observer(calendar,.true.)
      call sample_parallel_observer(calendar,.true.)
    end if
    call collect_conservation_balance(before);call collect_global_state(a)
    saved_u=up;saved_p=p(1:nx,1:ny,1:nz,:);saved_calendar=calendar
    if(trim(mode)=='io'.and.parallel_rank==0)then
      ! As in C01, apply the limit after PMIx creates its shared-memory files.
      call signal(25,1)
      call require(get_file_limit(1_c_int,previous_limit)==0,'read Linux file-size limit')
      injected_limit=previous_limit;injected_limit%soft=2000_c_long
      call require(set_file_limit(1_c_int,injected_limit)==0,'set root file-size limit')
    end if
    call save_observed_checkpoint(trim(path),layout,calendar,ios)
    if(trim(mode)=='io'.and.parallel_rank==0) &
      call require(set_file_limit(1_c_int,previous_limit)==0,'restore root file-size limit')
    if(trim(mode)=='empty-produce')then
      call require(ios==0,'empty observed generation saved');call report('empty_produced')
      call finalize_parallel_diagnostics();call finalize_halo_exchange();call parallel_finalize();stop
    end if
    if(trim(mode)=='io'.or.trim(mode)=='empty-produce')then
      call require(ios/=0,'real file limit must prevent observer publication')
      call report('io_rejected');call parallel_finalize();stop
    end if
    call require(ios==0,'observed generation saved')
    call save_observed_checkpoint(trim(path),layout,calendar,ios)
    call require(ios/=0,'immutable generation cannot be overwritten')
    call finalize_parallel_diagnostics()
    up=-12;u=-14;p=-13;numerical_events=0_int64;numerical_delta_u=0
    calendar=observer_calendar()
    call load_observed_checkpoint(trim(path),layout,calendar,ios)
    call require(ios==0,'observed generation recovered')
    call require(all(up==saved_u).and.all(u==saved_u).and.all(p(1:nx,1:ny,1:nz,:)==saved_p), &
      'roundtrip physical fields exact')
    call require(calendar%next_save==saved_calendar%next_save.and. &
      calendar%next_checkpoint==saved_calendar%next_checkpoint.and. &
      calendar%diagnostic_step==saved_calendar%diagnostic_step.and. &
      calendar%observable_step==saved_calendar%observable_step,'calendar preserved')
    call collect_conservation_balance(after);call collect_global_state(b)
    call require(all(before%initial==after%initial).and.all(before%flux==after%flux).and. &
      all(before%source==after%source).and.all(before%repair==after%repair).and. &
      all(before%perturbation==after%perturbation).and.all(before%residual==after%residual), &
      'budget history and residual exact at same partition')
    call require(all(a%events==b%events).and.all(a%raw_repair==b%raw_repair).and. &
      a%extra_interface_faces==b%extra_interface_faces,'counter history exact')
    call sample_parallel_observer(calendar,.true.)
    call report('roundtrip_exact')
  end if
  if(diagnostics_active)call finalize_parallel_diagnostics()
  call finalize_halo_exchange();call parallel_finalize()
contains
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    ! Called only on root in the output-cursor check, so do not use a collective here.
    if(.not.ok)call parallel_abort(label)
  end subroutine require
  subroutine report(label)
    character(len=*),intent(in)::label
    if(parallel_rank/=0)return
    open(newunit=unit,file='probe_result.txt',status='new')
    write(unit,'(A)')label;close(unit)
  end subroutine report
  subroutine setup()
    nx=37;ny=19;nz=8;nghost=3
    if(trim(mode)=='io')then
      nx=17;ny=1;nz=1
    end if
    case_name='Conv';case_id=6;geom_type='Cartesian';metric_type='Minkowski';bh_mass=0;a_spin=0
    r_min=0;r_max=1;y_min=0;y_max=1;z_min=0;z_max=1;use_log_r=.false.
    adb_idx=4.0d0/3.0d0;g1=adb_idx/(adb_idx-1);CFL=0.3d0
    rho_floor=2.0d-11;p_floor=7.0d-15;D_floor=rho_floor;tau_floor=p_floor*(g1-1)
    rec_method_id=REC_WENO5;riemann_solver_id=RS_HLLE;use_shock_sensor=.true.
    do_numerical_diagnostics=.true.;do_restart=.false.;do_radial_mass_budget=.false.;is_mhd=.false.
    do_gw_extraction=.false.;do_ppi_diagnostics=.false.;do_mdot_extraction=.false.
    apply_perturbation=.false.;perturbation_applied=.false.
    final_time=1;save_interval=0.1d0;checkpoint_interval=1;diagnostic_stride=3;extraction_stride=4
    output_folder='data';scheme_name='c03';output_prefix='c03'
    ios=0
    if(parallel_rank==0)call execute_command_line('mkdir -p data/c03',exitstat=ios)
    call parallel_require_collective(ios==0,'unit output directory')
    call build_domain([nx,ny,nz],nghost,parallel_rank,parallel_size,.true.,layout,message,grid,ny>1)
    call parallel_require_collective(len_trim(message)==0,message)
    call set_metric_type()
    call allocate_and_grid(layout);call precalculate_metric_cache();call init_wavespeed_solver()
    call set_initial_conditions()
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
        end do
      end do
    end do
    u=up;integration_time=0;n_steps=0;dt=0
    numerical_events=4294967296_int64+parallel_rank;numerical_delta_u=0
    call initialize_halo_exchange(layout);call fill_primitive_halos(p)
    call initialize_parallel_diagnostics(layout);call initialize_observer_calendar(calendar)
  end subroutine setup
end program test_parallel_observer
