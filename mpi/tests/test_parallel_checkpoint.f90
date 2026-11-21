program test_parallel_checkpoint
  use variables
  use parallel_runtime
  use parallel_checkpoint
  use parallel_diagnostics
  use parallel_perturbations, only: apply_parallel_perturbation
  use domain_decomposition
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use conditions
  use equations, only: prim_to_cons
  use fluxes, only: init_wavespeed_solver
  use metrics, only: set_metric_type
  use evolution, only: update_adaptive_dt
  use time_integrator, only: advance_rk_step
  use halo_exchange
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_c_binding, only: c_int,c_long
  implicit none
  type, bind(C) :: resource_limit
    integer(c_long) :: soft,hard
  end type resource_limit
  interface
    integer(c_int) function get_file_limit(resource,limits) bind(C,name='getrlimit')
      import c_int,resource_limit
      integer(c_int), value :: resource
      type(resource_limit), intent(out) :: limits
    end function get_file_limit
    integer(c_int) function set_file_limit(resource,limits) bind(C,name='setrlimit')
      import c_int,resource_limit
      integer(c_int), value :: resource
      type(resource_limit), intent(in) :: limits
    end function set_file_limit
  end interface
  type(domain_layout) :: layout
  integer :: grid(2)=[1,1],fixture=1,config=1,cases=0,checks=0,unit,ios,step_a,step_b,step_final,i
  integer(int64) :: events_a(NUM_EVENT_COUNT),events_b(NUM_EVENT_COUNT),events_final(NUM_EVENT_COUNT)
  real*8, allocatable :: ua(:,:,:,:),pa(:,:,:,:),ub(:,:,:,:),pb(:,:,:,:),uf(:,:,:,:),pf(:,:,:,:)
  real*8 :: time_a,time_b,time_final,dt_a,dt_b,dt_final,trace_reference(2),trace_resumed(2)
  real*8 :: delta_a(neq),delta_b(neq),delta_final(neq),max_error=0.0d0,budget_error=0.0d0,scale
  logical :: flag_a,flag_b,injected
  character(len=64) :: arg,mode,name
  character(len=512) :: target
  type(conservation_balance) :: balance

  call parallel_initialize()
  call get_command_argument(1,arg);if(len_trim(arg)>0)read(arg,*)grid(1)
  call get_command_argument(2,arg);if(len_trim(arg)>0)read(arg,*)grid(2)
  call get_command_argument(3,mode)
  call require(checkpoint_crc64(['1','2','3','4','5','6','7','8','9'])==int(z'6C40DF5F0B497347',int64), &
    'CRC64 standard check vector')
  if(trim(mode)/=''.and.trim(mode)/='report')call probe_mode()
  do fixture=1,13
    do config=1,merge(3,2,fixture>=6)
      call setup()
      call save_parallel_checkpoint(trim(output_folder)//'/generation_a',layout,ios)
      call require(ios==0,'initial generation committed')
      ua=up;pa=p;time_a=integration_time;step_a=n_steps;dt_a=dt
      events_a=numerical_events;delta_a=numerical_delta_u;flag_a=perturbation_applied
      call save_parallel_checkpoint(trim(output_folder)//'/generation_a',layout,ios)
      call require(ios/=0,'existing generation is preserved')
      do i=1,2
        call step_once()
      end do
      call save_parallel_checkpoint(trim(output_folder)//'/generation_b',layout,ios)
      call require(ios==0,'evolved generation committed')
      ub=up;pb=p;time_b=integration_time;step_b=n_steps;dt_b=dt
      events_b=numerical_events;delta_b=numerical_delta_u;flag_b=perturbation_applied
      do i=1,2
        call step_once();trace_reference(i)=dt
      end do
      uf=up;pf=p;time_final=integration_time;step_final=n_steps;dt_final=dt
      events_final=numerical_events;delta_final=numerical_delta_u
      call finalize_parallel_diagnostics()

      up=-12.0d0;u=up;p=-13.0d0;n_steps=99;integration_time=0.123d0;dt=0.001d0
      perturbation_applied=.not.flag_b;numerical_events=0_int64;numerical_delta_u=0
      call load_parallel_checkpoint(trim(output_folder)//'/generation_b',layout,ios)
      call require(ios==0,'same-partition generation recovered')
      call require(all(up==ub).and.all(u==ub).and.same_primitives(pb),'conserved, physical primitives and halos exact')
      call require(integration_time==time_b.and.n_steps==step_b.and.dt==dt_b,'epoch restored exactly')
      call require(all(numerical_events==events_b).and.all(numerical_delta_u==delta_b),'local counters restored exactly')
      call require(perturbation_applied.eqv.flag_b,'injection flag restored')
      call initialize_parallel_diagnostics(layout)
      call apply_parallel_perturbation(layout,injected)
      call require(.not.injected,'post-event recovery does not inject twice')
      do i=1,2
        call step_once();trace_resumed(i)=dt
      end do
      scale=max(1.0d0,maxval(abs(uf)))
      max_error=max(max_error,maxval(abs(up-uf))/scale)
      call require(all(up==uf).and.same_primitives(pf),'continued trajectory is bit-identical')
      call require(all(trace_reference==trace_resumed).and.integration_time==time_final.and. &
        n_steps==step_final.and.dt==dt_final,'continued step sequence is exact')
      call require(all(numerical_events==events_final),'continued physical counters are exact')
      call require(maxval(abs(numerical_delta_u-delta_final))<1.0d-11,'continued raw repair accumulator')
      call collect_conservation_balance(balance)
      scale=max(maxval(abs(balance%initial)),maxval(abs(balance%current)),1.0d-200)
      budget_error=max(budget_error,maxval(abs(balance%residual))/scale)
      call require(maxval(abs(balance%residual))/scale<1.0d-11,'new explicit D01 baseline closes')
      call finalize_parallel_diagnostics()

      call load_parallel_checkpoint(trim(output_folder)//'/generation_a',layout,ios)
      call require(ios==0.and.all(up==ua).and.same_primitives(pa),'older generation remains recoverable')
      call require(integration_time==time_a.and.n_steps==step_a.and.dt==dt_a,'older epoch restored')
      call require(all(numerical_events==events_a).and.all(numerical_delta_u==delta_a),'older counters intact')
      call require(perturbation_applied.eqv.flag_a,'pre-event flag restored')
      call teardown()
      deallocate(ua,pa,ub,pb,uf,pf)
      cases=cases+1
    end do
  end do
  write(name,'(A,I0,A)')'checkpoint_rank_',parallel_rank,'.txt'
  open(newunit=unit,file=trim(name),status='new')
  write(unit,'(4(I0,1X),2(ES24.16E3,1X))')parallel_rank,parallel_size,checks,cases,max_error,budget_error
  close(unit)
  if(parallel_rank==0)print *,'PARALLEL CHECKPOINT TESTS PASSED'
  call parallel_finalize()

contains
  logical function same_primitives(expected) result(same)
    real*8, intent(in) :: expected(:,:,:,:)
    ! Index -nghost is an unused extra allocation plane; the stencil starts at 1-nghost.
    same=all(p(1-nghost:,1-nghost:,1-nghost:,:)==expected(2:,2:,2:,:))
  end function same_primitives

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if(.not.ok)call parallel_abort('C01 test: '//label)
  end subroutine require

  subroutine setup()
    integer :: shape_global(3),parts(2),i,j,k
    logical :: periodic
    character(len=160) :: message
    shape_global=[37,19,8];parts=grid;periodic=.true.
    case_name='Conv';case_id=6;geom_type='Cartesian';metric_type='Minkowski';a_spin=0;bh_mass=1
    r_min=0;r_max=1;y_min=0;y_max=1;z_min=0;z_max=1;use_log_r=.false.
    select case(fixture)
    case(2)
      case_name='Sod';case_id=1;shape_global=[73,1,1];periodic=.false.
    case(3)
      case_name='KHI';case_id=4;shape_global=[37,1,19];z_min=-0.5d0;z_max=0.5d0
    case(4)
      case_name='Jet';case_id=5;shape_global=[37,1,19];periodic=.false.
      geom_type='Cylindrical';r_max=4;z_max=10
    case(5)
      case_name='MichelA';case_id=7;shape_global=[73,1,1];periodic=.false.
      geom_type='Spherical';metric_type='Eddington-Finkelstein';r_min=1.2d0;r_max=40;use_log_r=.true.
      michel_critical_radius=8;michel_critical_density=1
    case(6:13)
      periodic=.false.;case_name='FontDaigne3D';case_id=12
      if(mod((fixture-6)/2,2)==1)then
        case_name='FishMonc3D';case_id=15
      end if
      if(mod(fixture-6,2)==1)a_spin=0.9d0
      if(fixture>=10)then
        shape_global(3)=1
        if(case_id==12)then
          case_name='FontDaigneSag';case_id=11
        else
          case_name='FishMoncSag';case_id=14
        end if
      end if
      geom_type='Spheroidal';metric_type='Kerr-Schild';r_min=1.2d0;r_max=40
      y_max=pi;z_max=2*pi;use_log_r=.true.
    end select
    if(shape_global(2)==1)parts=[parallel_size,1]
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=3
    adb_idx=4.0d0/3.0d0;g1=adb_idx/(adb_idx-1);CFL=0.3d0
    rho_floor=2.0d-11;p_floor=7.0d-15;D_floor=rho_floor;tau_floor=p_floor*(g1-1)
    rec_method_id=REC_WENO5;riemann_solver_id=RS_HLLE;tvd_limiter_id=LIM_MC
    if(config==2)then
      rec_method_id=REC_MP5;riemann_solver_id=RS_HLLC
    end if
    use_shock_sensor=.true.;shock_sensor_thresholds=SENSOR_DEFAULTS
    do_numerical_diagnostics=.true.;do_restart=.false.;do_radial_mass_budget=.false.;is_mhd=.false.
    do_gw_extraction=.false.;do_ppi_diagnostics=.false.;do_mdot_extraction=.false.
    apply_perturbation=fixture>=6;perturbation_applied=.false.;perturbation_type=config
    perturbation_amplitude=0.01d0;perturbation_seed=3435;perturbation_mode=2;perturbation_time=0
    final_time=1;save_interval=1;checkpoint_interval=2;diagnostic_stride=5;extraction_stride=7;mass_monitor_stride=0
    walltime_limit_seconds=0;walltime_margin_seconds=300;walltime_checkpoint_seconds=60
    output_prefix='c01';scheme_name='c01'
    write(output_folder,'(A,I0,A,I0)')'case_',fixture,'_',config
    ios=0
    if(parallel_rank==0)call execute_command_line('mkdir -p '//trim(output_folder)//'/c01',exitstat=ios)
    call parallel_require_collective(ios==0,'C01 test directory prepared')
    call set_metric_type();call init_wavespeed_solver()
    call build_domain(shape_global,nghost,parallel_rank,parallel_size,periodic,layout,message,parts,periodic.and.ny>1)
    call parallel_require_collective(len_trim(message)==0,message)
    call allocate_and_grid(layout);call precalculate_metric_cache();call prepare_case_boundary_data()
    call set_initial_conditions()
    if(fixture>=6)then
      ! Some blocks contain only atmosphere; require the body only globally via the maximum.
      call require(maxval(p(1:nx,1:ny,1:nz,eq_de))>0,'positive torus or atmosphere')
    end if
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
        end do
      end do
    end do
    u=up;integration_time=0;n_steps=0;dt=0
    numerical_events=[(4294967296_int64+int(parallel_rank+i,int64),i=1,NUM_EVENT_COUNT)]
    numerical_delta_u=[(dble(parallel_rank+i)*1.0d-12,i=1,neq)]
    call initialize_halo_exchange(layout);call fill_primitive_halos(p)
    call initialize_parallel_diagnostics(layout)
  end subroutine setup

  subroutine step_once()
    call update_adaptive_dt();dt=min(dt,0.002d0)
    call apply_parallel_perturbation(layout,injected)
    n_steps=n_steps+1;integration_time=integration_time+dt
    call advance_rk_step(.true.)
  end subroutine step_once

  subroutine teardown()
    if(diagnostics_active)call finalize_parallel_diagnostics()
    call finalize_halo_exchange()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,s,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if(allocated(alpha_f_y))deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if(allocated(alpha_f_z))deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    if(allocated(michel_injector))deallocate(michel_injector)
  end subroutine teardown

  subroutine probe_mode()
    integer :: expected_status
    real*8, allocatable :: original(:,:,:,:),original_p(:,:,:,:)
    integer(int64) :: counters(NUM_EVENT_COUNT)
    real*8 :: original_delta(neq),original_time,original_dt
    integer :: original_step
    logical :: original_flag
    type(resource_limit) :: previous_limit,injected_limit
    fixture=1;config=1;call setup()
    call get_command_argument(4,target)
    select case(trim(mode))
    case('read-failure','read-success')
      call finalize_parallel_diagnostics()
      original=up;original_p=p;counters=numerical_events;original_delta=numerical_delta_u
      original_time=integration_time;original_dt=dt;original_step=n_steps;original_flag=perturbation_applied
      call load_parallel_checkpoint(trim(target),layout,ios)
      expected_status=1
      if(trim(mode)=='read-success')expected_status=0
      call require(merge(0,1,ios==0)==expected_status,'expected recovery status')
      if(expected_status/=0)then
        call require(all(up==original).and.all(p==original_p),'failed recovery preserves every physical cell and ghost')
        call require(all(numerical_events==counters).and.all(numerical_delta_u==original_delta),'failed recovery preserves counters')
        call require(integration_time==original_time.and.dt==original_dt.and.n_steps==original_step.and. &
          (perturbation_applied.eqv.original_flag),'failed recovery preserves epoch and flag')
      end if
      if(parallel_rank==0)print *,'CHECKPOINT RECOVERY PROBE PASSED'
    case('write-failure','write-crash')
      ! gfortran installs its own startup handler; ignore Linux SIGXFSZ only for this returned-error probe.
      if(trim(mode)=='write-failure')call signal(25,1)
      ! Set RLIMIT_FSIZE after MPI initializes its own shared-memory files.
      if(parallel_rank==parallel_size-1)then
        call require(get_file_limit(1_c_int,previous_limit)==0,'read Linux file-size resource limit')
        injected_limit=previous_limit;injected_limit%soft=4096_c_long
        call require(set_file_limit(1_c_int,injected_limit)==0,'inject child-only file-size limit')
      end if
      call save_parallel_checkpoint(trim(target),layout,ios)
      if(parallel_rank==parallel_size-1)call require(set_file_limit(1_c_int,previous_limit)==0,'restore resource limit')
      call require(ios/=0,'injected write failure is reported')
      if(trim(mode)=='write-crash')call parallel_abort('expected process death did not occur')
      if(parallel_rank==0)print *,'CHECKPOINT IO FAILURE PROBE PASSED'
    case default
      target='invalid_generation'
      select case(trim(mode))
      case('bad-layout');if(parallel_rank==parallel_size-1)layout%offset(1)=layout%offset(1)+1
      case('bad-path');if(parallel_rank==parallel_size-1)target='different_generation'
      case('bad-clock');if(parallel_rank==parallel_size-1)integration_time=0.1d0
      case('bad-seed');if(parallel_rank==parallel_size-1)perturbation_seed=99
      case('bad-flag');if(parallel_rank==parallel_size-1)perturbation_applied=.true.
      case('bad-state');if(parallel_rank==parallel_size-1)up(1,1,1,eq_de)=ieee_value(0.0d0,ieee_quiet_nan)
      case('bad-counters');if(parallel_rank==parallel_size-1)numerical_events(1)=-1_int64
      case('bad-eos');if(parallel_rank==parallel_size-1)g1=1.5d0
      case('bad-midstep');call begin_budget_step(0.01d0)
      case('bad-load-active');call load_parallel_checkpoint(trim(target),layout)
      case('bad-no-halos')
        call finalize_parallel_diagnostics();call finalize_halo_exchange()
        call load_parallel_checkpoint(trim(target),layout)
      case default;call parallel_abort('unknown C01 probe mode')
      end select
      call save_parallel_checkpoint(trim(target),layout)
      call parallel_abort('invalid C01 request accepted')
    end select
    call teardown();call parallel_finalize()
    stop
  end subroutine probe_mode
end program test_parallel_checkpoint
