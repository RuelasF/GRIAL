program test_parallel_restart
  use variables
  use parallel_runtime
  use parallel_checkpoint
  use parallel_diagnostics
  use parallel_perturbations, only: apply_parallel_perturbation
  use domain_decomposition
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use conditions
  use equations, only: prim_to_cons,cons_to_prim
  use fluxes, only: init_wavespeed_solver
  use metrics, only: set_metric_type
  use evolution, only: update_adaptive_dt
  use time_integrator, only: advance_rk_step
  use halo_exchange
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use output, only: save_checkpoint,read_checkpoint_metadata,read_checkpoint_state
  use initialization, only: choose_numerical_architecture
  use parameters, only: validate_parameters
  implicit none
  type(domain_layout) :: layout
  integer :: grid(2)=[1,1],fixture=1,config=1,cases=0,checks=0,unit,ios,i,j,k,parts(2),target_step=8,version
  integer :: saved_step,actual_step
  real*8 :: max_error=0.0d0,max_balance=0.0d0,actual_time,actual_dt,scale,legacy_time
  real*8, allocatable :: actual_u(:,:,:,:),actual_p(:,:,:,:)
  integer(int64) :: old_events(NUM_EVENT_COUNT)
  real*8 :: old_delta(neq)
  logical :: injected,actual_flag
  character(len=64) :: arg,mode,source_label,dest_label,reference_prefix,name
  character(len=256) :: case_dir
  character(len=512) :: source_path,dest_path,reference_path
  type(state_diagnostics) :: before_state,after_state,reference_state
  type(conservation_balance) :: balance

  call parallel_initialize()
  call get_command_argument(1,arg);if(len_trim(arg)>0)read(arg,*)grid(1)
  call get_command_argument(2,arg);if(len_trim(arg)>0)read(arg,*)grid(2)
  call get_command_argument(3,mode)
  call get_command_argument(4,source_label)
  call get_command_argument(5,dest_label)
  call get_command_argument(6,arg);if(len_trim(arg)>0)read(arg,*)target_step
  call get_command_argument(7,reference_prefix)
  if(len_trim(reference_prefix)==0)reference_prefix='ref'
  if(trim(mode)=='probe'.or.trim(mode)=='serial-probe')call probe()
  do fixture=1,13
    do config=1,merge(3,2,fixture>=6)
      write(case_dir,'(A,I0,A,I0)')'case_',fixture,'_',config
      parts=grid
      if(fixture>=2.and.fixture<=5)parts=[parallel_size,1]
      select case(trim(mode))
      case('produce')
        call require(parallel_size==1,'reference production uses one rank')
        call setup();perturbation_time=1.0d0
        output_folder=trim(case_dir)//'/events_ref'
        do while(n_steps<8)
          call step_once()
          if(n_steps==2)perturbation_time=integration_time
          if(mod(n_steps,2)==0)then
            write(name,'(A,I0)')'ref_',n_steps
            call save_parallel_checkpoint(trim(case_dir)//'/'//trim(name),layout,ios)
            call require(ios==0,'reference checkpoint')
            call collect_global_state(after_state)
            call write_history(trim(name),after_state,after_state)
            if(n_steps==2)call save_checkpoint(n_steps,integration_time,ios)
            call require(ios==0,'serial V9 source checkpoint')
          end if
        end do
      case('serial-reference')
        call require(parallel_size==1,'legacy reference uses one rank')
        source_path=trim(case_dir)//'/'//trim(source_label)
        call legacy_reference_setup(source_path)
        output_folder=trim(case_dir)//'/events_'//trim(reference_prefix)
        do_numerical_diagnostics=.true.;call initialize_parallel_diagnostics(layout)
        write(name,'(A,A,I0)')trim(reference_prefix),'_',n_steps
        call save_parallel_checkpoint(trim(case_dir)//'/'//trim(name),layout,ios)
        call require(ios==0,'legacy startup reference checkpoint')
        do while(n_steps<8)
          call step_once()
          if(mod(n_steps,2)==0)then
            write(name,'(A,A,I0)')trim(reference_prefix),'_',n_steps
            call save_parallel_checkpoint(trim(case_dir)//'/'//trim(name),layout,ios)
            call require(ios==0,'legacy reference checkpoint')
          end if
        end do
      case('resume','import')
        call poison_configuration()
        source_path=trim(case_dir)//'/'//trim(source_label)
        if(trim(mode)=='resume')then
          call read_parallel_checkpoint_configuration(source_path,layout,ios,parts)
        else
          call read_serial_checkpoint_configuration(source_path,layout,ios,parts)
        end if
        call require(ios==0,'checkpoint configuration loaded')
        call allocate_and_grid(layout);call precalculate_metric_cache();call prepare_case_boundary_data()
        call init_wavespeed_solver();call initialize_halo_exchange(layout)
        up=-12;u=-12;p=-13
        if(trim(mode)=='resume')then
          call load_parallel_checkpoint(source_path,layout,ios,.true.)
        else
          call import_serial_checkpoint(source_path,layout,ios)
        end if
        call require(ios==0,'checkpoint state loaded')
        do_numerical_diagnostics=.true.
        call initialize_parallel_diagnostics(layout);call collect_global_state(before_state)
        output_folder=trim(case_dir)//'/events_'//trim(dest_label)
        do while(n_steps<target_step)
          call step_once()
        end do
        call collect_global_state(after_state);call collect_conservation_balance(balance)
        scale=max(maxval(abs(balance%initial)),maxval(abs(balance%current)),1.0d-200)
        max_balance=max(max_balance,maxval(abs(balance%residual))/scale)
        call require(maxval(abs(balance%residual))/scale<1.0d-11,'new D01 baseline conservation')
        dest_path=trim(case_dir)//'/'//trim(dest_label)
        call save_parallel_checkpoint(dest_path,layout,ios)
        call require(ios==0,'repartitioned generation saved')
        call write_history(dest_label,before_state,after_state)
        actual_u=up;actual_p=p(1:nx,1:ny,1:nz,:)
        actual_time=integration_time;actual_dt=dt;actual_step=n_steps;actual_flag=perturbation_applied
        call finalize_parallel_diagnostics()
        output_folder=trim(case_dir)//'/events_'//trim(reference_prefix)
        write(reference_path,'(A,A,A,A,I0)')trim(case_dir),'/',trim(reference_prefix),'_',target_step
        call load_parallel_checkpoint(trim(reference_path),layout,ios,.true.)
        call require(ios==0,'independent full-grid reference loaded')
        scale=max(1.0d0,maxval(abs(up)),maxval(abs(p(1:nx,1:ny,1:nz,:))))
        max_error=max(max_error,maxval(abs(up-actual_u))/scale,maxval(abs(p(1:nx,1:ny,1:nz,:)-actual_p))/scale)
        call require(all(up==actual_u).and.all(p(1:nx,1:ny,1:nz,:)==actual_p),'repartitioned continuation matches reference exactly')
        call require(integration_time==actual_time.and.n_steps==actual_step.and.dt==actual_dt.and. &
          (perturbation_applied.eqv.actual_flag),'same endpoint, step and injection flag')
        call initialize_parallel_diagnostics(layout);call collect_global_state(reference_state)
        call require(all(after_state%events(8:14)==reference_state%events(8:14)),'physical counters equal serial reference')
        call require(maxval(abs(after_state%raw_repair-reference_state%raw_repair))<1.0d-11,'raw repair totals equal reference')
        deallocate(actual_u,actual_p)
      case default
        call parallel_abort('unknown C02 test mode')
      end select
      call teardown();cases=cases+1
    end do
  end do
  write(name,'(A,A,A,I0,A)')'restart_',trim(mode),'_rank_',parallel_rank,'.txt'
  open(newunit=unit,file=trim(name),status='new')
  write(unit,'(4(I0,1X),2(ES24.16E3,1X))')parallel_rank,parallel_size,checks,cases,max_error,max_balance
  close(unit)
  if(parallel_rank==0)print *,'PARALLEL RESTART TESTS PASSED'
  call parallel_finalize()
contains
  subroutine probe()
    real*8,allocatable::initial_u(:,:,:,:),initial_p(:,:,:,:)
    integer(int64)::initial_events(NUM_EVENT_COUNT)
    real*8::initial_delta(neq),initial_time,initial_dt
    integer::initial_step
    logical::initial_flag
    character::byte
    integer::file,close_status
    integer(int64)::length
    if(trim(mode)=='serial-probe')then
      call read_serial_checkpoint_configuration(trim(source_label),layout,ios,grid)
    else
      call read_parallel_checkpoint_configuration(trim(source_label),layout,ios,grid)
    end if
    if(trim(dest_label)=='config-failure')then
      call require(ios/=0.and..not.allocated(up),'bad configuration rejected before allocation')
    else
      call require(ios==0,'probe configuration loaded')
      call allocate_and_grid(layout);call precalculate_metric_cache();call prepare_case_boundary_data()
      call init_wavespeed_solver();call initialize_halo_exchange(layout);call set_initial_conditions()
      do k=1,nz
        do j=1,ny
          do i=1,nx
            call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
          end do
        end do
      end do
      call fill_primitive_halos(p)
      initial_u=up;initial_p=p;initial_events=numerical_events;initial_delta=numerical_delta_u
      initial_time=integration_time;initial_dt=dt;initial_step=n_steps;initial_flag=perturbation_applied
      if(trim(dest_label)=='settings-failure')rho_floor=2*rho_floor
      if(trim(mode)=='serial-probe')then
        ios=0
        if(parallel_rank==0)then
          open(newunit=file,file=trim(source_label),access='stream',form='unformatted',status='old',action='readwrite')
          inquire(unit=file,size=length)
          read(file,pos=length)byte
          write(file,pos=length)achar(ieor(iachar(byte),1))
          close(file,iostat=ios)
        end if
        call parallel_require_collective(ios==0,'serial mutation probe prepared')
        call import_serial_checkpoint(trim(source_label),layout,ios)
      else
        call load_parallel_checkpoint(trim(source_label),layout,ios,.true.)
      end if
      call require(ios/=0,'bad distributed source or settings rejected')
      call require(all(up==initial_u).and.all(p==initial_p),'failed redistribution preserves fields and halos')
      call require(all(numerical_events==initial_events).and.all(numerical_delta_u==initial_delta), &
        'failed redistribution preserves counters')
      call require(integration_time==initial_time.and.dt==initial_dt.and.n_steps==initial_step.and. &
        (perturbation_applied.eqv.initial_flag),'failed redistribution preserves epoch and flag')
      call teardown()
    end if
    if(parallel_rank==0)print *,'RESTART REJECTION PROBE PASSED'
    call parallel_finalize();stop
  end subroutine probe

  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    checks=checks+1
    if(.not.ok)call parallel_abort('C02 test: '//label)
  end subroutine require

  subroutine poison_configuration()
    nx=-7;ny=-5;nz=-6;nghost=-1;adb_idx=-4;g1=-2;bh_mass=-1;a_spin=-2
    rho_floor=-1;p_floor=-1;D_floor=-1;tau_floor=-1
    metric_type='invalid';geom_type='invalid';case_id=-1;case_name='invalid';scheme_name='invalid'
    output_folder='invalid';output_prefix='invalid';rec_method_id=-1;riemann_solver_id=-1
  end subroutine poison_configuration

  subroutine write_history(label,before,after)
    character(len=*),intent(in)::label
    type(state_diagnostics),intent(in)::before,after
    integer::file
    if(parallel_rank/=0)return
    open(newunit=file,file=trim(case_dir)//'/history_'//trim(label)//'.txt',status='new')
    write(file,'(*(I0,1X))')before%events
    write(file,'(*(ES24.16E3,1X))')before%raw_repair
    write(file,'(*(I0,1X))')after%events
    write(file,'(*(ES24.16E3,1X))')after%raw_repair
    close(file)
  end subroutine write_history

  subroutine legacy_reference_setup(path)
    character(len=*),intent(in)::path
    integer::file
    checkpoint_interval=100;extraction_stride=10;mass_monitor_stride=1000
    do_numerical_diagnostics=.false.;vtk_mapping_id=VTK_MAP_PHYSICAL
    walltime_limit_seconds=0;walltime_margin_seconds=900;walltime_checkpoint_seconds=0
    open(newunit=file,file=path,status='old',form='unformatted')
    call read_checkpoint_metadata(file)
    call choose_numerical_architecture();call validate_parameters(path)
    g1=adb_idx/(adb_idx-1);call set_metric_type()
    call build_domain([nx,ny,nz],nghost,0,1,case_id==4.or.case_id==6,layout,arg,[1,1], &
      (case_id==4.or.case_id==6).and.ny>1)
    call require(len_trim(arg)==0,'serial reference layout')
    call allocate_and_grid(layout);call precalculate_metric_cache();call init_wavespeed_solver()
    call read_checkpoint_state(file,n_steps,integration_time);close(file);dt=0
    call prepare_case_boundary_data();u=up
    old_events=numerical_events;old_delta=numerical_delta_u
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call cons_to_prim(u(i,j,k,:),p(i,j,k,:),i,j,k)
        end do
      end do
    end do
    numerical_events=old_events;numerical_delta_u=old_delta
    do_radial_mass_budget=.false.;do_restart=.false.
    call initialize_halo_exchange(layout);call fill_primitive_halos(p)
  end subroutine legacy_reference_setup

  subroutine setup()
    integer :: shape_global(3),parts(2),i,j,k
    logical :: periodic
    character(len=160) :: message
    shape_global=[37,19,8];parts=grid;periodic=.true.
    case_name='Conv';case_id=6;geom_type='Cartesian';metric_type='Minkowski';a_spin=0;bh_mass=0
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
      bh_mass=1
      geom_type='Spherical';metric_type='Eddington-Finkelstein';r_min=1.2d0;r_max=40;use_log_r=.true.
      michel_critical_radius=8;michel_critical_density=1
    case(6:13)
      bh_mass=1
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

end program test_parallel_restart
