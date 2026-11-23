! C03 CPU validation entry; the production grial MPI guard remains in place.
program mpi_segment_driver
  use variables
  use parallel_runtime
  use parallel_timing,only:any_time_request,clip_step_to_target
  use parallel_diagnostics,only:initialize_parallel_diagnostics,finalize_parallel_diagnostics
  use parallel_checkpoint,only:read_parallel_checkpoint_configuration
  use parallel_observer
  use parallel_visualization,only:save_parallel_vtk
  use parallel_perturbations,only:apply_parallel_perturbation
  use domain_decomposition
  use initialization,only:initialize_problem,allocate_and_grid,precalculate_metric_cache
  use parameters,only:read_parameter_header
  use conditions,only:set_initial_conditions,prepare_case_boundary_data
  use equations,only:prim_to_cons
  use fluxes,only:init_wavespeed_solver
  use halo_exchange
  use evolution,only:update_adaptive_dt
  use time_integrator,only:advance_rk_step
  use accelerator,only:wall_time_seconds
  use,intrinsic::iso_c_binding,only:c_char,c_int,c_null_char
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(domain_layout)::layout
  type(observer_calendar)::calendar
  character(len=512)::source,checkpoint,status_file,arg
  character(len=160)::message
  character(len=16)::mode
  character(len=24)::state
  integer::grid(2),ios,i,j,k,quota,session,start_step,last_checkpoint=-1,unit,close_status
  integer::last_vtk=-1
  character(len=512)::vtk_root,vtk_path
  real*8::started,budget,margin,step_started,max_step,checkpoint_cost,last_checkpoint_wall,measured,setup_done
  logical::periodic,injected,stop_requested,checkpoint_requested
  interface
    integer(c_int) function c_mkdir(path,permissions)bind(C,name='mkdir')
      import c_char,c_int
      character(kind=c_char),intent(in)::path(*)
      integer(c_int),value::permissions
    end function c_mkdir
  end interface

  started=wall_time_seconds();call parallel_initialize()
  call parallel_require_collective(command_argument_count()==8, &
    'expected cold/resume source grid_r grid_theta seconds margin test_step_quota session_id')
  call get_command_argument(1,mode);call get_command_argument(2,source)
  call get_command_argument(3,arg);read(arg,*,iostat=ios)grid(1)
  call parallel_require_collective(ios==0,'invalid radial process grid')
  call get_command_argument(4,arg);read(arg,*,iostat=ios)grid(2)
  call parallel_require_collective(ios==0,'invalid polar process grid')
  call get_command_argument(5,arg);read(arg,*,iostat=ios)budget
  call parallel_require_collective(ios==0,'invalid wall-clock budget')
  call get_command_argument(6,arg);read(arg,*,iostat=ios)margin
  call parallel_require_collective(ios==0,'invalid wall-clock margin')
  call get_command_argument(7,arg);read(arg,*,iostat=ios)quota
  call parallel_require_collective(ios==0.and.quota>=0,'invalid test step quota')
  call get_command_argument(8,arg);read(arg,*,iostat=ios)session
  call parallel_require_collective(ios==0.and.session>0,'invalid session number')
  call parallel_require_collective(all(ieee_is_finite([budget,margin])),'nonfinite wall-clock controls')
  call parallel_require_collective(budget>=0.0d0.and.margin>=0.0d0,'negative wall-clock controls')
  if(budget>0.0d0)call parallel_require_collective(margin>0.0d0.and.margin<budget,'require 0 < margin < budget')
  call parallel_require_collective(trim(mode)=='cold'.or.trim(mode)=='resume','invalid segment operation')
  if(trim(mode)=='cold')then
    call read_parameter_header(trim(source))
    call parallel_require_collective(.not.do_restart,'cold segment requires a cold parameter file')
    call initialize_problem(trim(source),.true.,distributed=.true.)
    periodic=case_id==4.or.case_id==6
    call build_domain([nx,ny,nz],nghost,parallel_rank,parallel_size,periodic,layout,message,grid,periodic.and.ny>1)
    call parallel_require_collective(len_trim(message)==0,message)
  else
    call read_parallel_checkpoint_configuration(trim(source),layout,ios,grid)
    call parallel_require_collective(ios==0,'cannot load segment configuration')
  end if
  call parallel_require_collective(do_numerical_diagnostics.and..not.do_radial_mass_budget.and. &
    .not.accelerator_enabled.and..not.is_mhd.and.diagnostic_stride>0.and.extraction_stride>0, &
    'C03 validation driver requires CPU hydro/D01 and positive diagnostic strides, without serial radial budgets')
  call allocate_and_grid(layout);call precalculate_metric_cache();call init_wavespeed_solver()
  if(trim(mode)=='cold')then
    call set_initial_conditions()
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
        end do
      end do
    end do
    u=up;integration_time=0.0d0;n_steps=0;dt=0.0d0
    numerical_events=0_int64;numerical_delta_u=0.0d0;perturbation_applied=.false.
  end if
  call prepare_case_boundary_data();call initialize_halo_exchange(layout)
  if(trim(mode)=='cold')then
    call fill_primitive_halos(p);call initialize_parallel_diagnostics(layout)
    call initialize_observer_calendar(calendar)
  else
    call load_observed_checkpoint(trim(source),layout,calendar,ios)
    call parallel_require_collective(ios==0,'observed segment recovery failed')
    last_checkpoint=n_steps;checkpoint=source
  end if
  call mkdirs(trim(output_folder)//'/'//trim(scheme_name));call mkdirs('checkpoints')
  vtk_root=trim(output_folder)//'/'//trim(scheme_name)//'/mpi_fields'
  call mkdirs(trim(vtk_root))
  call apply_parallel_perturbation(layout,injected)
  if(trim(mode)=='cold')call sample_parallel_observer(calendar,.true.)
  if(trim(mode)=='cold')then
    call vtk_now()
  else
    ! The recovered epoch was already published by the preceding C04 session.
    ! Older C03 checkpoints import fields but cannot recreate past visualizations.
    last_vtk=n_steps
  end if
  start_step=n_steps;max_step=0.0d0;checkpoint_cost=0.0d0
  last_checkpoint_wall=wall_time_seconds();setup_done=last_checkpoint_wall-started
  state='completed'
  do while(integration_time<final_time)
    ! Each rank checks its own elapsed clock; a request from any rank stops all.
    stop_requested=.false.
    if(budget>0.0d0)stop_requested=wall_time_seconds()-started+max_step>=budget-margin-checkpoint_cost
    if(quota>0.and.parallel_rank==parallel_size-1)stop_requested=stop_requested.or.n_steps-start_step>=quota
    if(any_time_request(stop_requested))then
      state='stopped_for_walltime';exit
    end if
    step_started=wall_time_seconds()
    call update_adaptive_dt();dt=clip_step_to_target(integration_time,dt,final_time)
    call parallel_require_collective(n_steps<huge(1),'step count exceeds integer range')
    n_steps=n_steps+1;integration_time=integration_time+dt
    call advance_rk_step(.true.)
    call apply_parallel_perturbation(layout,injected)
    call sample_parallel_observer(calendar)
    if(integration_time>=calendar%next_save)then
      call vtk_now()
      calendar%next_save=calendar%next_save+save_interval;calendar%save_index=calendar%save_index+1_int64
    end if
    checkpoint_requested=integration_time>=calendar%next_checkpoint
    if(walltime_checkpoint_seconds>0.0d0)checkpoint_requested=checkpoint_requested.or. &
      wall_time_seconds()-last_checkpoint_wall>=walltime_checkpoint_seconds
    ! Every rank must enter collective checkpoint IO in the same RK step.
    checkpoint_requested=any_time_request(checkpoint_requested)
    if(checkpoint_requested)then
      calendar%next_checkpoint=dble(int(integration_time/checkpoint_interval)+1)*checkpoint_interval
      call checkpoint_now()
    end if
    max_step=max(max_step,wall_time_seconds()-step_started)
  end do
  call sample_parallel_observer(calendar,.true.)
  call vtk_now()
  call checkpoint_now()
  write(status_file,'(A,I6.6,A)')'session_',session,'.dat'
  ios=0
  if(parallel_rank==0)then
    open(newunit=unit,file=trim(status_file)//'.tmp',status='new',action='write',iostat=ios)
    if(ios==0)then
      write(unit,'(A)',iostat=ios)'state='//trim(state)
      if(ios==0)write(unit,'(A)',iostat=ios)'checkpoint='//trim(checkpoint)
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'time=',integration_time
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'target=',final_time
      if(ios==0)write(unit,'(A,I0)',iostat=ios)'step=',n_steps
      if(ios==0)write(unit,'(A,I0)',iostat=ios)'ranks=',parallel_size
      if(ios==0)write(unit,'(A,I0)',iostat=ios)'grid_r=',grid(1)
      if(ios==0)write(unit,'(A,I0)',iostat=ios)'grid_theta=',grid(2)
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'elapsed=',wall_time_seconds()-started
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'setup_seconds=',setup_done
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'max_step_seconds=',max_step
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'checkpoint_seconds=',checkpoint_cost
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'wall_budget=',budget
      if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'margin=',margin
      close(unit,iostat=close_status)
      if(ios==0)ios=close_status
      if(ios==0)call rename(trim(status_file)//'.tmp',trim(status_file),ios)
    end if
  end if
  call parallel_require_collective(ios==0,'cannot publish segment receipt')
  call finalize_parallel_diagnostics();call finalize_halo_exchange();call parallel_finalize()
contains
  subroutine vtk_now()
    if(last_vtk==n_steps)return
    write(vtk_path,'(A,A,I10.10)')trim(vtk_root),'/step_',n_steps
    call save_parallel_vtk(trim(vtk_path),layout,ios)
    call parallel_require_collective(ios==0,'cannot publish distributed VTK frame')
    last_vtk=n_steps
  end subroutine vtk_now

  subroutine mkdirs(path)
    character(len=*),intent(in)::path
    integer::a,rc
    logical::exists
    rc=0
    if(parallel_rank==0)then
      do a=1,len_trim(path)
        if(path(a:a)/='/'.and.a/=len_trim(path))cycle
        if(a==1)cycle
        inquire(file=path(:a),exist=exists)
        if(.not.exists)then
          rc=c_mkdir(path(:a)//c_null_char,int(o'700',c_int))
          if(rc/=0)exit
        end if
      end do
    end if
    call parallel_require_collective(rc==0,'cannot create segment output directory')
  end subroutine mkdirs

  subroutine checkpoint_now()
    real*8::before
    if(last_checkpoint==n_steps)return
    call sample_parallel_observer(calendar,.true.)
    write(checkpoint,'(A,I9.9)')'checkpoints/step_',n_steps
    before=wall_time_seconds()
    call save_observed_checkpoint(trim(checkpoint),layout,calendar,ios)
    call parallel_require_collective(ios==0,'cannot publish observed checkpoint')
    measured=wall_time_seconds()-before
    checkpoint_cost=max(checkpoint_cost,measured)
    last_checkpoint=n_steps;last_checkpoint_wall=wall_time_seconds()
  end subroutine checkpoint_now
end program mpi_segment_driver
