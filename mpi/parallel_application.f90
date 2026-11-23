! Explicit CPU entry joining the tested distributed solver, observer and IO APIs.
module parallel_application
  use variables
  use parallel_runtime
  use parallel_timing, only: any_time_request,clip_step_to_target,periodic_output_due,next_periodic_index
  use parallel_diagnostics, only: initialize_parallel_diagnostics,finalize_parallel_diagnostics
  use parallel_checkpoint, only: read_parallel_checkpoint_configuration,checkpoint_file_signature
  use parallel_observer
  use parallel_visualization, only: save_parallel_vtk
  use parallel_perturbations, only: apply_parallel_perturbation
  use domain_decomposition
  use initialization, only: initialize_problem,allocate_and_grid,precalculate_metric_cache
  use parameters, only: read_parameter_header
  use conditions, only: set_initial_conditions,prepare_case_boundary_data
  use equations, only: prim_to_cons
  use fluxes, only: init_wavespeed_solver
  use halo_exchange
  use evolution, only: update_adaptive_dt
  use time_integrator, only: advance_rk_step
  use accelerator, only: wall_time_seconds
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: agree_command_line,run_distributed_application
  interface
    integer(c_int) function c_mkdir(path,permissions) bind(C,name='mkdir')
      import c_char,c_int
      character(kind=c_char),intent(in) :: path(*)
      integer(c_int),value :: permissions
    end function c_mkdir
    integer(c_int) function c_rmdir(path) bind(C,name='rmdir')
      import c_char,c_int
      character(kind=c_char),intent(in) :: path(*)
    end function c_rmdir
  end interface
contains
  subroutine agree_command_line()
    integer :: count,root_count,i,length,ios
    character(len=1024) :: value,root_value
    count=command_argument_count();root_count=count
#ifdef USE_MPI
    call MPI_Bcast(root_count,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
    call parallel_check_mpi(ios,'application argument count')
#endif
    call parallel_require_collective(count==root_count.and.count<=24,'command lines differ or have too many arguments')
    do i=1,count
      call get_command_argument(i,value,length=length,status=ios)
      call parallel_require_collective(ios==0.and.length>0.and.length<=len(value).and.index(value,achar(0))==0, &
        'empty, truncated or invalid command argument')
      root_value=value
#ifdef USE_MPI
      call MPI_Bcast(root_value,len(root_value),MPI_CHARACTER,0,MPI_COMM_WORLD,ios)
      call parallel_check_mpi(ios,'application argument agreement')
#endif
      call parallel_require_collective(value==root_value,'command arguments differ between ranks')
    end do
  end subroutine agree_command_line

  subroutine agree_input(path)
    character(len=*),intent(in) :: path
    integer(int64) :: signature(2),root_signature(2)
    integer :: ios
#ifdef USE_MPI
    type(MPI_Datatype) :: count_type
#endif
    call checkpoint_file_signature(path,signature,ios)
    call parallel_require_collective(ios==0,'cannot read distributed parameter file: '//trim(path))
    root_signature=signature
#ifdef USE_MPI
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(0_int64)/8,count_type,ios)
    call parallel_check_mpi(ios,'application signature datatype')
    call MPI_Bcast(root_signature,2,count_type,0,MPI_COMM_WORLD,ios)
    call parallel_check_mpi(ios,'application input agreement')
#endif
    call parallel_require_collective(all(signature==root_signature),'parameter contents differ between ranks')
  end subroutine agree_input

  subroutine mkdirs(path)
    character(len=*),intent(in) :: path
    integer :: a,rc
    logical :: exists
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
    call parallel_require_collective(rc==0,'cannot create distributed output directory')
  end subroutine mkdirs

  subroutine run_distributed_application(start_wall)
    real*8,intent(in) :: start_wall
    type(domain_layout) :: layout
    type(observer_calendar) :: calendar
    integer(int64) :: save_slot
    logical :: physical_checkpoint_due
    character(len=512) :: source,arg,root,checkpoint_root,checkpoint,vtk_root,path,lock_path,status_path
    character(len=160) :: message
    character(len=24) :: reason
    integer :: grid(2),ios,i,j,k,n,unit,close_status,last_checkpoint,last_vtk,start_step
    integer(int64) :: expected_size,actual_size
    real*8 :: started,budget,margin,max_step,checkpoint_cost,last_checkpoint_wall,step_started,before
    real*8 :: requested_budget,requested_margin,setup_seconds
#ifdef GRIAL_PROFILE
    real*8 :: loop_started
#endif
    logical :: resume,check_only,injected,stop_requested,checkpoint_requested,exists
    logical :: grid_seen,budget_seen,margin_seen

    started=start_wall;resume=.false.;check_only=.false.;source=''
    grid=[parallel_size,1];grid_seen=.false.;budget_seen=.false.;margin_seen=.false.
    requested_budget=0.0d0;requested_margin=0.0d0;i=2
    do while(i<=command_argument_count())
      call get_command_argument(i,arg)
      select case(trim(arg))
      case('--check')
        call parallel_require_collective(.not.check_only,'duplicate --check')
        check_only=.true.
      case('--resume')
        call parallel_require_collective(.not.resume,'duplicate --resume')
        resume=.true.
      case('--grid')
        call parallel_require_collective(.not.grid_seen.and.i+2<=command_argument_count(),'expected --grid radial polar')
        call get_command_argument(i+1,arg);read(arg,*,iostat=ios)grid(1)
        call parallel_require_collective(ios==0,'invalid radial process count')
        call get_command_argument(i+2,arg);read(arg,*,iostat=ios)grid(2)
        call parallel_require_collective(ios==0,'invalid polar process count')
        grid_seen=.true.;i=i+2
      case('--walltime','--margin')
        checkpoint_requested=trim(arg)=='--walltime'
        call parallel_require_collective(i<command_argument_count(),'missing wall-clock option value')
        call get_command_argument(i+1,arg)
        if(checkpoint_requested)then
          call parallel_require_collective(.not.budget_seen,'duplicate --walltime')
          read(arg,*,iostat=ios)requested_budget;budget_seen=.true.
        else
          call parallel_require_collective(.not.margin_seen,'duplicate --margin')
          read(arg,*,iostat=ios)requested_margin;margin_seen=.true.
        end if
        call parallel_require_collective(ios==0,'invalid wall-clock value');i=i+1
      case default
        call parallel_require_collective(len_trim(source)==0.and.index(arg,'--')/=1,'unexpected distributed argument')
        call parallel_require_collective(len_trim(arg)<=440,'distributed source path exceeds 440 characters')
        source=arg
      end select
      i=i+1
    end do
    call parallel_require_collective(len_trim(source)>0,'distributed entry requires a parameter file or --resume generation')
    if(resume)then
      call read_parallel_checkpoint_configuration(trim(source),layout,ios,grid)
      call parallel_require_collective(ios==0,'cannot load distributed restart configuration')
    else
      call agree_input(trim(source));call read_parameter_header(trim(source))
      call parallel_require_collective(.not.do_restart,'use --distributed --resume with an observed checkpoint directory')
      call initialize_problem(trim(source),.true.,distributed=.true.)
      call build_domain([nx,ny,nz],nghost,parallel_rank,parallel_size,case_id==4.or.case_id==6, &
        layout,message,grid,(case_id==4.or.case_id==6).and.ny>1)
      call parallel_require_collective(len_trim(message)==0,message)
    end if
    call parallel_require_collective(.not.accelerator_enabled.and..not.is_mhd.and..not.do_radial_mass_budget, &
      'distributed entry requires CPU hydro without serial radial budgets')
    call parallel_require_collective(do_numerical_diagnostics.and.diagnostic_stride>0.and.extraction_stride>0, &
      'distributed observer requires numerical_diagnostics=true and positive diagnostic/extraction strides')
    budget=walltime_limit_seconds;margin=walltime_margin_seconds
    if(budget_seen)budget=requested_budget
    if(margin_seen)margin=requested_margin
    call parallel_require_collective(all(ieee_is_finite([budget,margin])),'nonfinite wall-clock controls')
    call parallel_require_collective(budget>=0.0d0.and.margin>=0.0d0,'negative wall-clock controls')
    if(budget>0.0d0)call parallel_require_collective(margin>0.0d0.and.margin<budget,'require 0 < margin < walltime')
    root=trim(output_folder)//'/'//trim(scheme_name)
    call parallel_require_collective(len_trim(output_folder)+len_trim(scheme_name)+70<=440, &
      'distributed output path is too long')
    if(check_only.and..not.resume)then
      if(parallel_rank==0)write(*,'(A)')'DISTRIBUTED_CHECK passed; no output created'
      return
    end if
    call allocate_and_grid(layout);call precalculate_metric_cache();call init_wavespeed_solver()
    if(.not.resume)then
      call set_initial_conditions()
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      u=up;integration_time=0.0d0;n_steps=0;dt=0.0d0
      numerical_events=0_int64;numerical_delta_u=0.0d0;perturbation_applied=.false.
    end if
    call prepare_case_boundary_data();call initialize_halo_exchange(layout)
    if(resume)then
      call load_observed_checkpoint(trim(source),layout,calendar,ios)
      call parallel_require_collective(ios==0,'distributed observed checkpoint recovery failed')
    else
      call fill_primitive_halos(p);call initialize_parallel_diagnostics(layout)
      call initialize_observer_calendar(calendar)
    end if
    if(check_only)then
      if(parallel_rank==0)write(*,'(A)')'DISTRIBUTED_CHECK checkpoint passed; no output created'
      call finalize_parallel_diagnostics();call finalize_halo_exchange();return
    end if
    call mkdirs(trim(root))
    lock_path=trim(root)//'/.mpi_run_lock';status_path=trim(root)//'/run_status.dat'
    ios=0
    if(parallel_rank==0)then
      inquire(file=trim(status_path),exist=exists)
      if(exists.and..not.resume)then
        ios=1
      else
        ios=c_mkdir(trim(lock_path)//c_null_char,int(o'700',c_int))
      end if
    end if
    call parallel_require_collective(ios==0,'output already used or locked; inspect it and choose a fresh output_folder')
    checkpoint_root=trim(root)//'/mpi_checkpoints';vtk_root=trim(root)//'/mpi_fields'
    call mkdirs(trim(checkpoint_root));call mkdirs(trim(vtk_root))
    last_checkpoint=-1;last_vtk=-1;checkpoint='none'
    if(resume)then
      last_checkpoint=n_steps;last_vtk=n_steps;checkpoint=source
    end if
    ! The real deadlines remain on disk for compatibility with existing observers.
    save_slot=max(1_int64,nint(calendar%next_save/save_interval,kind=int64))
    if(resume.and.periodic_output_due(integration_time,calendar%next_checkpoint)) &
      calendar%next_checkpoint=dble(next_periodic_index(integration_time,checkpoint_interval))*checkpoint_interval
    call apply_parallel_perturbation(layout,injected)
    if(.not.resume)then
      call sample_parallel_observer(calendar,.true.);call vtk_now()
    end if
#ifdef GRIAL_PROFILE
    call profile_reset()
#endif
    start_step=n_steps;max_step=0.0d0;checkpoint_cost=0.0d0
    last_checkpoint_wall=wall_time_seconds();setup_seconds=last_checkpoint_wall-started
    reason='running';call publish_status()
#ifdef GRIAL_PROFILE
    loop_started=profile_clock()
#endif
    do while(integration_time<final_time)
      stop_requested=.false.
      if(budget>0.0d0)stop_requested=wall_time_seconds()-started+max_step>=budget-margin-checkpoint_cost
      if(any_time_request(stop_requested))exit
      step_started=wall_time_seconds()
#ifdef GRIAL_PROFILE
      block
        real*8 :: profile_start
        profile_start=profile_clock()
        call update_adaptive_dt()
        call profile_add(1,profile_start)
      end block
      dt=clip_step_to_target(integration_time,dt,final_time)
#else
      call update_adaptive_dt();dt=clip_step_to_target(integration_time,dt,final_time)
#endif
      call parallel_require_collective(n_steps<huge(1),'step count exceeds integer range')
      n_steps=n_steps+1;integration_time=integration_time+dt
      call advance_rk_step(.true.);call apply_parallel_perturbation(layout,injected)
#ifdef GRIAL_PROFILE
      block
        real*8 :: profile_start
        profile_start=profile_clock()
        call sample_parallel_observer(calendar)
        call profile_add(5,profile_start)
      end block
#else
      call sample_parallel_observer(calendar)
#endif
      if(periodic_output_due(integration_time,calendar%next_save))then
        call vtk_now()
        save_slot=save_slot+1_int64
        calendar%next_save=dble(save_slot)*save_interval
        calendar%save_index=calendar%save_index+1_int64
      end if
      physical_checkpoint_due=periodic_output_due(integration_time,calendar%next_checkpoint)
      checkpoint_requested=physical_checkpoint_due
      if(walltime_checkpoint_seconds>0.0d0)checkpoint_requested=checkpoint_requested.or. &
        wall_time_seconds()-last_checkpoint_wall>=walltime_checkpoint_seconds
      ! Wall clocks differ per rank: agree before entering collective checkpoint IO.
      checkpoint_requested=any_time_request(checkpoint_requested)
      if(checkpoint_requested)then
        if(physical_checkpoint_due) &
          calendar%next_checkpoint=dble(next_periodic_index(integration_time,checkpoint_interval))*checkpoint_interval
        call checkpoint_now();call publish_status()
      end if
      max_step=max(max_step,wall_time_seconds()-step_started)
    end do
#ifdef GRIAL_PROFILE
    write(*,'(A,I0,1X,ES24.16E3)') 'GRIAL_LOOP rank=',parallel_rank,profile_clock()-loop_started
    call profile_report(root,start_step)
#endif
    call sample_parallel_observer(calendar,.true.);call vtk_now();call checkpoint_now()
    reason='completed'
    if(integration_time<final_time)reason='stopped_for_walltime'
    call publish_status()
    ios=0
    if(parallel_rank==0)ios=c_rmdir(trim(lock_path)//c_null_char)
    call parallel_require_collective(ios==0,'cannot release distributed output lock')
    if(parallel_rank==0)write(*,'(A,A,A,ES24.16E3,A,I0)') &
      'DISTRIBUTED_RUN ',trim(reason),' time=',integration_time,' step=',n_steps
    call finalize_parallel_diagnostics();call finalize_halo_exchange()
  contains
    subroutine vtk_now()
      if(last_vtk==n_steps)return
      write(path,'(A,A,I10.10)')trim(vtk_root),'/step_',n_steps
      call save_parallel_vtk(trim(path),layout,ios)
      call parallel_require_collective(ios==0,'cannot publish distributed VTK frame')
      last_vtk=n_steps
    end subroutine vtk_now

    subroutine checkpoint_now()
      if(last_checkpoint==n_steps)return
      call sample_parallel_observer(calendar,.true.)
      write(checkpoint,'(A,A,I9.9)')trim(checkpoint_root),'/step_',n_steps
      before=wall_time_seconds()
      call save_observed_checkpoint(trim(checkpoint),layout,calendar,ios)
      call parallel_require_collective(ios==0,'cannot publish observed checkpoint')
      checkpoint_cost=max(checkpoint_cost,wall_time_seconds()-before)
      last_checkpoint=n_steps;last_checkpoint_wall=wall_time_seconds()
    end subroutine checkpoint_now

    subroutine publish_status()
      ios=0
      if(parallel_rank==0)then
        open(newunit=unit,file=trim(status_path)//'.tmp',status='new',access='stream', &
          form='formatted',action='write',iostat=ios)
        if(ios==0)then
          write(unit,'(A)',iostat=ios)'status='//trim(reason)
          if(ios==0)write(unit,'(A)',iostat=ios)'checkpoint='//trim(checkpoint)
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'physical_time=',integration_time
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'target_time=',final_time
          if(ios==0)write(unit,'(A,I0)',iostat=ios)'step=',n_steps
          if(ios==0)write(unit,'(A,I0)',iostat=ios)'ranks=',parallel_size
          if(ios==0)write(unit,'(A,I0)',iostat=ios)'grid_r=',grid(1)
          if(ios==0)write(unit,'(A,I0)',iostat=ios)'grid_theta=',grid(2)
          if(ios==0.and.polar_mesh_h/=1.0d0) &
            write(unit,'(A,ES24.16E3)',iostat=ios)'polar_mesh_h=',polar_mesh_h
          if(ios==0)write(unit,'(A,I0)',iostat=ios)'session_start_step=',start_step
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'elapsed_seconds=',wall_time_seconds()-started
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'setup_seconds=',setup_seconds
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'max_step_seconds=',max_step
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'checkpoint_seconds=',checkpoint_cost
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'wall_budget_seconds=',budget
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)'margin_seconds=',margin
          inquire(unit=unit,pos=expected_size)
          close(unit,iostat=close_status)
          if(ios==0)ios=close_status
          if(ios==0)then
            inquire(file=trim(status_path)//'.tmp',size=actual_size,iostat=ios)
            if(ios==0.and.actual_size/=expected_size-1_int64)ios=1
          end if
          if(ios==0)call rename(trim(status_path)//'.tmp',trim(status_path),ios)
        end if
      end if
      call parallel_require_collective(ios==0,'cannot publish distributed run status')
    end subroutine publish_status
  end subroutine run_distributed_application
end module parallel_application
