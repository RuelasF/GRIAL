! Validation entry point only; block dumps are analysis products, not restart files.
program srhd_long_driver
  use variables
  use parallel_runtime
  use parallel_timing, only: clip_step_to_target,any_time_request
  use parallel_diagnostics, only: initialize_parallel_diagnostics,finalize_parallel_diagnostics
  use domain_decomposition
  use parameters, only: read_parameter_header
  use initialization, only: initialize_problem,allocate_and_grid,precalculate_metric_cache
  use conditions, only: set_initial_conditions,set_boundary_conditions,prepare_case_boundary_data
  use equations, only: prim_to_cons
  use fluxes, only: init_wavespeed_solver
  use halo_exchange
  use evolution, only: update_adaptive_dt
  use time_integrator, only: advance_rk_step
  use output, only: write_distributed_diagnostics
  use accelerator, only: wall_time_seconds
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(domain_layout) :: layout
  character(len=1024) :: par,seed,arg
  character(len=160) :: message
  character(len=128) :: filename
  character(len=16), parameter :: magic='GRIAL_LONG_V1'
  character(len=16) :: marker
  integer :: shape_global(3),head(6),unit_file,ios,close_status,i,j,k,frame,trace_unit,last_sample=-1
  integer(int64) :: seed_bytes
  real*8 :: start_wall,budget,next_output,seed_clock(2)
  real*8, allocatable :: seed_u(:,:,:,:),seed_p(:,:,:,:)
  logical :: periodic,ok,stopped

  call parallel_initialize()
  call parallel_require_collective(command_argument_count()==3,'expected parameter file, seed file/native, wall budget')
  call get_command_argument(1,par);call get_command_argument(2,seed);call get_command_argument(3,arg)
  read(arg,*,iostat=ios)budget
  call parallel_require_collective(ios==0,'invalid validation wall budget')
  call parallel_require_collective(ieee_is_finite(budget).and.budget>0.0d0,'wall budget must be positive')
  start_wall=wall_time_seconds()
  ! The controller validates these frozen parameters with grial --check before launching MPI.
  call read_parameter_header(trim(par))
  call parallel_require_collective(.not.do_restart.and.(case_id==4.or.case_id==5),'driver supports cold KHI/Jet only')
  call initialize_problem(trim(par),.true.)
  call parallel_require_collective(ny==1.and..not.accelerator_enabled.and..not.apply_perturbation.and. &
    do_numerical_diagnostics.and..not.do_radial_mass_budget,'unsupported long-driver configuration')
  shape_global=[nx,ny,nz];periodic=case_id==4
  call build_domain(shape_global,nghost,parallel_rank,parallel_size,periodic,layout,message)
  call parallel_require_collective(len_trim(message)==0,message)
  call allocate_and_grid(layout);call precalculate_metric_cache();call init_wavespeed_solver()
  if (trim(seed)=='native') then
    call parallel_require_collective(.not.parallel_enabled.and.parallel_size==1,'native seed requires non-MPI reference')
    call set_initial_conditions();call set_boundary_conditions(p)
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  else
    allocate(seed_u(shape_global(1),1,nz,neq),seed_p(shape_global(1),1,nz,neq))
    open(newunit=unit_file,file=trim(seed),access='stream',form='unformatted',status='old',iostat=ios)
    ok=ios==0
    if (ok) then
      read(unit_file,iostat=ios)marker,head,seed_clock
      ok=ios==0
      if (ok)ok=marker==magic.and.all(head==[shape_global(1),nz,0,shape_global(1),0,neq]).and. &
        all(seed_clock==0.0d0)
      inquire(unit=unit_file,size=seed_bytes)
      ok=ok.and.seed_bytes==56_int64+16_int64*shape_global(1)*nz*neq
      if (ok) then
        read(unit_file,iostat=ios)seed_u,seed_p
        ok=ios==0
      end if
      close(unit_file,iostat=close_status);ok=ok.and.close_status==0
    end if
    if (ok)ok=all(ieee_is_finite(seed_u)).and.all(ieee_is_finite(seed_p))
    call parallel_require_collective(ok,'missing, invalid or incompatible canonical seed')
    up=seed_u(x_global_offset+1:x_global_offset+nx,:,:,:)
    p(1:nx,1:ny,1:nz,:)=seed_p(x_global_offset+1:x_global_offset+nx,:,:,:)
    deallocate(seed_u,seed_p)
    call prepare_case_boundary_data()
  end if
  integration_time=0.0d0;n_steps=0;dt=0.0d0
  numerical_events=0_int64;numerical_delta_u=0.0d0
  call initialize_halo_exchange(layout)
  if (parallel_enabled) then
    call fill_primitive_halos(p)
  else
    call set_boundary_conditions(p)
  end if
  call initialize_parallel_diagnostics(layout)
  ok=.true.
  if (parallel_rank==0) then
    open(newunit=trace_unit,file='trace.dat',status='new',action='write',iostat=ios)
    ok=ios==0
    if (ok)write(trace_unit,'(A)',iostat=ios)'# step time dt elapsed_wall_seconds'
    ok=ok.and.ios==0
  end if
  call parallel_require_collective(ok,'cannot create trace.dat')
  frame=0;call sample();next_output=save_interval;stopped=.false.
  do while(integration_time<final_time)
    if (any_time_request(wall_time_seconds()-start_wall>=budget)) then
      stopped=.true.;exit
    end if
    call update_adaptive_dt()
    dt=clip_step_to_target(integration_time,dt,final_time)
    n_steps=n_steps+1;integration_time=integration_time+dt
    call advance_rk_step(parallel_enabled)
    ok=.true.
    if (parallel_rank==0) then
      write(trace_unit,'(I0,3(1X,ES24.16E3))',iostat=ios)n_steps,integration_time,dt,wall_time_seconds()-start_wall
      ok=ios==0
      if (.not.ok)call parallel_abort('cannot write validation time trace')
    end if
    ! Reuse the common event cadence; do not alter dt to hit visualization times.
    if (integration_time>=next_output.or.integration_time>=final_time) then
      call parallel_require_collective(ok,'cannot write time trace')
      frame=frame+1;call sample();next_output=next_output+save_interval
      if (parallel_rank==0) then
        flush(trace_unit,iostat=ios)
        if (ios/=0)call parallel_abort('cannot flush validation time trace')
      end if
    end if
  end do
  if (stopped.and.n_steps>last_sample) then
    frame=frame+1;call sample()
  end if
  ok=.true.
  if (parallel_rank==0) then
    close(trace_unit,iostat=ios);ok=ios==0
    open(newunit=unit_file,file='driver_status.dat',status='new',action='write',iostat=ios)
    ok=ok.and.ios==0
    if (ios==0) then
      if (stopped) then
        write(unit_file,'(A)',iostat=ios)'status=incomplete_walltime'
      else
        write(unit_file,'(A)',iostat=ios)'status=completed'
      end if
      if (ios==0)write(unit_file,'(A,ES24.16E3)',iostat=ios)'physical_time=',integration_time
      if (ios==0)write(unit_file,'(A,I0)',iostat=ios)'steps=',n_steps
      if (ios==0)write(unit_file,'(A,I0)',iostat=ios)'ranks=',parallel_size
      if (ios==0)write(unit_file,'(A,ES24.16E3)',iostat=ios)'elapsed_seconds=',wall_time_seconds()-start_wall
      ok=ok.and.ios==0
      close(unit_file,iostat=close_status);ok=ok.and.close_status==0
    end if
  end if
  call parallel_require_collective(ok,'cannot publish driver status')
  call finalize_parallel_diagnostics();call finalize_halo_exchange();call parallel_finalize()

contains

  subroutine sample()
    logical :: valid
    valid=all(ieee_is_finite(up)).and.all(ieee_is_finite(p(1:nx,1:ny,1:nz,:)))
    if (valid)valid=minval(p(1:nx,1:ny,1:nz,eq_de))>0.0d0.and. &
      minval(p(1:nx,1:ny,1:nz,eq_pr))>0.0d0
    call parallel_require_collective(valid,'nonfinite or nonpositive long-run state')
    write(filename,'(A,I6.6,A,I4.4,A)')'frame_',frame,'_rank_',parallel_rank,'.bin'
    open(newunit=unit_file,file=trim(filename),access='stream',form='unformatted', &
      status='new',action='write',iostat=ios)
    valid=ios==0
    if (valid) then
      write(unit_file,iostat=ios)magic,[nx_global,nz_global,x_global_offset,nx,n_steps,neq], &
        [integration_time,dt],up,p(1:nx,1:ny,1:nz,:)
      valid=ios==0
      close(unit_file,iostat=close_status);valid=valid.and.close_status==0
    end if
    call parallel_require_collective(valid,'cannot write validation block')
    call write_distributed_diagnostics(integration_time)
    last_sample=n_steps
    if (parallel_rank==0) then
      write(*,'(A,I0,A,I0,A,ES12.4,A,ES12.4)')'LONG_PROGRESS frame=',frame, &
        ' step=',n_steps,' time=',integration_time,' elapsed=',wall_time_seconds()-start_wall
    end if
  end subroutine sample
end program srhd_long_driver
