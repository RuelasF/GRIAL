program test_parallel_runtime
  use parallel_runtime
  use omp_lib, only: omp_get_thread_num
  use, intrinsic :: ieee_exceptions, only: ieee_all,ieee_get_flag,ieee_get_halting_mode
#ifdef USE_MPI
  use mpi_f08, only: MPI_COMM_WORLD,MPI_Barrier,MPI_Send,MPI_INTEGER,MPI_SUCCESS,MPI_THREAD_FUNNELED
#endif
  implicit none
  character(len=32) :: mode
  integer :: workers
  logical :: worker_modes_ok,worker_traps(size(ieee_all))
  logical :: traps_before(size(ieee_all)),traps_after(size(ieee_all)),flags_before(size(ieee_all)),flags_after(size(ieee_all))
#ifdef USE_MPI
  integer :: status,value
#endif
  call ieee_get_flag(ieee_all,flags_before)
  call ieee_get_halting_mode(ieee_all,traps_before)
  call parallel_initialize()
  call ieee_get_flag(ieee_all,flags_after)
  call ieee_get_halting_mode(ieee_all,traps_after)
  call parallel_require_collective(all(traps_before.eqv.traps_after).and.all(flags_before.eqv.flags_after), &
    'MPI initialization changed the floating point environment')
  call get_command_argument(1,mode)
  if (len_trim(mode)==0) mode='ok'
#ifdef USE_MPI
  call parallel_require_collective(parallel_enabled.and.parallel_thread_level>=MPI_THREAD_FUNNELED, &
    'FUNNELED support not available')
#else
  call parallel_require_collective(.not.parallel_enabled.and.parallel_size==1.and.parallel_rank==0, &
    'serial topology invalid')
#endif
  select case(trim(mode))
  case('ok')
    workers=0
    worker_modes_ok=.true.
    !$omp parallel private(worker_traps) reduction(+:workers) reduction(.and.:worker_modes_ok)
    call ieee_get_halting_mode(ieee_all,worker_traps)
    worker_modes_ok=worker_modes_ok.and.all(worker_traps.eqv.traps_before)
    workers=workers+1
    !$omp end parallel
    call parallel_require_collective(worker_modes_ok,'OpenMP worker FP traps changed')
    call parallel_require_collective(workers>=1,'OpenMP workers missing')
    write(*,'(A,I0,A,I0,A,I0)') 'PARALLEL_TEST rank=',parallel_rank,' size=',parallel_size,' workers=',workers
  case('collective-fail')
    call parallel_require_collective(parallel_rank/=parallel_size-1,'deliberate collective failure')
    error stop 'Collective failure was not detected'
  case('abort')
    if (parallel_rank==parallel_size-1) call parallel_abort('deliberate local failure',7)
#ifdef USE_MPI
    call MPI_Barrier(MPI_COMM_WORLD,status)
    call parallel_check_mpi(status,'peer barrier')
#endif
    error stop 'Abort did not stop peers'
  case('mpi-error')
#ifdef USE_MPI
    if (parallel_rank==parallel_size-1) then
      value=1
      call MPI_Send(value,1,MPI_INTEGER,parallel_size,99,MPI_COMM_WORLD,status)
      if (status==MPI_SUCCESS) call parallel_abort('invalid rank unexpectedly accepted')
      call parallel_check_mpi(status,'deliberate invalid destination')
    else
      call MPI_Barrier(MPI_COMM_WORLD,status)
      call parallel_check_mpi(status,'peer barrier')
    end if
#else
    call parallel_check_mpi(1,'deliberate serial status failure')
#endif
    error stop 'MPI error was not detected'
  case('thread-fail')
    !$omp parallel num_threads(2)
    if (omp_get_thread_num()==1) call parallel_require_collective(.true.,'worker thread must not enter MPI')
    !$omp end parallel
    error stop 'Thread guard did not trigger'
  case default
    call parallel_abort('unknown runtime test mode')
  end select
  call ieee_get_flag(ieee_all,flags_before)
  call ieee_get_halting_mode(ieee_all,traps_before)
  call parallel_finalize()
  call ieee_get_flag(ieee_all,flags_after)
  call ieee_get_halting_mode(ieee_all,traps_after)
  if (.not.all(traps_before.eqv.traps_after).or..not.all(flags_before.eqv.flags_after)) &
    error stop 'MPI finalization changed the floating point environment'
  if (parallel_rank==0) print *, 'PARALLEL RUNTIME TESTS PASSED'
end program test_parallel_runtime
