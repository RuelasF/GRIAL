! Optional MPI lifecycle. No domain decomposition or hydrodynamic communication yet.
module parallel_runtime
  use, intrinsic :: iso_fortran_env, only: error_unit, int64
#ifdef _OPENMP
  use omp_lib, only: omp_in_parallel, omp_get_wtime
#endif
#ifdef USE_MPI
  use mpi_f08
  use, intrinsic :: ieee_exceptions, only: ieee_status_type,ieee_get_status,ieee_set_status,ieee_set_halting_mode,ieee_all
#endif
  implicit none
  private
  public :: parallel_initialize, parallel_finalize, parallel_require_collective, &
            parallel_abort, parallel_check_mpi, parallel_rank, parallel_size, &
            parallel_thread_level, parallel_enabled, parallel_require_active
#ifdef USE_MPI
  logical, parameter :: parallel_enabled = .true.
#else
  logical, parameter :: parallel_enabled = .false.
#endif
#ifdef GRIAL_PROFILE
  public :: profile_reset, profile_add, profile_report, profile_clock, profile_label_add
  real*8 :: profile_seconds(9)=0.0d0
  integer :: profile_calls(9)=0
  character(len=160) :: profile_labels(256)=""
  real*8 :: profile_label_seconds(256)=0.0d0
  integer :: profile_label_calls(256)=0
  ! Records are buffered outside the timed MPI call. All ranks on one host
  ! share the monotonic clock used by GNU Fortran's int64 SYSTEM_CLOCK.
  integer, parameter :: trace_capacity=131072
  integer, allocatable :: trace_site(:)
  real*8, allocatable :: trace_enter(:),trace_leave(:)
  integer :: trace_count=0
  logical :: trace_truncated=.false.
#endif
  integer, protected :: parallel_rank=0, parallel_size=1, parallel_thread_level=0
  logical :: active=.false.

contains

#ifdef GRIAL_PROFILE
  real*8 function profile_clock()
    integer(int64) :: ticks,rate
    call system_clock(ticks,rate)
    profile_clock=dble(ticks)/dble(rate)
  end function profile_clock

  subroutine profile_reset()
    character(len=8) :: setting
    integer :: status

    profile_seconds=0.0d0; profile_calls=0
    profile_labels=""; profile_label_seconds=0.0d0; profile_label_calls=0
    trace_count=0;trace_truncated=.false.
    call get_environment_variable('GRIAL_COLLECTIVE_TRACE',setting,status=status)
    if (allocated(trace_site)) deallocate(trace_site,trace_enter,trace_leave)
    if (status==0.and.setting=='1') allocate(trace_site(trace_capacity), &
      trace_enter(trace_capacity),trace_leave(trace_capacity))
  end subroutine profile_reset

  subroutine profile_add(slot,start)
    integer, intent(in) :: slot
    real*8, intent(in) :: start
    !$omp atomic update
    profile_seconds(slot)=profile_seconds(slot)+profile_clock()-start
    !$omp atomic update
    profile_calls(slot)=profile_calls(slot)+1
  end subroutine profile_add

  subroutine profile_label_add(label,start)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: start
    integer :: i
    real*8 :: finished

    finished=profile_clock()
    do i=1,size(profile_labels)
      if (profile_labels(i)==label.or.profile_labels(i)=="") then
        profile_labels(i)=label
        profile_label_seconds(i)=profile_label_seconds(i)+(finished-start)
        if (allocated(trace_site)) then
          if (trace_count<size(trace_site)) then
            trace_count=trace_count+1
            trace_site(trace_count)=i
            trace_enter(trace_count)=start
            trace_leave(trace_count)=finished
          else
            trace_truncated=.true.
          end if
        end if
        profile_label_calls(i)=profile_label_calls(i)+1
        return
      end if
    end do
    error stop "too many profiling labels"
  end subroutine profile_label_add

  subroutine profile_report(root,start_step)
    character(len=*), intent(in) :: root
    integer, intent(in) :: start_step
    integer :: i,unit,ios,close_status
    character(len=80) :: trace_name
    character(len=:), allocatable :: trace_path
    character(len=16),parameter :: names(9)=[character(len=16) :: 'cfl','rhs','recovery','halos', &
      'observer','validate_mpi','reconstruct','riemann_loop','source_assembly']
    do i=1,9
      write(*,'(A,I0,1X,A,1X,I0,1X,ES24.16E3)') 'GRIAL_PROFILE rank=',parallel_rank, &
        trim(names(i)),profile_calls(i),profile_seconds(i)
    end do
    do i=1,size(profile_labels)
      if (profile_label_calls(i)==0) cycle
      write(*,'(A,I0,1X,A,1X,I0,1X,ES24.16E3)') "GRIAL_PROFILE_COLLECTIVE rank=",parallel_rank, &
        trim(profile_labels(i)),profile_label_calls(i),profile_label_seconds(i)
    end do
    if (allocated(trace_site)) then
      write(trace_name,'(A,I0,A,I0,A)') 'grial_collectives_step_',start_step,'_rank_',parallel_rank,'.tsv'
      trace_path=trim(root)//'/'//trim(trace_name)
      open(newunit=unit,file=trace_path,status='new',action='write',iostat=ios)
      if (ios==0) then
        write(unit,'(A,L1)',iostat=ios) '# truncated=',trace_truncated
        do i=1,trace_count
          if (ios/=0) exit
          write(unit,'(I0,A,ES24.16E3,A,ES24.16E3,2A)',iostat=ios) i,achar(9),trace_enter(i),achar(9), &
            trace_leave(i),achar(9),trim(profile_labels(trace_site(i)))
        end do
        close(unit,iostat=close_status)
        if (ios==0) ios=close_status
      end if
      if (ios/=0) write(error_unit,'(A,I0,2A)') 'Collective profile unavailable on rank ',parallel_rank,': ',trace_path
    end if
  end subroutine profile_report
#endif

  subroutine require_main_thread()
    ! FUNNELED calls remain outside OpenMP regions, including the primary thread.
#ifdef _OPENMP
    if (omp_in_parallel()) then
      write(error_unit,'(A)') 'PARALLEL ERROR: MPI runtime called inside an OpenMP region'
      error stop 1
    end if
#endif
  end subroutine require_main_thread

  subroutine parallel_require_active()
    call require_main_thread()
    if (.not.active) call parallel_abort('parallel runtime is not active')
  end subroutine parallel_require_active

  subroutine parallel_initialize()
#ifdef USE_MPI
    integer :: status
    logical :: was_initialized
    type(ieee_status_type) :: fp_status
#endif
    call require_main_thread()
    if (active) call parallel_abort('parallel runtime initialized twice')
#ifdef USE_MPI
    was_initialized=.false.
    call MPI_Initialized(was_initialized,status)
    if (status/=MPI_SUCCESS.or.was_initialized) error stop 'MPI must be initialized once by GRIAL'
    ! Provider discovery can raise FP exceptions; preserve the solver's FP environment.
    call ieee_get_status(fp_status)
    call ieee_set_halting_mode(ieee_all,.false.)
    call MPI_Init_thread(MPI_THREAD_FUNNELED,parallel_thread_level,status)
    call ieee_set_status(fp_status)
    if (status/=MPI_SUCCESS) error stop 'MPI_Init_thread failed'
    active=.true.
    call MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,status)
    call parallel_check_mpi(status,'MPI_Comm_set_errhandler')
    call MPI_Comm_rank(MPI_COMM_WORLD,parallel_rank,status)
    call parallel_check_mpi(status,'MPI_Comm_rank')
    call MPI_Comm_size(MPI_COMM_WORLD,parallel_size,status)
    call parallel_check_mpi(status,'MPI_Comm_size')
    if (parallel_thread_level<MPI_THREAD_FUNNELED) call parallel_abort('MPI_THREAD_FUNNELED not provided')
#else
    call reject_nonmpi_launcher()
    active=.true.
#endif
  end subroutine parallel_initialize

#ifndef USE_MPI
  subroutine reject_nonmpi_launcher()
    character(len=20), parameter :: names(2)=[character(len=20)::'OMPI_COMM_WORLD_SIZE','PMI_SIZE']
    character(len=32) :: value
    integer :: i,status,count
    do i=1,size(names)
      call get_environment_variable(trim(names(i)),value,status=status)
      if (status==1) cycle
      if (status/=0) call parallel_abort('cannot read MPI launcher environment')
      read(value,*,iostat=status) count
      if (status/=0) call parallel_abort('invalid MPI launcher size')
      if (count>1) call parallel_abort('This executable has no MPI; rebuild with USE_MPI=1 and use grial_mpi')
    end do
  end subroutine reject_nonmpi_launcher
#endif

  subroutine parallel_require_collective(local_ok,message)
    logical, intent(in) :: local_ok
    character(len=*), intent(in) :: message
    logical :: all_ok
#ifdef GRIAL_PROFILE
    real*8 :: profile_start
#endif
#ifdef USE_MPI
    integer :: status
#endif
    call require_main_thread()
    if (.not.active) call parallel_abort('parallel runtime is not active')
#ifdef USE_MPI
#ifdef GRIAL_PROFILE
    profile_start=profile_clock()
#endif
    call MPI_Allreduce(local_ok,all_ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
    call profile_label_add(message,profile_start)
    call profile_add(6,profile_start)
#endif
    call parallel_check_mpi(status,'MPI_Allreduce validation')
#else
    all_ok=local_ok
#endif
    if (.not.all_ok) then
      if (.not.local_ok) write(error_unit,'(A,I0,2A)') 'PARALLEL ERROR rank=',parallel_rank,': ',trim(message)
      call parallel_abort('collective validation failed: '//trim(message))
    end if
  end subroutine parallel_require_collective

  subroutine parallel_check_mpi(status,operation)
    integer, intent(in) :: status
    character(len=*), intent(in) :: operation
#ifdef USE_MPI
    integer :: length,error_status
    character(len=MPI_MAX_ERROR_STRING) :: description
#endif
    call require_main_thread()
#ifdef USE_MPI
    if (status==MPI_SUCCESS) return
    if (.not.active) call parallel_abort('MPI error without active runtime')
    call MPI_Error_string(status,description,length,error_status)
    if (error_status==MPI_SUCCESS) then
      call parallel_abort(trim(operation)//': '//description(:length))
    else
      call parallel_abort(trim(operation)//': MPI error')
    end if
#else
    if (status/=0) call parallel_abort(trim(operation)//': nonzero status')
#endif
  end subroutine parallel_check_mpi

  subroutine parallel_abort(message,code)
    character(len=*), intent(in) :: message
    integer, intent(in), optional :: code
#ifdef USE_MPI
    integer :: status,exit_code
#endif
    call require_main_thread()
    write(error_unit,'(A,I0,2A)') 'PARALLEL ERROR rank=',parallel_rank,': ',trim(message)
#ifdef USE_MPI
    if (active) then
      exit_code=1
      if (present(code)) exit_code=max(1,code)
      call MPI_Abort(MPI_COMM_WORLD,exit_code,status)
    end if
#endif
    error stop 1
  end subroutine parallel_abort

  subroutine parallel_finalize()
#ifdef USE_MPI
    integer :: status
    type(ieee_status_type) :: fp_status
#endif
    call require_main_thread()
    if (.not.active) return
#ifdef USE_MPI
    call ieee_get_status(fp_status)
    call ieee_set_halting_mode(ieee_all,.false.)
    call MPI_Finalize(status)
    call ieee_set_status(fp_status)
    if (status/=MPI_SUCCESS) call parallel_abort('MPI_Finalize failed')
#endif
    active=.false.
  end subroutine parallel_finalize
end module parallel_runtime
