! CFL and temporal agreement; no fluid evolution, clocks or file writers.
module parallel_timing
  use parallel_runtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: validate_cfl_context,reduce_causal_maxima,cfl_step_from_maxima, &
            time_decisions,coordinate_time_actions,clip_step_to_target,any_time_request,validate_step_size

  type :: time_decisions
    logical :: perturb=.false.,save=.false.,checkpoint=.false.,stop=.false.,complete=.false.
  end type time_decisions

contains

  subroutine agree_metadata(values,valid,label)
    real*8, intent(in) :: values(:)
    logical, intent(in) :: valid
    character(len=*), intent(in) :: label
    real*8 :: signature(size(values)+1),minimum(size(values)+1),maximum(size(values)+1)
#ifdef USE_MPI
    integer :: status
#endif
    ! Never pass NaN to MPI_MIN/MAX, including a malformed parameter on one rank.
    signature=0.0d0
    if (valid) then
      signature(1:size(values))=values
    else
      signature(size(signature))=1.0d0
    end if
    minimum=signature
    maximum=signature
    if (parallel_size>1) then
      call parallel_require_active()
#ifdef USE_MPI
      call MPI_Allreduce(signature,minimum,size(signature),MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,status)
      call parallel_check_mpi(status,trim(label)//' minimum')
      call MPI_Allreduce(signature,maximum,size(signature),MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
      call parallel_check_mpi(status,trim(label)//' maximum')
#endif
    end if
    if (maximum(size(signature))>0.0d0) call parallel_abort('invalid '//trim(label))
    if (any(minimum/=maximum)) call parallel_abort('inconsistent '//trim(label)//' between ranks')
  end subroutine agree_metadata

  subroutine validate_cfl_context(spacing,active_axes,courant,cache_ok)
    real*8, intent(in) :: spacing(3),courant
    logical, intent(in) :: active_axes(3),cache_ok
    real*8 :: signature(7)
    logical :: valid
    integer :: direction
    signature=0.0d0
    valid=cache_ok.and.active_axes(1).and.ieee_is_finite(courant)
    if (ieee_is_finite(courant)) then
      valid=valid.and.courant>0.0d0
      signature(4)=courant
    end if
    do direction=1,3
      signature(4+direction)=merge(1.0d0,0.0d0,active_axes(direction))
      if (.not.active_axes(direction)) cycle
      if (.not.ieee_is_finite(spacing(direction))) then
        valid=.false.
      else
        valid=valid.and.spacing(direction)>0.0d0
        signature(direction)=spacing(direction)
      end if
    end do
    call agree_metadata(signature,valid,'CFL context or metric cache')
  end subroutine validate_cfl_context

  subroutine reduce_causal_maxima(local_maxima,local_ok,global_maxima)
    real*8, intent(in) :: local_maxima(3)
    logical, intent(in) :: local_ok
    real*8, intent(out) :: global_maxima(3)
    real*8 :: local(4),global(4)
    logical :: valid
#ifdef USE_MPI
    integer :: status
#endif
    valid=local_ok.and.all(ieee_is_finite(local_maxima))
    if (all(ieee_is_finite(local_maxima))) valid=valid.and.all(local_maxima>=0.0d0)
    local=0.0d0
    if (valid) then
      local(1:3)=local_maxima
    else
      local(4)=1.0d0
    end if
    global=local
    if (parallel_size>1) then
      call parallel_require_active()
#ifdef USE_MPI
      call MPI_Allreduce(local,global,4,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
      call parallel_check_mpi(status,'global directional causal maxima')
#endif
    end if
    if (global(4)>0.0d0) call parallel_abort('invalid metric data in CFL cells')
    global_maxima=global(1:3)
  end subroutine reduce_causal_maxima

  real*8 function cfl_step_from_maxima(maxima,spacing,active_axes,courant) result(step)
    real*8, intent(in) :: maxima(3),spacing(3),courant
    logical, intent(in) :: active_axes(3)
    real*8 :: inverse_dt,term
    integer :: direction
    if (.not.all(ieee_is_finite(maxima)).or..not.ieee_is_finite(courant)) &
      call parallel_abort('invalid CFL propagation data')
    if (any(maxima<0.0d0).or.courant<=0.0d0.or..not.active_axes(1)) &
      call parallel_abort('invalid CFL propagation data')
    inverse_dt=0.0d0
    do direction=1,3
      if (.not.active_axes(direction)) cycle
      if (.not.ieee_is_finite(spacing(direction))) call parallel_abort('invalid CFL spacing')
      if (spacing(direction)<=0.0d0) call parallel_abort('invalid CFL spacing')
      if (spacing(direction)<1.0d0) then
        if (maxima(direction)>huge(1.0d0)*spacing(direction)) call parallel_abort('CFL propagation rate overflow')
      end if
      term=maxima(direction)/spacing(direction)
      if (inverse_dt>huge(1.0d0)-term) call parallel_abort('CFL propagation sum overflow')
      inverse_dt=inverse_dt+term
    end do
    if (inverse_dt<=0.0d0) call parallel_abort('invalid multidimensional CFL propagation rate')
    if (inverse_dt<1.0d0) then
      if (courant>huge(1.0d0)*inverse_dt) call parallel_abort('CFL time step overflow')
    end if
    step=courant/inverse_dt
    if (.not.ieee_is_finite(step)) call parallel_abort('invalid CFL time step')
    if (step<=0.0d0) call parallel_abort('CFL time step does not advance')
  end function cfl_step_from_maxima

  subroutine coordinate_time_actions(time,step,target,event_times,event_enabled,local_requests,decisions)
    real*8, intent(in) :: time,target,event_times(3)
    integer, intent(in) :: step
    logical, intent(in) :: event_enabled(3),local_requests(2)
    type(time_decisions), intent(out) :: decisions
    real*8 :: signature(9),times(3)
    logical :: valid,due(3),requests(2)
#ifdef USE_MPI
    integer :: status
#endif
    call parallel_require_active()
    times=event_times
    where (.not.event_enabled) times=0.0d0
    valid=all(ieee_is_finite([time,target,times]))
    if (valid) valid=time>=0.0d0.and.target>=0.0d0.and.step>=0.and.all(times>=0.0d0)
    signature=[time,target,times,dble(step),merge(1.0d0,0.0d0,event_enabled)]
    call agree_metadata(signature,valid,'time epoch or event calendar')
    requests=local_requests
#ifdef USE_MPI
    if (parallel_size>1) then
      call MPI_Allreduce(local_requests,requests,2,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,status)
      call parallel_check_mpi(status,'checkpoint and stop requests')
    end if
#endif
    due=event_enabled.and.time>=times
    decisions%complete=time>=target
    decisions%perturb=due(1)
    decisions%save=due(2)
    decisions%stop=requests(1).and..not.decisions%complete
    decisions%checkpoint=due(3).or.requests(2).or.decisions%stop
  end subroutine coordinate_time_actions

  logical function any_time_request(local_request) result(request)
    logical, intent(in) :: local_request
#ifdef USE_MPI
    integer :: status
#endif
    call parallel_require_active()
    request=local_request
#ifdef USE_MPI
    if (parallel_size>1) then
      call MPI_Allreduce(local_request,request,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,status)
      call parallel_check_mpi(status,'walltime checkpoint request')
    end if
#endif
  end function any_time_request

  subroutine validate_step_size(step)
    real*8, intent(in) :: step
    logical :: valid
    call parallel_require_active()
    valid=ieee_is_finite(step)
    if (valid) valid=step>0.0d0
    call agree_metadata([step],valid,'RK time step')
  end subroutine validate_step_size

  real*8 function clip_step_to_target(time,proposal,target) result(step)
    real*8, intent(in) :: time,proposal,target
    if (.not.all(ieee_is_finite([time,proposal,target]))) call parallel_abort('invalid temporal step data')
    if (time<0.0d0.or.target<0.0d0.or.proposal<=0.0d0) call parallel_abort('invalid temporal step data')
    if (time>=target) then
      step=0.0d0
      return
    end if
    step=proposal
    if (proposal>huge(1.0d0)-time) then
      step=target-time
    else if (time+proposal>target) then
      step=target-time
    end if
    if (step<=0.0d0) call parallel_abort('time step does not advance physical time')
    if (time+step<=time) call parallel_abort('time step does not advance physical time')
  end function clip_step_to_target
end module parallel_timing
