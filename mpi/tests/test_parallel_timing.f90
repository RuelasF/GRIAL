program test_parallel_timing
  use parallel_runtime
  use parallel_timing
  use domain_decomposition
  use variables
  use evolution, only: update_adaptive_dt
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use metrics, only: set_metric_type
  use omp_lib, only: omp_get_thread_num
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  integer :: dimension,variant,geometry,checks=0,unit,timeline_steps=0
  character(len=48) :: mode,report_name
  real*8 :: cfl_error=0.0d0,wrong_min_2d=0.0d0,wrong_min_3d=0.0d0,temporary_step

  call parallel_initialize()
  call get_command_argument(1,mode)
  if (len_trim(mode)>0.and.trim(mode)/='report') call check_failure()
  do dimension=1,3
    do variant=0,1
      call check_directional_maxima(dimension,variant==1)
    end do
  end do
  do geometry=1,5
    call check_metric_partition(geometry)
  end do
  call check_requests()
  call check_periodic_boundary()
  call check_timeline(0.25d0,0.5d0)
  call check_timeline(0.3d0,0.6d0)
  if (trim(mode)=='report') then
    write(report_name,'(A,I0,A)') 'timing_rank_',parallel_rank,'.txt'
    open(newunit=unit,file=trim(report_name),status='new',action='write')
    write(unit,'(4(I0,1X),3(ES24.16E3,1X))') parallel_rank,parallel_size,checks,timeline_steps, &
      cfl_error,wrong_min_2d,wrong_min_3d
    close(unit)
  end if
  write(*,'(A,I0,A,I0,A,I0)') 'TIMING_TEST rank=',parallel_rank,' size=',parallel_size,' checks=',checks
  if (parallel_rank==0) write(*,'(A)') 'PARALLEL TIMING TESTS PASSED'
  call parallel_finalize()

contains

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if (.not.ok) call parallel_abort('timing test: '//label)
  end subroutine require

  subroutine synthetic_cache(dim,shifted)
    integer, intent(in) :: dim
    logical, intent(in) :: shifted
    type(domain_layout) :: block
    character(len=160) :: message
    integer :: global_shape(3),i,j,k,d,global_i
    global_shape=[37,1,1]
    if (dim>=2) global_shape(2)=3
    if (dim==3) global_shape(3)=2
    call build_domain(global_shape,3,parallel_rank,parallel_size,.false.,block,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=block%local_shape(1);ny=global_shape(2);nz=global_shape(3)
    nx_global=global_shape(1);ny_global=ny;nz_global=nz;x_global_offset=block%offset(1)
    dx=0.5d0;dy=0.25d0;dz=0.125d0;CFL=0.4d0
    if (dim==1) dy=0.0d0
    if (dim<=2) dz=0.0d0
    allocate(alpha_c(nx,ny),beta_c(3,nx,ny),gamma_inv_c(3,3,nx,ny))
    alpha_c=1.0d0;beta_c=0.0d0;gamma_inv_c=0.0d0
    if (shifted) then
      alpha_c=0.8d0
      beta_c(1,:,:)=0.1d0;beta_c(2,:,:)=-0.2d0;beta_c(3,:,:)=0.05d0
    end if
    do k=1,nz;do j=1,ny;do i=1,nx
      do d=1,3
        gamma_inv_c(d,d,i,j)=1.0d0
      end do
      global_i=x_global_offset+i
      if (global_i==1) gamma_inv_c(1,1,i,j)=100.0d0
      if (global_i==37) gamma_inv_c(2,2,i,j)=400.0d0
      if (global_i==19) gamma_inv_c(3,3,i,j)=900.0d0
    end do;end do;end do
    ! Synthetic geometry writes bypass precalculate_metric_cache. Publish their generation.
    metric_cache_generation=metric_cache_generation+1_int64
  end subroutine synthetic_cache

  subroutine check_directional_maxima(dim,shifted)
    integer, intent(in) :: dim
    logical, intent(in) :: shifted
    real*8 :: expected_maxima(3),local_maxima(3),rate,expected,wrong,wrong_global
    real*8 :: lapse,shift(3)
    integer :: i,j,k,d
    logical :: axes(3)
#ifdef USE_MPI
    integer :: status
#endif
    call synthetic_cache(dim,shifted)
    axes=[.true.,dim>=2,dim==3]
    lapse=1.0d0;shift=0.0d0
    if (shifted) then
      lapse=0.8d0;shift=[0.1d0,-0.2d0,0.05d0]
    end if
    expected_maxima=lapse*[10.0d0,20.0d0,30.0d0]+abs(shift)
    rate=expected_maxima(1)/dx
    if (dim>=2) rate=rate+expected_maxima(2)/dy
    if (dim==3) rate=rate+expected_maxima(3)/dz
    expected=CFL/rate
    local_maxima=0.0d0
    do k=1,nz;do j=1,ny;do i=1,nx
      do d=1,3
        if (axes(d)) local_maxima(d)=max(local_maxima(d),alpha_c(i,j)*sqrt(gamma_inv_c(d,d,i,j))+abs(beta_c(d,i,j)))
      end do
    end do;end do;end do
    wrong=cfl_step_from_maxima(local_maxima,[dx,dy,dz],axes,CFL)
    wrong_global=wrong
#ifdef USE_MPI
    call MPI_Allreduce(wrong,wrong_global,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,status)
    call parallel_check_mpi(status,'test minimum of local dt')
#endif
    call update_adaptive_dt()
    cfl_error=max(cfl_error,abs(dt-expected))
    call require(dt==expected,'directional global CFL equals serial formula')
    if (dim>=2.and.parallel_size>1) call require(wrong_global>dt,'minimum local dt misses separated directional maxima')
    if (.not.shifted.and.dim==2) wrong_min_2d=wrong_global
    if (.not.shifted.and.dim==3) wrong_min_3d=wrong_global
    deallocate(alpha_c,beta_c,gamma_inv_c)
  end subroutine check_directional_maxima

  subroutine check_metric_partition(which)
    integer, intent(in) :: which
    type(domain_layout) :: block
    character(len=160) :: message
    integer :: global_shape(3)
    real*8 :: reference_dt
    global_shape=[37,1,1]
    geom_type='Cartesian';metric_type='Minkowski';use_log_r=.false.
    r_min=0.0d0;r_max=1.0d0;y_min=0.0d0;y_max=pi;z_min=0.0d0;z_max=1.0d0
    bh_mass=1.0d0;a_spin=0.0d0;do_mdot_extraction=.false.;nghost=3;CFL=0.3d0
    select case(which)
    case(2)
      geom_type='Cylindrical';global_shape(3)=7
    case(3)
      geom_type='Spherical';metric_type='Eddington-Finkelstein';r_min=1.2d0;r_max=40.0d0
    case(4)
      geom_type='Spheroidal';metric_type='Kerr-Schild';r_min=1.2d0;r_max=40.0d0
      use_log_r=.true.;global_shape(2)=5
    case(5)
      geom_type='Spheroidal';metric_type='Kerr-Schild';r_min=1.2d0;r_max=40.0d0
      use_log_r=.true.;a_spin=0.9d0;global_shape=[37,4,6];z_max=2.0d0*pi
    end select
    call set_metric_type()
    nx=global_shape(1);ny=global_shape(2);nz=global_shape(3)
    call allocate_and_grid();call precalculate_metric_cache();call update_adaptive_dt()
    reference_dt=dt
    call release_grid()
    call build_domain(global_shape,3,parallel_rank,parallel_size,.false.,block,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=global_shape(1);ny=global_shape(2);nz=global_shape(3)
    call allocate_and_grid(block);call precalculate_metric_cache();call update_adaptive_dt()
    cfl_error=max(cfl_error,abs(dt-reference_dt))
    call require(dt==reference_dt,'partitioned metric CFL equals complete-grid CFL')
    call release_grid()
  end subroutine check_metric_partition

  subroutine check_periodic_boundary()
    real*8 :: interval,deadline,before_boundary
    type(time_decisions) :: decisions

    interval=0.1d0
    deadline=43.0d0*interval
    before_boundary=nearest(deadline,-1.0d0)
    call require(periodic_output_due(before_boundary,deadline),'one ULP below the deadline is due')
    call require(.not.periodic_output_due(deadline-1.0d-10,deadline),'an early sample is not due')
    call require(next_periodic_index(deadline-1.0d-10,interval)==43_int64,'retain a future slot')
    call require(next_periodic_index(before_boundary,interval)==44_int64,'advance past a rounded boundary')
    call require(next_periodic_index(deadline,interval)==44_int64,'advance past a rounded quotient')
    call require(next_periodic_index(0.0d0,interval)==1_int64,'first slot after initial output')
    call coordinate_time_actions(before_boundary,7,6.0d0,[0.0d0,deadline,deadline], &
      [.false.,.true.,.true.],[.false.,.false.],decisions)
    call require(decisions%save.and.decisions%checkpoint,'VTK and checkpoint agree at a rounded boundary')
  end subroutine check_periodic_boundary

  subroutine check_requests()
    type(time_decisions) :: decisions
    logical :: enabled(3),request
    real*8 :: events(3),nan
    enabled=[.true.,.true.,.true.];events=[0.5d0,0.25d0,0.75d0]
    call coordinate_time_actions(0.5d0,2,1.0d0,events,enabled,[parallel_rank==parallel_size-1,.false.],decisions)
    call require(decisions%stop.and.decisions%checkpoint,'one rank stop requests a complete checkpoint for all')
    call require(decisions%perturb.and.decisions%save.and..not.decisions%complete,'common physical events')
    call coordinate_time_actions(0.0d0,0,1.0d0,events,enabled,[.false.,parallel_rank==parallel_size-1],decisions)
    call require(decisions%checkpoint.and..not.decisions%stop,'one walltime checkpoint request propagates')
    call coordinate_time_actions(1.0d0,4,1.0d0,events,enabled,[parallel_rank==0,.false.],decisions)
    call require(decisions%complete.and..not.decisions%stop,'physical completion wins over walltime stop')
    nan=ieee_value(0.0d0,ieee_quiet_nan)
    call coordinate_time_actions(1.25d0,5,1.0d0,[nan,-1.0d0,nan],[.false.,.false.,.false.], &
      [.false.,.false.],decisions)
    call require(decisions%complete.and..not.decisions%perturb.and..not.decisions%save,'completed restart and disabled events')
    request=any_time_request(parallel_rank==parallel_size-1)
    call require(request,'walltime request sampled after output propagates')
    request=any_time_request(.false.)
    call require(.not.request,'no spurious walltime request')
    call require(clip_step_to_target(1.25d0,0.1d0,1.0d0)==0.0d0,'completed target does not move backward')
  end subroutine check_requests

  subroutine check_timeline(proposal,expected_injection_time)
    real*8, intent(in) :: proposal,expected_injection_time
    real*8 :: time,next_save,next_cp,chosen,old_step,injection_time
    integer :: step,perturbations,saves,checkpoints
    logical :: applied,request
    type(time_decisions) :: decisions
    time=0.0d0;step=0;next_save=0.125d0;next_cp=0.375d0
    applied=.false.;perturbations=0;saves=0;checkpoints=0;injection_time=-1.0d0
    do
      call coordinate_time_actions(time,step,1.0d0,[0.375d0,next_save,next_cp], &
        [.not.applied,.true.,.true.],[.false.,.false.],decisions)
      if (decisions%complete) exit
      call require(.not.decisions%stop,'no stop before physical completion')
      if (decisions%perturb) then
        applied=.true.;perturbations=perturbations+1;injection_time=time
      end if
      old_step=proposal
      if (time+old_step>1.0d0) old_step=1.0d0-time
      chosen=clip_step_to_target(time,proposal,1.0d0)
      call require(chosen==old_step,'final-time clipping preserves serial arithmetic')
      time=time+chosen;step=step+1
      call coordinate_time_actions(time,step,1.0d0,[0.375d0,next_save,next_cp], &
        [.not.applied,.true.,.true.],[.false.,.false.],decisions)
      call require(decisions%save.eqv.(time>=next_save),'save threshold unchanged')
      call require(decisions%checkpoint.eqv.(time>=next_cp),'checkpoint threshold unchanged')
      if (decisions%save) then
        saves=saves+1;next_save=next_save+0.125d0
      end if
      request=any_time_request(step==1.and.parallel_rank==parallel_size-1)
      if (decisions%checkpoint.or.request) then
        checkpoints=checkpoints+1
        next_cp=dble(int(time/0.375d0)+1)*0.375d0
      end if
    end do
    call require(time==1.0d0.and.step==4,'common exact endpoint')
    call require(perturbations==1.and.injection_time==expected_injection_time,'single start-of-step delayed injection')
    call require(saves==4.and.checkpoints==3,'serial event cadence and global walltime request')
    timeline_steps=timeline_steps+step
  end subroutine check_timeline

  subroutine check_failure()
    type(time_decisions) :: decisions
    real*8 :: nan,time,target,events(3),proposal
    integer :: step
    logical :: enabled(3),last
    nan=ieee_value(0.0d0,ieee_quiet_nan)
    last=parallel_rank==parallel_size-1
    if (index(trim(mode),'cfl-')==1) then
      call synthetic_cache(3,.false.)
      if (last) then
        select case(trim(mode))
        case('cfl-spacing');dx=0.0d0
        case('cfl-nan-spacing');dx=nan
        case('cfl-courant');CFL=0.3d0
        case('cfl-nan-courant');CFL=nan
        case('cfl-axis');ny_global=1
        case('cfl-nan-metric');alpha_c(1,1)=nan
        case('cfl-lapse');alpha_c(1,1)=0.0d0
        case('cfl-inverse');gamma_inv_c(1,1,1,1)=-1.0d0
        case('cfl-shift');beta_c(1,1,1)=ieee_value(0.0d0,ieee_positive_inf)
        case('cfl-cache');deallocate(alpha_c)
        case('cfl-shape');deallocate(beta_c);allocate(beta_c(2,nx,ny));beta_c=0.0d0
        end select
      end if
      if (trim(mode)=='cfl-thread') then
        !$omp parallel num_threads(2)
        if (omp_get_thread_num()==1) call update_adaptive_dt()
        !$omp end parallel
      else
        call update_adaptive_dt()
      end if
    else
      time=0.25d0;target=1.0d0;step=1;events=[0.5d0,0.5d0,0.75d0];enabled=.true.;proposal=0.25d0
      if (last) then
        select case(trim(mode))
        case('time-epoch');time=0.5d0
        case('time-step');step=2
        case('time-target');target=1.25d0
        case('time-calendar');events(2)=0.625d0
        case('time-enabled');enabled(1)=.false.
        case('time-nan');time=nan
        case('time-calendar-nan');events(2)=nan
        case('clip-stagnant');time=1.0d0;target=2.0d0;proposal=1.0d-300
        case('clip-nan');proposal=nan
        end select
      end if
      if (index(trim(mode),'clip-')==1) then
        temporary_step=clip_step_to_target(time,proposal,target)
      else if (trim(mode)=='rate-overflow') then
        if (last) temporary_step=cfl_step_from_maxima([0.5d0*huge(1.0d0),0.0d0,0.0d0], &
          [tiny(1.0d0),1.0d0,1.0d0],[.true.,.false.,.false.],0.4d0)
      else
        call coordinate_time_actions(time,step,target,events,enabled,[.false.,.false.],decisions)
      end if
    end if
    call parallel_require_collective(.true.,'waiting for invalid peer')
    call parallel_abort('invalid timing case was not rejected')
  end subroutine check_failure

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y)) deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z)) deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
  end subroutine release_grid
end program test_parallel_timing
