program test_numerical_diagnostics
  use variables
  use equations
  use reconstruction
  use fluxes
  implicit none
  integer, parameter :: nodes=8
  integer :: i,method
  real*8 :: line(neq,-3:nodes+3),left(neq,0:nodes),right(neq,0:nodes)
  real*8 :: conservative(neq),before(neq),primitive(neq),g(3,3),beta(3),flux(neq)

  if (accelerator_enabled) then
    print *, 'NUMERICAL DIAGNOSTICS TESTS SKIPPED: CPU counters only'
    stop
  end if

  nghost=3
  nx=nodes
  ny=1
  nz=1
  rho_floor=1.0d-10
  p_floor=1.0d-13
  adb_idx=1.4d0
  g1=adb_idx/(adb_idx-1.0d0)
  metric_type='Minkowski'
  geom_type='Cartesian'
  tvd_limiter_id=LIM_MC
  use_shock_sensor=.true.
  do_numerical_diagnostics=.true.
  numerical_events=0_int64
  numerical_delta_u=0.0d0
  line=0.0d0
  line(eq_de,:)=1.0d0
  line(eq_pr,:)=1.0d0
  do method=REC_GODUNOV,REC_WENO5
    rec_method_id=method
    call reconstruct_1d_core(line,nodes,left,right,DIR_X)
    if (numerical_events(method) /= nodes+1) error stop 'uniform reconstruction count'
  end do

  numerical_events=0_int64
  rec_method_id=REC_WENO5
  line(eq_de,4:)=0.01d0
  call reconstruct_1d_core(line,nodes,left,right,DIR_X)
  if (sum(numerical_events(1:5)) /= nodes+1) error stop 'one count per face'
  ! Six cells centered on i+1/2 straddle this jump at faces 1,...,5.
  if (numerical_events(NUM_RECON_GODUNOV) /= 5) error stop 'shock stencil count'

  numerical_events=0_int64
  use_shock_sensor=.false.
  rec_method_id=REC_GODUNOV
  line(eq_de,:)=-1.0d0
  call reconstruct_1d_core(line,nodes,left,right,DIR_X)
  if (numerical_events(NUM_FACE_RHO_FLOOR) /= 2*(nodes+1)) error stop 'face floor count'
  if (any(left(eq_de,:) /= rho_floor)) error stop 'floor behavior'

  allocate(gamma_c(3,3,1,1,1),gamma_inv_c(3,3,1,1,1),sqrt_gamma_c(1,1,1))
  g=0.0d0
  do i=1,3
    g(i,i)=1.0d0
  end do
  gamma_c(:,:,1,1,1)=g
  gamma_inv_c(:,:,1,1,1)=g
  sqrt_gamma_c=1.0d0
  conservative=0.0d0
  conservative(eq_de)=0.5d0*rho_floor
  before=conservative
  call cons_to_prim(conservative,primitive,1,1,1)
  if (numerical_events(NUM_C2P_CALLS) /= 1 .or. numerical_events(NUM_C2P_LOW_D) /= 1) &
    error stop 'C2P atmosphere count'
  if (maxval(abs(numerical_delta_u-(conservative-before))) > 1.0d-25) error stop 'repair delta'

  beta=0.0d0
  primitive=(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0/)
  riemann_solver_id=RS_HLLC
  call init_wavespeed_solver()
  call resolve_riemann_problem(primitive,primitive,DIR_X,1.0d0,beta,g,1.0d0,flux)
  if (numerical_events(NUM_HLLC_CALLS) /= 1 .or. numerical_events(NUM_HLLC_FALLBACK) /= 0) &
    error stop 'successful HLLC count'
  ! alpha=0 deliberately triggers the guarded fallback; it is not a physical test.
  call resolve_riemann_problem(primitive,primitive,DIR_X,0.0d0,beta,g,1.0d0,flux)
  if (numerical_events(NUM_HLLC_CALLS) /= 2 .or. numerical_events(NUM_HLLC_FALLBACK) /= 1) &
    error stop 'guarded HLLC fallback count'

  numerical_events=0_int64
  !$omp parallel do
  do i=1,10000
    call count_numerical_event(NUM_HLLC_CALLS)
  end do
  !$omp end parallel do
  if (numerical_events(NUM_HLLC_CALLS) /= 10000) error stop 'OpenMP atomic lost events'
  do_numerical_diagnostics=.false.
  call count_numerical_event(NUM_HLLC_CALLS)
  if (numerical_events(NUM_HLLC_CALLS) /= 10000) error stop 'disabled counters incremented'
  print *, 'NUMERICAL DIAGNOSTICS TESTS PASSED'
end program test_numerical_diagnostics
