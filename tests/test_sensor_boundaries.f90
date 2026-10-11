program test_sensor_boundaries
  use variables
  use initialization, only: choose_numerical_architecture
  use conditions, only: set_boundary_conditions
  use reconstruction, only: reconstruct_1d_core
  use fluxes, only: resolve_riemann_problem,init_wavespeed_solver
  implicit none
  integer :: method,sensor,solver,i
  real*8, allocatable :: state(:,:,:,:),line(:,:),left(:,:),right(:,:)
  real*8 :: g(3,3),beta(3),flux_left(neq),flux_right(neq)

  nx=32; ny=1; nz=1
  case_name='Conv'
  metric_type='Minkowski'; geom_type='Cartesian'
  rho_floor=1.0d-10; p_floor=1.0d-13
  adb_idx=4.0d0/3.0d0; g1=adb_idx/(adb_idx-1.0d0)
  tvd_limiter_id=LIM_MC
  g=0.0d0; beta=0.0d0
  do i=1,3
    g(i,i)=1.0d0
  end do
  do solver=RS_HLLE,RS_HLLC
    riemann_solver_id=solver
    call init_wavespeed_solver()
    do method=REC_GODUNOV,REC_WENO5
      rec_method_id=method
      do sensor=0,1
        use_shock_sensor=sensor == 1
        call choose_numerical_architecture()
        if ((use_shock_sensor .and. method /= REC_GODUNOV) .or. accelerator_enabled) then
          if (nghost < 3) error stop 'sensor/device needs three filled ghost layers'
        end if
        allocate(state(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq))
        allocate(line(neq,-nghost:nx+nghost),left(neq,0:nx),right(neq,0:nx))
        state=-123.0d0 ! Poison all ghosts: filling must not leave a sensor sample stale.
        do i=1,nx
          state(i,1,1,eq_de)=1.0d0+0.1d0*sin(2.0d0*pi*(i-0.5d0)/nx)
          state(i,1,1,eq_pr)=1.0d0
          state(i,1,1,eq_vx)=0.5d0
          state(i,1,1,eq_vy:eq_vz)=0.0d0
        end do
        call set_boundary_conditions(state)
        if (use_shock_sensor .and. method /= REC_GODUNOV) then
          if (any(state(-2,1,1,:) /= state(nx-2,1,1,:))) &
            error stop 'leftmost sensor ghost not periodic'
          if (any(state(nx+3,1,1,:) /= state(3,1,1,:))) &
            error stop 'rightmost sensor ghost not periodic'
        end if
        do i=-nghost,nx+nghost
          line(:,i)=state(i,1,1,:)
        end do
        call reconstruct_1d_core(line,nx,left,right,DIR_X)
        if (maxval(abs(left(:,0)-left(:,nx))) > 1.0d-14 .or. &
            maxval(abs(right(:,0)-right(:,nx))) > 1.0d-14) &
          error stop 'periodic reconstruction faces differ'
        call resolve_riemann_problem(left(:,0),right(:,0),DIR_X,1.0d0,beta,g,1.0d0,flux_left)
        call resolve_riemann_problem(left(:,nx),right(:,nx),DIR_X,1.0d0,beta,g,1.0d0,flux_right)
        if (maxval(abs(flux_left-flux_right)) > 1.0d-14) &
          error stop 'periodic fluxes do not telescope'
        deallocate(state,line,left,right)
      end do
    end do
  end do
  print *, 'SENSOR BOUNDARY TESTS PASSED'
end program test_sensor_boundaries
