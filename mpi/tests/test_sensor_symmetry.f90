program test_sensor_symmetry
  use variables
  use reconstruction, only: reconstruct_1d_core,reconstruct_face_state
  implicit none
  integer, parameter :: nodes=40
  real*8, parameter :: tolerance=5.0d-14
  real*8 :: line(neq,-3:nodes+3),mirror(neq,-3:nodes+3),signs(neq)
  real*8 :: left(neq,0:nodes),right(neq,0:nodes),mirror_left(neq,0:nodes),mirror_right(neq,0:nodes)
  real*8 :: stencil(neq,-2:3),face_left(neq),face_right(neq),position,error,worst_error
  integer(int64) :: counts(7)
  integer :: cell,wrapped,face,offset,method,sensor,profile,direction,checks

  nghost=3
  metric_type='Minkowski'; curved_metric=.false.
  rho_floor=1.0d-10; p_floor=1.0d-13
  tvd_limiter_id=LIM_MC
  do_numerical_diagnostics=.not. accelerator_enabled
  worst_error=0.0d0; checks=0

  ! Reflection swaps the two face states and reverses the normal velocity.
  ! Cover KHI, moderate/strong jumps, and the atmosphere fallback in all axes.
  do profile=1,4
    do cell=1,nodes
      position=-0.5d0+(dble(cell)-0.5d0)/dble(nodes)
      line(:,cell)=0.0d0
      line(eq_de,cell)=1.0d0
      line(eq_pr,cell)=2.5d0
      line(eq_vx,cell)=-0.5d0
      if (abs(position) < 0.25d0) then
        line(eq_de,cell)=2.0d0
        line(eq_vx,cell)=0.5d0
      end if
      line(eq_vz,cell)=0.01d0*exp(-0.5d0*((abs(position)-0.25d0)/0.05d0)**2)
      if (profile >= 2) then
        line(eq_de,cell)=1.0d0+0.1d0*sin(7.0d0*position)
        line(eq_pr,cell)=0.4d0+0.03d0*cos(9.0d0*position)
        line(eq_vx,cell)=0.05d0*sin(3.0d0*position)
        line(eq_vy,cell)=0.02d0*cos(5.0d0*position)
        if (cell > nodes/2) then
          select case(profile)
          case(2)
            line(eq_de,cell)=1.65d0*line(eq_de,cell)
            line(eq_pr,cell)=1.70d0*line(eq_pr,cell)
          case(3)
            line(eq_de,cell)=100.0d0*line(eq_de,cell)
            line(eq_pr,cell)=100.0d0*line(eq_pr,cell)
          case(4)
            line(eq_de,cell)=10.0d0*rho_floor
            line(eq_pr,cell)=10.0d0*p_floor
          end select
        end if
      end if
    end do
    do cell=-3,nodes+3
      wrapped=modulo(cell-1,nodes)+1
      if (cell < 1 .or. cell > nodes) line(:,cell)=line(:,wrapped)
    end do
    do direction=DIR_X,DIR_Z
      signs=1.0d0
      signs(eq_vx+direction-1)=-1.0d0
      do cell=-3,nodes+3
        wrapped=modulo(cell-1,nodes)+1
        mirror(:,cell)=signs*line(:,nodes+1-wrapped)
      end do
      do sensor=0,2
        use_shock_sensor=sensor /= 0
        shock_sensor_thresholds=SENSOR_DEFAULTS
        if (sensor == 2) shock_sensor_thresholds=[0.80d0,0.40d0,1.50d0,0.90d0,10.0d0,10.0d0,100.0d0]
        do method=REC_GODUNOV,REC_WENO5
          rec_method_id=method
          numerical_events=0_int64
          call reconstruct_1d_core(line,nodes,left,right,direction)
          counts=numerical_events(1:7)
          numerical_events=0_int64
          call reconstruct_1d_core(mirror,nodes,mirror_left,mirror_right,direction)
          if (do_numerical_diagnostics) then
            if (sum(counts(1:5)) /= nodes+1) error stop 'one reconstruction count per face'
            if (any(counts /= numerical_events(1:7))) error stop 'reflection changes reconstruction counters'
          end if
          error=0.0d0
          do face=0,nodes
            error=max(error,maxval(abs(left(:,face)-signs*mirror_right(:,nodes-face))), &
                            maxval(abs(right(:,face)-signs*mirror_left(:,nodes-face))))
            do offset=-2,3
              stencil(:,offset)=line(:,face+offset)
            end do
            call reconstruct_face_state(stencil,face,nodes,direction,method,LIM_MC,use_shock_sensor, &
                 .false.,rho_floor,p_floor,face_left,face_right,shock_sensor_thresholds)
            error=max(error,maxval(abs(left(:,face)-face_left)),maxval(abs(right(:,face)-face_right)))
            do offset=-2,3
              stencil(:,offset)=mirror(:,nodes-face+offset)
            end do
            call reconstruct_face_state(stencil,nodes-face,nodes,direction,method,LIM_MC,use_shock_sensor, &
                 .false.,rho_floor,p_floor,face_left,face_right,shock_sensor_thresholds)
            error=max(error,maxval(abs(left(:,face)-signs*face_right)), &
                            maxval(abs(right(:,face)-signs*face_left)))
          end do
          if (error > tolerance) then
            write(*,'(A,4(1X,I0),1X,ES14.6)') 'FAIL reflection profile/direction/sensor/method:', &
                 profile,direction,sensor,method,error
            error stop 'sensor breaks reflection symmetry'
          end if
          worst_error=max(worst_error,error)
          checks=checks+1
        end do
      end do
    end do
  end do
  write(*,'(A,I0,A,ES12.4)') 'SENSOR SYMMETRY TESTS PASSED: ',checks,' checks; max error ',worst_error
end program test_sensor_symmetry
