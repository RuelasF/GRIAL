program test_srhd_conservation
  use variables
  use parameters, only: read_parameter_header
  use initialization
  use metrics
  use conditions
  use equations
  use evolution
  use fluxes
  implicit none
  integer :: method,solver,i,rk,step,component,unit_file
  real*8 :: initial(neq),current(neq),expected(neq),reconverted(neq),v(3),w,h
  real*8 :: max_drift,max_residual,scale

  if (accelerator_enabled) then
    print *, 'SRHD CONSERVATION TESTS SKIPPED: CPU reference integration'
    stop
  end if
  open(newunit=unit_file,file='tests/srhd_conservation_test.par',status='replace')
  write(unit_file,'(A)') 'problem = convergence','nx = 64','ny = 1','nz = 1', &
    'reconstruction = mp5','riemann_solver = hlle','shock_sensor = true', &
    'numerical_diagnostics = true','final_time = 0.2','extract_gw = false', &
    'extract_mdot = false','output_folder = tests/conservation_output'
  close(unit_file)
  call read_parameter_header('tests/srhd_conservation_test.par')
  call initialize_problem('tests/srhd_conservation_test.par',.false.)
  dt=0.001d0
  v=(/0.35d0,0.2d0,-0.15d0/)
  w=1.0d0/sqrt(1.0d0-sum(v*v))
  h=1.0d0+g1
  expected(eq_de)=w
  expected(eq_pr)=h*w*w-1.0d0-w
  expected(eq_vx:eq_vz)=h*w*w*v
  scale=maxval(abs(expected))
  do method=REC_MP5,REC_WENO5
    rec_method_id=method
    do solver=RS_HLLE,RS_HLLC
      riemann_solver_id=solver
      call init_wavespeed_solver()
      numerical_events=0_int64
      numerical_delta_u=0.0d0
      do i=1,nx
        p(i,1,1,eq_de)=1.0d0+0.1d0*sin(2.0d0*pi*x(i))
        p(i,1,1,eq_pr)=1.0d0
        p(i,1,1,eq_vx:eq_vz)=v
        call prim_to_cons(p(i,1,1,:),up(i,1,1,:),i,1,1)
      end do
      call set_boundary_conditions(p)
      do component=1,neq
        initial(component)=sum(up(:,1,1,component))*dx
      end do
      if (maxval(abs(initial-expected))/scale>1.0d-13) error stop 'independent initial integrals'
      max_drift=0.0d0
      max_residual=0.0d0
      do step=1,200
        u=up
        do rk=1,3
          call calc_rhs(p,rhs)
          if (rk==1) then
            up=u+dt*rhs
          else if (rk==2) then
            up=0.75d0*u+0.25d0*(up+dt*rhs)
          else
            up=(u+2.0d0*(up+dt*rhs))/3.0d0
          end if
          do i=1,nx
            call cons_to_prim(up(i,1,1,:),p(i,1,1,:),i,1,1)
            call prim_to_cons(p(i,1,1,:),reconverted,i,1,1)
            max_residual=max(max_residual,maxval(abs(reconverted-up(i,1,1,:)))/scale)
          end do
          call set_boundary_conditions(p)
        end do
        do component=1,neq
          current(component)=sum(up(:,1,1,component))*dx
        end do
        max_drift=max(max_drift,maxval(abs(current-initial))/scale)
      end do
      write(*,'(A,2(I0,1X),2(ES24.16E3,1X))') '# conservation method/solver/drift/recovery ', &
        method,solver,max_drift,max_residual
      if (max_drift>1.0d-11.or.max_residual>1.0d-10) error stop 'periodic five-variable balance'
      if (sum(numerical_events(NUM_C2P_LOW_D:NUM_P2C_VELOCITY_CAP))/=0) error stop 'periodic state repaired'
      if (numerical_events(NUM_C2P_CALLS)/=3*200*nx) error stop 'periodic recovery count'
    end do
  end do
  print *, 'SRHD CONSERVATION TESTS PASSED: nonzero transverse momentum, 200 full RK3 steps'
end program test_srhd_conservation
