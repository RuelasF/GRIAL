program test_spatial_order
  use variables
  use reconstruction
  use fluxes
  implicit none
  integer, parameter :: levels(4)=(/16,32,64,128/)
  integer :: method,solver,level,i,n
  real*8, allocatable :: line(:,:),left(:,:),right(:,:),face_flux(:,:)
  real*8 :: g(3,3),beta(3),w,xcell,sinc,rhs_exact,error(4),order,required

  nghost=3
  metric_type='Minkowski'
  geom_type='Cartesian'
  use_shock_sensor=.false.
  do_numerical_diagnostics=.false.
  rho_floor=1.0d-10
  p_floor=1.0d-13
  adb_idx=4.0d0/3.0d0
  g1=adb_idx/(adb_idx-1.0d0)
  tvd_limiter_id=LIM_MC
  call init_wavespeed_solver()
  g=0.0d0
  do i=1,3
    g(i,i)=1.0d0
  end do
  beta=0.0d0
  w=1.0d0/sqrt(1.0d0-0.5d0**2)
  write(*,'(A)') '# method solver N L1_rhs_D observed_order'
  do method=REC_WENO3,REC_WENO5
    rec_method_id=method
    do solver=RS_HLLE,RS_HLLC
      riemann_solver_id=solver
      do level=1,size(levels)
        n=levels(level)
        allocate(line(neq,-3:n+3),left(neq,0:n),right(neq,0:n),face_flux(neq,0:n))
        sinc=sin(pi/dble(n))/(pi/dble(n))
        line=0.0d0
        do i=-3,n+3
          xcell=(dble(i)-0.5d0)/dble(n)
          line(eq_de,i)=1.0d0+0.1d0*sinc*sin(2.0d0*pi*xcell)
        end do
        line(eq_pr,:)=1.0d0
        line(eq_vx,:)=0.5d0
        call reconstruct_1d_core(line,n,left,right,DIR_X)
        do i=0,n
          call resolve_riemann_problem(left(:,i),right(:,i),DIR_X,1.0d0,beta,g,1.0d0,face_flux(:,i))
        end do
        error(level)=0.0d0
        do i=1,n
          xcell=(dble(i)-0.5d0)/dble(n)
          rhs_exact=-w*0.5d0*0.1d0*2.0d0*pi*sinc*cos(2.0d0*pi*xcell)
          error(level)=error(level)+abs(-dble(n)*(face_flux(eq_de,i)-face_flux(eq_de,i-1))-rhs_exact)/dble(n)
        end do
        order=0.0d0
        if (level>1) order=log(error(level-1)/error(level))/log(2.0d0)
        write(*,'(3(I0,1X),2(ES24.16E3,1X))') method,solver,n,error(level),order
        deallocate(line,left,right,face_flux)
      end do
      required=4.5d0
      if (method==REC_WENO3) required=2.5d0
      if (any(error(2:4)>=error(1:3)) .or. order<required) error stop 'spatial order below expected range'
    end do
  end do
  write(*,'(A)') '# SPATIAL ORDER TESTS PASSED: semidiscrete RHS, no time integrator'
end program test_spatial_order
