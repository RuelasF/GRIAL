program test_primitive_recovery
  use variables
  use equations
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
    ieee_positive_inf, ieee_negative_inf, ieee_quiet_nan, ieee_signaling_nan
  implicit none
  real*8, parameter :: gammas(4)=(/4.0d0/3.0d0,1.4d0,5.0d0/3.0d0,2.0d0/)
  real*8, parameter :: densities(9)=(/1.0d-9,1.0d-7,1.0d-5,1.0d-3,0.1d0,1.0d0,1.0d2,1.0d4,1.0d6/)
  real*8, parameter :: temperatures(7)=(/1.0d-6,1.0d-4,0.01d0,0.1d0,1.0d0,10.0d0,1.0d3/)
  real*8, parameter :: lorentz(7)=(/1.0d0,1.0001d0,1.05d0,2.0d0,10.0d0,100.0d0,300.0d0/)
  real*8 :: original(neq),recovered(neq),conserved(neq),before(neq),reconverted(neq),exact(neq)
  real*8 :: direction(3),metric_scale(3),speed,w,h,scale,eps,rho_error,p_error,v_error,residual
  real*8 :: rho_tol,p_tol,v_tol,max_errors(4),maximum_modification,change,conditioning,sound_speed_sq
  integer :: ig,ir,it,iw,id,im,component,cases,failures,unit_report,ios
  integer(int64) :: seed
  logical :: passed,report
  character(len=512) :: report_path

  if (accelerator_enabled) then
    print *, 'PRIMITIVE RECOVERY TESTS SKIPPED: CPU diagnostic counters required'
    stop
  end if
  rho_floor=1.0d-10
  p_floor=1.0d-13
  metric_type='Minkowski'
  geom_type='Cartesian'
  do_numerical_diagnostics=.true.
  numerical_events=0_int64
  numerical_delta_u=0.0d0
  allocate(gamma_c(3,3,1,1,1),gamma_inv_c(3,3,1,1,1),sqrt_gamma_c(1,1,1))
  call get_command_argument(1,report_path)
  report=len_trim(report_path)>0
  if (report) then
    open(newunit=unit_report,file=trim(report_path),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'cannot open recovery report'
    write(unit_report,'(A)') 'gamma metric rho p W direction rho_rel p_abs v_abs U_scaled U_modified passed'
  end if
  eps=epsilon(1.0d0)
  seed=20261011_int64
  cases=0
  failures=0
  max_errors=0.0d0
  maximum_modification=0.0d0
  do im=1,2
    metric_scale=1.0d0
    if (im==2) metric_scale=(/1.0d0,2.0d0,3.0d0/)
    gamma_c=0.0d0
    gamma_inv_c=0.0d0
    do component=1,3
      gamma_c(component,component,1,1,1)=metric_scale(component)**2
      gamma_inv_c(component,component,1,1,1)=1.0d0/metric_scale(component)**2
    end do
    sqrt_gamma_c=product(metric_scale)
    do ig=1,size(gammas)
      adb_idx=gammas(ig)
      g1=adb_idx/(adb_idx-1.0d0)
      do ir=1,size(densities)
        do it=1,size(temperatures)
          original(eq_de)=densities(ir)
          original(eq_pr)=densities(ir)*temperatures(it)
          if (original(eq_pr)<100.0d0*p_floor) cycle
          do iw=1,size(lorentz)
            speed=sqrt(1.0d0-1.0d0/lorentz(iw)**2)
            do id=1,7
              direction=0.0d0
              if (id<=3) then
                direction(id)=1.0d0
              else
                do component=1,3
                  seed=modulo(16807_int64*seed,2147483647_int64)
                  direction(component)=2.0d0*dble(seed)/2147483647.0d0-1.0d0
                end do
                direction=direction/sqrt(sum(direction**2))
              end if
              original(eq_vx:eq_vz)=speed*direction/metric_scale
              w=1.0d0/sqrt(1.0d0-sum((original(eq_vx:eq_vz)*metric_scale)**2))
              h=1.0d0+g1*original(eq_pr)/original(eq_de)
              exact(eq_de)=original(eq_de)*w
              exact(eq_pr)=original(eq_de)*h*w*w-original(eq_pr)-exact(eq_de)
              exact(eq_vx:eq_vz)=original(eq_de)*h*w*w*original(eq_vx:eq_vz)*metric_scale**2
              exact=exact*product(metric_scale)
              call prim_to_cons(original,conserved,1,1,1)
              scale=maxval(abs(exact))
              if (maxval(abs(conserved-exact))/scale>64.0d0*eps*w*w) error stop 'independent P2C formula'
              before=conserved
              call cons_to_prim(conserved,recovered,1,1,1)
              call prim_to_cons(recovered,reconverted,1,1,1)
              rho_error=abs(recovered(eq_de)-original(eq_de))/original(eq_de)
              p_error=abs(recovered(eq_pr)-original(eq_pr))
              v_error=maxval(abs((recovered(eq_vx:eq_vz)-original(eq_vx:eq_vz))*metric_scale))
              residual=maxval(abs(reconverted-before))/scale
              change=maxval(abs(conserved-before))/scale
              ! Inversion becomes ill-conditioned when the sound speed approaches light speed.
              sound_speed_sq=adb_idx*original(eq_pr)/(original(eq_de)*h)
              conditioning=w*w/(1.0d0-sound_speed_sq*(1.0d0-1.0d0/(w*w)))
              rho_tol=max(2.0d-9,64.0d0*eps*conditioning)
              p_tol=max(2.0d-9*original(eq_pr),64.0d0*eps*conditioning*(original(eq_de)+g1*original(eq_pr)))
              v_tol=max(2.0d-10,128.0d0*eps*w*w)
              passed=all(ieee_is_finite(recovered)).and.rho_error<=rho_tol.and.p_error<=p_tol.and. &
                v_error<=v_tol.and.residual<=2.0d-10.and.change==0.0d0
              cases=cases+1
              max_errors=max(max_errors,(/rho_error,p_error/max(p_tol,tiny(1.0d0)),v_error,residual/))
              maximum_modification=max(maximum_modification,change)
              if (.not.passed) then
                failures=failures+1
                if (failures<=12) write(*,'(A,4(ES12.4,1X),A,4(ES12.4,1X))') &
                  'FAIL Gamma/rho/p/W ',adb_idx,original(eq_de),original(eq_pr),w, &
                  ' errors rho/p_tol/v/U ',rho_error,p_error/p_tol,v_error,residual
              end if
              if (report) write(unit_report,'(ES24.16E3,1X,I0,1X,3(ES24.16E3,1X),I0,1X,5(ES24.16E3,1X),L1)') &
                adb_idx,im,original(eq_de),original(eq_pr),w,id,rho_error,p_error,v_error,residual,change,passed
            end do
          end do
        end do
      end do
    end do
  end do
  if (report) close(unit_report)
  write(*,'(A,I0,A,I0,A,5(ES12.4,1X))') 'RECOVERY states=',cases,' failures=',failures, &
    ' max rho_rel/p_tol/v/U/modified ',max_errors,maximum_modification
  if (numerical_events(NUM_C2P_CALLS)/=cases) error stop 'recovery calls not counted'
  if (sum(numerical_events(NUM_C2P_LOW_D:NUM_P2C_VELOCITY_CAP))/=0) then
    print *, 'Valid-state repair counters:',numerical_events(NUM_C2P_LOW_D:NUM_P2C_VELOCITY_CAP)
    failures=failures+1
  end if
  if (failures>0) error stop 'valid primitive recovery failed'

  ! Invalid low-density and excessive-momentum states exercise deliberate repairs separately.
  numerical_events=0_int64
  numerical_delta_u=0.0d0
  conserved=0.0d0
  conserved(eq_de)=0.5d0*rho_floor*sqrt_gamma_c(1,1,1)
  before=conserved
  call cons_to_prim(conserved,recovered,1,1,1)
  if (numerical_events(NUM_C2P_LOW_D)/=1) error stop 'low-D repair missing'
  if (maxval(abs(numerical_delta_u-(conserved-before)))>1.0d-25) error stop 'repair delta mismatch'
  numerical_events=0_int64
  conserved=(/1.0d0,1.0d0,100.0d0,0.0d0,0.0d0/)*sqrt_gamma_c(1,1,1)
  call cons_to_prim(conserved,recovered,1,1,1)
  if (numerical_events(NUM_C2P_MOMENTUM)/=1) error stop 'invalid momentum not repaired'
  if (.not.all(ieee_is_finite(recovered)).or.minval(recovered(1:2))<=0.0d0) error stop 'invalid repair output'
  call check_nonfinite_input()
  print *, 'PRIMITIVE RECOVERY TESTS PASSED; deliberate repairs tested separately'

contains

  subroutine check_nonfinite_input()
    real*8 :: invalid(4), initial(neq), atmosphere(neq), expected_delta(neq)
    integer :: field, kind

    invalid(1)=ieee_value(0.0d0,ieee_positive_inf)
    invalid(2)=ieee_value(0.0d0,ieee_negative_inf)
    invalid(3)=ieee_value(0.0d0,ieee_quiet_nan)
    invalid(4)=ieee_value(0.0d0,ieee_signaling_nan)
    original=[1.0d0,0.2d0,0.1d0,-0.05d0,0.025d0]
    call prim_to_cons(original,initial,1,1,1)
    atmosphere=0.0d0
    atmosphere(eq_de)=rho_floor*sqrt_gamma_c(1,1,1)
    atmosphere(eq_pr)=p_floor*(g1-1.0d0)*sqrt_gamma_c(1,1,1)
    do kind=1,size(invalid)
      do field=1,neq
        conserved=initial
        conserved(field)=invalid(kind)
        numerical_events=0_int64
        numerical_delta_u=0.0d0
        call cons_to_prim(conserved,recovered,1,1,1)
        if (.not.all(ieee_is_finite(conserved)).or..not.all(ieee_is_finite(recovered))) &
          error stop 'nonfinite input did not reach atmosphere'
        if (any(conserved/=atmosphere)) error stop 'nonfinite input: wrong conservative atmosphere'
        if (recovered(eq_de)/=rho_floor.or.recovered(eq_pr)/=p_floor.or.any(recovered(eq_vx:eq_vz)/=0.0d0)) &
          error stop 'nonfinite input: wrong primitive atmosphere'
        if (numerical_events(NUM_C2P_CALLS)/=1.or.numerical_events(NUM_C2P_SANITIZE)/=1.or. &
            sum(numerical_events)/=2) error stop 'nonfinite input: wrong recovery counters'
        expected_delta=atmosphere-initial
        expected_delta(field)=0.0d0
        if (.not.all(ieee_is_finite(numerical_delta_u))) error stop 'nonfinite input contaminated the budget'
        if (any(numerical_delta_u/=expected_delta)) error stop 'nonfinite input: wrong finite-component budget'
      end do
    end do
    print '(A,I0)', 'NONFINITE RECOVERY CASES PASSED: ',neq*size(invalid)
  end subroutine check_nonfinite_input
end program test_primitive_recovery
