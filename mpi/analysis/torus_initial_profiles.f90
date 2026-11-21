! Static initial profiles using the same FD/FM evaluators as GRIAL; no evolution.
program torus_initial_profiles
  use variables
  use metrics
  use conditions
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  character(len=512) :: arg, filename
  character(len=2) :: model
  integer :: nr, nt, i, j, unit_file, ios, eq_inside
  real*8 :: rin, rc, K_poly, rlo, rhi, dlr, dtheta, r, theta, angular, win
  real*8 :: weight, volume, mass, body_volume, values(7), eq(7), means(7), integrals(7)
  real*8 :: eq_h, eq_cs2, scratch_h, scratch_cs2
  logical :: valid, inside

  if (command_argument_count()/=11) &
    error stop 'usage: torus_initial_profiles FD/FM spin rin rc K Gamma nr ntheta rlo rhi output.csv'
  call get_command_argument(1,model)
  if(model/='FD'.and.model/='FM') error stop 'model must be FD or FM'
  call real_argument(2,a_spin);call real_argument(3,rin);call real_argument(4,rc)
  call real_argument(5,K_poly);call real_argument(6,adb_idx)
  call get_command_argument(7,arg);read(arg,*,iostat=ios) nr
  if(ios/=0) error stop 'invalid nr'
  call get_command_argument(8,arg);read(arg,*,iostat=ios) nt
  if(ios/=0) error stop 'invalid ntheta'
  call real_argument(9,rlo);call real_argument(10,rhi)
  call get_command_argument(11,filename)
  if(nr<16.or.nt<8.or.K_poly<=0.or.adb_idx<=1.or.adb_idx>2.or.abs(a_spin)>1.or. &
    rlo<=0.or.rlo>=rin.or.rhi<=rc) error stop 'invalid profile parameters'
  bh_mass=1.0d0;g1=adb_idx/(adb_idx-1.0d0)
  rho_floor=1.0d-10;p_floor=1.0d-13
  metric_type='Kerr-Schild';geom_type='Spheroidal';use_log_r=.true.
  call set_metric_type()
  if(model=='FD')then
    call prepare_font_daigne(rin,rc,angular,win,valid)
  else
    call prepare_fishbone_moncrief(rin,rc,angular,win,valid)
  end if
  if(.not.valid) error stop 'no closed torus for requested parameters'
  open(newunit=unit_file,file=trim(filename),status='new',action='write')
  write(unit_file,'(A)') 'r,inside_eq,rho_eq,p_eq,Omega_eq,uphi_eq,ell_eq,L_eq,E_eq,h_eq,cs2_eq,'// &
    'dMb_dlnr,dVbody_dlnr,rho_mean_M,p_mean_M,Omega_mean_M,uphi_mean_M,ell_mean_M,L_mean_M,E_mean_M,'// &
    'rho_integral_V,p_integral_V,Omega_integral_V,uphi_integral_V,ell_integral_V,L_integral_V,E_integral_V'
  dlr=log(rhi/rlo)/nr;dtheta=pi/nt
  do i=1,nr
    r=rlo*exp((i-.5d0)*dlr)
    mass=0.0d0;body_volume=0.0d0;means=0.0d0;integrals=0.0d0
    do j=1,nt
      theta=(j-.5d0)*dtheta
      call sample(r,theta,values,scratch_h,scratch_cs2,weight,volume,inside)
      if(.not.inside) cycle
      mass=mass+weight*dtheta*2*pi
      body_volume=body_volume+volume*dtheta*2*pi
      means=means+values*weight*dtheta*2*pi
      integrals=integrals+values*volume*dtheta*2*pi
    end do
    if(mass>0) means=means/mass
    call sample(r,pi/2,eq,eq_h,eq_cs2,weight,volume,inside)
    eq_inside=merge(1,0,inside)
    write(unit_file,'(ES25.16E3,",",I0,25(",",ES25.16E3))') &
      r,eq_inside,eq,eq_h,eq_cs2,mass,body_volume,means,integrals
  end do
  close(unit_file)
  print *, 'STATIC TORUS PROFILES PASSED: ',model,nr,nt

contains

  subroutine real_argument(index,value)
    integer, intent(in) :: index
    real*8, intent(out) :: value
    call get_command_argument(index,arg)
    read(arg,*,iostat=ios) value
    if(ios/=0) error stop 'invalid real parameter'
    if(.not.ieee_is_finite(value)) error stop 'nonfinite real parameter'
  end subroutine

  subroutine sample(radius,angle,q,h,cs2,dMb,dV,inside)
    real*8, intent(in) :: radius,angle
    real*8, intent(out) :: q(7),h,cs2,dMb,dV
    logical, intent(out) :: inside
    real*8 :: state(neq),alpha,beta(3),g(3,3),detg,v(3),v2,Wlor,ut,ucovt,ucovphi,omega
    call calculate_metric(log(radius),angle,alpha=alpha,beta=beta,gamma=g,det=detg)
    if(model=='FD')then
      call evaluate_font_daigne_state(radius,rin,K_poly,angular,win,alpha,beta,g,state,inside)
    else
      call evaluate_fishbone_moncrief_state(radius,rin,K_poly,angular,win,alpha,beta,g,state,inside)
    end if
    q=0.0d0;h=0.0d0;cs2=0.0d0;dMb=0.0d0;dV=0.0d0
    if(.not.inside) return
    v=state(eq_vx:eq_vz);v2=dot_product(v,matmul(g,v))
    if(v2<0.or.v2>=1) error stop 'nonphysical native velocity'
    Wlor=1/sqrt(1-v2);ut=Wlor/alpha
    ucovt=Wlor*(-alpha+dot_product(beta,matmul(g,v)))
    ucovphi=Wlor*dot_product(g(3,:),v)
    omega=alpha*v(3)-beta(3)
    h=1+g1*state(eq_pr)/state(eq_de)
    cs2=adb_idx*state(eq_pr)/(state(eq_de)*h)
    q=[state(eq_de),state(eq_pr),omega,ut*omega,-ucovphi/ucovt,ucovphi*ut,-h*ucovt]
    ! sqrt(det) includes dr/dln(r)=r in the logarithmic metric.
    dV=sqrt(detg);dMb=state(eq_de)*Wlor*dV
    if(.not.all(ieee_is_finite(q)).or.cs2<=0.or.cs2>=1.or.dMb<=0) &
      error stop 'invalid native profile state'
  end subroutine
end program
