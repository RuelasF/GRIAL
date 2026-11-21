! Prescribed rings exercise diagnostic routines only; no fluid evolution.
program test_ppi_gw_diagnostics
  use variables
  use metrics, only: calculate_metric,set_metric_type
  use equations, only: prim_to_cons
  use output, only: calc_ppi_modes,calc_gw_strain,gw_observer_distance
  implicit none
  integer :: i,k,n,coord,spin_index,ntimes,unit_ref,unit_modes,unit_gw,m,base,unit_actual
  real*8 :: r,alpha,beta(3),g(3,3),determinant,sample_time,omega,growth,a1,a2,phase1,phase2,rho_base
  real*8 :: c0_ref,mass_weight,mode_weight,v_r,v_phi,v_squared,lorentz,radial_jac
  real*8 :: expected(2+4*PPI_MAX_MODE),actual(2+4*PPI_MAX_MODE),wave(3),wave_ref(3),gw_scale
  complex*16 :: coefficients(PPI_MAX_MODE),wave_weight,polarizations
  character(len=512) :: folder,arg
  character(len=32) :: label
  logical :: exists

  call get_command_argument(1,folder)
  if(len_trim(folder)==0) error stop 'Provide a unique output directory'
  if(verify(trim(folder),'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789/_-.')/=0) &
    error stop 'Fixture output path requires ASCII letters, numbers, /, _, - or .'
  ntimes=65
  if(command_argument_count()>=2) then
    call get_command_argument(2,arg);read(arg,*) ntimes
  end if
  if(ntimes<8.or.ntimes>1025) error stop 'Invalid sample count'
  nx=2;ny=1;nz=64;bh_mass=1;g1=4;rho_floor=1.0d-10;p_floor=1.0d-13
  dy=pi;dz=2*pi/nz;omega=2*pi/256.0d0;growth=0.003d0
  output_folder=trim(folder);do_ppi_diagnostics=.true.;do_numerical_diagnostics=.false.
  allocate(x(nx),y(ny),z(nz),p(nx,ny,nz,neq),up(nx,ny,nz,neq))
  allocate(alpha_c(nx,ny,nz),beta_c(3,nx,ny,nz),gamma_c(3,3,nx,ny,nz),sqrt_gamma_c(nx,ny,nz))
  y=pi/2
  do k=1,nz
    z(k)=(k-0.5d0)*dz
  end do
  do spin_index=0,1
    a_spin=0.9d0*spin_index;metric_type='Kerr-Schild';geom_type='Spheroidal'
    do coord=0,1
      use_log_r=coord==1
      if(use_log_r) then
        x=log([6.0d0,8.0d0]);dx=log(8.0d0/6.0d0)
      else
        x=[6.0d0,8.0d0];dx=2
      end if
      write(label,'("a",I1,"_coord",I1)') spin_index,coord
      scheme_name=trim(label)
      inquire(file=trim(folder)//'/'//trim(label)//'/expected.dat',exist=exists)
      if(exists) error stop 'Refusing to overwrite diagnostic evidence'
      call execute_command_line('mkdir -p '//trim(folder)//'/'//trim(label))
      call set_metric_type()
      c0_ref=0;mode_weight=0;wave_weight=cmplx(0.0d0,0.0d0,kind=8)
      do i=1,nx
        r=6.0d0+2*(i-1);radial_jac=1
        if(use_log_r) radial_jac=r
        call calculate_metric(x(i),y(1),alpha=alpha,beta=beta,gamma=g,det=determinant)
        do k=1,nz
          alpha_c(i,1,k)=alpha;beta_c(:,i,1,k)=beta;gamma_c(:,:,i,1,k)=g
          sqrt_gamma_c(i,1,k)=sqrt(determinant)
        end do
        ! Independent equatorial KS formula for W and proper-volume weights.
        v_r=(2/r)/sqrt(1+2/r);v_phi=omega*sqrt(1+2/r)
        v_squared=(1+2/r)*v_r**2-2*a_spin*(1+2/r)*v_r*v_phi &
          +(r*r+a_spin*a_spin+2*a_spin*a_spin/r)*v_phi**2
        if(v_squared>=1) error stop 'Fixture is not causal'
        lorentz=1/sqrt(1-v_squared)
        mass_weight=r*r*sqrt(1+2/r)*radial_jac*dx*dy*lorentz
        rho_base=exp(-((r-7)/2)**2)
        c0_ref=c0_ref+2*pi*(rho_base+rho_floor)*mass_weight
        mode_weight=mode_weight+2*pi*rho_base*mass_weight
        ! Analytic azimuthal integral of the direct stress formula, not its cell helper.
        wave_weight=wave_weight-4*pi*rho_base*(omega**2+1/(r*r+a_spin*a_spin)**1.5d0) &
          *cmplx(r,a_spin,kind=8)**2*r*r*radial_jac*dx*dy/gw_observer_distance
      end do
      gw_scale=abs(wave_weight)*0.025d0
      open(newunit=unit_ref,file=trim(folder)//'/'//trim(label)//'/expected.dat',status='new')
      do n=0,ntimes-1
        sample_time=dble(n);a1=1.0d-4*exp(growth*sample_time);a2=0.025d0
        phase1=omega*sample_time+0.3d0;phase2=omega*sample_time+0.7d0
        do k=1,nz;do i=1,nx
          r=6.0d0+2*(i-1);rho_base=exp(-((r-7)/2)**2)
          p(i,1,k,eq_de)=rho_base*(1+2*a1*cos(z(k)-phase1)+2*a2*cos(2*(z(k)-phase2)))+rho_floor
          p(i,1,k,eq_pr)=1.0d-3+p_floor
          p(i,1,k,eq_vx:eq_vz)=beta_c(:,i,1,k)/alpha_c(i,1,k)
          p(i,1,k,eq_vz)=p(i,1,k,eq_vz)+omega/alpha_c(i,1,k)
          call prim_to_cons(p(i,1,k,:),up(i,1,k,:),i,1,k)
        end do;end do
        call calc_ppi_modes(sample_time);call calc_gw_strain(sample_time)
        coefficients=cmplx(0.0d0,0.0d0,kind=8)
        coefficients(1)=mode_weight*a1*exp(cmplx(0.0d0,-phase1,kind=8))
        coefficients(2)=mode_weight*a2*exp(cmplx(0.0d0,-2*phase2,kind=8))
        polarizations=wave_weight*a2*exp(cmplx(0.0d0,2*phase2,kind=8))
        expected=0;expected(1)=sample_time;expected(2)=c0_ref
        do m=1,PPI_MAX_MODE
          base=3+4*(m-1)
          expected(base)=real(coefficients(m),kind=8);expected(base+1)=aimag(coefficients(m))
          expected(base+2)=abs(coefficients(m))/c0_ref
          if(m<=2) expected(base+3)=-atan2(aimag(coefficients(m)),real(coefficients(m),kind=8))/m
        end do
        write(unit_ref,'(*(ES24.16E3,1X))') expected,real(polarizations,kind=8),aimag(polarizations),gw_scale
      end do
      close(unit_ref)
      open(newunit=unit_ref,file=trim(folder)//'/'//trim(label)//'/expected.dat',status='old')
      open(newunit=unit_modes,file=trim(folder)//'/'//trim(label)//'/ppi_modes.dat',status='old')
      open(newunit=unit_gw,file=trim(folder)//'/'//trim(label)//'/GW_signal.dat',status='old')
      do n=0,ntimes-1
        read(unit_ref,*) expected,wave_ref(2:3),gw_scale
        read(unit_modes,*) actual;read(unit_gw,*) wave
        if(abs(actual(2)-expected(2))/c0_ref>1.0d-11) error stop 'PPI mass weight failed'
        do m=1,PPI_MAX_MODE
          base=3+4*(m-1)
          if(maxval(abs(actual(base:base+1)-expected(base:base+1)))/c0_ref>1.0d-12) &
            error stop 'PPI coefficient sign/normalization failed'
          if(abs(actual(base+2)-expected(base+2))>1.0d-12) error stop 'PPI normalized amplitude failed'
          if(m<=2.and.2*abs(sin(0.5d0*m*(actual(base+3)-expected(base+3))))>1.0d-11) error stop 'PPI phase failed'
        end do
        if(maxval(abs(wave(2:3)-wave_ref(2:3)))/gw_scale>5.0d-7) error stop 'GW integral/serialization failed'
      end do
      close(unit_ref);close(unit_modes);close(unit_gw)
      scheme_name=trim(label)//'_axisym'
      call execute_command_line('mkdir -p '//trim(folder)//'/'//trim(scheme_name))
      do k=1,nz;do i=1,nx
        r=6.0d0+2*(i-1)
        p(i,1,k,eq_de)=exp(-((r-7)/2)**2)+rho_floor
        call prim_to_cons(p(i,1,k,:),up(i,1,k,:),i,1,k)
      end do;end do
      call calc_ppi_modes(0.0d0);call calc_gw_strain(0.0d0)
      call check_zero_wave(gw_scale*1.0d-12)
      open(newunit=unit_actual,file=trim(folder)//'/'//trim(scheme_name)//'/ppi_modes.dat',status='old')
      read(unit_actual,*) actual;close(unit_actual)
      do m=1,PPI_MAX_MODE
        base=3+4*(m-1)
        if(maxval(abs(actual(base:base+1)))/c0_ref>1.0d-12) error stop 'Axisymmetric mode control failed'
      end do
      ! Complete atmosphere must subtract to zero; all interior matter is excluded.
      scheme_name=trim(label)//'_atmosphere'
      call execute_command_line('mkdir -p '//trim(folder)//'/'//trim(scheme_name))
      p(:,:,:,eq_de)=rho_floor;p(:,:,:,eq_pr)=p_floor
      call calc_gw_strain(0.0d0);call calc_ppi_modes(0.0d0)
      call check_zero_wave()
      open(newunit=unit_actual,file=trim(folder)//'/'//trim(scheme_name)//'/ppi_modes.dat',status='old')
      read(unit_actual,*) actual;close(unit_actual)
      if(maxval(abs(actual))>0) error stop 'Empty PPI mask control failed'
      scheme_name=trim(label)//'_horizon'
      call execute_command_line('mkdir -p '//trim(folder)//'/'//trim(scheme_name))
      x=1.25d0
      if(use_log_r) x=log(1.25d0)
      p(:,:,:,eq_de)=1;p(:,:,:,eq_pr)=1.0d-3
      call calc_gw_strain(0.0d0)
      call check_zero_wave()
    end do
  end do
  print *, 'PPI GW DIAGNOSTIC TESTS PASSED'
contains
  subroutine check_zero_wave(tolerance)
    real*8, intent(in), optional :: tolerance
    real*8 :: limit
    limit=0
    if(present(tolerance)) limit=tolerance
    open(newunit=unit_actual,file=trim(folder)//'/'//trim(scheme_name)//'/GW_signal.dat',status='old')
    read(unit_actual,*) wave;close(unit_actual)
    if(maxval(abs(wave(2:3)))>limit) error stop 'GW symmetry/mask/atmosphere zero control failed'
  end subroutine check_zero_wave
end program test_ppi_gw_diagnostics
