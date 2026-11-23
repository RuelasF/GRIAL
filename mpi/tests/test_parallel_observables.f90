program test_parallel_observables
  use variables
  use parallel_runtime
  use parallel_timing, only: clip_step_to_target
  use parallel_diagnostics
  use domain_decomposition
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use metrics
  use conditions, only: set_initial_conditions,set_boundary_conditions,prepare_case_boundary_data
  use equations, only: prim_to_cons
  use fluxes, only: init_wavespeed_solver
  use evolution, only: update_adaptive_dt
  use halo_exchange
  use time_integrator, only: advance_rk_step
  use output, only: calc_ppi_modes,calc_gw_strain,calc_m_dot,ppi_mode_values, &
    write_distributed_observables,gw_observer_distance,write_distributed_diagnostics
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer :: grid(2)=[1,1],fixture=1,checks=0,cases=0,unit,n,shape_global(3),ios
  type(domain_layout) :: layout
  character(len=64) :: arg,mode,label
  real*8 :: field_error=0.0d0,observable_error=0.0d0,analytic_error=0.0d0
  real*8 :: coefficients(1+2*PPI_MAX_MODE),stress(3),mdot,reference(1+2*PPI_MAX_MODE+4,12)
  real*8 :: packet(1+2*PPI_MAX_MODE+4),values(2+4*PPI_MAX_MODE),expected(1+2*PPI_MAX_MODE+4)
  type(state_diagnostics) :: global_state,root_state
  type(conservation_balance) :: global_budget,root_budget

  call parallel_initialize()
  call get_command_argument(1,arg)
  if (len_trim(arg)>0) read(arg,*)grid(1)
  call get_command_argument(2,arg)
  if (len_trim(arg)>0) read(arg,*)grid(2)
  call get_command_argument(3,mode)
  if (len_trim(mode)>0.and.trim(mode)/='report') call invalid_call()
  do fixture=1,5
    call configure(.true.)
    call allocate_and_grid();call precalculate_metric_cache();call choose_probe()
    do n=1,12
      call prescribe(n,expected)
      call calc_ppi_modes(integration_time,coefficients)
      call calc_gw_strain(integration_time,stress);call calc_m_dot(integration_time,mdot)
      reference(:,n)=[coefficients,stress,mdot]
      call compare(reference(:,n),expected,'analytic prescribed coefficients, stress and accretion')
      if (parallel_rank==0) then
        write(label,'(A,I0)')'synthetic_',fixture
        call directory(trim(label)//'/serial')
        call calc_ppi_modes(integration_time);call calc_gw_strain(integration_time);call calc_m_dot(integration_time)
        call record_reference(reference(:,n))
      end if
    end do
    call release_grid()
    call configure(.true.);call make_layout(layout)
    call allocate_and_grid(layout);call precalculate_metric_cache();call choose_probe()
    call prescribe(1,expected);call initialize_parallel_diagnostics(layout)
    do n=1,12
      call prescribe(n,expected)
      ! Poison all ghosts: only physical cells belong in the observables.
      call poison_ghosts()
      call calc_ppi_modes(integration_time,coefficients)
      call calc_gw_strain(integration_time,stress);call calc_m_dot(integration_time,mdot)
      packet=[coefficients,stress,mdot];call reduce_observable_sum(packet)
      if (parallel_rank==0) call compare(packet,reference(:,n),'distributed prescribed observables')
      call collect_global_state(global_state)
      call collect_global_state(root_state,root_only=.true.)
      if (parallel_rank==0) then
        call require(maxval(abs(root_state%totals-global_state%totals))/maxval(abs(global_state%totals))<1.0d-11, &
          'root conservation totals equal all-rank totals')
        call require(maxval(abs(root_state%l2-global_state%l2))/maxval(global_state%l2)<1.0d-11, &
          'root L2 retains shared scale')
        call require(all(root_state%events==global_state%events).and.root_state%cells==global_state%cells, &
          'root exact int64 counters')
      else
        call require(root_state%cells==0.and.all(root_state%totals==0.0d0),'root-only result undefined elsewhere is zero')
      end if
      write(label,'(A,I0)')'synthetic_',fixture
      call directory(trim(label)//'/distributed')
      call write_distributed_observables(integration_time)
    end do
    call finalize_parallel_diagnostics();call release_grid();cases=cases+1
  end do
  do fixture=1,8
    call check_evolution()
  end do
  write(label,'(A,I0,A)')'observables_rank_',parallel_rank,'.txt'
  open(newunit=unit,file=trim(label),status='new')
  write(unit,'(4(I0,1X),3(ES24.16E3,1X))')parallel_rank,parallel_size,checks,cases,field_error,observable_error,analytic_error
  close(unit)
  if (parallel_rank==0) print *,'PARALLEL OBSERVABLES TESTS PASSED'
  call parallel_finalize()
contains
  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    checks=checks+1
    if (.not.ok) call parallel_abort('D02 test: '//message)
  end subroutine require

  subroutine configure(synthetic)
    logical, intent(in) :: synthetic
    shape_global=[37,19,32];case_name='FontDaigne3D';geom_type='Spheroidal';metric_type='Kerr-Schild'
    r_min=1.2d0;r_max=14.0d0;y_min=0;y_max=pi;z_min=0;z_max=2*pi
    bh_mass=1;a_spin=0;use_log_r=.false.
    if (synthetic) then
      if (fixture==2.or.fixture==4)a_spin=0.9d0
      use_log_r=fixture>=3
      if (fixture==5)then
        metric_type='Eddington-Finkelstein';geom_type='Spherical'
      end if
    else
      shape_global=[37,19,8];r_max=40.0d0;use_log_r=.true.
      if (mod(fixture,2)==0)a_spin=0.9d0
      if (mod((fixture-1)/2,2)==1)case_name='FishMonc3D'
      if (fixture>4)then
        shape_global(3)=1
        if (trim(case_name)=='FontDaigne3D')case_name='FontDaigneSag'
        if (trim(case_name)=='FishMonc3D')case_name='FishMoncSag'
      end if
    end if
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=3
    adb_idx=4.0d0/3.0d0;g1=4;CFL=0.3d0
    rho_floor=1.0d-10;p_floor=1.0d-13;D_floor=rho_floor;tau_floor=3*p_floor
    do_mdot_extraction=.true.;do_gw_extraction=.true.;do_ppi_diagnostics=.true.
    do_numerical_diagnostics=.true.;do_radial_mass_budget=.false.;do_restart=.false.
    apply_perturbation=.false.;perturbation_applied=.false.
    rec_method_id=REC_MP5
    if (mod(fixture,2)==0)rec_method_id=REC_WENO5
    riemann_solver_id=RS_HLLE
    if (mod((fixture-1)/2,2)==1)riemann_solver_id=RS_HLLC
    use_shock_sensor=.true.;shock_sensor_thresholds=SENSOR_DEFAULTS;tvd_limiter_id=LIM_MC
    integration_time=0;n_steps=0;numerical_events=0_int64;numerical_delta_u=0
    call set_metric_type();call init_wavespeed_solver()
  end subroutine configure

  subroutine make_layout(block)
    type(domain_layout), intent(out) :: block
    character(len=180) :: message
    call build_domain(shape_global,nghost,parallel_rank,parallel_size,.false.,block,message,grid,.false.)
    call parallel_require_collective(len_trim(message)==0,message)
  end subroutine make_layout

  subroutine choose_probe()
    idx_probe_global=shape_global(1)/2
    idx_probe=0
    if (idx_probe_global>x_global_offset.and.idx_probe_global<=x_global_offset+nx)idx_probe=idx_probe_global-x_global_offset
  end subroutine choose_probe

  subroutine directory(name)
    character(len=*), intent(in) :: name
    logical :: ok
    output_folder=trim(name);scheme_name='d02';ok=.true.
    if (parallel_rank==0) then
      call execute_command_line('mkdir -p '//trim(name)//'/d02',exitstat=ios);ok=ios==0
    end if
    ! Serial reference is called only on root; mkdir is synchronous there.
    if (.not.ok)call parallel_abort('D02 test output directory')
  end subroutine directory

  subroutine record_reference(integrals)
    real*8, intent(in) :: integrals(:)
    integer :: unit_ref
    character(len=512) :: path
    path=trim(output_folder)//'/d02/reference_integrals.dat'
    open(newunit=unit_ref,file=trim(path),status='unknown',position='append')
    write(unit_ref,'(*(ES24.16E3,1X))')integration_time,integrals
    close(unit_ref)
  end subroutine record_reference

  subroutine prescribe(sample,analytic)
    integer, intent(in) :: sample
    real*8, intent(out) :: analytic(:)
    integer :: i,j,k,m
    real*8 :: r,th,sigma,f,rho_base,omega,growth,a1,a2,phase1,phase2,radial_jac,vr,weight,cut,v2
    real*8 :: base_mass,disk_mass,mdot_ref,mode_mass,lorentz,horizon
    complex*16 :: wave_weight,qv,qx,polarization,cm
    integration_time=32.0d0*(sample-1);n_steps=sample-1
    omega=2*pi/256.0d0;growth=0.003d0;a1=1.0d-4*exp(growth*integration_time);a2=0.025d0
    phase1=omega*integration_time+0.3d0;phase2=omega*integration_time+0.7d0
    if (sample>=10)then
      a1=0;a2=0
    end if
    horizon=event_horizon_radius(bh_mass,a_spin)
    base_mass=0;disk_mass=0;mdot_ref=0;wave_weight=cmplx(0.0d0,0.0d0,kind=8)
    do j=1,ny;do i=1,nx
      r=x(i);if (use_log_r)r=exp(r)
      th=y(j);radial_jac=1;if (use_log_r)radial_jac=r
      rho_base=0
      if (r>4.and.r<10.and.th>0.2d0.and.th<pi-0.2d0)rho_base=1+0.2d0*sin(th)
      if (sample==11)rho_base=0
      if (sample==12)then
        rho_base=0
        if (r<horizon)rho_base=1
      end if
      vr=-0.004d0
      p(i,j,1,eq_vx:eq_vz)=beta_c(:,i,j)/alpha_c(i,j)
      p(i,j,1,eq_vx)=p(i,j,1,eq_vx)+vr/(radial_jac*alpha_c(i,j))
      p(i,j,1,eq_vz)=p(i,j,1,eq_vz)+omega/alpha_c(i,j)
      ! Use a causal Eulerian rest state in the strong-field interior.
      if (r<=3.0d0)then
        p(i,j,1,eq_vx:eq_vz)=0.0d0;vr=-beta_c(1,i,j)*radial_jac
      end if
      v2=dot_product(p(i,j,1,eq_vx:eq_vz),matmul(gamma_c(:,:,i,j),p(i,j,1,eq_vx:eq_vz)))
      call require(v2<1.and.v2>=0,'causal prescribed flow')
      lorentz=1/sqrt(1-v2)
      weight=sqrt_gamma_c(i,j)*lorentz*dx*dy*2*pi
      if (rho_base>0)then
        base_mass=base_mass+rho_base*weight
        disk_mass=disk_mass+(rho_base+rho_floor)*weight
      end if
      if (i==idx_probe)mdot_ref=mdot_ref-(rho_base+rho_floor)*weight*vr/(radial_jac*dx)
      sigma=r*r+a_spin*a_spin*cos(th)**2
      if (r>horizon)then
        qx=cmplx(r,a_spin,kind=8)*sin(th)
        qv=(vr+cmplx(0.0d0,omega,kind=8)*cmplx(r,a_spin,kind=8))*sin(th)
        wave_weight=wave_weight+4*pi*rho_base*(qv**2-qx**2/(r*r+a_spin*a_spin*sin(th)**2)**1.5d0) &
          *sigma*sin(th)*radial_jac*dx*dy/gw_observer_distance
      end if
      do k=1,nz
        p(i,j,k,eq_vx:eq_vz)=p(i,j,1,eq_vx:eq_vz)
        p(i,j,k,eq_de)=rho_base*(1+2*a1*cos(z(k)-phase1)+2*a2*cos(2*(z(k)-phase2)))+rho_floor
        p(i,j,k,eq_pr)=1.0d-3+p_floor
        call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
      end do
    end do;end do
    analytic=0;analytic(1)=disk_mass
    cm=base_mass*a1*exp(cmplx(0.0d0,-phase1,kind=8));analytic(2)=real(cm,kind=8);analytic(2+PPI_MAX_MODE)=aimag(cm)
    cm=base_mass*a2*exp(cmplx(0.0d0,-2*phase2,kind=8));analytic(3)=real(cm,kind=8);analytic(3+PPI_MAX_MODE)=aimag(cm)
    polarization=wave_weight*a2*exp(cmplx(0.0d0,2*phase2,kind=8))
    ! Compare TT combinations below; individual stress components include isotropic pressure.
    analytic(2+2*PPI_MAX_MODE)=real(polarization,kind=8)
    analytic(3+2*PPI_MAX_MODE)=aimag(polarization)
    analytic(5+2*PPI_MAX_MODE)=mdot_ref
  end subroutine prescribe

  subroutine compare(actual,reference_values,message)
    real*8, intent(in) :: actual(:),reference_values(:)
    character(len=*), intent(in) :: message
    real*8 :: error,scale,wave_actual(2)
    integer :: q
    q=1+2*PPI_MAX_MODE;scale=max(abs(reference_values(1)),1.0d-10)
    if (index(message,'analytic')>0)then
      wave_actual=[(actual(q+1)-actual(q+2))/gw_observer_distance,2*actual(q+3)/gw_observer_distance]
      error=max(maxval(abs(actual(1:q)-reference_values(1:q)))/scale, &
        maxval(abs(wave_actual-reference_values(q+1:q+2)))/scale,abs(actual(q+4)-reference_values(q+4))/scale)
      analytic_error=max(analytic_error,error)
    else
      error=maxval(abs(actual-reference_values))/max(maxval(abs(reference_values)),1.0d-10)
      observable_error=max(observable_error,error)
    end if
    if (error>=1.0d-11)write(*,*)trim(message),fixture,n,parallel_rank,error,actual,reference_values
    call require(error<1.0d-11,message)
  end subroutine compare

  subroutine poison_ghosts()
    integer :: i,j,k
    do k=lbound(p,3),ubound(p,3);do j=lbound(p,2),ubound(p,2);do i=lbound(p,1),ubound(p,1)
      if (i<1.or.i>nx.or.j<1.or.j>ny.or.k<1.or.k>nz)p(i,j,k,:)=ieee_value(0.0d0,ieee_quiet_nan)
    end do;end do;end do
  end subroutine poison_ghosts

  subroutine check_evolution()
    real*8, allocatable :: initial_p(:,:,:,:),initial_u(:,:,:,:),final_p(:,:,:,:),final_u(:,:,:,:)
    real*8 :: expected_packet(1+2*PPI_MAX_MODE+4),ref_dt(512),local_dt(512),error
    integer :: ref_steps,step
    call configure(.false.);call allocate_and_grid();call precalculate_metric_cache();call set_initial_conditions()
    call set_boundary_conditions(p)
    block
      integer :: i,j,k
      do k=1,nz;do j=1,ny;do i=1,nx
        call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
      end do;end do;end do
    end block
    initial_p=p(1:nx,1:ny,1:nz,:);initial_u=up
    final_time=0.1d0;integration_time=0;n_steps=0;numerical_events=0_int64;numerical_delta_u=0
    do while(integration_time<final_time)
      call update_adaptive_dt();dt=clip_step_to_target(integration_time,dt,final_time)
      n_steps=n_steps+1;ref_dt(n_steps)=dt;integration_time=integration_time+dt;call advance_rk_step()
    end do
    ref_steps=n_steps;final_p=p(1:nx,1:ny,1:nz,:);final_u=up
    call calc_ppi_modes(integration_time,coefficients);call calc_gw_strain(integration_time,stress)
    call calc_m_dot(integration_time,mdot);expected_packet=[coefficients,stress,mdot]
    if (parallel_rank==0)then
      write(label,'(A,I0)')'evolution_',fixture;call directory(trim(label)//'/serial')
      call calc_ppi_modes(integration_time);call calc_gw_strain(integration_time);call calc_m_dot(integration_time)
      call record_reference(expected_packet)
    end if
    call release_grid()
    call configure(.false.);call make_layout(layout)
    call allocate_and_grid(layout);call precalculate_metric_cache();call prepare_case_boundary_data()
    p=-12345;p(1:nx,1:ny,1:nz,:)=initial_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    up=initial_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    call initialize_halo_exchange(layout);call initialize_parallel_diagnostics(layout)
    final_time=0.1d0
    do while(integration_time<final_time)
      call update_adaptive_dt();dt=clip_step_to_target(integration_time,dt,final_time)
      n_steps=n_steps+1;local_dt(n_steps)=dt;integration_time=integration_time+dt;call advance_rk_step(.true.)
    end do
    call require(n_steps==ref_steps.and.all(local_dt(1:n_steps)==ref_dt(1:n_steps)),'FD/FM exact dt sequence')
    error=max(maxval(abs(up-final_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:))), &
      maxval(abs(p(1:nx,1:ny,1:nz,:)-final_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:))))
    field_error=max(field_error,error);call require(error==0,'FD/FM exact evolved fields')
    call calc_ppi_modes(integration_time,coefficients);call calc_gw_strain(integration_time,stress)
    call calc_m_dot(integration_time,mdot);packet=[coefficients,stress,mdot];call reduce_observable_sum(packet)
    if (parallel_rank==0)call compare(packet,expected_packet,'FD/FM distributed observables')
    call collect_conservation_balance(global_budget);call collect_conservation_balance(root_budget,root_only=.true.)
    if (parallel_rank==0)then
      call require(maxval(abs(global_budget%residual-root_budget%residual)) &
        /maxval(abs(global_budget%current))<1.0d-11,'root budget equals all-rank budget')
    end if
    write(label,'(A,I0)')'evolution_',fixture;call directory(trim(label)//'/distributed')
    call write_distributed_observables(integration_time);call write_distributed_diagnostics(integration_time)
    call finalize_parallel_diagnostics();call finalize_halo_exchange();call release_grid();cases=cases+1
  end subroutine check_evolution

  subroutine invalid_call()
    logical :: destination
    call configure(.true.);call make_layout(layout)
    call allocate_and_grid(layout);call precalculate_metric_cache();call choose_probe()
    call prescribe(1,expected);call initialize_parallel_diagnostics(layout)
    select case(trim(mode))
    case('bad-flags')
      if (parallel_rank==parallel_size-1)do_gw_extraction=.false.
    case('bad-spin')
      if (parallel_rank==parallel_size-1)a_spin=0.7d0
    case('bad-probe')
      if (parallel_rank==parallel_size-1)idx_probe=-1
    case('bad-state')
      if (parallel_rank==parallel_size-1)p(1,1,1,eq_de)=ieee_value(0.0d0,ieee_quiet_nan)
    case('bad-time')
      integration_time=1
    case('bad-midstep')
      call begin_budget_step(0.001d0)
    case('bad-root-destination')
      destination=parallel_rank==0
      call collect_global_state(root_state,root_only=destination)
      call parallel_abort('D02 invalid destination accepted')
    case('bad-duplicate','bad-existing','bad-io')
      call directory('invalid_output')
      if (parallel_rank==0.and.trim(mode)=='bad-existing')then
        open(newunit=unit,file='invalid_output/d02/mpi_ppi_modes.dat',status='new');write(unit,*)'preserve';close(unit)
      end if
      if (trim(mode)=='bad-io')output_folder='missing/parent/path'
      call write_distributed_observables(integration_time)
    case default
      call parallel_abort('unknown D02 invalid mode')
    end select
    call write_distributed_observables(0.0d0)
    call parallel_abort('D02 invalid call accepted')
  end subroutine invalid_call

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y))deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z))deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    if (allocated(michel_injector))deallocate(michel_injector)
  end subroutine release_grid
end program test_parallel_observables
