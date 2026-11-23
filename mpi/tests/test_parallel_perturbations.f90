program test_parallel_perturbations
  use variables
  use parallel_runtime
  use parallel_timing, only: clip_step_to_target
  use parallel_diagnostics
  use parallel_perturbations
  use domain_decomposition
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use metrics
  use conditions
  use equations, only: prim_to_cons
  use fluxes, only: init_wavespeed_solver
  use evolution, only: update_adaptive_dt
  use halo_exchange
  use time_integrator, only: advance_rk_step
  use output, only: calc_ppi_modes,calc_gw_strain,write_distributed_diagnostics,write_distributed_observables, &
    save_checkpoint,read_checkpoint_metadata,read_checkpoint_state
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer :: grid(2)=[1,1],fixture,kind,seed_id,cases=0,checks=0,unit,summary_unit,shape_global(3)
  integer :: seed_value=3435
  real*8 :: amplitude=0.01d0,max_field=0,max_budget=0,max_force=0,max_observable=0
  character(len=64) :: arg,mode,name
  type(domain_layout) :: layout

  call parallel_initialize()
  call get_command_argument(1,arg);if(len_trim(arg)>0)read(arg,*)grid(1)
  call get_command_argument(2,arg);if(len_trim(arg)>0)read(arg,*)grid(2)
  call get_command_argument(3,mode)
  if (trim(mode)/=''.and.trim(mode)/='report')call invalid_call()
  write(name,'(A,I0,A)')'perturbation_cases_rank_',parallel_rank,'.txt'
  open(newunit=summary_unit,file=trim(name),status='new')
  call check_extreme_seeds()
  do fixture=1,2
    do kind=1,3
      do seed_id=1,3
        seed_value=3435
        if(seed_id==2)seed_value=17
        if(seed_id==3)seed_value=huge(1)
        amplitude=0.025d0;call check_case(0)
      end do
      amplitude=0.0d0;seed_value=3435;call check_case(0)
    end do
  end do
  amplitude=0.01d0;seed_value=3435
  do fixture=3,10
    do kind=1,3
      call check_case(0)
    end do
  end do
  fixture=2
  do kind=1,3
    call check_case(1);call check_case(2)
  end do
  close(summary_unit)
  write(name,'(A,I0,A)')'perturbations_rank_',parallel_rank,'.txt'
  open(newunit=unit,file=trim(name),status='new')
  write(unit,'(4(I0,1X),4(ES24.16E3,1X))')parallel_rank,parallel_size,checks,cases, &
    max_field,max_budget,max_force,max_observable
  close(unit)
  if(parallel_rank==0)print *,'PARALLEL PERTURBATION TESTS PASSED'
  call parallel_finalize()

contains
  subroutine require(valid,label)
    logical,intent(in)::valid
    character(len=*),intent(in)::label
    checks=checks+1
    if(.not.valid)then
      write(*,*)'D03 fixture/type/seed/rank',fixture,kind,seed_value,parallel_rank
      call parallel_abort('D03 test: '//label)
    end if
  end subroutine require

  subroutine configure()
    shape_global=[37,19,8];case_name='Conv';case_id=6
    geom_type='Cartesian';metric_type='Minkowski'
    r_min=0;r_max=1;y_min=0;y_max=1;z_min=0;z_max=1
    use_log_r=.false.;bh_mass=1;a_spin=0
    if(fixture>=3)then
      case_name='FontDaigne3D';case_id=12
      if(mod((fixture-3)/2,2)==1)then
        case_name='FishMonc3D';case_id=15
      end if
      if(mod(fixture-3,2)==1)a_spin=0.9d0
      if(fixture>=7)then
        shape_global(3)=1
        if(case_id==12)then
          case_name='FontDaigneSag';case_id=11
        else
          case_name='FishMoncSag';case_id=14
        end if
      end if
      geom_type='Spheroidal';metric_type='Kerr-Schild';use_log_r=.true.
      r_min=1.2d0;r_max=40;y_max=pi;z_max=2*pi
    end if
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=3
    adb_idx=4.0d0/3.0d0;g1=adb_idx/(adb_idx-1);CFL=0.3d0
    rho_floor=1.0d-10;p_floor=1.0d-13;D_floor=rho_floor;tau_floor=p_floor*(g1-1)
    rec_method_id=REC_WENO5
    if(mod(fixture,2)==0)rec_method_id=REC_MP5
    riemann_solver_id=RS_HLLE
    if(kind==3)riemann_solver_id=RS_HLLC
    tvd_limiter_id=LIM_MC;use_shock_sensor=.true.;shock_sensor_thresholds=SENSOR_DEFAULTS
    do_numerical_diagnostics=.true.;do_radial_mass_budget=.false.;do_restart=.false.
    do_gw_extraction=.true.;do_mdot_extraction=.false.;do_ppi_diagnostics=fixture>=3.and.nz>1
    apply_perturbation=.true.;perturbation_applied=.false.;perturbation_type=kind
    perturbation_amplitude=amplitude;perturbation_seed=seed_value;perturbation_mode=2
    perturbation_time=0.012345d0;final_time=0.04d0;save_interval=1;scheme_name='d03'
    if(fixture<=2)perturbation_time=0.02d0
    call set_metric_type();call init_wavespeed_solver()
  end subroutine configure

  subroutine reset_epoch()
    integration_time=0;n_steps=0
    if(fixture<=2)then
      integration_time=0.02d0;n_steps=17
    end if
    numerical_events=0_int64;numerical_delta_u=0;perturbation_applied=.false.
  end subroutine reset_epoch

  subroutine seed_full()
    integer::i,j,k
    call set_initial_conditions()
    if(fixture<=2)then
      p(1:nx,1:ny,1:nz,eq_de)=0.2d0;p(1:nx,1:ny,1:nz,eq_pr)=0.8d0
      p(1:nx,1:ny,1:nz,eq_vx)=0.05d0;p(1:nx,1:ny,1:nz,eq_vy)=0.01d0
      p(1:nx,1:ny,1:nz,eq_vz)=0.03d0
      if(fixture==1)then
        p(nx-3:nx,ny-3:ny,:,eq_de)=2
        p(nx,ny,nz,eq_de)=20
        p(1,1,1,eq_de)=1 ! Exactly the strict 5% boundary: must remain unperturbed.
      else
        do k=1,nz
          do j=1,ny
            do i=1,nx
              p(i,j,k,eq_de)=2+0.1d0*i+0.01d0*j+0.001d0*k
            end do
          end do
        end do
      end if
    else
      call require(count(p(1:nx,1:ny,1:nz,eq_de)>100*rho_floor)>0,'torus contains fluid')
    end if
    call set_boundary_conditions(p)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
        end do
      end do
    end do
  end subroutine seed_full

  subroutine legacy_apply(change)
    real*8,intent(out)::change(neq)
    integer::i,j,k
    real*8::previous(neq)
    select case(kind)
    case(PERT_PRESSURE_NOISE);call inject_pressure_noise(p,amplitude,seed_value)
    case(PERT_DENSITY_NOISE);call inject_density_noise(p,amplitude,seed_value)
    case(PERT_DENSITY_MODE);call inject_density_mode(p,2.0d0,amplitude)
    end select
    change=0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          previous=up(i,j,k,:)
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
          change=change+up(i,j,k,:)-previous
        end do
      end do
    end do
    change=change*dx*dy*dz;perturbation_applied=.true.
    call set_boundary_conditions(p)
  end subroutine legacy_apply

  subroutine observables(values)
    real*8,intent(out)::values(12)
    values=0
    if(do_ppi_diagnostics)call calc_ppi_modes(integration_time,values(1:9))
    call calc_gw_strain(integration_time,values(10:12))
  end subroutine observables

  subroutine make_directory(label)
    character(len=*),intent(in)::label
    integer::ios
    logical::ok
    output_folder=label;ok=.true.
    if(parallel_rank==0)then
      call execute_command_line('mkdir -p '//trim(label)//'/d03',exitstat=ios);ok=ios==0
    end if
    call parallel_require_collective(ok,'create D03 test directory')
  end subroutine make_directory

  subroutine check_case(checkpoint_kind)
    integer,intent(in)::checkpoint_kind
    real*8,allocatable::initial_p(:,:,:,:),initial_u(:,:,:,:),final_p(:,:,:,:),final_u(:,:,:,:),copy(:,:,:,:)
    real*8::expected_force(neq),change(neq),expected_obs(12),actual_obs(12),reference_dt(4096),local_dt(4096)
    real*8::error,scale,force_error,obs_error,budget_error,step_error
    type(conservation_balance)::balance,before_balance
    type(state_diagnostics)::state
    integer::initial_steps,reference_steps,local_steps,i,j,k,m,events,ios,unit_file,read_step
    logical::injected,already_applied,exact
    character(len=160)::message,label,checkpoint
    call configure();call allocate_and_grid();call precalculate_metric_cache();call seed_full();call reset_epoch()
    already_applied=checkpoint_kind==2
    if(checkpoint_kind>0)then
      if(already_applied)call legacy_apply(change)
      write(label,'(A,I0)')'v9_',cases+1
      call make_directory(label)
      ios=0
      if(parallel_rank==0)call save_checkpoint(n_steps,integration_time,ios)
      call parallel_require_collective(ios==0,'V9 flag fixture saved')
      write(checkpoint,'(A,A,I9.9,A)')trim(label),'/checkpoint_d03_step_',n_steps,'.rst'
      perturbation_applied=.not.already_applied
      open(newunit=unit_file,file=trim(checkpoint),form='unformatted',status='old')
      call read_checkpoint_metadata(unit_file)
      call read_checkpoint_state(unit_file,read_step,integration_time);close(unit_file)
      call require(perturbation_applied.eqv.already_applied,'V9 restores before/after perturbation flag')
      call require(read_step==n_steps,'V9 preserves event step')
    end if
    initial_p=p(1:nx,1:ny,1:nz,:);initial_u=up;initial_steps=n_steps
    expected_force=0;reference_steps=0
    if(fixture<=2)then
      if(.not.already_applied)call legacy_apply(expected_force)
    else
      do while(integration_time<final_time)
        call update_adaptive_dt();dt=clip_step_to_target(integration_time,dt,final_time)
        if(integration_time>=perturbation_time.and..not.perturbation_applied)call legacy_apply(expected_force)
        reference_steps=reference_steps+1
        call require(reference_steps<=size(reference_dt),'bounded reference trajectory')
        reference_dt(reference_steps)=dt;n_steps=n_steps+1;integration_time=integration_time+dt
        call advance_rk_step(.false.)
      end do
    end if
    final_p=p(1:nx,1:ny,1:nz,:);final_u=up;call observables(expected_obs)
    call release_grid()
    call configure()
    call build_domain(shape_global,3,parallel_rank,parallel_size,fixture<=2,layout,message,grid,fixture<=2)
    call parallel_require_collective(len_trim(message)==0,message)
    call allocate_and_grid(layout);call precalculate_metric_cache();call prepare_case_boundary_data()
    p=-12345
    p(1:nx,1:ny,1:nz,:)=initial_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    up=initial_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    call reset_epoch();perturbation_applied=already_applied
    call initialize_halo_exchange(layout);call fill_primitive_halos(p);call initialize_parallel_diagnostics(layout)
    write(label,'(A,I0)')'case_',cases+1
    call make_directory(label);call write_distributed_diagnostics(integration_time)
    events=0;local_steps=0;step_error=0
    if(fixture<=2)then
      call apply_parallel_perturbation(layout,injected)
      call require(injected.eqv.(.not.already_applied),'one event or restored no-op')
      if(injected)events=events+1
    else
      do while(integration_time<final_time)
        call update_adaptive_dt();dt=clip_step_to_target(integration_time,dt,final_time)
        call apply_parallel_perturbation(layout,injected)
        if(injected)events=events+1
        local_steps=local_steps+1
        call require(local_steps<=reference_steps,'distributed trajectory length')
        local_dt(local_steps)=dt;n_steps=n_steps+1;integration_time=integration_time+dt
        call advance_rk_step(.true.)
      end do
      call require(local_steps==reference_steps,'same number of RK steps')
      step_error=maxval(abs(local_dt(1:local_steps)-reference_dt(1:local_steps)))/maxval(reference_dt(1:reference_steps))
      call require(step_error<1.0d-11,'matching time-step sequence')
      call write_distributed_diagnostics(integration_time)
    end if
    call require(events==merge(0,1,already_applied),'exactly one perturbation per evolution')
    copy=up
    call collect_conservation_balance(before_balance)
    call apply_parallel_perturbation(layout,injected)
    call require(.not.injected.and.all(up==copy),'repeated request does not inject again')
    call collect_conservation_balance(balance)
    call require(all(balance%perturbation==before_balance%perturbation),'repeated request leaves forcing unchanged')
    scale=max(1.0d0,maxval(abs(final_u)))
    error=maxval(abs(up-final_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)))/scale
    scale=max(1.0d0,maxval(abs(final_p)))
    error=max(error,maxval(abs(p(1:nx,1:ny,1:nz,:)-final_p(x_global_offset+1:x_global_offset+nx, &
      y_global_offset+1:y_global_offset+ny,:,:)))/scale)
    exact=all(up==final_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)).and. &
      all(p(1:nx,1:ny,1:nz,:)==final_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:))
    call require(error<1.0d-11,'global-mask/global-order perturbed fields')
    scale=max(maxval(abs(balance%initial)),maxval(abs(balance%current)),tiny(1.0d0))
    budget_error=maxval(abs(balance%residual))/scale
    force_error=maxval(abs(balance%perturbation-expected_force))/scale
    call require(budget_error<1.0d-11,'conservation budget includes the perturbation')
    call require(force_error<1.0d-11,'instantaneous forcing equals full-grid change')
    call collect_global_state(state)
    call require(state%events(NUM_C2P_CALLS)==3_int64*local_steps*product(int(shape_global,int64)), &
      'physical recoveries remain independent of injection')
    call require(state%invalid_cells==0,'physical post-injection state')
    call observables(actual_obs);call reduce_observable_sum(actual_obs)
    obs_error=0
    if(parallel_rank==0)then
      do i=1,2
        j=1;k=9
        if(i==2)then
          j=10;k=12
        end if
        obs_error=max(obs_error,maxval(abs(actual_obs(j:k)-expected_obs(j:k)))/ &
          max(maxval(abs(expected_obs(j:k))),1.0d-200))
      end do
      call require(obs_error<1.0d-11,'perturbed PPI coefficients and Finn-Evans tensor')
    end if
    call write_distributed_observables(integration_time)
    write(summary_unit,'(7(I0,1X),6(ES24.16E3,1X))')cases+1,fixture,kind,local_steps,events,checkpoint_kind, &
      merge(1,0,exact),integration_time,error,step_error,budget_error,force_error,obs_error
    max_field=max(max_field,error);max_budget=max(max_budget,budget_error)
    max_force=max(max_force,force_error);max_observable=max(max_observable,obs_error)
    call finalize_parallel_diagnostics();call finalize_halo_exchange();call release_grid()
    cases=cases+1
  end subroutine check_case

  subroutine check_extreme_seeds()
    integer::size_seed,i
    integer,allocatable::seeds(:)
    real*8::expected(32),actual(32)
    call random_seed(size=size_seed);allocate(seeds(size_seed))
    seeds=[(1+104729*(i-1),i=1,size_seed)]
    call random_seed(put=seeds);call random_number(expected)
    call set_reproducible_random_seed(huge(1));call random_number(actual)
    call require(all(actual==expected),'seed construction does not overflow at huge(integer)')
    seeds(1)=huge(1)-2
    do i=2,size_seed
      seeds(i)=104729*(i-1)-1
    end do
    call random_seed(put=seeds);call random_number(expected)
    call set_reproducible_random_seed(-huge(1));call random_number(actual)
    call require(all(actual==expected),'negative extreme seed construction')
  end subroutine check_extreme_seeds

  subroutine invalid_call()
    logical::injected
    integer::ios
    character(len=160)::message
    fixture=2;kind=1;call configure()
    call build_domain(shape_global,3,parallel_rank,parallel_size,.true.,layout,message,grid,.true.)
    call parallel_require_collective(len_trim(message)==0,message)
    call allocate_and_grid(layout);call precalculate_metric_cache()
    p(1:nx,1:ny,1:nz,eq_de)=2;p(1:nx,1:ny,1:nz,eq_pr)=0.8d0
    p(1:nx,1:ny,1:nz,eq_vx:eq_vz)=0
    up(:,:,:,eq_de)=2;up(:,:,:,eq_pr)=0.8d0*(g1-1);up(:,:,:,eq_vx:eq_vz)=0
    call reset_epoch();call initialize_halo_exchange(layout);call fill_primitive_halos(p)
    call initialize_parallel_diagnostics(layout);call make_directory('invalid')
    if(trim(mode)=='bad-midstep')call begin_budget_step(0.01d0)
    if(trim(mode)=='bad-halo')call finalize_halo_exchange()
    if(trim(mode)=='bad-io'.and.parallel_rank==0)then
      open(newunit=ios,file='invalid/d03/mpi_perturbation_event.dat',status='new')
      write(ios,'(A)')'preserve previous evidence';close(ios)
    end if
    if(parallel_rank==parallel_size-1)then
      select case(trim(mode))
      case('bad-amplitude');perturbation_amplitude=ieee_value(0.0d0,ieee_quiet_nan)
      case('bad-seed');perturbation_seed=perturbation_seed+1
      case('bad-flag');perturbation_applied=.true.
      case('bad-clock');integration_time=integration_time+0.001d0
      case('bad-type');perturbation_type=99
      case('bad-layout');layout%offset(1)=layout%offset(1)+1
      case('bad-state');p(1,1,1,eq_de)=0
      case('bad-up');up(1,1,1,1)=ieee_value(0.0d0,ieee_quiet_nan)
      end select
    end if
    call apply_parallel_perturbation(layout,injected)
    call parallel_abort('invalid perturbation probe unexpectedly returned')
  end subroutine invalid_call

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if(allocated(alpha_f_y))deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if(allocated(alpha_f_z))deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    if(allocated(michel_injector))deallocate(michel_injector)
  end subroutine release_grid
end program test_parallel_perturbations
