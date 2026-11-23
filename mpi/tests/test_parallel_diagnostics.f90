program test_parallel_diagnostics
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
  use output, only: write_distributed_diagnostics
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer :: grid(2)=[1,1],fixture=1,config=1,checks=0,cases=0,total_steps=0,unit
  character(len=64) :: mode,arg,name
  real*8 :: field_error=0.0d0,balance_error=0.0d0,state_error=0.0d0

  call parallel_initialize()
  call get_command_argument(1,arg)
  if (len_trim(arg)>0) read(arg,*)grid(1)
  call get_command_argument(2,arg)
  if (len_trim(arg)>0) read(arg,*)grid(2)
  call get_command_argument(3,mode)
  if (trim(mode)/=''.and.trim(mode)/='report') call invalid_call()
  call check_synthetic()
  call check_exterior_masks()
  do fixture=1,8
    do config=1,4
      call check_evolution()
    end do
  end do
  fixture=3
  do config=5,6
    call check_evolution()
  end do
  if (trim(mode)=='report') then
    write(name,'(A,I0,A)')'diagnostics_rank_',parallel_rank,'.txt'
    open(newunit=unit,file=trim(name),status='new')
    write(unit,'(5(I0,1X),3(ES24.16E3,1X))')parallel_rank,parallel_size,checks,cases,total_steps, &
      field_error,balance_error,state_error
    close(unit)
  end if
  if (parallel_rank==0) write(*,'(A)')'PARALLEL DIAGNOSTICS TESTS PASSED'
  call parallel_finalize()

contains

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if (.not.ok) then
      write(*,*)'D01 fixture/config/rank',fixture,config,parallel_rank
      call parallel_abort('D01 test: '//label)
    end if
  end subroutine require

  subroutine configure(shape_global,periodic)
    integer, intent(out) :: shape_global(3)
    logical, intent(out) :: periodic
    shape_global=[37,19,8]
    geom_type='Cartesian';metric_type='Minkowski';case_name='Conv'
    r_min=0.0d0;r_max=1.0d0;y_min=0.0d0;y_max=1.0d0;z_min=0.0d0;z_max=1.0d0
    a_spin=0.0d0;bh_mass=1.0d0;use_log_r=.false.;periodic=.true.;final_time=0.01d0
    select case(fixture)
    case(2,3)
      shape_global=[73,1,1];periodic=.false.;final_time=0.4d0;case_name='Sod'
      if (fixture==3)case_name='Strong'
      if (fixture==3.and.config>=5)shape_global(1)=400
    case(4)
      shape_global=[73,1,1];periodic=.false.;final_time=0.2d0;case_name='MichelA'
      geom_type='Spherical';metric_type='Eddington-Finkelstein'
      r_min=1.5d0;r_max=40.0d0;use_log_r=.true.
    case(5,6)
      periodic=.false.;final_time=0.03d0;case_name='FontDaigne3D'
      if (fixture==6) then
        case_name='FishMonc3D';a_spin=0.9d0
      end if
      geom_type='Spheroidal';metric_type='Kerr-Schild'
      r_min=1.2d0;r_max=40.0d0;use_log_r=.true.;y_max=pi;z_max=2.0d0*pi
    case(7)
      final_time=0.0001d0
    case(8)
      shape_global=[37,1,19];periodic=.false.;case_name='Jet';geom_type='Cylindrical'
      final_time=0.01d0
    end select
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=3
    adb_idx=4.0d0/3.0d0
    if (fixture==2) adb_idx=1.4d0
    if (fixture==3) adb_idx=5.0d0/3.0d0
    g1=adb_idx/(adb_idx-1.0d0);CFL=0.3d0
    rho_floor=1.0d-10;p_floor=1.0d-13
    if (fixture==7) then
      rho_floor=1.0d-6;p_floor=1.0d-8
    end if
    D_floor=rho_floor;tau_floor=p_floor*(g1-1.0d0)
    do_mdot_extraction=.false.;do_gw_extraction=.false.;do_ppi_diagnostics=.false.
    do_radial_mass_budget=.false.;do_numerical_diagnostics=.true.;do_restart=.false.
    shock_sensor_thresholds=SENSOR_DEFAULTS;use_shock_sensor=.true.
    rec_method_id=REC_MP5
    if (mod(config,2)==0)rec_method_id=REC_WENO5
    riemann_solver_id=RS_HLLE
    if (config>2)riemann_solver_id=RS_HLLC
    if (config>=5) then
      rec_method_id=REC_GODUNOV;riemann_solver_id=RS_HLLE
      if (config==6)riemann_solver_id=RS_HLLC
    end if
    tvd_limiter_id=LIM_MC
    call set_metric_type();call init_wavespeed_solver()
  end subroutine configure

  subroutine make_layout(shape_global,periodic,layout)
    integer, intent(in) :: shape_global(3)
    logical, intent(in) :: periodic
    type(domain_layout), intent(out) :: layout
    integer :: requested(2)
    character(len=160) :: message
    requested=grid
    if (shape_global(2)==1)requested=[parallel_size,1]
    call build_domain(shape_global,3,parallel_rank,parallel_size,periodic,layout,message,requested,periodic)
    call parallel_require_collective(len_trim(message)==0,message)
  end subroutine make_layout

  subroutine seed_problem()
    integer :: i,j,k
    call set_initial_conditions()
    if (fixture==5.or.fixture==6)call require(count(p(1:nx,1:ny,1:nz,eq_de)>100*rho_floor)>0,'nonempty torus')
    if (fixture==1.or.fixture==7) then
      do k=1,nz
        do j=1,ny
          do i=1,nx
            p(i,j,k,eq_de)=1.0d0+0.1d0*sin(2*pi*(x(i)+y(j)+z(k)))
            p(i,j,k,eq_pr)=1.0d0;p(i,j,k,eq_vx:eq_vz)=[0.2d0,0.1d0,-0.05d0]
            if (fixture==7) then
              p(i,j,k,eq_de)=rho_floor;p(i,j,k,eq_pr)=p_floor;p(i,j,k,eq_vx:eq_vz)=0.0d0
            end if
          end do
        end do
      end do
    end if
    call set_boundary_conditions(p)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          call prim_to_cons(p(i,j,k,:),up(i,j,k,:),i,j,k)
          if (fixture==7) then
            ! Deliberate inconsistent state exercises all three recovery corrections.
            up(i,j,k,:)=0.5d0*up(i,j,k,:)
            up(i,j,k,eq_vx:eq_vz)=[0.1d0,-0.2d0,0.3d0]*rho_floor
          end if
        end do
      end do
    end do
  end subroutine seed_problem

  subroutine reset_clock()
    integration_time=0.0d0;n_steps=0;numerical_events=0_int64;numerical_delta_u=0.0d0
  end subroutine reset_clock

  subroutine evolve(distributed,history,steps)
    logical, intent(in) :: distributed
    real*8, intent(out) :: history(:)
    integer, intent(out) :: steps
    do while(integration_time<final_time)
      call update_adaptive_dt();dt=clip_step_to_target(integration_time,dt,final_time)
      n_steps=n_steps+1
      call require(n_steps<=size(history),'bounded trajectory')
      history(n_steps)=dt;integration_time=integration_time+dt
      call advance_rk_step(distributed)
    end do
    steps=n_steps
  end subroutine evolve

  subroutine check_evolution()
    type(domain_layout) :: layout
    type(state_diagnostics) :: state
    type(conservation_balance) :: budget
    integer :: shape_global(3),reference_steps,local_steps,m,ios
    integer(int64) :: expected_calls
    logical :: periodic,io_ok
    real*8, allocatable :: initial_p(:,:,:,:),initial_u(:,:,:,:),final_p(:,:,:,:),final_u(:,:,:,:)
    real*8 :: reference_dt(4096),local_dt(4096),expected(neq),volume_global,scale,error,raw_repair(neq)
    real*8 :: l1(neq),l2(neq),linf(neq)
    call configure(shape_global,periodic)
    call allocate_and_grid();call precalculate_metric_cache();call seed_problem();call reset_clock()
    initial_p=p(1:nx,1:ny,1:nz,:);initial_u=up
    volume_global=dx*dy*dz
    call evolve(.false.,reference_dt,reference_steps)
    final_p=p(1:nx,1:ny,1:nz,:);final_u=up;raw_repair=numerical_delta_u*volume_global
    do m=1,neq
      expected(m)=sum(up(:,:,:,m))*volume_global
    end do
    call release_grid()
    call make_layout(shape_global,periodic,layout)
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
    call allocate_and_grid(layout);call precalculate_metric_cache();call prepare_case_boundary_data()
    p=-12345.0d0
    p(1:nx,1:ny,1:nz,:)=initial_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    up=initial_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:)
    call reset_clock();call initialize_halo_exchange(layout);call initialize_parallel_diagnostics(layout)
    if (trim(mode)=='report') then
      write(output_folder,'(A,I0,A,I0)')'case_',fixture,'_',config
      scheme_name='d01';io_ok=.true.
      if (parallel_rank==0) then
        call execute_command_line('mkdir -p '//trim(output_folder)//'/d01',exitstat=ios)
        io_ok=ios==0
      end if
      call parallel_require_collective(io_ok,'create D01 output directory')
      call write_distributed_diagnostics(integration_time)
    end if
    call evolve(.true.,local_dt,local_steps)
    call require(local_steps==reference_steps,'identical number of RK steps')
    call require(all(local_dt(1:local_steps)==reference_dt(1:reference_steps)),'identical dt sequence')
    error=max(maxval(abs(up-final_u(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:))), &
      maxval(abs(p(1:nx,1:ny,1:nz,:)-final_p(x_global_offset+1:x_global_offset+nx,y_global_offset+1:y_global_offset+ny,:,:))))
    field_error=max(field_error,error);call require(error==0.0d0,'diagnostics preserve evolved fields')
    call collect_global_state(state);call collect_conservation_balance(budget)
    scale=max(maxval(abs(budget%initial)),maxval(abs(budget%current)),tiny(1.0d0))
    error=maxval(abs(budget%residual))/scale
    balance_error=max(balance_error,error)
    if (error>=1.0d-11)write(*,*)'BUDGET residual/source/repair',budget%residual,budget%source,budget%repair
    call require(error<1.0d-11,'five-variable boundary/source/recovery balance')
    error=maxval(abs(state%totals-expected))/scale
    state_error=max(state_error,error);call require(error<1.0d-11,'global conserved integrals')
    call require(maxval(abs(state%raw_repair-raw_repair))/scale<1.0d-11,'raw recovery sum')
    expected_calls=3_int64*local_steps*int(shape_global(1),int64)*shape_global(2)*shape_global(3)
    call require(state%events(NUM_C2P_CALLS)==expected_calls,'global physical recovery calls')
    call require(sum(state%events(1:5))==3_int64*local_steps*state%work_faces_per_rhs,'face evaluation work')
    call require(state%extra_interface_faces==3_int64*local_steps* &
      (state%work_faces_per_rhs-state%unique_faces_per_rhs),'additional MPI face positions')
    call conserved_error_norms(final_u(x_global_offset+1:x_global_offset+nx, &
      y_global_offset+1:y_global_offset+ny,:,:),l1,l2,linf)
    call require(all(l1==0.0d0).and.all(l2==0.0d0).and.all(linf==0.0d0),'zero distributed error norms')
    if (fixture==1) then
      call require(all(budget%flux==0.0d0).and.all(budget%source==0.0d0),'periodic box without exterior flux/source')
      call require(all(budget%repair==0.0d0),'smooth box without repairs')
    end if
    if (fixture==7) then
      call require(maxval(abs(budget%repair-(budget%current-budget%initial)))/scale<1.0d-11, &
        'analytic forced recovery budget')
      call require(maxval(abs(state%raw_repair-budget%repair))>0.1d0*scale,'raw and RK-weighted recovery differ')
      call require(state%events(NUM_C2P_LOW_D)==expected_calls,'forced low-density recovery at every stage')
    end if
    if (trim(mode)=='report') call write_distributed_diagnostics(integration_time)
    call finalize_parallel_diagnostics();call finalize_halo_exchange();call release_grid()
    cases=cases+1;total_steps=total_steps+local_steps
  end subroutine check_evolution

  subroutine check_synthetic()
    type(domain_layout) :: layout
    type(state_diagnostics) :: state
    type(conservation_balance) :: budget
    integer :: shape_global(3),i,j,k,m,gi,gj
    integer(int64) :: expected_large
    logical :: periodic
    real*8 :: expected(neq),value(neq),l1(neq),l2(neq),linf(neq),error
    real*8, allocatable :: reference(:,:,:,:)
    call configure(shape_global,periodic);call make_layout(shape_global,periodic,layout)
    call allocate_and_grid(layout);call precalculate_metric_cache();call reset_clock()
    p=ieee_value(0.0d0,ieee_quiet_nan)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          gi=i+x_global_offset;gj=j+y_global_offset
          up(i,j,k,:)=[1.0d0+0.125d0*gi,2.0d0+0.25d0*gj,0.5d0*k,-0.25d0*gi,0.125d0*gj]
          p(i,j,k,:)=[1.0d0,2.0d0,0.1d0,0.2d0,-0.1d0]
        end do
      end do
    end do
    p(1,1,1,eq_pr)=-1.0d0
    numerical_events(NUM_HLLC_CALLS)=2_int64**31+parallel_rank
    call initialize_parallel_diagnostics(layout);call collect_global_state(state)
    expected=0.0d0
    do k=1,shape_global(3)
      do j=1,shape_global(2)
        do i=1,shape_global(1)
          value=[1.0d0+0.125d0*i,2.0d0+0.25d0*j,0.5d0*k,-0.25d0*i,0.125d0*j]
          expected=expected+value*dx*dy*dz
        end do
      end do
    end do
    error=maxval(abs(state%totals-expected))/maxval(abs(expected))
    state_error=max(state_error,error);call require(error<1.0d-11,'synthetic global mass and momenta')
    expected_large=int(parallel_size,int64)*2_int64**31+int(parallel_size,int64)*(parallel_size-1)/2
    call require(state%events(NUM_HLLC_CALLS)==expected_large,'MPI counters retain int64 values')
    call require(state%cells==37_int64*19*8.and.state%invalid_cells==parallel_size,'physical cells counted once')
    call require(state%p_min==-1.0d0.and.state%p_max==2.0d0,'global pressure extrema')
    call require(abs(state%max_vsq-0.06d0)<1.0d-14,'metric velocity extremum')
    call require(state%unique_faces_per_rhs==18023_int64,'unique global face positions')
    reference=up+1.0d0
    call conserved_error_norms(reference,l1,l2,linf)
    call require(all(l1==1.0d0).and.all(l2==1.0d0).and.all(linf==1.0d0),'analytic global L1 L2 Linf errors')
    call collect_conservation_balance(budget)
    call require(all(budget%residual==0.0d0),'initial balance is zero')
    call finalize_parallel_diagnostics();call release_grid()
  end subroutine check_synthetic

  subroutine check_exterior_masks()
    type(domain_layout) :: layout
    type(conservation_balance) :: budget
    integer :: shape_global(3),variant,stage,i,j,k,m,d
    logical :: periodic,periodic_axes(3)
    character(len=160) :: message
    real*8 :: left(neq),right(neq),zero(neq),expected_flux(neq,2,3)
    fixture=1;config=1
    do variant=1,3
      call configure(shape_global,periodic)
      periodic_axes=.true.
      if (variant==2) then
        case_name='Sod';periodic_axes(1)=.false.
      end if
      if (variant==3) then
        case_name='Jet';periodic_axes=.false.
      end if
      call build_domain(shape_global,3,parallel_rank,parallel_size,periodic_axes(1), &
        layout,message,grid,periodic_axes(2))
      call parallel_require_collective(len_trim(message)==0,message)
      call allocate_and_grid(layout);call precalculate_metric_cache();call reset_clock()
      up=1.0d0;p=1.0d0
      p(:,:,:,eq_vx:eq_vz)=0.0d0
      call initialize_parallel_diagnostics(layout)
      left=[(dble(m),m=1,neq)];right=2.0d0*left;zero=0.0d0
      call begin_budget_step(1.0d0)
      do stage=1,3
        call reset_rhs_budget()
        do k=1,nz
          do j=1,ny
            call capture_boundary_fluxes(DIR_X,left,right)
          end do
          do i=1,nx
            call capture_boundary_fluxes(DIR_Y,left,right)
          end do
        end do
        do j=1,ny
          do i=1,nx
            call capture_boundary_fluxes(DIR_Z,left,right)
          end do
        end do
        call capture_source_sum(zero);call accumulate_rhs_budget(stage)
        call accumulate_recovery_budget(stage,zero)
      end do
      call end_budget_step();call collect_conservation_balance(budget)
      expected_flux=0.0d0
      do d=1,3
        if (periodic_axes(d))cycle
        expected_flux(:,1,d)=left;expected_flux(:,2,d)=right
      end do
      call require(maxval(abs(budget%flux-expected_flux))<1.0d-11,'physical exterior masks in three axes')
      call require(maxval(abs(budget%residual-count(.not.periodic_axes)*left))<1.0d-11, &
        'outward flux sign; MPI interfaces excluded')
      call finalize_parallel_diagnostics();call release_grid()
    end do
  end subroutine check_exterior_masks

  subroutine invalid_call()
    type(domain_layout) :: layout
    type(state_diagnostics) :: state
    type(conservation_balance) :: budget
    integer :: shape_global(3),ios
    logical :: periodic
    real*8 :: l1(neq),l2(neq),linf(neq)
    real*8, allocatable :: reference(:,:,:,:)
    call configure(shape_global,periodic);call make_layout(shape_global,periodic,layout)
    call allocate_and_grid(layout);call precalculate_metric_cache();call seed_problem();call reset_clock()
    if (parallel_rank==parallel_size-1) then
      if (trim(mode)=='bad-layout')y_global_offset=y_global_offset+1
      if (trim(mode)=='bad-counters-disabled')do_numerical_diagnostics=.false.
      if (trim(mode)=='bad-restart')do_restart=.true.
      if (trim(mode)=='bad-geometry')metric_type='Kerr-Schild'
      if (trim(mode)=='bad-floors')rho_floor=2*rho_floor
    end if
    if (trim(mode)=='bad-inactive') call collect_global_state(state)
    call initialize_parallel_diagnostics(layout)
    if (trim(mode)=='bad-stage') then
      call begin_budget_step(0.001d0);call reset_rhs_budget();call accumulate_rhs_budget(2)
    end if
    if (trim(mode)=='bad-midstep') then
      call begin_budget_step(0.001d0);call collect_conservation_balance(budget)
    end if
    if (trim(mode)=='bad-step') then
      dt=0.001d0
      if (parallel_rank==parallel_size-1)dt=0.002d0
      call begin_budget_step(dt)
    end if
    if (trim(mode)=='bad-shape') then
      reference=up
      if (parallel_rank==parallel_size-1) then
        deallocate(reference);allocate(reference(nx-1,ny,nz,neq));reference=0.0d0
      end if
      call conserved_error_norms(reference,l1,l2,linf)
    end if
    if (trim(mode)=='bad-duplicate') then
      output_folder='duplicate_output';scheme_name='d01'
      if (parallel_rank==0)call execute_command_line('mkdir -p duplicate_output/d01',exitstat=ios)
      call write_distributed_diagnostics(integration_time)
      call write_distributed_diagnostics(integration_time)
    end if
    if (parallel_rank==parallel_size-1) then
      select case(trim(mode))
      case('bad-epoch');integration_time=0.001d0
      case('bad-grid');y_global_offset=y_global_offset+1
      case('bad-state');up(1,1,1,1)=ieee_value(0.0d0,ieee_quiet_nan)
      case('bad-count-range');numerical_events(1)=huge(0_int64)
      end select
    end if
    call collect_global_state(state)
    call parallel_abort('invalid diagnostic call was not rejected')
  end subroutine invalid_call

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y)) deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z)) deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    if (allocated(michel_injector))deallocate(michel_injector)
  end subroutine release_grid
end program test_parallel_diagnostics
