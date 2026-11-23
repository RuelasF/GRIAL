! CPU global diagnostics and finite-volume budgets; conserved fields are densitized.
module parallel_diagnostics
  use variables
  use domain_decomposition, only: domain_layout,validate_domain
  use parallel_runtime
  use parallel_timing, only: validate_step_size
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: state_diagnostics,conservation_balance,diagnostics_active
  public :: initialize_parallel_diagnostics,finalize_parallel_diagnostics
  public :: collect_global_state,conserved_error_norms,collect_conservation_balance
  public :: begin_budget_step,reset_rhs_budget,capture_boundary_fluxes,capture_source_sum
  public :: accumulate_rhs_budget,accumulate_recovery_budget,end_budget_step
  public :: accumulate_perturbation_budget
  public :: reduce_observable_sum,validate_observable_sample
  public :: validate_checkpoint_boundary
  public :: diagnostic_history,export_diagnostic_history,valid_diagnostic_history

  type :: state_diagnostics
    real*8 :: totals(neq)=0.0d0,disk_mass=0.0d0
    real*8 :: rho_min=0.0d0,p_min=0.0d0,rho_max=0.0d0,p_max=0.0d0,max_vsq=0.0d0
    real*8 :: l1(neq)=0.0d0,l2(neq)=0.0d0,linf(neq)=0.0d0,raw_repair(neq)=0.0d0
    integer(int64) :: cells=0,atmosphere_cells=0,invalid_cells=0,events(NUM_EVENT_COUNT)=0
    integer(int64) :: unique_faces_per_rhs=0,work_faces_per_rhs=0,extra_interface_faces=0
  end type state_diagnostics

  type :: conservation_balance
    real*8 :: initial(neq)=0.0d0,current(neq)=0.0d0,flux(neq,2,3)=0.0d0
    real*8 :: source(neq)=0.0d0,repair(neq)=0.0d0,perturbation(neq)=0.0d0,residual(neq)=0.0d0
    integer(int64) :: steps=0
  end type conservation_balance

  type :: diagnostic_history
    real*8 :: sums(10*neq)=0.0d0
    integer(int64) :: steps=0_int64,extra_faces=0_int64
  end type diagnostic_history

  logical, protected :: diagnostics_active=.false.
  type(domain_layout) :: domain
  logical :: exterior_lower(3),exterior_upper(3),inside_step=.false.,rhs_ready=.false.
  character(len=len(case_name)) :: boundary_case=''
  character(len=len(geom_type)) :: geometry_name=''
  character(len=len(metric_type)) :: metric_name=''
  real*8 :: context_metadata(5)
  real*8 :: spacing(3),face_area(3),volume,budget_dt,initial_local(neq)
  real*8 :: flux_total(neq,2,3),source_total(neq),repair_total(neq)
  real*8 :: perturbation_total(neq)
  real*8 :: rhs_flux(neq,2,3),rhs_source(neq)
  integer :: next_rhs=0,next_recovery=0
  integer(int64) :: budget_steps=0,global_cells=0,unique_faces=0,work_faces=0
  integer(int64) :: history_start_steps=0_int64,history_extra_faces=0_int64
  real*8, parameter :: rhs_weights(3)=[1.0d0/6.0d0,1.0d0/6.0d0,2.0d0/3.0d0]
  real*8, parameter :: repair_weights(3)=[1.0d0/6.0d0,2.0d0/3.0d0,1.0d0]
  integer, parameter :: OP_SUM=1,OP_MIN=2,OP_MAX=3
#ifdef USE_MPI
  type(MPI_Datatype) :: count_type
#endif

contains

  subroutine validate_checkpoint_boundary()
    real*8 :: minimum(1),maximum(1)
    call parallel_require_active()
    minimum=merge(1.0d0,0.0d0,diagnostics_active);maximum=minimum
    call reduce_real(minimum,OP_MIN,'validate_checkpoint_boundary:minimum:OP_MIN');call reduce_real(maximum,OP_MAX,'validate_checkpoint_boundary:maximum:OP_MAX')
    call parallel_require_collective(all(minimum==maximum),'checkpoint diagnostic activation differs between ranks')
    if (diagnostics_active) then
      call validate_context();call agree_epoch()
      call parallel_require_collective(.not.inside_step,'cannot checkpoint inside an RK step')
    end if
  end subroutine validate_checkpoint_boundary

  subroutine reduce_real(values,operation,site,root_only)
    real*8, intent(inout) :: values(:)
    integer, intent(in) :: operation
    character(len=*), intent(in) :: site
    logical, intent(in), optional :: root_only
    logical :: to_root
    real*8 :: received(size(values))
#ifdef USE_MPI
    integer :: status
    type(MPI_Op) :: op
    to_root=.false.
    if (present(root_only)) to_root=root_only
    if (parallel_size==1) return
    select case(operation)
    case(OP_SUM);op=MPI_SUM
    case(OP_MIN);op=MPI_MIN
    case(OP_MAX);op=MPI_MAX
    end select
    if (to_root) then
#ifdef GRIAL_PROFILE
      block
        real*8 :: collective_start
        collective_start=profile_clock()
#endif
      call MPI_Reduce(values,received,size(values),MPI_DOUBLE_PRECISION,op,0,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
        call profile_label_add(site,collective_start)
      end block
#endif
      if (parallel_rank==0) values=received
    else
#ifdef GRIAL_PROFILE
      block
        real*8 :: collective_start
        collective_start=profile_clock()
#endif
      call MPI_Allreduce(MPI_IN_PLACE,values,size(values),MPI_DOUBLE_PRECISION,op,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
        call profile_label_add(site,collective_start)
      end block
#endif
    end if
    call parallel_check_mpi(status,'global diagnostic real reduction')
#endif
  end subroutine reduce_real

  subroutine reduce_counts(values,root_only)
    integer(int64), intent(inout) :: values(:)
    logical, intent(in) :: root_only
    integer(int64) :: received(size(values))
#ifdef USE_MPI
    integer :: status
#endif
    call parallel_require_collective(all(values>=0_int64).and.all(values<=huge(0_int64)/parallel_size), &
      'diagnostic counters exceed the int64 reduction range')
#ifdef USE_MPI
    if (parallel_size>1) then
      if (root_only) then
#ifdef GRIAL_PROFILE
        block
          real*8 :: collective_start
          collective_start=profile_clock()
#endif
        call MPI_Reduce(values,received,size(values),count_type,MPI_SUM,0,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
          call profile_label_add('parallel_diagnostics:reduce_counts_root:125',collective_start)
        end block
#endif
        if (parallel_rank==0) values=received
      else
#ifdef GRIAL_PROFILE
        block
          real*8 :: collective_start
          collective_start=profile_clock()
#endif
        call MPI_Allreduce(MPI_IN_PLACE,values,size(values),count_type,MPI_SUM,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
          call profile_label_add('parallel_diagnostics:reduce_counts:2',collective_start)
        end block
#endif
      end if
      call parallel_check_mpi(status,'global diagnostic int64 reduction')
    end if
#endif
  end subroutine reduce_counts

  logical function root_destination(root_only) result(to_root)
    logical, intent(in), optional :: root_only
    real*8 :: lower(1),upper(1)
    to_root=.false.
    if (present(root_only)) to_root=root_only
    lower=merge(1.0d0,0.0d0,to_root);upper=lower
    call reduce_real(lower,OP_MIN,'reduce_counts:lower:OP_MIN');call reduce_real(upper,OP_MAX,'reduce_counts:upper:OP_MAX')
    call parallel_require_collective(all(lower==upper),'diagnostic reduction destination differs between ranks')
  end function root_destination

  subroutine reduce_observable_sum(values)
    real*8, intent(inout) :: values(:)
    call parallel_require_active()
    call parallel_require_collective(all(ieee_is_finite(values)),'nonfinite local observable contribution')
    call reduce_real(values,OP_SUM,'reduce_observable_sum:values:OP_SUM',.true.)
    call parallel_require_collective(parallel_rank/=0.or.all(ieee_is_finite(values)), &
      'nonfinite global observable sum')
  end subroutine reduce_observable_sum

  subroutine validate_observable_sample(current_time)
    real*8, intent(in) :: current_time
    real*8 :: lower(9),upper(9)
    logical :: valid,owns_probe
    call validate_context();call agree_epoch()
    call parallel_require_collective(.not.inside_step,'cannot sample observables inside an RK step')
    valid=ieee_is_finite(current_time).and.ieee_is_finite(bh_mass).and.ieee_is_finite(a_spin)
    call parallel_require_collective(valid,'nonfinite observable time or metric parameters')
    call parallel_require_collective(current_time==integration_time,'observable time differs from RK epoch')
    lower=[merge(1.0d0,0.0d0,do_ppi_diagnostics),merge(1.0d0,0.0d0,do_gw_extraction), &
      merge(1.0d0,0.0d0,do_mdot_extraction),bh_mass,a_spin,merge(1.0d0,0.0d0,use_log_r), &
      dble(idx_probe_global),dble(PPI_MAX_MODE),0.0d0];upper=lower
    call reduce_real(lower,OP_MIN,'validate_observable_sample:lower:OP_MIN');call reduce_real(upper,OP_MAX,'validate_observable_sample:upper:OP_MAX')
    call parallel_require_collective(all(lower==upper),'observable controls differ between ranks')
    valid=all(ieee_is_finite(p(1:nx,1:ny,1:nz,:))).and.all(ieee_is_finite(up)).and. &
      all(ieee_is_finite(alpha_c)).and.all(ieee_is_finite(beta_c))
    call parallel_require_collective(valid,'nonfinite physical observable state or metric')
    if ((do_gw_extraction.or.do_mdot_extraction).and.trim(metric_type)/='Minkowski'.and.bh_mass>=0.1d0) then
      call parallel_require_collective((trim(geom_type)=='Spheroidal'.or.trim(geom_type)=='Spherical').and. &
        (trim(metric_type)=='Kerr-Schild'.or.trim(metric_type)=='Eddington-Finkelstein'), &
        'curved observables require spherical EF or spheroidal KS coordinates')
    end if
    if (do_mdot_extraction.and.trim(metric_type)/='Minkowski'.and.bh_mass>=0.1d0) then
      valid=idx_probe_global>=1.and.idx_probe_global<=nx_global
      owns_probe=idx_probe_global>x_global_offset.and.idx_probe_global<=x_global_offset+nx
      if (owns_probe) then
        valid=valid.and.idx_probe==idx_probe_global-x_global_offset
      else
        valid=valid.and.idx_probe==0
      end if
      call parallel_require_collective(valid,'invalid distributed accretion-probe ownership')
    end if
  end subroutine validate_observable_sample

  subroutine validate_context()
    logical :: valid
    call parallel_require_active()
    call parallel_require_collective(diagnostics_active,'parallel diagnostics are not initialized')
    valid=all([nx,ny,nz]==domain%local_shape).and.nghost==domain%halo.and.case_name==boundary_case
    valid=valid.and.all([nx_global,ny_global,nz_global]==domain%global_shape).and. &
      all([x_global_offset,y_global_offset]==domain%offset(1:2))
    valid=valid.and.geom_type==geometry_name.and.metric_type==metric_name
    if (all(ieee_is_finite([dx,dy,dz,rho_floor,p_floor]))) then
      valid=valid.and.all([dx,dy,dz,rho_floor,p_floor]==context_metadata)
    else
      valid=.false.
    end if
    valid=valid.and.do_numerical_diagnostics.and..not.accelerator_enabled
    call parallel_require_collective(valid,'parallel diagnostic grid or controls changed')
  end subroutine validate_context

  subroutine agree_epoch()
    real*8 :: lower(3),upper(3)
    logical :: valid
    valid=ieee_is_finite(integration_time)
    if (valid) valid=integration_time>=0.0d0.and.n_steps>=0
    call parallel_require_collective(valid,'invalid diagnostic time epoch')
    lower=[integration_time,dble(n_steps),dble(budget_steps)];upper=lower
    call reduce_real(lower,OP_MIN,'agree_epoch:lower:OP_MIN');call reduce_real(upper,OP_MAX,'agree_epoch:upper:OP_MAX')
    if (any(lower/=upper)) call parallel_abort('inconsistent diagnostic time epoch')
  end subroutine agree_epoch

  integer(int64) function cell_count(shape_global) result(count)
    integer, intent(in) :: shape_global(3)
    integer :: d
    count=1_int64
    do d=1,3
      if (shape_global(d)<1) then
        count=-1_int64;return
      end if
      if (count>huge(0_int64)/int(shape_global(d),int64)) then
        count=-1_int64;return
      end if
      count=count*int(shape_global(d),int64)
    end do
  end function cell_count

  subroutine initialize_parallel_diagnostics(layout,history)
    type(domain_layout), intent(in) :: layout
    type(diagnostic_history), intent(in), optional :: history
    character(len=160) :: message
    integer :: i,d,signature(8+len(case_name)+len(geom_type)+len(metric_type)),face_shape(3)
    character(len=len(case_name)+len(geom_type)+len(metric_type)) :: labels
    real*8 :: metadata(5),lower(5),upper(5)
    logical :: valid,periodic(3),supported
#ifdef USE_MPI
    integer :: status,minimum(size(signature)),maximum(size(signature))
#endif
    call parallel_require_active()
    call parallel_require_collective(.not.diagnostics_active,'parallel diagnostics already initialized')
    call validate_domain(layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    call parallel_require_collective(.not.accelerator_enabled.and.do_numerical_diagnostics, &
      'parallel diagnostics require CPU numerical counters')
    call parallel_require_collective(.not.do_restart.or.present(history),'distributed diagnostic restart requires history')
    if(present(history))call parallel_require_collective(valid_diagnostic_history(history,n_steps), &
      'invalid persistent diagnostic history')
    valid=layout%rank==parallel_rank.and.layout%ranks==parallel_size.and.layout%halo==nghost
    valid=valid.and.all(layout%local_shape==[nx,ny,nz]).and. &
      all(layout%global_shape==[nx_global,ny_global,nz_global]).and. &
      all(layout%offset(1:2)==[x_global_offset,y_global_offset])
    call parallel_require_collective(valid,'diagnostic layout does not match the allocated grid')
    valid=allocated(up).and.allocated(p).and.allocated(gamma_c)
    if (valid) valid=all(shape(up)==[nx,ny,nz,neq]).and.all(shape(gamma_c)==[3,3,nx,ny]).and. &
      all(shape(p)==[nx+2*nghost+1,ny+2*nghost+1,nz+2*nghost+1,neq])
    call parallel_require_collective(valid,'diagnostic state or metric arrays do not match the grid')
    metadata=[dx,dy,dz,rho_floor,p_floor]
    valid=all(ieee_is_finite(metadata))
    if (valid) valid=all(metadata>0.0d0)
    call parallel_require_collective(valid,'invalid diagnostic spacing or atmosphere floors')
    lower=metadata;upper=metadata
    call reduce_real(lower,OP_MIN,'initialize_parallel_diagnostics:lower:OP_MIN');call reduce_real(upper,OP_MAX,'initialize_parallel_diagnostics:upper:OP_MAX')
    if (any(lower/=upper)) call parallel_abort('inconsistent diagnostic spacing or floors')
    supported=.true.;periodic=[.false.,.true.,.true.]
    select case(trim(case_name))
    case('Conv','6','KHI','4');periodic=.true.
    case('Jet','5');periodic=.false.
    case('FontDaigneSag','11','FishMoncSag','14','FontDaigne3D','12','FishMonc3D','15')
      periodic(2)=.false.
    case('Sod','1','Strong','2','Shu','3','MichelA','7','Dust','8','OffAxis','9', &
         'FontDaigneEqu','10','FishMoncEqu','13')
    case default;supported=.false.
    end select
    call parallel_require_collective(supported,'unsupported diagnostic boundary case')
    call parallel_require_collective(layout%periodic_x.eqv.periodic(1),'diagnostic radial periodicity mismatch')
    if (layout%process_shape(2)>1) call parallel_require_collective(layout%periodic_y.eqv.periodic(2), &
      'diagnostic polar periodicity mismatch')
    signature(1:8)=[layout%global_shape,layout%process_shape,merge(1,0,periodic)]
    labels=case_name//geom_type//metric_type
    do i=1,len(labels)
      signature(8+i)=iachar(labels(i:i))
    end do
#ifdef USE_MPI
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(signature,minimum,size(signature),MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_diagnostics:initialize_parallel_diagnostics:3',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'diagnostic configuration minimum')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(signature,maximum,size(signature),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_diagnostics:initialize_parallel_diagnostics:4',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'diagnostic configuration maximum')
    call parallel_require_collective(all(minimum==maximum),'inconsistent diagnostic configuration')
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(0_int64)/8,count_type,status)
    call parallel_check_mpi(status,'diagnostic int64 MPI datatype')
#endif
    domain=layout;boundary_case=case_name;geometry_name=geom_type;metric_name=metric_type
    context_metadata=metadata;spacing=[dx,dy,dz];volume=dx*dy*dz
    face_area=[dy*dz,dx*dz,dx*dy]
    exterior_lower=layout%owns_lower.and.layout%active_axes.and..not.periodic
    exterior_upper=layout%owns_upper.and.layout%active_axes.and..not.periodic
    global_cells=cell_count(layout%global_shape)
    call parallel_require_collective(global_cells>0_int64,'diagnostic global cell count exceeds int64')
    unique_faces=0_int64;work_faces=0_int64
    do d=1,3
      if (.not.layout%active_axes(d)) cycle
      face_shape=layout%global_shape;face_shape(d)=face_shape(d)+1
      call add_face_count(unique_faces,cell_count(face_shape))
      face_shape(d)=layout%global_shape(d)
      if (d<=2) then
        call parallel_require_collective(face_shape(d)<=huge(1)-layout%process_shape(d), &
          'diagnostic face index exceeds the integer range')
        face_shape(d)=face_shape(d)+layout%process_shape(d)
      else
        face_shape(d)=face_shape(d)+1
      end if
      call add_face_count(work_faces,cell_count(face_shape))
    end do
    call parallel_require_collective(all(ieee_is_finite(up)),'nonfinite initial diagnostic conserved state')
    do i=1,neq
      initial_local(i)=sum(up(:,:,:,i))*volume
    end do
    flux_total=0.0d0;source_total=0.0d0;repair_total=0.0d0
    perturbation_total=0.0d0
    budget_steps=0_int64;inside_step=.false.;rhs_ready=.false.
    history_start_steps=0_int64;history_extra_faces=0_int64
    if(present(history))then
      initial_local=history%sums(1:neq)
      flux_total=reshape(history%sums(neq+1:7*neq),[neq,2,3])
      source_total=history%sums(7*neq+1:8*neq)
      repair_total=history%sums(8*neq+1:9*neq)
      perturbation_total=history%sums(9*neq+1:10*neq)
      budget_steps=history%steps;history_start_steps=history%steps
      history_extra_faces=history%extra_faces
    end if
    diagnostics_active=.true.
    call agree_epoch()
  end subroutine initialize_parallel_diagnostics

  logical function valid_diagnostic_history(history,step) result(valid)
    type(diagnostic_history), intent(in) :: history
    integer, intent(in) :: step
    valid=all(ieee_is_finite(history%sums)).and.history%steps>=0_int64.and. &
      history%steps<=int(step,int64).and.history%extra_faces>=0_int64
  end function valid_diagnostic_history

  integer(int64) function extra_face_work() result(work)
    integer(int64) :: since_restart,per_step
    per_step=3_int64*(work_faces-unique_faces)
    since_restart=budget_steps-history_start_steps
    call parallel_require_collective(since_restart>=0_int64,'diagnostic history step went backwards')
    call parallel_require_collective(since_restart<=(huge(0_int64)-history_extra_faces)/max(1_int64,per_step), &
      'diagnostic interface work exceeds int64')
    work=history_extra_faces+since_restart*per_step
  end function extra_face_work

  subroutine export_diagnostic_history(history)
    type(diagnostic_history), intent(out) :: history
    call validate_checkpoint_boundary()
    call parallel_require_collective(diagnostics_active,'cannot export an inactive diagnostic history')
    history%sums(1:neq)=initial_local
    history%sums(neq+1:7*neq)=reshape(flux_total,[6*neq])
    history%sums(7*neq+1:8*neq)=source_total
    history%sums(8*neq+1:9*neq)=repair_total
    history%sums(9*neq+1:10*neq)=perturbation_total
    history%steps=budget_steps;history%extra_faces=extra_face_work()
    call parallel_require_collective(valid_diagnostic_history(history,n_steps),'nonfinite diagnostic history')
  end subroutine export_diagnostic_history

  subroutine add_face_count(total,addition)
    integer(int64), intent(inout) :: total
    integer(int64), intent(in) :: addition
    call parallel_require_collective(addition>0_int64,'diagnostic face count exceeds int64')
    call parallel_require_collective(total<=huge(0_int64)-addition,'diagnostic face sum exceeds int64')
    total=total+addition
  end subroutine add_face_count

  subroutine begin_budget_step(step)
    real*8, intent(in) :: step
    call validate_context()
    call validate_step_size(step)
    call parallel_require_collective(.not.inside_step,'invalid conservation budget step')
    budget_dt=step;inside_step=.true.;next_rhs=1;next_recovery=1;rhs_ready=.false.
  end subroutine begin_budget_step

  subroutine reset_rhs_budget()
    call validate_context()
    call parallel_require_collective(inside_step.and.next_rhs==next_recovery.and.next_rhs<=3.and..not.rhs_ready, &
      'RHS budget capture out of sequence')
    rhs_flux=0.0d0;rhs_source=0.0d0;rhs_ready=.true.
  end subroutine reset_rhs_budget

  subroutine capture_boundary_fluxes(direction,left,right)
    integer, intent(in) :: direction
    real*8, intent(in) :: left(neq),right(neq)
    integer :: m
    real*8 :: area
    area=face_area(direction)
    do m=1,neq
      if (exterior_lower(direction)) then
        !$omp atomic update
        rhs_flux(m,1,direction)=rhs_flux(m,1,direction)+left(m)*area
      end if
      if (exterior_upper(direction)) then
        !$omp atomic update
        rhs_flux(m,2,direction)=rhs_flux(m,2,direction)+right(m)*area
      end if
    end do
  end subroutine capture_boundary_fluxes

  subroutine capture_source_sum(source_sum)
    real*8, intent(in) :: source_sum(neq)
    rhs_source=source_sum*volume
  end subroutine capture_source_sum

  subroutine accumulate_rhs_budget(stage)
    integer, intent(in) :: stage
    call parallel_require_collective(stage>=1.and.stage<=3,'invalid RK budget stage')
    call parallel_require_collective(inside_step.and.rhs_ready.and.stage==next_rhs.and.stage==next_recovery, &
      'RK RHS budget stage out of sequence')
    call parallel_require_collective(all(ieee_is_finite(rhs_flux)).and.all(ieee_is_finite(rhs_source)), &
      'nonfinite RK budget flux or source')
    flux_total=flux_total+budget_dt*rhs_weights(stage)*rhs_flux
    source_total=source_total+budget_dt*rhs_weights(stage)*rhs_source
    next_rhs=stage+1;rhs_ready=.false.
  end subroutine accumulate_rhs_budget

  subroutine accumulate_recovery_budget(stage,change)
    integer, intent(in) :: stage
    real*8, intent(in) :: change(neq)
    call parallel_require_collective(stage>=1.and.stage<=3,'invalid RK budget stage')
    call parallel_require_collective(inside_step.and.stage==next_recovery.and.next_rhs==stage+1, &
      'RK recovery budget stage out of sequence')
    call parallel_require_collective(all(ieee_is_finite(change)),'nonfinite RK recovery budget')
    repair_total=repair_total+repair_weights(stage)*change*volume
    next_recovery=stage+1
  end subroutine accumulate_recovery_budget

  subroutine end_budget_step()
    call parallel_require_collective(inside_step.and.next_rhs==4.and.next_recovery==4.and..not.rhs_ready, &
      'incomplete RK conservation budget')
    call parallel_require_collective(budget_steps<(huge(0_int64)/3_int64)/max(1_int64,work_faces), &
      'diagnostic work counter exceeds int64')
    budget_steps=budget_steps+1_int64;inside_step=.false.
  end subroutine end_budget_step

  subroutine accumulate_perturbation_budget(change)
    real*8, intent(in) :: change(neq)
    call validate_context();call agree_epoch()
    call parallel_require_collective(.not.inside_step,'cannot inject a perturbation inside an RK step')
    call parallel_require_collective(all(ieee_is_finite(change)),'nonfinite perturbation budget')
    perturbation_total=perturbation_total+change*volume
    call parallel_require_collective(all(ieee_is_finite(perturbation_total)),'nonfinite accumulated perturbation budget')
  end subroutine accumulate_perturbation_budget

  subroutine collect_conservation_balance(balance,root_only)
    type(conservation_balance), intent(out) :: balance
    logical, intent(in), optional :: root_only
    logical :: to_root
    real*8 :: packet(11*neq)
    integer :: m
    call validate_context();call agree_epoch()
    to_root=root_destination(root_only)
    call parallel_require_collective(.not.inside_step,'cannot collect a partial RK budget')
    call parallel_require_collective(all(ieee_is_finite(up)),'nonfinite diagnostic conserved state')
    packet(1:neq)=initial_local
    do m=1,neq
      packet(neq+m)=sum(up(:,:,:,m))*volume
    end do
    packet(2*neq+1:8*neq)=reshape(flux_total,[6*neq])
    packet(8*neq+1:9*neq)=source_total
    packet(9*neq+1:10*neq)=repair_total
    packet(10*neq+1:11*neq)=perturbation_total
    call reduce_real(packet,OP_SUM,'collect_conservation_balance:packet:OP_SUM',to_root)
    if (to_root.and.parallel_rank/=0) return
    balance%initial=packet(1:neq);balance%current=packet(neq+1:2*neq)
    balance%flux=reshape(packet(2*neq+1:8*neq),[neq,2,3])
    balance%source=packet(8*neq+1:9*neq);balance%repair=packet(9*neq+1:10*neq)
    balance%perturbation=packet(10*neq+1:11*neq)
    balance%residual=balance%current-balance%initial- &
      sum(balance%flux(:,1,:)-balance%flux(:,2,:),dim=2)-balance%source-balance%repair-balance%perturbation
    balance%steps=budget_steps
  end subroutine collect_conservation_balance

  subroutine collect_global_state(state,root_only)
    type(state_diagnostics), intent(out) :: state
    logical, intent(in), optional :: root_only
    logical :: to_root
    real*8 :: sums(3*neq+1),minimum(2),maximum(3+neq),squares(neq),v2,rho,pres
    integer(int64) :: counts(NUM_EVENT_COUNT+3),extra_work
    integer :: i,j,k,a,b,m
    call validate_context();call agree_epoch()
    to_root=root_destination(root_only)
    call parallel_require_collective(.not.inside_step,'cannot collect state inside an RK step')
    call parallel_require_collective(all(ieee_is_finite(up)).and.all(ieee_is_finite(p(1:nx,1:ny,1:nz,:))).and. &
      all(ieee_is_finite(gamma_c)).and.all(ieee_is_finite(numerical_delta_u)), &
      'nonfinite diagnostic state or metric')
    sums=0.0d0;minimum=huge(1.0d0);maximum=0.0d0;counts=0_int64
    maximum(1:2)=-huge(1.0d0)
    counts(1:NUM_EVENT_COUNT)=numerical_events
    counts(NUM_EVENT_COUNT+1)=cell_count([nx,ny,nz])
    do m=1,neq
      sums(m)=sum(up(:,:,:,m))*volume
      sums(neq+m)=sum(abs(up(:,:,:,m)))
      sums(2*neq+m)=numerical_delta_u(m)*volume
      maximum(3+m)=maxval(abs(up(:,:,:,m)))
    end do
    do k=1,nz
      do j=1,ny
        do i=1,nx
          rho=p(i,j,k,eq_de);pres=p(i,j,k,eq_pr);v2=0.0d0
          do a=1,3
            do b=1,3
              v2=v2+gamma_c(a,b,i,j)*p(i,j,k,eq_vx+a-1)*p(i,j,k,eq_vx+b-1)
            end do
          end do
          minimum=min(minimum,[rho,pres]);maximum(1:3)=max(maximum(1:3),[rho,pres,v2])
          if (rho>1.0d3*rho_floor) then
            sums(3*neq+1)=sums(3*neq+1)+up(i,j,k,eq_de)*volume
          else
            counts(NUM_EVENT_COUNT+2)=counts(NUM_EVENT_COUNT+2)+1_int64
          end if
          if (rho<=0.0d0.or.pres<=0.0d0.or.v2<0.0d0.or.v2>=1.0d0) &
            counts(NUM_EVENT_COUNT+3)=counts(NUM_EVENT_COUNT+3)+1_int64
        end do
      end do
    end do
    call reduce_real(sums,OP_SUM,'collect_global_state:sums:OP_SUM',to_root);call reduce_real(minimum,OP_MIN,'collect_global_state:minimum:OP_MIN',to_root)
    ! Every rank needs this scale before forming its local L2 sum.
    call reduce_real(maximum,OP_MAX,'collect_global_state:maximum:OP_MAX')
    call reduce_counts(counts,to_root)
    call parallel_require_collective((to_root.and.parallel_rank/=0).or.counts(NUM_EVENT_COUNT+1)==global_cells, &
      'incomplete diagnostic cell coverage')
    squares=0.0d0
    do m=1,neq
      if (maximum(3+m)>0.0d0) squares(m)=sum((up(:,:,:,m)/maximum(3+m))**2)
    end do
    call reduce_real(squares,OP_SUM,'collect_global_state:squares:OP_SUM',to_root)
    extra_work=extra_face_work()
    if (to_root.and.parallel_rank/=0) return
    state%totals=sums(1:neq);state%l1=sums(neq+1:2*neq)/dble(global_cells)
    state%raw_repair=sums(2*neq+1:3*neq);state%disk_mass=sums(3*neq+1)
    state%rho_min=minimum(1);state%p_min=minimum(2)
    state%rho_max=maximum(1);state%p_max=maximum(2);state%max_vsq=maximum(3)
    state%linf=maximum(4:3+neq);state%l2=state%linf*sqrt(squares/dble(global_cells))
    state%cells=counts(NUM_EVENT_COUNT+1);state%atmosphere_cells=counts(NUM_EVENT_COUNT+2)
    state%invalid_cells=counts(NUM_EVENT_COUNT+3);state%events=counts(1:NUM_EVENT_COUNT)
    state%unique_faces_per_rhs=unique_faces;state%work_faces_per_rhs=work_faces
    state%extra_interface_faces=extra_work
  end subroutine collect_global_state

  subroutine conserved_error_norms(reference,l1,l2,linf)
    real*8, intent(in) :: reference(:,:,:,:)
    real*8, intent(out) :: l1(neq),l2(neq),linf(neq)
    real*8 :: sums(2*neq)
    integer :: m
    call validate_context();call agree_epoch()
    call parallel_require_collective(.not.inside_step,'cannot collect norms inside an RK step')
    call parallel_require_collective(all(shape(reference)==[nx,ny,nz,neq]),'diagnostic reference has wrong shape')
    call parallel_require_collective(all(ieee_is_finite(reference)).and.all(ieee_is_finite(up)), &
      'nonfinite diagnostic reference or conserved state')
    do m=1,neq
      linf(m)=maxval(abs(up(:,:,:,m)-reference(:,:,:,m)))
    end do
    call reduce_real(linf,OP_MAX,'conserved_error_norms:linf:OP_MAX')
    sums=0.0d0
    do m=1,neq
      sums(m)=sum(abs(up(:,:,:,m)-reference(:,:,:,m)))
      if (linf(m)>0.0d0) sums(neq+m)=sum(((up(:,:,:,m)-reference(:,:,:,m))/linf(m))**2)
    end do
    call reduce_real(sums,OP_SUM,'conserved_error_norms:sums:OP_SUM')
    l1=sums(1:neq)/dble(global_cells);l2=linf*sqrt(sums(neq+1:2*neq)/dble(global_cells))
  end subroutine conserved_error_norms

  subroutine finalize_parallel_diagnostics()
    call parallel_require_active()
    call parallel_require_collective(diagnostics_active.and..not.inside_step,'cannot finalize an inactive or partial budget')
    diagnostics_active=.false.;boundary_case='';budget_steps=0_int64
  end subroutine finalize_parallel_diagnostics
end module parallel_diagnostics
