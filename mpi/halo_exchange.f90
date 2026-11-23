! CPU primitive slabs along x/r. Transverse boundaries follow the receive.
module halo_exchange
  use parallel_runtime
  use domain_decomposition, only: domain_layout, validate_domain
  use variables, only: nx,ny,nz,nghost,neq,int64,case_name,accelerator_enabled, &
                       alpha_c,beta_c,michel_injector,nx_global,ny_global,nz_global,x_global_offset,y_global_offset
  use conditions, only: set_boundary_conditions
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: initialize_halo_exchange, exchange_primitive_halos, fill_primitive_halos, finalize_halo_exchange
  public :: halo_layout_ready
  public :: begin_radial_overlap, end_radial_overlap, finish_polar_overlap
  type(domain_layout) :: domain
  logical :: initialized=.false.
  character(len=len(case_name)) :: boundary_name=''
  real*8, allocatable :: send_buffer(:,:,:,:),receive_buffer(:,:,:,:)
  real*8, allocatable :: polar_send(:,:,:,:),polar_receive(:,:,:,:)
  integer :: buffer_count=0,polar_count=0
#ifdef USE_MPI
  real*8, allocatable, asynchronous :: overlap_send(:,:,:,:,:),overlap_recv(:,:,:,:,:)
  real*8, allocatable, asynchronous :: overlap_polar_send(:,:,:,:,:),overlap_polar_recv(:,:,:,:,:)
  type(MPI_Request) :: overlap_requests(4),polar_requests(4)
  logical :: polar_active=.false.
  logical :: overlap_active=.false.
#endif
#ifdef USE_MPI
  integer :: lower_peer,upper_peer,lower_polar,upper_polar
#endif

contains

  logical function halo_layout_ready(layout) result(ready)
    type(domain_layout), intent(in) :: layout
    ready=initialized.and.case_name==boundary_name.and.nghost==domain%halo.and. &
      all([nx,ny,nz]==domain%local_shape).and.all(layout%global_shape==domain%global_shape).and. &
      all(layout%local_shape==domain%local_shape).and.all(layout%offset==domain%offset).and. &
      all(layout%process_shape==domain%process_shape).and.layout%rank==domain%rank.and. &
      layout%ranks==domain%ranks
  end function halo_layout_ready

  subroutine initialize_halo_exchange(layout)
    type(domain_layout), intent(in) :: layout
    character(len=160) :: message
    logical :: supported,periodic,cache_ok,injector_ok,polar_supported
    integer :: signature(9+len(case_name)),i,allocation_status
    integer(int64) :: count64
#ifdef USE_MPI
    integer :: minimum(size(signature)),maximum(size(signature)),status
#endif

    call parallel_require_active()
    call validate_domain(layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    call parallel_require_collective(.not.initialized,'halo exchange already initialized')
    call parallel_require_collective(.not.accelerator_enabled,'halo exchange supports CPU only')
    call parallel_require_collective(layout%rank==parallel_rank.and.layout%ranks==parallel_size, &
                                     'halo layout must match MPI_COMM_WORLD')
    call parallel_require_collective(all(layout%local_shape==[nx,ny,nz]).and.layout%halo==nghost, &
                                     'halo layout must match allocated local grid')
    call parallel_require_collective(all(layout%global_shape==[nx_global,ny_global,nz_global]).and. &
                                     all(layout%offset(1:2)==[x_global_offset,y_global_offset]),'halo layout must match global grid metadata')
    supported=.true.
    select case(trim(case_name))
    case('Sod','1','Strong','2','Shu','3','KHI','4','Jet','5','Conv','6','MichelA','7','Dust','8', &
         'OffAxis','9','FontDaigneEqu','10','FontDaigneSag','11','FontDaigne3D','12', &
         'FishMoncEqu','13','FishMoncSag','14','FishMonc3D','15')
    case default
      supported=.false.
    end select
    call parallel_require_collective(supported,'unsupported boundary case for halo exchange')
    periodic=trim(case_name)=='KHI'.or.trim(case_name)=='4'.or.trim(case_name)=='Conv'.or.trim(case_name)=='6'
    call parallel_require_collective(layout%periodic_x.eqv.periodic,'halo periodicity does not match the boundary case')
    if (layout%process_shape(2)>1) then
      polar_supported=periodic
      select case(trim(case_name))
      case('FontDaigneSag','11','FishMoncSag','14','FontDaigne3D','12','FishMonc3D','15')
        polar_supported=.true.
      end select
      call parallel_require_collective(polar_supported,'polar decomposition is unsupported for this boundary case')
      call parallel_require_collective(layout%periodic_y.eqv.periodic,'polar periodicity does not match boundary case')
    end if
    signature(1:9)=[layout%global_shape,layout%halo,neq,merge(1,0,periodic), &
      merge(1,0,layout%periodic_y),layout%process_shape]
    do i=1,len(case_name)
      signature(9+i)=iachar(case_name(i:i))
    end do
#ifdef USE_MPI
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(signature,minimum,size(signature),MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('halo_exchange:initialize_halo_exchange:1',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'halo global configuration minimum')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(signature,maximum,size(signature),MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('halo_exchange:initialize_halo_exchange:2',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'halo global configuration maximum')
    call parallel_require_collective(all(minimum==maximum),'inconsistent halo configuration between ranks')
#endif
    cache_ok=.true.
    select case(trim(case_name))
    case('OffAxis','9','FontDaigneEqu','10','FontDaigneSag','11','FontDaigne3D','12', &
         'FishMoncEqu','13','FishMoncSag','14','FishMonc3D','15')
      cache_ok=allocated(alpha_c).and.allocated(beta_c)
      if (cache_ok) cache_ok=all(shape(alpha_c)==[nx,ny]).and.all(shape(beta_c)==[3,nx,ny])
    end select
    call parallel_require_collective(cache_ok,'radial physical boundaries require the local metric cache')
    injector_ok=.true.
    if ((trim(case_name)=='MichelA'.or.trim(case_name)=='7').and.layout%physical_upper(1)) then
      injector_ok=allocated(michel_injector)
      if (injector_ok) injector_ok=all(shape(michel_injector)==[nghost,ny,nz,neq])
    end if
    call parallel_require_collective(injector_ok,'outer Michel owner requires its injector')
    select case(trim(case_name))
    case('FontDaigne3D','12','FishMonc3D','15')
      call parallel_require_collective(mod(nz,2)==0,'polar boundaries require even global nz')
    end select
    domain=layout
    boundary_name=case_name
    if (parallel_size>1) then
      ! Check each multiplication before it can overflow a default MPI count.
      count64=int(nghost,int64)
      do i=1,3
        if (i==1) allocation_status=ny
        if (i==2) allocation_status=nz
        if (i==3) allocation_status=neq
        call parallel_require_collective(count64<=int(huge(buffer_count),int64)/int(allocation_status,int64), &
                                         'halo slab exceeds the MPI count range')
        count64=count64*int(allocation_status,int64)
      end do
      buffer_count=int(count64)
      allocate(send_buffer(nghost,ny,nz,neq),receive_buffer(nghost,ny,nz,neq),stat=allocation_status)
      call parallel_require_collective(allocation_status==0,'halo buffer allocation failed')
    end if
    if (layout%process_shape(2)>1) then
      count64=int(nx,int64)+2_int64*nghost
      do i=1,3
        if (i==1) allocation_status=nghost
        if (i==2) allocation_status=nz
        if (i==3) allocation_status=neq
        call parallel_require_collective(count64<=int(huge(polar_count),int64)/int(allocation_status,int64), &
          'polar halo slab exceeds the MPI count range')
        count64=count64*int(allocation_status,int64)
      end do
      polar_count=int(count64)
      allocate(polar_send(nx+2*nghost,nghost,nz,neq),polar_receive(nx+2*nghost,nghost,nz,neq),stat=allocation_status)
      call parallel_require_collective(allocation_status==0,'polar halo buffer allocation failed')
    end if
#ifdef USE_MPI
    lower_peer=MPI_PROC_NULL
    upper_peer=MPI_PROC_NULL
    lower_polar=MPI_PROC_NULL
    upper_polar=MPI_PROC_NULL
    if (domain%lower_neighbor(2)>=0) lower_polar=domain%lower_neighbor(2)
    if (domain%upper_neighbor(2)>=0) upper_polar=domain%upper_neighbor(2)
    if (domain%lower_rank>=0) lower_peer=domain%lower_rank
    if (domain%upper_rank>=0) upper_peer=domain%upper_rank
#endif
    initialized=.true.
  end subroutine initialize_halo_exchange

  subroutine exchange_primitive_halos(q_state)
    real*8, intent(inout) :: q_state(-nghost:,-nghost:,-nghost:,:)
    integer :: g
#ifdef USE_MPI
    integer :: status,received_count
    type(MPI_Status) :: receipt
#endif

    call parallel_require_active()
    if (.not.initialized) call parallel_abort('halo exchange is not initialized')
    if (any([nx,ny,nz]/=domain%local_shape).or.nghost/=domain%halo.or.case_name/=boundary_name) &
      call parallel_abort('halo grid or boundary case changed after initialization')
    if (any([nx_global,ny_global,nz_global]/=domain%global_shape).or. &
        any([x_global_offset,y_global_offset]/=domain%offset(1:2))) &
      call parallel_abort('global grid metadata changed after halo initialization')
    if (any(shape(q_state)/=[nx+2*nghost+1,ny+2*nghost+1,nz+2*nghost+1,neq])) &
      call parallel_abort('primitive halo array has the wrong shape')
    if (parallel_size==1) then
      if (domain%periodic_x) then
        do g=1,nghost
          q_state(1-g,1:ny,1:nz,:)=q_state(modulo(-g,nx)+1,1:ny,1:nz,:)
          q_state(nx+g,1:ny,1:nz,:)=q_state(modulo(g-1,nx)+1,1:ny,1:nz,:)
        end do
      end if
      return
    end if
#ifdef USE_MPI
    send_buffer=q_state(1:nghost,1:ny,1:nz,:)
    call MPI_Sendrecv(send_buffer,buffer_count,MPI_DOUBLE_PRECISION,lower_peer,201, &
                      receive_buffer,buffer_count,MPI_DOUBLE_PRECISION,upper_peer,201,MPI_COMM_WORLD,receipt,status)
    call parallel_check_mpi(status,'halo left send/right receive')
    if (upper_peer/=MPI_PROC_NULL) then
      call MPI_Get_count(receipt,MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'halo right received count')
      if (received_count/=buffer_count) call parallel_abort('incomplete right halo slab')
      q_state(nx+1:nx+nghost,1:ny,1:nz,:)=receive_buffer
    end if
    send_buffer=q_state(nx-nghost+1:nx,1:ny,1:nz,:)
    call MPI_Sendrecv(send_buffer,buffer_count,MPI_DOUBLE_PRECISION,upper_peer,202, &
                      receive_buffer,buffer_count,MPI_DOUBLE_PRECISION,lower_peer,202,MPI_COMM_WORLD,receipt,status)
    call parallel_check_mpi(status,'halo right send/left receive')
    if (lower_peer/=MPI_PROC_NULL) then
      call MPI_Get_count(receipt,MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'halo left received count')
      if (received_count/=buffer_count) call parallel_abort('incomplete left halo slab')
      q_state(1-nghost:0,1:ny,1:nz,:)=receive_buffer
    end if
#endif
  end subroutine exchange_primitive_halos

  ! Experimental radial-only exchange. Buffers stay alive until the matching wait.
  subroutine begin_radial_overlap(q_state,started)
    real*8, intent(in) :: q_state(-nghost:,-nghost:,-nghost:,:)
    logical, intent(out) :: started
#ifdef USE_MPI
    integer :: status
    call parallel_require_active()
    if (.not.initialized) call parallel_abort('halo exchange is not initialized')
    if (overlap_active.or.polar_active) call parallel_abort('overlap exchange already active')
    started=parallel_size>1.and. &
      (trim(case_name)=='Jet'.or.trim(case_name)=='5'.or. &
       trim(case_name)=='FontDaigne3D'.or.trim(case_name)=='12'.or. &
       trim(case_name)=='FishMonc3D'.or.trim(case_name)=='15')
    if (.not.started) return
    if (.not.allocated(overlap_send)) then
      allocate(overlap_send(nghost,ny,nz,neq,2),overlap_recv(nghost,ny,nz,neq,2),stat=status)
      if (status/=0) call parallel_abort('overlap buffer allocation failed')
    end if
    overlap_send(:,:,:,:,1)=q_state(1:nghost,1:ny,1:nz,:)
    overlap_send(:,:,:,:,2)=q_state(nx-nghost+1:nx,1:ny,1:nz,:)
    call MPI_Irecv(overlap_recv(:,:,:,:,1),buffer_count,MPI_DOUBLE_PRECISION,upper_peer,201, &
      MPI_COMM_WORLD,overlap_requests(1),status)
    call parallel_check_mpi(status,'overlap receive right')
    call MPI_Irecv(overlap_recv(:,:,:,:,2),buffer_count,MPI_DOUBLE_PRECISION,lower_peer,202, &
      MPI_COMM_WORLD,overlap_requests(2),status)
    call parallel_check_mpi(status,'overlap receive left')
    call MPI_Isend(overlap_send(:,:,:,:,1),buffer_count,MPI_DOUBLE_PRECISION,lower_peer,201, &
      MPI_COMM_WORLD,overlap_requests(3),status)
    call parallel_check_mpi(status,'overlap send left')
    call MPI_Isend(overlap_send(:,:,:,:,2),buffer_count,MPI_DOUBLE_PRECISION,upper_peer,202, &
      MPI_COMM_WORLD,overlap_requests(4),status)
    call parallel_check_mpi(status,'overlap send right')
    overlap_active=.true.
#else
    started=.false.
#endif
  end subroutine begin_radial_overlap

  subroutine end_radial_overlap(q_state,polar_started)
    real*8, intent(inout) :: q_state(-nghost:,-nghost:,-nghost:,:)
    logical, intent(out), optional :: polar_started
#ifdef USE_MPI
    integer :: status,received_count
    type(MPI_Status) :: receipts(4)
    if (.not.overlap_active) call parallel_abort('overlap exchange is not active')
    call MPI_Waitall(4,overlap_requests,receipts,status)
    call parallel_check_mpi(status,'overlap wait')
    if (upper_peer/=MPI_PROC_NULL) then
      call MPI_Get_count(receipts(1),MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'overlap right count')
      if (received_count/=buffer_count) call parallel_abort('incomplete right overlap halo')
      q_state(nx+1:nx+nghost,1:ny,1:nz,:)=overlap_recv(:,:,:,:,1)
    end if
    if (lower_peer/=MPI_PROC_NULL) then
      call MPI_Get_count(receipts(2),MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'overlap left count')
      if (received_count/=buffer_count) call parallel_abort('incomplete left overlap halo')
      q_state(1-nghost:0,1:ny,1:nz,:)=overlap_recv(:,:,:,:,2)
    end if
    overlap_active=.false.
    if (present(polar_started)) polar_started=.false.
    if (domain%process_shape(2)==1) then
      call set_boundary_conditions(q_state,domain%physical_lower(1),domain%physical_upper(1))
    else
      call set_boundary_conditions(q_state,domain%physical_lower(1),domain%physical_upper(1),.false.,.false.)
      if (.not.allocated(overlap_polar_send)) then
        allocate(overlap_polar_send(nx+2*nghost,nghost,nz,neq,2), &
                 overlap_polar_recv(nx+2*nghost,nghost,nz,neq,2),stat=status)
        if (status/=0) call parallel_abort('polar overlap buffer allocation failed')
      end if
      overlap_polar_send(:,:,:,:,1)=q_state(1-nghost:nx+nghost,1:nghost,1:nz,:)
      overlap_polar_send(:,:,:,:,2)=q_state(1-nghost:nx+nghost,ny-nghost+1:ny,1:nz,:)
      call MPI_Irecv(overlap_polar_recv(:,:,:,:,1),polar_count,MPI_DOUBLE_PRECISION,upper_polar,211, &
                     MPI_COMM_WORLD,polar_requests(1),status)
      call parallel_check_mpi(status,'overlap polar receive upper')
      call MPI_Irecv(overlap_polar_recv(:,:,:,:,2),polar_count,MPI_DOUBLE_PRECISION,lower_polar,212, &
                     MPI_COMM_WORLD,polar_requests(2),status)
      call parallel_check_mpi(status,'overlap polar receive lower')
      call MPI_Isend(overlap_polar_send(:,:,:,:,1),polar_count,MPI_DOUBLE_PRECISION,lower_polar,211, &
                     MPI_COMM_WORLD,polar_requests(3),status)
      call parallel_check_mpi(status,'overlap polar send lower')
      call MPI_Isend(overlap_polar_send(:,:,:,:,2),polar_count,MPI_DOUBLE_PRECISION,upper_polar,212, &
                     MPI_COMM_WORLD,polar_requests(4),status)
      call parallel_check_mpi(status,'overlap polar send upper')
      polar_active=.true.
      if (present(polar_started)) then
        polar_started=.true.
      else
        call finish_polar_overlap(q_state)
      end if
    end if
#else
    call parallel_abort('overlap requires MPI')
#endif
  end subroutine end_radial_overlap

  subroutine finish_polar_overlap(q_state)
    real*8, intent(inout) :: q_state(-nghost:,-nghost:,-nghost:,:)
#ifdef USE_MPI
    integer :: status,received_count
    type(MPI_Status) :: receipts(4)
    if (.not.polar_active) call parallel_abort('polar overlap is not active')
    call MPI_Waitall(4,polar_requests,receipts,status)
    call parallel_check_mpi(status,'polar overlap wait')
    if (upper_polar/=MPI_PROC_NULL) then
      call MPI_Get_count(receipts(1),MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'polar overlap upper count')
      if (received_count/=polar_count) call parallel_abort('incomplete upper polar overlap')
      q_state(1-nghost:nx+nghost,ny+1:ny+nghost,1:nz,:)=overlap_polar_recv(:,:,:,:,1)
    end if
    if (lower_polar/=MPI_PROC_NULL) then
      call MPI_Get_count(receipts(2),MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'polar overlap lower count')
      if (received_count/=polar_count) call parallel_abort('incomplete lower polar overlap')
      q_state(1-nghost:nx+nghost,1-nghost:0,1:nz,:)=overlap_polar_recv(:,:,:,:,2)
    end if
    polar_active=.false.
    call set_boundary_conditions(q_state,.false.,.false.,domain%physical_lower(2),domain%physical_upper(2))
#else
    call parallel_abort('polar overlap requires MPI')
#endif
  end subroutine finish_polar_overlap

  subroutine fill_primitive_halos(q_state)
    real*8, intent(inout) :: q_state(-nghost:,-nghost:,-nghost:,:)
    call exchange_primitive_halos(q_state)
    if (domain%process_shape(2)==1) then
      call set_boundary_conditions(q_state,domain%physical_lower(1),domain%physical_upper(1))
    else
      ! Radial physical ghosts travel with theta slabs, so corners need no diagonal messages.
      call set_boundary_conditions(q_state,domain%physical_lower(1),domain%physical_upper(1),.false.,.false.)
      call exchange_polar_halos(q_state)
      call set_boundary_conditions(q_state,.false.,.false.,domain%physical_lower(2),domain%physical_upper(2))
    end if
  end subroutine fill_primitive_halos

  subroutine exchange_polar_halos(q_state)
    real*8, intent(inout) :: q_state(-nghost:,-nghost:,-nghost:,:)
#ifdef USE_MPI
    integer :: status,received_count
    type(MPI_Status) :: receipt
    polar_send=q_state(1-nghost:nx+nghost,1:nghost,1:nz,:)
    call MPI_Sendrecv(polar_send,polar_count,MPI_DOUBLE_PRECISION,lower_polar,211, &
      polar_receive,polar_count,MPI_DOUBLE_PRECISION,upper_polar,211,MPI_COMM_WORLD,receipt,status)
    call parallel_check_mpi(status,'polar lower send/upper receive')
    if (upper_polar/=MPI_PROC_NULL) then
      call MPI_Get_count(receipt,MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'polar upper received count')
      if (received_count/=polar_count) call parallel_abort('incomplete upper polar halo')
      q_state(1-nghost:nx+nghost,ny+1:ny+nghost,1:nz,:)=polar_receive
    end if
    polar_send=q_state(1-nghost:nx+nghost,ny-nghost+1:ny,1:nz,:)
    call MPI_Sendrecv(polar_send,polar_count,MPI_DOUBLE_PRECISION,upper_polar,212, &
      polar_receive,polar_count,MPI_DOUBLE_PRECISION,lower_polar,212,MPI_COMM_WORLD,receipt,status)
    call parallel_check_mpi(status,'polar upper send/lower receive')
    if (lower_polar/=MPI_PROC_NULL) then
      call MPI_Get_count(receipt,MPI_DOUBLE_PRECISION,received_count,status)
      call parallel_check_mpi(status,'polar lower received count')
      if (received_count/=polar_count) call parallel_abort('incomplete lower polar halo')
      q_state(1-nghost:nx+nghost,1-nghost:0,1:nz,:)=polar_receive
    end if
#endif
  end subroutine exchange_polar_halos

  subroutine finalize_halo_exchange()
    call parallel_require_active()
#ifdef USE_MPI
    if (overlap_active.or.polar_active) call parallel_abort('cannot finalize active overlap')
    if (allocated(overlap_send)) deallocate(overlap_send,overlap_recv)
    if (allocated(overlap_polar_send)) deallocate(overlap_polar_send,overlap_polar_recv)
#endif
    if (allocated(send_buffer)) deallocate(send_buffer)
    if (allocated(receive_buffer)) deallocate(receive_buffer)
    if (allocated(polar_send)) deallocate(polar_send,polar_receive)
    initialized=.false.
    buffer_count=0
    polar_count=0
  end subroutine finalize_halo_exchange
end module halo_exchange
