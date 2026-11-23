! CPU checkpoint generations: immutable blocks, a manifest, then one commit marker.
module parallel_checkpoint
  use variables
  use parallel_runtime
  use domain_decomposition, only: domain_layout,validate_domain,build_domain
  use parallel_diagnostics, only: diagnostics_active,validate_checkpoint_boundary
  use halo_exchange, only: halo_layout_ready,fill_primitive_halos
  use output, only: write_checkpoint_metadata,read_checkpoint_metadata
  use parameters, only: validate_parameters
  use initialization, only: choose_numerical_architecture
  use metrics, only: set_metric_type
  use equations, only: cons_to_prim
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: save_parallel_checkpoint,load_parallel_checkpoint,checkpoint_crc64
  public :: read_parallel_checkpoint_configuration,read_serial_checkpoint_configuration,import_serial_checkpoint
  public :: checkpoint_file_signature
  character(len=16), parameter :: BLOCK_MAGIC='GRIAL_CP_BLOCK1',MANIFEST_MAGIC='GRIAL_MPI_CP_V1'
  character(len=8), parameter :: POLAR_CONFIG_MAGIC='GRPOLAR1'
  character(len=16), parameter :: COMMIT_MAGIC='GRIAL_CP_COMMIT1',END_MAGIC='GRIAL_CP_END1'
  integer, parameter :: ABI(5)=[storage_size(0),storage_size(0_int64),storage_size(1.0d0), &
    storage_size(.true.),int(z'01020304')]
  type :: checkpoint_manifest
    integer :: abi_tag(5)=0,global_shape(3)=0,process_shape(2)=0,halo=0,ranks=0,equations=0,step=0
    logical :: periodic(2)=.false.,applied=.false.
    real*8 :: time=0.0d0,step_size=0.0d0
    integer(int64) :: configuration(2)=0_int64
    integer(int64) :: signature(2)=0_int64
    integer(int64), allocatable :: blocks(:,:)
  end type checkpoint_manifest
  integer(int64), save :: crc_table(0:255)
  logical, save :: table_ready=.false.
  type :: serial_checkpoint_index
    integer :: version=0,shape(3)=0,step=0
    integer(int64) :: payload=0_int64,signature(2)=0_int64,events(NUM_EVENT_COUNT)=0_int64
    real*8 :: time=0.0d0,delta(neq)=0.0d0
    logical :: applied=.false.
  end type serial_checkpoint_index
  type(serial_checkpoint_index), save :: serial_source
  character(len=512), save :: serial_source_path=''
  character, allocatable, save :: serial_configuration(:)
#ifdef USE_MPI
  type(MPI_Datatype) :: count_type
#endif
  interface
    integer(c_int) function posix_mkdir(path,mode) bind(C,name='mkdir')
      import c_char,c_int
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value :: mode
    end function posix_mkdir
  end interface
contains

  subroutine checkpoint_file_signature(path,signature,status)
    character(len=*), intent(in) :: path
    integer(int64), intent(out) :: signature(2)
    integer, intent(out) :: status
    call file_signature(path,signature,status)
  end subroutine checkpoint_file_signature

  subroutine prepare_crc()
    integer :: i,bit
    integer(int64) :: value
    if (table_ready) return
    do i=0,255
      value=shiftl(int(i,int64),56)
      do bit=1,8
        if (btest(value,63)) then
          value=ieor(shiftl(value,1),int(z'42F0E1EBA9EA3693',int64))
        else
          value=shiftl(value,1)
        end if
      end do
      crc_table(i)=value
    end do
    table_ready=.true.
  end subroutine prepare_crc

  integer(int64) function checkpoint_crc64(bytes,initial) result(crc)
    character, intent(in) :: bytes(:)
    integer(int64), intent(in), optional :: initial
    integer :: i,index
    call prepare_crc()
    crc=0_int64
    if (present(initial)) crc=initial
    do i=1,size(bytes)
      index=int(ieor(shiftr(crc,56),int(iachar(bytes(i)),int64)))
      crc=ieor(shiftl(crc,8),crc_table(index))
    end do
  end function checkpoint_crc64

  subroutine file_signature(path,signature,status)
    character(len=*), intent(in) :: path
    integer(int64), intent(out) :: signature(2)
    integer, intent(out) :: status
    character :: buffer(65536)
    integer :: unit,amount,close_status
    integer(int64) :: remaining
    signature=0_int64
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=status)
    if (status/=0) return
    inquire(unit=unit,size=remaining,iostat=status)
    if (status==0) then
      signature(1)=remaining
      do while(remaining>0_int64)
        amount=int(min(remaining,int(size(buffer),int64)))
        read(unit,iostat=status) buffer(:amount)
        if (status/=0) exit
        signature(2)=checkpoint_crc64(buffer(:amount),signature(2))
        remaining=remaining-amount
      end do
    end if
    close(unit,iostat=close_status)
    if (status==0) status=close_status
  end subroutine file_signature

  logical function io_agreement(local_ok) result(ok)
    logical, intent(in) :: local_ok
#ifdef USE_MPI
    integer :: status
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(local_ok,ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_checkpoint:io_agreement:1',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'checkpoint IO agreement')
#else
    ok=local_ok
#endif
  end function io_agreement

  subroutine finish_io(ok,status,operation)
    logical, intent(in) :: ok
    integer, intent(out), optional :: status
    character(len=*), intent(in) :: operation
    if (parallel_rank==0.and..not.ok)write(*,'(2A)')'CHECKPOINT ERROR: ',operation
    if (present(status)) then
      status=merge(0,1,ok)
    else
      call parallel_require_collective(ok,operation)
    end if
  end subroutine finish_io

  subroutine agree_values(values,message)
    integer(int64), intent(in) :: values(:)
    character(len=*), intent(in) :: message
#ifdef USE_MPI
    integer(int64) :: minimum(size(values)),maximum(size(values))
    integer :: status
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(values,minimum,size(values),count_type,MPI_MIN,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_checkpoint:agree_values:2',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'checkpoint controls minimum')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(values,maximum,size(values),count_type,MPI_MAX,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_checkpoint:agree_values:3',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'checkpoint controls maximum')
    call parallel_require_collective(all(minimum==maximum),message)
#endif
  end subroutine agree_values

  subroutine configuration_bytes(layout,bytes,status)
    type(domain_layout), intent(in) :: layout
    character, allocatable, intent(out) :: bytes(:)
    integer, intent(out) :: status
    integer :: unit,length,close_status
    open(newunit=unit,status='scratch',access='stream',form='unformatted',iostat=status)
    if (status/=0) return
    call write_checkpoint_metadata(unit,status,layout%global_shape,.true.)
    ! The identity map keeps the original configuration record byte for byte.
    if (status==0.and.polar_mesh_h/=1.0d0) write(unit,iostat=status) POLAR_CONFIG_MAGIC,polar_mesh_h
    if (status==0) inquire(unit=unit,size=length,iostat=status)
    if (status==0) allocate(bytes(length),stat=status)
    if (status==0) read(unit,pos=1,iostat=status) bytes
    close(unit,iostat=close_status)
    if (status==0) status=close_status
  end subroutine configuration_bytes

  subroutine agree_path(path,has_status,operation)
    character(len=*), intent(in) :: path
    logical, intent(in) :: has_status
    integer, intent(in) :: operation
    character(len=512) :: root_path
#ifdef USE_MPI
    integer :: status
#endif
    call parallel_require_active()
    call parallel_require_collective(len_trim(path)>0.and.len_trim(path)<=440.and.index(path,c_null_char)==0, &
      'invalid checkpoint generation path')
    root_path=path
#ifdef USE_MPI
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(0_int64)/8,count_type,status)
    call parallel_check_mpi(status,'checkpoint int64 MPI datatype')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(root_path,len(root_path),MPI_CHARACTER,0,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_checkpoint:agree_path:216',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'checkpoint path agreement')
#endif
    call parallel_require_collective(path==root_path,'checkpoint path differs between ranks')
    call agree_values([merge(1_int64,0_int64,has_status),int(operation,int64)],'checkpoint call differs between ranks')
  end subroutine agree_path

  integer(int64) function payload_bytes(shape) result(bytes)
    integer, intent(in) :: shape(3)
    integer :: d
    bytes=2_int64*neq*(storage_size(1.0d0)/8)
    do d=1,3
      if(shape(d)<1)then
        bytes=-1;return
      end if
      if(bytes>(huge(0_int64)-4096_int64)/shape(d))then
        bytes=-1;return
      end if
      bytes=bytes*shape(d)
    end do
  end function payload_bytes

  subroutine validate_call(path,layout,has_status,loading)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    logical, intent(in) :: has_status,loading
    character(len=512) :: root_path
    character(len=160) :: message
    logical :: valid
#ifdef USE_MPI
    integer :: status
#endif
    call parallel_require_active()
    call validate_domain(layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    valid=len_trim(path)>0.and.len_trim(path)<=440.and.index(path,c_null_char)==0
    call parallel_require_collective(valid,'invalid checkpoint generation path')
    root_path=path
#ifdef USE_MPI
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(0_int64)/8,count_type,status)
    call parallel_check_mpi(status,'checkpoint int64 MPI datatype')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(root_path,len(root_path),MPI_CHARACTER,0,MPI_COMM_WORLD,status)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_checkpoint:agree_path:257',collective_start)
    end block
#endif
    call parallel_check_mpi(status,'checkpoint path agreement')
#endif
    call parallel_require_collective(path==root_path,'checkpoint path differs between ranks')
    valid=layout%rank==parallel_rank.and.layout%ranks==parallel_size.and.layout%halo==nghost.and. &
      all(layout%local_shape==[nx,ny,nz]).and.all(layout%global_shape==[nx_global,ny_global,nz_global]).and. &
      all(layout%offset(1:2)==[x_global_offset,y_global_offset])
    call parallel_require_collective(valid,'checkpoint layout does not match allocated grid')
    valid=allocated(up).and.allocated(u).and.allocated(p).and.allocated(gamma_c)
    if (valid) valid=all(shape(up)==[nx,ny,nz,neq]).and.all(shape(u)==[nx,ny,nz,neq]).and. &
      all(shape(p)==[nx+2*nghost+1,ny+2*nghost+1,nz+2*nghost+1,neq]).and. &
      all(shape(gamma_c)==[3,3,nx,ny])
    call parallel_require_collective(valid,'checkpoint arrays or metric cache are not prepared')
    call parallel_require_collective(.not.accelerator_enabled.and..not.is_mhd,'checkpoint supports CPU hydrodynamics only')
    valid=all(ieee_is_finite([integration_time,dt,adb_idx,g1]))
    if(valid) valid=integration_time>=0.0d0.and.dt>=0.0d0.and.n_steps>=0.and.adb_idx>1.0d0
    if(valid) valid=g1==adb_idx/(adb_idx-1.0d0)
    call parallel_require_collective(valid,'invalid checkpoint epoch or equation of state')
    call agree_values([int(layout%global_shape,int64),int(layout%process_shape,int64),int(nghost,int64), &
      merge(1_int64,0_int64,layout%periodic_x),merge(1_int64,0_int64,layout%periodic_y), &
      merge(1_int64,0_int64,has_status),merge(1_int64,0_int64,loading)],'checkpoint call differs between ranks')
    call validate_checkpoint_boundary()
    if (loading) then
      call parallel_require_collective(.not.diagnostics_active,'checkpoint observer restoration requires C03; finalize D01 first')
      call parallel_require_collective(halo_layout_ready(layout),'checkpoint recovery requires prepared halos')
    end if
  end subroutine validate_call

  logical function physical_state(primitives,conserved) result(valid)
    real*8, intent(in) :: primitives(:,:,:,:),conserved(:,:,:,:)
    integer :: i,j,k
    real*8 :: velocity(3),vsq
    valid=all(ieee_is_finite(primitives)).and.all(ieee_is_finite(conserved)).and.all(ieee_is_finite(gamma_c))
    if (.not.valid) return
    valid=all(primitives(:,:,:,eq_de)>0.0d0).and.all(primitives(:,:,:,eq_pr)>0.0d0).and. &
      all(conserved(:,:,:,eq_de)>0.0d0)
    if (.not.valid) return
    do k=1,nz
      do j=1,ny
        do i=1,nx
          velocity=primitives(i,j,k,eq_vx:eq_vz)
          vsq=dot_product(velocity,matmul(gamma_c(:,:,i,j),velocity))
          valid=valid.and.ieee_is_finite(vsq)
          if (ieee_is_finite(vsq)) valid=valid.and.vsq>=0.0d0.and.vsq<1.0d0
        end do
      end do
    end do
  end function physical_state

  function block_path(path,rank) result(name)
    character(len=*), intent(in) :: path
    integer, intent(in) :: rank
    character(len=512) :: name
    write(name,'(A,A,I9.9,A)')trim(path),'/block_',rank,'.bin'
  end function block_path

  integer(int64) function expected_block_bytes(layout) result(bytes)
    type(domain_layout), intent(in) :: layout
    integer :: header
    inquire(iolength=header) BLOCK_MAGIC,ABI,layout%rank,layout%ranks,layout%global_shape, &
      layout%local_shape,layout%offset,layout%process_shape,layout%halo, &
      [layout%periodic_x,layout%periodic_y],integration_time,dt,n_steps,perturbation_applied, &
      numerical_events,numerical_delta_u,END_MAGIC
    bytes=payload_bytes(layout%local_shape)
    if(bytes>=0)bytes=bytes+int(header,int64)
  end function expected_block_bytes

  subroutine report_checkpoint_io(path,operation,status,message)
    character(len=*), intent(in) :: path,operation,message
    integer, intent(in) :: status
    if(status/=0)write(error_unit,'(A,I0,5A,I0,2A)') &
      'CHECKPOINT I/O rank=',parallel_rank,' ',operation,' ',trim(path),' iostat=',status,' ',trim(message)
  end subroutine report_checkpoint_io

  subroutine write_block(path,layout,status)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    integer, intent(out) :: status
    integer :: unit,close_status
    character(len=512) :: message
    open(newunit=unit,file=path,access='stream',form='unformatted',status='new',action='write',iostat=status,iomsg=message)
    if(status/=0)then
      call report_checkpoint_io(path,'open',status,message)
      return
    end if
    write(unit,iostat=status,iomsg=message) BLOCK_MAGIC,ABI,layout%rank,layout%ranks,layout%global_shape, &
      layout%local_shape,layout%offset,layout%process_shape,layout%halo, &
      [layout%periodic_x,layout%periodic_y],integration_time,dt,n_steps,perturbation_applied
    if(status==0)write(unit,iostat=status,iomsg=message) numerical_events,numerical_delta_u
    if(status==0)write(unit,iostat=status,iomsg=message) up,p(1:nx,1:ny,1:nz,:)
    if(status==0)write(unit,iostat=status,iomsg=message) END_MAGIC
    if(status/=0)call report_checkpoint_io(path,'write',status,message)
    close(unit,iostat=close_status,iomsg=message)
    if(close_status/=0)call report_checkpoint_io(path,'close',close_status,message)
    if(status==0)status=close_status
  end subroutine write_block

  subroutine write_manifest(path,manifest,status)
    character(len=*), intent(in) :: path
    type(checkpoint_manifest), intent(in) :: manifest
    integer, intent(out) :: status
    integer :: unit,close_status
    character(len=512) :: message
    open(newunit=unit,file=path,access='stream',form='unformatted',status='new',action='write',iostat=status,iomsg=message)
    if(status/=0)then
      call report_checkpoint_io(path,'open',status,message)
      return
    end if
    write(unit,iostat=status,iomsg=message) MANIFEST_MAGIC,ABI,manifest%global_shape,manifest%process_shape, &
      manifest%halo,manifest%ranks,neq,manifest%periodic,manifest%time,manifest%step_size, &
      manifest%step,manifest%applied,manifest%configuration
    if(status==0)write(unit,iostat=status,iomsg=message) manifest%blocks
    if(status==0)write(unit,iostat=status,iomsg=message) END_MAGIC
    if(status/=0)call report_checkpoint_io(path,'write',status,message)
    close(unit,iostat=close_status,iomsg=message)
    if(close_status/=0)call report_checkpoint_io(path,'close',close_status,message)
    if(status==0)status=close_status
  end subroutine write_manifest

  subroutine save_parallel_checkpoint(path,layout,status)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    integer, intent(out), optional :: status
    type(checkpoint_manifest) :: manifest
    character, allocatable :: configuration(:)
    real*8, allocatable :: checked_u(:,:,:,:),checked_p(:,:,:,:)
    integer(int64) :: signature(2),descriptor(9),manifest_signature(2)
    integer(int64) :: checked_events(NUM_EVENT_COUNT),checked_signature(2)
    real*8 :: checked_delta(neq)
    character(len=512) :: name,message
    integer :: ios,unit,close_status,manifest_header
    logical :: ok
#ifdef USE_MPI
    integer :: mpi_result
#endif
    call validate_call(path,layout,present(status),.false.)
    call parallel_require_collective(physical_state(p(1:nx,1:ny,1:nz,:),up),'invalid physical checkpoint state')
    call parallel_require_collective(all(numerical_events>=0_int64).and.all(ieee_is_finite(numerical_delta_u)), &
      'invalid checkpoint counters')
    call configuration_bytes(layout,configuration,ios)
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    signature=[int(size(configuration),int64),checkpoint_crc64(configuration)]
    call agree_values([signature,int(n_steps,int64),transfer(integration_time,0_int64),transfer(dt,0_int64), &
      merge(1_int64,0_int64,perturbation_applied)],'checkpoint configuration or epoch differs between ranks')
    ios=0
    if(parallel_rank==0)ios=posix_mkdir(trim(path)//c_null_char,int(o'700',c_int))
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    if(parallel_rank==0)then
      name=trim(path)//'/configuration.bin'
      open(newunit=unit,file=trim(name),status='new',access='stream', &
        form='unformatted',action='write',iostat=ios,iomsg=message)
      if(ios==0)then
        write(unit,iostat=ios,iomsg=message)configuration
        if(ios/=0)call report_checkpoint_io(name,'write',ios,message)
        close(unit,iostat=close_status,iomsg=message)
        if(close_status/=0)call report_checkpoint_io(name,'close',close_status,message)
        if(ios==0)ios=close_status
      else
        call report_checkpoint_io(name,'open',ios,message)
      end if
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    name=trim(path)//'/configuration.bin'
    call file_signature(trim(name),signature,ios)
    if(ios==0)then
      if(any(signature/=[int(size(configuration),int64),checkpoint_crc64(configuration)]))ios=1
    end if
    if(ios/=0)call report_checkpoint_io(name,'verify',ios,'configuration size or checksum could not be verified')
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    name=block_path(path,parallel_rank)
    call write_block(trim(name),layout,ios)
    if(ios==0)then
      call file_signature(trim(name),signature,ios)
      if(ios==0.and.signature(1)/=expected_block_bytes(layout))then
        ios=1
        write(message,'(A,I0,A,I0)') 'block size ',signature(1),', expected ',expected_block_bytes(layout)
      else
        message='could not read block signature'
      end if
      if(ios/=0)call report_checkpoint_io(name,'verify',ios,message)
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    descriptor=[int(parallel_rank,int64),int(layout%local_shape,int64),int(layout%offset,int64),signature]
    allocate(manifest%blocks(9,parallel_size),stat=ios)
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    manifest%blocks=0_int64
    manifest%global_shape=layout%global_shape;manifest%process_shape=layout%process_shape
    manifest%halo=layout%halo;manifest%ranks=parallel_size
    manifest%periodic=[layout%periodic_x,layout%periodic_y]
    manifest%time=integration_time;manifest%step_size=dt;manifest%step=n_steps;manifest%applied=perturbation_applied
    manifest%blocks(:,parallel_rank+1)=descriptor
    call read_block(trim(name),layout,manifest,checked_u,checked_p,checked_events,checked_delta,ios)
    if(ios==0)then
      if(.not.all(checked_u==up).or..not.all(checked_p==p(1:nx,1:ny,1:nz,:)).or. &
        any(checked_events/=numerical_events).or.any(checked_delta/=numerical_delta_u))ios=1
    end if
    if(ios/=0)call report_checkpoint_io(name,'verify',ios,'block readback failed')
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    deallocate(checked_u,checked_p)
#ifdef USE_MPI
    call MPI_Gather(descriptor,9,count_type,manifest%blocks,9,count_type,0,MPI_COMM_WORLD,mpi_result)
    call parallel_check_mpi(mpi_result,'checkpoint block descriptors')
#else
    manifest%blocks(:,1)=descriptor
#endif
    ios=0
    if(parallel_rank==0)then
      manifest%configuration=[int(size(configuration),int64),checkpoint_crc64(configuration)]
      call write_manifest(trim(path)//'/manifest.bin',manifest,ios)
      if(ios==0)call file_signature(trim(path)//'/manifest.bin',manifest_signature,ios)
      if(ios==0)then
        inquire(iolength=manifest_header) MANIFEST_MAGIC,ABI,manifest%global_shape,manifest%process_shape, &
          manifest%halo,manifest%ranks,neq,manifest%periodic,manifest%time,manifest%step_size, &
          manifest%step,manifest%applied,manifest%configuration,END_MAGIC
        if(manifest_signature(1)/=int(manifest_header,int64)+ &
          int(storage_size(0_int64)/8,int64)*size(manifest%blocks,kind=int64))ios=1
      end if
      if(ios==0)then
        name=trim(path)//'/COMMITTED.tmp'
        open(newunit=unit,file=trim(name),status='new',access='stream', &
          form='unformatted',action='write',iostat=ios,iomsg=message)
        if(ios==0)then
          write(unit,iostat=ios,iomsg=message) COMMIT_MAGIC,manifest_signature,END_MAGIC
          if(ios/=0)call report_checkpoint_io(name,'write',ios,message)
          close(unit,iostat=close_status,iomsg=message)
          if(close_status/=0)call report_checkpoint_io(name,'close',close_status,message)
          if(ios==0)ios=close_status
        else
          call report_checkpoint_io(name,'open',ios,message)
        end if
        if(ios==0)call read_commit(trim(path)//'/COMMITTED.tmp',checked_signature,ios)
        if(ios==0)then
          if(any(checked_signature/=manifest_signature))ios=1
        end if
        if(ios==0)call rename(trim(path)//'/COMMITTED.tmp',trim(path)//'/COMMITTED',ios)
      end if
    end if
    ok=io_agreement(ios==0)
900 continue
    call finish_io(ok,status,'checkpoint generation not published: '//trim(path))
  end subroutine save_parallel_checkpoint

  subroutine read_commit(path,signature,status)
    character(len=*), intent(in) :: path
    integer(int64), intent(out) :: signature(2)
    integer, intent(out) :: status
    character(len=16) :: magic,footer
    character :: trailing
    integer :: unit,close_status,tail_status
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)magic,signature,footer
    if(status==0)then
      if(magic/=COMMIT_MAGIC.or.footer/=END_MAGIC.or.signature(1)<=0)status=1
    end if
    if(status==0)then
      read(unit,iostat=tail_status)trailing
      if(tail_status>=0)status=1
    end if
    close(unit,iostat=close_status)
    if(status==0)status=close_status
  end subroutine read_commit

  subroutine read_manifest_data(path,manifest,status)
    character(len=*), intent(in) :: path
    type(checkpoint_manifest), intent(out) :: manifest
    integer, intent(out) :: status
    integer(int64) :: expected(2),actual(2)
    character(len=16) :: magic,footer
    character(len=160) :: message
    character :: trailing
    type(domain_layout) :: block
    integer :: unit,close_status,tail_status,r,header
    call read_commit(trim(path)//'/COMMITTED',expected,status)
    if(status/=0)return
    manifest%signature=expected
    call file_signature(trim(path)//'/manifest.bin',actual,status)
    if(status/=0)return
    if(any(actual/=expected))then
      status=1;return
    end if
    open(newunit=unit,file=trim(path)//'/manifest.bin',access='stream',form='unformatted', &
      status='old',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)magic,manifest%abi_tag,manifest%global_shape,manifest%process_shape, &
      manifest%halo,manifest%ranks,manifest%equations,manifest%periodic,manifest%time,manifest%step_size, &
      manifest%step,manifest%applied,manifest%configuration
    if(status==0)then
      if(magic/=MANIFEST_MAGIC.or.any(manifest%abi_tag/=ABI).or.manifest%equations/=neq.or. &
        manifest%ranks<1.or.manifest%configuration(1)<=0)status=1
    end if
    if(status==0)then
      call build_domain(manifest%global_shape,manifest%halo,0,manifest%ranks,manifest%periodic(1),block,message, &
        manifest%process_shape,manifest%periodic(2))
      if(len_trim(message)/=0.or.payload_bytes(manifest%global_shape)<0)status=1
    end if
    if(status==0)then
      if(.not.all(ieee_is_finite([manifest%time,manifest%step_size])))status=1
    end if
    if(status==0)then
      if(manifest%time<0.or.manifest%step_size<0.or.manifest%step<0)status=1
    end if
    if(status==0)then
      inquire(iolength=header) MANIFEST_MAGIC,ABI,manifest%global_shape,manifest%process_shape, &
        manifest%halo,manifest%ranks,neq,manifest%periodic,manifest%time,manifest%step_size, &
        manifest%step,manifest%applied,manifest%configuration,END_MAGIC
      if(actual(1)/=int(header,int64)+72_int64*manifest%ranks)status=1
    end if
    if(status==0)allocate(manifest%blocks(9,manifest%ranks),stat=status)
    if(status==0)read(unit,iostat=status)manifest%blocks
    if(status==0)read(unit,iostat=status)footer
    if(status==0)then
      read(unit,iostat=tail_status)trailing
      if(footer/=END_MAGIC.or.tail_status>=0)status=1
    end if
    close(unit,iostat=close_status)
    if(status==0)status=close_status
    if(status/=0)return
    do r=0,manifest%ranks-1
      call build_domain(manifest%global_shape,manifest%halo,r,manifest%ranks,manifest%periodic(1),block,message, &
        manifest%process_shape,manifest%periodic(2))
      if(any(manifest%blocks(1:7,r+1)/=[int(r,int64),int(block%local_shape,int64),int(block%offset,int64)]).or. &
        manifest%blocks(8,r+1)/=expected_block_bytes(block))status=1
    end do
  end subroutine read_manifest_data

  subroutine read_manifest(path,layout,manifest,status,redistribute)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    type(checkpoint_manifest), intent(out) :: manifest
    integer, intent(out) :: status
    logical, intent(in), optional :: redistribute
    logical :: same_only
    same_only=.true.
    if(present(redistribute))same_only=.not.redistribute
    call read_manifest_data(path,manifest,status)
    if(status/=0)return
    if(any(manifest%global_shape/=layout%global_shape).or.manifest%halo/=nghost.or. &
      any(manifest%periodic.neqv.[layout%periodic_x,layout%periodic_y]))status=1
    if(same_only)then
      if(any(manifest%process_shape/=layout%process_shape).or.manifest%ranks/=parallel_size)status=1
    end if
  end subroutine read_manifest

  subroutine read_parallel_checkpoint_configuration(path,layout,status,process_shape)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(out) :: layout
    integer, intent(out), optional :: status
    integer, intent(in), optional :: process_shape(2)
    type(checkpoint_manifest) :: manifest
    integer(int64) :: signature(2)
    integer :: unit,ios,close_status,parts(2)
    integer(int64) :: position
    character(len=len(POLAR_CONFIG_MAGIC)) :: polar_magic
    real*8 :: saved_h
    logical :: ok
    character(len=160) :: message
    layout=domain_layout()
    call agree_path(path,present(status),3)
    call parallel_require_collective(.not.allocated(up).and..not.allocated(p).and..not.allocated(u).and. &
      .not.diagnostics_active,'configuration loading requires an empty grid and inactive diagnostics')
    call read_manifest_data(path,manifest,ios)
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    call agree_values(manifest%signature,'checkpoint generation differs between ranks')
    call file_signature(trim(path)//'/configuration.bin',signature,ios)
    if(ios==0)then
      if(any(signature/=manifest%configuration))ios=1
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    open(newunit=unit,file=trim(path)//'/configuration.bin',access='stream',form='unformatted', &
      status='old',action='read',iostat=ios)
    if(ios==0)then
      call read_checkpoint_metadata(unit,ios)
      saved_h=1.0d0
      if (ios==0) inquire(unit=unit,pos=position,iostat=ios)
      if (ios==0) then
        select case(manifest%configuration(1)-position+1)
        case(0_int64)
        case(len(POLAR_CONFIG_MAGIC)+storage_size(saved_h)/8)
          read(unit,iostat=ios) polar_magic,saved_h
          if (ios==0) then
            if (polar_magic/=POLAR_CONFIG_MAGIC.or..not.ieee_is_finite(saved_h)) ios=1
          end if
          if (ios==0) then
            if (saved_h<=0.0d0.or.saved_h>=1.0d0) ios=1
          end if
        case default
          ios=1
        end select
      end if
      close(unit,iostat=close_status)
      if(ios==0)ios=close_status
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    ok=io_agreement(all([nx,ny,nz]==manifest%global_shape).and.nghost==manifest%halo)
    if(.not.ok)goto 900
    polar_mesh_h=saved_h
    call validate_parameters(path)
    g1=adb_idx/(adb_idx-1.0d0);call set_metric_type()
    do_restart=.false.;do_radial_mass_budget=.false.
    integration_time=manifest%time;dt=manifest%step_size;n_steps=manifest%step
    perturbation_applied=manifest%applied
    parts=[parallel_size,1]
    if(present(process_shape))parts=process_shape
    call agree_values(int(parts,int64),'checkpoint process grid differs between ranks')
    call build_domain([nx,ny,nz],nghost,parallel_rank,parallel_size,manifest%periodic(1),layout,message, &
      parts,manifest%periodic(2))
    ok=io_agreement(len_trim(message)==0)
900 continue
    call finish_io(ok,status,'checkpoint configuration rejected: '//trim(path))
  end subroutine read_parallel_checkpoint_configuration

  subroutine read_block(path,layout,manifest,saved_u,saved_p,events,delta,status)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    type(checkpoint_manifest), intent(in) :: manifest
    real*8, allocatable, intent(out) :: saved_u(:,:,:,:),saved_p(:,:,:,:)
    integer(int64), intent(out) :: events(NUM_EVENT_COUNT)
    real*8, intent(out) :: delta(neq)
    integer, intent(out) :: status
    character(len=16) :: magic,footer
    character :: trailing
    integer :: unit,close_status,tail_status,abi_tag(5),rank,ranks,global_shape(3),local_shape(3), &
      offset(3),process_shape(2),halo,step
    integer(int64) :: signature(2)
    logical :: periodic(2),applied
    real*8 :: time,step_size
    call file_signature(path,signature,status)
    if(status/=0)return
    if(any(signature/=manifest%blocks(8:9,parallel_rank+1)))then
      status=1;return
    end if
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)magic,abi_tag,rank,ranks,global_shape,local_shape,offset,process_shape,halo, &
      periodic,time,step_size,step,applied
    if(status==0)then
      if(magic/=BLOCK_MAGIC.or.any(abi_tag/=ABI).or.rank/=parallel_rank.or.ranks/=parallel_size.or. &
        any(global_shape/=layout%global_shape).or.any(local_shape/=layout%local_shape).or. &
        any(offset/=layout%offset).or.any(process_shape/=layout%process_shape).or.halo/=nghost.or. &
        any(periodic.neqv.manifest%periodic).or.(applied.neqv.manifest%applied))status=1
    end if
    if(status==0)then
      if(.not.all(ieee_is_finite([time,step_size])))status=1
    end if
    if(status==0)then
      if(time/=manifest%time.or.step_size/=manifest%step_size.or.step/=manifest%step)status=1
    end if
    if(status==0)read(unit,iostat=status)events,delta
    if(status==0)allocate(saved_u(nx,ny,nz,neq),saved_p(nx,ny,nz,neq),stat=status)
    if(status==0)read(unit,iostat=status)saved_u,saved_p
    if(status==0)read(unit,iostat=status)footer
    if(status==0)then
      read(unit,iostat=tail_status)trailing
      if(footer/=END_MAGIC.or.tail_status>=0)status=1
    end if
    close(unit,iostat=close_status)
    if(status==0)status=close_status
    if(status/=0)return
    if(any(events<0_int64).or..not.all(ieee_is_finite(delta)))status=1
    if(.not.physical_state(saved_p,saved_u))status=1
  end subroutine read_block

  subroutine open_block_header(path,layout,manifest,unit,payload,events,delta,status)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    type(checkpoint_manifest), intent(in) :: manifest
    integer, intent(out) :: unit,status
    integer(int64), intent(out) :: payload,events(NUM_EVENT_COUNT)
    real*8, intent(out) :: delta(neq)
    character(len=16) :: magic,footer
    integer :: abi_tag(5),rank,ranks,global_shape(3),local_shape(3),offset(3),parts(2),halo,step,close_status
    integer(int64) :: bytes
    logical :: periodic(2),applied
    real*8 :: time,step_size
    open(newunit=unit,file=path,access='stream',form='unformatted',status='old',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)magic,abi_tag,rank,ranks,global_shape,local_shape,offset,parts,halo, &
      periodic,time,step_size,step,applied
    if(status==0)then
      if(magic/=BLOCK_MAGIC.or.any(abi_tag/=ABI).or.rank/=layout%rank.or.ranks/=layout%ranks.or. &
        any(global_shape/=layout%global_shape).or.any(local_shape/=layout%local_shape).or. &
        any(offset/=layout%offset).or.any(parts/=layout%process_shape).or.halo/=layout%halo.or. &
        any(periodic.neqv.manifest%periodic).or.(applied.neqv.manifest%applied))status=1
    end if
    if(status==0)then
      if(.not.all(ieee_is_finite([time,step_size])))status=1
    end if
    if(status==0)then
      if(time/=manifest%time.or.step_size/=manifest%step_size.or.step/=manifest%step)status=1
    end if
    if(status==0)read(unit,iostat=status)events,delta
    if(status==0)inquire(unit=unit,pos=payload,size=bytes,iostat=status)
    if(status==0)then
      if(bytes/=expected_block_bytes(layout).or.any(events<0_int64).or..not.all(ieee_is_finite(delta)))status=1
    end if
    if(status==0)read(unit,pos=bytes-15_int64,iostat=status)footer
    if(status==0)then
      if(footer/=END_MAGIC)status=1
    end if
    if(status/=0)close(unit,iostat=close_status)
  end subroutine open_block_header

  subroutine add_counter_history(total,delta_total,events,delta,status)
    integer(int64), intent(inout) :: total(NUM_EVENT_COUNT)
    real*8, intent(inout) :: delta_total(neq)
    integer(int64), intent(in) :: events(NUM_EVENT_COUNT)
    real*8, intent(in) :: delta(neq)
    integer, intent(out) :: status
    integer :: q
    status=1
    if(any(events<0_int64).or.any(events>huge(0_int64)-total))return
    do q=1,neq
      if(delta(q)>=0.0d0)then
        if(delta_total(q)>huge(1.0d0)-delta(q))return
      else
        if(delta_total(q)<-huge(1.0d0)-delta(q))return
      end if
    end do
    total=total+events;delta_total=delta_total+delta;status=0
  end subroutine add_counter_history

  subroutine redistribute_blocks(path,layout,manifest,saved_u,saved_p,events,delta,status)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    type(checkpoint_manifest), intent(in) :: manifest
    real*8, allocatable, intent(out) :: saved_u(:,:,:,:),saved_p(:,:,:,:)
    integer(int64), intent(out) :: events(NUM_EVENT_COUNT)
    real*8, intent(out) :: delta(neq)
    integer, intent(out) :: status
    type(domain_layout) :: block
    integer, allocatable :: covered(:,:)
    integer(int64) :: signature(2),payload,offset,old_events(NUM_EVENT_COUNT)
    real*8 :: old_delta(neq)
    integer :: r,unit,close_status,lo(3),hi(3),q,k,j
    character(len=160) :: message
    status=0;events=0_int64;delta=0.0d0
    ! One verifier per stored file, independent of how many new blocks intersect it.
    do r=parallel_rank,manifest%ranks-1,parallel_size
      call file_signature(trim(block_path(path,r)),signature,status)
      if(status/=0)exit
      if(any(signature/=manifest%blocks(8:9,r+1)))then
        status=1;exit
      end if
    end do
    if(.not.io_agreement(status==0))then
      status=1;return
    end if
    allocate(saved_u(nx,ny,nz,neq),saved_p(nx,ny,nz,neq),covered(nx,ny),stat=status)
    if(.not.io_agreement(status==0))then
      status=1;return
    end if
    covered=0
    do r=0,manifest%ranks-1
      call build_domain(manifest%global_shape,manifest%halo,r,manifest%ranks,manifest%periodic(1),block,message, &
        manifest%process_shape,manifest%periodic(2))
      lo=max(layout%offset,block%offset)+1
      hi=min(layout%offset+layout%local_shape,block%offset+block%local_shape)
      if(any(lo>hi).and.parallel_rank/=0)cycle
      call open_block_header(trim(block_path(path,r)),block,manifest,unit,payload,old_events,old_delta,status)
      if(status/=0)exit
      if(parallel_rank==0)call add_counter_history(events,delta,old_events,old_delta,status)
      if(status==0.and.all(lo<=hi))then
        do q=1,neq
          do k=lo(3),hi(3)
            do j=lo(2),hi(2)
              offset=(((int(q-1,int64)*block%local_shape(3)+k-block%offset(3)-1)*block%local_shape(2)+ &
                j-block%offset(2)-1)*block%local_shape(1)+lo(1)-block%offset(1)-1)*8_int64
              read(unit,pos=payload+offset,iostat=status) &
                saved_u(lo(1)-layout%offset(1):hi(1)-layout%offset(1),j-layout%offset(2),k-layout%offset(3),q)
              if(status/=0)exit
              read(unit,pos=payload+payload_bytes(block%local_shape)/2+offset,iostat=status) &
                saved_p(lo(1)-layout%offset(1):hi(1)-layout%offset(1),j-layout%offset(2),k-layout%offset(3),q)
              if(status/=0)exit
            end do
            if(status/=0)exit
          end do
          if(status/=0)exit
        end do
        covered(lo(1)-layout%offset(1):hi(1)-layout%offset(1),lo(2)-layout%offset(2):hi(2)-layout%offset(2))= &
          covered(lo(1)-layout%offset(1):hi(1)-layout%offset(1),lo(2)-layout%offset(2):hi(2)-layout%offset(2))+1
      end if
      close(unit,iostat=close_status)
      if(status==0)status=close_status
      if(status/=0)exit
    end do
    if(status==0)then
      if(any(covered/=1))status=1
    end if
    if(status==0)then
      if(.not.physical_state(saved_p,saved_u))status=1
    end if
  end subroutine redistribute_blocks

  subroutine load_parallel_checkpoint(path,layout,status,redistribute)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    integer, intent(out), optional :: status
    logical, intent(in), optional :: redistribute
    type(checkpoint_manifest) :: manifest
    character, allocatable :: configuration(:),stored_configuration(:)
    real*8, allocatable :: saved_u(:,:,:,:),saved_p(:,:,:,:)
    integer(int64) :: signature(2),events(NUM_EVENT_COUNT)
    real*8 :: delta(neq)
    integer :: ios,unit,close_status
    logical :: ok,allow_redistribution,changed
    call validate_call(path,layout,present(status),.true.)
    allow_redistribution=.false.
    if(present(redistribute))allow_redistribution=redistribute
    call agree_values([merge(1_int64,0_int64,allow_redistribution)],'checkpoint redistribution option differs between ranks')
    call read_manifest(path,layout,manifest,ios,allow_redistribution)
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    call agree_values(manifest%signature,'checkpoint generation differs between ranks')
    call configuration_bytes(layout,configuration,ios)
    if(ios==0)then
      signature=[int(size(configuration),int64),checkpoint_crc64(configuration)]
      if(any(signature/=manifest%configuration))ios=1
    end if
    if(ios==0)call file_signature(trim(path)//'/configuration.bin',signature,ios)
    if(ios==0)then
      if(any(signature/=manifest%configuration))ios=1
    end if
    if(ios==0)allocate(stored_configuration(size(configuration)),stat=ios)
    if(ios==0)then
      open(newunit=unit,file=trim(path)//'/configuration.bin',status='old',access='stream', &
        form='unformatted',action='read',iostat=ios)
      if(ios==0)then
        read(unit,iostat=ios)stored_configuration
        close(unit,iostat=close_status)
        if(ios==0)ios=close_status
        if(ios==0)then
          if(any(stored_configuration/=configuration))ios=1
        end if
      end if
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    changed=manifest%ranks/=parallel_size.or.any(manifest%process_shape/=layout%process_shape)
    if(changed)then
      call redistribute_blocks(path,layout,manifest,saved_u,saved_p,events,delta,ios)
    else
      call read_block(trim(block_path(path,parallel_rank)),layout,manifest,saved_u,saved_p,events,delta,ios)
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    ! Commit memory only after every rank has validated its own complete block.
    up=saved_u;u=saved_u;p(1:nx,1:ny,1:nz,:)=saved_p
    integration_time=manifest%time;dt=manifest%step_size;n_steps=manifest%step
    perturbation_applied=manifest%applied;numerical_events=events;numerical_delta_u=delta
    if(changed.and.parallel_rank==0)write(*,'(A)') &
      'CHECKPOINT: changed partition; historical counter totals belong to rank 0; other ranks start at zero.'
    call fill_primitive_halos(p)
900 continue
    call finish_io(ok,status,'checkpoint recovery rejected: '//trim(path))
  end subroutine load_parallel_checkpoint

  subroutine root_file_signature(path,signature,status)
    character(len=*), intent(in) :: path
    integer(int64), intent(out) :: signature(2)
    integer, intent(out) :: status
#ifdef USE_MPI
    integer :: mpi_result
#endif
    status=0;signature=0_int64
    if(parallel_rank==0)call file_signature(path,signature,status)
    if(.not.io_agreement(status==0))then
      status=1;return
    end if
#ifdef USE_MPI
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(signature,2,count_type,0,MPI_COMM_WORLD,mpi_result)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_checkpoint:restore_history_signature:916',collective_start)
    end block
#endif
    call parallel_check_mpi(mpi_result,'serial checkpoint signature')
#endif
  end subroutine root_file_signature

  subroutine index_serial_file(path,index,status)
    character(len=*), intent(in) :: path
    type(serial_checkpoint_index), intent(out) :: index
    integer, intent(out) :: status
    integer(int64) :: positions(33),lengths(33),position,bytes
    integer :: unit,n,marker,trailer,close_status,metadata_records,expected_records,saved_step,dims(4)
    character(len=12) :: magic,footer
    real*8 :: saved_time
    index=serial_checkpoint_index()
    open(newunit=unit,file=path,status='old',access='stream',form='unformatted',action='read',iostat=status)
    if(status/=0)return
    inquire(unit=unit,size=bytes,iostat=status)
    position=1;n=0
    do while(status==0.and.position<=bytes)
      if(n==size(positions).or.bytes-position<7)then
        status=1;exit
      end if
      read(unit,pos=position,iostat=status)marker
      if(status/=0)exit
      ! Native 4-byte record markers; segmented records are explicitly unsupported.
      if(marker<0.or.int(marker,int64)>bytes-position-7)then
        status=1;exit
      end if
      n=n+1;positions(n)=position+4;lengths(n)=marker
      read(unit,pos=position+4_int64+marker,iostat=status)trailer
      if(status/=0)exit
      if(trailer/=marker)then
        status=1;exit
      end if
      position=position+8_int64+marker
    end do
    if(status==0)then
      if(n<1.or.position/=bytes+1)status=1
    end if
    if(status==0)then
      if(lengths(1)/=12)status=1
    end if
    if(status==0)read(unit,pos=positions(1),iostat=status)magic
    if(status==0)then
      select case(magic)
      case('GRHDCP_V7');index%version=7;metadata_records=24
      case('GRHDCP_V8');index%version=8;metadata_records=27
      case('GRHDCP_V9');index%version=9;metadata_records=29
      case default;status=1
      end select
    end if
    if(status==0)then
      expected_records=metadata_records+3
      if(index%version>=8)expected_records=expected_records+1
      if(n/=expected_records.or.lengths(2)/=16)status=1
    end if
    if(status==0)then
      read(unit,pos=positions(2),iostat=status)trailer,index%shape
      if(status==0)then
        if(trailer<1.or.payload_bytes(index%shape)<0)status=1
      end if
    end if
    if(status==0)then
      if(lengths(metadata_records+1)/=8.or.lengths(metadata_records+2)/=4.or. &
        lengths(metadata_records+3)/=payload_bytes(index%shape)/2)status=1
    end if
    if(status==0)read(unit,pos=positions(metadata_records+1),iostat=status)index%time
    if(status==0)read(unit,pos=positions(metadata_records+2),iostat=status)index%step
    if(status==0)then
      if(.not.ieee_is_finite(index%time))status=1
    end if
    if(status==0)then
      if(index%time<0.0d0.or.index%step<0)status=1
    end if
    if(status==0.and.index%version>=8)then
      if(lengths(n)/=40)then
        status=1
      else
        read(unit,pos=positions(n),iostat=status)footer,saved_step,saved_time,dims
        if(status==0)then
          if(.not.ieee_is_finite(saved_time))status=1
        end if
        if(status==0)then
          if((index%version==8.and.footer/='GRIAL_END_V8').or. &
            (index%version==9.and.footer/='GRIAL_END_V9').or.saved_step/=index%step.or. &
            saved_time/=index%time.or.any(dims/=[index%shape,neq]))status=1
        end if
      end if
    end if
    if(status==0)index%payload=positions(metadata_records+3)
    close(unit,iostat=close_status)
    if(status==0)status=close_status
  end subroutine index_serial_file

  subroutine read_serial_checkpoint_configuration(path,layout,status,process_shape)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(out) :: layout
    integer, intent(out), optional :: status
    integer, intent(in), optional :: process_shape(2)
    type(serial_checkpoint_index) :: index
    integer :: unit,ios,close_status,parts(2)
    logical :: ok,periodic
    character(len=160) :: message
    layout=domain_layout();serial_source_path=''
    polar_mesh_h=1.0d0
    if(allocated(serial_configuration))deallocate(serial_configuration)
    call agree_path(path,present(status),4)
    call parallel_require_collective(.not.allocated(up).and..not.allocated(p).and..not.allocated(u).and. &
      .not.diagnostics_active,'configuration loading requires an empty grid and inactive diagnostics')
    call root_file_signature(path,index%signature,ios)
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    block
      integer(int64) :: signature(2)
      signature=index%signature
      call index_serial_file(path,index,ios)
      index%signature=signature
    end block
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    open(newunit=unit,file=path,status='old',form='unformatted',action='read',iostat=ios)
    if(ios==0)then
      call read_checkpoint_metadata(unit,ios)
      close(unit,iostat=close_status)
      if(ios==0)ios=close_status
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    if(index%version==7)then
      checkpoint_interval=100.0d0;extraction_stride=10;mass_monitor_stride=1000
      do_numerical_diagnostics=.false.;vtk_mapping_id=VTK_MAP_PHYSICAL
      walltime_limit_seconds=0;walltime_margin_seconds=900;walltime_checkpoint_seconds=0
    end if
    ok=io_agreement(all([nx,ny,nz]==index%shape))
    if(.not.ok)goto 900
    call choose_numerical_architecture();call validate_parameters(path)
    g1=adb_idx/(adb_idx-1.0d0);call set_metric_type()
    integration_time=index%time;n_steps=index%step;dt=0.0d0
    do_restart=.false.;do_radial_mass_budget=.false.
    index%events=numerical_events;index%delta=numerical_delta_u;index%applied=perturbation_applied
    parts=[parallel_size,1]
    if(present(process_shape))parts=process_shape
    call agree_values(int(parts,int64),'checkpoint process grid differs between ranks')
    periodic=case_id==4.or.case_id==6
    call build_domain([nx,ny,nz],nghost,parallel_rank,parallel_size,periodic,layout,message,parts,periodic.and.ny>1)
    ok=io_agreement(len_trim(message)==0)
    if(.not.ok)goto 900
    call configuration_bytes(layout,serial_configuration,ios)
    ok=io_agreement(ios==0)
    if(ok)then
      serial_source=index;serial_source_path=path
      if(parallel_rank==0)write(*,'(A,I0,A)')'CHECKPOINT: importing serial V',index%version, &
        '; source has no saved primitives or last dt. Historical counter totals belong to rank 0.'
    end if
900 continue
    call finish_io(ok,status,'serial checkpoint configuration rejected: '//trim(path))
  end subroutine read_serial_checkpoint_configuration

  subroutine import_serial_checkpoint(path,layout,status)
    character(len=*), intent(in) :: path
    type(domain_layout), intent(in) :: layout
    integer, intent(out), optional :: status
    real*8, allocatable :: saved_u(:,:,:,:),saved_p(:,:,:,:)
    character, allocatable :: configuration(:)
    integer(int64) :: signature(2),offset
    real*8 :: work(neq)
    integer :: ios,unit,close_status,q,k,j,i
    logical :: ok,diagnostics_enabled
    call validate_call(path,layout,present(status),.true.)
    call agree_path(path,present(status),5)
    ok=serial_source_path==path.and.allocated(serial_configuration)
    ok=io_agreement(ok)
    if(.not.ok)goto 900
    call configuration_bytes(layout,configuration,ios)
    if(ios==0)then
      if(size(configuration)/=size(serial_configuration))then
        ios=1
      else if(any(configuration/=serial_configuration))then
        ios=1
      end if
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    call root_file_signature(path,signature,ios)
    ok=io_agreement(ios==0.and.all(signature==serial_source%signature))
    if(.not.ok)goto 900
    allocate(saved_u(nx,ny,nz,neq),saved_p(nx,ny,nz,neq),stat=ios)
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    open(newunit=unit,file=path,status='old',access='stream',form='unformatted',action='read',iostat=ios)
    if(ios==0)then
      do q=1,neq
        do k=1,nz
          do j=1,ny
            offset=(((int(q-1,int64)*serial_source%shape(3)+k+layout%offset(3)-1)*serial_source%shape(2)+ &
              j+layout%offset(2)-1)*serial_source%shape(1)+layout%offset(1))*8_int64
            read(unit,pos=serial_source%payload+offset,iostat=ios)saved_u(:,j,k,q)
            if(ios/=0)exit
          end do
          if(ios/=0)exit
        end do
        if(ios/=0)exit
      end do
      close(unit,iostat=close_status)
      if(ios==0)ios=close_status
    end if
    if(ios==0)then
      if(.not.all(ieee_is_finite(saved_u)))ios=1
    end if
    ok=io_agreement(ios==0)
    if(.not.ok)goto 900
    diagnostics_enabled=do_numerical_diagnostics;do_numerical_diagnostics=.false.
    do k=1,nz
      do j=1,ny
        do i=1,nx
          work=saved_u(i,j,k,:)
          call cons_to_prim(work,saved_p(i,j,k,:),i,j,k)
        end do
      end do
    end do
    do_numerical_diagnostics=diagnostics_enabled
    ok=io_agreement(physical_state(saved_p,saved_u))
    if(.not.ok)goto 900
    call root_file_signature(path,signature,ios)
    ok=io_agreement(ios==0.and.all(signature==serial_source%signature))
    if(.not.ok)goto 900
    up=saved_u;u=saved_u;p(1:nx,1:ny,1:nz,:)=saved_p
    integration_time=serial_source%time;n_steps=serial_source%step;dt=0.0d0
    perturbation_applied=serial_source%applied
    numerical_events=0_int64;numerical_delta_u=0.0d0
    if(parallel_rank==0)then
      numerical_events=serial_source%events;numerical_delta_u=serial_source%delta
    end if
    call fill_primitive_halos(p)
900 continue
    call finish_io(ok,status,'serial checkpoint import rejected: '//trim(path))
  end subroutine import_serial_checkpoint
end module parallel_checkpoint
