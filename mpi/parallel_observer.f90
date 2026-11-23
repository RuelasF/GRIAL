! C03: bind C01 fields, persistent D01 integrals and append-only diagnostic histories.
module parallel_observer
  use variables
  use parallel_runtime
  use domain_decomposition, only: domain_layout
  use parallel_diagnostics, only: diagnostic_history,export_diagnostic_history, &
    valid_diagnostic_history,initialize_parallel_diagnostics,diagnostics_active
  use parallel_checkpoint, only: save_parallel_checkpoint,load_parallel_checkpoint, &
    checkpoint_file_signature,checkpoint_crc64
  use output, only: distributed_output_cursor,export_distributed_output_cursor, &
    restore_distributed_output_cursor,valid_distributed_output_cursor, &
    write_distributed_diagnostics,write_distributed_observables
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: observer_calendar,initialize_observer_calendar,sample_parallel_observer
  public :: save_observed_checkpoint,load_observed_checkpoint
  type :: observer_calendar
    real*8 :: next_save=0.0d0,next_checkpoint=0.0d0
    integer(int64) :: save_index=0_int64
    integer :: diagnostic_step=-1,observable_step=-1
  end type observer_calendar
  character(len=16),parameter :: MAGIC='GRIAL_OBSERVER1',MARKER='GRIAL_OBS_COM1',FOOTER='GRIAL_OBS_END1'
  integer,parameter :: ABI(5)=[storage_size(0),storage_size(0_int64),storage_size(1.0d0), &
    storage_size(.true.),int(z'01020304')]
  character(len=40),parameter :: NAMES(7)=[character(len=40)::'mpi_global_diagnostics.dat', &
    'mpi_conservation_budget.dat','mpi_numerical_diagnostics.dat','mpi_ppi_modes.dat', &
    'mpi_GW_signal.dat','mpi_m_dot.dat','mpi_perturbation_event.dat']
contains

  logical function agreed(local) result(ok)
    logical,intent(in)::local
#ifdef USE_MPI
    integer::ios
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(local,ok,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ios)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_observer:agreed:1',collective_start)
    end block
#endif
    call parallel_check_mpi(ios,'observer IO agreement')
#else
    ok=local
#endif
  end function agreed

  subroutine finish(ok,status)
    logical,intent(in)::ok
    integer,intent(out),optional::status
    if(present(status))then
      status=merge(0,1,ok)
    else
      call parallel_require_collective(ok,'observed checkpoint rejected or not published')
    end if
  end subroutine finish

  subroutine agree_call(path,has_status,operation,calendar)
    character(len=*),intent(in)::path
    logical,intent(in)::has_status
    integer,intent(in)::operation
    type(observer_calendar),intent(in),optional::calendar
    character(len=512)::root_path
    integer(int64)::values(8),lower(8),upper(8)
#ifdef USE_MPI
    type(MPI_Datatype)::count_type
    integer::ios
#endif
    call parallel_require_active()
    call parallel_require_collective(len_trim(path)>0.and.len_trim(path)<=440.and.index(path,achar(0))==0, &
      'invalid observed checkpoint path')
    root_path=path;values=0_int64
    values(1:3)=[merge(1_int64,0_int64,has_status),int(operation,int64),merge(1_int64,0_int64,present(calendar))]
    if(present(calendar))values(4:8)=[transfer(calendar%next_save,0_int64), &
      transfer(calendar%next_checkpoint,0_int64),calendar%save_index, &
      int(calendar%diagnostic_step,int64),int(calendar%observable_step,int64)]
#ifdef USE_MPI
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(root_path,len(root_path),MPI_CHARACTER,0,MPI_COMM_WORLD,ios)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_observer:agree_call_path:84',collective_start)
    end block
#endif
    call parallel_check_mpi(ios,'observer path agreement')
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(0_int64)/8,count_type,ios)
    call parallel_check_mpi(ios,'observer integer datatype')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(values,lower,size(values),count_type,MPI_MIN,MPI_COMM_WORLD,ios)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_observer:agree_call:2',collective_start)
    end block
#endif
    call parallel_check_mpi(ios,'observer controls minimum')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(values,upper,size(values),count_type,MPI_MAX,MPI_COMM_WORLD,ios)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_observer:agree_call:3',collective_start)
    end block
#endif
    call parallel_check_mpi(ios,'observer controls maximum')
    call parallel_require_collective(all(lower==upper),'observer call or calendar differs between ranks')
#endif
    call parallel_require_collective(path==root_path,'observed checkpoint paths differ between ranks')
  end subroutine agree_call

  logical function valid_calendar(calendar,step) result(ok)
    type(observer_calendar),intent(in)::calendar
    integer,intent(in)::step
    ok=all(ieee_is_finite([calendar%next_save,calendar%next_checkpoint]))
    if(.not.ok)return
    ok=calendar%next_save>0.0d0.and.calendar%next_checkpoint>0.0d0.and.calendar%save_index>=0_int64.and. &
      calendar%diagnostic_step>=-1.and.calendar%diagnostic_step<=step.and. &
      calendar%observable_step>=-1.and.calendar%observable_step<=step
  end function valid_calendar

  subroutine initialize_observer_calendar(calendar)
    type(observer_calendar),intent(out)::calendar
    calendar=observer_calendar(next_save=save_interval,next_checkpoint=checkpoint_interval)
    call parallel_require_collective(valid_calendar(calendar,n_steps),'invalid initial observer calendar')
  end subroutine initialize_observer_calendar

  subroutine sample_parallel_observer(calendar,force)
    type(observer_calendar),intent(inout)::calendar
    logical,intent(in),optional::force
    logical::forced
    forced=.false.
    if(present(force))forced=force
    call agree_call('sample',present(force),merge(5,4,forced),calendar)
    call parallel_require_collective(valid_calendar(calendar,n_steps).and.diagnostic_stride>0.and.extraction_stride>0, &
      'invalid observer sampling calendar or strides')
    if((forced.or.mod(n_steps,diagnostic_stride)==0).and.n_steps>calendar%diagnostic_step)then
      call write_distributed_diagnostics(integration_time);calendar%diagnostic_step=n_steps
    end if
    if((forced.or.mod(n_steps,extraction_stride)==0).and.n_steps>calendar%observable_step)then
      call write_distributed_observables(integration_time);calendar%observable_step=n_steps
    end if
  end subroutine sample_parallel_observer

  ! Root checks both CRC/size and semantic tails. An older checkpoint cannot append
  ! to histories that advanced after it; rollback requires an explicit separate copy.
  subroutine history_signatures(cursor,time,step,applied,signatures,status,expected)
    type(distributed_output_cursor),intent(in)::cursor
    real*8,intent(in)::time
    integer,intent(in)::step
    logical,intent(in)::applied
    integer(int64),intent(out)::signatures(2,7)
    integer,intent(out)::status
    integer(int64),intent(in),optional::expected(2,7)
    character(len=512)::path,root
    character(len=8192)::line
    real*8::previous,value,last
    integer::i,unit,ios,close_status,row_step,last_step,rows
    logical::wanted,exists
    status=0;signatures=0_int64;root=trim(output_folder)//'/'//trim(scheme_name)
    do i=1,7
      path=trim(root)//'/'//trim(NAMES(i))
      wanted=len_trim(cursor%diagnostic_root)>0
      if(i>=4.and.i<=6)wanted=len_trim(cursor%observable_path(i-3))>0
      if(i==7)wanted=applied
      inquire(file=trim(path),exist=exists,iostat=ios)
      if(ios/=0.or.(exists.neqv.wanted))then
        status=1;return
      end if
      if(.not.exists)then
        signatures(:,i)=[-1_int64,0_int64]
        if(present(expected))then
          if(any(signatures(:,i)/=expected(:,i)))then
            status=1;return
          end if
        end if
        cycle
      end if
      call checkpoint_file_signature(trim(path),signatures(:,i),status)
      if(status/=0)return
      if(present(expected))then
        if(any(signatures(:,i)/=expected(:,i)))then
          status=1;return
        end if
      end if
      open(newunit=unit,file=trim(path),status='old',action='read',iostat=status)
      if(status/=0)return
      previous=-1.0d0;last_step=-1;rows=0;last=-1.0d0
      do
        read(unit,'(A)',iostat=ios)line
        if(ios<0)exit
        if(ios/=0)then
          status=1;exit
        end if
        if(len_trim(line)==0)cycle
        if(line(1:1)=='#')cycle
        row_step=-1
        if(i<=3.or.i==7)then
          read(line,*,iostat=ios)value,row_step
        else
          read(line,*,iostat=ios)value
        end if
        if(ios/=0)then
          status=1;exit
        end if
        if(.not.ieee_is_finite(value))then
          status=1;exit
        end if
        if(value<=previous.or.value>time.or.value<0.0d0)then
          status=1;exit
        end if
        if(i<=3.or.i==7)then
          if(row_step<=last_step.or.row_step>step)then
            status=1;exit
          end if
        end if
        rows=rows+1;previous=value;last=value;last_step=row_step
      end do
      close(unit,iostat=close_status)
      if(close_status/=0.or.rows<1)status=1
      if(i<=3)then
        if(last/=cursor%diagnostic_time.or.last_step/=cursor%diagnostic_step)status=1
      else if(i<=6)then
        if(last/=cursor%observable_time(i-3))status=1
      else
        if(rows/=1)status=1
      end if
      if(status/=0)return
    end do
  end subroutine history_signatures

  subroutine write_bytes(path,bytes,status)
    character(len=*),intent(in)::path
    character,intent(in)::bytes(:)
    integer,intent(out)::status
    integer::unit,close_status
    integer(int64)::actual(2)
    open(newunit=unit,file=path,status='new',access='stream',form='unformatted',action='write',iostat=status)
    if(status/=0)return
    write(unit,iostat=status)bytes
    close(unit,iostat=close_status)
    if(status==0)status=close_status
    if(status==0)call checkpoint_file_signature(path,actual,status)
    if(status==0)then
      if(any(actual/=[int(size(bytes),int64),checkpoint_crc64(bytes)]))status=1
    end if
  end subroutine write_bytes

  subroutine scratch_bytes(unit,bytes,status)
    integer,intent(in)::unit
    character,allocatable,intent(out)::bytes(:)
    integer,intent(out)::status
    integer::length
    inquire(unit=unit,size=length,iostat=status)
    if(status==0)allocate(bytes(length),stat=status)
    if(status==0)read(unit,pos=1,iostat=status)bytes
  end subroutine scratch_bytes

  subroutine save_observed_checkpoint(path,layout,calendar,status)
    character(len=*),intent(in)::path
    type(domain_layout),intent(in)::layout
    type(observer_calendar),intent(in)::calendar
    integer,intent(out),optional::status
    type(diagnostic_history)::history
    type(distributed_output_cursor)::cursor
    real*8,allocatable::packets(:,:)
    character,allocatable::bytes(:),marker_bytes(:)
    integer(int64)::physical(2),sidecar(2),signatures(2,7)
    integer::ios,unit,close_status
    logical::ok
#ifdef USE_MPI
    integer::mpi_result
#endif
    call agree_call(path,present(status),1,calendar)
    call parallel_require_collective(valid_calendar(calendar,n_steps),'invalid checkpoint observer calendar')
    call export_diagnostic_history(history);call export_distributed_output_cursor(cursor)
    ios=0
    if(parallel_rank==0)then
      if(.not.valid_distributed_output_cursor(cursor,integration_time,n_steps))ios=1
      if(calendar%diagnostic_step/=cursor%diagnostic_step)ios=1
      if(ios==0)call history_signatures(cursor,integration_time,n_steps,perturbation_applied,signatures,ios)
    end if
    ok=agreed(ios==0)
    if(.not.ok)goto 900
    allocate(packets(10*neq,parallel_size),stat=ios)
    ok=agreed(ios==0)
    if(.not.ok)goto 900
#ifdef USE_MPI
    call MPI_Gather(history%sums,10*neq,MPI_DOUBLE_PRECISION,packets,10*neq,MPI_DOUBLE_PRECISION, &
      0,MPI_COMM_WORLD,mpi_result)
    call parallel_check_mpi(mpi_result,'observer integral histories')
#else
    packets(:,1)=history%sums
#endif
    call save_parallel_checkpoint(path,layout,ios)
    ok=ios==0
    if(.not.ok)goto 900
    if(parallel_rank==0)then
      call checkpoint_file_signature(trim(path)//'/COMMITTED',physical,ios)
      if(ios==0)open(newunit=unit,status='scratch',access='stream',form='unformatted',iostat=ios)
      if(ios==0)then
        write(unit,iostat=ios)MAGIC,ABI,layout%global_shape,layout%process_shape,parallel_size, &
          n_steps,integration_time,dt,perturbation_applied,calendar,cursor, &
          history%steps,history%extra_faces,physical,signatures,packets,FOOTER
        if(ios==0)call scratch_bytes(unit,bytes,ios)
        close(unit,iostat=close_status)
        if(ios==0)ios=close_status
      end if
      if(ios==0)call write_bytes(trim(path)//'/observer.bin',bytes,ios)
      if(ios==0)then
        sidecar=[int(size(bytes),int64),checkpoint_crc64(bytes)]
        open(newunit=unit,status='scratch',access='stream',form='unformatted',iostat=ios)
        if(ios==0)then
          write(unit,iostat=ios)MARKER,physical,sidecar,FOOTER
          if(ios==0)call scratch_bytes(unit,marker_bytes,ios)
          close(unit,iostat=close_status)
          if(ios==0)ios=close_status
        end if
      end if
      if(ios==0)call write_bytes(trim(path)//'/OBSERVED.tmp',marker_bytes,ios)
      if(ios==0)call rename(trim(path)//'/OBSERVED.tmp',trim(path)//'/OBSERVED',ios)
    end if
    ok=agreed(ios==0)
900 continue
    call finish(ok,status)
  end subroutine save_observed_checkpoint

  subroutine read_observer(path,layout,calendar,cursor,history,status)
    character(len=*),intent(in)::path
    type(domain_layout),intent(in)::layout
    type(observer_calendar),intent(out)::calendar
    type(distributed_output_cursor),intent(out)::cursor
    type(diagnostic_history),intent(out)::history
    integer,intent(out)::status
    integer::unit,close_status,length,abi_tag(5),shape_global(3),parts(2),ranks,step,header,i
    integer(int64)::physical(2),sidecar(2),actual(2),linked(2),signatures(2,7),checked(2,7),file_size
    character(len=16)::label,ending
    real*8::time,step_size
    real*8,allocatable::packets(:,:)
    logical::applied
    history=diagnostic_history();calendar=observer_calendar();cursor=distributed_output_cursor()
    open(newunit=unit,file=trim(path)//'/OBSERVED',status='old',access='stream', &
      form='unformatted',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)label,physical,sidecar,ending
    inquire(iolength=length)MARKER,physical,sidecar,FOOTER
    inquire(unit=unit,size=file_size)
    close(unit,iostat=close_status)
    if(status/=0)return
    if(close_status/=0.or.file_size/=length.or.label/=MARKER.or.ending/=FOOTER)then
      status=1;return
    end if
    call checkpoint_file_signature(trim(path)//'/COMMITTED',actual,status)
    if(status/=0)return
    if(any(actual/=physical))then
      status=1;return
    end if
    call checkpoint_file_signature(trim(path)//'/observer.bin',actual,status)
    if(status/=0)return
    if(any(actual/=sidecar))then
      status=1;return
    end if
    open(newunit=unit,file=trim(path)//'/observer.bin',status='old',access='stream', &
      form='unformatted',action='read',iostat=status)
    if(status/=0)return
    read(unit,iostat=status)label,abi_tag,shape_global,parts,ranks,step,time,step_size,applied,calendar,cursor, &
      history%steps,history%extra_faces,linked,signatures
    inquire(iolength=header)MAGIC,ABI,shape_global,parts,ranks,step,time,step_size,applied,calendar,cursor, &
      history%steps,history%extra_faces,linked,signatures
    if(status==0)then
      if(label/=MAGIC.or.any(abi_tag/=ABI).or.any(shape_global/=layout%global_shape).or. &
        ranks<1.or.any(parts<1).or.step<0.or.any(linked/=physical))status=1
      if(int(parts(1),int64)*int(parts(2),int64)/=int(ranks,int64))status=1
      if(.not.all(ieee_is_finite([time,step_size])))status=1
    end if
    if(status==0)then
      if(time/=integration_time.or.step/=n_steps.or.step_size/=dt.or.(applied.neqv.perturbation_applied))status=1
      if(.not.valid_calendar(calendar,step).or..not.valid_diagnostic_history(history,step))status=1
      if(.not.valid_distributed_output_cursor(cursor,time,step))status=1
      if(calendar%diagnostic_step/=cursor%diagnostic_step)status=1
      if(sidecar(1)/=int(header,int64)+int(ranks,int64)*int(10*neq*8,int64)+16_int64)status=1
    end if
    if(status==0)allocate(packets(10*neq,ranks),stat=status)
    if(status==0)read(unit,iostat=status)packets,ending
    close(unit,iostat=close_status)
    if(status==0)status=close_status
    if(status/=0)return
    if(ending/=FOOTER.or..not.all(ieee_is_finite(packets)))then
      status=1;return
    end if
    if(all(parts==layout%process_shape))then
      history%sums=packets(:,parallel_rank+1)
    else
      history%sums=0.0d0
      ! Validate on every rank so a bad sum cannot become a post-field-load failure.
      do i=1,ranks
        if(any(abs(history%sums)>huge(1.0d0)-abs(packets(:,i))))then
          status=1;return
        end if
        history%sums=history%sums+packets(:,i)
      end do
      if(parallel_rank/=0)history%sums=0.0d0
    end if
    if(parallel_rank==0)call history_signatures(cursor,time,step,applied,checked,status,signatures)
  end subroutine read_observer

  subroutine load_observed_checkpoint(path,layout,calendar,status)
    character(len=*),intent(in)::path
    type(domain_layout),intent(in)::layout
    type(observer_calendar),intent(inout)::calendar
    integer,intent(out),optional::status
    type(observer_calendar)::saved_calendar
    type(distributed_output_cursor)::cursor
    type(diagnostic_history)::history
    integer::ios
    logical::ok
    call agree_call(path,present(status),2)
    call parallel_require_collective(.not.diagnostics_active,'finalize D01 before observer recovery')
    call read_observer(path,layout,saved_calendar,cursor,history,ios)
    ok=agreed(ios==0)
    if(.not.ok)goto 900
    call load_parallel_checkpoint(path,layout,ios,redistribute=.true.)
    ok=ios==0
    if(.not.ok)goto 900
    call initialize_parallel_diagnostics(layout,history)
    call restore_distributed_output_cursor(cursor);calendar=saved_calendar
900 continue
    call finish(ok,status)
  end subroutine load_observed_checkpoint
end module parallel_observer
