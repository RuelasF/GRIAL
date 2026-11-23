! C04: local physical cells in VTU pieces; root publishes only small XML indices.
module parallel_visualization
  use variables
  use parallel_runtime
  use domain_decomposition, only: domain_layout,validate_domain
  use parallel_diagnostics, only: validate_checkpoint_boundary
  use output, only: vtk_polar_position
  use metrics, only: polar_angle,polar_jacobian
  use, intrinsic :: iso_fortran_env, only: int8
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  private
  public :: save_parallel_vtk
  character(len=24), parameter :: names(7)=[character(len=24):: &
    'Density','Pressure','Velocity_X','Velocity_Y','Velocity_Z','Lapse_alpha','Shift_beta_r']
  interface
    integer(c_int) function posix_mkdir(path,mode) bind(C,name='mkdir')
      import c_char,c_int
      character(kind=c_char),intent(in)::path(*)
      integer(c_int),value::mode
    end function posix_mkdir
  end interface
contains
  subroutine agree_status(status)
    integer,intent(inout)::status
#ifdef USE_MPI
    integer::rc
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Allreduce(MPI_IN_PLACE,status,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,rc)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_visualization:publish_status:36',collective_start)
    end block
#endif
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_visualization:agree_status:1',collective_start)
    end block
#endif
    call parallel_check_mpi(rc,'VTK writer status')
#endif
  end subroutine agree_status

  function piece_name(rank) result(name)
    integer,intent(in)::rank
    character(len=32)::name
    write(name,'(A,I10.10,A)')'piece_',rank,'.vtu'
  end function piece_name

  subroutine validate_call(path,layout)
    character(len=*),intent(in)::path
    type(domain_layout),intent(in)::layout
    character(len=512)::root_path
    character(len=len(geom_type))::root_geometry
    character(len=len(metric_type))::root_metric
    character(len=160)::message
    integer(int64)::values(27),root_values(27)
    logical::valid
#ifdef USE_MPI
    integer::rc
    type(MPI_Datatype)::count_type
#endif
    call parallel_require_active()
    call validate_domain(layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    valid=len_trim(path)>0.and.len_trim(path)<=440.and.index(path,c_null_char)==0
    call parallel_require_collective(valid,'invalid VTK generation path')
    root_path=path
    root_geometry=geom_type;root_metric=metric_type
#ifdef USE_MPI
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(root_path,len(root_path),MPI_CHARACTER,0,MPI_COMM_WORLD,rc)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_visualization:agree_frame_path:72',collective_start)
    end block
#endif
    call parallel_check_mpi(rc,'VTK path agreement')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(root_geometry,len(root_geometry),MPI_CHARACTER,0,MPI_COMM_WORLD,rc)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_visualization:agree_frame_geometry:74',collective_start)
    end block
#endif
    call parallel_check_mpi(rc,'VTK geometry agreement')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(root_metric,len(root_metric),MPI_CHARACTER,0,MPI_COMM_WORLD,rc)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_visualization:agree_frame_metric:76',collective_start)
    end block
#endif
    call parallel_check_mpi(rc,'VTK metric agreement')
#endif
    call parallel_require_collective(path==root_path,'VTK path differs between ranks')
    call parallel_require_collective(geom_type==root_geometry.and.metric_type==root_metric, &
      'VTK geometry or metric differs between ranks')
    valid=layout%rank==parallel_rank.and.layout%ranks==parallel_size.and.layout%halo==nghost.and. &
      all(layout%local_shape==[nx,ny,nz]).and.all(layout%global_shape==[nx_global,ny_global,nz_global]).and. &
      all(layout%offset==[x_global_offset,y_global_offset,0])
    call parallel_require_collective(valid,'VTK layout does not match allocated grid')
    valid=allocated(p).and.allocated(alpha_c).and.allocated(beta_c).and. &
      allocated(x_face).and.allocated(y_face).and.allocated(z_face)
    call parallel_require_collective(valid,'VTK fields or coordinates are not allocated')
    valid=all(shape(p)==[nx+2*nghost+1,ny+2*nghost+1,nz+2*nghost+1,neq]).and. &
      all(shape(alpha_c)==[nx,ny]).and.all(shape(beta_c)==[3,nx,ny])
    call parallel_require_collective(valid,'VTK field shapes do not match layout')
    call parallel_require_collective(.not.accelerator_enabled.and..not.is_mhd.and.neq==5, &
      'VTK pieces require CPU hydrodynamics')
    valid=all(ieee_is_finite(p(1:nx,1:ny,1:nz,:))).and.all(ieee_is_finite(alpha_c)).and. &
      all(ieee_is_finite(beta_c)).and.all(ieee_is_finite(x_face)).and. &
      all(ieee_is_finite(y_face)).and.all(ieee_is_finite(z_face))
    valid=valid.and.ieee_is_finite(integration_time).and.integration_time>=0.0d0.and.n_steps>=0
    call parallel_require_collective(valid,'nonfinite VTK state or epoch')
    values=[int(layout%global_shape,int64),int(layout%process_shape,int64),int(n_steps,int64), &
      transfer(integration_time,0_int64),int(vtk_mapping_id,int64),merge(1_int64,0_int64,use_log_r), &
      transfer(a_spin,0_int64),transfer(bh_mass,0_int64), &
      transfer(r_min,0_int64),transfer(r_max,0_int64),transfer(y_min,0_int64),transfer(y_max,0_int64), &
      transfer(z_min,0_int64),transfer(z_max,0_int64), &
      int(iachar(geom_type(1:1)),int64),int(iachar(metric_type(1:1)),int64), &
      int(nghost,int64),int(neq,int64),int(case_id,int64), &
      transfer(dx,0_int64),transfer(dy,0_int64),transfer(dz,0_int64),int(parallel_size,int64),transfer(polar_mesh_h,0_int64)]
    root_values=values
#ifdef USE_MPI
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER,storage_size(0_int64)/8,count_type,rc)
    call parallel_check_mpi(rc,'VTK int64 type')
#ifdef GRIAL_PROFILE
    block
      real*8 :: collective_start
      collective_start=profile_clock()
#endif
    call MPI_Bcast(root_values,size(values),count_type,0,MPI_COMM_WORLD,rc)
#ifdef GRIAL_PROFILE
      call profile_label_add('parallel_visualization:agree_frame_controls:111',collective_start)
    end block
#endif
    call parallel_check_mpi(rc,'VTK epoch agreement')
#endif
    call parallel_require_collective(all(values==root_values),'VTK epoch or configuration differs between ranks')
    call validate_checkpoint_boundary()
  end subroutine validate_call

  subroutine save_parallel_vtk(path,layout,status)
    character(len=*),intent(in)::path
    type(domain_layout),intent(in)::layout
    integer,intent(out)::status
    call validate_call(path,layout)
    status=0
    if(parallel_rank==0)then
      ! A failed generation is retained; retry requires another directory.
      if(posix_mkdir(trim(path)//c_null_char,int(o'700',c_int))/=0)status=1
    end if
    call agree_status(status)
    if(status/=0)return
    call write_piece(trim(path)//'/'//trim(piece_name(parallel_rank)),layout,status)
    call agree_status(status)
    if(status/=0)return
    if(parallel_rank==0)call write_indices(path,status)
    call agree_status(status)
  end subroutine save_parallel_vtk

  subroutine make_mesh(points,connectivity,offsets,types)
    real*8,allocatable,intent(out)::points(:,:)
    integer(int64),allocatable,intent(out)::connectivity(:),offsets(:)
    integer(int8),allocatable,intent(out)::types(:)
    integer(int64),allocatable::ids(:,:,:)
    integer::i,j,k,jp,kp,npoints,ncells,nconnect,cell,vertices(8),kind
    real*8::r,theta,phi,spin,tolerance
    logical::polar,closed_phi,north,south
    tolerance=128.0d0*epsilon(1.0d0)
    polar=trim(geom_type)=='Spherical'.or.trim(geom_type)=='Spheroidal'
    jp=0;kp=0
    if(ny_global>1)jp=ny
    if(nz_global>1)kp=nz
    closed_phi=polar.and.jp>0.and.kp>0.and.abs(z_max-z_min-2.0d0*pi)<tolerance*2.0d0*pi
    north=polar.and.jp>0.and.kp>0.and.abs(y_face(0))<tolerance
    south=polar.and.jp>0.and.kp>0.and.abs(y_face(ny)-pi)<tolerance*pi
    spin=0.0d0
    if(trim(metric_type)=='Kerr-Schild')spin=a_spin
    allocate(ids(0:nx,0:jp,0:kp),points(3,(nx+1)*(jp+1)*(kp+1)))
    npoints=0
    do k=0,kp
      phi=z(1)
      if(kp>0)phi=z_face(k)
      do j=0,jp
        theta=y(1)
        if(jp>0)theta=y_face(j)
        if(polar)theta=polar_angle(theta)
        do i=0,nx
          if(closed_phi.and.k==kp)then
            ids(i,j,k)=ids(i,j,0);cycle
          end if
          if(k>0.and.((north.and.j==0).or.(south.and.j==jp)))then
            ids(i,j,k)=ids(i,j,0);cycle
          end if
          r=x_face(i)
          if(use_log_r)r=exp(r)
          ids(i,j,k)=int(npoints,int64);npoints=npoints+1
          select case(trim(geom_type))
          case('Spherical','Spheroidal')
            call vtk_polar_position(r,theta,phi,spin,points(:,npoints))
          case('Cylindrical')
            points(:,npoints)=[r*cos(theta),r*sin(theta),phi]
          case default
            points(:,npoints)=[r,theta,phi]
          end select
        end do
      end do
    end do
    points=points(:,:npoints)
    ncells=nx*ny*nz
    allocate(connectivity(8*ncells),offsets(ncells),types(ncells))
    nconnect=0;cell=0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          cell=cell+1
          if(jp==0.and.kp==0)then
            vertices(1:2)=[ids(i-1,0,0),ids(i,0,0)];kind=3
            call add_cell(2)
          else if(jp==0)then
            vertices(1:4)=[ids(i-1,0,k-1),ids(i,0,k-1),ids(i,0,k),ids(i-1,0,k)];kind=9
            call add_cell(4)
          else if(kp==0)then
            vertices(1:4)=[ids(i-1,j-1,0),ids(i,j-1,0),ids(i,j,0),ids(i-1,j,0)];kind=9
            call add_cell(4)
          else if(north.and.j==1)then
            vertices(1:6)=[ids(i-1,0,k-1),ids(i-1,1,k-1),ids(i-1,1,k), &
              ids(i,0,k-1),ids(i,1,k-1),ids(i,1,k)];kind=13
            call add_cell(6)
          else if(south.and.j==ny)then
            vertices(1:6)=[ids(i-1,ny,k-1),ids(i-1,ny-1,k),ids(i-1,ny-1,k-1), &
              ids(i,ny,k-1),ids(i,ny-1,k),ids(i,ny-1,k-1)];kind=13
            call add_cell(6)
          else
            vertices=[ids(i-1,j-1,k-1),ids(i,j-1,k-1),ids(i,j,k-1),ids(i-1,j,k-1), &
              ids(i-1,j-1,k),ids(i,j-1,k),ids(i,j,k),ids(i-1,j,k)];kind=12
            call add_cell(8)
          end if
        end do
      end do
    end do
    connectivity=connectivity(:nconnect)
  contains
    subroutine add_cell(count)
      integer,intent(in)::count
      connectivity(nconnect+1:nconnect+count)=int(vertices(:count),int64)
      nconnect=nconnect+count;offsets(cell)=int(nconnect,int64);types(cell)=int(kind,int8)
    end subroutine add_cell
  end subroutine make_mesh

  subroutine write_piece(path,layout,status)
    character(len=*),intent(in)::path
    type(domain_layout),intent(in)::layout
    integer,intent(out)::status
    real*8,allocatable::points(:,:),field(:)
    integer(int64),allocatable::connectivity(:),offsets(:),ids(:)
    integer(int8),allocatable::types(:)
    integer::unit,rc,i,j,k,c,q,ncells
    integer(int64)::next_offset,expected,actual,bytes
    character(len=256)::line
    status=0;expected=-2_int64;actual=-1_int64
    ! Arrays and VTK IDs are local; reject index overflow before allocation.
    if(int(nx+1,int64)*int(ny+1,int64)*int(nz+1,int64)>huge(1)/3.or. &
       int(nx,int64)*int(ny,int64)*int(nz,int64)>huge(1)/8)then
      status=1;return
    end if
    call make_mesh(points,connectivity,offsets,types)
    ncells=nx*ny*nz;allocate(field(ncells),ids(ncells))
    open(newunit=unit,file=trim(path)//'.tmp',status='new',access='stream',form='unformatted', &
      action='write',convert='big_endian',iostat=status)
    if(status/=0)then
      status=1;return
    end if
    call xml('<?xml version="1.0"?>')
    call xml('<VTKFile type="UnstructuredGrid" version="1.0" byte_order="BigEndian" header_type="UInt64">')
    call xml('<UnstructuredGrid><FieldData>')
    write(line,'(ES24.16E3)')integration_time
    call metadata('Float64','TimeValue',trim(line))
    write(line,'(I0)')n_steps;call metadata('Int64','Step',trim(line))
    write(line,'(3(I0,1X))')layout%global_shape;call metadata('Int64','GlobalShape',trim(line),3)
    write(line,'(3(I0,1X))')layout%local_shape;call metadata('Int64','LocalShape',trim(line),3)
    write(line,'(3(I0,1X))')layout%offset;call metadata('Int64','GlobalOffset',trim(line),3)
    write(line,'(2(I0,1X))')layout%process_shape;call metadata('Int64','ProcessGrid',trim(line),2)
    write(line,'(I0)')parallel_rank;call metadata('Int64','Rank',trim(line))
    write(line,'(I0)')vtk_mapping_id;call metadata('Int64','MappingId',trim(line))
    write(line,'(I0)')merge(1,0,use_log_r);call metadata('Int64','LogRadius',trim(line))
    write(line,'(ES24.16E3)')a_spin;call metadata('Float64','Spin',trim(line))
    if(polar_mesh_h/=1.0d0)then
      write(line,'(ES24.16E3)')polar_mesh_h;call metadata('Float64','PolarMeshH',trim(line))
    end if
    call xml('</FieldData>')
    write(line,'(A,I0,A,I0,A)')'<Piece NumberOfPoints="',size(points,2),'" NumberOfCells="',ncells,'">'
    call xml(trim(line));call xml('<PointData/><CellData Scalars="Density">')
    next_offset=0_int64;bytes=8_int64*int(ncells,int64)
    do q=1,7
      call array_header('Float64',trim(names(q)),1,bytes)
    end do
    call array_header('Int64','GlobalCellId',1,bytes)
    call array_header('Int64','RankId',1,bytes)
    call xml('</CellData><Points>')
    call array_header('Float64','Points',3,8_int64*int(size(points),int64))
    call xml('</Points><Cells>')
    call array_header('Int64','connectivity',1,8_int64*int(size(connectivity),int64))
    call array_header('Int64','offsets',1,bytes)
    call array_header('UInt8','types',1,int(ncells,int64))
    call xml('</Cells></Piece></UnstructuredGrid>')
    call xml('<AppendedData encoding="raw">',.false.)
    if(status==0)write(unit,iostat=status)'_'
    do q=1,7
      c=0
      do k=1,nz
        do j=1,ny
          do i=1,nx
            c=c+1
            select case(q)
            case(1:5);field(c)=p(i,j,k,q)
            case(6);field(c)=alpha_c(i,j)
            case(7);field(c)=beta_c(1,i,j)
            end select
            if(use_log_r.and.(q==eq_vx.or.q==7))field(c)=field(c)*exp(x(i))
            if(polar_mesh_h/=1.0d0.and.q==eq_vy)field(c)=field(c)*polar_jacobian(y(j))
            ids(c)=int(i-1+x_global_offset,int64)+int(nx_global,int64)* &
              (int(j-1+y_global_offset,int64)+int(ny_global,int64)*int(k-1,int64))
          end do
        end do
      end do
      if(status==0)write(unit,iostat=status)bytes,field
    end do
    if(status==0)write(unit,iostat=status)bytes,ids
    ids=int(parallel_rank,int64)
    if(status==0)write(unit,iostat=status)bytes,ids
    if(status==0)write(unit,iostat=status)8_int64*int(size(points),int64),points
    if(status==0)write(unit,iostat=status)8_int64*int(size(connectivity),int64),connectivity
    if(status==0)write(unit,iostat=status)bytes,offsets
    if(status==0)write(unit,iostat=status)int(ncells,int64),types
    call xml('</AppendedData></VTKFile>')
    inquire(unit=unit,pos=expected,iostat=rc)
    if(rc/=0)status=1
    close(unit,iostat=rc)
    if(rc/=0)status=1
    if(status==0)then
      inquire(file=trim(path)//'.tmp',size=actual,iostat=rc)
      if(rc/=0.or.actual/=expected-1_int64)status=1
    end if
    if(status==0)call rename(trim(path)//'.tmp',trim(path),status)
    if(status/=0)status=1
  contains
    subroutine xml(text,newline)
      character(len=*),intent(in)::text
      logical,intent(in),optional::newline
      if(status/=0)return
      if(present(newline))then
        if(.not.newline)then
          write(unit,iostat=status)text;return
        end if
      end if
      write(unit,iostat=status)text//achar(10)
    end subroutine xml
    subroutine metadata(kind,name,value,components)
      character(len=*),intent(in)::kind,name,value
      integer,intent(in),optional::components
      character(len=12)::count
      count='1'
      if(present(components))write(count,'(I0)')components
      call xml('<DataArray type="'//kind//'" Name="'//name//'" NumberOfTuples="1" NumberOfComponents="'// &
        trim(count)//'" format="ascii">'//value//'</DataArray>')
    end subroutine metadata
    subroutine array_header(kind,name,components,length)
      character(len=*),intent(in)::kind,name
      integer,intent(in)::components
      integer(int64),intent(in)::length
      character(len=32)::a,b
      write(a,'(I0)')components;write(b,'(I0)')next_offset
      call xml('<DataArray type="'//kind//'" Name="'//name//'" NumberOfComponents="'//trim(a)// &
        '" format="appended" offset="'//trim(b)//'"/>')
      next_offset=next_offset+8_int64+length
    end subroutine array_header
  end subroutine write_piece

  subroutine write_indices(path,status)
    character(len=*),intent(in)::path
    integer,intent(out)::status
    integer::unit,rc,rank,q
    integer(int64)::expected,actual
    character(len=32)::time
    status=0;expected=-2_int64;actual=-1_int64
    open(newunit=unit,file=trim(path)//'/frame.pvtu.tmp',status='new',access='stream', &
      form='unformatted',action='write',iostat=status)
    if(status/=0)then
      status=1;return
    end if
    call line('<?xml version="1.0"?>')
    call line('<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="BigEndian">')
    call line('<PUnstructuredGrid GhostLevel="0"><PPointData/><PCellData Scalars="Density">')
    do q=1,7
      call line('<PDataArray type="Float64" Name="'//trim(names(q))//'" NumberOfComponents="1"/>')
    end do
    call line('<PDataArray type="Int64" Name="GlobalCellId" NumberOfComponents="1"/>')
    call line('<PDataArray type="Int64" Name="RankId" NumberOfComponents="1"/>')
    call line('</PCellData><PPoints><PDataArray type="Float64" NumberOfComponents="3"/></PPoints>')
    do rank=0,parallel_size-1
      call line('<Piece Source="'//trim(piece_name(rank))//'"/>')
    end do
    call line('</PUnstructuredGrid></VTKFile>')
    call close_index('frame.pvtu.tmp')
    if(status/=0)return
    open(newunit=unit,file=trim(path)//'/frame.pvd',status='new',access='stream', &
      form='unformatted',action='write',iostat=status)
    if(status/=0)then
      status=1;return
    end if
    write(time,'(ES24.16E3)')integration_time
    call line('<?xml version="1.0"?>')
    call line('<VTKFile type="Collection" version="0.1" byte_order="BigEndian"><Collection>')
    call line('<DataSet timestep="'//trim(adjustl(time))//'" group="" part="0" file="frame.pvtu"/>')
    call line('</Collection></VTKFile>')
    call close_index('frame.pvd')
    ! The PVTU index is the final publication marker for this immutable frame.
    if(status==0)call rename(trim(path)//'/frame.pvtu.tmp',trim(path)//'/frame.pvtu',status)
    if(status/=0)status=1
  contains
    subroutine line(text)
      character(len=*),intent(in)::text
      if(status==0)write(unit,iostat=status)text//achar(10)
    end subroutine line
    subroutine close_index(name)
      character(len=*),intent(in)::name
      inquire(unit=unit,pos=expected,iostat=rc)
      if(rc/=0)status=1
      close(unit,iostat=rc)
      if(rc/=0)status=1
      if(status/=0)return
      inquire(file=trim(path)//'/'//name,size=actual,iostat=rc)
      if(rc/=0.or.actual/=expected-1_int64)status=1
    end subroutine close_index
  end subroutine write_indices
end module parallel_visualization
