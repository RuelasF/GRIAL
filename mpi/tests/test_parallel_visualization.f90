program test_parallel_visualization
  use variables
  use parallel_runtime
  use parallel_visualization
  use domain_decomposition
  use initialization,only:allocate_and_grid,precalculate_metric_cache
  use metrics,only:set_metric_type
  use output,only:save_vtk_data
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use,intrinsic::iso_c_binding,only:c_int,c_long
  implicit none
  type,bind(C)::resource_limit
    integer(c_long)::soft,hard
  end type resource_limit
  interface
    integer(c_int)function get_file_limit(resource,limits)bind(C,name='getrlimit')
      import c_int,resource_limit
      integer(c_int),value::resource
      type(resource_limit),intent(out)::limits
    end function get_file_limit
    integer(c_int)function set_file_limit(resource,limits)bind(C,name='setrlimit')
      import c_int,resource_limit
      integer(c_int),value::resource
      type(resource_limit),intent(in)::limits
    end function set_file_limit
    function ignore_signal(signal,handler)bind(C,name='signal')result(previous)
      import c_int,c_long
      integer(c_int),value::signal
      integer(c_long),value::handler
      integer(c_long)::previous
    end function ignore_signal
  end interface
  type(resource_limit)::limits,limited
  type(domain_layout)::layout
  integer::grid(2),fixture,shape_global(3),i,j,k,id,ios,rc
  integer(c_long)::previous_handler
  real*8,allocatable::saved(:,:,:,:)
  character(len=160)::message
  character(len=512)::arg,path
  character(len=32)::mode
  call parallel_initialize()
  call get_command_argument(1,arg);read(arg,*)grid(1)
  call get_command_argument(2,arg);read(arg,*)grid(2)
  call get_command_argument(3,arg);read(arg,*)fixture
  call get_command_argument(4,mode)
  shape_global=[37,19,8]
  geom_type='Cartesian';metric_type='Minkowski';use_log_r=.false.;a_spin=0.0d0;bh_mass=1.0d0
  r_min=0.0d0;r_max=1.0d0;y_min=0.0d0;y_max=1.0d0;z_min=0.0d0;z_max=1.0d0
  select case(fixture)
  case(1);shape_global=[37,1,1]
  case(2);shape_global=[37,1,19]
  case(3);shape_global=[37,19,1]
  case(5);shape_global=[37,1,19];geom_type='Cylindrical';r_min=1.0d0;r_max=4.0d0
  case(6:14)
    geom_type='Spheroidal';metric_type='Kerr-Schild';a_spin=0.9d0
    r_min=1.2d0;r_max=40.0d0;y_min=0.0d0;y_max=pi;z_min=0.0d0;z_max=2.0d0*pi
    if(fixture==6.or.fixture==14)then
      geom_type='Spherical';metric_type='Eddington-Finkelstein';a_spin=0.0d0
    end if
    if(fixture>=8)use_log_r=.true.
    if(fixture==10)shape_global=[37,19,1]
    if(fixture==11)shape_global=[37,1,8]
    if(fixture==12)then
      y_min=0.2d0;y_max=pi-0.2d0
    end if
    if(fixture==13)z_max=pi
    if(fixture==14)shape_global=[37,1,1]
  case(15);shape_global=[1,1,1]
  end select
  vtk_mapping_id=VTK_MAP_PHYSICAL
  if(fixture==9)vtk_mapping_id=VTK_MAP_UNTWISTED
  nghost=3;case_id=4;case_name='Conv';is_mhd=.false.
  adb_idx=4.0d0/3.0d0;g1=adb_idx/(adb_idx-1.0d0)
  call set_metric_type()
  call build_domain(shape_global,nghost,parallel_rank,parallel_size,.false.,layout,message,grid)
  call parallel_require_collective(len_trim(message)==0,message)
  nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
  call allocate_and_grid(layout);call precalculate_metric_cache()
  p=ieee_value(0.0d0,ieee_quiet_nan)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        id=(i-1+x_global_offset)+nx_global*((j-1+y_global_offset)+ny_global*(k-1))
        p(i,j,k,:)=[1.0d0+0.001d0*dble(id),0.1d0+0.0001d0*dble(id), &
          0.01d0+1.0d-8*dble(id),0.002d0,0.003d0]
      end do
    end do
  end do
  saved=p(1:nx,1:ny,1:nz,:);n_steps=7;integration_time=0.125d0
  numerical_events=29_int64;numerical_delta_u=3.0d0
  call save_parallel_vtk('good',layout,ios)
  call require(ios==0,'valid generation')
  call require(all(saved==p(1:nx,1:ny,1:nz,:)),'writer must not change fields')
  call require(all(numerical_events==29_int64).and.all(numerical_delta_u==3.0d0),'writer must not change counters')
  if(parallel_size==1)then
    output_folder='.';scheme_name='good'
    call save_vtk_data('legacy.vtk')
  end if
  path='bad'
  select case(trim(mode))
  case('duplicate')
    call save_parallel_vtk('good',layout,ios)
    call require(ios/=0,'refuse existing generation')
  case('io')
    call save_parallel_vtk('missing/child',layout,ios)
    call require(ios/=0,'missing parent is a collective IO failure')
  case('limit')
    rc=0
    if(parallel_rank==parallel_size-1)then
      previous_handler=ignore_signal(25_c_int,1_c_long)
      rc=get_file_limit(1_c_int,limits);limited=limits;limited%soft=2000_c_long
      if(rc==0)rc=set_file_limit(1_c_int,limited)
    end if
    call require(rc==0,'establish real file limit after MPI startup')
    call save_parallel_vtk('truncated',layout,ios)
    call require(ios/=0,'truncated piece must not publish index')
    if(parallel_rank==parallel_size-1)rc=set_file_limit(1_c_int,limits)
    call require(rc==0,'restore file limit')
  case('path-disagree')
    if(parallel_rank==parallel_size-1)path='other'
    call save_parallel_vtk(trim(path),layout,ios)
    call parallel_abort('different paths accepted')
  case('epoch-disagree')
    if(parallel_rank==parallel_size-1)integration_time=0.25d0
    call save_parallel_vtk(trim(path),layout,ios)
    call parallel_abort('different epochs accepted')
  case('grid-disagree')
    if(parallel_rank==parallel_size-1)layout%offset(1)=layout%offset(1)+1
    call save_parallel_vtk(trim(path),layout,ios)
    call parallel_abort('invalid layout accepted')
  case('nan')
    if(parallel_rank==parallel_size-1)p(1,1,1,1)=ieee_value(0.0d0,ieee_quiet_nan)
    call save_parallel_vtk(trim(path),layout,ios)
    call parallel_abort('nonfinite field accepted')
  case('spin-disagree')
    if(parallel_rank==parallel_size-1)a_spin=0.1d0
    call save_parallel_vtk(trim(path),layout,ios)
    call parallel_abort('different geometry accepted')
  end select
  if(parallel_rank==0)write(*,'(A,I0,A,I0)')'VISUALIZATION PASSED fixture=',fixture,' ranks=',parallel_size
  call parallel_finalize()
contains
  subroutine require(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    call parallel_require_collective(ok,'visualization test: '//label)
  end subroutine require
end program test_parallel_visualization
