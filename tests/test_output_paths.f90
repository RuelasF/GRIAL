program test_output_paths
  use variables
  use output
  implicit none

  integer :: i,unit_file,time_index
  real*8 :: stored_time
  character(len=120) :: title
  logical :: file_exists
  character(len=*), parameter :: test_directory = &
    '/tmp/grhd_output_path_test/123456789012345678901234567890123456789012345678901234567890'
  character(len=*), parameter :: long_prefix = &
    '12345678901234567890123456789012345678901234567890'
  character(len=*), parameter :: expected_file = test_directory // &
    '/weno5_superbee_hlle/' // long_prefix // '_weno5_superbee_hlle_0000.vtk'

  nghost = 1
  nx = 2
  ny = 1
  nz = 2
  use_log_r = .false.
  geom_type = 'Cartesian'
  metric_type = 'Minkowski'
  a_spin = 0.0d0
  save_interval = 1.0d0
  integration_time = 0.125d0
  output_prefix = long_prefix
  output_folder = test_directory
  scheme_name = 'weno5_superbee_hlle'

  allocate(x(-nghost:nx+nghost), y(-nghost:ny+nghost), z(-nghost:nz+nghost))
  allocate(x_face(0:nx), y_face(0:ny), z_face(0:nz))
  allocate(p(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq))
  allocate(alpha_c(nx,ny,nz), beta_c(3,nx,ny,nz))
  allocate(var_names(neq))

  do i = -nghost, nx+nghost
    x(i) = dble(i) - 0.5d0
  end do
  do i = -nghost, ny+nghost
    y(i) = dble(i) - 0.5d0
  end do
  do i = -nghost, nz+nghost
    z(i) = dble(i) - 0.5d0
  end do
  x_face = [0.0d0, 1.0d0, 2.0d0]
  y_face = [0.0d0, 1.0d0]
  z_face = [0.0d0, 1.0d0, 2.0d0]
  p = 0.0d0
  p(:,:,:,eq_de) = 1.0d0
  p(:,:,:,eq_pr) = 0.1d0
  alpha_c = 1.0d0
  beta_c = 0.0d0
  var_names(eq_de) = 'Density'
  var_names(eq_pr) = 'Pressure'
  var_names(eq_vx) = 'Velocity_X'
  var_names(eq_vy) = 'Velocity_Y'
  var_names(eq_vz) = 'Velocity_Z'

  call execute_command_line('mkdir -p '//test_directory//'/'//trim(scheme_name))
  call save_vtk_at_time(0.0d0)
  inquire(file=expected_file, exist=file_exists)
  if (.not. file_exists) then
    write(*,'(A)') 'OUTPUT PATH TEST FAILED: long VTK name was not created'
    error stop 1
  end if
  open(newunit=unit_file,file=expected_file,status='old',access='stream',form='formatted')
  read(unit_file,'(A)') title
  read(unit_file,'(A)') title
  close(unit_file)
  time_index=index(title,'t=')
  if (time_index == 0) error stop 'VTK actual physical time missing'
  read(title(time_index+2:),*) stored_time
  if (abs(stored_time-integration_time) > 1.0d-14) error stop 'VTK physical time differs from integration_time'

  call execute_command_line('rm -f '//expected_file)
  call execute_command_line('rmdir '//test_directory//'/'//trim(scheme_name))
  call execute_command_line('rmdir '//test_directory)
  deallocate(x, y, z, x_face, y_face, z_face, p, alpha_c, beta_c, var_names)
  write(*,'(A)') 'OUTPUT PATH TEST PASSED'
end program test_output_paths
