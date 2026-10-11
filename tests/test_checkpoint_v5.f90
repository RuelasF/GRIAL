program test_checkpoint_v5
  use variables
  use output
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer :: failures, unit_file, step_read, i, j, k, component
  real*8 :: time_read
  real*8, allocatable :: expected_state(:,:,:,:)
  character(len=*), parameter :: test_directory = '/tmp/grhd_checkpoint_v5_test'
  character(len=*), parameter :: checkpoint_file = &
    '/tmp/grhd_checkpoint_v5_test/checkpoint_checkpoint_test_step_000000073.rst'

  failures = 0
  call execute_command_line('mkdir -p '//test_directory)

  nghost = 1
  nx = 2
  ny = 1
  nz = 2
  allocate(up(nx,ny,nz,neq), expected_state(nx,ny,nz,neq))
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        do component = 1, neq
          up(i,j,k,component) = dble(1000*component + 100*i + 10*j + k)
        end do
      end do
    end do
  end do
  expected_state = up

  r_min = 1.125d0
  r_max = 32.5d0
  y_min = 0.0d0
  y_max = pi
  z_min = 0.0d0
  z_max = 2.0d0*pi
  bh_mass = 1.25d0
  adb_idx = 1.4d0
  a_spin = 0.35d0
  use_log_r = .true.
  use_shock_sensor = .false.
  do_gw_extraction = .true.
  do_mdot_extraction = .true.
  do_ppi_diagnostics = .true.
  metric_type = 'Kerr-Schild'
  geom_type = 'Spheroidal'
  final_time = 25.0d0
  CFL = 0.31d0
  save_interval = 2.0d0
  output_prefix = 'checkpoint_v4'
  output_folder = test_directory
  case_id = 7
  case_name = 'MichelA'
  scheme_name = 'checkpoint_test'
  rec_method_id = REC_WENO3
  tvd_limiter_id = LIM_MC
  riemann_solver_id = RS_HLLE
  apply_perturbation = .true.
  perturbation_applied = .true.
  perturbation_type = PERT_DENSITY_MODE
  perturbation_seed = 9876
  diagnostic_stride = 37
  perturbation_time = 4.5d0
  perturbation_amplitude = 0.023d0
  perturbation_mode = 3.0d0

  sod_interface = 0.41d0
  strong_pressure_left = 777.0d0
  shu_wave_number = 23.0d0
  khi_perturbation_width = 0.037d0
  jet_density = 0.234d0
  jet_velocity = 0.876d0
  advected_wave_amplitude = 0.043d0
  michel_critical_radius = 17.25d0
  michel_critical_density = 0.2345d0
  dust_accretion_constant = -0.765d0
  offaxis_radial_center = 12.75d0
  fm_polytropic_constant = 0.00234d0
  shock_sensor_thresholds = [0.80d0,0.40d0,1.50d0,0.90d0,10.0d0,10.0d0,100.0d0]

  call save_checkpoint(73, 12.5d0)

  up = -1.0d0
  r_min = -1.0d0
  bh_mass = -1.0d0
  a_spin = -1.0d0
  case_id = -1
  rec_method_id = -1
  perturbation_seed = -1
  sod_interface = -1.0d0
  strong_pressure_left = -1.0d0
  shu_wave_number = -1.0d0
  khi_perturbation_width = -1.0d0
  jet_density = -1.0d0
  jet_velocity = -1.0d0
  advected_wave_amplitude = -1.0d0
  michel_critical_radius = -1.0d0
  michel_critical_density = -1.0d0
  dust_accretion_constant = 1.0d0
  offaxis_radial_center = -1.0d0
  fm_polytropic_constant = -1.0d0
  shock_sensor_thresholds = SENSOR_DEFAULTS

  open(newunit=unit_file,file=checkpoint_file,status='old',form='unformatted',action='read')
  call read_checkpoint_metadata(unit_file)
  call read_checkpoint_state(unit_file,step_read,time_read)
  close(unit_file)

  call assert_close('r_min', r_min, 1.125d0, failures)
  call assert_close('bh_mass', bh_mass, 1.25d0, failures)
  call assert_close('spin', a_spin, 0.35d0, failures)
  call assert_integer('case id', case_id, 7, failures)
  call assert_integer('step', step_read, 73, failures)
  call assert_integer('reconstruction', rec_method_id, REC_WENO3, failures)
  call assert_integer('perturbation seed', perturbation_seed, 9876, failures)
  call assert_close('time', time_read, 12.5d0, failures)
  call assert_close('Sod parameter', sod_interface, 0.41d0, failures)
  call assert_close('Strong parameter', strong_pressure_left, 777.0d0, failures)
  call assert_close('Shu parameter', shu_wave_number, 23.0d0, failures)
  call assert_close('KHI parameter', khi_perturbation_width, 0.037d0, failures)
  call assert_close('Jet density', jet_density, 0.234d0, failures)
  call assert_close('Jet velocity', jet_velocity, 0.876d0, failures)
  call assert_close('convergence parameter', advected_wave_amplitude, 0.043d0, failures)
  call assert_close('Michel radius', michel_critical_radius, 17.25d0, failures)
  call assert_close('Michel density', michel_critical_density, 0.2345d0, failures)
  call assert_close('Dust constant', dust_accretion_constant, -0.765d0, failures)
  call assert_close('OffAxis center', offaxis_radial_center, 12.75d0, failures)
  call assert_close('FM parameter', fm_polytropic_constant, 0.00234d0, failures)
  call assert_close('sensor B WENO3 pressure', shock_sensor_thresholds(1), 0.80d0, failures)
  call assert_close('sensor B TVD density', shock_sensor_thresholds(4), 0.90d0, failures)
  if (any(up /= expected_state)) then
    failures = failures + 1
    write(*,'(A)') 'FAIL conserved state round trip'
  end if

  ! Sintetizar V5 con su secuencia original de registros: compatible con
  ! estados SoA antiguos y umbrales A implícitos, sin leer un registro de más.
  open(newunit=unit_file,file=checkpoint_file,status='replace',form='unformatted')
  write(unit_file) 'GRHDCP_V5   '
  write(unit_file) nghost,nx,ny,nz
  write(unit_file) r_min,r_max,y_min,y_max,z_min,z_max
  write(unit_file) bh_mass,adb_idx,a_spin
  write(unit_file) use_log_r,use_shock_sensor,do_gw_extraction,do_mdot_extraction,do_ppi_diagnostics
  write(unit_file) metric_type,geom_type
  write(unit_file) final_time,CFL,save_interval
  write(unit_file) output_prefix,output_folder
  write(unit_file) case_id,case_name,scheme_name
  write(unit_file) rec_method_id,tvd_limiter_id,riemann_solver_id
  call write_checkpoint_initial_parameters(unit_file)
  write(unit_file) apply_perturbation,perturbation_applied
  write(unit_file) perturbation_type,perturbation_seed,diagnostic_stride
  write(unit_file) perturbation_time,perturbation_amplitude,perturbation_mode
  write(unit_file) 12.5d0
  write(unit_file) 73
  write(unit_file) expected_state
  close(unit_file)
  up=-1.0d0
  open(newunit=unit_file,file=checkpoint_file,status='old',form='unformatted')
  call read_checkpoint_metadata(unit_file)
  call read_checkpoint_state(unit_file,step_read,time_read)
  close(unit_file,status='delete')
  if (any(shock_sensor_thresholds /= SENSOR_DEFAULTS) .or. any(up /= expected_state)) &
    error stop 'legacy V5 sensor defaults/state mismatch'
  call assert_close('V5 time',time_read,12.5d0,failures)
  call assert_integer('V5 steps',step_read,73,failures)
  call execute_command_line('rmdir '//test_directory)
  deallocate(up,expected_state)

  if (failures /= 0) then
    write(*,'(A,I0)') 'CHECKPOINT V5 TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,'(A)') 'CHECKPOINT V5/V6 TESTS PASSED'

contains

  subroutine assert_close(label, actual, expected, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual, expected
    integer, intent(inout) :: failure_count

    if (.not. ieee_is_finite(actual) .or. &
        abs(actual-expected) > 64.0d0*epsilon(1.0d0)*max(1.0d0,abs(expected))) then
      failure_count = failure_count + 1
      write(*,'(A,1X,A,2(1X,ES16.8))') 'FAIL',trim(label),actual,expected
    end if
  end subroutine assert_close


  subroutine assert_integer(label, actual, expected, failure_count)
    character(len=*), intent(in) :: label
    integer, intent(in) :: actual, expected
    integer, intent(inout) :: failure_count

    if (actual /= expected) then
      failure_count = failure_count + 1
      write(*,'(A,1X,A,2(1X,I0))') 'FAIL',trim(label),actual,expected
    end if
  end subroutine assert_integer

end program test_checkpoint_v5
