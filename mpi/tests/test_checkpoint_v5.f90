program test_checkpoint_v5
  use variables
  use output
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  character(len=1024) :: argument, checkpoint_path
  integer :: io_status
  integer :: failures, unit_file, step_read, i, j, k, component, version
  real*8 :: time_read
  real*8, allocatable :: expected_state(:,:,:,:)
  character(len=*), parameter :: test_directory = '/tmp/grial_mpi_checkpoint_v8_test'
  character(len=*), parameter :: checkpoint_file = &
    '/tmp/grial_mpi_checkpoint_v8_test/checkpoint_checkpoint_test_step_000000073.rst'

  ! The Python control supplies RLIMIT_FSIZE. GNU Fortran installs a startup
  ! handler for Linux SIGXFSZ, so ignore it here to exercise the returned error.
  if(command_argument_count()==2) then
    call get_command_argument(1,argument)
    if(trim(argument)/='--rewrite-failure') error stop 'Unknown checkpoint probe mode'
    call get_command_argument(2,checkpoint_path)
    open(newunit=unit_file,file=trim(checkpoint_path),status='old',form='unformatted',action='read')
    call read_checkpoint_metadata(unit_file)
    allocate(up(nx,ny,nz,neq))
    call read_checkpoint_state(unit_file,step_read,time_read)
    close(unit_file)
    call signal(25,1)
    call save_checkpoint(step_read,time_read,io_status)
    if(io_status==0) error stop 'Expected limited-file write to fail'
    print *, 'CHECKPOINT WRITE FAILURE HANDLED'
    stop
  end if
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
  fd_polytropic_constant = 0.00234d0
  fm_inner_radius=8.5d0;fm_pressure_max_radius=13.25d0;fm_polytropic_constant=0.00321d0
  shock_sensor_thresholds = [0.80d0,0.40d0,1.50d0,0.90d0,10.0d0,10.0d0,100.0d0]

  rho_floor=2.0d-11; p_floor=7.0d-15; D_floor=3.0d-11; tau_floor=9.0d-15
  checkpoint_interval=7.5d0; extraction_stride=17; mass_monitor_stride=0
  do_numerical_diagnostics=.true.; vtk_mapping_id=VTK_MAP_UNTWISTED
  walltime_limit_seconds=86400.0d0; walltime_margin_seconds=900.0d0
  walltime_checkpoint_seconds=3600.0d0
  numerical_events=[(int(i,int64),i=1,NUM_EVENT_COUNT)]
  numerical_delta_u=[1.0d-9,2.0d-9,3.0d-9,4.0d-9,5.0d-9]
  do_radial_mass_budget=.true.;radial_mass_initial=987.125d0;radial_mass_flux_integral=[-0.125d0,-0.25d0]
  call save_checkpoint(73, 12.5d0)

  rho_floor=-1.0d0; p_floor=-1.0d0; D_floor=-1.0d0; tau_floor=-1.0d0
  checkpoint_interval=-1.0d0; extraction_stride=-1; mass_monitor_stride=-1
  do_numerical_diagnostics=.false.; vtk_mapping_id=-1
  walltime_limit_seconds=-1.0d0; walltime_margin_seconds=-1.0d0; walltime_checkpoint_seconds=-1.0d0
  numerical_events=0_int64;numerical_delta_u=0.0d0
  radial_mass_initial=0.0d0;radial_mass_flux_integral=0.0d0;do_radial_mass_budget=.false.
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
  fd_polytropic_constant = -1.0d0
  fm_inner_radius=-1.0d0;fm_pressure_max_radius=-1.0d0;fm_polytropic_constant=-1.0d0
  shock_sensor_thresholds = SENSOR_DEFAULTS

  open(newunit=unit_file,file=checkpoint_file,status='old',form='unformatted',action='read')
  call read_checkpoint_metadata(unit_file)
  call read_checkpoint_state(unit_file,step_read,time_read)
  close(unit_file)

  if(any(numerical_events/=[(int(i,int64),i=1,NUM_EVENT_COUNT)])) error stop 'V9 counters changed'
  if(any(numerical_delta_u/=[1.0d-9,2.0d-9,3.0d-9,4.0d-9,5.0d-9])) error stop 'V9 repair budget changed'
  if(.not.do_radial_mass_budget.or.radial_mass_initial/=987.125d0.or. &
      any(radial_mass_flux_integral/=[-0.125d0,-0.25d0])) error stop 'V9 flux budget changed'
  if (rho_floor/=2.0d-11.or.p_floor/=7.0d-15.or.D_floor/=3.0d-11.or.tau_floor/=9.0d-15) &
    error stop 'V8 effective floors changed'
  if (checkpoint_interval/=7.5d0.or.extraction_stride/=17.or.mass_monitor_stride/=0.or. &
      .not.do_numerical_diagnostics.or.vtk_mapping_id/=VTK_MAP_UNTWISTED) error stop 'V8 controls changed'
  if (walltime_limit_seconds/=86400.0d0.or.walltime_margin_seconds/=900.0d0.or. &
      walltime_checkpoint_seconds/=3600.0d0) error stop 'V8 walltime controls changed'

  ! Force a failed temporary open; the published checkpoint must remain readable.
  call execute_command_line('mkdir -p '//checkpoint_file//'.tmp')
  call save_checkpoint(73,13.0d0,io_status)
  if(io_status==0) error stop 'Checkpoint write failure was hidden'
  call execute_command_line('rmdir '//checkpoint_file//'.tmp')
  ! A failed publication must also retain the previous generation.
  call execute_command_line('mkdir -p '//test_directory//'/checkpoint_checkpoint_test_step_000000074.rst')
  call save_checkpoint(74,13.0d0,io_status)
  if(io_status==0) error stop 'Checkpoint rename failure was hidden'
  call execute_command_line('rmdir '//test_directory//'/checkpoint_checkpoint_test_step_000000074.rst')
  open(newunit=unit_file,file=test_directory//'/checkpoint_checkpoint_test_step_000000074.rst.tmp',status='old')
  close(unit_file,status='delete')
  open(newunit=unit_file,file=checkpoint_file,status='old',form='unformatted',action='read')
  call read_checkpoint_metadata(unit_file)
  call read_checkpoint_state(unit_file,step_read,time_read)
  close(unit_file)
  if(time_read/=12.5d0.or.any(up/=expected_state)) error stop 'Failed write damaged valid generation'

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
  call assert_close('FD parameter', fd_polytropic_constant, 0.00234d0, failures)
  call assert_close('FM inner radius',fm_inner_radius,8.5d0,failures)
  call assert_close('FM pressure maximum',fm_pressure_max_radius,13.25d0,failures)
  call assert_close('FM polytropic constant',fm_polytropic_constant,0.00321d0,failures)
  call assert_close('sensor B WENO3 pressure', shock_sensor_thresholds(1), 0.80d0, failures)
  call assert_close('sensor B TVD density', shock_sensor_thresholds(4), 0.90d0, failures)
  if (any(up /= expected_state)) then
    failures = failures + 1
    write(*,'(A)') 'FAIL conserved state round trip'
  end if

  case_id=13;case_name='FishMoncEqu'
  call save_checkpoint(73,12.5d0)
  case_id=-1;case_name='invalid';fm_polytropic_constant=-1.0d0
  open(newunit=unit_file,file=checkpoint_file,status='old',form='unformatted')
  call read_checkpoint_metadata(unit_file)
  call read_checkpoint_state(unit_file,step_read,time_read)
  close(unit_file)
  call assert_integer('V7 FM problem id',case_id,13,failures)
  call assert_close('V7 FM parameter',fm_polytropic_constant,0.00321d0,failures)
  if(trim(case_name)/='FishMoncEqu'.or.any(up/=expected_state)) error stop 'V7 FM restart changed identity/state'

  ! Sintetizar V5/V6 con su secuencia original de registros: compatible con
  ! estados SoA antiguos y umbrales A implícitos, sin leer un registro de más.
  do version=5,6
  case_id=11;case_name='FishMoncSag'
  shock_sensor_thresholds = [0.80d0,0.40d0,1.50d0,0.90d0,10.0d0,10.0d0,100.0d0]
  open(newunit=unit_file,file=checkpoint_file,status='replace',form='unformatted')
  if(version==5) then
    write(unit_file) 'GRHDCP_V5   '
  else
    write(unit_file) 'GRHDCP_V6   '
  end if
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
  if(version==6) write(unit_file) shock_sensor_thresholds
  write(unit_file) 12.5d0
  write(unit_file) 73
  write(unit_file) expected_state
  close(unit_file)
  up=-1.0d0
  open(newunit=unit_file,file=checkpoint_file,status='old',form='unformatted')
  call read_checkpoint_metadata(unit_file)
  call read_checkpoint_state(unit_file,step_read,time_read)
  close(unit_file)
  if (any(up /= expected_state)) error stop 'legacy FD state mismatch'
  if(version==5.and.any(shock_sensor_thresholds/=SENSOR_DEFAULTS)) error stop 'V5 sensor defaults mismatch'
  if(version==6.and.shock_sensor_thresholds(1)/=0.80d0) error stop 'V6 sensor metadata mismatch'
  if(case_id/=11.or.trim(case_name)/='FontDaigneSag') error stop 'Legacy FM name was not identified as FD'
  call assert_close('legacy FM defaults independent from FD',fm_inner_radius,8.0d0,failures)
  call assert_close('V5 time',time_read,12.5d0,failures)
  call assert_integer('V5 steps',step_read,73,failures)
  end do
  open(newunit=unit_file,file=checkpoint_file,status='old',form='unformatted')
  close(unit_file,status='delete')
  call execute_command_line('rmdir '//test_directory)
  deallocate(up,expected_state)

  if (failures /= 0) then
    write(*,'(A,I0)') 'CHECKPOINT V5 TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,'(A)') 'CHECKPOINT V5/V6/V9 FD/FM, ACCUMULATORS AND ATOMIC WRITE TESTS PASSED'

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
