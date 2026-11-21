! =======================================================================================
! PROGRAMA PRINCIPAL: GRIAL (General Relativistic Instability Accretion Laboratory)
! ---------------------------------------------------------------------------------------
! Orquestador central de la simulación. Controla el bucle temporal principal, la
! integración espacial mediante el Método de Líneas (Method of Lines) delegada al
! módulo 'evolution', la actualización temporal usando Runge-Kutta SSP de 3er orden,
! y el paralelismo de memoria compartida masivo (OpenMP).
! =======================================================================================
program grial
  use parallel_runtime
  use parallel_application, only: agree_command_line,run_distributed_application
  use parallel_timing, only: time_decisions,coordinate_time_actions,clip_step_to_target,any_time_request,next_periodic_index
  use variables
  use parameters
  use metrics
  use initialization
  use conditions
  use equations
  use fluxes
  use reconstruction
  use evolution
  use time_integrator, only: advance_rk_step
  use accelerator
  use output
  implicit none

  real*8 :: start_time, end_time
  integer :: i, j, k, argument_count, io_status
  character(len=1024) :: first_argument, parameter_file
  type(time_decisions) :: decisions
  logical :: wall_stop_request,wall_checkpoint_request
  logical :: check_only, stopped_for_walltime
  integer :: last_checkpoint_step
  integer(int64) :: save_slot
  real*8 :: last_wall_checkpoint, step_wall_start, max_step_wall, now_wall
  character(len=512) :: last_checkpoint_file

  ! ---  VARIABLES PARA EL CRONÓMETRO Y ETA ---
  real*8 :: session_start_time, session_progress
  real*8 :: current_wtime, elapsed_time, eta
  integer :: current_percent, last_percent
  integer :: eta_h, eta_m, eta_s
  integer :: elap_h, elap_m, elap_s

  real*8 :: total_mass
  integer(int64) :: restart_events(NUM_EVENT_COUNT)
  real*8 :: restart_delta_u(neq)

  start_time = wall_time_seconds()
  call parallel_initialize()
  stopped_for_walltime = .false.
  last_checkpoint_step = -1
  last_checkpoint_file = 'none'
  max_step_wall = 0.0d0
  last_wall_checkpoint = start_time
  argument_count = command_argument_count()
  check_only = .false.
  first_argument = ''
  parameter_file = ''

  call agree_command_line()
  if (argument_count>0) call get_command_argument(1,first_argument)
  if (trim(first_argument)=='--distributed') then
    call run_distributed_application(start_time)
    call parallel_finalize()
    stop
  end if
  if (argument_count==1) then
    call get_command_argument(1,first_argument)
    if (trim(first_argument)=='--parallel-check') then
      call parallel_require_collective(parallel_size>=1.and.parallel_rank>=0.and. &
        parallel_rank<parallel_size,'invalid MPI topology')
      if (parallel_rank==0) write(*,'(A,L1,A,I0,A,I0)') 'PARALLEL_CHECK mpi=',parallel_enabled, &
        ' size=',parallel_size,' thread_level=',parallel_thread_level
      call parallel_finalize()
      stop
    end if
  end if
  call parallel_require_collective(parallel_size==1, &
    'Distributed evolution is not enabled yet; use --parallel-check or one process')
  if (argument_count == 1) then
    call get_command_argument(1, parameter_file)
  else if (argument_count == 2) then
    call get_command_argument(1, first_argument)
    call get_command_argument(2, parameter_file)
    if (trim(first_argument) /= '--check') then
      write(*,*) 'Uso: ./grial [--check] archivo.par'
      call parallel_abort('expected --check before the parameter file')
    end if
    check_only = .true.
  else
    write(*,*) 'Uso: ./grial [--check] archivo.par'
    call parallel_abort('invalid command line')
  end if

  write(*,'(A)') 'GRIAL: General Relativistic Instability Accretion Laboratory'
  call read_parameter_header(trim(parameter_file))

  ! =========================================================================
  ! 1. ARRANQUE DEL SISTEMA (Cold Start vs Checkpoint Restart)
  ! =========================================================================
  if (do_restart) then
    print *, '=========================================='
    print *, '>>> REANUDANDO SIMULACION DESDE BINARIO'
    print *, '>>> Archivo: ', trim(restart_file)
    print *, '=========================================='

    open(20, file=trim(restart_file), status='old', form='unformatted', iostat=io_status)
    if (io_status /= 0) then
      write(*,*) 'ERROR: no se pudo abrir el checkpoint: ', trim(restart_file)
      call parallel_abort('checkpoint open failed: '//trim(restart_file))
    end if
    call read_checkpoint_metadata(20)
    call read_parameter_file(trim(parameter_file), .true.)
    ! Checkpoints store only physical conserved cells. Re-derive ghost depth
    ! so a historical TVD/WENO3 checkpoint cannot restore an insufficient
    ! filled stencil after the five-cell sensor boundary correction.
    call choose_numerical_architecture()
    call validate_parameters(trim(parameter_file))
    call print_parameter_summary(trim(parameter_file))
    if (check_only) then
      allocate(up(nx,ny,nz,neq))
      call read_checkpoint_state(20,n_steps,integration_time)
      deallocate(up)
      close(20)
      write(*,'(A,I0,A,ES24.16E3)') 'CHECKPOINT_STATE step=',n_steps,' time=',integration_time
      print *, 'Archivo de parametros y checkpoint validos; no se inicio la simulacion.'
      call parallel_finalize()
      stop
    end if

    ! Recalculamos parámetros termodinámicos auxiliares
    g1 = adb_idx / (adb_idx - 1.0d0)

    call set_metric_type()
    call allocate_and_grid()
    call precalculate_metric_cache()
    call init_wavespeed_solver()
    call read_checkpoint_state(20, n_steps, integration_time)
    close(20)

    ! Algunos bordes necesitan datos analíticos auxiliares que no están en U.
    ! Se reconstruyen desde los parámetros guardados antes de recuperar P.
    call prepare_case_boundary_data()

    ! Sincronizamos y recuperamos primitivas para el integrador RK3
    u = up
    restart_events=numerical_events
    restart_delta_u=numerical_delta_u

    ! Propagación de y(j) para la inversión termodinámica 3D
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call cons_to_prim(u(i,j,k,:), p(i,j,k,:), i, j, k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call set_boundary_conditions(p)
    numerical_events=restart_events
    numerical_delta_u=restart_delta_u
    next_checkpoint = dble(next_periodic_index(integration_time,checkpoint_interval))*checkpoint_interval
  else
    ! --- ARRANQUE EN FRÍO NORMAL (T=0) ---
    integration_time = 0.0d0
    call initialize_problem(trim(parameter_file), check_only)
    if (check_only) then
      call parallel_finalize()
      stop
    end if
    next_checkpoint = checkpoint_interval

    ! Llenado inicial del vector de conservativas
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call prim_to_cons(p(i,j,k,:), up(i,j,k,:), i, j, k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call save_vtk_at_time(integration_time)
  end if

  call initialize_accelerator()
  call initialize_rhs_accelerator()

  call initialize_auxiliary_output(.not. do_restart)
  call initialize_radial_mass_output()
  call initialize_numerical_output()
  if (do_gw_extraction) call calc_gw_strain(integration_time)
  if (do_mdot_extraction) call calc_m_dot(integration_time)
  if (do_ppi_diagnostics) then
    call calc_ppi_modes(integration_time)
    call calc_global_diagnostics(integration_time)
  end if

  ! -------------------------------------------
  ! Mantiene el calendario absoluto de salidas al reiniciar; desplazarlo desde
  ! integration_time omitiría los frames nominales posteriores al checkpoint.
  save_slot = next_periodic_index(integration_time,save_interval)
  next_save_time = dble(save_slot)*save_interval
  call update_adaptive_dt() ! Inicialización del primer dt

  session_start_time = integration_time
  if (final_time > 0.0d0) then
    last_percent = int((integration_time / final_time) * 100.0d0)
  else
    last_percent = 100
  end if

    total_mass = 0.0d0

    ! Integral de volumen lógico 3D Coordenado
    ! (El Jacobiano dr/dx ya está matemáticamente embebido en el sqrt(gamma) dentro de u(eq_de))
    !$OMP PARALLEL DO REDUCTION(+:total_mass) PRIVATE(i, j, k)
    !$ACC PARALLEL LOOP COLLAPSE(3) REDUCTION(+:total_mass) PRESENT(up)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          total_mass = total_mass + up(i,j,k,eq_de) * dx * dy * dz
        end do
      end do
    end do
    !$ACC END PARALLEL LOOP
    !$OMP END PARALLEL DO

    write(*,*) 'Initial Mass in domain: ', total_mass

  ! =========================================================================
  ! 2. BUCLE PRINCIPAL DE INTEGRACIÓN TEMPORAL (SSP-RK3)
  ! =========================================================================
  call publish_run_status('running')
  do
    now_wall=wall_time_seconds()
    wall_stop_request=.false.
    if (walltime_limit_seconds>0.0d0) &
      wall_stop_request=now_wall-start_time+max_step_wall>=walltime_limit_seconds-walltime_margin_seconds
    call coordinate_time_actions(integration_time,n_steps,final_time, &
      [perturbation_time,next_save_time,next_checkpoint], &
      [apply_perturbation.and..not.perturbation_applied,.true.,.true.], &
      [wall_stop_request,.false.],decisions)
    if (decisions%complete) exit
    if (decisions%stop) then
      stopped_for_walltime=.true.
      exit
    end if
    step_wall_start=now_wall
    call update_adaptive_dt()
    dt=clip_step_to_target(integration_time,dt,final_time)

    ! The delayed perturbation keeps the existing start-of-step threshold rule.
    if (decisions%perturb) then
      call accelerator_update_host()
      print *, ""
      print *, "=========================================================="
      print *, ">>> TIEMPO DE PERTURBACION ALCANZADO: ", integration_time
      print *, "=========================================================="

      select case(perturbation_type)
      case(PERT_PRESSURE_NOISE)
        call inject_pressure_noise(p, perturbation_amplitude, perturbation_seed)
      case(PERT_DENSITY_NOISE)
        call inject_density_noise(p, perturbation_amplitude, perturbation_seed)
      case(PERT_DENSITY_MODE)
        call inject_density_mode(p, perturbation_mode, perturbation_amplitude)
      case default
        write(*,*) 'Perturbation disabled: type = PERT_NONE.'
        apply_perturbation = .false.
      end select

      if (perturbation_type /= PERT_NONE) then
        ! Re-sincronizamos U y UP tras modificar P.
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              call prim_to_cons(p(i,j,k,:), u(i,j,k,:), i, j, k)
              up(i,j,k,:) = u(i,j,k,:)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        call set_boundary_conditions(p)
        call accelerator_update_device()

        perturbation_applied = .true.
        call write_perturbation_event(integration_time)
        if (do_ppi_diagnostics) then
          call calc_ppi_modes(integration_time)
          call calc_global_diagnostics(integration_time)
        end if
      end if
    end if

    n_steps = n_steps + 1
    integration_time = integration_time + dt

    call advance_rk_step()

    call coordinate_time_actions(integration_time,n_steps,final_time, &
      [perturbation_time,next_save_time,next_checkpoint], &
      [apply_perturbation.and..not.perturbation_applied,.true.,.true.],[.false.,.false.],decisions)

    ! =====================================================================
    ! 3. DIAGNÓSTICOS Y SALIDA DE DATOS
    ! =====================================================================
    ! Corrección: Sustitución de 'is_flat' por la validación arquitectónica global
    if (trim(metric_type) /= 'Minkowski' .and. bh_mass > 0.0d0) then
      if (mod(n_steps, extraction_stride) == 0) then
        call accelerator_update_host()
        if (do_gw_extraction) call calc_gw_strain(integration_time)
        if (do_mdot_extraction) call calc_m_dot(integration_time)
      end if
    end if
    if (do_ppi_diagnostics .and. mod(n_steps, diagnostic_stride) == 0) then
      call accelerator_update_host()
      call calc_ppi_modes(integration_time)
      call calc_global_diagnostics(integration_time)
    end if
    if (do_numerical_diagnostics .and. mod(n_steps,diagnostic_stride) == 0) &
      call write_numerical_diagnostics()

    ! --- CONSERVACIÓN DE MASA GLOBAL ---
    if (mass_monitor_stride > 0) then
      if (mod(n_steps, mass_monitor_stride) == 0) then
        total_mass = 0.0d0

      ! Integral de volumen lógico 3D Coordenado
      ! (El Jacobiano dr/dx ya está matemáticamente embebido en el sqrt(gamma) dentro de u(eq_de))
        !$OMP PARALLEL DO REDUCTION(+:total_mass) PRIVATE(i, j, k)
        !$ACC PARALLEL LOOP COLLAPSE(3) REDUCTION(+:total_mass) PRESENT(up)
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              total_mass = total_mass + up(i,j,k,eq_de) * dx * dy * dz
            end do
          end do
        end do
        !$ACC END PARALLEL LOOP
        !$OMP END PARALLEL DO

        write(*,*) 'Total Mass in domain: ', total_mass
      end if
    end if

    ! --- CRONÓMETRO INTELIGENTE (ETA) ---
    current_percent = int((integration_time / final_time) * 100.0d0)
    if (current_percent > last_percent) then
      current_wtime = wall_time_seconds()
      elapsed_time = current_wtime - start_time

      session_progress = (integration_time - session_start_time) / (final_time - session_start_time)

      if (session_progress > 0.0d0) then
        eta = (elapsed_time / session_progress) - elapsed_time

        eta_h = int(eta) / 3600
        eta_m = mod(int(eta), 3600) / 60
        eta_s = mod(int(eta), 60)

        elap_h = int(elapsed_time) / 3600
        elap_m = mod(int(elapsed_time), 3600) / 60
        elap_s = mod(int(elapsed_time), 60)

        write(*, '(A, I3, A, I3, A, I2.2, A, I2.2, A, I3, A, I2.2, A, I2.2, A, A, I8, A, F10.3)') &
          '>>> PROGRESO: ', current_percent, '% | Transcurrido: ', &
          elap_h, 'h ', elap_m, 'm ', elap_s, 's | Faltan aprox: ', &
          eta_h, 'h ', eta_m, 'm ', eta_s, 's  <<< | ', &
          'Paso ', n_steps, ' | Tiempo: ', integration_time
      end if

      last_percent = current_percent
    end if

    ! --- EJECUCIÓN DEL GUARDADO A DISCO ---
    if (decisions%save) then
      call accelerator_update_host()
      call save_vtk_at_time(next_save_time)
      save_slot = save_slot+1_int64
      next_save_time = dble(save_slot)*save_interval
    end if

    now_wall=wall_time_seconds()
    wall_checkpoint_request=walltime_checkpoint_seconds>0.0d0.and. &
      now_wall-last_wall_checkpoint>=walltime_checkpoint_seconds
    ! Evaluate the collective before the conditional; logical evaluation may short-circuit.
    wall_checkpoint_request=any_time_request(wall_checkpoint_request)
    if (decisions%checkpoint.or.wall_checkpoint_request) then
      call checkpoint_now()
      if (decisions%checkpoint) &
        next_checkpoint = dble(next_periodic_index(integration_time,checkpoint_interval))*checkpoint_interval
    end if
    max_step_wall=max(max_step_wall,wall_time_seconds()-step_wall_start)
  end do

  ! Garantiza una muestra diagnóstica exactamente en el tiempo final aunque
  ! el último paso no coincida con el stride periódico.
  if (do_ppi_diagnostics .and. mod(n_steps, diagnostic_stride) /= 0) then
    call accelerator_update_host()
    call calc_ppi_modes(integration_time)
    call calc_global_diagnostics(integration_time)
  end if
  if (trim(metric_type) /= 'Minkowski' .and. bh_mass > 0.0d0 .and. &
      mod(n_steps, extraction_stride) /= 0) then
    call accelerator_update_host()
    if (do_gw_extraction) call calc_gw_strain(integration_time)
    if (do_mdot_extraction) call calc_m_dot(integration_time)
  end if

  ! Evaluación final de errores L1/L2 si corresponde al test
  if (do_numerical_diagnostics .and. mod(n_steps,diagnostic_stride) /= 0) &
    call write_numerical_diagnostics()
  call accelerator_update_host()
  call calc_convergence_norms()

  ! A final complete generation also covers a stop during initialization (step zero).
  if (last_checkpoint_step /= n_steps) call checkpoint_now()
  if (stopped_for_walltime) then
    call publish_run_status('stopped_for_walltime')
  else
    call publish_run_status('completed')
  end if

  ! =========================================================================
  ! 4. FINALIZACIÓN Y LIMPIEZA
  ! =========================================================================
  end_time = wall_time_seconds()
  elapsed_time = end_time - start_time

  elap_h = int(elapsed_time) / 3600
  elap_m = mod(int(elapsed_time), 3600) / 60
  elap_s = mod(int(elapsed_time), 60)

  write(*, '(/, A)') '========================================================='
  if (stopped_for_walltime) then
    write(*,'(A,F14.6)') ' PARADA PLANIFICADA POR WALLTIME; REINICIO NECESARIO. T = ', integration_time
  else
    write(*, '(A, F10.3)') ' SIMULACIÓN FINALIZADA EXITOSAMENTE EN T = ', integration_time
  end if
  write(*, '(A, I8)')    ' TOTAL DE PASOS TEMPORALES: ', n_steps
  write(*, '(A, I3, A, I2.2, A, I2.2, A)') &
    ' TIEMPO TOTAL DE EJECUCIÓN: ', elap_h, 'h ', elap_m, 'm ', elap_s, 's'
  write(*, '(A, /)') '========================================================='

  call finalize_rhs_accelerator()
  call finalize_accelerator()

  deallocate(x,y,z,var_names)
  deallocate(u,up,s,rhs,p)
  deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)

  deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)

  if (ny > 1) then
    deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
  end if

  if (nz > 1) then
    deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
  end if

  call parallel_finalize()

contains

  subroutine checkpoint_now()
    call accelerator_update_host()
    call save_checkpoint(n_steps,integration_time)
    last_checkpoint_step=n_steps
    last_wall_checkpoint=wall_time_seconds()
    write(last_checkpoint_file,'(A,"/checkpoint_",A,"_step_",I9.9,".rst")') &
      trim(output_folder),trim(scheme_name),n_steps
    call publish_run_status('running')
  end subroutine checkpoint_now

  subroutine publish_run_status(reason)
    character(len=*), intent(in) :: reason
    character(len=512) :: path, temporary_path
    integer :: unit_file, status, close_status
    path=trim(output_folder)//'/run_status.dat'
    temporary_path=trim(path)//'.tmp'
    open(newunit=unit_file,file=trim(temporary_path),status='replace',action='write',iostat=status)
    if (status /= 0) call parallel_abort('Cannot create run_status.dat')
    write(unit_file,'(A)',iostat=status) 'status='//trim(reason)
    if (status == 0) write(unit_file,'(A,ES24.16E3)',iostat=status) 'physical_time=',integration_time
    if (status == 0) write(unit_file,'(A,I0)',iostat=status) 'step=',n_steps
    if (status == 0) write(unit_file,'(A,ES24.16E3)',iostat=status) 'target_time=',final_time
    if (status == 0) write(unit_file,'(A)',iostat=status) 'checkpoint='//trim(last_checkpoint_file)
    if (status == 0) write(unit_file,'(A,ES24.16E3)',iostat=status) &
      'elapsed_wall_seconds=',wall_time_seconds()-start_time
    close(unit_file,iostat=close_status)
    if (status /= 0 .or. close_status /= 0) call parallel_abort('Cannot write run_status.dat')
    call rename(trim(temporary_path),trim(path),status)
    if (status /= 0) call parallel_abort('Cannot publish run_status.dat')
  end subroutine publish_run_status

end program grial
