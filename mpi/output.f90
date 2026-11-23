! =======================================================================================
! Módulo: output (Visualización, Checkpoints y Extracción de Observables)
! ---------------------------------------------------------------------------------------
! Propósito: Gestiona toda la salida de datos del simulador hacia el disco duro.
! Características Clave:
! 1. Exportación VTK Binario: Máxima velocidad de escritura (I/O) minimizando el
!    cuello de botella del disco duro. Genera archivos compatibles con VisIt y ParaView.
! 2. Tolerancia a Fallos: Sistema de Checkpoints que guarda un volcado binario exacto
!    de la memoria RAM para pausar y reanudar simulaciones de semanas de duración.
! 3. Física Multimensajero: Extracción al vuelo de la señal de Ondas Gravitacionales
!    (polarizaciones h+ y hx) y la tasa de acreción de masa al agujero negro.
! =======================================================================================
module output
  use variables
  use metrics
  use conditions
  use parallel_diagnostics, only: state_diagnostics,conservation_balance,collect_global_state, &
    collect_conservation_balance,reduce_observable_sum,validate_observable_sample
  use parallel_runtime, only: parallel_rank,parallel_require_collective
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  public :: save_vtk_at_time, save_checkpoint, read_checkpoint_metadata, &
            read_checkpoint_state, calc_gw_strain, calc_m_dot, calc_convergence_norms, &
            initialize_auxiliary_output, calc_ppi_modes, calc_global_diagnostics, &
            write_perturbation_event, finn_evans_cell_stress, gw_cartesian_kinematics, &
            kerr_schild_spheroidal_position
  public :: write_distributed_diagnostics,write_distributed_observables,write_parallel_perturbation_event
  public :: initialize_numerical_output, write_numerical_diagnostics
  public :: initialize_radial_mass_output, accumulate_radial_mass_budget

  ! La extracción Finn--Evans es un proxy de campo débil evaluado sobre el
  ! fluido GRHD. La distancia está expresada en las unidades geométricas del
  ! código y se registra en el manifiesto de cada corrida.
  real*8, parameter :: gw_observer_distance = 1000.0d0
  logical, parameter :: gw_exclude_inside_horizon = .true.
  logical, parameter :: gw_subtract_numerical_atmosphere = .true.
  logical :: checkpoint_restores_diagnostics = .false.
  character(len=12) :: checkpoint_footer_magic = 'GRIAL_END_V8'
  integer, parameter :: SERIES_GW=1, SERIES_MDOT=2, SERIES_PPI=3, SERIES_GLOBAL=4
  integer, parameter :: SERIES_EVENTS=5, SERIES_KHI=6, SERIES_NUMERICAL=7, SERIES_RADIAL=8
  character(len=512) :: series_paths(8) = ''
  real*8 :: series_last_time(8) = -huge(1.0d0)
  logical :: checkpoint_has_footer = .false.
  logical :: checkpoint_uses_soa_layout = .false.
  character(len=512) :: d01_series_root=''
  real*8 :: d01_last_time=-1.0d0
  integer :: d01_last_step=-1
  character(len=512) :: observable_paths(3)=''
  real*8 :: observable_last_time(3)=-1.0d0

  type :: distributed_output_cursor
    character(len=512) :: diagnostic_root='',observable_path(3)=''
    real*8 :: diagnostic_time=-1.0d0,observable_time(3)=-1.0d0
    integer :: diagnostic_step=-1
  end type distributed_output_cursor

contains

  subroutine export_distributed_output_cursor(cursor)
    type(distributed_output_cursor), intent(out) :: cursor
    cursor=distributed_output_cursor()
    if(parallel_rank/=0)return
    cursor%diagnostic_root=d01_series_root;cursor%diagnostic_time=d01_last_time
    cursor%diagnostic_step=d01_last_step
    cursor%observable_path=observable_paths;cursor%observable_time=observable_last_time
  end subroutine export_distributed_output_cursor

  logical function valid_distributed_output_cursor(cursor,time,step) result(valid)
    type(distributed_output_cursor), intent(in) :: cursor
    real*8, intent(in) :: time
    integer, intent(in) :: step
    character(len=512) :: root
    character(len=24), parameter :: names(3)=[character(len=24):: &
      'mpi_ppi_modes.dat','mpi_GW_signal.dat','mpi_m_dot.dat']
    integer :: i
    root=trim(output_folder)//'/'//trim(scheme_name)
    valid=all(ieee_is_finite([time,cursor%diagnostic_time,cursor%observable_time])).and.step>=0
    if(.not.valid)return
    if(len_trim(cursor%diagnostic_root)==0)then
      valid=cursor%diagnostic_time==-1.0d0.and.cursor%diagnostic_step==-1
    else
      valid=cursor%diagnostic_root==root.and.cursor%diagnostic_time>=0.0d0.and. &
        cursor%diagnostic_time<=time.and.cursor%diagnostic_step>=0.and.cursor%diagnostic_step<=step
    end if
    do i=1,3
      if(len_trim(cursor%observable_path(i))==0)then
        valid=valid.and.cursor%observable_time(i)==-1.0d0
      else
        valid=valid.and.cursor%observable_path(i)==trim(root)//'/'//trim(names(i)).and. &
          cursor%observable_time(i)>=0.0d0.and.cursor%observable_time(i)<=time
      end if
    end do
  end function valid_distributed_output_cursor

  subroutine restore_distributed_output_cursor(cursor)
    type(distributed_output_cursor), intent(in) :: cursor
    ! C03 validates histories and cursor before committing any fields or observer state.
    if(parallel_rank/=0)return
    d01_series_root=cursor%diagnostic_root;d01_last_time=cursor%diagnostic_time
    d01_last_step=cursor%diagnostic_step
    observable_paths=cursor%observable_path;observable_last_time=cursor%observable_time
  end subroutine restore_distributed_output_cursor

  subroutine write_distributed_observables(current_time)
    real*8, intent(in) :: current_time
    real*8 :: coefficients(1+2*PPI_MAX_MODE),values(2+4*PPI_MAX_MODE),stress(3),mdot(1),wave(3)
    values=0.0d0;wave=0.0d0
    call validate_observable_sample(current_time)
    if (do_ppi_diagnostics.and.nz>=2) then
      call calc_ppi_modes(current_time,coefficients)
      call reduce_observable_sum(coefficients)
      if (parallel_rank==0) call ppi_mode_values(current_time,coefficients,values)
      call write_observable_series(1,'mpi_ppi_modes.dat',current_time,values)
    end if
    if (trim(metric_type)=='Minkowski'.or.bh_mass<0.1d0) return
    if (do_gw_extraction) then
      call calc_gw_strain(current_time,stress)
      call reduce_observable_sum(stress)
      if (parallel_rank==0) wave=[current_time,(stress(1)-stress(2))/gw_observer_distance, &
        2.0d0*stress(3)/gw_observer_distance]
      call write_observable_series(2,'mpi_GW_signal.dat',current_time,wave)
    end if
    if (do_mdot_extraction) then
      call calc_m_dot(current_time,mdot(1))
      call reduce_observable_sum(mdot)
      call write_observable_series(3,'mpi_m_dot.dat',current_time,[current_time,mdot(1)])
    end if
  end subroutine write_distributed_observables

  subroutine write_observable_series(series_id,name,current_time,values)
    integer, intent(in) :: series_id
    character(len=*), intent(in) :: name
    real*8, intent(in) :: current_time,values(:)
    integer :: unit_file,status,close_status
    logical :: io_ok,new_series
    character(len=512) :: path
    path=trim(output_folder)//'/'//trim(scheme_name)//'/'//name
    io_ok=.true.
    if (parallel_rank==0) then
      new_series=path/=observable_paths(series_id)
      if (.not.new_series) io_ok=current_time>observable_last_time(series_id)
      io_ok=io_ok.and.all(ieee_is_finite(values))
      if (io_ok) then
        if (new_series) then
          open(newunit=unit_file,file=trim(path),status='new',action='write',iostat=status)
        else
          open(newunit=unit_file,file=trim(path),status='old',position='append',action='write',iostat=status)
        end if
        io_ok=status==0
        if (io_ok) then
          write(unit_file,'(*(ES24.16E3,1X))',iostat=status) values
          close(unit_file,iostat=close_status)
          io_ok=status==0.and.close_status==0
        end if
      end if
      if (io_ok) then
        observable_paths(series_id)=path;observable_last_time(series_id)=current_time
      end if
    end if
    call parallel_require_collective(io_ok,'cannot publish root observable: '//trim(name))
  end subroutine write_observable_series

  subroutine write_parallel_perturbation_event(current_time,rho_peak,selected,change)
    real*8, intent(in) :: current_time,rho_peak,change(neq)
    integer(int64), intent(in) :: selected
    integer :: unit_file,status,close_status
    logical :: io_ok
    character(len=512) :: path
    io_ok=.true.
    if (parallel_rank==0) then
      path=trim(output_folder)//'/'//trim(scheme_name)//'/mpi_perturbation_event.dat'
      open(newunit=unit_file,file=trim(path),status='new',action='write',iostat=status)
      if (status==0) then
        write(unit_file,'(A)',iostat=status) '# D03 legacy intrinsic RNG, global k/j/i order, global 5% body mask'
        if (status==0)write(unit_file,'(A)',iostat=status) &
          '# time step type seed amplitude mode rho_peak selected_cells delta_integral_U(5)'
        if (status==0)write(unit_file,'(ES24.16E3,3(1X,I0),3(1X,ES24.16E3),1X,I0,5(1X,ES24.16E3))',iostat=status) &
          current_time,n_steps,perturbation_type,perturbation_seed,perturbation_amplitude,perturbation_mode, &
          rho_peak,selected,change
        io_ok=status==0
        close(unit_file,iostat=close_status);io_ok=io_ok.and.close_status==0
      else
        io_ok=.false.
      end if
    end if
    call parallel_require_collective(io_ok,'cannot publish unique distributed perturbation event')
  end subroutine write_parallel_perturbation_event

  subroutine write_distributed_diagnostics(current_time)
    real*8, intent(in) :: current_time
    type(state_diagnostics) :: state
    type(conservation_balance) :: balance
    character(len=512) :: root,path
    character(len=40), parameter :: names(3)=[character(len=40):: &
      'mpi_global_diagnostics.dat','mpi_conservation_budget.dat','mpi_numerical_diagnostics.dat']
    logical :: valid,new_series,io_ok
    integer :: unit_file,status,close_status,file_id
    valid=ieee_is_finite(current_time)
    if (valid) valid=current_time==integration_time
    call parallel_require_collective(valid,'invalid distributed diagnostic sample time')
    call collect_global_state(state,root_only=.true.)
    call collect_conservation_balance(balance,root_only=.true.)
    root=trim(output_folder)//'/'//trim(scheme_name)
    io_ok=.true.
    if (parallel_rank==0) then
      new_series=root/=d01_series_root
      if (.not.new_series) io_ok=current_time>d01_last_time.and.n_steps>d01_last_step
      if (io_ok) then
        do file_id=1,3
          path=trim(root)//'/'//trim(names(file_id))
          if (new_series) then
            open(newunit=unit_file,file=trim(path),status='new',action='write',iostat=status)
          else
            open(newunit=unit_file,file=trim(path),status='old',position='append',action='write',iostat=status)
          end if
          if (status/=0) then
            io_ok=.false.;exit
          end if
          if (new_series) then
            write(unit_file,'(A)',iostat=status) '# D01 CPU global diagnostics; C03 restores validated append cursors'
            if (status==0) then
              select case(file_id)
              case(1)
                write(unit_file,'(A)',iostat=status) &
                  '# time step mass disk_mass S3 rho_min p_min rho_max p_max vmax2 '// &
                  'l1_D l1_tau l1_S1 l1_S2 l1_S3 l2_D l2_tau l2_S1 l2_S2 l2_S3 '// &
                  'linf_D linf_tau linf_S1 linf_S2 linf_S3 cells atmosphere_cells invalid_cells'
              case(2)
                write(unit_file,'(A)',iostat=status) '# time step Q_initial(5) Q_current(5) '// &
                  'integral_flux(5,2,3) integral_source(5) weighted_repair(5) perturbation(5) residual(5)'
                if (status==0) write(unit_file,'(A)',iostat=status) &
                  '# variable order D,tau,S1,S2,S3; flux lower/upper then x/y/z; outward net=upper-lower'
              case(3)
                write(unit_file,'(A)',iostat=status) '# time step events(16) raw_repair(5) '// &
                  'unique_face_positions_per_rhs work_faces_per_rhs extra_interface_evaluations'
                if (status==0) write(unit_file,'(A)',iostat=status) &
                  '# C2P events 8:13 count physical recovery calls; face/reconstruction/HLLC counts include MPI duplicates'
              end select
            end if
          end if
          if (status==0) then
            select case(file_id)
            case(1)
              write(unit_file,'(ES24.16E3,1X,I0,23(1X,ES24.16E3),3(1X,I0))',iostat=status) &
                current_time,n_steps,state%totals(eq_de),state%disk_mass,state%totals(eq_vz), &
                state%rho_min,state%p_min,state%rho_max,state%p_max,state%max_vsq, &
                state%l1,state%l2,state%linf,state%cells,state%atmosphere_cells,state%invalid_cells
            case(2)
              write(unit_file,'(ES24.16E3,1X,I0,*(1X,ES24.16E3))',iostat=status) &
                current_time,n_steps,balance%initial,balance%current,balance%flux, &
                balance%source,balance%repair,balance%perturbation,balance%residual
            case(3)
              write(unit_file,'(ES24.16E3,1X,I0,16(1X,I0),5(1X,ES24.16E3),3(1X,I0))',iostat=status) &
                current_time,n_steps,state%events,state%raw_repair,state%unique_faces_per_rhs, &
                state%work_faces_per_rhs,state%extra_interface_faces
            end select
          end if
          close(unit_file,iostat=close_status)
          if (status/=0.or.close_status/=0) then
            io_ok=.false.;exit
          end if
        end do
      end if
      if (io_ok) then
        d01_series_root=root;d01_last_time=current_time;d01_last_step=n_steps
      end if
    end if
    call parallel_require_collective(io_ok,'distributed diagnostic output failed, already exists, or sample is not increasing')
  end subroutine write_distributed_diagnostics


  subroutine initialize_radial_mass_output()
    integer :: unit_file
    logical :: resume_budget
    character(len=512) :: path
    resume_budget=do_restart.and.checkpoint_restores_diagnostics.and.do_radial_mass_budget
    do_radial_mass_budget = do_numerical_diagnostics .and. .not. accelerator_enabled .and. &
      trim(case_name) == 'MichelA' .and. ny == 1 .and. nz == 1
    if (.not. do_radial_mass_budget) return
    if (.not.resume_budget) then
      radial_mass_initial = sum(up(:,:,:,eq_de))*dx*dy*dz
      radial_mass_flux_integral = 0.0d0
    end if
    path=trim(output_folder)//'/'//trim(scheme_name)//'/radial_mass_budget.dat'
    call initialize_text_file(path, &
      '# time steps mass initial_mass integral_F_left integral_F_right residual',.not. do_restart,SERIES_RADIAL)
    open(newunit=unit_file,file=trim(path),position='append',action='write')
    write(unit_file,'(A,ES24.16E3)') '# radial flux budget initialized/restored at t=',integration_time
    write(unit_file,'(A)') '# F points towards increasing logical radius; residual=M-M0-(I_left-I_right)'
    close(unit_file)
  end subroutine initialize_radial_mass_output

  subroutine accumulate_radial_mass_budget(rk_stage)
    integer, intent(in) :: rk_stage
    real*8, parameter :: weights(3) = [1.0d0/6.0d0,1.0d0/6.0d0,2.0d0/3.0d0]
    radial_mass_flux_integral = radial_mass_flux_integral+dt*weights(rk_stage)*radial_mass_flux
  end subroutine accumulate_radial_mass_budget

  subroutine write_radial_mass_budget()
    integer :: unit_file
    real*8 :: mass, residual
    character(len=512) :: path
    mass = sum(up(:,:,:,eq_de))*dx*dy*dz
    residual = mass-radial_mass_initial-(radial_mass_flux_integral(1)-radial_mass_flux_integral(2))
    path=trim(output_folder)//'/'//trim(scheme_name)//'/radial_mass_budget.dat'
    if (.not.claim_series_sample(SERIES_RADIAL,path,integration_time)) return
    open(newunit=unit_file,file=trim(path),position='append',action='write')
    write(unit_file,'(ES24.16E3,1X,I0,5(1X,ES24.16E3))') integration_time,n_steps, &
      mass,radial_mass_initial,radial_mass_flux_integral,residual
    close(unit_file)
  end subroutine write_radial_mass_budget

  subroutine initialize_numerical_output()
    integer :: unit_file
    character(len=512) :: path
    if (.not. do_numerical_diagnostics) return
    path=trim(output_folder)//'/'//trim(scheme_name)//'/numerical_diagnostics.dat'
    call initialize_text_file(path, &
      '# time steps recon_godunov recon_tvd recon_weno3 recon_mp5 recon_weno5 face_rho_floor face_p_floor '// &
      'c2p_calls c2p_low_D c2p_momentum_clip c2p_root_fail c2p_bad_V c2p_sanitize p2c_velocity_cap '// &
      'hllc_calls hllc_fallback raw_delta_D raw_delta_tau raw_delta_Sx raw_delta_Sy raw_delta_Sz',.not. do_restart,SERIES_NUMERICAL)
    open(newunit=unit_file,file=trim(path),position='append',action='write')
    write(unit_file,'(A,ES24.16E3)') '# CPU counters cumulative with V9 checkpoints; startup recovery excluded; segment t=',integration_time
    write(unit_file,'(A)') '# recon counts include boundary faces and all RK stages; floors count states L/R'
    write(unit_file,'(A)') '# raw_delta sums C2P changes times logical cell volume; no RK weighting or boundary-flux budget'
    close(unit_file)
    call write_numerical_diagnostics()
  end subroutine initialize_numerical_output

  subroutine write_numerical_diagnostics()
    integer :: unit_file
    character(len=512) :: path
    if (.not. do_numerical_diagnostics) return
    path=trim(output_folder)//'/'//trim(scheme_name)//'/numerical_diagnostics.dat'
    if (.not.claim_series_sample(SERIES_NUMERICAL,path,integration_time)) return
    open(newunit=unit_file,file=trim(path),position='append',action='write')
    write(unit_file,'(ES24.16E3,1X,I0,16(1X,I0),5(1X,ES24.16E3))') &
      integration_time,n_steps,numerical_events,numerical_delta_u*dx*dy*dz
    close(unit_file)
    if (do_radial_mass_budget) call write_radial_mass_budget()
  end subroutine write_numerical_diagnostics

  subroutine initialize_text_file(file_name, header, replace_existing, series_id)
    character(len=*), intent(in) :: file_name, header
    logical, intent(in) :: replace_existing
    integer, intent(in), optional :: series_id
    logical :: file_exists
    integer :: unit_file

    inquire(file=trim(file_name), exist=file_exists)
    if (file_exists .and. .not. replace_existing) then
      if (present(series_id)) call load_series_tail(series_id,file_name)
      return
    end if

    open(newunit=unit_file, file=trim(file_name), status='replace', action='write')
    write(unit_file, '(A)') trim(header)
    close(unit_file)
    if (present(series_id)) then
      series_paths(series_id)=''
      call load_series_tail(series_id,file_name)
    end if
  end subroutine initialize_text_file


  subroutine load_series_tail(series_id,file_name)
    integer, intent(in) :: series_id
    character(len=*), intent(in) :: file_name
    integer :: unit_file, status
    real*8 :: sample_time, previous_time
    character(len=4096) :: line
    logical :: exists
    previous_time=-huge(1.0d0)
    inquire(file=trim(file_name),exist=exists)
    if (exists) then
      open(newunit=unit_file,file=trim(file_name),status='old',action='read')
      do
        read(unit_file,'(A)',iostat=status) line
        if (status < 0) exit
        if (status /= 0) error stop 'Cannot read diagnostic history'
        line=adjustl(line)
        if (len_trim(line)==0.or.line(1:1)=='#') cycle
        read(line,*,iostat=status) sample_time
        if (status /= 0) error stop 'Malformed diagnostic history'
        if (.not.ieee_is_finite(sample_time)) error stop 'Nonfinite diagnostic time'
        if (sample_time <= previous_time) error stop 'Diagnostic history is not strictly increasing'
        previous_time=sample_time
      end do
      close(unit_file)
    end if
    if (do_restart) then
      if (previous_time>integration_time) then
        write(*,'(A,1X,A)') 'ERROR history is ahead of checkpoint; preserve it and use another output folder:',trim(file_name)
        error stop 1
      end if
    end if
    series_paths(series_id)=file_name
    series_last_time(series_id)=previous_time
  end subroutine load_series_tail

  logical function claim_series_sample(series_id,file_name,sample_time)
    integer, intent(in) :: series_id
    character(len=*), intent(in) :: file_name
    real*8, intent(in) :: sample_time
    if (series_paths(series_id)/=file_name) call load_series_tail(series_id,file_name)
    if (sample_time < series_last_time(series_id)) error stop 'Diagnostic time went backwards'
    claim_series_sample=sample_time>series_last_time(series_id)
    if (claim_series_sample) series_last_time(series_id)=sample_time
  end function claim_series_sample

  subroutine write_run_manifest()
    integer :: unit_file
    character(len=512) :: file_name
    real*8 :: manifest_spin

    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/run_manifest.txt'
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(A)') '[run]'
    write(unit_file, '(A,ES24.16E3)') 'start_time=', integration_time
    write(unit_file, '(A,A)') 'case=', trim(case_name)
    write(unit_file, '(A,A)') 'metric=', trim(metric_type)
    write(unit_file, '(A,A)') 'geometry=', trim(geom_type)
    if (vtk_mapping_id == VTK_MAP_PHYSICAL) then
      write(unit_file, '(A)') 'vtk_mapping=physical'
    else
      write(unit_file, '(A)') 'vtk_mapping=untwisted'
    end if
    write(unit_file, '(A,ES24.16E3)') 'bh_mass=', bh_mass
    write(unit_file, '(A,ES24.16E3)') 'a_spin=', a_spin
    if (trim(metric_type) /= 'Minkowski' .and. bh_mass > 0.0d0) then
      manifest_spin = 0.0d0
      if (trim(metric_type) == 'Kerr-Schild') manifest_spin = a_spin
      write(unit_file, '(A,ES24.16E3)') 'event_horizon_radius=', &
                                       event_horizon_radius(bh_mass, manifest_spin)
    end if
    write(unit_file, '(A,A)') 'scheme=', trim(scheme_name)
    write(unit_file, '(A,3(I0,1X))') 'grid_dimensions=', nx, ny, nz
    write(unit_file, '(A,2(ES24.16E3,1X))') 'radial_domain=', r_min, r_max
    write(unit_file, '(A,L1)') 'logarithmic_radial_grid=', use_log_r
    write(unit_file, '(A,ES24.16E3)') 'final_time=', final_time
    write(unit_file, '(A,ES24.16E3)') 'cfl=', CFL
    write(unit_file, '(A,4(ES24.16E3,1X))') 'effective_floors_rho_p_D_tau=',rho_floor,p_floor,D_floor,tau_floor
    write(unit_file, '(A,3(ES24.16E3,1X))') 'walltime_limit_margin_checkpoint_seconds=', &
      walltime_limit_seconds,walltime_margin_seconds,walltime_checkpoint_seconds
    write(unit_file, '(A,L1)') 'ppi_diagnostics=', do_ppi_diagnostics
    write(unit_file, '(A,L1)') 'numerical_diagnostics=', do_numerical_diagnostics
    write(unit_file, '(A,L1)') 'shock_sensor=', use_shock_sensor
    write(unit_file, '(A)') 'shock_sensor_stencil=face_centered_i_minus_2_to_i_plus_3'
    write(unit_file, '(A,7(ES24.16E3,1X))') 'sensor_weno3_p_rho_tvd_p_rho_godunov_p_rho_atmosphere=', &
                                         shock_sensor_thresholds
    write(unit_file, '(A,I0)') 'diagnostic_stride=', diagnostic_stride
    write(unit_file, '(A,I0)') 'extraction_stride=', extraction_stride
    write(unit_file, '(A,I0)') 'mass_monitor_stride=', mass_monitor_stride
    if (do_gw_extraction) then
      write(unit_file, '(A)') 'gw_method=Finn-Evans weak-field stress proxy'
      write(unit_file, '(A)') 'gw_observer_axis=positive-z'
      write(unit_file, '(A,ES24.16E3)') 'gw_observer_distance=', gw_observer_distance
      write(unit_file, '(A,L1)') 'gw_exclude_inside_horizon=', gw_exclude_inside_horizon
      write(unit_file, '(A,L1)') 'gw_subtract_numerical_atmosphere=', &
                                 gw_subtract_numerical_atmosphere
      write(unit_file, '(A)') 'gw_volume_measure=flat Cartesian Kerr-Schild Jacobian'
      write(unit_file, '(A)') 'gw_density_measure=Newtonian rest-mass density proxy'
      write(unit_file, '(A,I0)') 'gw_polar_cells=', ny
    end if
    write(unit_file, '(A,L1)') 'apply_perturbation=', apply_perturbation
    write(unit_file, '(A,L1)') 'perturbation_already_applied=', perturbation_applied
    write(unit_file, '(A,I0)') 'perturbation_type=', perturbation_type
    write(unit_file, '(A,I0)') 'perturbation_seed=', perturbation_seed
    write(unit_file, '(A,ES24.16E3)') 'perturbation_time=', perturbation_time
    write(unit_file, '(A,ES24.16E3)') 'perturbation_amplitude=', perturbation_amplitude
    write(unit_file, '(A,ES24.16E3)') 'perturbation_mode=', perturbation_mode
    call write_initial_parameters_manifest(unit_file)
    write(unit_file, '(A)') '[/run]'
    close(unit_file)
  end subroutine write_run_manifest


  subroutine initialize_auxiliary_output(replace_existing)
    logical, intent(in) :: replace_existing
    integer :: exit_status
    character(len=512) :: directory, file_name

    series_paths=''
    series_last_time=-huge(1.0d0)
    directory = trim(output_folder) // '/' // trim(scheme_name)
    call execute_command_line('mkdir -p "' // trim(directory) // '"', &
                              exitstat=exit_status)
    if (exit_status /= 0) then
      write(*,*) 'CRITICAL ERROR: unable to create output directory: ', trim(directory)
      error stop 1
    end if

    if (do_gw_extraction) then
      file_name = trim(directory) // '/GW_signal.dat'
      call initialize_text_file(file_name, '# Time h_plus h_cross', replace_existing,SERIES_GW)
      if (ny == 1) then
        write(*,*) 'WARNING: GW amplitudes with ny=1 retain the single-cell polar quadrature.'
        write(*,*) '         Frequencies/phases are meaningful proxies; absolute 3D normalization is not fixed.'
      end if
    end if
    if (do_mdot_extraction) then
      file_name = trim(directory) // '/m_dot.dat'
      call initialize_text_file(file_name, '# Time M_dot', replace_existing,SERIES_MDOT)
    end if
    if (do_ppi_diagnostics) then
      file_name = trim(directory) // '/ppi_modes.dat'
      call initialize_text_file(file_name, &
        '# t C0 ReC1 ImC1 A1 phase1 ReC2 ImC2 A2 phase2 ReC3 ImC3 A3 phase3 ReC4 ImC4 A4 phase4', &
        replace_existing,SERIES_PPI)
      file_name = trim(directory) // '/global_diagnostics.dat'
      call initialize_text_file(file_name, &
        '# t M_total M_disk J_phi rho_min p_min rho_max p_max v2_max N_atmosphere N_invalid', &
        replace_existing,SERIES_GLOBAL)
      file_name = trim(directory) // '/perturbation_events.dat'
      call initialize_text_file(file_name, &
        '# t type seed amplitude mode', replace_existing,SERIES_EVENTS)
      if (trim(case_name) == 'KHI') then
        file_name = trim(directory) // '/khi_diagnostics.dat'
        call initialize_text_file(file_name, &
          '# t D_total tau_total Sx_total Sy_total Sz_total mode_plus mode_minus Kz_proxy', replace_existing,SERIES_KHI)
      end if
    end if

    call write_run_manifest()
  end subroutine initialize_auxiliary_output

  ! ===================================================================================
  ! 1. EXPORTADOR VTK (VisIt/ParaView) - FORMATO BINARIO
  ! ===================================================================================
  subroutine save_vtk_data(file_name)
    character(len=*), intent(in) :: file_name
    integer :: i, j, k, unit_file, io_status
    integer :: nx_points, ny_points, nz_points, n_points
    character(len=512) :: full_path
    character(len=100) :: header_line
    real*8, allocatable :: pts_buffer(:,:,:,:)
    real*8 :: r_phys, theta_point, phi_point
    real*8 :: geometry_spin, cartesian_position(3)
    real*8 :: topology_tolerance
    logical :: use_polar_unstructured_grid

    geometry_spin = 0.0d0
    if (trim(metric_type) == 'Kerr-Schild') geometry_spin = a_spin

    ! Un dominio angular polar completo no puede representarse con hexaedros
    ! estructurados regulares en los polos: allí distintos índices azimutales
    ! ocupan el mismo punto y VisIt puede producir nudos al cortar esas celdas.
    topology_tolerance = 128.0d0 * epsilon(1.0d0)
    use_polar_unstructured_grid = is_polar_geometry(geom_type) .and. &
      ny > 1 .and. nz > 1 .and. &
      abs(y_face(0)) <= topology_tolerance .and. &
      abs(y_face(ny) - pi) <= topology_tolerance * max(1.0d0, abs(pi)) .and. &
      abs((z_face(nz) - z_face(0)) - 2.0d0*pi) <= &
        topology_tolerance * max(1.0d0, abs(2.0d0*pi))

    if (use_polar_unstructured_grid) then
      call save_vtk_polar_unstructured(file_name, geometry_spin)
      return
    end if

    ! Las variables hidrodinámicas viven en celdas. La malla VTK se construye
    ! con sus caras verdaderas para que cada celda numérica corresponda a una
    ! celda geométrica y los polos no se formen con anillos artificiales.
    nx_points = nx + 1
    ny_points = 1
    nz_points = 1
    if (ny > 1) ny_points = ny + 1
    if (nz > 1) nz_points = nz + 1
    n_points = nx_points * ny_points * nz_points

    full_path = trim(output_folder) //  '/' // trim(scheme_name) // '/' // trim(file_name)
    print *, '>>> Guardando frame VTK Binario: ', trim(full_path)

    unit_file = 10
    open(unit_file, file=full_path, status='replace', access='stream', form='unformatted', convert='big_endian', iostat=io_status)
    if (io_status /= 0) then
      print *, 'CRITICAL ERROR: No se pudo crear el archivo VTK: ', trim(full_path)
      return
    end if

    write(unit_file) '# vtk DataFile Version 3.0' // char(10)
    write(header_line,'(A,ES24.16E3)') 'GRHD Simulation Data t=',integration_time
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) 'BINARY' // char(10)
    write(unit_file) 'DATASET STRUCTURED_GRID' // char(10)

    write(header_line, '(a, i0, 1x, i0, 1x, i0)') &
      'DIMENSIONS ', nx_points, ny_points, nz_points
    write(unit_file) trim(header_line) // char(10)

    write(header_line, '(a, i0, a)') 'POINTS ', n_points, ' double'
    write(unit_file) trim(header_line) // char(10)

    allocate(pts_buffer(3, 1:nx_points, 1:ny_points, 1:nz_points))

    ! -------------------------------------------------------------
    ! BUCLE DE MAPEO GEOMÉTRICO 3D (X, Y, Z)
    ! -------------------------------------------------------------
    do k = 1, nz_points
      if (nz > 1) then
        phi_point = z_face(k - 1)
      else
        phi_point = z(1)
      end if

      do j = 1, ny_points
        if (ny > 1) then
          theta_point = y_face(j - 1)
        else
          theta_point = y(1)
        end if

        do i = 1, nx_points
          if (use_log_r) then
            r_phys = exp(x_face(i - 1))
          else
            r_phys = x_face(i - 1)
          end if

          select case(trim(geom_type))
            case('Spherical', 'Spheroidal')
              call vtk_polar_position(r_phys, theta_point, phi_point, &
                                      geometry_spin, cartesian_position)
              pts_buffer(:, i, j, k) = cartesian_position

            case('Cylindrical')
              pts_buffer(1, i, j, k) = r_phys * cos(theta_point)
              pts_buffer(2, i, j, k) = r_phys * sin(theta_point)
              pts_buffer(3, i, j, k) = phi_point

            case default
              pts_buffer(1, i, j, k) = r_phys
              pts_buffer(2, i, j, k) = theta_point
              pts_buffer(3, i, j, k) = phi_point
          end select
        end do
      end do
    end do

    write(unit_file) pts_buffer
    write(unit_file) char(10)
    deallocate(pts_buffer)

    call write_vtk_cell_data(unit_file)
    close(unit_file)
  end subroutine save_vtk_data


  subroutine save_vtk_polar_unstructured(file_name, geometry_spin)
    character(len=*), intent(in) :: file_name
    real*8, intent(in) :: geometry_spin
    integer :: i, j, k, k_next, point_id, cell_id, offset
    integer :: unit_file, io_status, n_points, n_cells, cell_list_size
    integer(kind=4), allocatable :: ring_id(:,:,:), north_id(:), south_id(:)
    integer(kind=4), allocatable :: cell_buffer(:), cell_types(:)
    real*8, allocatable :: point_buffer(:,:)
    real*8 :: r_phys, position(3)
    character(len=512) :: full_path
    character(len=100) :: header_line

    ! Hay nz puntos por anillo: la conectividad periódica une el último con el
    ! primero. Cada cara radial tiene un único punto norte y uno sur.
    n_points = (nx + 1) * ((ny - 1) * nz + 2)
    n_cells = nx * ny * nz
    cell_list_size = nx * nz * (9 * ny - 4)

    allocate(point_buffer(3, n_points))
    allocate(ring_id(0:nx, 1:ny-1, 0:nz-1))
    allocate(north_id(0:nx), south_id(0:nx))

    point_id = 0
    do i = 0, nx
      if (use_log_r) then
        r_phys = exp(x_face(i))
      else
        r_phys = x_face(i)
      end if

      north_id(i) = point_id
      call vtk_polar_position(r_phys, 0.0d0, z_face(0), geometry_spin, position)
      point_buffer(:, point_id + 1) = position
      point_id = point_id + 1

      do k = 0, nz - 1
        do j = 1, ny - 1
          ring_id(i, j, k) = point_id
          call vtk_polar_position(r_phys, y_face(j), z_face(k), &
                                  geometry_spin, position)
          point_buffer(:, point_id + 1) = position
          point_id = point_id + 1
        end do
      end do

      south_id(i) = point_id
      call vtk_polar_position(r_phys, pi, z_face(0), geometry_spin, position)
      point_buffer(:, point_id + 1) = position
      point_id = point_id + 1
    end do

    if (point_id /= n_points) error stop 'Internal VTK point-count mismatch.'

    allocate(cell_buffer(cell_list_size), cell_types(n_cells))
    offset = 1
    cell_id = 0
    do k = 0, nz - 1
      k_next = mod(k + 1, nz)
      do j = 1, ny
        do i = 0, nx - 1
          cell_id = cell_id + 1
          if (j == 1) then
            ! Prisma polar norte: caras triangulares a radio constante.
            cell_buffer(offset:offset+6) = [6, north_id(i), &
              ring_id(i, 1, k), ring_id(i, 1, k_next), north_id(i+1), &
              ring_id(i+1, 1, k), ring_id(i+1, 1, k_next)]
            cell_types(cell_id) = 13       ! VTK_WEDGE
            offset = offset + 7
          else if (j == ny) then
            ! En el polo sur se invierte el orden azimutal para conservar la
            ! orientación de la celda.
            cell_buffer(offset:offset+6) = [6, south_id(i), &
              ring_id(i, ny-1, k_next), ring_id(i, ny-1, k), south_id(i+1), &
              ring_id(i+1, ny-1, k_next), ring_id(i+1, ny-1, k)]
            cell_types(cell_id) = 13       ! VTK_WEDGE
            offset = offset + 7
          else
            cell_buffer(offset:offset+8) = [8, ring_id(i, j-1, k), &
              ring_id(i+1, j-1, k), ring_id(i+1, j, k), ring_id(i, j, k), &
              ring_id(i, j-1, k_next), ring_id(i+1, j-1, k_next), &
              ring_id(i+1, j, k_next), ring_id(i, j, k_next)]
            cell_types(cell_id) = 12       ! VTK_HEXAHEDRON
            offset = offset + 9
          end if
        end do
      end do
    end do

    if (cell_id /= n_cells .or. offset - 1 /= cell_list_size) then
      error stop 'Internal VTK cell-count mismatch.'
    end if

    full_path = trim(output_folder) // '/' // trim(scheme_name) // '/' // trim(file_name)
    print *, '>>> Guardando frame VTK Binario: ', trim(full_path)
    unit_file = 10
    open(unit_file, file=full_path, status='replace', access='stream', &
         form='unformatted', convert='big_endian', iostat=io_status)
    if (io_status /= 0) then
      print *, 'CRITICAL ERROR: No se pudo crear el archivo VTK: ', trim(full_path)
      return
    end if

    write(unit_file) '# vtk DataFile Version 3.0' // char(10)
    write(header_line,'(A,ES24.16E3)') 'GRHD Simulation Data t=',integration_time
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) 'BINARY' // char(10)
    write(unit_file) 'DATASET UNSTRUCTURED_GRID' // char(10)
    write(header_line, '(a, i0, a)') 'POINTS ', n_points, ' double'
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) point_buffer
    write(unit_file) char(10)

    write(header_line, '(a, i0, 1x, i0)') 'CELLS ', n_cells, cell_list_size
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) cell_buffer
    write(unit_file) char(10)
    write(header_line, '(a, i0)') 'CELL_TYPES ', n_cells
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) cell_types
    write(unit_file) char(10)

    call write_vtk_cell_data(unit_file)
    close(unit_file)
    deallocate(point_buffer, ring_id, north_id, south_id, cell_buffer, cell_types)
  end subroutine save_vtk_polar_unstructured


  subroutine write_vtk_cell_data(unit_file)
    integer, intent(in) :: unit_file
    integer :: eq, i, j, k, n_cells
    character(len=100) :: header_line
    real*8, allocatable :: var_buffer(:,:,:)

    n_cells = nx * ny * nz
    write(header_line, '(a, i0)') 'CELL_DATA ', n_cells
    write(unit_file) trim(header_line) // char(10)
    allocate(var_buffer(1:nx, 1:ny, 1:nz))

    do eq = 1, neq
      write(header_line, '(a, 1x, a, a)') 'SCALARS', trim(var_names(eq)), ' double 1'
      write(unit_file) trim(header_line) // char(10)
      write(unit_file) 'LOOKUP_TABLE default' // char(10)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            var_buffer(i, j, k) = p(i,j,k,eq)
            if (use_log_r .and. eq == eq_vx) then
              var_buffer(i, j, k) = var_buffer(i, j, k) * exp(x(i))
            end if
          end do
        end do
      end do
      write(unit_file) var_buffer
      write(unit_file) char(10)
    end do

    write(unit_file) 'SCALARS Lapse_alpha double 1' // char(10)
    write(unit_file) 'LOOKUP_TABLE default' // char(10)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          var_buffer(i, j, k) = alpha_c(i, j)
        end do
      end do
    end do
    write(unit_file) var_buffer
    write(unit_file) char(10)

    write(unit_file) 'SCALARS Shift_beta_r double 1' // char(10)
    write(unit_file) 'LOOKUP_TABLE default' // char(10)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          var_buffer(i, j, k) = beta_c(1, i, j)
          if (use_log_r) var_buffer(i, j, k) = var_buffer(i, j, k) * exp(x(i))
        end do
      end do
    end do
    write(unit_file) var_buffer
    write(unit_file) char(10)

    deallocate(var_buffer)
  end subroutine write_vtk_cell_data


  subroutine save_vtk_at_time(current_time)
    real*8, intent(in) :: current_time
    character(len=512) :: file_name
    integer :: frame_id
    logical :: exists

    if (save_interval > 0.0d0) then
      frame_id = nint(current_time / save_interval)
    else
      frame_id = 0
    end if

    write(file_name, '(A, "_", A, "_", I4.4, ".vtk")') trim(output_prefix), trim(scheme_name), frame_id
    if (do_restart) then
      inquire(file=trim(output_folder)//'/'//trim(scheme_name)//'/'//trim(file_name),exist=exists)
      if (exists) error stop 'Restart would overwrite an existing VTK frame; use another output folder'
    end if
    call save_vtk_data(file_name)
  end subroutine save_vtk_at_time

  ! ===================================================================================
  ! 2. SISTEMA DE CHECKPOINTS (Gestión de memoria binaria directa)
  ! ===================================================================================

  ! The default path preserves every V9 metadata record; C01 requests canonical global controls.
  subroutine write_checkpoint_metadata(unit_file,io_status,global_shape,canonical)
    integer, intent(in) :: unit_file
    integer, intent(out) :: io_status
    integer, intent(in), optional :: global_shape(3)
    logical, intent(in), optional :: canonical
    character(len=12), parameter :: checkpoint_magic='GRHDCP_V9'
    integer :: sizes(3)
    integer(int64) :: events(NUM_EVENT_COUNT)
    real*8 :: delta_u(neq),radial_initial,radial_integral(2)
    logical :: saved_applied,radial_enabled
    sizes=[nx,ny,nz]
    if (present(global_shape)) sizes=global_shape
    saved_applied=perturbation_applied
    events=numerical_events;delta_u=numerical_delta_u
    radial_enabled=do_radial_mass_budget;radial_initial=radial_mass_initial;radial_integral=radial_mass_flux_integral
    if (present(canonical)) then
      if (canonical) then
        saved_applied=.false.;events=0_int64;delta_u=0.0d0
        radial_enabled=.false.;radial_initial=0.0d0;radial_integral=0.0d0
      end if
    end if
    write(unit_file,iostat=io_status) checkpoint_magic
    if (io_status /= 0) return
    ! Topología de malla 3D y profundidad de fantasmas
    write(unit_file,iostat=io_status) nghost, sizes
    if (io_status /= 0) return

    ! Fronteras lógicas (Actualizado a r_min/r_max)
    write(unit_file,iostat=io_status) r_min, r_max, y_min, y_max, z_min, z_max
    if (io_status /= 0) return

    ! Física y Banderas Arquitectónicas (Actualizado con banderas de extracción)
    write(unit_file,iostat=io_status) bh_mass, adb_idx, a_spin
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) use_log_r, use_shock_sensor, do_gw_extraction, &
                     do_mdot_extraction, do_ppi_diagnostics
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) metric_type, geom_type
    if (io_status /= 0) return

    ! Control de tiempo y salida
    write(unit_file,iostat=io_status) final_time, CFL, save_interval
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) output_prefix, output_folder
    if (io_status /= 0) return

    ! Identificadores de esquemas
    write(unit_file,iostat=io_status) case_id, case_name, scheme_name
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) rec_method_id, tvd_limiter_id, riemann_solver_id
    if (io_status /= 0) return

    ! V5 conserva los parámetros específicos de V4 y registra mediante su
    ! versión que el estado se almacena como [x,y,z,variable].
    ! fronteras. Son necesarios para reinicios reproducibles de Jet, Michel
    ! y Dust, y para documentar exactamente cualquier otro caso.
    call write_checkpoint_initial_parameters(unit_file,.true.,io_status)
    if (io_status /= 0) return

    ! Estado reproducible de la perturbación y frecuencia de diagnóstico
    write(unit_file,iostat=io_status) apply_perturbation, saved_applied
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) perturbation_type, perturbation_seed, diagnostic_stride
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) perturbation_time, perturbation_amplitude, perturbation_mode
    if (io_status /= 0) return
    ! V6 añade umbrales; evita reanudar una corrida B silenciosamente con A.
    write(unit_file,iostat=io_status) shock_sensor_thresholds
    if (io_status /= 0) return

    ! V8: effective atmosphere and runtime controls, independent of physical t.
    write(unit_file,iostat=io_status) rho_floor, p_floor, D_floor, tau_floor
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) checkpoint_interval, extraction_stride, mass_monitor_stride, &
      do_numerical_diagnostics, vtk_mapping_id
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) walltime_limit_seconds, walltime_margin_seconds, walltime_checkpoint_seconds
    if (io_status /= 0) return

    write(unit_file,iostat=io_status) events,delta_u
    if (io_status /= 0) return
    write(unit_file,iostat=io_status) radial_enabled,radial_initial,radial_integral
    if (io_status /= 0) return

  end subroutine write_checkpoint_metadata

  subroutine save_checkpoint(step_num, current_time, status)
    integer, intent(in) :: step_num
    real*8, intent(in) :: current_time
    integer, intent(out), optional :: status
    integer :: unit_file, io_status, close_status
    character(len=512) :: file_name, temporary_name, message
    character(len=24) :: operation

    if (present(status)) status = 0
    write(file_name, '(A, "/checkpoint_", A, "_step_", I9.9, ".rst")') &
      trim(output_folder), trim(scheme_name), step_num
    temporary_name=trim(file_name)//'.tmp'
    print *, '>>> Guardando Checkpoint de Respaldo: ', trim(file_name)
    operation='open'
    message=''
    open(newunit=unit_file, file=trim(temporary_name), status='replace', action='write', &
      form='unformatted', iostat=io_status,iomsg=message)
    if (io_status /= 0) goto 900
    operation='write metadata'
    call write_checkpoint_metadata(unit_file,io_status)
    if (io_status /= 0) goto 800

    operation='write data'
    ! Estado temporal de la simulación
    write(unit_file,iostat=io_status,iomsg=message) current_time
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status,iomsg=message) step_num
    if (io_status /= 0) goto 800

    ! Vuelco binario masivo del arreglo de estado conservado 4D
    write(unit_file,iostat=io_status,iomsg=message) up
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status,iomsg=message) 'GRIAL_END_V9', step_num, current_time, nx, ny, nz, neq
    if (io_status /= 0) goto 800
800 continue
    if (io_status == 0) then
      operation='close'
      close(unit_file,iostat=io_status,iomsg=message)
    else
      close(unit_file,iostat=close_status)
    end if
    if (io_status /= 0) goto 900
    operation='verify footer'
    call verify_checkpoint_footer(trim(temporary_name),'GRIAL_END_V9',step_num,current_time,io_status,message)
    if (io_status /= 0) goto 900
    ! Same-directory rename publishes a complete generation; older steps survive.
    operation='rename'
    message=''
    call rename(trim(temporary_name),trim(file_name),io_status)
    if (io_status == 0) return
900 continue
    write(error_unit,'(5A,I0,2A)') 'ERROR checkpoint not published: ',trim(operation),' ', &
      trim(temporary_name),' iostat=',io_status,' ',trim(message)
    if (present(status)) then
      status=io_status
    else
      error stop 1
    end if
  end subroutine save_checkpoint

  subroutine verify_checkpoint_footer(path, expected_footer, step, time, io_status, message)
    character(len=*), intent(in) :: path, expected_footer
    integer, intent(in) :: step
    real*8, intent(in) :: time
    integer, intent(out) :: io_status
    character(len=*), intent(out) :: message
    integer :: unit_file, saved_step, saved_shape(4), close_status
    real*8 :: saved_time
    character(len=12) :: footer

    ! Verify the final record from physical EOF. A forward sequential read
    ! may accept a truncated trailing record marker with GNU Fortran.
    message=''
    open(newunit=unit_file,file=path,status='old',action='read', &
      form='unformatted',position='append',iostat=io_status,iomsg=message)
    if (io_status /= 0) return
    backspace(unit_file,iostat=io_status,iomsg=message)
    if (io_status == 0) read(unit_file,iostat=io_status,iomsg=message) footer,saved_step,saved_time,saved_shape
    if (io_status == 0) then
      if (footer/=expected_footer.or.saved_step/=step.or. &
          transfer(saved_time,0_int64)/=transfer(time,0_int64).or.any(saved_shape/=[nx,ny,nz,neq])) then
        io_status=1
        message='checkpoint footer disagrees with written state'
      end if
    end if
    if (io_status == 0) then
      close(unit_file,iostat=io_status,iomsg=message)
    else
      close(unit_file,iostat=close_status)
    end if
  end subroutine verify_checkpoint_footer

  subroutine read_checkpoint_metadata(unit_file, status)
    integer, intent(in) :: unit_file
    integer, intent(out), optional :: status
    character(len=160) :: failure_message
    integer :: io_status
    logical :: is_v2, is_v3, is_v4, is_v5, is_v6, is_v7, is_v8, is_v9, is_versioned
    character(len=12) :: checkpoint_magic
    character(len=12), parameter :: checkpoint_magic_v2 = 'GRHDCP_V2'
    character(len=12), parameter :: checkpoint_magic_v3 = 'GRHDCP_V3'
    character(len=12), parameter :: checkpoint_magic_v4 = 'GRHDCP_V4'
    character(len=12), parameter :: checkpoint_magic_v5 = 'GRHDCP_V5'
    character(len=12), parameter :: checkpoint_magic_v6 = 'GRHDCP_V6'
    character(len=12), parameter :: checkpoint_magic_v7 = 'GRHDCP_V7'

    character(len=12), parameter :: checkpoint_magic_v8 = 'GRHDCP_V8'

    character(len=12), parameter :: checkpoint_magic_v9 = 'GRHDCP_V9'

    if (present(status)) status=0
    checkpoint_magic = ''
    read(unit_file, iostat=io_status) checkpoint_magic
    is_v2 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v2)
    is_v3 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v3)
    is_v4 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v4)
    is_v5 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v5)
    is_v6 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v6)
    is_v7 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v7)
    is_v8 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v8)
    is_v9 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v9)
    checkpoint_restores_diagnostics=is_v9
    checkpoint_has_footer = is_v8.or.is_v9
    checkpoint_footer_magic='GRIAL_END_V8'
    if (is_v9) checkpoint_footer_magic='GRIAL_END_V9'
    if (index(checkpoint_magic,'GRHDCP_') == 1 .and. .not.  (is_v2.or.is_v3.or.is_v4.or.is_v5.or.is_v6.or.is_v7.or.is_v8.or.is_v9)) then
      failure_message='Unsupported checkpoint version';goto 880
    end if
    is_versioned = is_v2 .or. is_v3 .or. is_v4 .or. is_v5 .or. is_v6 .or. is_v7 .or. is_v8 .or. is_v9
    checkpoint_uses_soa_layout = is_v5 .or. is_v6 .or. is_v7 .or. is_v8 .or. is_v9
    if (.not. is_versioned) rewind(unit_file)

    ! Topología de malla
    read(unit_file,iostat=io_status) nghost, nx, ny, nz
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
    end if
    if (min(nx,ny,nz,nghost) < 1) then
      failure_message='Invalid checkpoint topology';goto 880

    end if
    ! Fronteras lógicas
    read(unit_file,iostat=io_status) r_min, r_max, y_min, y_max, z_min, z_max
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
    end if
    ! Física y Banderas Arquitectónicas (Actualizado con banderas de extracción)
    if (is_versioned) then
      read(unit_file,iostat=io_status) bh_mass, adb_idx, a_spin
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
    else
      read(unit_file,iostat=io_status) bh_mass, adb_idx
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
      a_spin = 0.0d0
      print *, 'WARNING: legacy checkpoint has no Kerr spin; assuming a = 0.'
    end if
    if (is_v3 .or. is_v4 .or. is_v5 .or. is_v6 .or. is_v7 .or. is_v8 .or. is_v9) then
      read(unit_file,iostat=io_status) use_log_r, use_shock_sensor, do_gw_extraction,  do_mdot_extraction, do_ppi_diagnostics
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
    else
      read(unit_file,iostat=io_status) use_log_r, use_shock_sensor, do_gw_extraction, do_mdot_extraction
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
      do_ppi_diagnostics = .false.
    end if
    read(unit_file,iostat=io_status) metric_type, geom_type
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
    end if
    ! Control
    read(unit_file,iostat=io_status) final_time, CFL, save_interval
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
    end if
    read(unit_file,iostat=io_status) output_prefix, output_folder
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
    end if
    ! Eschemas
    read(unit_file,iostat=io_status) case_id, case_name, scheme_name
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
    end if
    ! IDs historicos 10--12 siempre representan FD, incluso con nombre FM antiguo.
    select case(case_id)
    case(10); case_name='FontDaigneEqu'
    case(11); case_name='FontDaigneSag'
    case(12); case_name='FontDaigne3D'
    case(13:15)
      if (.not.(is_v7.or.is_v8.or.is_v9)) then
      failure_message='FM requires a V7 checkpoint with separate FD/FM parameters';goto 880
      end if
    end select

    read(unit_file,iostat=io_status) rec_method_id, tvd_limiter_id, riemann_solver_id
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
    end if
    if (is_v4 .or. is_v5 .or. is_v6 .or. is_v7 .or. is_v8 .or. is_v9) then
      call read_checkpoint_initial_parameters(unit_file,is_v7.or.is_v8.or.is_v9,io_status)
      if (io_status/=0) then
        failure_message='Incomplete checkpoint initial parameters';goto 880
      end if
    else if (case_id == 5 .or. case_id == 7 .or. case_id == 8) then
      write(*,'(A)') 'WARNING: legacy checkpoint lacks case-specific boundary parameters; using defaults.'
    end if

    if (is_v3 .or. is_v4 .or. is_v5 .or. is_v6 .or. is_v7 .or. is_v8 .or. is_v9) then
      read(unit_file,iostat=io_status) apply_perturbation, perturbation_applied
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
      read(unit_file,iostat=io_status) perturbation_type, perturbation_seed, diagnostic_stride
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
      read(unit_file,iostat=io_status) perturbation_time, perturbation_amplitude, perturbation_mode
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
    else
      perturbation_applied = .false.
    end if
    shock_sensor_thresholds = SENSOR_DEFAULTS
    if (is_v6 .or. is_v7 .or. is_v8 .or. is_v9) then
      read(unit_file,iostat=io_status) shock_sensor_thresholds
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
    end if
    if (is_v8.or.is_v9) then
      read(unit_file,iostat=io_status) rho_floor, p_floor, D_floor, tau_floor
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
      read(unit_file,iostat=io_status) checkpoint_interval, extraction_stride, mass_monitor_stride,  do_numerical_diagnostics, vtk_mapping_id
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
      read(unit_file,iostat=io_status) walltime_limit_seconds, walltime_margin_seconds, walltime_checkpoint_seconds
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint metadata';goto 880
      end if
    else
      rho_floor=1.0d-10; p_floor=1.0d-13
      if (.not.ieee_is_finite(adb_idx)) then
      failure_message='Invalid checkpoint adiabatic index';goto 880
      end if
      if (adb_idx <= 1.0d0) then
      failure_message='Invalid checkpoint adiabatic index';goto 880
      end if
      D_floor=rho_floor; tau_floor=p_floor*(adb_idx/(adb_idx-1.0d0)-1.0d0)
      write(*,'(A)') 'WARNING: legacy checkpoint lacks floors: assuming rho=1e-10, p=1e-13; override explicitly if needed.'
    end if
    numerical_events=0_int64; numerical_delta_u=0.0d0
    do_radial_mass_budget=.false.; radial_mass_initial=0.0d0; radial_mass_flux_integral=0.0d0
    if (is_v9) then
      read(unit_file,iostat=io_status) numerical_events,numerical_delta_u
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint diagnostic counters';goto 880
      end if
      read(unit_file,iostat=io_status) do_radial_mass_budget,radial_mass_initial,radial_mass_flux_integral
      if (io_status /= 0) then
      failure_message='Incomplete checkpoint radial budget';goto 880
      end if
      if (any(numerical_events < 0_int64)) then
      failure_message='Invalid checkpoint counters';goto 880
      end if
      if (.not.all(ieee_is_finite(numerical_delta_u))) then
      failure_message='Invalid checkpoint repair budget';goto 880
      end if
      if (.not.all(ieee_is_finite([radial_mass_initial,radial_mass_flux_integral]))) then
      failure_message='Invalid checkpoint radial budget';goto 880
      end if
    else
      write(*,'(A)') 'WARNING: pre-V9 checkpoint lacks diagnostic accumulators; starting a new diagnostic baseline.'
    end if
    return
880 continue
    if (present(status)) then
      status=1
    else
      write(*,'(A)')trim(failure_message)
      error stop 1
    end if
  end subroutine read_checkpoint_metadata

  subroutine read_checkpoint_state(unit_file, step_num, current_time)
    integer, intent(in) :: unit_file
    integer, intent(out) :: step_num
    real*8, intent(out) :: current_time
    real*8, allocatable :: legacy_up(:,:,:,:)
    integer :: component, i, j, k, io_status, saved_step, dims(4)
    real*8 :: saved_time
    character(len=12) :: footer
    character(len=1024) :: path
    character(len=512) :: message

    read(unit_file,iostat=io_status) current_time
    if (io_status /= 0) error stop 'Incomplete checkpoint state'
    read(unit_file,iostat=io_status) step_num
    if (io_status /= 0) error stop 'Incomplete checkpoint state'

    ! El programa principal debe asegurar que 'up' ya esté alojado en memoria
    ! dinámicamente con nx, ny, nz leídos en la metadata antes de llegar aquí.
    if (checkpoint_uses_soa_layout) then
      read(unit_file,iostat=io_status) up
      if (io_status /= 0) error stop 'Incomplete checkpoint state'
    else
      allocate(legacy_up(neq,nx,ny,nz))
      read(unit_file,iostat=io_status) legacy_up
      if (io_status /= 0) error stop 'Incomplete checkpoint state'
      do component = 1, neq
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              up(i,j,k,component) = legacy_up(component,i,j,k)
            end do
          end do
        end do
      end do
      deallocate(legacy_up)
    end if
    if (.not.ieee_is_finite(current_time)) error stop 'Invalid checkpoint time'
    if (current_time < 0.0d0 .or. step_num < 0) error stop 'Invalid checkpoint time/step'
    if (.not.all(ieee_is_finite(up))) error stop 'Nonfinite checkpoint state'
    if (checkpoint_has_footer) then
      read(unit_file,iostat=io_status) footer,saved_step,saved_time,dims
      if (io_status /= 0) error stop 'Incomplete checkpoint footer'
      if (footer /= checkpoint_footer_magic .or. saved_step /= step_num .or. saved_time /= current_time .or. &
          any(dims /= [nx,ny,nz,neq])) error stop 'Inconsistent checkpoint footer'
      read(unit_file,iostat=io_status)
      if (io_status >= 0) error stop 'Unexpected checkpoint trailing record'
      inquire(unit=unit_file,name=path,iostat=io_status)
      if (io_status /= 0) error stop 'Cannot identify checkpoint file'
      call verify_checkpoint_footer(trim(path),checkpoint_footer_magic,step_num,current_time,io_status,message)
      if (io_status /= 0) then
        write(error_unit,'(4A)') 'Invalid checkpoint footer: ',trim(path),' ',trim(message)
        error stop 1
      end if
    end if
  end subroutine read_checkpoint_state

  subroutine write_checkpoint_initial_parameters(unit_file,include_fm,status)
    integer, intent(in) :: unit_file
    logical, intent(in), optional :: include_fm

    integer, intent(out), optional :: status
    integer :: io_status

    if (present(status)) status=0
    write(unit_file,iostat=io_status) sod_interface, sod_rho_left, sod_rho_right, &
      sod_pressure_left, sod_pressure_right
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) strong_interface, strong_rho_left, strong_rho_right, &
      strong_pressure_left, strong_pressure_right
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) shu_interface, shu_rho_left, shu_pressure_left, &
      shu_rho_right, shu_rho_amplitude, shu_wave_number, shu_pressure_right
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) khi_half_width, khi_rho_inner, khi_rho_outer, &
      khi_vx_inner, khi_vx_outer, khi_pressure, khi_perturbation_amplitude, &
      khi_perturbation_wave_number, khi_perturbation_width
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) jet_ambient_density, jet_ambient_pressure, &
      jet_ambient_velocity, jet_nozzle_radius, jet_nozzle_length, &
      jet_density, jet_pressure, jet_velocity
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) advected_wave_density, advected_wave_amplitude, &
      advected_wave_speed, advected_wave_pressure, advected_wave_number
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) michel_critical_radius, michel_critical_density, &
      dust_accretion_constant
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) offaxis_radial_center, offaxis_phi_center, offaxis_width, &
      offaxis_background_density, offaxis_background_pressure, &
      offaxis_density_amplitude, offaxis_pressure_amplitude
    if (io_status /= 0) goto 800
    write(unit_file,iostat=io_status) fd_inner_radius, fd_pressure_max_radius, &
      fd_polytropic_constant
    if (io_status /= 0) goto 800

    if (present(include_fm)) then
      if (include_fm) then
        write(unit_file,iostat=io_status) fm_inner_radius,fm_pressure_max_radius,fm_polytropic_constant
        if (io_status /= 0) goto 800
      end if
    end if
    return
800 continue
    if (present(status)) then
      status=io_status
    else
      error stop 'Checkpoint initial parameters write failed'
    end if
  end subroutine write_checkpoint_initial_parameters


  subroutine read_checkpoint_initial_parameters(unit_file, include_fm, status)
    integer, intent(in) :: unit_file
    logical, intent(in), optional :: include_fm

    integer, intent(out), optional :: status
    character(len=160) :: failure_message
    integer :: io_status

    if (present(status)) status=0
    read(unit_file,iostat=io_status) sod_interface, sod_rho_left, sod_rho_right,  sod_pressure_left, sod_pressure_right
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) strong_interface, strong_rho_left, strong_rho_right,  strong_pressure_left, strong_pressure_right
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) shu_interface, shu_rho_left, shu_pressure_left,  shu_rho_right, shu_rho_amplitude, shu_wave_number, shu_pressure_right
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) khi_half_width, khi_rho_inner, khi_rho_outer,  khi_vx_inner, khi_vx_outer, khi_pressure, khi_perturbation_amplitude,  khi_perturbation_wave_number, khi_perturbation_width
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) jet_ambient_density, jet_ambient_pressure,  jet_ambient_velocity, jet_nozzle_radius, jet_nozzle_length,  jet_density, jet_pressure, jet_velocity
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) advected_wave_density, advected_wave_amplitude,  advected_wave_speed, advected_wave_pressure, advected_wave_number
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) michel_critical_radius, michel_critical_density,  dust_accretion_constant
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) offaxis_radial_center, offaxis_phi_center, offaxis_width,  offaxis_background_density, offaxis_background_pressure,  offaxis_density_amplitude, offaxis_pressure_amplitude
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    read(unit_file,iostat=io_status) fd_inner_radius, fd_pressure_max_radius,  fd_polytropic_constant
    if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
    end if
    fm_inner_radius=8.0d0;fm_pressure_max_radius=12.0d0;fm_polytropic_constant=0.0015d0
    if (present(include_fm)) then
      if (include_fm) then
        read(unit_file,iostat=io_status) fm_inner_radius,fm_pressure_max_radius,fm_polytropic_constant
        if (io_status /= 0) then
      failure_message='Incomplete checkpoint initial parameters';goto 880
        end if
      end if
    end if
    return
880 continue
    if (present(status)) then
      status=1
    else
      write(*,'(A)')trim(failure_message)
      error stop 1
    end if
  end subroutine read_checkpoint_initial_parameters


  subroutine write_initial_parameters_manifest(unit_file)
    integer, intent(in) :: unit_file

    select case(case_id)
    case(1)
      call manifest_real(unit_file, 'sod_interface', sod_interface)
      call manifest_real(unit_file, 'sod_density_left', sod_rho_left)
      call manifest_real(unit_file, 'sod_density_right', sod_rho_right)
      call manifest_real(unit_file, 'sod_pressure_left', sod_pressure_left)
      call manifest_real(unit_file, 'sod_pressure_right', sod_pressure_right)
    case(2)
      call manifest_real(unit_file, 'strong_interface', strong_interface)
      call manifest_real(unit_file, 'strong_density_left', strong_rho_left)
      call manifest_real(unit_file, 'strong_density_right', strong_rho_right)
      call manifest_real(unit_file, 'strong_pressure_left', strong_pressure_left)
      call manifest_real(unit_file, 'strong_pressure_right', strong_pressure_right)
    case(3)
      call manifest_real(unit_file, 'shu_interface', shu_interface)
      call manifest_real(unit_file, 'shu_density_left', shu_rho_left)
      call manifest_real(unit_file, 'shu_pressure_left', shu_pressure_left)
      call manifest_real(unit_file, 'shu_density_right', shu_rho_right)
      call manifest_real(unit_file, 'shu_density_amplitude', shu_rho_amplitude)
      call manifest_real(unit_file, 'shu_wave_number', shu_wave_number)
      call manifest_real(unit_file, 'shu_pressure_right', shu_pressure_right)
    case(4)
      call manifest_real(unit_file, 'khi_half_width', khi_half_width)
      call manifest_real(unit_file, 'khi_density_inner', khi_rho_inner)
      call manifest_real(unit_file, 'khi_density_outer', khi_rho_outer)
      call manifest_real(unit_file, 'khi_velocity_inner', khi_vx_inner)
      call manifest_real(unit_file, 'khi_velocity_outer', khi_vx_outer)
      call manifest_real(unit_file, 'khi_pressure', khi_pressure)
      call manifest_real(unit_file, 'khi_perturbation_amplitude', khi_perturbation_amplitude)
      call manifest_real(unit_file, 'khi_perturbation_wave_number', khi_perturbation_wave_number)
      call manifest_real(unit_file, 'khi_perturbation_width', khi_perturbation_width)
    case(5)
      call manifest_real(unit_file, 'jet_ambient_density', jet_ambient_density)
      call manifest_real(unit_file, 'jet_ambient_pressure', jet_ambient_pressure)
      call manifest_real(unit_file, 'jet_ambient_velocity', jet_ambient_velocity)
      call manifest_real(unit_file, 'jet_nozzle_radius', jet_nozzle_radius)
      call manifest_real(unit_file, 'jet_nozzle_length', jet_nozzle_length)
      call manifest_real(unit_file, 'jet_density', jet_density)
      call manifest_real(unit_file, 'jet_pressure', jet_pressure)
      call manifest_real(unit_file, 'jet_velocity', jet_velocity)
    case(6)
      call manifest_real(unit_file, 'convergence_density', advected_wave_density)
      call manifest_real(unit_file, 'convergence_amplitude', advected_wave_amplitude)
      call manifest_real(unit_file, 'convergence_speed', advected_wave_speed)
      call manifest_real(unit_file, 'convergence_pressure', advected_wave_pressure)
      call manifest_real(unit_file, 'convergence_wave_number', advected_wave_number)
    case(7)
      call manifest_real(unit_file, 'michel_critical_radius', michel_critical_radius)
      call manifest_real(unit_file, 'michel_critical_density', michel_critical_density)
    case(8)
      call manifest_real(unit_file, 'dust_accretion_constant', dust_accretion_constant)
    case(9)
      call manifest_real(unit_file, 'offaxis_radial_center', offaxis_radial_center)
      call manifest_real(unit_file, 'offaxis_phi_center', offaxis_phi_center)
      call manifest_real(unit_file, 'offaxis_width', offaxis_width)
      call manifest_real(unit_file, 'offaxis_background_density', offaxis_background_density)
      call manifest_real(unit_file, 'offaxis_background_pressure', offaxis_background_pressure)
      call manifest_real(unit_file, 'offaxis_density_amplitude', offaxis_density_amplitude)
      call manifest_real(unit_file, 'offaxis_pressure_amplitude', offaxis_pressure_amplitude)
    case(10:12)
      write(unit_file,'(A)') 'torus_model=FD'
      write(unit_file,'(A)') 'torus_rotation_law=constant_minus_u_phi_over_u_t'
      call manifest_real(unit_file, 'fd_inner_radius', fd_inner_radius)
      call manifest_real(unit_file, 'fd_pressure_max_radius', fd_pressure_max_radius)
      call manifest_real(unit_file, 'fd_polytropic_constant', fd_polytropic_constant)
    case(13:15)
      write(unit_file,'(A)') 'torus_model=FM'
      write(unit_file,'(A)') 'torus_rotation_law=constant_u_phi_times_u_t_contravariant'
      call manifest_real(unit_file,'fm_inner_radius',fm_inner_radius)
      call manifest_real(unit_file,'fm_pressure_max_radius',fm_pressure_max_radius)
      call manifest_real(unit_file,'fm_polytropic_constant',fm_polytropic_constant)
    end select
  end subroutine write_initial_parameters_manifest


  subroutine manifest_real(unit_file, key, value)
    integer, intent(in) :: unit_file
    character(len=*), intent(in) :: key
    real*8, intent(in) :: value

    write(unit_file,'(A,A,ES24.16E3)') trim(key), '=', value
  end subroutine manifest_real

  ! ===================================================================================
  ! 3. EXTRACCIÓN DE OBSERVABLES FÍSICOS (Física de Agujeros Negros)
  ! ===================================================================================

  pure subroutine kerr_schild_spheroidal_position(r_phys, theta, phi, spin, position)
    real*8, intent(in) :: r_phys, theta, phi, spin
    real*8, intent(out) :: position(3)
    real*8 :: sin_theta, cos_theta, sin_phi, cos_phi

    sin_theta = sin(theta)
    cos_theta = cos(theta)

    ! sin(pi) no es exactamente cero en punto flotante. En Kerr-Schild ese
    ! residuo genera un anillo diminuto cuya orientación depende de r y spin;
    ! VisIt puede mostrarlo como un nudo en el polo. Colapsamos ambos ejes de
    ! coordenadas de manera exacta.
    if (abs(sin_theta) <= 64.0d0 * epsilon(1.0d0)) then
      position(1) = 0.0d0
      position(2) = 0.0d0
      position(3) = sign(r_phys, cos_theta)
      return
    end if

    sin_phi = sin(phi)
    cos_phi = cos(phi)

    position(1) = (r_phys*cos_phi - spin*sin_phi) * sin_theta
    position(2) = (r_phys*sin_phi + spin*cos_phi) * sin_theta
    position(3) = r_phys*cos_theta
  end subroutine kerr_schild_spheroidal_position


  subroutine vtk_polar_position(r_phys, theta, phi, spin, position)
    real*8, intent(in) :: r_phys, theta, phi, spin
    real*8, intent(out) :: position(3)
    real*8 :: sin_theta, cos_theta, cylindrical_radius

    if (vtk_mapping_id == VTK_MAP_PHYSICAL) then
      call kerr_schild_spheroidal_position(r_phys, theta, phi, spin, position)
      return
    end if

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    if (abs(sin_theta) <= 64.0d0 * epsilon(1.0d0)) then
      position = [0.0d0, 0.0d0, sign(r_phys, cos_theta)]
      return
    end if

    ! Visualizacion oblata sin la rotacion radial atan(a/r). Conserva las
    ! superficies r=constante y colapsa exactamente ambos polos, pero no cambia
    ! las coordenadas ni los datos usados durante la evolucion.
    cylindrical_radius = sqrt(r_phys**2 + spin**2) * sin_theta
    position(1) = cylindrical_radius * cos(phi)
    position(2) = cylindrical_radius * sin(phi)
    position(3) = r_phys * cos_theta
  end subroutine vtk_polar_position

  pure subroutine gw_cartesian_kinematics(r_phys, theta, phi, spin, alpha, beta, &
                                          eulerian_velocity, logarithmic_r, &
                                          position, velocity, flat_jacobian)
    real*8, intent(in) :: r_phys, theta, phi, spin, alpha
    real*8, intent(in) :: beta(3), eulerian_velocity(3)
    logical, intent(in) :: logarithmic_r
    real*8, intent(out) :: position(3), velocity(3), flat_jacobian
    real*8 :: coordinate_rates(3)
    real*8 :: sin_theta, cos_theta, sin_phi, cos_phi, sigma
    real*8 :: r_dot, theta_dot, phi_dot

    ! Las primitivas son velocidades eulerianas. Finn--Evans necesita las
    ! tasas coordenadas para construir velocidades cartesianas.
    coordinate_rates = alpha * eulerian_velocity - beta
    if (logarithmic_r) coordinate_rates(1) = r_phys * coordinate_rates(1)

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_phi = sin(phi)
    cos_phi = cos(phi)
    r_dot = coordinate_rates(1)
    theta_dot = coordinate_rates(2)
    phi_dot = coordinate_rates(3)

    ! Para spin=0 se recupera la transformación esférica usada por EF.
    call kerr_schild_spheroidal_position(r_phys, theta, phi, spin, position)

    velocity(1) = r_dot * cos_phi * sin_theta &
                + theta_dot * (r_phys * cos_phi - spin * sin_phi) * cos_theta &
                - phi_dot * (r_phys * sin_phi + spin * cos_phi) * sin_theta
    velocity(2) = r_dot * sin_phi * sin_theta &
                + theta_dot * (r_phys * sin_phi + spin * cos_phi) * cos_theta &
                + phi_dot * (r_phys * cos_phi - spin * sin_phi) * sin_theta
    velocity(3) = r_dot * cos_theta - r_phys * theta_dot * sin_theta

    sigma = r_phys**2 + spin**2 * cos_theta**2
    flat_jacobian = sigma * sin_theta
  end subroutine gw_cartesian_kinematics


  pure subroutine finn_evans_cell_stress(rho, pressure, central_mass_radius3, &
                                         position, velocity, cell_volume, &
                                         Ixx_cell, Iyy_cell, Ixy_cell)
    real*8, intent(in) :: rho, pressure, central_mass_radius3, cell_volume
    real*8, intent(in) :: position(3), velocity(3)
    real*8, intent(out) :: Ixx_cell, Iyy_cell, Ixy_cell

    Ixx_cell = (2.0d0 * rho * velocity(1)**2 + 2.0d0 * pressure &
               -2.0d0 * rho * central_mass_radius3 * position(1)**2) * cell_volume
    Iyy_cell = (2.0d0 * rho * velocity(2)**2 + 2.0d0 * pressure &
               -2.0d0 * rho * central_mass_radius3 * position(2)**2) * cell_volume
    Ixy_cell = (2.0d0 * rho * velocity(1) * velocity(2) &
               -2.0d0 * rho * central_mass_radius3 * position(1) * position(2)) * cell_volume
  end subroutine finn_evans_cell_stress

  subroutine calc_gw_strain(current_time,local_stress)
    real*8, intent(in) :: current_time
    real*8, intent(out), optional :: local_stress(3)
    integer :: i, j, k
    real*8 :: Ixx_ddot, Iyy_ddot, Ixy_ddot
    real*8 :: r, th, phi, rho, pres, dV, mass_radius3
    real*8 :: h_plus, h_cross
    real*8 :: alpha, beta(3), eulerian_velocity(3)
    real*8 :: position(3), cartesian_velocity(3), flat_jacobian
    real*8 :: cartesian_radius2, radial_grid_jacobian, polar_grid_jacobian
    real*8 :: horizon_radius, geometry_spin
    real*8 :: Ixx_cell, Iyy_cell, Ixy_cell

    if (present(local_stress)) local_stress=0.0d0
    if (trim(metric_type) == 'Minkowski' .or. bh_mass < 0.1d0) return

    ! EF es el caso Schwarzschild a=0. KS usa las coordenadas oblatas
    ! asociadas con el spin configurado en la métrica.
    geometry_spin = 0.0d0
    if (trim(metric_type) == 'Kerr-Schild') geometry_spin = a_spin
    horizon_radius = event_horizon_radius(bh_mass, geometry_spin)

    Ixx_ddot = 0.0d0
    Iyy_ddot = 0.0d0
    Ixy_ddot = 0.0d0

    !$OMP PARALLEL DO PRIVATE(i, j, k, r, th, phi, rho, pres, dV, mass_radius3, &
    !$OMP alpha, beta, eulerian_velocity, position, &
    !$OMP cartesian_velocity, flat_jacobian, cartesian_radius2, radial_grid_jacobian) &
    !$OMP PRIVATE(Ixx_cell, Iyy_cell, Ixy_cell, polar_grid_jacobian) &
    !$OMP REDUCTION(+:Ixx_ddot, Iyy_ddot, Ixy_ddot)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (use_log_r) then
            r = exp(x(i))
          else
            r = x(i)
          end if

          ! El proxy representa únicamente la materia del disco exterior. La
          ! stress formula newtoniana no modela la desconexión causal del
          ! horizonte; incluir celdas interiores daría una fuente instantánea
          ! manifiestamente no física.
          if (gw_exclude_inside_horizon .and. r <= horizon_radius) cycle

          th = polar_angle(y(j))
          phi = z(k)

          rho = p(i,j,k,eq_de)
          pres = p(i,j,k,eq_pr)
          if (gw_subtract_numerical_atmosphere) then
            rho = max(0.0d0, rho - rho_floor)
            pres = max(0.0d0, pres - p_floor)
          end if
          if (rho <= 0.0d0 .and. pres <= 0.0d0) cycle

          alpha = alpha_c(i,j)
          beta = beta_c(:,i,j)
          eulerian_velocity = p(i,j,k,eq_vx:eq_vz)
          polar_grid_jacobian=polar_jacobian(y(j))
          if (polar_mesh_h/=1.0d0) then
            beta(2)=beta(2)*polar_grid_jacobian
            eulerian_velocity(2)=eulerian_velocity(2)*polar_grid_jacobian
          end if
          call gw_cartesian_kinematics(r, th, phi, geometry_spin, alpha, beta, &
                                       eulerian_velocity, use_log_r, position, &
                                       cartesian_velocity, flat_jacobian)

          ! d^3X usa el Jacobiano cartesiano plano de las coordenadas KS,
          ! tal como requiere la aproximación cuadrupolar de campo débil.
          radial_grid_jacobian = 1.0d0
          if (use_log_r) radial_grid_jacobian = r
          dV = flat_jacobian * radial_grid_jacobian * dx * dy * dz
          if (polar_mesh_h/=1.0d0) dV=dV*polar_grid_jacobian

          ! A orden newtoniano el potencial exterior sólo contiene el
          ! monopolo de masa. El spin afecta la dinámica GRHD y la geometría
          ! de coordenadas, pero no se inventa aquí un potencial escalar de
          ! frame dragging.
          cartesian_radius2 = sum(position**2)
          mass_radius3 = bh_mass / (cartesian_radius2 * sqrt(cartesian_radius2))

          call finn_evans_cell_stress(rho, pres, mass_radius3, position, &
                                      cartesian_velocity, dV, Ixx_cell, &
                                      Iyy_cell, Ixy_cell)
          Ixx_ddot = Ixx_ddot + Ixx_cell
          Iyy_ddot = Iyy_ddot + Iyy_cell
          Ixy_ddot = Ixy_ddot + Ixy_cell
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    if (present(local_stress)) then
      local_stress=[Ixx_ddot,Iyy_ddot,Ixy_ddot]
      return
    end if

    h_plus  = (Ixx_ddot - Iyy_ddot) / gw_observer_distance
    h_cross = (2.0d0 * Ixy_ddot) / gw_observer_distance

    block
      character(len=512) :: gw_filename
      gw_filename = trim(output_folder) //  '/' // trim(scheme_name) // '/GW_signal.dat'
      if (.not.claim_series_sample(SERIES_GW,gw_filename,current_time)) return
      open(30, file=gw_filename, position='append', status='unknown')
      write(30, '(3(ES24.16E3,1X))') current_time, h_plus, h_cross
      close(30)
    end block
  end subroutine calc_gw_strain


  subroutine calc_m_dot(current_time,local_mdot)
    implicit none
    real*8, intent(in) :: current_time
    real*8, intent(out), optional :: local_mdot
    real*8 :: m_dot
    real*8 :: alpha, beta(3), v_eff
    integer :: j, k

    if (present(local_mdot)) local_mdot=0.0d0
    ! Corrección: Actualización de bandera
    if (trim(metric_type) == 'Minkowski' .or. bh_mass < 0.1d0) return
    if (present(local_mdot).and.idx_probe==0) return
    if (idx_probe < 1 .or. idx_probe > nx) then
      write(*,*) 'CRITICAL ERROR: invalid accretion-probe index: ', idx_probe
      error stop 1
    end if

    m_dot = 0.0d0
    ! Tasa de acreción integrada sobre una esfera a un radio fijo (idx_probe)
    ! Se propaga la coordenada polar j
    !$OMP PARALLEL DO REDUCTION(+:m_dot) PRIVATE(j, k, alpha, beta, v_eff)
    do k = 1, nz
      do j = 1, ny
        ! call calculate_metric(x(idx_probe), y(j), alpha=alpha, beta=beta)
        alpha = alpha_c(idx_probe, j)
        beta(:) = beta_c(:, idx_probe, j)

        v_eff = alpha * p(idx_probe,j,k,eq_vx) - beta(1)

        ! 'up' contiene el estado conservado densificado. up(eq_de) es D = rho * W * sqrt(gamma)
        m_dot = m_dot + up(idx_probe,j,k,eq_de) * v_eff * dy * dz
      end do
    end do
    !$OMP END PARALLEL DO

    m_dot = - m_dot
    if (present(local_mdot)) then
      local_mdot=m_dot
      return
    end if

    block
      character(len=512) :: m_dot_filename
      m_dot_filename = trim(output_folder) //  '/' // trim(scheme_name) // '/m_dot.dat'
      if (.not.claim_series_sample(SERIES_MDOT,m_dot_filename,current_time)) return
      open(88, file=m_dot_filename, position='append', status='unknown')
      write(88, '(2(ES24.16E3,1X))') current_time, m_dot
      close(88)
    end block
  end subroutine calc_m_dot


  subroutine calc_ppi_modes(current_time,local_coefficients)
    implicit none
    real*8, intent(in) :: current_time
    real*8, intent(out), optional :: local_coefficients(1+2*PPI_MAX_MODE)
    integer :: i, j, k, mode_m, unit_file
    real*8 :: cell_volume, density_cut, azimuthal_mass
    real*8 :: c0, c_real(PPI_MAX_MODE), c_imag(PPI_MAX_MODE)
    real*8 :: values(2 + 4 * PPI_MAX_MODE)
    character(len=512) :: file_name

    if (present(local_coefficients)) local_coefficients=0.0d0
    if (.not. do_ppi_diagnostics .or. nz < 2) return

    cell_volume = dx * dy * dz
    density_cut = 1.0d3 * rho_floor
    c0 = 0.0d0
    c_real = 0.0d0
    c_imag = 0.0d0

    ! Primero se integra la masa de cada sector azimutal. Esto evita evaluar
    ! funciones trigonométricas dentro del bucle tridimensional completo.
    do k = 1, nz
      azimuthal_mass = 0.0d0
      do j = 1, ny
        do i = 1, nx
          if (p(i,j,k,eq_de) > density_cut) then
            azimuthal_mass = azimuthal_mass + up(i,j,k,eq_de) * cell_volume
          end if
        end do
      end do
      c0 = c0 + azimuthal_mass
      do mode_m = 1, PPI_MAX_MODE
        c_real(mode_m) = c_real(mode_m) + azimuthal_mass * cos(dble(mode_m) * z(k))
        c_imag(mode_m) = c_imag(mode_m) - azimuthal_mass * sin(dble(mode_m) * z(k))
      end do
    end do

    if (present(local_coefficients)) then
      local_coefficients=[c0,c_real,c_imag]
      return
    end if
    call ppi_mode_values(current_time,[c0,c_real,c_imag],values)

    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/ppi_modes.dat'
    if (.not.claim_series_sample(SERIES_PPI,file_name,current_time)) return
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(*(ES24.16E3,1X))') values
    close(unit_file)
  end subroutine calc_ppi_modes


  pure subroutine ppi_mode_values(current_time,coefficients,values)
    real*8, intent(in) :: current_time,coefficients(1+2*PPI_MAX_MODE)
    real*8, intent(out) :: values(2+4*PPI_MAX_MODE)
    real*8 :: c0,c_real,c_imag
    integer :: m,base
    c0=coefficients(1);values=0.0d0;values(1)=current_time;values(2)=c0
    do m=1,PPI_MAX_MODE
      base=3+4*(m-1);c_real=coefficients(1+m);c_imag=coefficients(1+PPI_MAX_MODE+m)
      values(base)=c_real;values(base+1)=c_imag
      if (c0>tiny(1.0d0)) then
        values(base+2)=sqrt(c_real**2+c_imag**2)/c0
        values(base+3)=-atan2(c_imag,c_real)/dble(m)
      end if
    end do
  end subroutine ppi_mode_values

  subroutine calc_global_diagnostics(current_time)
    implicit none
    real*8, intent(in) :: current_time
    integer :: i, j, k, ii, jj, unit_file
    integer(kind=8) :: atmosphere_cells, invalid_cells
    real*8 :: cell_volume, density_cut, rho, pres, v_sq
    real*8 :: mass_total, disk_mass, angular_momentum
    real*8 :: rho_minimum, pressure_minimum, rho_maximum, pressure_maximum
    real*8 :: max_velocity_squared
    character(len=512) :: file_name

    if (.not. do_ppi_diagnostics) return

    cell_volume = dx * dy * dz
    density_cut = 1.0d3 * rho_floor
    mass_total = 0.0d0
    disk_mass = 0.0d0
    angular_momentum = 0.0d0
    rho_minimum = huge(1.0d0)
    pressure_minimum = huge(1.0d0)
    rho_maximum = -huge(1.0d0)
    pressure_maximum = -huge(1.0d0)
    max_velocity_squared = 0.0d0
    atmosphere_cells = 0_8
    invalid_cells = 0_8

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho = p(i,j,k,eq_de)
          pres = p(i,j,k,eq_pr)
          v_sq = 0.0d0
          do ii = 1, 3
            do jj = 1, 3
              v_sq = v_sq + gamma_c(ii, jj, i, j) * &
                     p(i,j,k,eq_vx + ii - 1) * p(i,j,k,eq_vx + jj - 1)
            end do
          end do

          if (.not. ieee_is_finite(rho) .or. .not. ieee_is_finite(pres) .or. &
              .not. ieee_is_finite(v_sq) .or. rho <= 0.0d0 .or. pres <= 0.0d0 .or. &
              v_sq < 0.0d0 .or. v_sq >= 1.0d0) then
            invalid_cells = invalid_cells + 1_8
          end if

          if (ieee_is_finite(rho)) then
            rho_minimum = min(rho_minimum, rho)
            rho_maximum = max(rho_maximum, rho)
          end if
          if (ieee_is_finite(pres)) then
            pressure_minimum = min(pressure_minimum, pres)
            pressure_maximum = max(pressure_maximum, pres)
          end if
          if (ieee_is_finite(v_sq)) max_velocity_squared = max(max_velocity_squared, v_sq)

          mass_total = mass_total + up(i,j,k,eq_de) * cell_volume
          angular_momentum = angular_momentum + up(i,j,k,eq_vz) * cell_volume
          if (rho > density_cut) then
            disk_mass = disk_mass + up(i,j,k,eq_de) * cell_volume
          else
            atmosphere_cells = atmosphere_cells + 1_8
          end if
        end do
      end do
    end do

    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/global_diagnostics.dat'
    if (.not.claim_series_sample(SERIES_GLOBAL,file_name,current_time)) return
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(9(ES24.16E3,1X),2(I0,1X))') current_time, mass_total, disk_mass, &
      angular_momentum, rho_minimum, pressure_minimum, rho_maximum, pressure_maximum, &
      max_velocity_squared, atmosphere_cells, invalid_cells
    close(unit_file)
    if (trim(case_name) == 'KHI') call calc_khi_diagnostics(current_time)
  end subroutine calc_global_diagnostics

  ! Caja periódica SRHD: balances de las variables realmente evolucionadas
  ! y amplitud del modo sembrado en cada interfaz original. Kz es un proxy
  ! transversal, no una descomposición covariante de la energía cinética.
  subroutine calc_khi_diagnostics(current_time)
    real*8, intent(in) :: current_time
    integer :: i,k,m,unit_file
    real*8 :: totals(neq),volume,phase,cx(nx),sx(nx),wp(nz),wm(nz)
    real*8 :: rp,ip,rm,im,vz,rho,pres,v2,Kz,enthalpy,amp_plus,amp_minus
    character(len=512) :: file_name
    if (trim(case_name) /= 'KHI' .or. ny /= 1 .or. trim(metric_type) /= 'Minkowski') return
    volume=dx*dy*dz
    totals=0.0d0; rp=0.0d0; ip=0.0d0; rm=0.0d0; im=0.0d0; Kz=0.0d0
    do i=1,nx
      phase=khi_perturbation_wave_number*pi*x(i)
      cx(i)=cos(phase); sx(i)=sin(phase)
    end do
    do k=1,nz
      wp(k)=exp(-0.5d0*((z(k)-khi_half_width)/khi_perturbation_width)**2)
      wm(k)=exp(-0.5d0*((z(k)+khi_half_width)/khi_perturbation_width)**2)
      do i=1,nx
        do m=1,neq
          totals(m)=totals(m)+up(i,1,k,m)*volume
        end do
        vz=p(i,1,k,eq_vz)
        rp=rp+vz*cx(i)*wp(k); ip=ip-vz*sx(i)*wp(k)
        rm=rm+vz*cx(i)*wm(k); im=im-vz*sx(i)*wm(k)
        rho=p(i,1,k,eq_de); pres=p(i,1,k,eq_pr)
        v2=sum(p(i,1,k,eq_vx:eq_vz)**2)
        enthalpy=rho+g1*pres
        Kz=Kz+0.5d0*enthalpy*vz*vz/max(1.0d0-v2,tiny(1.0d0))*volume
      end do
    end do
    amp_plus=2.0d0*sqrt(rp*rp+ip*ip)/(dble(nx)*sum(wp))
    amp_minus=2.0d0*sqrt(rm*rm+im*im)/(dble(nx)*sum(wm))
    file_name=trim(output_folder)//'/'//trim(scheme_name)//'/khi_diagnostics.dat'
    if (.not.claim_series_sample(SERIES_KHI,file_name,current_time)) return
    open(newunit=unit_file,file=trim(file_name),status='unknown',position='append',action='write')
    write(unit_file,'(9(ES24.16E3,1X))') current_time,totals(eq_de),totals(eq_pr), &
      totals(eq_vx),totals(eq_vy),totals(eq_vz),amp_plus,amp_minus,Kz
    close(unit_file)
  end subroutine calc_khi_diagnostics


  subroutine write_perturbation_event(current_time)
    implicit none
    real*8, intent(in) :: current_time
    integer :: unit_file
    character(len=512) :: file_name

    if (.not. do_ppi_diagnostics) return
    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/perturbation_events.dat'
    if (.not.claim_series_sample(SERIES_EVENTS,file_name,current_time)) return
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(ES24.16E3,1X,2(I0,1X),2(ES24.16E3,1X))') current_time, &
      perturbation_type, perturbation_seed, perturbation_amplitude, perturbation_mode
    close(unit_file)
  end subroutine write_perturbation_event


  subroutine calc_convergence_norms()
    implicit none
    integer :: i, i_start, i_end, j_idx, unit_out, idx_max_err
    real*8 :: rho_exact, dust_exact(neq)
    real*8 :: err_L1, err_L2, sum_L1, sum_L2
    real*8 :: err_Linf, diff, x_max_err
    real*8 :: sinc_factor, sinc_argument, celdas_validas, advected_coordinate

    real*8, allocatable :: p_conv(:,:,:,:)
    character(len=512) :: file_name
    logical :: file_exists

    if (trim(case_name) /= 'Conv' .and. trim(case_name) /= 'MichelA' .and. &
        trim(case_name) /= 'Dust') return

    sum_L1 = 0.0d0
    sum_L2 = 0.0d0
    err_Linf = 0.0d0
    idx_max_err = 1
    j_idx = max(1, ny / 2)
    x_max_err = x(1)

    if (trim(case_name) == 'Conv') then
      ! Solución exacta en promedios celulares. Las fronteras son periódicas.
      sinc_argument = 0.5d0 * advected_wave_number * pi * dx
      if (abs(sinc_argument) > epsilon(1.0d0)) then
        sinc_factor = sin(sinc_argument) / sinc_argument
      else
        sinc_factor = 1.0d0
      end if
      celdas_validas = dble(nx)
      do i = 1, nx
        advected_coordinate = x(i) - advected_wave_speed * integration_time
        rho_exact = advected_wave_density + advected_wave_amplitude * sinc_factor * &
                    sin(advected_wave_number * pi * advected_coordinate)
        diff = abs(p(i,j_idx,1,eq_de) - rho_exact)
        sum_L1 = sum_L1 + diff
        sum_L2 = sum_L2 + diff**2
        if (diff > err_Linf) then
          err_Linf = diff
          idx_max_err = i
          x_max_err = x(i)
        end if
      end do
      file_name = trim(output_folder) // '/convergence_' // trim(scheme_name) // '.dat'
    else if (trim(case_name) == 'MichelA') then
      ! Michel conserva la comparación histórica en la región interior segura.
      allocate(p_conv(1:nx,1:ny,1:nz,neq))
      call setup_michel_accretion_initial(p_conv, michel_critical_radius, &
                                           michel_critical_density)
      i_start = max(1, int(0.10d0 * dble(nx)))
      i_end = min(nx, max(i_start, int(0.95d0 * dble(nx))))
      celdas_validas = dble(i_end - i_start + 1)
      do i = i_start, i_end
        rho_exact = p_conv(i,j_idx,1,eq_de)
        diff = abs(p(i,j_idx,1,eq_de) - rho_exact)
        sum_L1 = sum_L1 + diff
        sum_L2 = sum_L2 + diff**2
        if (diff > err_Linf) then
          err_Linf = diff
          idx_max_err = i
          if (use_log_r) then
            x_max_err = exp(x(i))
          else
            x_max_err = x(i)
          end if
        end if
      end do
      deallocate(p_conv)
      file_name = trim(output_folder) // '/michel_convergence_' // trim(scheme_name) // '.dat'
    else
      ! La solución fría de caída libre también es estacionaria. Se excluyen
      ! franjas junto a ambos bordes para medir el error del volumen interior.
      i_start = max(1, int(0.10d0 * dble(nx)))
      i_end = min(nx, max(i_start, int(0.95d0 * dble(nx))))
      celdas_validas = dble(i_end - i_start + 1)
      do i = i_start, i_end
        call evaluate_dust_accretion_state(x(i), dust_exact)
        rho_exact = dust_exact(eq_de)
        diff = abs(p(i,j_idx,1,eq_de) - rho_exact)
        sum_L1 = sum_L1 + diff
        sum_L2 = sum_L2 + diff**2
        if (diff > err_Linf) then
          err_Linf = diff
          idx_max_err = i
          if (use_log_r) then
            x_max_err = exp(x(i))
          else
            x_max_err = x(i)
          end if
        end if
      end do
      file_name = trim(output_folder) // '/dust_convergence_' // trim(scheme_name) // '.dat'
    end if

    err_L1 = sum_L1 / celdas_validas
    err_L2 = sqrt(sum_L2 / celdas_validas)

    print *, "N =", nx, "| L1 =", err_L1, "| L2 =", err_L2, &
             "| Linf =", err_Linf, " localizado en x =", x_max_err, &
             " indice = ", idx_max_err

    unit_out = 20

    ! ====================================================================
    ! LÓGICA DE ESCRITURA INTELIGENTE
    ! ====================================================================
    inquire(file=trim(file_name), exist=file_exists)
    open(unit=unit_out, file=trim(file_name), status='unknown', position='append')
    if (.not. file_exists) write(unit_out, '(A)') '# N L1 L2 Linf'

    write(unit_out, '(I6, 3(2X, E20.12))') nx, err_L1, err_L2, err_Linf
    close(unit_out)

  end subroutine calc_convergence_norms

end module output
