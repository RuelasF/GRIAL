! Shared SSP-RK3. Distributed callers initialize the CPU halo exchange first.
module time_integrator
  use variables
  use equations, only: cons_to_prim
  use evolution, only: calc_rhs
#ifdef GRIAL_OVERLAP
  use evolution, only: prepare_overlap_sources,prepare_overlap_polar_fluxes
#endif
  use conditions, only: set_boundary_conditions,set_boundary_conditions_accelerator
  use output, only: accumulate_radial_mass_budget
  use halo_exchange, only: fill_primitive_halos
#ifdef GRIAL_OVERLAP
  use halo_exchange, only: begin_radial_overlap,end_radial_overlap,finish_polar_overlap
#endif
  use parallel_diagnostics, only: diagnostics_active,begin_budget_step,accumulate_rhs_budget, &
    accumulate_recovery_budget,end_budget_step
  use parallel_runtime, only: parallel_require_active,parallel_require_collective,parallel_size
#ifdef GRIAL_PROFILE
  use parallel_runtime, only: profile_clock,profile_add
#endif
  use parallel_timing, only: validate_step_size
  implicit none
  private
  public :: advance_rk_step

contains

  subroutine advance_rk_step(distributed)
    logical, intent(in), optional :: distributed
    logical :: exchange_halos,state_ok
#ifdef GRIAL_OVERLAP
    logical :: overlap_started,polar_started
#endif
    integer :: i,j,k,rk
    real*8 :: recovery_before(neq),recovery_change(neq)
#ifdef GRIAL_PROFILE
    real*8 :: profile_start
#endif

    exchange_halos=.false.
    if (present(distributed)) exchange_halos=distributed
    if (exchange_halos) then
      call parallel_require_active()
      call parallel_require_collective(.not.accelerator_enabled,'distributed RK supports CPU only')
      call validate_step_size(dt)
      state_ok=allocated(u).and.allocated(up).and.allocated(rhs).and.allocated(p)
      if (state_ok) state_ok=all(shape(u)==[nx,ny,nz,neq]).and.all(shape(up)==[nx,ny,nz,neq]).and. &
        all(shape(rhs)==[nx,ny,nz,neq]).and. &
        all(shape(p)==[nx+2*nghost+1,ny+2*nghost+1,nz+2*nghost+1,neq])
      call parallel_require_collective(state_ok,'distributed RK state arrays do not match the local grid')
      ! Local endpoint fluxes cannot be published as a global Michel budget.
      call parallel_require_collective(.not.do_radial_mass_budget.or.parallel_size==1, &
        'distributed radial mass budget requires D01 collector rather than serial output')
    end if

    if (diagnostics_active) then
      call parallel_require_collective(exchange_halos.or.parallel_size==1,'distributed diagnostics require distributed RK')
      call begin_budget_step(dt)
    end if

#ifdef GRIAL_OVERLAP
    overlap_started=.false.
    if (exchange_halos) call begin_radial_overlap(p,overlap_started)
#endif
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(u,up)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          u(i,j,k,:)=up(i,j,k,:)
        end do
      end do
    end do
    !$ACC END PARALLEL LOOP
    !$OMP END PARALLEL DO
#ifdef GRIAL_OVERLAP
    if (overlap_started) then
      call prepare_overlap_sources(p)
      call end_radial_overlap(p,polar_started)
      if (polar_started) then
        call prepare_overlap_polar_fluxes(p)
        call finish_polar_overlap(p)
      end if
    else
      call refresh_stage_primitives(exchange_halos)
    end if
#else
    call refresh_stage_primitives(exchange_halos)
#endif

    ! Integrador Runge-Kutta Fuerte de Preservación de Estabilidad (SSP)
    do rk = 1, 3

      ! A. Integración Espacial: Reconstrucción, Flujos y Términos Fuente
#ifdef GRIAL_PROFILE
      profile_start=profile_clock()
#endif
      call calc_rhs(p, rhs)
#ifdef GRIAL_PROFILE
      call profile_add(2,profile_start)
#endif
      if (diagnostics_active) call accumulate_rhs_budget(rk)
      if (do_radial_mass_budget) call accumulate_radial_mass_budget(rk)

      ! B. Avance Temporal Low-Storage
      ! (Las operaciones matriciales vectorizadas funcionan sin bucles for)
      !$OMP PARALLEL DO PRIVATE(i,j,k)
      !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(u,up,rhs)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            if (rk == 1) then
              up(i,j,k,:)=u(i,j,k,:)+rhs(i,j,k,:)*dt
            else if (rk == 2) then
              up(i,j,k,:)=0.75d0*u(i,j,k,:)+0.25d0*(up(i,j,k,:)+rhs(i,j,k,:)*dt)
            else
              up(i,j,k,:)=(u(i,j,k,:)+2.0d0*(up(i,j,k,:)+rhs(i,j,k,:)*dt))/3.0d0
            end if
          end do
        end do
      end do
      !$ACC END PARALLEL LOOP
      !$OMP END PARALLEL DO

      ! C. Recuperación Termodinámica (Conservativas -> Primitivas)
#ifdef GRIAL_PROFILE
      profile_start=profile_clock()
#endif
      recovery_change=0.0d0
      !$OMP PARALLEL DO PRIVATE(i,j,k,recovery_before) REDUCTION(+:recovery_change)
      !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(up,p,gamma_c,gamma_inv_c,sqrt_gamma_c)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
#ifndef USE_OPENACC
            if (diagnostics_active) recovery_before=up(i,j,k,:)
#endif
            call cons_to_prim(up(i,j,k,:), p(i,j,k,:), i, j, k)
#ifndef USE_OPENACC
            if (diagnostics_active) recovery_change=recovery_change+(up(i,j,k,:)-recovery_before)
#endif
          end do
        end do
      end do
      !$ACC END PARALLEL LOOP
      !$OMP END PARALLEL DO

#ifdef GRIAL_PROFILE
      call profile_add(3,profile_start)
#endif
      if (diagnostics_active) call accumulate_recovery_budget(rk,recovery_change)

      ! D. Actualización de Topología de Red
#ifdef GRIAL_OVERLAP
      overlap_started=.false.
      if (exchange_halos.and.rk<3) call begin_radial_overlap(p,overlap_started)
      if (overlap_started) then
        call prepare_overlap_sources(p)
        call end_radial_overlap(p,polar_started)
      if (polar_started) then
        call prepare_overlap_polar_fluxes(p)
        call finish_polar_overlap(p)
      end if
      else
        call refresh_stage_primitives(exchange_halos)
      end if
#else
      call refresh_stage_primitives(exchange_halos)
#endif

    end do

    if (diagnostics_active) call end_budget_step()

  end subroutine advance_rk_step

  subroutine refresh_stage_primitives(exchange_halos)
    logical, intent(in) :: exchange_halos
#ifdef GRIAL_PROFILE
    real*8 :: profile_start
    profile_start=profile_clock()
#endif
    if (exchange_halos) then
      call fill_primitive_halos(p)
    else if (accelerator_enabled) then
      call set_boundary_conditions_accelerator(p)
    else
      call set_boundary_conditions(p)
    end if
#ifdef GRIAL_PROFILE
    call profile_add(4,profile_start)
#endif
  end subroutine refresh_stage_primitives
end module time_integrator
