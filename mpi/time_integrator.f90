! Shared SSP-RK3. Distributed callers initialize the CPU halo exchange first.
module time_integrator
  use variables
  use equations, only: cons_to_prim
  use evolution, only: calc_rhs
  use conditions, only: set_boundary_conditions,set_boundary_conditions_accelerator
  use output, only: accumulate_radial_mass_budget
  use halo_exchange, only: fill_primitive_halos
  use parallel_diagnostics, only: diagnostics_active,begin_budget_step,accumulate_rhs_budget, &
    accumulate_recovery_budget,end_budget_step
  use parallel_runtime, only: parallel_require_active,parallel_require_collective,parallel_size
  use parallel_timing, only: validate_step_size
  implicit none
  private
  public :: advance_rk_step

contains

  subroutine advance_rk_step(distributed)
    logical, intent(in), optional :: distributed
    logical :: exchange_halos,state_ok
    integer :: i,j,k,rk
    real*8 :: recovery_before(neq),recovery_change(neq)

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
    call refresh_stage_primitives(exchange_halos)

    ! Integrador Runge-Kutta Fuerte de Preservación de Estabilidad (SSP)
    do rk = 1, 3

      ! A. Integración Espacial: Reconstrucción, Flujos y Términos Fuente
      call calc_rhs(p, rhs)
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

      if (diagnostics_active) call accumulate_recovery_budget(rk,recovery_change)

      ! D. Actualización de Topología de Red
      call refresh_stage_primitives(exchange_halos)

    end do

    if (diagnostics_active) call end_budget_step()

  end subroutine advance_rk_step

  subroutine refresh_stage_primitives(exchange_halos)
    logical, intent(in) :: exchange_halos
    if (exchange_halos) then
      call fill_primitive_halos(p)
    else if (accelerator_enabled) then
      call set_boundary_conditions_accelerator(p)
    else
      call set_boundary_conditions(p)
    end if
  end subroutine refresh_stage_primitives
end module time_integrator
