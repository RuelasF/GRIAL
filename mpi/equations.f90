! =======================================================================================
! Módulo: equations (Física y Transformaciones de Estado GRHD)
! ---------------------------------------------------------------------------------------
! Este módulo es el corazón hidrodinámico del código. En GRHD, a diferencia de la física
! Newtoniana, las variables primitivas (rho, p, v) y conservativas (D, S, tau) están
! acopladas de forma altamente no lineal por el factor de Lorentz (W) y la entalpía (h).
! Este módulo gestiona estas conversiones, calcula los flujos a través de las interfaces
! numéricas y evalúa los términos fuente debidos a la geometría y curvatura espacial.
! =======================================================================================
module equations
  use variables
  use metrics
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: prim_to_cons, cons_to_prim, calc_fluxes, calc_sources

contains

  ! ===================================================================================
  ! TRANSFORMACIONES DIRECTAS (Primitivas -> Conservativas)
  ! ===================================================================================

  ! Subrutina: prim_to_cons
  ! Propósito: Mapeo directo P -> U. Es puramente algebraico (sin iteraciones).
  ! Fórmulas analíticas (Densidades Tensoriales conservadas):
  ! D   = rho * W * sqrt(gamma)
  ! S_i = rho * h * W^2 * v_i * sqrt(gamma)
  ! tau = (rho * h * W^2 - p - rho * W) * sqrt(gamma)

  subroutine prim_to_cons(prim_state, cons_state, i, j, k)
    !$acc routine seq
    implicit none
    real*8, intent(in)  :: prim_state(:)
    integer, intent(in) :: i, j, k
    real*8, intent(out) :: cons_state(:)

    real*8 :: v_sq, h, W
    real*8 :: v(3), v_cov(3)
    real*8 :: g(3,3), sqg
    integer :: ii, jj

    ! 1. LEER MÉTRICA DE LA CACHÉ
    g(:,:) = gamma_c(:,:,i,j)
    sqg    = sqrt_gamma_c(i,j)

    v(1) = prim_state(eq_vx)
    v(2) = prim_state(eq_vy)
    v(3) = prim_state(eq_vz)

    ! Cuadrado de la velocidad Euleriana tridimensional general
    v_sq = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_sq = v_sq + g(ii,jj) * v(ii) * v(jj)
      end do
    end do

    ! --- ESCUDO DE CAUSALIDAD NUMÉRICA ---
    if (v_sq >= v_max) then
#ifndef USE_OPENACC
      call count_numerical_event(NUM_P2C_VELOCITY_CAP)
#endif
      v(1) = v(1) * sqrt(v_max / v_sq)
      v(2) = v(2) * sqrt(v_max / v_sq)
      v(3) = v(3) * sqrt(v_max / v_sq)
      v_sq = v_max
    end if

    ! Calcular la velocidad covariante (v_i = gamma_ij * v^j)
    v_cov = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_cov(ii) = v_cov(ii) + g(ii,jj) * v(jj)
      end do
    end do

    ! Entalpía específica (gas ideal) y Factor de Lorentz
    h = 1.0d0 + g1 * prim_state(eq_pr) / prim_state(eq_de)
    W = 1.0d0 / sqrt(1.0d0 - v_sq)

    ! Cálculo de variables conservativas
    cons_state(eq_de) = prim_state(eq_de) * W
    cons_state(eq_pr) = prim_state(eq_de) * h * W**2 - prim_state(eq_pr) - prim_state(eq_de) * W

    ! Momentos espaciales covariantes S_i = rho * h * W^2 * v_i
    cons_state(eq_vx) = prim_state(eq_de) * h * W**2 * v_cov(1)
    cons_state(eq_vy) = prim_state(eq_de) * h * W**2 * v_cov(2)
    cons_state(eq_vz) = prim_state(eq_de) * h * W**2 * v_cov(3)

    ! Transformación a densidades respecto a la malla coordenada
    cons_state = cons_state * sqg
  end subroutine prim_to_cons


  ! ===================================================================================
  ! TRANSFORMACIONES INVERSAS (Conservativas -> Primitivas)
  ! ===================================================================================
  subroutine cons_to_prim(cons_state, prim_state, i, j, k)
    !$acc routine seq
    implicit none
    real*8, intent(inout) :: cons_state(:)
    integer, intent(in)   :: i, j, k
    real*8, intent(out)   :: prim_state(:)

    real*8 :: W, V, dV, v_sq_final
    real*8 :: E_tot, S_sq, S_sq_max, scale_factor
    real*8 :: local_D_floor, local_tau_floor
    real*8 :: g(3,3), ginv(3,3), sqg
    real*8 :: S_cov(3)
    integer :: ii, jj, repair_event
    logical :: finite_input(neq)
    ! Binary64 exponent mask also detects signaling NaNs without raising INVALID.
    integer(int64), parameter :: exponent_mask=int(z'7FF0000000000000',int64)
#ifndef USE_OPENACC
    real*8 :: diagnostic_before(neq)
    if (do_numerical_diagnostics) then
      diagnostic_before=cons_state
      call count_numerical_event(NUM_C2P_CALLS)
    end if
#endif

    ! 1. LEER MÉTRICA DE LA CACHÉ
    g(:,:)    = gamma_c(:,:,i,j)
    ginv(:,:) = gamma_inv_c(:,:,i,j)
    sqg       = sqrt_gamma_c(i,j)

    local_D_floor   = rho_floor * sqg
    local_tau_floor = p_floor * (g1 - 1.0d0) * sqg

    repair_event=NUM_C2P_SANITIZE
    do ii=1,neq
      finite_input(ii)=iand(transfer(cons_state(ii),0_int64),exponent_mask)/=exponent_mask
    end do
    if (.not.all(finite_input)) goto 800

    if (cons_state(eq_de) <= local_D_floor) then
      repair_event=NUM_C2P_LOW_D
      goto 800
    end if

    ! Momentos Covariantes Iniciales
    S_cov(1) = cons_state(eq_vx)
    S_cov(2) = cons_state(eq_vy)
    S_cov(3) = cons_state(eq_vz)

    ! Cuadrado del Momento: S^2 = S_i * g^ij * S_j
    S_sq = S_cov(1) * (ginv(1,1)*S_cov(1) + ginv(1,2)*S_cov(2) + ginv(1,3)*S_cov(3)) &
         + S_cov(2) * (ginv(2,1)*S_cov(1) + ginv(2,2)*S_cov(2) + ginv(2,3)*S_cov(3)) &
         + S_cov(3) * (ginv(3,1)*S_cov(1) + ginv(3,2)*S_cov(2) + ginv(3,3)*S_cov(3))

    ! --- ESCUDO DE CAUSALIDAD (LIMÍTE DEL CONO DE LUZ) ---
    E_tot = cons_state(eq_pr) + cons_state(eq_de)
    if (.not.ieee_is_finite(E_tot).or..not.ieee_is_finite(S_sq)) goto 800
    if (S_sq<0.0d0) goto 800
    ! E >= sqrt(D^2+S^2) is the admissibility bound, not a velocity cap.
    S_sq_max = max(0.0d0, (E_tot-cons_state(eq_de))*(E_tot+cons_state(eq_de)))
    if (.not.ieee_is_finite(S_sq_max)) goto 800

    if (S_sq > S_sq_max .and. S_sq > 1.0d-20) then
#ifndef USE_OPENACC
      call count_numerical_event(NUM_C2P_MOMENTUM)
#endif
      scale_factor = sqrt(v_max * S_sq_max / S_sq)
      cons_state(eq_vx) = cons_state(eq_vx) * scale_factor
      cons_state(eq_vy) = cons_state(eq_vy) * scale_factor
      cons_state(eq_vz) = cons_state(eq_vz) * scale_factor
    end if

    ! --- SOLUCIONADOR DE RAÍCES ---
    ! Asumimos que root_safe_newton recalcula S_sq internamente o usa cons_state actualizado
    W = root_safe_newton(cons_state, i, j, k)

    repair_event=NUM_C2P_ROOT_FAIL
    if (.not.ieee_is_finite(W)) goto 800
    if (W<1.0d0) goto 800

    call eval_V_and_dV(W, cons_state, V, dV)
    repair_event=NUM_C2P_BAD_V
    if (.not.ieee_is_finite(V).or..not.ieee_is_finite(dV)) goto 800
    if (V<=0.0d0) goto 800

    prim_state(eq_de) = cons_state(eq_de) / (W * sqg)
    ! Avoid subtracting two O(rho*W^2) quantities to recover small pressure.
    prim_state(eq_pr) = (cons_state(eq_pr)-cons_state(eq_de)*(W-1.0d0)) / ((g1*W**2-1.0d0)*sqg)

    ! EXTRACCIÓN DE VELOCIDAD CONTRAVARIANTE (Sana y libre de desincronización)
    ! v^i = (gamma^ij * S_j) / V. Utilizamos siempre el cons_state final.
    prim_state(eq_vx) = (ginv(1,1)*cons_state(eq_vx) + ginv(1,2)*cons_state(eq_vy) + ginv(1,3)*cons_state(eq_vz)) / V
    prim_state(eq_vy) = (ginv(2,1)*cons_state(eq_vx) + ginv(2,2)*cons_state(eq_vy) + ginv(2,3)*cons_state(eq_vz)) / V
    prim_state(eq_vz) = (ginv(3,1)*cons_state(eq_vx) + ginv(3,2)*cons_state(eq_vy) + ginv(3,3)*cons_state(eq_vz)) / V

    repair_event=NUM_C2P_SANITIZE
    if (.not.all(ieee_is_finite(prim_state))) goto 800
    if (prim_state(eq_de)<rho_floor.or.prim_state(eq_pr)<p_floor) goto 800

    v_sq_final = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_sq_final = v_sq_final + g(ii,jj) * prim_state(eq_vx + ii - 1) * prim_state(eq_vx + jj - 1)
      end do
    end do

    if (.not.ieee_is_finite(v_sq_final)) goto 800
    if (v_sq_final>=v_max) goto 800
#ifndef USE_OPENACC
    if (do_numerical_diagnostics) call record_numerical_change(diagnostic_before,cons_state)
#endif
    return

800 continue
    prim_state(eq_de) = rho_floor
    prim_state(eq_pr) = p_floor
    prim_state(eq_vx:eq_vz) = 0.0d0
    cons_state(eq_de) = local_D_floor
    cons_state(eq_pr) = local_tau_floor
    cons_state(eq_vx:eq_vz) = 0.0d0
#ifndef USE_OPENACC
    call count_numerical_event(repair_event)
    if (do_numerical_diagnostics) then
      ! A nonfinite input has no defined conservative delta. Count the repair
      ! and accumulate the change only for the finite input components.
      where (.not.finite_input) diagnostic_before=cons_state
      call record_numerical_change(diagnostic_before,cons_state)
    end if
#endif
  end subroutine cons_to_prim

  ! ===================================================================================
  ! HERRAMIENTAS MATEMÁTICAS PARA INVERSIÓN (Root-Finding)
  ! ===================================================================================

  ! Función auxiliar V(W) y su derivada dV/dW para Newton-Raphson.
  subroutine eval_V_and_dV(W, cons_state, V, dV)
    !$acc routine seq
    real*8, intent(in)  :: W, cons_state(:)
    real*8, intent(out) :: V, dV
    real*8 :: denom

    denom = g1 * W**2 - 1.0d0
    V  = (g1 * W**2 * (cons_state(eq_pr) + cons_state(eq_de)) - cons_state(eq_de) * W) / denom
    dV = (g1 * W * (cons_state(eq_de) * W - 2.0d0 * (cons_state(eq_pr) + cons_state(eq_de))) + cons_state(eq_de)) / denom**2
  end subroutine eval_V_and_dV

  ! Residual adimensional fn(W) y su derivada, a estado conservado fijo.
  subroutine eval_residual(W, cons_state, state_scale, scaled_momentum_sq, residual, d_residual)
    !$acc routine seq
    real*8, intent(in) :: W, cons_state(:), state_scale, scaled_momentum_sq
    real*8, intent(out) :: residual, d_residual
    real*8 :: V, dV

    call eval_V_and_dV(W, cons_state, V, dV)
    V=V/state_scale
    dV=dV/state_scale
    residual   = V**2 * (1.0d0 - 1.0d0/W**2) - scaled_momentum_sq
    d_residual = 2.0d0 * V * (dV * (W**3 - W) + V) / (W**3)
  end subroutine eval_residual

  ! Solucionador híbrido Newton-Raphson + Bisección
  function root_safe_newton(cons_state, i, j, k, x1_in, x2_in) result(root)
    !$acc routine seq
    implicit none
    real*8, intent(in) :: cons_state(:)
    integer, intent(in) :: i, j, k
    real*8, intent(in), optional :: x1_in, x2_in
    real*8 :: x1, x2
    real*8, parameter :: x_acc = 1.0d-12
    real*8 :: root
    integer, parameter :: MAX_ITER = 100
    integer :: iter_j, ii, jj
    real*8 :: df, diff_x, diff_x_old, f, f_high, f_low, temp, x_high, x_low
    logical :: converged

    real*8 :: momentum_sq, state_scale, scaled_momentum_sq, W_max_guess
    real*8 :: ginv(3,3)

    ! El estado conservado permanece fijo durante Newton y Brent.
    ginv(:,:) = gamma_inv_c(:,:,i,j)
    momentum_sq = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        momentum_sq = momentum_sq + ginv(ii,jj) * cons_state(eq_vx+ii-1) * cons_state(eq_vx+jj-1)
      end do
    end do
    state_scale = max(abs(cons_state(eq_pr)+cons_state(eq_de)),cons_state(eq_de),tiny(1.0d0))
    scaled_momentum_sq = momentum_sq / state_scale**2

    if (present(x1_in)) then ; x1 = x1_in ; else ; x1 = 1.0d0 ; end if

    if (present(x2_in)) then
      x2 = x2_in
    else
      W_max_guess = (sqrt(momentum_sq) + cons_state(eq_de)) / max(cons_state(eq_de), rho_floor)
      x2 = max(1.0d3, W_max_guess * 2.0d0)
    end if

    x_low = x1; x_high = x2
    converged = .false.

    call eval_residual(x1, cons_state, state_scale, scaled_momentum_sq, f_low, df)
    call eval_residual(x2, cons_state, state_scale, scaled_momentum_sq, f_high, df)

    if (abs(f_low) < 1.0d-14) then
      root = x_low ; return
    else if (abs(f_high) < 1.0d-14) then
      root = x_high ; return
    else if (sign(1.0d0, f_low) * sign(1.0d0, f_high) > 0.0d0) then
      root = -1.0d0 ; return
    end if

    if (f_low < 0.0d0) then
      x_low = x1 ; x_high = x2
    else
      x_high = x1 ; x_low = x2
    end if

    root = 0.5d0 * (x1 + x2)
    diff_x_old = abs(x2 - x1)
    diff_x = diff_x_old
    call eval_residual(root, cons_state, state_scale, scaled_momentum_sq, f, df)

    do iter_j = 1, MAX_ITER
      if (((root - x_high)*df - f)*((root - x_low)*df - f) > 0.0d0 .or. &
          abs(2.0d0*f) > abs(diff_x_old*df)) then
        diff_x_old = diff_x
        diff_x = 0.5d0 * (x_high - x_low)
        root = x_low + diff_x
        if (abs(x_low - root) < 1.0d-14) return
      else
        diff_x_old = diff_x
        diff_x = f / df
        temp = root
        root = root - diff_x
        if (abs(temp - root) < 1.0d-14) return
      end if

      if (abs(diff_x) < x_acc) then
        converged = .true.
        return
      end if

      call eval_residual(root, cons_state, state_scale, scaled_momentum_sq, f, df)

      if (f < 0.0d0) then
        x_low = root
      else
        x_high = root
      end if
    end do

    if (.not. converged .or. abs(f) > 1.0d-6) then
      root = root_brent(cons_state, momentum_sq, state_scale, scaled_momentum_sq)
    end if
  end function root_safe_newton

  ! Método de Brent
  function root_brent(cons_state, momentum_sq, state_scale, scaled_momentum_sq) result(root)
    !$acc routine seq
    implicit none
    real*8, intent(in) :: cons_state(:)
    real*8, intent(in) :: momentum_sq, state_scale, scaled_momentum_sq
    real*8 :: x2
    real*8, parameter :: x1 = 1.0d0, x_acc = 1.0d-12
    real*8 :: root
    integer, parameter :: MAX_ITER = 100
    real*8, parameter :: EPS = epsilon(1.0d0)
    integer :: iter
    real*8 :: a, b, c, d, e, fa, fb, fc, o, q, r, y1, tol1, xm, df

    real*8 :: W_max_guess

    W_max_guess = (sqrt(momentum_sq) + cons_state(eq_de)) / max(cons_state(eq_de), rho_floor)
    x2 = max(1.0d4, W_max_guess * 2.0d0) ! Límite dinámico blindado para Brent

    call eval_residual(x1, cons_state, state_scale, scaled_momentum_sq, fa, df)
    call eval_residual(x2, cons_state, state_scale, scaled_momentum_sq, fb, df)

    if (sign(1.0d0, fa) * sign(1.0d0, fb) > 0.0d0) then
      root = -1.0d0 ; return
    end if

    a = x1; b = x2; c = b; fc = fb
    d = 0.0d0; e = 0.0d0

    do iter = 1, MAX_ITER
      if (sign(1.0d0, fb) * sign(1.0d0, fc) > 0.0d0) then
        c = a; fc = fa; d = b - a; e = d
      end if
      if (abs(fc) < abs(fb)) then
        a = b; b = c; c = a; fa = fb; fb = fc; fc = fa
      end if

      tol1 = 2.0d0*EPS*abs(b) + 0.5d0*x_acc
      xm = 0.5d0*(c - b)
      if (abs(xm) <= tol1 .or. abs(fb) < 1.0d-14) then
        root = b ; return
      end if

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        y1 = fb/fa
        if (abs(a - c) < 1.0d-14) then
          o = 2.0d0*xm*y1; q = 1.0d0 - y1
        else
          q = fa/fc; r = fb/fc
          o = y1*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
          q = (q - 1.0d0)*(r - 1.0d0)*(y1 - 1.0d0)
        end if
        if (o > 0.0d0) q = -q
        o = abs(o)
        if (2.0d0*o < min(3.0d0*xm*q - abs(tol1*q), abs(e*q))) then
          e = d; d = o/q
        else
          d = xm; e = d
        end if
      else
        d = xm; e = d
      end if

      a = b; fa = fb
      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1, xm)
      end if
      call eval_residual(b, cons_state, state_scale, scaled_momentum_sq, fb, df)
    end do

    root = -1.0d0 ! Time-out irrecuperable
    return
  end function root_brent

  ! ===================================================================================
  ! FÍSICA EULERIANA (Flujos y Términos Fuente)
  ! ===================================================================================

  ! Subrutina: calc_fluxes
  ! Propósito: Calcula los Flujos Numéricos F(U) a través de las interfaces de celda.
  ! Cuantifica cuánta masa, momento y energía cruza dinámicamente la pared.

  subroutine calc_fluxes(prim_state, f_out, scan_dir, alpha, beta, g, sqg)
    !$acc routine seq
    implicit none
    real*8, intent(in)  :: prim_state(:)
    integer, intent(in) :: scan_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: f_out(:)

    real*8 :: v_sq, h, W, v_eff
    real*8 :: v_norm, beta_norm
    real*8 :: v(3), v_cov(3)
    integer :: ii, jj

    v(1) = prim_state(eq_vx)
    v(2) = prim_state(eq_vy)
    v(3) = prim_state(eq_vz)

    ! Cuadrado de la velocidad Euleriana tridimensional (Tensorial Universal)
    v_sq = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_sq = v_sq + g(ii,jj) * v(ii) * v(jj)
      end do
    end do

    ! --- ESCUDO DE CAUSALIDAD NUMÉRICA ---
    if (v_sq >= v_max) then
      v(1) = v(1) * sqrt(v_max / v_sq)
      v(2) = v(2) * sqrt(v_max / v_sq)
      v(3) = v(3) * sqrt(v_max / v_sq)
      v_sq = v_max
    end if

    ! Calcular el vector de velocidad covariante (v_i = gamma_ij * v^j)
    v_cov = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_cov(ii) = v_cov(ii) + g(ii,jj) * v(jj)
      end do
    end do

    ! Extrae la velocidad alineada con la dirección espacial del barrido actual
    if (scan_dir == DIR_X) then
      v_norm = v(1); beta_norm = beta(1)
    else if (scan_dir == DIR_Y) then
      v_norm = v(2); beta_norm = beta(2)
    else
      v_norm = v(3); beta_norm = beta(3)
    end if

    h = 1.0d0 + g1 * prim_state(eq_pr) / prim_state(eq_de)
    W = 1.0d0 / sqrt(1.0d0 - v_sq)

    ! --- VELOCIDAD DE ADVECCIÓN EFECTIVA ---
    v_eff = alpha * v_norm - beta_norm

    ! Ecuación de Continuidad: F(D)
    f_out(eq_de) = prim_state(eq_de) * W * v_eff

    ! Ecuación de Energía: F(tau)
    f_out(eq_pr) = prim_state(eq_de) * h * W**2 * v_eff + prim_state(eq_pr) * beta_norm - f_out(eq_de)

    ! Ecuaciones de Momento: F(S_i) (Usando el vector v_cov pre-calculado)
    f_out(eq_vx) = prim_state(eq_de) * h * W**2 * v_cov(1) * v_eff + merge(alpha*prim_state(eq_pr), 0.0d0, scan_dir == DIR_X)
    f_out(eq_vy) = prim_state(eq_de) * h * W**2 * v_cov(2) * v_eff + merge(alpha*prim_state(eq_pr), 0.0d0, scan_dir == DIR_Y)
    f_out(eq_vz) = prim_state(eq_de) * h * W**2 * v_cov(3) * v_eff + merge(alpha*prim_state(eq_pr), 0.0d0, scan_dir == DIR_Z)

    ! Transformamos de flujo propio a flujo coordenado usando sqrt_gamma (sqg)
    f_out = f_out * sqg
  end subroutine calc_fluxes


  ! Subrutina: calc_sources
  ! Propósito: Calcula los Términos Fuente S(U). Aportan los cambios en el momento y
  ! energía causados exclusivamente por la geometría del espacio curvo (Gravedad).

  subroutine calc_sources(prim_state, src_out, i, j, k)
    !$acc routine seq
    implicit none
    real*8, intent(in)  :: prim_state(:)
    integer, intent(in) :: i, j, k
    real*8, intent(out) :: src_out(:)

    real*8 :: T_mn(0:3,0:3), chris(0:3,0:3,0:3), dg(0:3,0:3,1:3)
    real*8 :: alpha, sqg, dlna(0:3)
    integer :: mu, nu, dir_idx

    src_out(:) = 0.0d0

    if (flat_cartesian_sources) return

    ! 1. LEER TODO DE LA CACHÉ GLOBAL (Gravedad a coste computacional CERO)
    alpha      = alpha_c(i,j)
    sqg        = sqrt_gamma_c(i,j)
    chris(:,:,:) = chris_c(:,:,:,i,j)
    dg(:,:,:)    = dg_c(:,:,:,i,j)
    dlna(:)      = dlna_c(:,i,j)

    ! 2. EVALUAR LA MATERIA ON-THE-FLY
    call calculate_stress_energy_tensor(prim_state, i, j, k, T_mn)

    ! --- FUERZA GRAVITACIONAL SOBRE EL MOMENTO ESPACIAL (S_j) ---
    do dir_idx = 1, 3
      do mu = 0, 3
        do nu = 0, 3
          src_out(eq_vx + dir_idx - 1) = src_out(eq_vx + dir_idx - 1) + 0.5d0 * T_mn(mu,nu) * dg(mu,nu,dir_idx)
        end do
      end do
    end do

    ! --- TRABAJO GRAVITACIONAL SOBRE LA ENERGÍA (tau) ---
    do mu = 0, 3
      src_out(eq_pr) = src_out(eq_pr) + T_mn(mu,0) * dlna(mu)
      do nu = 0, 3
        src_out(eq_pr) = src_out(eq_pr) - T_mn(mu,nu) * chris(0,mu,nu)
      end do
    end do

    src_out = src_out * alpha * sqg
    src_out(eq_pr) = src_out(eq_pr) * alpha
  end subroutine calc_sources

end module equations
