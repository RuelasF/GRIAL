program test_michel_initial
  use variables
  use metrics
  use conditions, only: setup_michel_accretion_initial, prepare_case_boundary_data, set_boundary_conditions
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: failures, coordinate, background, i
  real*8 :: radii(8), expected_rho(8), expected_ur(8), state(8,1,1,neq)
  real*8 :: alpha, beta(3), metric(3,3), velocity(3), w, ur, ut, h, energy, cs2
  real*8 :: reference_energy, physical_ur, normalization, q, tolerance

  failures = 0
  tolerance = 3.0d-10
  ! Fixtures from an independent SciPy Brent solve of the two Bernoulli branches.
  ! docs/grhd_1d/michel_reference.py; Aguayo-Ortiz et al. (2021), eqs. 2.14--2.16.
  radii = [1.5d0,2.0d0,3.0d0,7.9d0,8.0d0,8.1d0,40.0d0,100.0d0]
  expected_rho = [7.5947756314630588d-1,5.1861252366202615d-1,3.0861489319038488d-1, &
    1.0130314469189540d-1,1.0d-1,9.8734016261155549d-2,3.0327634167769221d-2,2.2948921848588506d-2]
  expected_ur = [-9.3631615417995240d-1,-7.7128874014750060d-1,-5.7605054616760343d-1, &
    -2.5307128656947275d-1,-2.5d-1,-2.4699214483118206d-1,-3.2973228128119300d-2,-6.9720050926854737d-3]
  adb_idx = 4.0d0/3.0d0
  g1 = adb_idx/(adb_idx-1.0d0)
  bh_mass = 1.0d0
  a_spin = 0.0d0
  nx = 8
  ny = 1
  nz = 1
  nghost = 3
  allocate(x(-nghost:nx+nghost),y(-nghost:ny+nghost),z(-nghost:nz+nghost))
  y = pi/2.0d0
  z = 0.0d0
  reference_energy = 1.1718041645257966d0
  do background = 1, 2
    metric_type = 'Eddington-Finkelstein'
    geom_type = 'Spherical'
    if (background == 2) then
      metric_type = 'Kerr-Schild'
      geom_type = 'Spheroidal'
    end if
    do coordinate = 1, 2
      use_log_r = coordinate == 2
      call set_metric_type()
      x(1:nx) = radii
      if (use_log_r) x(1:nx) = log(radii)
      call setup_michel_accretion_initial(state,8.0d0,0.1d0)
      call require('finite state across sonic point and horizon',all(ieee_is_finite(state)))
      do i = 1, nx
        call calculate_metric(x(i),y(1),alpha=alpha,beta=beta,gamma=metric)
        velocity = state(i,1,1,eq_vx:eq_vz)
        w = 1.0d0/sqrt(1.0d0-dot_product(velocity,matmul(metric,velocity)))
        ur = w*(velocity(1)-beta(1)/alpha)
        if (use_log_r) ur = ur*radii(i)
        call require('density independent fixture',abs(state(i,1,1,eq_de)/expected_rho(i)-1.0d0)<tolerance)
        call require('radial four-velocity independent fixture',abs(ur/expected_ur(i)-1.0d0)<tolerance)
        call require('constant mass flux',abs(-radii(i)**2*state(i,1,1,eq_de)*ur/1.6d0-1.0d0)<tolerance)
        h = 1.0d0+g1*state(i,1,1,eq_pr)/state(i,1,1,eq_de)
        energy = h*sqrt(1.0d0-2.0d0/radii(i)+ur**2)
        call require('constant Bernoulli energy',abs(energy/reference_energy-1.0d0)<tolerance)
        ut = w/alpha
        q = 2.0d0/radii(i)
        normalization = -(1.0d0-q)*ut**2+2.0d0*q*ut*ur+(1.0d0+q)*ur**2
        call require('future directed normalized four-velocity',ut>0.0d0.and.abs(normalization+1.0d0)<tolerance)
        cs2 = adb_idx*state(i,1,1,eq_pr)/(state(i,1,1,eq_de)*h)
        physical_ur = ur**2/(1.0d0-2.0d0/radii(i)+ur**2)
        if (radii(i)>2.0d0.and.radii(i)<8.0d0) call require('inner supersonic branch',physical_ur>cs2)
        if (radii(i)>8.0d0) call require('outer subsonic branch',physical_ur<cs2)
        if (radii(i)==8.0d0) call require('sonic regularity',abs(physical_ur-cs2)<tolerance)
      end do
    end do
  end do

  ! The outer injector must also select the physical branch, including restart preparation.
  metric_type = 'Eddington-Finkelstein'
  geom_type = 'Spherical'
  use_log_r = .false.
  call set_metric_type()
  x(1:nx) = radii
  x(nx+1:nx+nghost) = [101.0d0,102.0d0,103.0d0]
  case_name = 'MichelA'
  michel_critical_radius = 8.0d0
  michel_critical_density = 0.1d0
  allocate(p(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq))
  p = 0.0d0
  call prepare_case_boundary_data()
  state = p(1:nx,1:1,1:1,:)
  call set_boundary_conditions(p)
  call require('outer injector physical density',p(nx+1,1,1,eq_de)>0.02d0)
  call require('inner outflow ghost',all(p(0,1,1,:)==p(1,1,1,:)))
  call require('boundary preparation preserves reference',maxval(abs(state(:,1,1,eq_de)/expected_rho-1.0d0))<tolerance)
  if (failures /= 0) then
    write(*,*) 'MICHEL INITIAL TESTS FAILED: ',failures
    error stop 1
  end if
  write(*,*) 'MICHEL INITIAL TESTS PASSED'
contains
  subroutine require(label,passed)
    character(len=*), intent(in) :: label
    logical, intent(in) :: passed
    if (.not.passed) then
      failures = failures+1
      write(*,*) 'FAIL: ',label
    end if
  end subroutine require
end program test_michel_initial
