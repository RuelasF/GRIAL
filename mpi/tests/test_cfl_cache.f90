! Checks invalidation on one rank, changed parameters, geometry and clipping.
program test_cfl_cache
  use variables
  use evolution, only: update_adaptive_dt
  use parallel_runtime
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  character(len=32) :: mode
  real*8 :: expected
  call parallel_initialize()
  call get_command_argument(1,mode)
  nx=4;ny=1;nz=1;nx_global=4*parallel_size;ny_global=1;nz_global=1
  x_global_offset=4*parallel_rank;y_global_offset=0
  dx=0.5d0;dy=1.0d0;dz=1.0d0;CFL=0.4d0
  allocate(alpha_c(nx,ny,nz),beta_c(3,nx,ny,nz),gamma_inv_c(3,3,nx,ny,nz))
  alpha_c=1.0d0;beta_c=0.0d0;gamma_inv_c=0.0d0
  gamma_inv_c(1,1,:,:,:)=1.0d0
  metric_cache_generation=metric_cache_generation+1_int64
  expected=0.2d0
  call check()
  dt=1.0d-10 ! An event-clipped step must not poison the cached causal bound.
  call check()
  if (parallel_rank==parallel_size-1) then
    alpha_c=2.0d0
    metric_cache_generation=metric_cache_generation+1_int64
  end if
  expected=0.1d0
  call check()
  call check()
  CFL=0.2d0;expected=0.05d0;call check()
  dx=0.25d0;expected=0.025d0;call check()
  if (trim(mode)=='bad_metric') then
    if (parallel_rank==parallel_size-1) then
      alpha_c=-1.0d0
      metric_cache_generation=metric_cache_generation+1_int64
    end if
    call update_adaptive_dt()
    error stop 'invalid metric was accepted'
  else if (trim(mode)=='bad_cfl') then
    if (parallel_rank==parallel_size-1) CFL=ieee_value(1.0d0,ieee_quiet_nan)
    call update_adaptive_dt()
    error stop 'invalid CFL was accepted'
  end if
  deallocate(alpha_c,beta_c,gamma_inv_c)
  nx=5
  allocate(alpha_c(nx,ny,nz),beta_c(3,nx,ny,nz),gamma_inv_c(3,3,nx,ny,nz))
  alpha_c=1.0d0;beta_c=0.0d0;gamma_inv_c=0.0d0
  gamma_inv_c(1,1,:,:,:)=1.0d0
  expected=0.05d0;call check()
  print *, 'CFL CACHE CHECKS PASSED rank=',parallel_rank
  call parallel_finalize()
contains
  subroutine check()
    call update_adaptive_dt()
    if (abs(dt-expected)>8*epsilon(dt)) error stop 'wrong cached time step'
  end subroutine check
end program test_cfl_cache
