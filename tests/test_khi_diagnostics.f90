program test_khi_diagnostics
  use variables
  use output, only: calc_khi_diagnostics
  implicit none
  integer :: i,k,m,unit_file
  real*8 :: row(9),expected(neq),distance
  character(len=*), parameter :: folder='/tmp/grhd_khi_diagnostics_test'
  nx=40; ny=1; nz=40
  dx=1.0d0/nx; dy=1.0d0; dz=1.0d0/nz
  case_name='KHI'; metric_type='Minkowski'; scheme_name='test'
  output_folder=folder
  g1=2.5d0
  allocate(x(nx),z(nz),p(nx,ny,nz,neq),up(nx,ny,nz,neq))
  p=0.0d0
  do i=1,nx
    x(i)=-0.5d0+(i-0.5d0)*dx
  end do
  do k=1,nz
    z(k)=-0.5d0+(k-0.5d0)*dz
    distance=abs(abs(z(k))-khi_half_width)
    do i=1,nx
      p(i,1,k,eq_de)=1.0d0
      p(i,1,k,eq_pr)=2.5d0
      p(i,1,k,eq_vz)=0.01d0*sin(10.0d0*pi*x(i))*exp(-0.5d0*(distance/0.05d0)**2)
      do m=1,neq
        up(i,1,k,m)=dble(m)+0.01d0*dble(i-k)
      end do
    end do
  end do
  expected=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0]
  call execute_command_line('mkdir -p '//folder//'/test')
  call calc_khi_diagnostics(1.234d0)
  open(newunit=unit_file,file=folder//'/test/khi_diagnostics.dat',status='old')
  read(unit_file,*) row
  close(unit_file,status='delete')
  if (abs(row(1)-1.234d0) > 1.0d-14 .or. maxval(abs(row(2:6)-expected)) > 1.0d-12) &
    error stop 'KHI conserved integrals/order mismatch'
  if (maxval(abs(row(7:8)-0.01d0/sqrt(2.0d0))) > 1.0d-7 .or. row(9) <= 0.0d0) &
    error stop 'KHI seeded-mode normalization/proxy mismatch'
  call execute_command_line('rmdir '//folder//'/test '//folder)
  print *, 'KHI DIAGNOSTIC TESTS PASSED'
end program test_khi_diagnostics
