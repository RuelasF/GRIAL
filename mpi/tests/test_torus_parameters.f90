program test_torus_parameters
  use variables
  use parameters
  implicit none
  integer :: unit_file
  character(len=*), parameter :: path='/tmp/grial_torus_parameters_test.par'

  fd_inner_radius=6.25d0;fd_pressure_max_radius=9.25d0;fd_polytropic_constant=.0015d0
  fm_inner_radius=8.0d0;fm_pressure_max_radius=12.0d0;fm_polytropic_constant=.002d0
  call write_case('fishbone_equatorial','fm_',6.5d0)
  call read_parameter_header(path)
  call read_parameter_file(path,.false.)
  if(case_id/=10.or.fd_inner_radius/=6.5d0.or.fm_inner_radius/=8.0d0) &
    error stop 'Legacy fishbone/FM keys must still select FD without changing FM'

  call write_case('fd_equatorial','fd_',7.0d0)
  call read_parameter_header(path)
  call read_parameter_file(path,.false.)
  if(case_id/=10.or.fd_inner_radius/=7.0d0.or.fm_inner_radius/=8.0d0) &
    error stop 'Canonical FD keys/model mismatch'

  call write_case('fm_equatorial','fm_',8.5d0)
  call read_parameter_header(path)
  call read_parameter_file(path,.false.)
  if(case_id/=13.or.fm_inner_radius/=8.5d0.or.fd_inner_radius/=7.0d0) &
    error stop 'Canonical FM must have independent parameters'

  open(newunit=unit_file,file=path,status='old')
  close(unit_file,status='delete')
  print *, 'FD/FM PARAMETER TESTS PASSED'
contains
  subroutine write_case(problem_name,prefix,inner)
    character(len=*), intent(in) :: problem_name,prefix
    real*8, intent(in) :: inner
    open(newunit=unit_file,file=path,status='replace')
    write(unit_file,'(A)') 'problem='//problem_name
    write(unit_file,'(A)') 'restart=false'
    write(unit_file,'(A,ES24.16E3)') prefix//'inner_radius=',inner
    write(unit_file,'(A)') prefix//'pressure_max_radius=13'
    write(unit_file,'(A)') prefix//'polytropic_constant=0.002'
    close(unit_file)
  end subroutine write_case
end program test_torus_parameters
