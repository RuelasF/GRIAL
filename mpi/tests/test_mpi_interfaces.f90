program test_mpi_interfaces
  use parallel_runtime
  use domain_decomposition
  use halo_exchange
  use variables
  use metrics
  use initialization, only: allocate_and_grid,precalculate_metric_cache
  use conditions, only: set_boundary_conditions
  use reconstruction
  use fluxes
  use evolution, only: calc_rhs_cpu,calc_rhs_accelerator,initialize_rhs_accelerator,finalize_rhs_accelerator
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none
  integer, parameter :: config_count=14
  integer :: fixture,profile,config,solver,checks=0,fixtures=0,unit
  character(len=64) :: mode,report_name
  real*8 :: state_error=0.0d0,flux_error=0.0d0,rhs_error=0.0d0,shared_error=0.0d0
  real*8 :: face_error=0.0d0,face_flux_error=0.0d0,face_rhs_error=0.0d0

  call parallel_initialize()
  call get_command_argument(1,mode)
  if (trim(mode)=='bad-offset'.or.trim(mode)=='bad-total') call check_invalid_metadata()
  call check_policy_windows()
  do fixture=1,5
    do profile=1,5
      call check_fixture()
    end do
  end do
  if (trim(mode)=='report') then
    write(report_name,'(A,I0,A)') 'interface_rank_',parallel_rank,'.txt'
    open(newunit=unit,file=trim(report_name),status='new',action='write')
    write(unit,'(4(I0,1X),7(ES24.16E3,1X))') parallel_rank,parallel_size,checks,fixtures, &
      state_error,flux_error,rhs_error,shared_error,face_error,face_flux_error,face_rhs_error
    close(unit)
  end if
  write(*,'(A,I0,A,I0,A,I0)') 'INTERFACE_TEST rank=',parallel_rank,' size=',parallel_size,' checks=',checks
  if (parallel_rank==0) write(*,'(A)') 'MPI INTERFACE TESTS PASSED'
  call parallel_finalize()

contains

  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    checks=checks+1
    if (.not.ok) then
      write(*,*) 'INTERFACE fixture/profile/config/solver/rank:',fixture,profile,config,solver,parallel_rank
      call parallel_abort('interface test: '//label)
    end if
  end subroutine require

  subroutine choose_config(number)
    integer, intent(in) :: number
    integer, parameter :: methods(7)=[REC_GODUNOV,REC_TVD,REC_TVD,REC_TVD,REC_WENO3,REC_MP5,REC_WENO5]
    integer, parameter :: limiters(7)=[LIM_MC,LIM_MINMOD,LIM_SUPERBEE,LIM_MC,LIM_MC,LIM_MC,LIM_MC]
    integer :: choice
    choice=mod(number-1,7)+1
    rec_method_id=methods(choice)
    tvd_limiter_id=limiters(choice)
    use_shock_sensor=number>7
  end subroutine choose_config

  subroutine check_policy_windows()
    integer, parameter :: n=37,nlocal=9,origins(7)=[0,1,2,3,10,26,28]
    integer :: i,direction,geometry,variant,origin,selected,offset,position,method
    integer(int64) :: expected(7)
    real*8 :: full(neq,-3:n+3),part(neq,-3:nlocal+3),ql(neq,0:n),qr(neq,0:n)
    real*8 :: pl(neq,0:nlocal),pr(neq,0:nlocal),stencil(neq,-2:3),fl(neq),fr(neq),error,scale
    nghost=3;rho_floor=1.0d-10;p_floor=1.0d-13
    do i=-3,n+3
      full(eq_de,i)=1.0d0+dble(i*i*i)/1048576.0d0
      full(eq_pr,i)=0.5d0+dble(i*i)/65536.0d0
      full(eq_vx:eq_vz,i)=0.0d0
    end do
    do geometry=0,1
      metric_type='Minkowski'
      if (geometry==1) metric_type='Kerr-Schild'
      curved_metric=geometry==1
      do variant=0,1
        shock_sensor_thresholds=SENSOR_DEFAULTS
        if (variant==1) shock_sensor_thresholds=[0.8d0,0.4d0,1.5d0,0.9d0,10.0d0,10.0d0,100.0d0]
        do config=1,config_count
          call choose_config(config)
          do direction=DIR_X,DIR_Z
            do_numerical_diagnostics=.false.
            call reconstruct_1d_core(full,n,ql,qr,direction)
            do position=1,size(origins)
              origin=origins(position)
              part=full(:,origin-3:origin+nlocal+3)
              numerical_events=0_int64
              do_numerical_diagnostics=.true.
              call reconstruct_1d_core(part,nlocal,pl,pr,direction,origin,n)
              error=max(maxval(abs(pl-ql(:,origin:origin+nlocal))),maxval(abs(pr-qr(:,origin:origin+nlocal))))
              state_error=max(state_error,error)
              call require(error==0.0d0,'global-position reconstruction window')
              expected=0_int64
              do i=0,nlocal
                method=rec_method_id
                selected=origin+i
                if (curved_metric) then
                  if (direction==DIR_X.and.selected<=0.and.(use_shock_sensor.or.rec_method_id>=REC_WENO3)) then
                    method=REC_GODUNOV
                  else if (use_shock_sensor.and.rec_method_id==REC_MP5) then
                    if ((direction==DIR_X.and.selected<=3).or. &
                        (direction==DIR_Y.and.(selected<=2.or.selected>=n-2))) method=REC_WENO3
                  end if
                end if
                expected(method)=expected(method)+1_int64
                do offset=-2,3
                  stencil(:,offset)=part(:,i+offset)
                end do
                call reconstruct_face_state(stencil,selected,n,direction,rec_method_id,tvd_limiter_id, &
                  use_shock_sensor,curved_metric,rho_floor,p_floor,fl,fr,shock_sensor_thresholds)
                scale=max(1.0d0,maxval(abs(pl(:,i))),maxval(abs(pr(:,i))))
                error=max(maxval(abs(fl-pl(:,i))),maxval(abs(fr-pr(:,i))))/scale
                face_error=max(face_error,error)
                call require(error<=5.0d-15,'face/strip policy window')
              end do
              call require(all(numerical_events(1:7)==expected),'global method audit counters')
            end do
          end do
        end do
      end do
    end do
    do_numerical_diagnostics=.false.;numerical_events=0_int64;numerical_delta_u=0.0d0
    shock_sensor_thresholds=SENSOR_DEFAULTS
  end subroutine check_policy_windows

  subroutine configure(shape_global,periodic)
    integer, intent(out) :: shape_global(3)
    logical, intent(out) :: periodic
    shape_global=[37,1,1]
    geom_type='Cartesian';metric_type='Minkowski';case_name='Sod'
    r_min=0.0d0;r_max=1.0d0;y_min=0.0d0;y_max=pi;z_min=-0.5d0;z_max=0.5d0
    use_log_r=.false.;a_spin=0.0d0;bh_mass=1.0d0
    periodic=.false.
    select case(fixture)
    case(2)
      shape_global(3)=7;case_name='KHI';periodic=.true.;r_min=-0.5d0;r_max=0.5d0
    case(3)
      shape_global(3)=7;case_name='Jet';geom_type='Cylindrical';z_min=0.0d0;z_max=1.0d0
    case(4)
      shape_global(2)=5;case_name='FontDaigneSag';geom_type='Spheroidal';metric_type='Kerr-Schild'
      r_min=1.2d0;r_max=12.0d0;use_log_r=.true.
    case(5)
      shape_global=[37,4,6];case_name='FishMonc3D';geom_type='Spheroidal';metric_type='Kerr-Schild'
      r_min=1.2d0;r_max=12.0d0;use_log_r=.true.;a_spin=0.9d0;z_min=0.0d0;z_max=2.0d0*pi
    end select
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3);nghost=3
    adb_idx=4.0d0/3.0d0;g1=adb_idx/(adb_idx-1.0d0)
    rho_floor=1.0d-10;p_floor=1.0d-13
    do_mdot_extraction=.false.;do_radial_mass_budget=.false.;do_numerical_diagnostics=.false.
    shock_sensor_thresholds=SENSOR_DEFAULTS
    call set_metric_type()
    call init_wavespeed_solver()
  end subroutine configure

  subroutine seed_profile()
    integer :: i,j,k,global_i
    real*8 :: factor
    p=-12345.0d0
    do k=1,nz
      do j=1,ny
        do i=1,nx
          global_i=x_global_offset+i
          p(i,j,k,eq_de)=1.0d0+dble(global_i*global_i)/65536.0d0+dble(j+k)/4096.0d0
          p(i,j,k,eq_pr)=0.5d0+dble(global_i*global_i*global_i)/1048576.0d0
          p(i,j,k,eq_vx)=0.001d0+dble(global_i)/1048576.0d0
          p(i,j,k,eq_vy)=dble(j)/4194304.0d0
          p(i,j,k,eq_vz)=0.001d0+dble(k)/1048576.0d0
          factor=1.0d0
          if (global_i>=20) then
            if (profile==2) factor=1.25d0
            if (profile==3) factor=1.75d0
            if (profile==4) factor=32.0d0
          end if
          p(i,j,k,eq_de:eq_pr)=factor*p(i,j,k,eq_de:eq_pr)
          if (profile==5.and.global_i>=18.and.global_i<=22) then
            p(i,j,k,eq_de)=10.0d0*rho_floor
            p(i,j,k,eq_pr)=10.0d0*p_floor
          end if
        end do
      end do
    end do
  end subroutine seed_profile

  subroutine capture_x_packet(packet)
    real*8, intent(out) :: packet(neq,4,0:nx,ny,nz)
    real*8 :: line(neq,-nghost:nx+nghost),ql(neq,0:nx),qr(neq,0:nx)
    real*8 :: stencil(neq,-2:3),fl(neq),fr(neq),computed(neq),error,scale
    integer :: i,j,k,offset,rs
    do k=1,nz
      do j=1,ny
        do i=-nghost,nx+nghost
          line(:,i)=p(i,j,k,:)
        end do
        call reconstruct_1d_core(line,nx,ql,qr,DIR_X,x_global_offset,nx_global)
        packet(:,1,:,j,k)=ql;packet(:,2,:,j,k)=qr
        do i=0,nx
          do offset=-2,3
            stencil(:,offset)=line(:,i+offset)
          end do
          call reconstruct_face_state(stencil,x_global_offset+i,nx_global,DIR_X,rec_method_id,tvd_limiter_id, &
            use_shock_sensor,curved_metric,rho_floor,p_floor,fl,fr,shock_sensor_thresholds)
          scale=max(1.0d0,maxval(abs(ql(:,i))),maxval(abs(qr(:,i))))
          error=max(maxval(abs(fl-ql(:,i))),maxval(abs(fr-qr(:,i))))/scale
          face_error=max(face_error,error)
          call require(error<=5.0d-15,'radial face/strip states')
          do rs=RS_HLLE,RS_HLLC
            riemann_solver_id=rs
            if (trim(geom_type)=='Cylindrical'.and.abs(x_face(i))<=64.0d0*epsilon(1.0d0)) then
              packet(:,2+rs,i,j,k)=0.0d0;computed=0.0d0
            else
              call resolve_riemann_problem(ql(:,i),qr(:,i),DIR_X,alpha_f_x(i+1,j),beta_f_x(:,i+1,j), &
                gamma_f_x(:,:,i+1,j),sqrt_gamma_f_x(i+1,j),packet(:,2+rs,i,j,k))
              call resolve_riemann_problem(fl,fr,DIR_X,alpha_f_x(i+1,j),beta_f_x(:,i+1,j), &
                gamma_f_x(:,:,i+1,j),sqrt_gamma_f_x(i+1,j),computed)
            end if
            scale=max(1.0d0,maxval(abs(packet(:,2+rs,i,j,k))))
            error=maxval(abs(computed-packet(:,2+rs,i,j,k)))/scale
            face_flux_error=max(face_flux_error,error)
            call require(error<=2.0d-13,'radial face/strip flux')
          end do
        end do
      end do
    end do
  end subroutine capture_x_packet

  subroutine check_shared_face(packet,layout)
    real*8, intent(in) :: packet(neq,4,0:nx,ny,nz)
    type(domain_layout), intent(in) :: layout
#ifdef USE_MPI
    integer :: lower,upper,status
    type(MPI_Status) :: receipt
    real*8 :: sent(neq,4,ny,nz),received(neq,4,ny,nz),error
    lower=MPI_PROC_NULL;upper=MPI_PROC_NULL
    if (layout%lower_rank>=0) lower=layout%lower_rank
    if (layout%upper_rank>=0) upper=layout%upper_rank
    sent=packet(:,:,nx,:,:)
    call MPI_Sendrecv(sent,size(sent),MPI_DOUBLE_PRECISION,upper,401,received,size(received), &
      MPI_DOUBLE_PRECISION,lower,401,MPI_COMM_WORLD,receipt,status)
    call parallel_check_mpi(status,'test shared interface')
    if (lower/=MPI_PROC_NULL) then
      error=maxval(abs(received-packet(:,:,0,:,:)))
      shared_error=max(shared_error,error)
      call require(error==0.0d0,'two owners of the same face agree')
    end if
#endif
  end subroutine check_shared_face

  subroutine check_fixture()
    real*8, allocatable :: accelerator_rhs(:,:,:,:)
    integer :: shape_global(3)
    logical :: periodic
    type(domain_layout) :: layout
    character(len=160) :: message
    real*8, allocatable :: reference_x(:,:,:,:,:,:),reference_rhs(:,:,:,:,:,:),packet(:,:,:,:,:)
    real*8 :: error,scale
    call configure(shape_global,periodic)
    call allocate_and_grid();call precalculate_metric_cache()
    call seed_profile();call set_boundary_conditions(p)
    allocate(reference_x(neq,4,0:nx,ny,nz,config_count),reference_rhs(nx,ny,nz,neq,2,config_count))
    do config=1,config_count
      call choose_config(config)
      call capture_x_packet(reference_x(:,:,:,:,:,config))
      do solver=RS_HLLE,RS_HLLC
        riemann_solver_id=solver
        call calc_rhs_cpu(p,reference_rhs(:,:,:,:,solver,config))
      end do
    end do
    call release_grid()

    call build_domain(shape_global,3,parallel_rank,parallel_size,periodic,layout,message)
    call parallel_require_collective(len_trim(message)==0,message)
    nx=shape_global(1);ny=shape_global(2);nz=shape_global(3)
    call allocate_and_grid(layout);call precalculate_metric_cache()
    call seed_profile();call initialize_halo_exchange(layout);call fill_primitive_halos(p)
    allocate(accelerator_rhs,mold=rhs)
    call initialize_rhs_accelerator(.true.)
    allocate(packet(neq,4,0:nx,ny,nz))
    do config=1,config_count
      call choose_config(config)
      call capture_x_packet(packet)
      error=maxval(abs(packet(:,1:2,:,:,:)-reference_x(:,1:2,x_global_offset:x_global_offset+nx,:,:,config)))
      state_error=max(state_error,error)
      call require(error==0.0d0,'local/full radial reconstructed states')
      error=maxval(abs(packet(:,3:4,:,:,:)-reference_x(:,3:4,x_global_offset:x_global_offset+nx,:,:,config)))
      flux_error=max(flux_error,error)
      call require(error==0.0d0,'local/full radial numerical fluxes')
      call check_shared_face(packet,layout)
      do solver=RS_HLLE,RS_HLLC
        riemann_solver_id=solver
        call calc_rhs_cpu(p,rhs)
        error=maxval(abs(rhs-reference_rhs(x_global_offset+1:x_global_offset+nx,:,:,:,solver,config)))
        rhs_error=max(rhs_error,error)
        call require(error==0.0d0,'local/full static RHS')
        call calc_rhs_accelerator(p,accelerator_rhs)
        scale=max(1.0d0,maxval(abs(rhs)))
        error=maxval(abs(rhs-accelerator_rhs))/scale
        face_rhs_error=max(face_rhs_error,error)
        call require(error<=2.0d-13,'host face/strip static RHS')
      end do
    end do
    call finalize_rhs_accelerator(.true.)
    call finalize_halo_exchange();call release_grid()
    fixtures=fixtures+1
  end subroutine check_fixture

  subroutine check_invalid_metadata()
    real*8 :: line(neq,-3:11),ql(neq,0:8),qr(neq,0:8)
    integer :: offset,total
    nghost=3;line=1.0d0;rec_method_id=REC_GODUNOV
    offset=0;total=8
    if (parallel_rank==parallel_size-1) then
      if (trim(mode)=='bad-offset') offset=-1
      if (trim(mode)=='bad-total') total=7
    end if
    call reconstruct_1d_core(line,8,ql,qr,DIR_X,offset,total)
    call parallel_require_collective(.true.,'wait for metadata rejection')
    call parallel_abort('invalid reconstruction metadata not rejected')
  end subroutine check_invalid_metadata

  subroutine release_grid()
    deallocate(x,x_face,y,y_face,z,z_face,u,up,rhs,p)
    deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (allocated(alpha_f_y)) deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    if (allocated(alpha_f_z)) deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
  end subroutine release_grid
end program test_mpi_interfaces
