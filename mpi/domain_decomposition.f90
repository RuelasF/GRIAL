! CPU x/r by y/theta blocks, with the third direction complete.
module domain_decomposition
  implicit none
  private
  public :: domain_layout, build_domain, validate_domain, domain_owner

  type :: domain_layout
    integer :: global_shape(3)=0, local_shape(3)=0, offset(3)=0
    integer :: rank=0, ranks=1, halo=0
    integer :: process_shape(2)=[1,1], process_coords(2)=0
    integer :: lower_neighbor(2)=-1, upper_neighbor(2)=-1
    integer :: lower_rank=-1, upper_rank=-1 ! Legacy aliases for the first direction.
    logical :: periodic_x=.false., periodic_y=.false., active_axes(3)=.false.
    logical :: owns_lower(3)=.false., owns_upper(3)=.false.
    logical :: physical_lower(3)=.false., physical_upper(3)=.false.
  end type domain_layout

contains

  pure subroutine build_domain(global_shape,halo,rank,ranks,periodic_x,layout,message,process_shape,periodic_y)
    integer, intent(in) :: global_shape(3),halo,rank,ranks
    logical, intent(in) :: periodic_x
    type(domain_layout), intent(out) :: layout
    character(len=*), intent(out) :: message
    integer, intent(in), optional :: process_shape(2)
    logical, intent(in), optional :: periodic_y
    integer :: grid(2),coords(2),base,remainder,d,stride,lower,upper
    logical :: periodic(2)

    layout=domain_layout()
    message=''
    if (any(global_shape<1).or.halo<1) then
      message='global cell counts and halo must be positive'
      return
    end if
    if (any(global_shape>huge(1)-halo)) then
      message='grid plus halo exceeds integer index range'
      return
    end if
    if (ranks<1.or.rank<0.or.rank>=ranks) then
      message='invalid rank or process count'
      return
    end if
    grid=[ranks,1]
    if (present(process_shape)) grid=process_shape
    if (any(grid<1).or.any(grid>ranks)) then
      message='invalid process grid'
      return
    end if
    if (mod(ranks,grid(1))/=0.or.ranks/grid(1)/=grid(2)) then
      message='process grid must cover the communicator'
      return
    end if
    coords=[mod(rank,grid(1)),rank/grid(1)]
    periodic=[periodic_x,.false.]
    if (present(periodic_y)) periodic(2)=periodic_y
    layout%global_shape=global_shape
    layout%local_shape=global_shape
    layout%owns_lower=.true.;layout%owns_upper=.true.
    layout%physical_lower=.true.;layout%physical_upper=.true.
    do d=1,2
      base=global_shape(d)/grid(d)
      remainder=mod(global_shape(d),grid(d))
      ! Local dimension switches remain valid: never collapse an active split axis.
      if (grid(d)>1.and.base<max(2,halo)) then
        message='first-axis block too small for halo or active dimension'
        if (d==2) message='polar block too small for halo or active dimension'
        return
      end if
      layout%local_shape(d)=base
      if (coords(d)<remainder) layout%local_shape(d)=base+1
      layout%offset(d)=coords(d)*base+min(coords(d),remainder)
      layout%owns_lower(d)=coords(d)==0
      layout%owns_upper(d)=coords(d)==grid(d)-1
      layout%physical_lower(d)=layout%owns_lower(d).and..not.periodic(d)
      layout%physical_upper(d)=layout%owns_upper(d).and..not.periodic(d)
      lower=coords(d)-1;upper=coords(d)+1
      if (periodic(d)) then
        lower=modulo(lower,grid(d));upper=modulo(upper,grid(d))
      end if
      stride=1
      if (d==2) stride=grid(1)
      if (lower>=0) layout%lower_neighbor(d)=rank+(lower-coords(d))*stride
      if (upper<grid(d)) layout%upper_neighbor(d)=rank+(upper-coords(d))*stride
    end do
    layout%rank=rank;layout%ranks=ranks;layout%halo=halo
    layout%process_shape=grid;layout%process_coords=coords
    layout%periodic_x=periodic(1);layout%periodic_y=periodic(2)
    layout%active_axes=global_shape>1
    layout%lower_rank=layout%lower_neighbor(1);layout%upper_rank=layout%upper_neighbor(1)
  end subroutine build_domain

  pure subroutine validate_domain(layout,message)
    type(domain_layout), intent(in) :: layout
    character(len=*), intent(out) :: message
    type(domain_layout) :: expected
    call build_domain(layout%global_shape,layout%halo,layout%rank,layout%ranks, &
      layout%periodic_x,expected,message,layout%process_shape,layout%periodic_y)
    if (len_trim(message)>0) return
    if (any(layout%local_shape/=expected%local_shape).or.any(layout%offset/=expected%offset).or. &
        layout%lower_rank/=expected%lower_rank.or.layout%upper_rank/=expected%upper_rank.or. &
        any(layout%process_coords/=expected%process_coords).or. &
        any(layout%lower_neighbor/=expected%lower_neighbor).or.any(layout%upper_neighbor/=expected%upper_neighbor).or. &
        any(layout%active_axes.neqv.expected%active_axes).or. &
        any(layout%owns_lower.neqv.expected%owns_lower).or.any(layout%owns_upper.neqv.expected%owns_upper).or. &
        any(layout%physical_lower.neqv.expected%physical_lower).or. &
        any(layout%physical_upper.neqv.expected%physical_upper)) message='inconsistent domain layout'
  end subroutine validate_domain

  pure integer function axis_owner(cells,parts,index) result(owner)
    integer, intent(in) :: cells,parts,index
    integer :: base,remainder,wide
    owner=-1
    if (index<1.or.index>cells.or.parts<1) return
    base=cells/parts;remainder=mod(cells,parts)
    if (base<1) return
    wide=base*remainder+remainder
    if (index<=wide) then
      owner=(index-1)/(base+1)
    else
      owner=remainder+(index-wide-1)/base
    end if
  end function axis_owner

  pure integer function domain_owner(layout,global_i,global_j) result(owner)
    type(domain_layout), intent(in) :: layout
    integer, intent(in) :: global_i
    integer, intent(in), optional :: global_j
    integer :: radial,polar
    owner=-1
    radial=axis_owner(layout%global_shape(1),layout%process_shape(1),global_i)
    ! Without j, the radial probe belongs to the caller's polar column.
    polar=layout%process_coords(2)
    if (present(global_j)) polar=axis_owner(layout%global_shape(2),layout%process_shape(2),global_j)
    if (radial<0.or.polar<0) return
    owner=radial+layout%process_shape(1)*polar
  end function domain_owner
end module domain_decomposition
