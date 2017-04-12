!< Wenoof interpolations factory.
module wenoof_interpolations_factory
!< Wenoof interpolations factory.

use penf, only : I_P, R_P
use wenoof_interpolations_object, only : interpolations_object, interpolations_object_constructor
use wenoof_interpolations_rec_js, only : interpolations_rec_js, interpolations_rec_js_constructor
use wenoof_interpolations_int_js, only : interpolations_int_js, interpolations_int_js_constructor

implicit none
private
public :: interpolations_factory

type :: interpolations_factory
  !< Factory, create an instance of concrete extension of [[interpolations_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create                                          !< Create a concrete instance of [[interpolations_object]].
    procedure, nopass :: create_constructor_rec
    procedure, nopass :: create_constructor_int
    generic           :: create_constructor => create_constructor_rec, & !< Create a concrete instance
                                               create_constructor_int    !< of [[interpolations_object_constructor]].
endtype interpolations_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[interpolations_object]] given its constructor.
  class(interpolations_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolations_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(interpolations_int_js_constructor)
    allocate(interpolations_int_js :: object)
  type is(interpolations_rec_js_constructor)
    allocate(interpolations_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor_rec(interpolator_type, S, constructor)
  !< Create an instance of concrete extension of [[beta_object_constructor]].
  character(*),                                          intent(in)  :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                          intent(in)  :: S                 !< Stencils dimension.
  class(interpolations_object_constructor), allocatable, intent(out) :: constructor       !< Constructor.

  allocate(interpolations_rec_js_constructor :: constructor)
  call constructor%create(S=S)
  endsubroutine create_constructor_rec

  subroutine create_constructor_int(interpolator_type, S, stencil, x_target, constructor)
  !< Create an instance of concrete extension of [[beta_object_constructor]].
  character(*),                                          intent(in)  :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                          intent(in)  :: S                 !< Stencils dimension.
  real(R_P),                                             intent(in)  :: stencil(1-S:)     !< Stencil used for inter, [1-S:-1+S].
  real(R_P),                                             intent(in)  :: x_target          !< Coordinate of the interp point.
  class(interpolations_object_constructor), allocatable, intent(out) :: constructor       !< Constructor.

  allocate(interpolations_int_js_constructor :: constructor)
  allocate(constructor%stencil(1-S:S-1))
  constructor%stencil = stencil
  constructor%x_target = x_target
  call constructor%create(S=S)
  endsubroutine create_constructor_int
endmodule wenoof_interpolations_factory
