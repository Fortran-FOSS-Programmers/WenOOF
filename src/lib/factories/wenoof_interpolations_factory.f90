!< Wenoof interpolations factory.
module wenoof_interpolations_factory
!< Wenoof interpolations factory.

use penf, only: I_P
use wenoof_interpolations_object
use wenoof_interpolations_rec_js

implicit none
private
public :: interpolations_factory

type :: interpolations_factory
  !< Factory, create an instance of concrete extension of [[interpolations_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create             !< Create a concrete instance of [[interpolations_object]].
    procedure, nopass :: create_constructor !< Create a concrete instance of [[interpolations_object_constructor]].
endtype interpolations_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[interpolations_object]] given its constructor.
  class(interpolations_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolations_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(interpolations_rec_js_constructor)
    allocate(interpolations_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor(interpolator_type, S, constructor, face_left, face_right)
  !< Create an instance of concrete extension of [[beta_object_constructor]].
  character(*),                                          intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                          intent(in)           :: S                 !< Stencils dimension.
  class(interpolations_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  logical,                                               intent(in), optional :: face_left         !< Activate left-face interp.
  logical,                                               intent(in), optional :: face_right        !< Activate right-face interp.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    ! @TODO implement this
    error stop 'interpolator-JS to be implemented'
  case('reconstructor-JS')
    allocate(interpolations_rec_js_constructor :: constructor)
  case('reconstructor-M-JS')
    allocate(interpolations_rec_js_constructor :: constructor)
  case('reconstructor-M-Z')
    allocate(interpolations_rec_js_constructor :: constructor)
  case('reconstructor-Z')
    allocate(interpolations_rec_js_constructor :: constructor)
  endselect
  call constructor%create(S=S, face_left=face_left, face_right=face_right)
  endsubroutine create_constructor
endmodule wenoof_interpolations_factory
