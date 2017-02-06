!< Wenoof beta factory.
module wenoof_beta_factory
!< Wenoof beta factory.

use penf, only: I_P
use wenoof_beta_object
use wenoof_beta_rec_js
use wenoof_beta_int_js

implicit none
private
public :: beta_factory

type :: beta_factory
  !< Factory, create an instance of concrete extension of [[beta_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create             !< Create a concrete instance of [[beta_object]].
    procedure, nopass :: create_constructor !< Create a concrete instance of [[beta_object_constructor]].
endtype beta_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[beta_object]] given its constructor.
  class(beta_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(beta_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(beta_rec_js_constructor)
    allocate(beta_rec_js :: object)
  type is(beta_int_js_constructor)
    allocate(beta_int_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor(interpolator_type, S, constructor, face_left, face_right)
  !< Create an instance of concrete extension of [[beta_object_constructor]].
  character(*),                                intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                intent(in)           :: S                 !< Stencils dimension.
  class(beta_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  logical,                                     intent(in), optional :: face_left         !< Activate left-face interpolations.
  logical,                                     intent(in), optional :: face_right        !< Activate right-face interpolations.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    allocate(beta_int_js_constructor :: constructor)
  case('interpolator-M-JS')
    allocate(beta_int_js_constructor :: constructor)
  case('interpolator-M-Z')
    allocate(beta_int_js_constructor :: constructor)
  case('interpolator-Z')
    allocate(beta_int_js_constructor :: constructor)
  case('reconstructor-JS')
    allocate(beta_rec_js_constructor :: constructor)
  case('reconstructor-M-JS')
    allocate(beta_rec_js_constructor :: constructor)
  case('reconstructor-M-Z')
    allocate(beta_rec_js_constructor :: constructor)
  case('reconstructor-Z')
    allocate(beta_rec_js_constructor :: constructor)
  endselect
  call constructor%create(S=S, face_left=face_left, face_right=face_right)
  endsubroutine create_constructor
endmodule wenoof_beta_factory
