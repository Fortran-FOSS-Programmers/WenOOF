!< Wenoof weights factory.
module wenoof_weights_factory
!< Wenoof weights factory.

use penf, only: I_P
use wenoof_alpha_object
use wenoof_beta_object
use wenoof_kappa_object
use wenoof_weights_object
use wenoof_weights_js

implicit none
private
public :: weights_factory

type :: weights_factory
  !< Factory, create an instance of concrete extension of [[weights_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create             !< Create a concrete instance of [[weights_object]].
    procedure, nopass :: create_constructor !< Create a concrete instance of [[weights_object_constructor]].
endtype weights_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[weights_object]] given its constructor.
  class(weights_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(weights_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(weights_js_constructor)
    allocate(weights_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor(interpolator_type, S, alpha_constructor, beta_constructor, kappa_constructor, &
                                constructor, face_left, face_right)
  !< Create an instance of concrete extension of [[weights_object_constructor]].
  character(*),                                   intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                   intent(in)           :: S                 !< Stencils dimension.
  class(alpha_object_constructor),                intent(in)           :: alpha_constructor !< Alpha constructor.
  class(beta_object_constructor),                 intent(in)           :: beta_constructor  !< Beta constructor.
  class(kappa_object_constructor),                intent(in)           :: kappa_constructor !< kappa constructor.
  class(weights_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  logical,                                        intent(in), optional :: face_left         !< Activate left-face interpolations.
  logical,                                        intent(in), optional :: face_right        !< Activate right-face interpolations.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    ! @TODO implement this
    error stop 'interpolator-JS to be implemented'
  case('reconstructor-JS')
    allocate(weights_js_constructor :: constructor)
  case('reconstructor-M-JS')
    allocate(weights_js_constructor :: constructor)
  case('reconstructor-M-Z')
    allocate(weights_js_constructor :: constructor)
  case('reconstructor-Z')
    allocate(weights_js_constructor :: constructor)
  endselect
  call constructor%create(S=S, face_left=face_left, face_right=face_right)
  select type(constructor)
  type is(weights_js_constructor)
    allocate(constructor%alpha_constructor, source=alpha_constructor)
    allocate(constructor%beta_constructor, source=beta_constructor)
    allocate(constructor%kappa_constructor, source=kappa_constructor)
  endselect
  endsubroutine create_constructor
endmodule wenoof_weights_factory
