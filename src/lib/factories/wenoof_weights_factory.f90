!< Wenoof weights factory.
module wenoof_weights_factory
!< Wenoof weights factory.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only: I_P
use wenoof_alpha_object  , only : alpha_object_constructor
use wenoof_beta_object   , only : beta_object_constructor
use wenoof_kappa_object  , only : kappa_object_constructor
use wenoof_weights_object, only : weights_object, weights_object_constructor
use wenoof_weights_int_js, only : weights_int_js, weights_int_js_constructor
use wenoof_weights_rec_js, only : weights_rec_js, weights_rec_js_constructor

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
  type is(weights_int_js_constructor)
    allocate(weights_int_js :: object)
  type is(weights_rec_js_constructor)
    allocate(weights_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor(interpolator_type, S, alpha_constructor, beta_constructor, kappa_constructor, &
                                constructor)
  !< Create an instance of concrete extension of [[weights_object_constructor]].
  character(*),                                   intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                   intent(in)           :: S                 !< Stencils dimension.
  class(alpha_object_constructor),                intent(in)           :: alpha_constructor !< Alpha constructor.
  class(beta_object_constructor),                 intent(in)           :: beta_constructor  !< Beta constructor.
  class(kappa_object_constructor),                intent(in)           :: kappa_constructor !< kappa constructor.
  class(weights_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    allocate(weights_int_js_constructor :: constructor)
  case('interpolator-M-JS')
    allocate(weights_int_js_constructor :: constructor)
  case('interpolator-M-Z')
    allocate(weights_int_js_constructor :: constructor)
  case('interpolator-Z')
    allocate(weights_int_js_constructor :: constructor)
  case('reconstructor-JS')
    allocate(weights_rec_js_constructor :: constructor)
  case('reconstructor-M-JS')
    allocate(weights_rec_js_constructor :: constructor)
  case('reconstructor-M-Z')
    allocate(weights_rec_js_constructor :: constructor)
  case('reconstructor-Z')
    allocate(weights_rec_js_constructor :: constructor)
  case default
    write(stderr, '(A)') 'error: interpolator type "'//trim(adjustl(interpolator_type))//'" is unknown!'
    stop
  endselect
  call constructor%create(S=S)
  select type(constructor)
  type is(weights_int_js_constructor)
    allocate(constructor%alpha_constructor, source=alpha_constructor)
    allocate(constructor%beta_constructor, source=beta_constructor)
    allocate(constructor%kappa_constructor, source=kappa_constructor)
  type is(weights_rec_js_constructor)
    allocate(constructor%alpha_constructor, source=alpha_constructor)
    allocate(constructor%beta_constructor, source=beta_constructor)
    allocate(constructor%kappa_constructor, source=kappa_constructor)
  endselect
  endsubroutine create_constructor
endmodule wenoof_weights_factory
