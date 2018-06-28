!< Wenoof interpolator factory.
module wenoof_interpolator_factory
!< Wenoof interpolator factory.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P
use wenoof_interpolations_object, only : interpolations_object_constructor
use wenoof_interpolator_object, only : interpolator_object, interpolator_object_constructor
use wenoof_interpolator_js, only : interpolator_js, interpolator_js_constructor
use wenoof_reconstructor_js, only : reconstructor_js, reconstructor_js_constructor
use wenoof_weights_object, only : weights_object_constructor

implicit none
private
public :: interpolator_factory

type :: interpolator_factory
  !< Factory, create an instance of concrete extension of [[interpolator_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create             !< Create a concrete instance of [[interpolator_object]].
    procedure, nopass :: create_constructor !< Create a concrete instance of [[interpolator_object_constructor]].
endtype interpolator_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[interpolator_object]] given its constructor.
  class(interpolator_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolator_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(interpolator_js_constructor)
    allocate(interpolator_js :: object)
  type is(reconstructor_js_constructor)
    allocate(reconstructor_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor(interpolator_type, S, interpolations_constructor, weights_constructor, &
                                constructor)
  !< Create an instance of concrete extension of [[interpolator_object_constructor]].
  character(*),                                        intent(in)  :: interpolator_type          !< Type of interpolator.
  integer(I_P),                                        intent(in)  :: S                          !< Stencils dimension.
  class(interpolations_object_constructor),            intent(in)  :: interpolations_constructor !< Interpolations const.
  class(weights_object_constructor),                   intent(in)  :: weights_constructor        !< Weights constructor.
  class(interpolator_object_constructor), allocatable, intent(out) :: constructor                !< Constructor.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS', 'interpolator-M-JS', 'interpolator-M-Z', 'interpolator-Z')
    allocate(interpolator_js_constructor :: constructor)
  case('reconstructor-JS', 'reconstructor-M-JS', 'reconstructor-M-Z', 'reconstructor-Z')
    allocate(reconstructor_js_constructor :: constructor)
  case default
    write(stderr, '(A)') 'error: interpolator type "'//trim(adjustl(interpolator_type))//'" is unknown!'
    stop
  endselect
  call constructor%create(S=S)
  select type(constructor)
  type is(interpolator_js_constructor)
    allocate(constructor%interpolations_constructor, mold=interpolations_constructor)
    constructor%interpolations_constructor = interpolations_constructor
    allocate(constructor%weights_constructor, mold=weights_constructor)
    constructor%weights_constructor = weights_constructor
  type is(reconstructor_js_constructor)
    allocate(constructor%interpolations_constructor, mold=interpolations_constructor)
    constructor%interpolations_constructor = interpolations_constructor
    allocate(constructor%weights_constructor, mold=weights_constructor)
    constructor%weights_constructor = weights_constructor
  endselect
  endsubroutine create_constructor
endmodule wenoof_interpolator_factory
