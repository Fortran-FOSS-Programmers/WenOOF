!< Wenoof interpolator factory.
module wenoof_interpolator_factory
!< Wenoof interpolator factory.

use penf, only: I_P
use wenoof_interpolations_object
use wenoof_interpolator_object
! use wenoof_interpolator_js
use wenoof_reconstructor_js
use wenoof_weights_object

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
  ! type is(interpolator_js_constructor)
    ! allocate(interpolator_js :: object)
  type is(reconstructor_js_constructor)
    allocate(reconstructor_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor(interpolator_type, S, interpolations_constructor, weights_constructor, &
                                constructor, face_left, face_right)
  !< Create an instance of concrete extension of [[weights_object_constructor]].
  character(*),                                        intent(in)           :: interpolator_type          !< Type of interpolator.
  integer(I_P),                                        intent(in)           :: S                          !< Stencils dimension.
  class(interpolations_object_constructor),            intent(in)           :: interpolations_constructor !< Interpolations const.
  class(weights_object_constructor),                   intent(in)           :: weights_constructor        !< Weights constructor.
  class(interpolator_object_constructor), allocatable, intent(out)          :: constructor                !< Constructor.
  logical,                                             intent(in), optional :: face_left                  !< Activate left interp.
  logical,                                             intent(in), optional :: face_right                 !< Activate right interp.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    ! @TODO implement this
    error stop 'interpolator-JS to be implemented'
  case('reconstructor-JS')
    allocate(reconstructor_js_constructor :: constructor)
  endselect
  call constructor%create(S=S, face_left=face_left, face_right=face_right)
  select type(constructor)
  type is(reconstructor_js_constructor)
    allocate(constructor%interpolations_constructor, source=interpolations_constructor)
    allocate(constructor%weights_constructor, source=weights_constructor)
  endselect
  endsubroutine create_constructor
endmodule wenoof_interpolator_factory
