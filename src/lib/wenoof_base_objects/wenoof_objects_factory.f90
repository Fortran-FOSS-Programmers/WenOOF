!< Wenoof objects factory.
module wenoof_objects_factory
!< Wenoof factory.

use wenoof_alpha_coefficients
use wenoof_alpha_coefficients_js
use wenoof_alpha_coefficients_m
use wenoof_alpha_coefficients_z
use wenoof_base_object
use wenoof_optimal_weights
use wenoof_optimal_weights_js
use wenoof_polynomials
use wenoof_polynomials_js
use wenoof_smoothness_indicators
use wenoof_smoothness_indicators_js

implicit none
private
public :: objects_factory

type :: objects_factory
  !< Factory, create an instance of concrete extension of [[base_object]] given its constructor.
  contains
    ! public methods
    generic :: create => create_alpha_coefficients, &
                         create_optimal_weights,    &
                         create_polynomials,        &
                         create_smoothness_indicators !< Create a concrete instance of [[alpha_coefficients]] or
                                                      !< [[optimal_weights]] or [[polynomials]] or [[smoothness_indicators]].
    procedure, nopass :: create_base_object           !< Create a concrete instance of [[base_object]].
    ! private methods
    procedure, nopass, private :: create_alpha_coefficients    !< Create a concrete instance of [[alpha_coefficients]].
    procedure, nopass, private :: create_optimal_weights       !< Create a concrete instance of [[optimal_weights]].
    procedure, nopass, private :: create_polynomials           !< Create a concrete instance of [[polynomials]].
    procedure, nopass, private :: create_smoothness_indicators !< Create a concrete instance of [[smoothness_indicators]].
endtype objects_factory

contains
  subroutine create_base_object(constructor, object)
  !< Create an instance of concrete extension of [[base_object]] given its constructor.
  class(base_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(base_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(alpha_coefficients_js_constructor)
    allocate(alpha_coefficients_js :: object)
  type is(alpha_coefficients_m_constructor)
    allocate(alpha_coefficients_m :: object)
  type is(alpha_coefficients_z_constructor)
    allocate(alpha_coefficients_z :: object)
  type is(optimal_weights_js_constructor)
    allocate(optimal_weights_js :: object)
  type is(polynomials_js_constructor)
    allocate(polynomials_js :: object)
  type is(smoothness_indicators_js_constructor)
    allocate(smoothness_indicators_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_base_object

  subroutine create_alpha_coefficients(constructor, object)
  !< Create an instance of concrete extension of [[alpha_coefficients]] given its constructor.
  class(alpha_coefficients_constructor),  intent(in)  :: constructor !< Constructor.
  class(alpha_coefficients), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(alpha_coefficients_js_constructor)
    allocate(alpha_coefficients_js :: object)
  type is(alpha_coefficients_m_constructor)
    allocate(alpha_coefficients_m :: object)
  type is(alpha_coefficients_z_constructor)
    allocate(alpha_coefficients_z :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_alpha_coefficients

  subroutine create_optimal_weights(constructor, object)
  !< Create an instance of concrete extension of [[optimal_weights]] given its constructor.
  class(optimal_weights_constructor),  intent(in)  :: constructor !< Constructor.
  class(optimal_weights), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(optimal_weights_js_constructor)
    allocate(optimal_weights_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_optimal_weights

  subroutine create_polynomials(constructor, object)
  !< Create an instance of concrete extension of [[polynomials]] given its constructor.
  class(polynomials_constructor),  intent(in)  :: constructor !< Constructor.
  class(polynomials), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(polynomials_js_constructor)
    allocate(polynomials_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_polynomials

  subroutine create_smoothness_indicators(constructor, object)
  !< Create an instance of concrete extension of [[smoothness_indicators]] given its constructor.
  class(smoothness_indicators_constructor),  intent(in)  :: constructor !< Constructor.
  class(smoothness_indicators), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(smoothness_indicators_js_constructor)
    allocate(smoothness_indicators_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_smoothness_indicators
endmodule wenoof_objects_factory
