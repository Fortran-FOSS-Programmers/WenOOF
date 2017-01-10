!< Wenoof objects factory.
module wenoof_objects_factory
!< Wenoof factory.

use wenoof_alpha_coefficients_js
use wenoof_alpha_coefficients_m
use wenoof_alpha_coefficients_z
use wenoof_base_object
use wenoof_optimal_weights_js
use wenoof_polynomials_js
use wenoof_smoothness_indicators_js

implicit none
private
public :: objects_factory

type :: objects_factory
  !< Factory, create an instance of concrete extension of [base_object] given its constructor.
  contains
    procedure, nopass :: create !< Create an instance of concrete extension of [base_object].
endtype objects_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [base_object] given its constructor.
  class(base_object_constructor),  intent(in)  :: constructor !< Constructor
  class(base_object), allocatable, intent(out) :: object      !< Object

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
  endsubroutine create
endmodule wenoof_objects_factory
