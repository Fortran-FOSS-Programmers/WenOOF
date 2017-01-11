!< WenOOF, WENO interpolation Object Oriented Fortran library
module wenoof
!< WenOOF, WENO interpolation Object Oriented Fortran library

use penf
use wenoof_alpha_coefficients
use wenoof_alpha_coefficients_js
use wenoof_alpha_coefficients_z
use wenoof_alpha_coefficients_m
use wenoof_interpolator
use wenoof_interpolator_js
use wenoof_optimal_weights
use wenoof_optimal_weights_js
use wenoof_polynomials
use wenoof_polynomials_js
use wenoof_smoothness_indicators
use wenoof_smoothness_indicators_js

implicit none
private
public :: interpolator
public :: wenoof_create

contains
  subroutine wenoof_create(interpolator_type, S, wenoof_interpolator, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                     intent(in)           :: interpolator_type   !< Type of the interpolator.
  integer(I_P),                     intent(in)           :: S                   !< Stencil dimension.
  class(interpolator), allocatable, intent(out)          :: wenoof_interpolator !< The concrete WENO interpolator.
  real(R_P),                        intent(in), optional :: eps                 !< Parameter for avoiding division by zero.
  class(smoothness_indicators_constructor), allocatable  :: is_constructor      !< Smoothness indicators constructor.
  class(alpha_coefficients_constructor),    allocatable  :: alpha_constructor   !< Alpha coefficients constructor.
  class(optimal_weights_constructor),       allocatable  :: weights_constructor !< Optimal weights constructor.
  class(polynomials_constructor),           allocatable  :: polynom_constructor !< Polynomials constructor.
  class(interpolator_constructor),          allocatable  :: interp_constructor  !< Interpolator constructor.

  select case(trim(adjustl(interpolator_type)))
  case('JS')
    call create_smoothness_indicators_js_constructor(S=S, constructor=is_constructor)
    call create_alpha_coefficients_js_constructor(S=S, constructor=alpha_constructor)
    call create_optimal_weights_js_constructor(S=S, constructor=weights_constructor)
    call create_polynomials_js_constructor(S=S, constructor=polynom_constructor)
    call create_interpolator_js_constructor(is=is_constructor,           &
                                            alpha=alpha_constructor,     &
                                            weights=weights_constructor, &
                                            polynom=polynom_constructor, &
                                            S=S, eps=eps,                &
                                            constructor=interp_constructor)
    allocate(interpolator_js :: wenoof_interpolator)
    call wenoof_interpolator%create(constructor=interp_constructor)
  case('JS-Z')
    ! @TODO add Z support
  case('JS-M')
    ! @TODO add M support
  case default
    ! @TODO add error handling
  endselect
  endsubroutine wenoof_create
endmodule wenoof
