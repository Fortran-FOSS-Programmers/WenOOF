!< WenOOF, WENO interpolation Object Oriented Fortran library
module wenoof
!< WenOOF, WENO interpolation Object Oriented Fortran library

use penf
use wenoof_interpolator_abstract
use wenoof_interpolator_js
use wenoof_alpha_coefficient_abstract
use wenoof_smoothness_indicators_abstract
use wenoof_smoothness_indicators_js
use wenoof_alpha_coefficient_js
use wenoof_alpha_coefficient_z
use wenoof_alpha_coefficient_m
use wenoof_optimal_weights_abstract
use wenoof_optimal_weights_js
use wenoof_polynomials_abstract
use wenoof_polynomials_js

implicit none
private
public :: wenoof_factory
public :: wenoof_constructor
public :: wenoof_interpolator
public :: wenoof_constructor_upwind
public :: wenoof_interpolator_upwind

type :: wenoof_factory
  !< WENO factory, create and return a concrete WENO interpolator.
  contains
    procedure, nopass :: create !< Create and return a concrete WENO interpolator object.
endtype
contains
  subroutine create(constructor, IS_type, alpha_type, weights_opt_type, polynomial_type, interpolator, alpha_base_type)
  !< Create and return a concrete WENO interpolator object.
  class(wenoof_constructor),               intent(in)           :: constructor       !< The concrete WENO constructor selected by
                                                                                     !< client code.
  character(*),                            intent(in)           :: IS_type           !< The concrete WENO smoothness indicator.
  character(*),                            intent(in)           :: alpha_type        !< The concrete WENO alpha coefficient.
  character(*),                            intent(in)           :: weights_opt_type  !< The concrete WENO optimal weights.
  character(*),                            intent(in)           :: polynomial_type   !< The concrete WENO polynomial.
  class(wenoof_interpolator), allocatable, intent(out)          :: interpolator      !< The concrete WENO interpolator.
  character(*),                            intent(in), optional :: alpha_base_type   !< The WENO alpha coefficient base for WENO
                                                                                     !< Mapped.
  class(IS),                allocatable                         :: IS_type_          !< The concrete WENO smoothness indicator.
  class(alpha_coefficient), allocatable                         :: alpha_type_       !< The concrete WENO alpha coefficient.
  class(optimal_weights),   allocatable                         :: weights_opt_type_ !< The concrete WENO optimal weights.
  class(polynomials),       allocatable                         :: polynomial_type_  !< The concrete WENO polynomial.

  select type(constructor)
  type is(wenoof_constructor_upwind)
    allocate(wenoof_interpolator_upwind :: interpolator)
    ! instantiate WENO smoothness indicators
    select case(IS_type)
    case('JS')
      allocate(IS_js :: IS_type_)
    case default
      call interpolator%init_error(error_code = 3_I_P)
    endselect
    ! instantiate WENO alpha coefficients
    select case(alpha_type)
    case('JS')
      allocate(alpha_coefficient_js :: alpha_type_)
    case('Z')
      allocate(alpha_coefficient_z :: alpha_type_)
    case('M')
      allocate(alpha_coefficient_m :: alpha_type_)
      if (.not.present(alpha_base_type)) call interpolator%init_error(error_code = 0_I_P)
      associate(alpha => alpha_type_)
        select type(alpha)
        type is(alpha_coefficient_m)
          if (alpha_base_type/='JS'.or.alpha_base_type/='Z') call interpolator%init_error(error_code = 1_I_P)
          call alpha%initialize(alpha_base = alpha_base_type)
        endselect
      endassociate
    case default
      call interpolator%init_error(error_code = 2_I_P)
    endselect
    ! instantiate WENO optimal weights
    select case(weights_opt_type)
    case('JS')
      allocate(optimal_weights_js :: weights_opt_type_)
    case default
      call interpolator%init_error(error_code = 4_I_P)
    endselect
    ! instantiate WENO polynomials
    select case(polynomial_type)
    case('JS')
      allocate(polynomials_js :: polynomial_type_)
    case default
      call interpolator%init_error(error_code = 5_I_P)
    endselect
      call interpolator%create(constructor=constructor,            &
                               IS_type=IS_type_,                   &
                               alpha_type=alpha_type_,             &
                               weights_opt_type=weights_opt_type_, &
                               polynomial_type=polynomial_type_)
  endselect
  endsubroutine create
endmodule wenoof
