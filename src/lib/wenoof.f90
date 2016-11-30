module wenoof
!-----------------------------------------------------------------------------------------------------------------------------------
!< WenOOF, WENO interpolation Object Oriented Fortran library
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : R_P, I_P
use wenoof_interpolator_abstract, only : wenoof_constructor, wenoof_interpolator
use wenoof_interpolator_js,       only : wenoof_constructor_upwind, wenoof_interpolator_upwind
use wenoof_alpha_coefficient_abstract
use wenoof_smoothness_indicators_js
use wenoof_alpha_coefficient_js
use wenoof_alpha_coefficient_z
use wenoof_alpha_coefficient_m
use wenoof_optimal_weights_js
use wenoof_polynomials_js
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: wenoof_factory, wenoof_constructor, wenoof_interpolator
public :: wenoof_constructor_upwind, wenoof_interpolator_upwind
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: wenoof_factory
  !< WENO factory aimed to create and return a concrete WENO interpolator to the client code without exposing the concrete
  !< interpolators classes.
  contains
    procedure, nopass :: create
endtype
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine create(constructor, IS_type, alpha_type, alpha_base_type, weights_opt_type, polynomial_type, interpolator)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create and return a concrete WENO interpolator object being an extension of the abstract *weno_interpolator* type.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(wenoof_constructor),  intent(IN)               :: constructor       !< The concrete WENO constructor selected by client code.
  character(*),               intent(IN)               :: IS_type           !< The concrete WENO smoothness indicator.
  character(*),               intent(IN)               :: alpha_type        !< The concrete WENO alpha coefficient.
  character(*),               intent(IN), optional     :: alpha_base_type   !< The WENO alpha coefficient base for WENO Mapped.
  character(*),               intent(IN)               :: weights_opt_type  !< The concrete WENO optimal weights.
  character(*),               intent(IN)               :: polynomial_type   !< The concrete WENO polynomial.
  class(wenoof_interpolator), allocatable, intent(OUT) :: interpolator      !< The concrete WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(constructor)
  type is(wenoof_constructor_upwind)
    allocate(wenoof_interpolator_upwind :: interpolator)
    ! Instantiate WENO smoothness indicators
    select case(IS_type)
    case('JS')
      allocate(IS_js :: interpolator%IS)
    case default
      call interpolator%init_error(error_code = 3_I_P)
    endselect
    ! Instantiate WENO alpha coefficients
    select case(alpha_type)
    case('JS')
      allocate(alpha_coefficient_js :: interpolator%alpha)
    case('Z')
      allocate(alpha_coefficient_z :: interpolator%alpha)
    case('M')
      allocate(alpha_coefficient_m :: interpolator%alpha)
      if (.not.present(alpha_base_type)) call interpolator%init_error(error_code = 0_I_P)
      associate(alpha => interpolator%alpha)
        select type(alpha)
        type is(alpha_coefficient_m)
          if (alpha_base_type.NE.'JS'.or.alpha_base_type.NE.'Z') call interpolator%init_error(error_code = 1_I_P)
          call alpha%initialize(alpha_base = alpha_base_type)
        endselect
      endassociate
    case default
      call interpolator%init_error(error_code = 2_I_P)
    endselect
    ! Instantiate WENO optimal weights
    select case(weights_opt_type)
    case('JS')
      allocate(optimal_weights_js :: interpolator%weights)
    case default
      call interpolator%init_error(error_code = 4_I_P)
    endselect
    ! Instantiate WENO polynomials
    select case(polynomial_type)
    case('JS')
      allocate(polynomials_js :: interpolator%polynom)
    case default
      call interpolator%init_error(error_code = 5_I_P)
    endselect
      call interpolator%create(constructor=constructor, IS_type=interpolator%IS, alpha_type=interpolator%alpha,&
                               weights_opt_type = interpolator%weights, polynomial_type = interpolator%polynom)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create
endmodule wenoof
