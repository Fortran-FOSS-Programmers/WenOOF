module wenoof
!-----------------------------------------------------------------------------------------------------------------------------------
!< WenOOF, WENO interpolation Object Oriented Fortran library
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use type_weno_interpolator, only : weno_constructor, weno_interpolator
use type_weno_interpolator_js, only : weno_constructor_js, weno_interpolator_js
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_factory, weno_constructor, weno_interpolator
public :: weno_constructor_upwind, weno_interpolator_upwind
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type :: weno_factory
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
  class(weno_constructor),  intent(IN)               :: constructor       !< The concrete WENO constructor selected by client code.
  character(*),             intent(IN)               :: IS_type           !< The concrete WENO smoothness indicator.
  character(*),             intent(IN)               :: alpha_type        !< The concrete WENO alpha coefficient.
  character(*),             intent(IN), optional     :: alpha_base_type   !< The WENO alpha coefficient base for WENO Mapped.
  character(*),             intent(IN)               :: weights_opt_type  !< The concrete WENO optimal weights.
  character(*),             intent(IN)               :: polynomial_type   !< The concrete WENO polynomial.
  class(weno_interpolator), allocatable, intent(OUT) :: interpolator      !< The concrete WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(constructor)
  type is(weno_constructor_upwind)
    allocate(weno_interpolator_upwind :: interpolator)
    ! Instantiate WENO smoothness indicators
    select case(IS_type)
    case('JS')
      allocate(weno_IS_js :: interpolator%IS)
    case default
      call interpolator%init_error(error_code = 3_I_P)
    endselect
    ! Instantiate WENO alpha coefficients
    select case(alpha_type)
    case('JS')
      allocate(weno_alpha_coefficient_js :: interpolator%alpha)
    case('Z')
      allocate(weno_alpha_coefficient_z :: interpolator%alpha)
    case('M')
      allocate(weno_alpha_coefficient_m :: interpolator%alpha)
      if (.not.present(alpha_base_type)) call interpolator%init_error(error_code = 0_I_P)
      select case(alpha_base_type)
      case('JS')
        allocate(weno_alpha_coefficient_js :: interpolator%alpha%alpha_base)
      case('Z')
        allocate(weno_alpha_coefficient_z :: interpolator%alpha%alpha_base)
      case default
        call interpolator%init_error(error_code = 1_I_P)
      endselect
    case default
      call interpolator%init_error(error_code = 2_I_P)
    endselect
    ! Instantiate WENO optimal weights
    select case(weights_opt_type)
    case('JS')
      allocate(weno_optimal_weights_js :: interpolator%weights)
    case default
      call interpolator%init_error(error_code = 4_I_P)
    endselect
    ! Instantiate WENO polynomials
    select case(polynomial_type)
    case('JS')
      allocate(weno_polynomials_js :: interpolator%polynom)
    case default
      call interpolator%init_error(error_code = 5_I_P)
    endselect
    call interpolator%create(constructor=constructor)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create
endmodule wenoof
