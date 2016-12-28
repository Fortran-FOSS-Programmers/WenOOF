module wenoof_interpolator_abstract
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO interpolator object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
use wenoof_alpha_coefficient_abstract
use wenoof_optimal_weights_abstract
use wenoof_smoothness_indicators_abstract
use wenoof_polynomials_abstract
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: wenoof_interpolator, wenoof_constructor
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: wenoof_constructor
  !< Abstract type used for create new concrete WENO interpolators.
  !<
  !< @note Every concrete WENO interpolator implementations must define their own constructor type.
  private
endtype wenoof_constructor

type, abstract :: wenoof_interpolator
  !< WENO interpolator object.
  !<
  !< @note Do not implement any real interpolator: provide the interface for the different interpolators implemented.
  class(IS),                allocatable :: IS      !< Pointer to the WENO smoothness indicators.
  class(alpha_coefficient), allocatable :: alpha   !< Pointer to the WENO alpha coefficients.
  class(optimal_weights),   allocatable :: weights !< Pointer to the WENO optimal weights.
  class(polynomials),       allocatable :: polynom !< Pointer to the WENO polynomilas.
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
    procedure(description_interface), pass(self), deferred, public :: description
    procedure(init_error_interface),  pass(self), deferred, public :: init_error
    procedure(interpolate_interface), pass(self), deferred, public :: interpolate
endtype wenoof_interpolator

abstract interface
  elemental subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: wenoof_interpolator
  class(wenoof_interpolator), intent(inout) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  subroutine constructor_interface(self, constructor, IS_type, alpha_type, weights_opt_type, polynomial_type)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create a WENO interpolator.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: wenoof_constructor, wenoof_interpolator, IS, alpha_coefficient, optimal_weights, polynomials
  class(wenoof_interpolator), intent(inout) :: self              !< WENO interpolator.
  class(wenoof_constructor),  intent(in)    :: constructor       !< WENO constructor.
  class(IS),                  intent(in)    :: IS_type           !< The concrete WENO smoothness indicator.
  class(alpha_coefficient),   intent(in)    :: alpha_type        !< The concrete WENO alpha coefficient.
  class(optimal_weights),     intent(in)    :: weights_opt_type  !< The concrete WENO optimal weights.
  class(polynomials),         intent(in)    :: polynomial_type   !< The concrete WENO polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: wenoof_interpolator
  class(wenoof_interpolator),    intent(in)  :: self              !< WENO interpolator.
  character(len=:), allocatable, intent(out) :: string            !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  subroutine init_error_interface(self, error_code)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing the error in the choose of smoothness indicators, alpha coefficients, optimal weights, polynomials
  !<  for the selected WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: I_P, wenoof_interpolator
  class(wenoof_interpolator), intent(inout) :: self              !< WENO interpolator.
  integer(I_P),               intent(in)    :: error_code        !< Error code.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_error_interface
endinterface

abstract interface
  pure subroutine interpolate_interface(self, S, stencil, location, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stencil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: I_P, R_P, wenoof_interpolator, alpha_coefficient, optimal_weights, IS, polynomials
  class(wenoof_interpolator), intent(inout) :: self                 !< WENO interpolator.
  integer,                    intent(IN)    :: S                    !< Number of stencils actually used.
  real(R_P),                  intent(IN)    :: stencil(1:, 1 - S:)  !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  character(*),               intent(in)    :: location             !< Location of interpolated values: left, right, both.
  real(R_P),                  intent(out)   :: interpolation(1:)    !< Result of the interpolation, [1:2].
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interpolate_interface
endinterface

!-----------------------------------------------------------------------------------------------------------------------------------
endmodule wenoof_interpolator_abstract
