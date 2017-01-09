!< Abstract WENO interpolator object.
module wenoof_interpolator_abstract
!< Abstract WENO interpolator object.

use penf, only : I_P, R_P
use wenoof_alpha_coefficient_abstract
use wenoof_optimal_weights_abstract
use wenoof_smoothness_indicators_abstract
use wenoof_polynomials_abstract

implicit none
private
public :: wenoof_constructor
public :: wenoof_interpolator

type, abstract :: wenoof_constructor
  !< Abstract constructor of new concrete WENO interpolators.
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
    procedure(constructor_interface), pass(self), deferred :: create      !< Create interpolator.
    procedure(description_interface), pass(self), deferred :: description !< Return string-description of interpolator.
    procedure(destructor_interface),  pass(self), deferred :: destroy     !< Destroy interpolator.
    procedure(init_error_interface),  pass(self), deferred :: init_error  !< Return error code for wrong intialization.
    procedure(interpolate_interface), pass(self), deferred :: interpolate !< Interpolate values.
endtype wenoof_interpolator

abstract interface
  !< Create interpolator.
  subroutine constructor_interface(self, constructor, IS_type, alpha_type, weights_opt_type, polynomial_type)
  !< Create interpolator.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  import :: wenoof_constructor, wenoof_interpolator, IS, alpha_coefficient, optimal_weights, polynomials
  class(wenoof_interpolator), intent(inout) :: self             !< WENO interpolator.
  class(wenoof_constructor),  intent(in)    :: constructor      !< WENO constructor.
  class(IS),                  intent(in)    :: IS_type          !< Concrete WENO smoothness indicator.
  class(alpha_coefficient),   intent(in)    :: alpha_type       !< Concrete WENO alpha coefficient.
  class(optimal_weights),     intent(in)    :: weights_opt_type !< Concrete WENO optimal weights.
  class(polynomials),         intent(in)    :: polynomial_type  !< Concrete WENO polynomial.
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Destoy interpolator.
  elemental subroutine destructor_interface(self)
  !< Destoy interpolator.
  import :: wenoof_interpolator
  class(wenoof_interpolator), intent(inout) :: self !< WENO interpolator.
  endsubroutine destructor_interface
endinterface

abstract interface
  !< Return string-description of interpolator.
  pure subroutine description_interface(self, string)
  !< Return string-description of interpolator.
  import :: wenoof_interpolator
  class(wenoof_interpolator),    intent(in)  :: self   !< WENO interpolator.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  endsubroutine description_interface
endinterface

abstract interface
  !< Return error code for wrong intialization.
  subroutine init_error_interface(self, error_code)
  !< Return error code for wrong intialization.
  import :: I_P, wenoof_interpolator
  class(wenoof_interpolator), intent(inout) :: self       !< WENO interpolator.
  integer(I_P),               intent(in)    :: error_code !< Error code.
  endsubroutine init_error_interface
endinterface

abstract interface
  !< Interpolate values.
  pure subroutine interpolate_interface(self, S, stencil, location, interpolation)
  !< Interpolate values.
  import :: I_P, R_P, wenoof_interpolator, alpha_coefficient, optimal_weights, IS, polynomials
  class(wenoof_interpolator), intent(inout) :: self                !< WENO interpolator.
  integer,                    intent(in)    :: S                   !< Number of stencils actually used.
  real(R_P),                  intent(in)    :: stencil(1:, 1 - S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  character(*),               intent(in)    :: location            !< Location of interpolated values: left, right, both.
  real(R_P),                  intent(out)   :: interpolation(1:)   !< Result of the interpolation, [1:2].
  endsubroutine interpolate_interface
endinterface
endmodule wenoof_interpolator_abstract
