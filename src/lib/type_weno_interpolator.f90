module type_weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO interpolator object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
use type_weno_alpha_coefficient
use type_weno_optimal_weights
use type_weno_smoothness_indicators
use type_weno_polynomials
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_interpolator, weno_constructor
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_constructor
  !< Abstract type used for create new concrete WENO interpolators.
  !<
  !< @note Every concrete WENO interpolator implementations must define their own constructor type.
  private
endtype weno_constructor

type, abstract :: weno_interpolator
  !< WENO interpolator object.
  !<
  !< @note Do not implement any real interpolator: provide the interface for the different interpolators implemented.
  class(weno_IS), pointer                :: IS      => null() !< Pointer to the WENO smoothness indicators.
  class(weno_alpha_coefficient), pointer :: alpha   => null() !< Pointer to the WENO alpha coefficients.
  class(weno_optimal_weights), pointer   :: weights => null() !< Pointer to the WENO optimal weights.
  class(weno_polynomials), pointer       :: polynom => null() !< Pointer to the WENO polynomilas.
  private
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
    procedure(description_interface), pass(self), deferred, public :: description
    procedure(init_error_interface),  pass(self), deferred, public :: init_error
    procedure(interpolate_interface), pass(self), deferred, public :: interpolate
endtype weno_interpolator

abstract interface
  elemental subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  class(weno_interpolator), intent(inout) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  subroutine constructor_interface(self, constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create a WENO interpolator.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_constructor, weno_interpolator
  class(weno_interpolator), intent(inout) :: self              !< WENO interpolator.
  class(weno_constructor),  intent(in)    :: constructor       !< WENO constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  class(weno_interpolator),      intent(in)  :: self              !< WENO interpolator.
  character(len=:), allocatable, intent(out) :: string            !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  pure subroutine init_error_interface(self, error_code, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing the error in the choose of smoothness indicators, alpha coefficients, optimal weights, polynomials
  !<  for the selected WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: I_P, weno_interpolator
  class(weno_interpolator),      intent(in)  :: self              !< WENO interpolator.
  integer(I_P),                  intent(in)  :: error_code        !< Error code.
  character(len=:), allocatable, intent(out) :: string            !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_error_interface
endinterface

abstract interface
  pure subroutine interpolate_interface(self, S, stencil, location, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stencil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: I_P, R_P, weno_interpolator, weno_alpha_coefficient, weno_optimal_weights, weno_IS, weno_polynomials
  class(weno_interpolator), intent(in)  :: self                 !< WENO interpolator.
  integer,                  intent(IN)  :: S                    !< Number of stencils actually used.
  real(R_P),                intent(IN)  :: stencil(1:, 1 - S:)  !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(R_P),                intent(in)  :: location             !< Location of the interpolation.
  real(R_P),                intent(out) :: interpolation(1:)    !< Result of the interpolation, [1:2].
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interpolate_interface
endinterface

!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_interpolator
