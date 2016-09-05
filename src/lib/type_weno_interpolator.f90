module type_weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO interpolator object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
use type_weno_alpha_coefficient
use type_weno_optimal_weights
use type_weno_IS
use type_weno_poly_coefficients
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
  private
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
    procedure(description_interface), pass(self), deferred, public :: description
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
  subroutine constructor_interface(self, constructor, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create a WENO interpolator.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_constructor, weno_interpolator, weno_optimal_weights, weno_IS, weno_poly_coefficients
  class(weno_interpolator),      intent(inout) :: self              !< WENO interpolator.
  class(weno_constructor),       intent(in)    :: constructor       !< WENO constructor.
  class(weno_optimal_weights),   intent(in)    :: optimal_weights   !< WENO optimal weights.
  class(weno_IS),                intent(in)    :: IS                !< WENO smoothness indicators coefficients.
  class(weno_poly_coefficients), intent(in)    :: poly_coefficients !< WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  pure subroutine description_interface(self, alfa, weno_IS, poly_coefficients, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator, weno_alpha_coefficient, weno_optimal_weights, weno_IS, weno_poly_coefficients
  class(weno_interpolator),      intent(in)  :: self              !< WENO interpolator.
  class(weno_alpha_coefficient), intent(in)  :: alpha             !< WENO weights.
  class(weno_IS),                intent(in)  :: IS                !< WENO smoothness indicators.
  class(weno_poly_coefficients), intent(in)  :: poly_coefficients !< WENO smoothness indicators.
  character(len=:), allocatable, intent(out) :: string            !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  pure subroutine interpolate_interface(self, alpha_coefficient, S, stencil, location, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stencil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: I_P, R_P, weno_interpolator, weno_alpha_coefficient, weno_optimal_weights, weno_IS, weno_poly_coefficients
  class(weno_interpolator),      intent(in)  :: self               !< WENO interpolator.
  class(weno_alpha_coefficient), intent(in)  :: alpha              !< WENO weights.
  class(weno_optimal_weights),   intent(in)  :: optimal_weights    !< WENO optimal weights.
  class(weno_IS),                intent(in)  :: IS                 !< WENO IS coefficients and values.
  class(weno_poly_coefficients), intent(in)  :: poly_coefficients  !< WENO polynomial coefficients.
  real(R_P),                     intent(in)  :: location           !< Location of the interpolation.
  real(R_P),                     intent(out) :: interpolation(1:)  !< Result of the interpolation, [1:2].
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interpolate_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_interpolator
