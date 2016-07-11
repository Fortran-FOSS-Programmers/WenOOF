module type_weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO interpolator object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
use type_weno_weights
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
    procedure(abstract_destructor),  pass(self), deferred, public :: destroy
    procedure(abstract_constructor), pass(self), deferred, public :: create
    procedure(abstract_description), pass(self), deferred, public :: description
    procedure(abstract_interpolate), pass(self), deferred, public :: interpolate
endtype weno_interpolator

abstract interface
  elemental subroutine abstract_destructor(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  class(weno_interpolator), intent(INOUT) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_destructor

  subroutine abstract_constructor(self, constructor, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create a WENO interpolator.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_constructor, weno_interpolator, weno_weights, weno_optimal_weights, weno_IS, weno_poly_coefficients
  class(weno_interpolator),      intent(INOUT) :: self              !< WENO interpolator.
  class(weno_constructor),       intent(IN)    :: constructor       !< WENO constructor.
  class(weno_weights),           intent(IN)    :: weights           !< WENO weights.
  class(weno_optimal_weights),   intent(IN)    :: optimal_weights   !< WENO optimal weights.
  class(weno_IS),                intent(IN)    :: IS_coefficients   !< WENO smoothness indicators coefficients.
  class(weno_poly_coefficients), intent(IN)    :: poly_coefficients !< WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_constructor

  pure subroutine abstract_description(self, weights, weno_IS, poly_coefficients, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator, weno_weights, weno_optimal_weights, weno_IS, weno_poly_coefficients
  class(weno_interpolator),      intent(IN)  :: self                      !< WENO interpolator.
  class(weno_weights),           intent(IN)  :: weights_string            !< WENO weights.
  class(weno_IS),                intent(IN)  :: IS_string                 !< WENO smoothness indicators.
  class(weno_poly_coefficients), intent(IN)  :: poly_coefficients_string  !< WENO smoothness indicators.
  character(len=:), allocatable, intent(OUT) :: string                    !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_description

  pure subroutine abstract_interpolate(self, weights, S, stencil, location, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stencil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: I_P, R_P, weno_interpolator, weno_weights, weno_optimal_weights, weno_IS, weno_poly_coefficients
  class(weno_interpolator),      intent(IN)  :: self                !< WENO interpolator.
  class(weno_weights),           intent(IN)  :: weights             !< WENO weights.
  class(weno_optimal_weights),   intent(IN)  :: optimal_weights     !< WENO optimal weights.
  class(weno_IS),                intent(IN)  :: IS_coefficients, IS !< WENO IS coefficients and values.
  class(weno_poly_coefficients), intent(IN)  :: poly_coefficients   !< WENO polynomial coefficients.
  real(R_P),                     intent(IN)  :: location            !< Location of the interpolation.
  real(R_P),                     intent(OUT) :: interpolation(1:)   !< Result of the interpolation, [1:2].
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_interpolate
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_interpolator
