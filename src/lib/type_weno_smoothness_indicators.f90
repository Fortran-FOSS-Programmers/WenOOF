module type_weno_smoothness_indicators
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO smoothness indicators object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_IS
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_IS
  !< WENO smoothness indicators.
  !<
  !< @note Do not implement any real smoothness indicators: provide the interface for the different smoothness_indicators implemented.
  private
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destructor
    procedure(constructor_interface), pass(self), deferred, public :: constructor
    procedure(description_interface), pass(self), deferred, public :: description
    procedure(compute_interface),     pass(self), deferred, public :: compute
endtype weno_IS

abstract interface
  !< Destroy WENO polynomial coefficients.
  pure subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS
  class(weno_IS), intent(inout) :: self   !< WENO smoothenss indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  !< Create WENO polynomial coefficients.
  pure subroutine constructor_interface(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO smoothness indicators coefficients.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS, I_P
  class(weno_IS), intent(inout) :: self        !< WENO smoothness indicators.
  integer(I_P),   intent(in)    :: S           !< Number of stencils used.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Return a string describing WENO smoothness_indicators.
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO smoothness_indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS
  class(weno_IS),                intent(in)  :: self   !< WENO smoothness indicator.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  !< Compute the smoothness indicators of the WENO interpolating polynomial.
  pure function compute_interface(self, smooth_coef, v1, v2) result(IS)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the partial value of the smoothness indicator of a single WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS, I_P, R_P
  class(weno_IS), intent(in) :: self        !< WENO smoothness indicator.
  real(R_P),      intent(in) :: smooth_coef !< Coefficient of the smoothness indicator.
  real(R_P),                 :: v1          !< First (pivotal) value from the stencil used for the interpolation.
  real(R_P),                 :: v2          !< Second value from the stencil used for the interpolation.
  real(R_P),                 :: IS          !< Partial value of the smoothness indicator of polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_smoothness_indicators
