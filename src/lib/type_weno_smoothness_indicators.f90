module type_weno_IS
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
  !< @note Do not implement any real smoothness indicator: provide the interface for the different smoothness_indicators implemented.
  private
  contains
    procedure(IS_abstract_description),      pass(self), deferred, public :: IS_description
    procedure(IS_abstract_set_coefficients), pass(self), deferred, public :: IS_set_coefficients
    procedure(IS_abstract_compute),          pass(self), deferred, public :: IS_compute
endtype weno_IS

abstract interface

  pure subroutine IS_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO smoothness_indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS
  class(weno_IS),                intent(IN)  :: self   !< WENO smoothness indicator.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IS_abstract_description

  pure subroutine IS_abstract_set_coefficients(self, IS_coefficients)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the smoothness_indicators of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS, I_P, R_P
  class(weno_IS), intent(IN)  :: self                !< WENO smoothness_indicator.
  real(R_P),      intent(OUT) :: IS_coefficients(:)  !< smoothness indicators coefficients,
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IS_abstract_set_coefficients

  pure subroutine IS_abstract_compute(self, IS)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the smoothness_indicators of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS, I_P, R_P
  class(weno_IS), intent(IN)  :: self   !< WENO smoothness_indicator.
  real(R_P),      intent(OUT) :: IS(:)  !< smoothness indicators of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IS_abstract_compute
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_smoothness_indicators
