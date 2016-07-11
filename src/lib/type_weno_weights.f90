module type_weno_weights
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO weights object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_weights
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_weights
  !< WENO weights.
  !<
  !< @note Do not implement any real weight: provide the interface for the different weights implemented.
  private
  contains
    procedure(weights_abstract_description), pass(self), deferred, public :: weights_description
    procedure(weights_abstract_compute),     pass(self), deferred, public :: weights_compute
endtype weno_weights

abstract interface

  pure subroutine weights_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_weights
  class(weno_weights),           intent(IN)  :: self   !< WENO smoothness indicator.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine weights_abstract_description

  pure subroutine weights_abstract_compute(self, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the weights of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_weights, I_P, R_P
  class(weno_weights), intent(IN)  :: self        !< WENO smoothness_indicator.
  real(R_P),           intent(OUT) :: weights(:)  !< Weights of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine weights_abstract_compute
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_weights
