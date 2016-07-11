module type_weno_optimal_weights
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO optimal weights object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_optimal_weights
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_optimal_weights
  !< WENO optimal weights.
  !<
  !< @note Do not implement any real optimal weight: provide the interface for the different optimal weights implemented.
  private
  contains
    procedure(optimal_weights_abstract_description), pass(self), deferred, public :: optimal_weights_description
    procedure(optimal_weights_abstract_set),         pass(self), deferred, public :: optimal_weights_set
endtype weno_optimal_weights

abstract interface

  pure subroutine optimal_weights_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_weights
  class(weno_weights),           intent(IN)  :: self   !< WENO smoothness indicator.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine optimal_weights_abstract_description

  pure subroutine optimal_weights_abstract_set(self, optimal_weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the optimal weights of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_weights, I_P, R_P
  class(weno_weights), intent(IN)  :: self                !< WENO smoothness_indicator.
  real(R_P),           intent(OUT) :: optimal_weights(:)  !< Optimal weights of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine optimal_weights_abstract_set
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_optimal_weights
