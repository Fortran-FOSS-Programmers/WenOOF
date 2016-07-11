module type_weno_poly_coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO polynomial coefficients object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_poly_coefficients
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_poly_coefficients
  !< WENO polynomial coefficients.
  !<
  !< @note Do not implement any real polynomial coefficient: provide the interface for the different polynomial coefficients implemented.
  private
  contains
    procedure(poly_coefficients_abstract_description), pass(self), deferred, public :: poly_coefficients_description
    procedure(poly_coefficients_abstract_set),         pass(self), deferred, public :: poly_coefficients_set
endtype weno_poly_coefficients

abstract interface

  pure subroutine poly_coefficients_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_poly_coefficients
  class(weno_poly_coefficients), intent(IN)  :: self   !< WENO polynomial coefficients.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine poly_coefficients_abstract_description

  pure subroutine poly_coefficients_abstract_set(self, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_poly _coefficients, I_P, R_P
  class(weno_poly_coefficients), intent(IN)  :: self                  !< WENO polynomial coefficients
  real(R_P),                     intent(OUT) :: poly_coefficients(:)  !< WENO polynomial coefficients for each stencil.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine poly_coefficients_abstract_set
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_poly_coefficients
