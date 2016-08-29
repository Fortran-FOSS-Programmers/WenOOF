module type_weno_alpha_coefficient
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
public :: weno_alpha_coefficient
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_alpha_coefficient_constructor
  !< Abstract type used for create new concrete WENO alpha coefficient.
  !<
  !< @note Every concrete WENO alpha coefficient implementation must define its own constructor type.
  private
endtype weno_alpha_coefficient_constructor

type, abstract :: weno_alpha_coefficient
  !< WENO weights.
  !<
  !< @note Do not implement any real alpha coefficient provide the interface for the different alpha coefficient implemented.
  private
  contains
    procedure(alpha_coefficient_abstract_description), pass(self), deferred, public :: alpha_coefficient_description
    procedure(alpha_coefficient_abstract_compute),     pass(self), deferred, public :: alpha_coefficient_compute
endtype weno_alpha_coefficient

abstract interface

  pure subroutine alpha_coefficient_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_alpha_coefficient
  class(weno_alpha_coefficient),           intent(IN)  :: self   !< WENO alpha coefficient.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine alpha_coefficient_abstract_description

  pure function alpha_coefficient_abstract_compute(self, weights)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_alpha_coefficient, I_P, R_P
  class(weno_alpha_coefficient), intent(IN)  :: self        !< WENO alpha coefficient.
  real(R_P),           intent(OUT) :: weights(:)  !< Alpha coefficient of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction alpha_coefficient_abstract_compute
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_alpha_coefficient
