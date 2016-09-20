module type_weno_alpha_coefficient
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO alpha coefficient object.
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
type, abstract :: weno_alpha_coefficient
  !< WENO weights.
  !<
  !< @note Do not implement any real alpha coefficient provide the interface for the different alpha coefficient implemented.
  private
  contains
    procedure(description_interface), pass(self), deferred, public :: description
    procedure(compute_interface),     pass(self), deferred, public :: compute
endtype weno_alpha_coefficient

abstract interface
  !< Return a string describing WENO alpha coefficient.
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_alpha_coefficient
  class(weno_alpha_coefficient), intent(in)  :: self   !< WENO alpha coefficient.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  pure function compute_interface(self, S, weight_opt, IS, eps) result(alpha)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_alpha_coefficient, I_P, R_P
  class(weno_alpha_coefficient), intent(in) :: self        !< WENO alpha coefficient.
  integer(I_P),                  intent(in) :: S           !< Number of stencils used.
  real(R_P),                     intent(in) :: weight_opt  !< Optimal weight of the stencil.
  real(R_P),                     intent(in) :: IS          !< Smoothness indicator of the stencil.
  real(R_P),                     intent(in) :: eps         !< Parameter for avoiding divided by zero.
  real(R_P)                                 :: alpha       !< Alpha coefficient of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_alpha_coefficient
