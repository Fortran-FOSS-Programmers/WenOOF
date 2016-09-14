module type_weno_alpha_coefficient_js
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Jiang-Shu alpha coefficient for WENO schemes.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
use type_weno_alpha_coefficient
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_alpha_coefficient_js
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_alpha_coefficient) :: weno_alpha_coefficient_js
  !< Jiang-Shu WENO alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  private
  contains
    ! deferred public methods
    procedure, pass(self), public :: description
    procedure, pass(self), public :: compute
endtype weno_alpha_coefficient_js
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! deferred public methods
  pure subroutine description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_js), intent(in)  :: self   !< WENO alpha coefficient.
  character(len=:), allocatable,    intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO alpha coefficient'//nl
  string = string//'  Based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    alpha(S,weigt_opt,IS,eps)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(IN), the number of the stencils used'//nl
  string = string//'    weight_opt: real(R_P), intent(IN), the optimal weight of the actual stencil'//nl
  string = string//'    IS: real(R_P), intent(IN), the smoothness indicator of the actual stencil'//nl
  string = string//'    eps: real(R_P), intent(IN), the coefficient to avoid zero division used'//nl
  string = string//'    alpha: realR_P, intent(OUT), the alpha value'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure function compute(self, S, weight_opt, IS, eps) result(alpha)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_js), intent(in) :: self       !< WENO alpha coefficient.
  integer(I_P),                     intent(in) :: S          !< Number of stencils used.
  real(R_P),                        intent(in) :: weight_opt !< Optimal weight of the stencil.
  real(R_P),                        intent(in) :: IS         !< Smoothness indicator of the stencil.
  real(R_P),                        intent(in) :: eps        !< Parameter for avoiding divided by zero.
  real(R_P),                                   :: alpha      !< Alpha coefficient of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  alpha = weight_opt * (1._R_P/(eps + IS)**S)
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute

endmodule type_weno_alpha_coefficient_js
