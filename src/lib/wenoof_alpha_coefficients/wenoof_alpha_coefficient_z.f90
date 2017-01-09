!< Borges WENO alpha coefficient object.
module wenoof_alpha_coefficient_z
!< Borges WENO alpha coefficient object.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *An improved weighted essentially
!< non-oscillatory scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP, 2008,
!< vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.

use penf, only : I_P, R_P
use wenoof_alpha_coefficient_abstract
use wenoof_alpha_coefficient_js

implicit none
private
public :: alpha_coefficient_z

type, extends(alpha_coefficient_js) :: alpha_coefficient_z
  !< Borges WENO alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *An improved weighted essentially
  !< non-oscillatory scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP,
  !< 2008, vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.
  private
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute coefficients.
    procedure, nopass     :: description !< Return string-description of coefficients.
    ! public methods
    procedure, nopass :: tau      !< Compute the tau coefficient used in the WENO-Z alpha coefficient.
    procedure, nopass :: weno_exp !< Compute the exponent used in the alpha function.
    procedure, nopass :: weno_odd !< Compute the distinguisher between odd and even number of stencils.
endtype alpha_coefficient_z
contains
  ! deferred public methods
  pure subroutine description(string)
  !< Return string-descripition of coefficients.
  !<
  !< @TODO change to function.
  character(len=:), allocatable, intent(out) :: string           !< String returned.
  character(len=1), parameter                :: nl=new_line('a') !< New line character.

  string = 'WENO alpha coefficient'//nl
  string = string//'  Based on the work by Borges, Carmona, Costa and Don "An improved weighted essentially non-oscillatory '// &
           'scheme for hyperbolic conservation laws", see '// &
           'JCP, 2008, vol. 227, pp. 3191--3211, doi:10.1016/j.jcp.2007.11.038'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    compute(S, weigt_opt, IS, eps, f1, f2)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(in), the number of the stencils used'//nl
  string = string//'    weight_opt: real(R_P), intent(in), the optimal weight of the actual stencil'//nl
  string = string//'    IS: real(R_P), intent(in), the smoothness indicator of the actual stencil'//nl
  string = string//'    eps: real(R_P), intent(in), the coefficient to avoid zero division used'//nl
  string = string//'    f1, f2: integer(I_P), intent(in), the faces to be computed (1 => left interface, 2 => right interface)'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute coefficients.
  class(alpha_coefficient_z), intent(inout) :: self                        !< WENO alpha coefficient.
  integer(I_P),               intent(in)    :: S                           !< Number of stencils used.
  real(R_P),                  intent(in)    :: weight_opt(1: 2, 0 : S - 1) !< Optimal weight of the stencil.
  real(R_P),                  intent(in)    :: IS(1: 2, 0 : S - 1)         !< Smoothness indicators of the stencils.
  real(R_P),                  intent(in)    :: eps                         !< Parameter for avoiding divided by zero.
  integer(I_P),               intent(in)    :: f1, f2                      !< Faces to be computed.
  integer(I_P)                              :: f, s1                       !< Counters.

  self%alpha_tot = 0._R_P
  do s1=0, S - 1 ! stencil loops
    do f=f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%alpha_coef(f, s1) = weight_opt(f, s1) * ((1._R_P + (tau(S,IS)/(eps+IS(f, s1)))) ** (weno_exp(S)))
      self%alpha_tot(f) = self%alpha_tot(f) + self%alpha_coef(f, s1)
    enddo
  enddo
  endsubroutine compute

  pure function tau(S, IS) result(w_tau)
  !< Compute the tau coefficient used in the WENO-Z alpha coefficient.
  integer(I_P), intent(in) :: S           !< Number of stencils used.
  real(R_P),    intent(in) :: IS(0:S - 1) !< Smoothness indicators.
  real(R_P)                :: w_tau       !< Tau coefficient.

  w_tau = abs(IS(0) - (1-weno_odd(S))*IS(1) - (1-weno_odd(S))*IS(S-2_I_P) + (1-2*weno_odd(S))*IS(S-1_I_P))
  endfunction tau

  pure function weno_exp(S) result(w_exp)
  !< Compute the exponent used in the alpha function.
  integer(I_P), intent(in) :: S     !< Number of stencils used.
  integer(I_P)             :: w_exp !< Exponent used in the alpha function.

  w_exp = int(S, I_P)
  endfunction weno_exp

  pure function weno_odd(S) result(w_odd)
  !< Compute the distinguisher between odd and even number of stencils.
  integer(I_P), intent(in) :: S     !< Number of stencils used.
  integer(I_P)             :: w_odd !< Distinguishing between odd and even number of stencils.

  w_odd = int(mod(S, 2_I_P), I_P)
  endfunction weno_odd
endmodule wenoof_alpha_coefficient_z
