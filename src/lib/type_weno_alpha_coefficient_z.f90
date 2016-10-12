module type_weno_alpha_coefficient_z
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Borges et al. alpha coefficient for WENO schemes.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *An improved weighted essentially
!< non-oscillatory scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP, 2008,
!< vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P
use type_weno_alpha_coefficient
use type_weno_alpha_coefficient_js
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_alpha_coefficient_z, associate_WENO_alpha_z
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_alpha_coefficient_js) :: weno_alpha_coefficient_z
  !< Borges et al. WENO alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *An improved weighted essentially
  !< non-oscillatory scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP, 2008,
  !< vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.
  private
  contains
    ! deferred public methods
    procedure, pass(self), public :: description
    procedure, pass(self), public :: compute
    ! public methods
    procedure, nopass,     public :: tau
    procedure, nopass,     public :: weno_exp
    procedure, nopass,     public :: weno_odd
endtype weno_alpha_coefficient_z
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public, non TBP
  function associate_WENO_alpha_z(alpha_input) result(alpha_pointer)
    !< Check the type of the alpha coefficient passed as input and return a WENO Z alpha coefficient associated to the alpha coefficient.
    class(weno_alpha_coefficient), intent(in), target  :: alpha_input   !< Input alpha coefficient.
    class(weno_alpha_coefficient_z),           pointer :: alpha_pointer !< WENO Z alpha coefficients.

    select type(alpha_input)
      type is(weno_alpha_coefficient_z)
        alpha_pointer => alpha_input
      class default
        write(stderr, '(A)')'error: wrong alpha coefficient type chosen'
        stop
    end select
  end function associate_WENO_alpha_z

  ! deferred public methods
  pure subroutine description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_z), intent(in)  :: self   !< WENO alpha coefficient.
  character(len=:), allocatable,   intent(out) :: string !< String returned.
  character(len=1), parameter                  :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO alpha coefficient'//nl
  string = string//'  Based on the work by Borges, Carmona, Costa and Don "An improved weighted essentially non-oscillatory '// &
           'scheme for hyperbolic conservation laws", see '// &
           'JCP, 2008, vol. 227, pp. 3191--3211, doi:10.1016/j.jcp.2007.11.038'//nl
  string = string//'  The "alpha" method has the following public API'//nl
  string = string//'    alpha(S,weigt_opt,IS,eps)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(IN), the number of the stencils used'//nl
  string = string//'    weight_opt: real(R_P), intent(IN), the optimal weight of the actual stencil'//nl
  string = string//'    IS: real(R_P), intent(IN), the smoothness indicator of the actual stencil'//nl
  string = string//'    eps: real(R_P), intent(IN), the coefficient to avoid zero division used'//nl
  string = string//'    alpha: realR_P, intent(OUT), the alpha value'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure function compute(self, S, weight_opt, IS, IS_i, eps) result(alpha)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_z), intent(in)           :: self         !< WENO alpha coefficient.
  integer(I_P),                    intent(in)           :: S            !< Number of stencils used.
  real(R_P),                       intent(in)           :: weight_opt   !< Optimal weight of the stencil.
  real(R_P),                       intent(in), optional :: IS(0:S - 1)  !< Smoothness indicators of the stencils.
  real(R_P),                       intent(in)           :: IS_i         !< Smoothness indicator of the i-th stencil.
  real(R_P),                       intent(in)           :: eps          !< Parameter for avoiding divided by zero.
  real(R_P)                                             :: alpha        !< Alpha coefficient of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  alpha = weight_opt * ((1._R_P + (tau(S,IS)/(eps+IS_i))) ** (weno_exp(S)))
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute

  pure function weno_exp(S) result(w_exp)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the exponent used in the alpha function.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(in) :: S     !< Number of stencils used.
  integer(I_P)             :: w_exp !< Exponent used in the alpha function.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  w_exp = int(S, I_P)
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction weno_exp

  pure function weno_odd(S) result(w_odd)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the distinguisher between odd and even number of stencils.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(in) :: S     !< Number of stencils used.
  integer(I_P)             :: w_odd !< Distinguishing between odd and even number of stencils.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  w_odd = int(mod(S, 2_I_P), I_P)
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction weno_odd

  pure function tau(S, IS) result(w_tau)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the tau coefficient used in the WENO-Z alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(in) :: S           !< Number of stencils used.
  real(R_P),    intent(in) :: IS(0:S - 1) !< Smoothness indicators.
  real(R_P)                :: w_tau       !< Tau coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  w_tau = abs(IS(0) - (1-weno_odd(S))*IS(1) - (1-weno_odd(S))*IS(S-2_I_P) + (1-2*weno_odd(S))*IS(S-1_I_P))
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction tau
endmodule type_weno_alpha_coefficient_z
