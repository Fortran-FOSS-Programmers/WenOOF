!< Borges alpha coefficients (non linear weights) object.
module wenoof_alpha_rec_z
!< Borges alpha coefficients (non linear weights) object.
!<
!< @note The provided WENO alpha coefficients implements the alpha coefficients defined in *An improved weighted essentially
!< non-oscillatory scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP, 2008,
!< vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.

use penf, only : I_P, R_P
use wenoof_alpha_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_rec_z
public :: alpha_rec_z_constructor

type, extends(alpha_object_constructor) :: alpha_rec_z_constructor
  !< Borges alpha coefficients (non linear weights) object constructor.
endtype alpha_rec_z_constructor

type, extends(alpha_object) :: alpha_rec_z
  !< Borges alpha coefficients (non linear weights) object.
  !<
  !< @note The provided alpha coefficients implements the alpha coefficients defined in *An improved weighted essentially
  !< non-oscillatory scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP,
  !< 2008, vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.
  contains
    ! public deferred methods
    procedure, pass(self) :: compute     !< Compute coefficients.
    procedure, nopass     :: description !< Return coefficients string-description.
endtype alpha_rec_z
contains
  ! public deferred methods
  pure subroutine compute(self, beta, kappa)
  !< Compute alpha coefficients.
  class(alpha_rec_z),  intent(inout) :: self  !< Alpha coefficients.
  class(beta_object),  intent(in)    :: beta  !< Beta coefficients.
  class(kappa_object), intent(in)    :: kappa !< Kappa coefficients.
  integer(I_P)                       :: f, s1 !< Counters.

  self%values_sum = 0._R_P
  do s1=0, self%S - 1 ! stencil loops
    do f=self%f1, self%f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%values(f, s1) = kappa%values(f, s1) * &
                           ((1._R_P + (tau(S=self%S, beta=beta%values) / (self%eps + beta%values(f, s1)))) ** (weno_exp(self%S)))
      self%values_sum(f) = self%values_sum(f) + self%values(f, s1)
    enddo
  enddo
  endsubroutine compute

  pure function description(self) result(string)
  !< Return alpha coefficients string-descripition.
  class(alpha_rec_z), intent(in) :: self             !< Alpha coefficients.
  character(len=:), allocatable  :: string           !< String-description.
  character(len=1), parameter    :: nl=new_line('a') !< New line character.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'alpha_rec_z%description to be implemented, do not use!'
#endif
  endfunction description

  ! private non TBP
  pure function tau(S, beta) result(w_tau)
  !< Compute the tau coefficient used in the WENO-Z alpha coefficients.
  integer(I_P), intent(in) :: S           !< Number of stencils used.
  real(R_P),    intent(in) :: beta(0:S-1) !< Smoothness indicators.
  real(R_P)                :: w_tau       !< Tau coefficient.

  w_tau = abs(beta(0) -                           &
              (1_I_P - weno_odd(S)) * beta(1) -   &
              (1_I_P - weno_odd(S)) * beta(S-2) + &
              (1_I_P - 2_I_P * weno_odd(S)) * beta(S-1))
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
endmodule wenoof_alpha_rec_z
