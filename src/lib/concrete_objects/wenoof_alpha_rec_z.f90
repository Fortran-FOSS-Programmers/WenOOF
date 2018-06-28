!< Borges alpha (non linear weights) object.
module wenoof_alpha_rec_z
!< Borges alpha (non linear weights) object.
!<
!< @note The provided WENO alpha implements the alpha coefficients defined in *An improved weighted essentially non-oscillatory
!< scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP, 2008,
!< vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.

use penf, only : I_P, R_P, str
use wenoof_alpha_object, only : alpha_object, alpha_object_constructor
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: alpha_rec_z
public :: alpha_rec_z_constructor

type, extends(alpha_object_constructor) :: alpha_rec_z_constructor
  !< Borges alpha (non linear weights) object constructor.
  contains
    ! public deferred methods
    procedure, pass(lhs) :: constr_assign_constr !< `=` operator.
endtype alpha_rec_z_constructor

type, extends(alpha_object) :: alpha_rec_z
  !< Borges alpha (non linear weights) object.
  !<
  !< @note The provided alpha implements the alpha coefficients defined in *An improved weighted essentially non-oscillatory
  !< scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP,
  !< 2008, vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.
  contains
    ! public deferred methods
    procedure, pass(self) :: create               !< Create alpha.
    procedure, pass(self) :: compute_int          !< Compute alpha (interpolate).
    procedure, pass(self) :: compute_rec          !< Compute alpha (reconstruct).
    procedure, pass(self) :: description          !< Return alpha string-description.
    procedure, pass(self) :: destroy              !< Destroy alpha.
    procedure, pass(lhs)  :: object_assign_object !< `=` operator.
endtype alpha_rec_z

contains
  ! constructor

  ! deferred public methods
  subroutine constr_assign_constr(lhs, rhs)
  !< `=` operator.
  class(alpha_rec_z_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor), intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  endsubroutine constr_assign_constr

  ! public deferred methods
  subroutine create(self, constructor)
  !< Create alpha.
  class(alpha_rec_z),             intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  endsubroutine create

  pure subroutine compute_int(self, ord, beta, kappa, values)
  !< Compute alpha (interpolate).
  class(alpha_rec_z), intent(in)  :: self       !< Alpha.
  integer(I_P),       intent(in)  :: ord        !< Order of interpolation.
  real(R_P),          intent(in)  :: beta(0:)   !< Beta [0:S-1].
  real(R_P),          intent(in)  :: kappa(0:)  !< Kappa [0:S-1].
  real(R_P),          intent(out) :: values(0:) !< Alpha values [0:S-1].
  ! empty procedure
  endsubroutine compute_int

  pure subroutine compute_rec(self, ord, beta, kappa, values)
  !< Compute alpha.
  class(alpha_rec_z), intent(in)  :: self          !< Alpha.
  integer(I_P),       intent(in)  :: ord           !< Order of reconstruction.
  real(R_P),          intent(in)  :: beta(1:,0:)   !< Beta [1:2,0:S-1].
  real(R_P),          intent(in)  :: kappa(1:,0:)  !< Kappa [1:2,0:S-1].
  real(R_P),          intent(out) :: values(1:,0:) !< Alpha values [1:2,0:S-1].
  integer(I_P)                    :: f, s1         !< Counters.

  do s1=0, ord - 1 ! stencil loops
    do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      values(f, s1) = kappa(f, s1) * ((1._R_P + (tau(S=ord, beta=beta) / (self%eps + beta(f, s1)))) ** (weno_exp(ord)))
    enddo
  enddo
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(alpha_rec_z), intent(in)           :: self             !< Alpha coefficient.
  character(len=*),   intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable            :: string           !< String-description.
  character(len=:), allocatable            :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter              :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Borges alpha coefficients object for reconstruction:'//NL
  string = string//prefix_//'  - S   = '//trim(str(self%S))//NL
  string = string//prefix_//'  - eps = '//trim(str(self%eps))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy alpha.
  class(alpha_rec_z), intent(inout) :: self !< Alpha.

  call self%destroy_
  endsubroutine destroy

  pure subroutine object_assign_object(lhs, rhs)
  !< `=` operator.
  class(alpha_rec_z), intent(inout) :: lhs !< Left hand side.
  class(base_object), intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  endsubroutine object_assign_object

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
