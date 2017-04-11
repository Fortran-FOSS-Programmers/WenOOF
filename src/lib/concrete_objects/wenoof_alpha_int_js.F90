!< Jiang-Shu alpha (non linear weights) object.
module wenoof_alpha_int_js
!< Jiang-Shu alpha (non linear weights) object.
!<
!< @note The provided alpha implements the alpha coefficients defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.

#ifdef r16p
use penf, only: I_P, RPP=>R16P, str
#else
use penf, only: I_P, RPP=>R8P, str
#endif
use wenoof_alpha_object, only : alpha_object, alpha_object_constructor
use wenoof_base_object, only : base_object_constructor

implicit none
private
public :: alpha_int_js
public :: alpha_int_js_constructor

type, extends(alpha_object_constructor) :: alpha_int_js_constructor
  !< Jiang-Shu alpha object constructor.
endtype alpha_int_js_constructor

type, extends(alpha_object) :: alpha_int_js
  !< Jiang-Shu alpha object.
  !<
  !< @note The provided WENO alpha implements the alpha coefficients defined in *Efficient Implementation of Weighted
  !< ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  contains
    ! public deferred methods
    procedure, pass(self) :: create      !< Create alpha.
    procedure, pass(self) :: compute_int !< Compute alpha (interpolate).
    procedure, pass(self) :: compute_rec !< Compute alpha (reconstruct).
    procedure, pass(self) :: description !< Return object string-description.
    procedure, pass(self) :: destroy     !< Destroy alpha.
endtype alpha_int_js

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create alpha.
  class(alpha_int_js),            intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  endsubroutine create

  pure subroutine compute_int(self, beta, kappa, values)
  !< Compute alpha (interpolate).
  class(alpha_int_js), intent(in)  :: self       !< Alpha coefficient.
  real(RPP),           intent(in)  :: beta(0:)   !< Beta [0:S-1].
  real(RPP),           intent(in)  :: kappa(0:)  !< Kappa [0:S-1].
  real(RPP),           intent(out) :: values(0:) !< Alpha values [0:S-1].
  integer(I_P)                     :: s1         !< Counter.

  do s1=0, self%S - 1 ! stencil loops
    values(s1) = kappa(s1) / (self%eps + beta(s1)) ** self%S
  enddo
  endsubroutine compute_int

  pure subroutine compute_rec(self, beta, kappa, values)
  !< Compute alpha (reconstruct).
  class(alpha_int_js), intent(in)  :: self          !< Alpha coefficient.
  real(RPP),           intent(in)  :: beta(1:,0:)   !< Beta [1:2,0:S-1].
  real(RPP),           intent(in)  :: kappa(1:,0:)  !< Kappa [1:2,0:S-1].
  real(RPP),           intent(out) :: values(1:,0:) !< Alpha values [1:2,0:S-1].
  ! empty procedure
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(alpha_int_js), intent(in)           :: self             !< Alpha coefficient.
  character(len=*),    intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable             :: string           !< String-description.
  character(len=:), allocatable             :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter               :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu alpha coefficients object for interpolation:'//NL
  string = string//prefix_//'  - S   = '//trim(str(self%S))//NL
  string = string//prefix_//'  - eps = '//trim(str(self%eps))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy alpha.
  class(alpha_int_js), intent(inout) :: self !< Alpha.

  call self%destroy_
  endsubroutine destroy
endmodule wenoof_alpha_int_js
