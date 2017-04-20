!< Henrick alpha (non linear weights) object.
module wenoof_alpha_int_m
!< Henrick alpha (non linear weights) object.
!<
!< @note The provided alpha implements the alpha coefficients defined in *Mapped weighted essentially non-oscillatory schemes:
!< Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers, JCP,
!< 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023

use penf, only : I_P, R_P, str
use wenoof_alpha_object, only : alpha_object, alpha_object_constructor
use wenoof_alpha_rec_js, only : alpha_rec_js, alpha_rec_js_constructor
use wenoof_alpha_rec_z, only : alpha_rec_z, alpha_rec_z_constructor
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: alpha_int_m
public :: alpha_int_m_constructor

type, extends(alpha_object_constructor) :: alpha_int_m_constructor
  !< Henrick alpha (non linear weights) object constructor.
  character(len=:), allocatable :: base_type !< Base alpha coefficient type.
  contains
    ! public deferred methods
    procedure, pass(lhs) :: constr_assign_constr !< `=` operator.
endtype alpha_int_m_constructor

type, extends(alpha_object) :: alpha_int_m
  !< Henrick alpha (non linear weights) object.
  !<
  !< @note The provided alpha implements the alpha coefficients defined in *Mapped weighted essentially non-oscillatory schemes:
  !< Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
  !< JCP, 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023.
  class(alpha_object), allocatable :: alpha_base !< Base alpha to be re-mapped.
  contains
    ! public deferred methods
    procedure, pass(self) :: create               !< Create alpha.
    procedure, pass(self) :: compute_int          !< Compute alpha (interpolate).
    procedure, pass(self) :: compute_rec          !< Compute alpha (reconstruct).
    procedure, pass(self) :: description          !< Return object string-description.
    procedure, pass(self) :: destroy              !< Destroy alpha.
    procedure, pass(lhs)  :: object_assign_object !< `=` operator.
endtype alpha_int_m

contains
  ! constructor

  ! deferred public methods
  subroutine constr_assign_constr(lhs, rhs)
  !< `=` operator.
  class(alpha_int_m_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor), intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(alpha_int_m_constructor)
     if (allocated(rhs%base_type)) lhs%base_type = rhs%base_type
  endselect
  endsubroutine constr_assign_constr

  ! deferred public methods
  subroutine create(self, constructor)
  !< Create alpha.
  class(alpha_int_m),             intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  select type(constructor)
  type is(alpha_int_m_constructor)
    if (allocated(constructor%base_type)) then
      select case(constructor%base_type)
      case('JS')
        if (allocated(self%alpha_base)) deallocate(self%alpha_base)
        allocate(alpha_rec_js :: self%alpha_base)
        call self%alpha_base%create(constructor=constructor)
      case('Z')
        if (allocated(self%alpha_base)) deallocate(self%alpha_base)
        allocate(alpha_rec_z :: self%alpha_base)
        call self%alpha_base%create(constructor=constructor)
      endselect
    endif
  class default
    ! TODO add error handling
  endselect
  endsubroutine create

  pure subroutine compute_int(self, beta, kappa, values)
  !< Compute alpha (interpolate).
  class(alpha_int_m),  intent(in)  :: self                   !< Alpha.
  real(R_P),           intent(in)  :: beta(0:)               !< Beta [0:S-1].
  real(R_P),           intent(in)  :: kappa(0:)              !< Kappa [0:S-1].
  real(R_P),           intent(out) :: values(0:)             !< Alpha values [0:S-1].
  real(R_P)                        :: alpha_base(0:self%S-1) !< Alpha base coefficients.
  real(R_P)                        :: alpha_base_sum         !< Sum of alpha base coefficients.
  real(R_P)                        :: kappa_base             !< Kappa base coefficient.
  integer(I_P)                     :: s1                     !< Counter.

  call self%alpha_base%compute(beta=beta, kappa=kappa, values=alpha_base)
  alpha_base_sum = sum(alpha_base)
  do s1=0, self%S - 1 ! stencil loops
    kappa_base = alpha_base(s1) / alpha_base_sum
    values(s1) =                                                                                                       &
      (kappa_base * (kappa(s1) + kappa(s1) * kappa(s1) - 3._R_P * kappa(s1) * kappa_base + kappa_base * kappa_base)) / &
      (kappa(s1) * kappa(s1) + kappa_base * (1._R_P - 2._R_P * kappa(s1)))
  enddo
  endsubroutine compute_int

  pure subroutine compute_rec(self, beta, kappa, values)
  !< Compute alpha (reconstruct).
  class(alpha_int_m), intent(in)  :: self          !< Alpha coefficient.
  real(R_P),          intent(in)  :: beta(1:,0:)   !< Beta [1:2,0:S-1].
  real(R_P),          intent(in)  :: kappa(1:,0:)  !< Kappa [1:2,0:S-1].
  real(R_P),          intent(out) :: values(1:,0:) !< Alpha values [1:2,0:S-1].
  ! empty procedure
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(alpha_int_m), intent(in)           :: self             !< Alpha coefficient.
  character(len=*),   intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable            :: string           !< String-description.
  character(len=:), allocatable            :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter              :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu alpha coefficients object for interpolation:'//NL
  string = string//prefix_//'  - S   = '//trim(str(self%S))//NL
  string = string//prefix_//'  - eps = '//trim(str(self%eps))//NL
  associate(alpha_base=>self%alpha_base)
    select type(alpha_base)
    type is(alpha_rec_js)
      string = string//prefix_//'  - base-mapped-alpha type = Jiang-Shu'
    type is(alpha_rec_z)
      string = string//prefix_//'  - base-mapped-alpha type = Bogeg'
    endselect
  endassociate
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy alpha.
  class(alpha_int_m), intent(inout) :: self !< Alpha.

  call self%destroy_
  if (allocated(self%alpha_base)) deallocate(self%alpha_base)
  endsubroutine destroy

  subroutine object_assign_object(lhs, rhs)
  !< `=` operator.
  class(alpha_int_m), intent(inout) :: lhs !< Left hand side.
  class(base_object), intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(alpha_int_m)
     if (allocated(rhs%alpha_base)) then
        if (.not.allocated(lhs%alpha_base)) allocate(lhs%alpha_base, mold=rhs%alpha_base)
        lhs%alpha_base = rhs%alpha_base
     else
        if (allocated(lhs%alpha_base)) deallocate(lhs%alpha_base)
     endif
  endselect
  endsubroutine object_assign_object
endmodule wenoof_alpha_int_m
