!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
module wenoof_weights_int_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
!<
!< @note The provided WENO weights implements the weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P, str
use wenoof_alpha_factory,  only : alpha_factory
use wenoof_alpha_object, only : alpha_object, alpha_object_constructor
use wenoof_base_object, only : base_object, base_object_constructor
use wenoof_beta_factory, only : beta_factory
use wenoof_beta_object, only : beta_object, beta_object_constructor
use wenoof_kappa_factory, only : kappa_factory
use wenoof_kappa_int_js, only : kappa_int_js
use wenoof_kappa_object, only : kappa_object, kappa_object_constructor
use wenoof_weights_object, only : weights_object, weights_object_constructor

implicit none
private
public :: weights_int_js
public :: weights_int_js_constructor

type, extends(weights_object_constructor) :: weights_int_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights object constructor.
  class(alpha_object_constructor), allocatable :: alpha_constructor !< Alpha coefficients (non linear weights) constructor.
  class(beta_object_constructor),  allocatable :: beta_constructor  !< Beta coefficients (smoothness indicators) constructor.
  class(kappa_object_constructor), allocatable :: kappa_constructor !< kappa coefficients (optimal, linear weights) constructor.
  contains
    ! public deferred methods
    procedure, pass(lhs) :: constr_assign_constr !< `=` operator.
endtype weights_int_js_constructor

type, extends(weights_object):: weights_int_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet weights object.
  !<
  !< @note The provided WENO weights implements the weights defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  class(alpha_object), allocatable :: alpha !< Alpha coefficients (non linear weights).
  class(beta_object),  allocatable :: beta  !< Beta coefficients (smoothness indicators).
  class(kappa_object), allocatable :: kappa !< kappa coefficients (optimal, linear weights).
  contains
    ! deferred public methods
    procedure, pass(self) :: create                    !< Create weights.
    procedure, pass(self) :: compute_int               !< Compute weights (interpolate).
    procedure, pass(self) :: compute_rec               !< Compute weights (reconstruct).
    procedure, pass(self) :: description               !< Return object string-description.
    procedure, pass(self) :: destroy                   !< Destroy weights.
    procedure, pass(self) :: smoothness_indicators_int !< Return smoothness indicators (interpolate).
    procedure, pass(self) :: smoothness_indicators_rec !< Return smoothness indicators (reconstrcut).
    procedure, pass(lhs)  :: object_assign_object      !< `=` operator.
endtype weights_int_js

contains
  ! constructor

  ! deferred public methods
  subroutine constr_assign_constr(lhs, rhs)
  !< `=` operator.
  class(weights_int_js_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(weights_int_js_constructor)
     if (allocated(rhs%alpha_constructor)) then
        if (.not.allocated(lhs%alpha_constructor)) allocate(lhs%alpha_constructor, mold=rhs%alpha_constructor)
           lhs%alpha_constructor = rhs%alpha_constructor
     else
        if (allocated(lhs%alpha_constructor)) deallocate(lhs%alpha_constructor)
     endif
     if (allocated(rhs%beta_constructor)) then
        if (.not.allocated(lhs%beta_constructor)) allocate(lhs%beta_constructor, mold=rhs%beta_constructor)
           lhs%beta_constructor = rhs%beta_constructor
     else
        if (allocated(lhs%beta_constructor)) deallocate(lhs%beta_constructor)
     endif
     if (allocated(rhs%kappa_constructor)) then
        if (.not.allocated(lhs%kappa_constructor)) allocate(lhs%kappa_constructor, mold=rhs%kappa_constructor)
           lhs%kappa_constructor = rhs%kappa_constructor
     else
        if (allocated(lhs%kappa_constructor)) deallocate(lhs%kappa_constructor)
     endif
  endselect
  endsubroutine constr_assign_constr

  ! deferred public methods
  subroutine create(self, constructor)
  !< Create reconstructor.
  class(weights_int_js),          intent(inout) :: self        !< Weights.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.
  type(alpha_factory)                           :: a_factory   !< Alpha factory.
  type(beta_factory)                            :: b_factory   !< Beta factory.
  type(kappa_factory)                           :: k_factory   !< Kappa factory.

  call self%destroy
  call self%create_(constructor=constructor)
  select type(constructor)
  type is(weights_int_js_constructor)
    associate(alpha_constructor=>constructor%alpha_constructor, &
              beta_constructor=>constructor%beta_constructor,   &
              kappa_constructor=>constructor%kappa_constructor)
      call a_factory%create(constructor=alpha_constructor, object=self%alpha)
      call b_factory%create(constructor=beta_constructor, object=self%beta)
      call k_factory%create(constructor=kappa_constructor, object=self%kappa)
    endassociate
  endselect
  endsubroutine create

  pure subroutine compute_int(self, stencil, values)
  !< Compute weights.
  class(weights_int_js), intent(in)  :: self               !< Weights.
  real(R_P),             intent(in)  :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(R_P),             intent(out) :: values(0:)         !< Weights values.
  real(R_P)                          :: alpha(0:self%S-1)  !< Aplha values.
  real(R_P)                          :: beta(0:self%S-1)   !< Beta values.
  real(R_P)                          :: alpha_sum          !< Sum of aplha values.
  integer(I_P)                       :: s                  !< Counters.

  call self%beta%compute(stencil=stencil, values=beta)
  select type(kappa => self%kappa)
  class is(kappa_int_js)
    call self%alpha%compute(beta=beta, kappa=kappa%values, values=alpha)
  endselect
  alpha_sum = sum(alpha)
  do s=0, self%S - 1 ! stencils loop
    values(s) = alpha(s) / alpha_sum
  enddo
  endsubroutine compute_int

  pure subroutine compute_rec(self, stencil, values)
  !< Compute weights.
  class(weights_int_js), intent(in)  :: self                  !< Weights.
  real(R_P),             intent(in)  :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(R_P),             intent(out) :: values(1:,0:)         !< Weights values of stencil interpolations.
  ! empty procedure
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(weights_int_js), intent(in)           :: self             !< Weights.
  character(len=*),      intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable               :: string           !< String-description.
  character(len=:), allocatable               :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                 :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu weights object for interpolation:'//NL
  string = string//prefix_//'  - S   = '//trim(str(self%S))//NL
  string = string//prefix_//self%alpha%description(prefix=prefix_//'  ')
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy weights.
  class(weights_int_js), intent(inout) :: self !< Weights.

  call self%destroy_
  if (allocated(self%alpha)) deallocate(self%alpha)
  if (allocated(self%beta)) deallocate(self%beta)
  if (allocated(self%kappa)) deallocate(self%kappa)
  endsubroutine destroy

  pure subroutine smoothness_indicators_int(self, si)
  !< Return smoothness indicators (interpolate).
  class(weights_int_js),  intent(in)  :: self  !< Weights.
  real(R_P),              intent(out) :: si(:) !< Smoothness indicators.
  ! TODO implement this
  endsubroutine smoothness_indicators_int

  pure subroutine smoothness_indicators_rec(self, si)
  !< Return smoothness indicators (reconstruct).
  class(weights_int_js),  intent(in)  :: self    !< Weights.
  real(R_P),              intent(out) :: si(:,:) !< Smoothness indicators.
  ! empty procedure
  endsubroutine smoothness_indicators_rec

  pure subroutine object_assign_object(lhs, rhs)
  !< `=` operator.
  class(weights_int_js), intent(inout) :: lhs !< Left hand side.
  class(base_object),    intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(weights_int_js)
     if (allocated(rhs%alpha)) then
        if (.not.allocated(lhs%alpha)) allocate(lhs%alpha, mold=rhs%alpha)
        lhs%alpha = rhs%alpha
     else
        if (allocated(lhs%alpha)) deallocate(lhs%alpha)
     endif
     if (allocated(rhs%beta)) then
        if (.not.allocated(lhs%beta)) allocate(lhs%beta, mold=rhs%beta)
        lhs%beta = rhs%beta
     else
        if (allocated(lhs%beta)) deallocate(lhs%beta)
     endif
     if (allocated(rhs%kappa)) then
        if (.not.allocated(lhs%kappa)) allocate(lhs%kappa, mold=rhs%kappa)
        lhs%kappa = rhs%kappa
     else
        if (allocated(lhs%kappa)) deallocate(lhs%kappa)
     endif
  endselect
  endsubroutine object_assign_object
endmodule wenoof_weights_int_js
