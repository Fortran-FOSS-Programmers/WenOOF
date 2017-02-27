!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
module wenoof_weights_int_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
!<
!< @note The provided WENO weights implements the weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

#ifdef r16p
use penf, only: I_P, RPP=>R16P, str
#else
use penf, only: I_P, RPP=>R8P, str
#endif
use wenoof_alpha_factory
use wenoof_alpha_object
use wenoof_alpha_rec_js
use wenoof_alpha_rec_m
use wenoof_alpha_rec_z
use wenoof_base_object
use wenoof_beta_factory
use wenoof_beta_object
use wenoof_beta_rec_js
use wenoof_kappa_factory
use wenoof_kappa_object
use wenoof_kappa_rec_js
use wenoof_weights_object

implicit none
private
public :: weights_int_js
public :: weights_int_js_constructor

type, extends(weights_object_constructor) :: weights_int_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights object constructor.
  class(alpha_object_constructor), allocatable :: alpha_constructor !< Alpha coefficients (non linear weights) constructor.
  class(beta_object_constructor),  allocatable :: beta_constructor  !< Beta coefficients (smoothness indicators) constructor.
  class(kappa_object_constructor), allocatable :: kappa_constructor !< kappa coefficients (optimal, linear weights) constructor.
endtype weights_int_js_constructor

type, extends(weights_object):: weights_int_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet weights object.
  !<
  !< @note The provided WENO weights implements the weights defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  class(alpha_object), allocatable :: alpha     !< Alpha coefficients (non linear weights).
  class(beta_object),  allocatable :: beta      !< Beta coefficients (smoothness indicators).
  class(kappa_object), allocatable :: kappa     !< kappa coefficients (optimal, linear weights).
  contains
    ! deferred public methods
    procedure, pass(self) :: create                          !< Create weights.
    procedure, pass(self) :: compute_with_stencil_of_rank_1  !< Compute weights.
    procedure, pass(self) :: compute_with_stencil_of_rank_2  !< Compute weights.
    procedure, pass(self) :: description                     !< Return weights string-description.
    procedure, pass(self) :: destroy                         !< Destroy weights.
    procedure, pass(self) :: smoothness_indicators_of_rank_1 !< Return smoothness indicators.
    procedure, pass(self) :: smoothness_indicators_of_rank_2 !< Return smoothness indicators.
endtype weights_int_js

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create reconstructor.
  class(weights_int_js),           intent(inout) :: self        !< Weights.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.
  type(alpha_factory)                           :: a_factory   !< Alpha factory.
  type(beta_factory)                            :: b_factory   !< Beta factory.
  type(kappa_factory)                           :: k_factory   !< Kappa factory.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values_rank_1(0:self%S - 1))
  self%values_rank_1 = 0._RPP
  select type(constructor)
  type is(weights_int_js_constructor)
    associate(alpha_constructor=>constructor%alpha_constructor, &
              beta_constructor=>constructor%beta_constructor,   &
              kappa_constructor=>constructor%kappa_constructor)

      call a_factory%create(constructor=alpha_constructor, object=self%alpha)
      ! select type(alpha_constructor)
      ! type is(alpha_rec_js_constructor)
      !   call factory%create(constructor=alpha_constructor, object=self%alpha)
      ! type is(alpha_rec_m_constructor)
      !   call factory%create(constructor=alpha_constructor, object=self%alpha)
      ! type is(alpha_rec_z_constructor)
      !   call factory%create(constructor=alpha_constructor, object=self%alpha)
      ! endselect

      call b_factory%create(constructor=beta_constructor, object=self%beta)
      ! select type(beta_constructor)
      ! type is(beta_rec_js_constructor)
      !   allocate(beta_rec_js :: self%beta)
      !   call self%beta%create(constructor=beta_constructor)
      ! endselect

      call k_factory%create(constructor=kappa_constructor, object=self%kappa)
      ! select type(kappa_constructor)
      ! type is(kappa_rec_js_constructor)
      !   allocate(kappa_rec_js :: self%kappa)
      !   call self%kappa%create(constructor=kappa_constructor)
      ! endselect
    endassociate
  endselect
  endsubroutine create

  pure subroutine compute_with_stencil_of_rank_1(self, stencil)
  !< Compute weights.
  class(weights_int_js), intent(inout) :: self               !< Weights.
  real(RPP),            intent(in)    :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  integer(I_P)                        :: s                  !< Counters.

  call self%beta%compute(stencil=stencil)
  call self%alpha%compute(beta=self%beta, kappa=self%kappa)
  do s=0, self%S - 1 ! stencils loop
    self%values_rank_1(s) = self%alpha%values_rank_1(s) / self%alpha%values_sum_rank_1
  enddo
  endsubroutine compute_with_stencil_of_rank_1

  pure subroutine compute_with_stencil_of_rank_2(self, stencil)
  !< Compute weights.
  class(weights_int_js), intent(inout) :: self               !< Weights.
  real(RPP),         intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].

  ! Empty routine.
  endsubroutine compute_with_stencil_of_rank_2

  pure function description(self) result(string)
  !< Return string-description of weights.
  class(weights_int_js), intent(in) :: self             !< Weights.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line char.

  string = '  Jiang-Shu weights:'//nl
  string = string//'    - S   = '//trim(str(self%S))//nl
  string = string//self%alpha%description()
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy weights.
  class(weights_int_js), intent(inout) :: self !< Weights.

  call self%destroy_
  if (allocated(self%values_rank_1)) deallocate(self%values_rank_1)
  if (allocated(self%alpha)) deallocate(self%alpha)
  if (allocated(self%beta)) deallocate(self%beta)
  if (allocated(self%kappa)) deallocate(self%kappa)
  endsubroutine destroy

  pure subroutine smoothness_indicators_of_rank_1(self, si)
  !< Return smoothness indicators..
  class(weights_int_js),  intent(in)  :: self  !< Weights.
  real(RPP),              intent(out) :: si(:) !< Smoothness indicators.

  if (allocated(self%beta)) then
    if (allocated(self%beta%values_rank_1)) then
      si = self%beta%values_rank_1
    endif
  endif
  endsubroutine smoothness_indicators_of_rank_1

  pure subroutine smoothness_indicators_of_rank_2(self, si)
  !< Return smoothness indicators..
  class(weights_int_js),  intent(in)  :: self    !< Weights.
  real(RPP),              intent(out) :: si(:,:) !< Smoothness indicators.

  ! Empty routine
  endsubroutine smoothness_indicators_of_rank_2
endmodule wenoof_weights_int_js
