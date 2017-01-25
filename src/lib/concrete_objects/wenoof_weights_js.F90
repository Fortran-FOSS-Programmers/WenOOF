!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
module wenoof_weights_js
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
use wenoof_beta_object
use wenoof_beta_rec_js
use wenoof_kappa_object
use wenoof_kappa_rec_js
use wenoof_weights_object

implicit none
private
public :: weights_js
public :: weights_js_constructor

type, extends(weights_object_constructor) :: weights_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights object constructor.
  class(alpha_object_constructor), allocatable :: alpha_constructor !< Alpha coefficients (non linear weights) constructor.
  class(beta_object_constructor),  allocatable :: beta_constructor  !< Beta coefficients (smoothness indicators) constructor.
  class(kappa_object_constructor), allocatable :: kappa_constructor !< kappa coefficients (optimal, linear weights) constructor.
endtype weights_js_constructor

type, extends(weights_object):: weights_js
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
    procedure, pass(self) :: create                !< Create weights.
    procedure, pass(self) :: compute               !< Compute weights.
    procedure, pass(self) :: description           !< Return weights string-description.
    procedure, pass(self) :: destroy               !< Destroy weights.
    procedure, pass(self) :: smoothness_indicators !< Return smoothness indicators.
endtype weights_js

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create reconstructor.
  class(weights_js),              intent(inout) :: self        !< Weights.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.
  type(alpha_factory)                           :: factory     !< Objects factory.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._RPP
  select type(constructor)
  type is(weights_js_constructor)
    associate(alpha_constructor=>constructor%alpha_constructor, &
              beta_constructor=>constructor%beta_constructor,   &
              kappa_constructor=>constructor%kappa_constructor)

      select type(alpha_constructor)
      type is(alpha_rec_js_constructor)
        ! allocate(alpha_rec_js :: self%alpha)
        ! call self%alpha%create(constructor=alpha_constructor)
        call factory%create(constructor=alpha_constructor, object=self%alpha)
      type is(alpha_rec_m_constructor)
        ! @TODO implement this
        error stop 'alpha_rec_m to be implemented'
      type is(alpha_rec_z_constructor)
        ! @TODO implement this
        error stop 'alpha_rec_z to be implemented'
      endselect

      select type(beta_constructor)
      type is(beta_rec_js_constructor)
        allocate(beta_rec_js :: self%beta)
        call self%beta%create(constructor=beta_constructor)
      endselect

      select type(kappa_constructor)
      type is(kappa_rec_js_constructor)
        allocate(kappa_rec_js :: self%kappa)
        call self%kappa%create(constructor=kappa_constructor)
      endselect
    endassociate
  endselect
  endsubroutine create

  pure subroutine compute(self, stencil)
  !< Compute weights.
  class(weights_js), intent(inout) :: self                  !< Weights.
  real(RPP),         intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P)                     :: f, s                  !< Counters.

  call self%beta%compute(stencil=stencil)
  call self%alpha%compute(beta=self%beta, kappa=self%kappa)
  do s=0, self%S - 1 ! stencils loop
    do f=self%f1, self%f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%values(f + self%ff, s) = self%alpha%values(f, s) / self%alpha%values_sum(f)
    enddo
  enddo
  endsubroutine compute

  pure function description(self) result(string)
  !< Return string-description of weights.
  class(weights_js), intent(in) :: self             !< Weights.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line char.

  string = '  Jiang-Shu weights:'//nl
  string = string//'    - S   = '//trim(str(self%S))//nl
  string = string//'    - f1  = '//trim(str(self%f1))//nl
  string = string//'    - f2  = '//trim(str(self%f2))//nl
  string = string//'    - ff  = '//trim(str(self%ff))//nl
  string = string//self%alpha%description()
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy weights.
  class(weights_js), intent(inout) :: self !< Weights.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  if (allocated(self%alpha)) deallocate(self%alpha)
  if (allocated(self%beta)) deallocate(self%beta)
  if (allocated(self%kappa)) deallocate(self%kappa)
  endsubroutine destroy

  pure function smoothness_indicators(self) result(si)
  !< Return smoothness indicators..
  class(weights_js), intent(in) :: self    !< Weights.
  real(RPP), allocatable        :: si(:,:) !< Smoothness indicators.

  if (allocated(self%beta)) then
    if (allocated(self%beta%values)) then
      si = self%beta%values
    endif
  endif
  endfunction smoothness_indicators
endmodule wenoof_weights_js
