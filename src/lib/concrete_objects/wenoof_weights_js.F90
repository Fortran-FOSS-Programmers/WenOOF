!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
module wenoof_weights_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
!<
!< @note The provided WENO weights implements the weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P
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
public :: create_weights_js_constructor

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
    procedure, pass(self) :: create      !< Create weights.
    procedure, pass(self) :: compute     !< Compute weights.
    procedure, pass(self) :: description !< Return weights string-description.
    procedure, pass(self) :: destroy     !< Destroy weights.
endtype weights_js

contains
  ! public non TBP
  subroutine create_weights_js_constructor(S, alpha_constructor, beta_constructor, kappa_constructor, constructor, &
                                           face_left, face_right, eps)
  !< Create weights constructor.
  integer(I_P),                                   intent(in)           :: S                 !< Stencils dimension.
  class(alpha_object_constructor),                intent(in)           :: alpha_constructor !< Alpha constructor.
  class(beta_object_constructor),                 intent(in)           :: beta_constructor  !< Beta constructor.
  class(kappa_object_constructor),                intent(in)           :: kappa_constructor !< kappa constructor.
  class(weights_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  logical,                                        intent(in), optional :: face_left         !< Activate left-face interpolations.
  logical,                                        intent(in), optional :: face_right        !< Activate right-face interpolations.
  real(R_P),                                      intent(in), optional :: eps               !< Small epsilon to avoid zero-div.

  allocate(weights_js_constructor :: constructor)
  constructor%S = S
  if (present(face_left)) constructor%face_left = face_left
  if (present(face_right)) constructor%face_right = face_right
  if (present(eps)) constructor%eps = eps
  select type(constructor)
  type is(weights_js_constructor)
    allocate(constructor%alpha_constructor, source=alpha_constructor)
    allocate(constructor%beta_constructor, source=beta_constructor)
    allocate(constructor%kappa_constructor, source=kappa_constructor)
  endselect
  endsubroutine create_weights_js_constructor

  ! deferred public methods
  subroutine create(self, constructor)
  !< Create reconstructor.
  class(weights_js),              intent(inout) :: self        !< Weights.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._R_P
  select type(constructor)
  type is(weights_js_constructor)
    associate(alpha_constructor=>constructor%alpha_constructor, &
              beta_constructor=>constructor%beta_constructor,   &
              kappa_constructor=>constructor%kappa_constructor)

      select type(alpha_constructor)
      type is(alpha_rec_js_constructor)
        allocate(alpha_rec_js :: self%alpha)
        call self%alpha%create(constructor=alpha_constructor)
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
  real(R_P),         intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
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
  class(weights_js), intent(in) :: self   !< Weights.
  character(len=:), allocatable :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'weights_js%description to be implemented, do not use!'
#endif
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
endmodule wenoof_weights_js
