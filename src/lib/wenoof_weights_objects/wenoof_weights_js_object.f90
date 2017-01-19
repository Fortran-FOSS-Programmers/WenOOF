!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
module wenoof_weights_js_object
!< Jiang-Shu and Gerolymos-Senechal-Vallet weights.
!<
!< @note The provided WENO weights implements the weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P
use wenoof_alpha_object
use wenoof_beta_object
use wenoof_kappa_object
use wenoof_weights_object

implicit none
private
public :: weights_js_object
public :: weights_js_object_constructor
public :: create_weights_js_object_constructor

type, extends(weights_object_constructor) :: weights_js_object_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights object constructor.
  class(alpha_object_constructor), allocatable :: alpha_constructor !< Alpha coefficients (non linear weights) constructor.
  class(beta_object_constructor),  allocatable :: beta_constructor  !< Beta coefficients (smoothness indicators) constructor.
  class(kappa_object_constructor), allocatable :: kappa_constructor !< kappa coefficients (optimal, linear weights) constructor.
endtype weights_js_object_constructor

type, extends(weights_object):: weights_js_object
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
    procedure, pass(self) :: compute     !< Compute weights.
    procedure, pass(self) :: description !< Return weights string-description.
endtype weights_js_object

contains
  ! public non TBP
  subroutine create_weights_js_object_constructor(S, constructor)
  !< Create weights constructor.
  integer(I_P),                                   intent(in)  :: S           !< Stencils dimension.
  class(weights_object_constructor), allocatable, intent(out) :: constructor !< Weights constructor.

  allocate(weights_js_object_constructor :: constructor)
  constructor%S = S
  endsubroutine create_weights_js_object_constructor

  ! deferred public methods
  pure subroutine compute(self, stencil)
  !< Compute weights.
  class(weights_js_object), intent(inout) :: self                  !< Weights.
  real(R_P),                intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P)                            :: f, s                  !< Counters.

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
  class(weights_js_object), intent(in) :: self             !< Weights.
  character(len=:), allocatable        :: string           !< String-description.
  character(len=1), parameter          :: nl=new_line('a') !< New line character.

  string = 'WENO optimal weights,'//nl
  string = string//'  based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and'//nl
  string = string//'  on the work by Gerolymos, Senechal and Vallet "Very-high-order weno schemes", see '// &
           'JCP, 2009, vol. 228, pp. 8481--8524, doi:10.1016/j.jcp.2009.07.039'//nl
  string = string//'    The optimal weights are allocated in a two-dimensional array, in which the first index'//nl
  string = string//'    is the face selected (1 => i-1/2, 2 => i+1/2) and the second index is the number of the stencil '//nl
  string = string//'    (from 0 to S-1)'
  endfunction description

  ! overridden methods
  subroutine create(self, constructor)
  !< Create reconstructor.
  class(weights_js_object),       intent(inout) :: self        !< Weights.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.

  call self%destroy
  call self%weights_object%create(constructor=constructor)
  select type(constructor)
  type is(weights_js_object_constructor)
  endselect
  endsubroutine create

endmodule wenoof_optimal_weights_js
