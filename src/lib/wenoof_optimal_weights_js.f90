!< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights.
module wenoof_optimal_weights_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights.
!<
!< @note The provided WENO optimal weights implements the optimal weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P
use wenoof_optimal_weights

implicit none
private
public :: optimal_weights_js
public :: optimal_weights_js_constructor
public :: create_optimal_weights_js_constructor

type, extends(optimal_weights_constructor) :: optimal_weights_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights object constructor.
endtype optimal_weights_js_constructor

type, extends(optimal_weights):: optimal_weights_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal weights object.
  !<
  !< @note The provided WENO optimal weights implements the optimal weights defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute weights.
    procedure, nopass     :: description !< Return weights string-description.
endtype optimal_weights_js

contains
  ! public non TBP
  subroutine create_optimal_weights_js_constructor(S, constructor)
  !< Create optimal weights constructor.
  integer(I_P),                                    intent(in)  :: S           !< Stencils dimension.
  class(optimal_weights_constructor), allocatable, intent(out) :: constructor !< Optimal weights constructor.

  allocate(optimal_weights_js_constructor :: constructor)
  constructor%S = S
  endsubroutine create_optimal_weights_js_constructor

  ! deferred public methods
  pure subroutine compute(self, S)
  !< Compute weights.
  class(optimal_weights_js), intent(inout) :: self !< Optimal weights.
  integer(I_P),              intent(in)    :: S    !< Number of stencils used.

  associate(opt => self%opt)
    select case(S)
      case(2) ! 3rd order
        ! 1 => left interface (i-1/2)
        opt(1, 0) = 2._R_P/3._R_P ! stencil 0
        opt(1, 1) = 1._R_P/3._R_P ! stencil 1
        ! 2 => right interface (i+1/2)
        opt(2, 0) = 1._R_P/3._R_P ! stencil 0
        opt(2, 1) = 2._R_P/3._R_P ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) = 0.3_R_P ! stencil 0
        opt(1, 1) = 0.6_R_P ! stencil 1
        opt(1, 2) = 0.1_R_P ! stencil 2
        ! 2 => right interface (i+1/2)
        opt(2, 0) = 0.1_R_P ! stencil 0
        opt(2, 1) = 0.6_R_P ! stencil 1
        opt(2, 2) = 0.3_R_P ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =  4._R_P/35._R_P ! stencil 0
        opt(1, 1) = 18._R_P/35._R_P ! stencil 1
        opt(1, 2) = 12._R_P/35._R_P ! stencil 2
        opt(1, 3) =  1._R_P/35._R_P ! stencil 3
        ! 2 => right interface (i+1/2)
        opt(2, 0) =  1._R_P/35._R_P ! stencil 0
        opt(2, 1) = 12._R_P/35._R_P ! stencil 1
        opt(2, 2) = 18._R_P/35._R_P ! stencil 2
        opt(2, 3) =  4._R_P/35._R_P ! stencil 3
      case(5) ! 9th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =  5._R_P/126._R_P ! stencil 0
        opt(1, 1) = 20._R_P/63._R_P  ! stencil 1
        opt(1, 2) = 10._R_P/21._R_P  ! stencil 2
        opt(1, 3) = 10._R_P/63._R_P  ! stencil 3
        opt(1, 4) =  1._R_P/126._R_P ! stencil 4
        ! 2 => right interface (i+1/2)
        opt(2, 0) =  1._R_P/126._R_P ! stencil 0
        opt(2, 1) = 10._R_P/63._R_P  ! stencil 1
        opt(2, 2) = 10._R_P/21._R_P  ! stencil 2
        opt(2, 3) = 20._R_P/63._R_P  ! stencil 3
        opt(2, 4) =  5._R_P/126._R_P ! stencil 4
      case(6) ! 11th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =   1._R_P/77._R_P  ! stencil 0
        opt(1, 1) =  25._R_P/154._R_P ! stencil 1
        opt(1, 2) = 100._R_P/231._R_P ! stencil 2
        opt(1, 3) =  25._R_P/77._R_P  ! stencil 3
        opt(1, 4) =   5._R_P/77._R_P  ! stencil 4
        opt(1, 5) =   1._R_P/462._R_P ! stencil 5
        ! 2 => right interface (i+1/2)
        opt(2, 0) =   1._R_P/462._R_P ! stencil 0
        opt(2, 1) =   5._R_P/77._R_P  ! stencil 1
        opt(2, 2) =  25._R_P/77._R_P  ! stencil 2
        opt(2, 3) = 100._R_P/231._R_P ! stencil 3
        opt(2, 4) =  25._R_P/154._R_P ! stencil 4
        opt(2, 5) =   1._R_P/77._R_P  ! stencil 5
      case(7) ! 13th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =   7._R_P/1716._R_P ! stencil 0
        opt(1, 1) =  21._R_P/286._R_P  ! stencil 1
        opt(1, 2) = 175._R_P/572._R_P  ! stencil 2
        opt(1, 3) = 175._R_P/429._R_P  ! stencil 3
        opt(1, 4) = 105._R_P/572._R_P  ! stencil 4
        opt(1, 5) =   7._R_P/286._R_P  ! stencil 5
        opt(1, 6) =   1._R_P/1716._R_P ! stencil 6
        ! 2 => right interface (i+1/2)
        opt(2, 0) =   1._R_P/1716._R_P ! stencil 0
        opt(2, 1) =   7._R_P/286._R_P  ! stencil 1
        opt(2, 2) = 105._R_P/572._R_P  ! stencil 2
        opt(2, 3) = 175._R_P/429._R_P  ! stencil 3
        opt(2, 4) = 175._R_P/572._R_P  ! stencil 4
        opt(2, 5) =  21._R_P/286._R_P  ! stencil 5
        opt(2, 6) =   7._R_P/1716._R_P ! stencil 6
      case(8) ! 15th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =   8._R_P/6435._R_P ! stencil 0
        opt(1, 1) = 196._R_P/6435._R_P ! stencil 1
        opt(1, 2) = 392._R_P/2145._R_P ! stencil 2
        opt(1, 3) = 490._R_P/1287._R_P ! stencil 3
        opt(1, 4) = 392._R_P/1287._R_P ! stencil 4
        opt(1, 5) = 196._R_P/2145._R_P ! stencil 5
        opt(1, 6) =  56._R_P/6435._R_P ! stencil 6
        opt(1, 7) =   1._R_P/6435._R_P ! stencil 7
        ! 2 => right interface (i+1/2)
        opt(2, 0) =   1._R_P/6435._R_P ! stencil 0
        opt(2, 1) =  56._R_P/6435._R_P ! stencil 1
        opt(2, 2) = 196._R_P/2145._R_P ! stencil 2
        opt(2, 3) = 392._R_P/1287._R_P ! stencil 3
        opt(2, 4) = 490._R_P/1287._R_P ! stencil 4
        opt(2, 5) = 392._R_P/2145._R_P ! stencil 5
        opt(2, 6) = 196._R_P/6435._R_P ! stencil 6
        opt(2, 7) =   8._R_P/6435._R_P ! stencil 7
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =    9._R_P/24310._R_P ! stencil 0
        opt(1, 1) =  144._R_P/12155._R_P ! stencil 1
        opt(1, 2) = 1176._R_P/12155._R_P ! stencil 2
        opt(1, 3) = 3528._R_P/12155._R_P ! stencil 3
        opt(1, 4) =  882._R_P/2431._R_P  ! stencil 4
        opt(1, 5) = 2352._R_P/12155._R_P ! stencil 5
        opt(1, 6) =  504._R_P/12155._R_P ! stencil 6
        opt(1, 7) =   36._R_P/12155._R_P ! stencil 7
        opt(1, 8) =    1._R_P/24310._R_P ! stencil 8
        ! 2 => right interface (i+1/2)
        opt(2, 0) =    1._R_P/24310._R_P ! stencil 0
        opt(2, 1) =   36._R_P/12155._R_P ! stencil 1
        opt(2, 2) =  504._R_P/12155._R_P ! stencil 2
        opt(2, 3) = 2352._R_P/12155._R_P ! stencil 3
        opt(2, 4) =  882._R_P/2431._R_P  ! stencil 4
        opt(2, 5) = 3528._R_P/12155._R_P ! stencil 5
        opt(2, 6) = 1176._R_P/12155._R_P ! stencil 6
        opt(2, 7) =  144._R_P/12155._R_P ! stencil 7
        opt(2, 8) =    9._R_P/24310._R_P ! stencil 8
    endselect
  endassociate
  endsubroutine compute

  pure function description() result(string)
  !< Return string-description of weights.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line character.

  string = 'WENO optimal weights,'//nl
  string = string//'  based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and'//nl
  string = string//'  on the work by Gerolymos, Senechal and Vallet "Very-high-order weno schemes", see '// &
           'JCP, 2009, vol. 228, pp. 8481--8524, doi:10.1016/j.jcp.2009.07.039'//nl
  string = string//'    The optimal weights are allocated in a two-dimensional array, in which the first index'//nl
  string = string//'    is the face selected (1 => i-1/2, 2 => i+1/2) and the second index is the number of the stencil '//nl
  string = string//'    (from 0 to S-1)'
  endfunction description
endmodule wenoof_optimal_weights_js
