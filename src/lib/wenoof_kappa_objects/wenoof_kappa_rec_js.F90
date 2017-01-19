!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for reconstruction.
module wenoof_kappa_rec_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for reconstruction.
!<
!< @note The provided WENO kappa implements the linear weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P
use wenoof_kappa_object

implicit none
private
public :: kappa_rec_js
public :: kappa_rec_js_constructor

type, extends(kappa_object_constructor) :: kappa_rec_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal kappa object constructor.
endtype kappa_rec_js_constructor

type, extends(kappa_object):: kappa_rec_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet kappa object.
  !<
  !< @note The provided WENO kappa implements the weights defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  contains
    ! public deferred methods
    procedure, pass(self) :: compute     !< Compute kappa.
    procedure, pass(self) :: description !< Return kappa string-description.
    ! public overridden methods
    procedure, pass(self) :: create !< Create interpolator.
endtype kappa_rec_js

contains
  ! deferred public methods
  pure subroutine compute(self)
  !< Compute kappa.
  class(kappa_rec_js), intent(inout) :: self !< Kappa.

  associate(val => self%values)
    select case(self%S)
      case(2) ! 3rd order
        ! 1 => left interface (i-1/2)
        val(1, 0) = 2._R_P/3._R_P ! stencil 0
        val(1, 1) = 1._R_P/3._R_P ! stencil 1
        ! 2 => right interface (i+1/2)
        val(2, 0) = 1._R_P/3._R_P ! stencil 0
        val(2, 1) = 2._R_P/3._R_P ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        val(1, 0) = 0.3_R_P ! stencil 0
        val(1, 1) = 0.6_R_P ! stencil 1
        val(1, 2) = 0.1_R_P ! stencil 2
        ! 2 => right interface (i+1/2)
        val(2, 0) = 0.1_R_P ! stencil 0
        val(2, 1) = 0.6_R_P ! stencil 1
        val(2, 2) = 0.3_R_P ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =  4._R_P/35._R_P ! stencil 0
        val(1, 1) = 18._R_P/35._R_P ! stencil 1
        val(1, 2) = 12._R_P/35._R_P ! stencil 2
        val(1, 3) =  1._R_P/35._R_P ! stencil 3
        ! 2 => right interface (i+1/2)
        val(2, 0) =  1._R_P/35._R_P ! stencil 0
        val(2, 1) = 12._R_P/35._R_P ! stencil 1
        val(2, 2) = 18._R_P/35._R_P ! stencil 2
        val(2, 3) =  4._R_P/35._R_P ! stencil 3
      case(5) ! 9th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =  5._R_P/126._R_P ! stencil 0
        val(1, 1) = 20._R_P/63._R_P  ! stencil 1
        val(1, 2) = 10._R_P/21._R_P  ! stencil 2
        val(1, 3) = 10._R_P/63._R_P  ! stencil 3
        val(1, 4) =  1._R_P/126._R_P ! stencil 4
        ! 2 => right interface (i+1/2)
        val(2, 0) =  1._R_P/126._R_P ! stencil 0
        val(2, 1) = 10._R_P/63._R_P  ! stencil 1
        val(2, 2) = 10._R_P/21._R_P  ! stencil 2
        val(2, 3) = 20._R_P/63._R_P  ! stencil 3
        val(2, 4) =  5._R_P/126._R_P ! stencil 4
      case(6) ! 11th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   1._R_P/77._R_P  ! stencil 0
        val(1, 1) =  25._R_P/154._R_P ! stencil 1
        val(1, 2) = 100._R_P/231._R_P ! stencil 2
        val(1, 3) =  25._R_P/77._R_P  ! stencil 3
        val(1, 4) =   5._R_P/77._R_P  ! stencil 4
        val(1, 5) =   1._R_P/462._R_P ! stencil 5
        ! 2 => right interface (i+1/2)
        val(2, 0) =   1._R_P/462._R_P ! stencil 0
        val(2, 1) =   5._R_P/77._R_P  ! stencil 1
        val(2, 2) =  25._R_P/77._R_P  ! stencil 2
        val(2, 3) = 100._R_P/231._R_P ! stencil 3
        val(2, 4) =  25._R_P/154._R_P ! stencil 4
        val(2, 5) =   1._R_P/77._R_P  ! stencil 5
      case(7) ! 13th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   7._R_P/1716._R_P ! stencil 0
        val(1, 1) =  21._R_P/286._R_P  ! stencil 1
        val(1, 2) = 175._R_P/572._R_P  ! stencil 2
        val(1, 3) = 175._R_P/429._R_P  ! stencil 3
        val(1, 4) = 105._R_P/572._R_P  ! stencil 4
        val(1, 5) =   7._R_P/286._R_P  ! stencil 5
        val(1, 6) =   1._R_P/1716._R_P ! stencil 6
        ! 2 => right interface (i+1/2)
        val(2, 0) =   1._R_P/1716._R_P ! stencil 0
        val(2, 1) =   7._R_P/286._R_P  ! stencil 1
        val(2, 2) = 105._R_P/572._R_P  ! stencil 2
        val(2, 3) = 175._R_P/429._R_P  ! stencil 3
        val(2, 4) = 175._R_P/572._R_P  ! stencil 4
        val(2, 5) =  21._R_P/286._R_P  ! stencil 5
        val(2, 6) =   7._R_P/1716._R_P ! stencil 6
      case(8) ! 15th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   8._R_P/6435._R_P ! stencil 0
        val(1, 1) = 196._R_P/6435._R_P ! stencil 1
        val(1, 2) = 392._R_P/2145._R_P ! stencil 2
        val(1, 3) = 490._R_P/1287._R_P ! stencil 3
        val(1, 4) = 392._R_P/1287._R_P ! stencil 4
        val(1, 5) = 196._R_P/2145._R_P ! stencil 5
        val(1, 6) =  56._R_P/6435._R_P ! stencil 6
        val(1, 7) =   1._R_P/6435._R_P ! stencil 7
        ! 2 => right interface (i+1/2)
        val(2, 0) =   1._R_P/6435._R_P ! stencil 0
        val(2, 1) =  56._R_P/6435._R_P ! stencil 1
        val(2, 2) = 196._R_P/2145._R_P ! stencil 2
        val(2, 3) = 392._R_P/1287._R_P ! stencil 3
        val(2, 4) = 490._R_P/1287._R_P ! stencil 4
        val(2, 5) = 392._R_P/2145._R_P ! stencil 5
        val(2, 6) = 196._R_P/6435._R_P ! stencil 6
        val(2, 7) =   8._R_P/6435._R_P ! stencil 7
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =    9._R_P/24310._R_P ! stencil 0
        val(1, 1) =  144._R_P/12155._R_P ! stencil 1
        val(1, 2) = 1176._R_P/12155._R_P ! stencil 2
        val(1, 3) = 3528._R_P/12155._R_P ! stencil 3
        val(1, 4) =  882._R_P/2431._R_P  ! stencil 4
        val(1, 5) = 2352._R_P/12155._R_P ! stencil 5
        val(1, 6) =  504._R_P/12155._R_P ! stencil 6
        val(1, 7) =   36._R_P/12155._R_P ! stencil 7
        val(1, 8) =    1._R_P/24310._R_P ! stencil 8
        ! 2 => right interface (i+1/2)
        val(2, 0) =    1._R_P/24310._R_P ! stencil 0
        val(2, 1) =   36._R_P/12155._R_P ! stencil 1
        val(2, 2) =  504._R_P/12155._R_P ! stencil 2
        val(2, 3) = 2352._R_P/12155._R_P ! stencil 3
        val(2, 4) =  882._R_P/2431._R_P  ! stencil 4
        val(2, 5) = 3528._R_P/12155._R_P ! stencil 5
        val(2, 6) = 1176._R_P/12155._R_P ! stencil 6
        val(2, 7) =  144._R_P/12155._R_P ! stencil 7
        val(2, 8) =    9._R_P/24310._R_P ! stencil 8
    endselect
  endassociate
  endsubroutine compute

  pure function description(self) result(string)
  !< Return string-description of kappa.
  class(kappa_rec_js), intent(in) :: self             !< Kappa.
  character(len=:), allocatable   :: string           !< String-description.
  character(len=1), parameter     :: nl=new_line('a') !< New line character.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'kappa_rec_js%description to be implemented, do not use!'
#endif
  endfunction description

  ! overridden methods
  subroutine create(self, constructor)
  !< Create interpolator.
  !<
  !< @note The kappa coefficients are also computed, they being constants.
  class(kappa_rec_js),            intent(inout) :: self        !< Kappa.
  class(base_object_constructor), intent(in)    :: constructor !< Kappa constructor.

  call self%destroy
  call self%kappa_object%create(constructor=constructor)
  call self%compute
  endsubroutine create
endmodule wenoof_kappa_rec_js
