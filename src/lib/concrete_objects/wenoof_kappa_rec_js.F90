!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for reconstruction.
module wenoof_kappa_rec_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for reconstruction.
!<
!< @note The provided WENO kappa implements the linear weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
use wenoof_base_object
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
  real(RPP), allocatable :: values(:,:) !< Kappa coefficients values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: create             !< Create kappa.
    procedure, pass(self) :: compute_kappa_int  !< Compute kappa.
    procedure, pass(self) :: compute_kappa_rec  !< Compute kappa.
    procedure, pass(self) :: description        !< Return kappa string-description.
    procedure, pass(self) :: destroy            !< Destroy kappa.
endtype kappa_rec_js

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create kappa.
  !<
  !< @note The kappa coefficients are also computed, they being constants.
  class(kappa_rec_js),            intent(inout) :: self        !< Kappa.
  class(base_object_constructor), intent(in)    :: constructor !< Kappa constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._RPP
  associate(val => self%values)
    select case(self%S)
      case(2) ! 3rd order
        ! 1 => left interface (i-1/2)
        val(1, 0) = 2._RPP/3._RPP ! stencil 0
        val(1, 1) = 1._RPP/3._RPP ! stencil 1
        ! 2 => right interface (i+1/2)
        val(2, 0) = 1._RPP/3._RPP ! stencil 0
        val(2, 1) = 2._RPP/3._RPP ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        val(1, 0) = 0.3_RPP ! stencil 0
        val(1, 1) = 0.6_RPP ! stencil 1
        val(1, 2) = 0.1_RPP ! stencil 2
        ! 2 => right interface (i+1/2)
        val(2, 0) = 0.1_RPP ! stencil 0
        val(2, 1) = 0.6_RPP ! stencil 1
        val(2, 2) = 0.3_RPP ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =  4._RPP/35._RPP ! stencil 0
        val(1, 1) = 18._RPP/35._RPP ! stencil 1
        val(1, 2) = 12._RPP/35._RPP ! stencil 2
        val(1, 3) =  1._RPP/35._RPP ! stencil 3
        ! 2 => right interface (i+1/2)
        val(2, 0) =  1._RPP/35._RPP ! stencil 0
        val(2, 1) = 12._RPP/35._RPP ! stencil 1
        val(2, 2) = 18._RPP/35._RPP ! stencil 2
        val(2, 3) =  4._RPP/35._RPP ! stencil 3
      case(5) ! 9th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =  5._RPP/126._RPP ! stencil 0
        val(1, 1) = 20._RPP/63._RPP  ! stencil 1
        val(1, 2) = 10._RPP/21._RPP  ! stencil 2
        val(1, 3) = 10._RPP/63._RPP  ! stencil 3
        val(1, 4) =  1._RPP/126._RPP ! stencil 4
        ! 2 => right interface (i+1/2)
        val(2, 0) =  1._RPP/126._RPP ! stencil 0
        val(2, 1) = 10._RPP/63._RPP  ! stencil 1
        val(2, 2) = 10._RPP/21._RPP  ! stencil 2
        val(2, 3) = 20._RPP/63._RPP  ! stencil 3
        val(2, 4) =  5._RPP/126._RPP ! stencil 4
      case(6) ! 11th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   1._RPP/77._RPP  ! stencil 0
        val(1, 1) =  25._RPP/154._RPP ! stencil 1
        val(1, 2) = 100._RPP/231._RPP ! stencil 2
        val(1, 3) =  25._RPP/77._RPP  ! stencil 3
        val(1, 4) =   5._RPP/77._RPP  ! stencil 4
        val(1, 5) =   1._RPP/462._RPP ! stencil 5
        ! 2 => right interface (i+1/2)
        val(2, 0) =   1._RPP/462._RPP ! stencil 0
        val(2, 1) =   5._RPP/77._RPP  ! stencil 1
        val(2, 2) =  25._RPP/77._RPP  ! stencil 2
        val(2, 3) = 100._RPP/231._RPP ! stencil 3
        val(2, 4) =  25._RPP/154._RPP ! stencil 4
        val(2, 5) =   1._RPP/77._RPP  ! stencil 5
      case(7) ! 13th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   7._RPP/1716._RPP ! stencil 0
        val(1, 1) =  21._RPP/286._RPP  ! stencil 1
        val(1, 2) = 175._RPP/572._RPP  ! stencil 2
        val(1, 3) = 175._RPP/429._RPP  ! stencil 3
        val(1, 4) = 105._RPP/572._RPP  ! stencil 4
        val(1, 5) =   7._RPP/286._RPP  ! stencil 5
        val(1, 6) =   1._RPP/1716._RPP ! stencil 6
        ! 2 => right interface (i+1/2)
        val(2, 0) =   1._RPP/1716._RPP ! stencil 0
        val(2, 1) =   7._RPP/286._RPP  ! stencil 1
        val(2, 2) = 105._RPP/572._RPP  ! stencil 2
        val(2, 3) = 175._RPP/429._RPP  ! stencil 3
        val(2, 4) = 175._RPP/572._RPP  ! stencil 4
        val(2, 5) =  21._RPP/286._RPP  ! stencil 5
        val(2, 6) =   7._RPP/1716._RPP ! stencil 6
      case(8) ! 15th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   8._RPP/6435._RPP ! stencil 0
        val(1, 1) = 196._RPP/6435._RPP ! stencil 1
        val(1, 2) = 392._RPP/2145._RPP ! stencil 2
        val(1, 3) = 490._RPP/1287._RPP ! stencil 3
        val(1, 4) = 392._RPP/1287._RPP ! stencil 4
        val(1, 5) = 196._RPP/2145._RPP ! stencil 5
        val(1, 6) =  56._RPP/6435._RPP ! stencil 6
        val(1, 7) =   1._RPP/6435._RPP ! stencil 7
        ! 2 => right interface (i+1/2)
        val(2, 0) =   1._RPP/6435._RPP ! stencil 0
        val(2, 1) =  56._RPP/6435._RPP ! stencil 1
        val(2, 2) = 196._RPP/2145._RPP ! stencil 2
        val(2, 3) = 392._RPP/1287._RPP ! stencil 3
        val(2, 4) = 490._RPP/1287._RPP ! stencil 4
        val(2, 5) = 392._RPP/2145._RPP ! stencil 5
        val(2, 6) = 196._RPP/6435._RPP ! stencil 6
        val(2, 7) =   8._RPP/6435._RPP ! stencil 7
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =    9._RPP/24310._RPP ! stencil 0
        val(1, 1) =  144._RPP/12155._RPP ! stencil 1
        val(1, 2) = 1176._RPP/12155._RPP ! stencil 2
        val(1, 3) = 3528._RPP/12155._RPP ! stencil 3
        val(1, 4) =  882._RPP/2431._RPP  ! stencil 4
        val(1, 5) = 2352._RPP/12155._RPP ! stencil 5
        val(1, 6) =  504._RPP/12155._RPP ! stencil 6
        val(1, 7) =   36._RPP/12155._RPP ! stencil 7
        val(1, 8) =    1._RPP/24310._RPP ! stencil 8
        ! 2 => right interface (i+1/2)
        val(2, 0) =    1._RPP/24310._RPP ! stencil 0
        val(2, 1) =   36._RPP/12155._RPP ! stencil 1
        val(2, 2) =  504._RPP/12155._RPP ! stencil 2
        val(2, 3) = 2352._RPP/12155._RPP ! stencil 3
        val(2, 4) =  882._RPP/2431._RPP  ! stencil 4
        val(2, 5) = 3528._RPP/12155._RPP ! stencil 5
        val(2, 6) = 1176._RPP/12155._RPP ! stencil 6
        val(2, 7) =  144._RPP/12155._RPP ! stencil 7
        val(2, 8) =    9._RPP/24310._RPP ! stencil 8
    endselect
  endassociate
  endsubroutine create

  pure function description(self) result(string)
  !< Return string-description of kappa.
  class(kappa_rec_js), intent(in) :: self   !< Kappa.
  character(len=:), allocatable   :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'kappa_rec_js%description to be implemented, do not use!'
#endif
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy kappa.
  class(kappa_rec_js), intent(inout) :: self !< Kappa.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_kappa_rec_js
