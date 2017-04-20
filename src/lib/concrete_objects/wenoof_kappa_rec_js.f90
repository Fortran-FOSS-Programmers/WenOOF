!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for reconstruction.
module wenoof_kappa_rec_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for reconstruction.
!<
!< @note The provided WENO kappa implements the linear weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P, str
use wenoof_base_object, only : base_object, base_object_constructor
use wenoof_kappa_object, only : kappa_object, kappa_object_constructor

implicit none
private
public :: kappa_rec_js
public :: kappa_rec_js_constructor

type, extends(kappa_object_constructor) :: kappa_rec_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal kappa object constructor.
  contains
    ! public deferred methods
    procedure, pass(lhs) :: constr_assign_constr !< `=` operator.
endtype kappa_rec_js_constructor

type, extends(kappa_object):: kappa_rec_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet kappa object.
  !<
  !< @note The provided WENO kappa implements the weights defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  real(R_P), allocatable :: values(:,:) !< Kappa coefficients values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: create               !< Create kappa.
    procedure, pass(self) :: compute_int          !< Compute kappa (interpolate).
    procedure, pass(self) :: compute_rec          !< Compute kappa (reconstruct).
    procedure, pass(self) :: description          !< Return object string-description.
    procedure, pass(self) :: destroy              !< Destroy kappa.
    procedure, pass(lhs)  :: object_assign_object !< `=` operator.
endtype kappa_rec_js

contains
  ! constructor

  ! deferred public methods
  subroutine constr_assign_constr(lhs, rhs)
  !< `=` operator.
  class(kappa_rec_js_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  endsubroutine constr_assign_constr

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
  call self%compute_rec(values=self%values)
  endsubroutine create

  pure subroutine compute_int(self, stencil, x_target, values)
  !< Compute kappa (interpolate).
  class(kappa_rec_js), intent(in)  :: self               !< Kappa.
  real(R_P),           intent(in)  :: stencil(1-self%S:) !< Stencil used for interpolation, [1-S:S-1].
  real(R_P),           intent(in)  :: x_target           !< Coordinate of the interpolation point.
  real(R_P),           intent(out) :: values(0:)         !< Kappa values.
  ! empty procedure
  endsubroutine compute_int

  pure subroutine compute_rec(self, values)
  !< Compute kappa (reconstruct).
  class(kappa_rec_js), intent(in)  :: self          !< Kappa.
  real(R_P),           intent(out) :: values(1:,0:) !< Kappa values.

  select case(self%S)
    case(2) ! 3rd order
      ! 1 => left interface (i-1/2)
      values(1, 0) = 2._R_P/3._R_P ! stencil 0
      values(1, 1) = 1._R_P/3._R_P ! stencil 1
      ! 2 => right interface (i+1/2)
      values(2, 0) = 1._R_P/3._R_P ! stencil 0
      values(2, 1) = 2._R_P/3._R_P ! stencil 1
    case(3) ! 5th order
      ! 1 => left interface (i-1/2)
      values(1, 0) = 0.3_R_P ! stencil 0
      values(1, 1) = 0.6_R_P ! stencil 1
      values(1, 2) = 0.1_R_P ! stencil 2
      ! 2 => right interface (i+1/2)
      values(2, 0) = 0.1_R_P ! stencil 0
      values(2, 1) = 0.6_R_P ! stencil 1
      values(2, 2) = 0.3_R_P ! stencil 2
    case(4) ! 7th order
      ! 1 => left interface (i-1/2)
      values(1, 0) =  4._R_P/35._R_P ! stencil 0
      values(1, 1) = 18._R_P/35._R_P ! stencil 1
      values(1, 2) = 12._R_P/35._R_P ! stencil 2
      values(1, 3) =  1._R_P/35._R_P ! stencil 3
      ! 2 => right interface (i+1/2)
      values(2, 0) =  1._R_P/35._R_P ! stencil 0
      values(2, 1) = 12._R_P/35._R_P ! stencil 1
      values(2, 2) = 18._R_P/35._R_P ! stencil 2
      values(2, 3) =  4._R_P/35._R_P ! stencil 3
    case(5) ! 9th order
      ! 1 => left interface (i-1/2)
      values(1, 0) =  5._R_P/126._R_P ! stencil 0
      values(1, 1) = 20._R_P/63._R_P  ! stencil 1
      values(1, 2) = 10._R_P/21._R_P  ! stencil 2
      values(1, 3) = 10._R_P/63._R_P  ! stencil 3
      values(1, 4) =  1._R_P/126._R_P ! stencil 4
      ! 2 => right interface (i+1/2)
      values(2, 0) =  1._R_P/126._R_P ! stencil 0
      values(2, 1) = 10._R_P/63._R_P  ! stencil 1
      values(2, 2) = 10._R_P/21._R_P  ! stencil 2
      values(2, 3) = 20._R_P/63._R_P  ! stencil 3
      values(2, 4) =  5._R_P/126._R_P ! stencil 4
    case(6) ! 11th order
      ! 1 => left interface (i-1/2)
      values(1, 0) =   1._R_P/77._R_P  ! stencil 0
      values(1, 1) =  25._R_P/154._R_P ! stencil 1
      values(1, 2) = 100._R_P/231._R_P ! stencil 2
      values(1, 3) =  25._R_P/77._R_P  ! stencil 3
      values(1, 4) =   5._R_P/77._R_P  ! stencil 4
      values(1, 5) =   1._R_P/462._R_P ! stencil 5
      ! 2 => right interface (i+1/2)
      values(2, 0) =   1._R_P/462._R_P ! stencil 0
      values(2, 1) =   5._R_P/77._R_P  ! stencil 1
      values(2, 2) =  25._R_P/77._R_P  ! stencil 2
      values(2, 3) = 100._R_P/231._R_P ! stencil 3
      values(2, 4) =  25._R_P/154._R_P ! stencil 4
      values(2, 5) =   1._R_P/77._R_P  ! stencil 5
    case(7) ! 13th order
      ! 1 => left interface (i-1/2)
      values(1, 0) =   7._R_P/1716._R_P ! stencil 0
      values(1, 1) =  21._R_P/286._R_P  ! stencil 1
      values(1, 2) = 175._R_P/572._R_P  ! stencil 2
      values(1, 3) = 175._R_P/429._R_P  ! stencil 3
      values(1, 4) = 105._R_P/572._R_P  ! stencil 4
      values(1, 5) =   7._R_P/286._R_P  ! stencil 5
      values(1, 6) =   1._R_P/1716._R_P ! stencil 6
      ! 2 => right interface (i+1/2)
      values(2, 0) =   1._R_P/1716._R_P ! stencil 0
      values(2, 1) =   7._R_P/286._R_P  ! stencil 1
      values(2, 2) = 105._R_P/572._R_P  ! stencil 2
      values(2, 3) = 175._R_P/429._R_P  ! stencil 3
      values(2, 4) = 175._R_P/572._R_P  ! stencil 4
      values(2, 5) =  21._R_P/286._R_P  ! stencil 5
      values(2, 6) =   7._R_P/1716._R_P ! stencil 6
    case(8) ! 15th order
      ! 1 => left interface (i-1/2)
      values(1, 0) =   8._R_P/6435._R_P ! stencil 0
      values(1, 1) = 196._R_P/6435._R_P ! stencil 1
      values(1, 2) = 392._R_P/2145._R_P ! stencil 2
      values(1, 3) = 490._R_P/1287._R_P ! stencil 3
      values(1, 4) = 392._R_P/1287._R_P ! stencil 4
      values(1, 5) = 196._R_P/2145._R_P ! stencil 5
      values(1, 6) =  56._R_P/6435._R_P ! stencil 6
      values(1, 7) =   1._R_P/6435._R_P ! stencil 7
      ! 2 => right interface (i+1/2)
      values(2, 0) =   1._R_P/6435._R_P ! stencil 0
      values(2, 1) =  56._R_P/6435._R_P ! stencil 1
      values(2, 2) = 196._R_P/2145._R_P ! stencil 2
      values(2, 3) = 392._R_P/1287._R_P ! stencil 3
      values(2, 4) = 490._R_P/1287._R_P ! stencil 4
      values(2, 5) = 392._R_P/2145._R_P ! stencil 5
      values(2, 6) = 196._R_P/6435._R_P ! stencil 6
      values(2, 7) =   8._R_P/6435._R_P ! stencil 7
    case(9) ! 17th order
      ! 1 => left interface (i-1/2)
      values(1, 0) =    9._R_P/24310._R_P ! stencil 0
      values(1, 1) =  144._R_P/12155._R_P ! stencil 1
      values(1, 2) = 1176._R_P/12155._R_P ! stencil 2
      values(1, 3) = 3528._R_P/12155._R_P ! stencil 3
      values(1, 4) =  882._R_P/2431._R_P  ! stencil 4
      values(1, 5) = 2352._R_P/12155._R_P ! stencil 5
      values(1, 6) =  504._R_P/12155._R_P ! stencil 6
      values(1, 7) =   36._R_P/12155._R_P ! stencil 7
      values(1, 8) =    1._R_P/24310._R_P ! stencil 8
      ! 2 => right interface (i+1/2)
      values(2, 0) =    1._R_P/24310._R_P ! stencil 0
      values(2, 1) =   36._R_P/12155._R_P ! stencil 1
      values(2, 2) =  504._R_P/12155._R_P ! stencil 2
      values(2, 3) = 2352._R_P/12155._R_P ! stencil 3
      values(2, 4) =  882._R_P/2431._R_P  ! stencil 4
      values(2, 5) = 3528._R_P/12155._R_P ! stencil 5
      values(2, 6) = 1176._R_P/12155._R_P ! stencil 6
      values(2, 7) =  144._R_P/12155._R_P ! stencil 7
      values(2, 8) =    9._R_P/24310._R_P ! stencil 8
  endselect
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(kappa_rec_js), intent(in)           :: self             !< Kappa coefficient.
  character(len=*),    intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable             :: string           !< String-description.
  character(len=:), allocatable             :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter               :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu kappa coefficients object for reconstruction:'//NL
  string = string//prefix_//'  - S   = '//trim(str(self%S))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy kappa.
  class(kappa_rec_js), intent(inout) :: self !< Kappa.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy

  subroutine object_assign_object(lhs, rhs)
  !< `=` operator.
  class(kappa_rec_js), intent(inout) :: lhs !< Left hand side.
  class(base_object),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(kappa_rec_js)
     if (allocated(rhs%values)) then
        lhs%values = rhs%values
     else
        if (allocated(lhs%values)) deallocate(lhs%values)
     endif
  endselect
  endsubroutine object_assign_object
endmodule wenoof_kappa_rec_js
