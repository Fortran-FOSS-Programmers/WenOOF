!< Jiang-Shu (Lagrange) interpolations object for derivative reconstruction.
module wenoof_interpolations_rec_js
!< Jiang-Shu (Lagrange) interpolations object for derivative reconstruction.
!<
!< @note The provided interpolations implement the Lagrange interpolations defined in *Efficient Implementation
!< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

#ifdef r16p
use penf, only: I_P, RPP=>R16P, str
#else
use penf, only: I_P, RPP=>R8P, str
#endif
use wenoof_base_object, only : base_object_constructor
use wenoof_interpolations_object, only : interpolations_object, interpolations_object_constructor

implicit none
private
public :: interpolations_rec_js
public :: interpolations_rec_js_constructor

type, extends(interpolations_object_constructor) :: interpolations_rec_js_constructor
  !< Jiang-Shu (Lagrange) interpolations object for derivative reconstruction constructor.
endtype interpolations_rec_js_constructor

type, extends(interpolations_object) :: interpolations_rec_js
  !< Jiang-Shu (Lagrange) interpolations object for derivative reconstruction.
  !<
  !< @note The provided interpolations implement the Lagrange interpolations defined in *Efficient Implementation
  !< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  private
  real(RPP), allocatable :: coef(:,:,:) !< Polynomial coefficients [1:2,0:S-1,0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: create      !< Create interpolations.
    procedure, pass(self) :: compute_int !< Compute interpolations (interpolate).
    procedure, pass(self) :: compute_rec !< Compute interpolations (reconstruct).
    procedure, pass(self) :: description !< Return object string-description.
    procedure, pass(self) :: destroy     !< Destroy interpolations.
endtype interpolations_rec_js

contains
  ! public deferred methods
  subroutine create(self, constructor)
  !< Create interpolations.
  class(interpolations_rec_js),   intent(inout) :: self        !< Interpolations.
  class(base_object_constructor), intent(in)    :: constructor !< Interpolations constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%coef(1:2, 0:self%S - 1, 0:self%S - 1))
  associate(c => self%coef)
    select case(self%S)
      case(2) ! 3rd order
        ! 1 => left interface (i-1/2)
        !  cell  0           ;    cell  1
        c(1,0,0) =  0.5_RPP; c(1,1,0) =  0.5_RPP ! stencil 0
        c(1,0,1) = -0.5_RPP; c(1,1,1) =  1.5_RPP ! stencil 1
        ! 2 => right interface (i+1/2)
        !  cell  0           ;    cell  1
        c(2,0,0) =  1.5_RPP; c(2,1,0) = -0.5_RPP ! stencil 0
        c(2,0,1) =  0.5_RPP; c(2,1,1) =  0.5_RPP ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        !  cell  0                 ;    cell  1                 ;    cell  2
        c(1,0,0) =  1._RPP/3._RPP; c(1,1,0) =  5._RPP/6._RPP; c(1,2,0) = -1._RPP/6._RPP ! stencil 0
        c(1,0,1) = -1._RPP/6._RPP; c(1,1,1) =  5._RPP/6._RPP; c(1,2,1) =  1._RPP/3._RPP ! stencil 1
        c(1,0,2) =  1._RPP/3._RPP; c(1,1,2) = -7._RPP/6._RPP; c(1,2,2) = 11._RPP/6._RPP ! stencil 2
        ! 2 => right interface (i+1/2)
        !  cell  0                 ;    cell  1                 ;    cell  2
        c(2,0,0) = 11._RPP/6._RPP; c(2,1,0) = -7._RPP/6._RPP; c(2,2,0) =  1._RPP/3._RPP ! stencil 0
        c(2,0,1) =  1._RPP/3._RPP; c(2,1,1) =  5._RPP/6._RPP; c(2,2,1) = -1._RPP/6._RPP ! stencil 1
        c(2,0,2) = -1._RPP/6._RPP; c(2,1,2) =  5._RPP/6._RPP; c(2,2,2) =  1._RPP/3._RPP ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        !  cell  0               ;    cell  1               ;    cell  2               ;    cell  3
        c(1,0,0)=  1._RPP/4._RPP ; c(1,1,0)= 13._RPP/12._RPP; c(1,2,0)= -5._RPP/12._RPP; c(1,3,0)=  1._RPP/12._RPP ! stencil 0
        c(1,0,1)= -1._RPP/12._RPP; c(1,1,1)=  7._RPP/12._RPP; c(1,2,1)=  7._RPP/12._RPP; c(1,3,1)= -1._RPP/12._RPP ! stencil 1
        c(1,0,2)=  1._RPP/12._RPP; c(1,1,2)= -5._RPP/12._RPP; c(1,2,2)= 13._RPP/12._RPP; c(1,3,2)=  1._RPP/4._RPP  ! stencil 2
        c(1,0,3)= -1._RPP/4._RPP ; c(1,1,3)= 13._RPP/12._RPP; c(1,2,3)=-23._RPP/12._RPP; c(1,3,3)= 25._RPP/12._RPP ! stencil 3
        ! 2 => right interface (i+1/2)
        !  cell  0               ;    cell  1               ;   cell  2                ;    cell  3
        c(2,0,0)= 25._RPP/12._RPP; c(2,1,0)=-23._RPP/12._RPP; c(2,2,0)= 13._RPP/12._RPP; c(2,3,0)= -1._RPP/4._RPP  ! stencil 0
        c(2,0,1)=  1._RPP/4._RPP ; c(2,1,1)= 13._RPP/12._RPP; c(2,2,1)= -5._RPP/12._RPP; c(2,3,1)=  1._RPP/12._RPP ! stencil 1
        c(2,0,2)= -1._RPP/12._RPP; c(2,1,2)=  7._RPP/12._RPP; c(2,2,2)=  7._RPP/12._RPP; c(2,3,2)= -1._RPP/12._RPP ! stencil 2
        c(2,0,3)=  1._RPP/12._RPP; c(2,1,3)= -5._RPP/12._RPP; c(2,2,3)= 13._RPP/12._RPP; c(2,3,3)=  1._RPP/4._RPP  ! stencil 3
      case(5) ! 9th order
        ! 1 => left interface (i-1/2)
        !  cell  0                ;    cell  1                ;    cell  2                ;    cell  3
        c(1,0,0)=   1._RPP/5._RPP ; c(1,1,0)=  77._RPP/60._RPP; c(1,2,0)= -43._RPP/60._RPP; c(1,3,0)=  17._RPP/60._RPP  ! stencil 0
        c(1,0,1)=  -1._RPP/20._RPP; c(1,1,1)=   9._RPP/20._RPP; c(1,2,1)=  47._RPP/60._RPP; c(1,3,1)= -13._RPP/60._RPP  ! stencil 1
        c(1,0,2)=   1._RPP/30._RPP; c(1,1,2)= -13._RPP/60._RPP; c(1,2,2)=  47._RPP/60._RPP; c(1,3,2)=   9._RPP/20._RPP  ! stencil 2
        c(1,0,3)=  -1._RPP/20._RPP; c(1,1,3)=  17._RPP/60._RPP; c(1,2,3)= -43._RPP/60._RPP; c(1,3,3)=  77._RPP/60._RPP  ! stencil 3
        c(1,0,4)=   1._RPP/5._RPP ; c(1,1,4)= -21._RPP/20._RPP; c(1,2,4)= 137._RPP/60._RPP; c(1,3,4)=-163._RPP/60._RPP  ! stencil 4
        !  cell  4
        c(1,4,0)=  -1._RPP/20._RPP  ! stencil 0
        c(1,4,1)=   1._RPP/30._RPP  ! stencil 1
        c(1,4,2)=  -1._RPP/20._RPP  ! stencil 2
        c(1,4,3)=   1._RPP/5._RPP   ! stencil 3
        c(1,4,4)= 137._RPP/60._RPP  ! stencil 4
        ! 2 => right interface (i+1/2)
        !  cell  0               ;    cell  1               ;   cell  2                ;    cell  3
        c(2,0,0)= 137._RPP/60._RPP; c(2,1,0)=-163._RPP/60._RPP; c(2,2,0)= 137._RPP/60._RPP; c(2,3,0)= -21._RPP/20._RPP  ! stencil 0
        c(2,0,1)=   1._RPP/5._RPP ; c(2,1,1)=  77._RPP/60._RPP; c(2,2,1)= -43._RPP/60._RPP; c(2,3,1)=  17._RPP/60._RPP  ! stencil 1
        c(2,0,2)=  -1._RPP/20._RPP; c(2,1,2)=   9._RPP/20._RPP; c(2,2,2)=  47._RPP/60._RPP; c(2,3,2)= -13._RPP/60._RPP  ! stencil 2
        c(2,0,3)=   1._RPP/30._RPP; c(2,1,3)= -13._RPP/60._RPP; c(2,2,3)=  47._RPP/60._RPP; c(2,3,3)=   9._RPP/20._RPP  ! stencil 3
        c(2,0,4)=  -1._RPP/20._RPP; c(2,1,4)=  17._RPP/60._RPP; c(2,2,4)= -43._RPP/60._RPP; c(2,3,4)=  77._RPP/60._RPP  ! stencil 4
        !  cell  4
        c(2,4,0)=   1._RPP/5._RPP  ! stencil 0
        c(2,4,1)=  -1._RPP/20._RPP ! stencil 1
        c(2,4,2)=   1._RPP/30._RPP ! stencil 2
        c(2,4,3)=  -1._RPP/20._RPP ! stencil 3
        c(2,4,4)=   1._RPP/5._RPP  ! stencil 4
      case(6) ! 11th order
        ! 1 => left interface (i-1/2)
        !  cell  0                ;    cell  1                ;    cell  2                ;    cell  3
        c(1,0,0)=   1._RPP/6._RPP ; c(1,1,0)=  29._RPP/20._RPP; c(1,2,0)= -21._RPP/20._RPP; c(1,3,0)=  37._RPP/60._RPP  ! stencil 0
        c(1,0,1)=  -1._RPP/30._RPP; c(1,1,1)=  11._RPP/30._RPP; c(1,2,1)=  19._RPP/20._RPP; c(1,3,1)= -23._RPP/60._RPP  ! stencil 1
        c(1,0,2)=   1._RPP/60._RPP; c(1,1,2)=  -2._RPP/15._RPP; c(1,2,2)=  37._RPP/60._RPP; c(1,3,2)=  37._RPP/60._RPP  ! stencil 2
        c(1,0,3)=  -1._RPP/60._RPP; c(1,1,3)=   7._RPP/60._RPP; c(1,2,3)= -23._RPP/60._RPP; c(1,3,3)=  19._RPP/20._RPP  ! stencil 3
        c(1,0,4)=   1._RPP/30._RPP; c(1,1,4)= -13._RPP/60._RPP; c(1,2,4)=  37._RPP/60._RPP; c(1,3,4)= -21._RPP/20._RPP  ! stencil 4
        c(1,0,5)=  -1._RPP/6._RPP ; c(1,1,5)=  31._RPP/30._RPP; c(1,2,5)=-163._RPP/60._RPP; c(1,3,5)=  79._RPP/20._RPP  ! stencil 5
        !  cell  4                ;    cell  5
        c(1,4,0)= -13._RPP/60._RPP; c(1,5,0)=   1._RPP/30._RPP  ! stencil 0
        c(1,4,1)=   7._RPP/60._RPP; c(1,5,1)=  -1._RPP/60._RPP  ! stencil 1
        c(1,4,2)=  -2._RPP/15._RPP; c(1,5,2)=   1._RPP/60._RPP  ! stencil 2
        c(1,4,3)=  11._RPP/30._RPP; c(1,5,3)=  -1._RPP/30._RPP  ! stencil 3
        c(1,4,4)=  29._RPP/20._RPP; c(1,5,4)=   1._RPP/6._RPP   ! stencil 4
        c(1,4,5)= -71._RPP/20._RPP; c(1,5,5)=  49._RPP/20._RPP  ! stencil 5
        ! 2 => right interface (i+1/2)
        !  cell  0                ;    cell  1                ;   cell  2                 ;    cell  3
        c(2,0,0)=  49._RPP/20._RPP; c(2,1,0)= -71._RPP/20._RPP; c(2,2,0)=  79._RPP/20._RPP; c(2,3,0)=-163._RPP/60._RPP  ! stencil 0
        c(2,0,1)=   1._RPP/6._RPP ; c(2,1,1)=  29._RPP/20._RPP; c(2,2,1)= -21._RPP/20._RPP; c(2,3,1)=  37._RPP/60._RPP  ! stencil 1
        c(2,0,2)=  -1._RPP/30._RPP; c(2,1,2)=  11._RPP/30._RPP; c(2,2,2)=  19._RPP/20._RPP; c(2,3,2)= -23._RPP/60._RPP  ! stencil 2
        c(2,0,3)=   1._RPP/60._RPP; c(2,1,3)=  -2._RPP/15._RPP; c(2,2,3)=  37._RPP/60._RPP; c(2,3,3)=  37._RPP/60._RPP  ! stencil 3
        c(2,0,4)=  -1._RPP/60._RPP; c(2,1,4)=   7._RPP/60._RPP; c(2,2,4)= -23._RPP/60._RPP; c(2,3,4)=  19._RPP/20._RPP  ! stencil 4
        c(2,0,5)=   1._RPP/30._RPP; c(2,1,5)= -13._RPP/60._RPP; c(2,2,5)=  37._RPP/60._RPP; c(2,3,5)= -21._RPP/20._RPP  ! stencil 5
        !  cell  4                ;    cell  5
        c(2,4,0)=  31._RPP/30._RPP; c(2,5,0)=  -1._RPP/6._RPP   ! stencil 0
        c(2,4,1)= -13._RPP/60._RPP; c(2,5,1)=   1._RPP/30._RPP  ! stencil 1
        c(2,4,2)=   7._RPP/60._RPP; c(2,5,2)=  -1._RPP/60._RPP  ! stencil 2
        c(2,4,3)=  -2._RPP/15._RPP; c(2,5,3)=   1._RPP/60._RPP  ! stencil 3
        c(2,4,4)=  11._RPP/30._RPP; c(2,5,4)=  -1._RPP/30._RPP  ! stencil 4
        c(2,4,5)=  29._RPP/20._RPP; c(2,5,5)=   1._RPP/6._RPP   ! stencil 5
      case(7) ! 13th order
        !! 1 => left interface (i-1/2)
        !  cell  0                  ;    cell  1                  ;    cell  2
        !c(1,0,0)=    1._RPP/7._RPP  ; c(1,1,0)=  223._RPP/140._RPP; c(1,2,0)= -197._RPP/140._RPP  ! stencil 0
        !c(1,0,1)=   -1._RPP/42._RPP ; c(1,1,1)=   13._RPP/42._RPP ; c(1,2,1)=  153._RPP/140._RPP  ! stencil 1
        !c(1,0,2)=    1._RPP/105._RPP; c(1,1,2)=  -19._RPP/210._RPP; c(1,2,2)=  107._RPP/210._RPP  ! stencil 2
        !c(1,0,3)=   -1._RPP/140._RPP; c(1,1,3)=    5._RPP/84._RPP ; c(1,2,3)= -101._RPP/420._RPP  ! stencil 3
        !c(1,0,4)=    1._RPP/105._RPP; c(1,1,4)=  -31._RPP/420._RPP; c(1,2,4)=  109._RPP/420._RPP  ! stencil 4
        !c(1,0,5)=   -1._RPP/42._RPP ; c(1,1,5)=  -37._RPP/210._RPP; c(1,2,5)= -241._RPP/420._RPP  ! stencil 5
        !c(1,0,6)=   -1._RPP/7._RPP  ; c(1,1,6)=  -43._RPP/42._RPP ; c(1,2,6)=  667._RPP/210._RPP  ! stencil 6
        !!  cell  3                       cell  4                  ;    cell  5
        !c(1,3,0)=  153._RPP/140._RPP; c(1,4,0)= -241._RPP/420._RPP; c(1,5,0)=   37._RPP/210._RPP  ! stencil 0
        !c(1,3,1)= -241._RPP/420._RPP; c(1,4,1)=  109._RPP/420._RPP; c(1,5,1)=  -31._RPP/420._RPP  ! stencil 1
        !c(1,3,2)=  319._RPP/420._RPP; c(1,4,2)= -101._RPP/420._RPP; c(1,5,2)=    5._RPP/84._RPP   ! stencil 2
        !c(1,3,3)=  319._RPP/420._RPP; c(1,4,3)=  107._RPP/210._RPP; c(1,5,3)=  -19._RPP/210._RPP  ! stencil 3
        !c(1,3,4)= -241._RPP/420._RPP; c(1,4,4)=  153._RPP/140._RPP; c(1,5,4)=   13._RPP/42._RPP   ! stencil 4
        !c(1,3,5)=  153._RPP/140._RPP; c(1,4,5)= -197._RPP/140._RPP; c(1,5,5)=  223._RPP/140._RPP  ! stencil 5
        !c(1,3,6)=-2341._RPP/420._RPP; c(1,4,6)=  853._RPP/140._RPP; c(1,5,6)= -617._RPP/140._RPP  ! stencil 6
        !!  cell  6
        !c(1,6,0)=   -1._RPP/42._RPP   ! stencil 0
        !c(1,6,1)=    1._RPP/105._RPP  ! stencil 1
        !c(1,6,2)=   -1._RPP/140._RPP  ! stencil 2
        !c(1,6,3)=    1._RPP/105._RPP  ! stencil 3
        !c(1,6,4)=   -1._RPP/42._RPP   ! stencil 4
        !c(1,6,5)=    1._RPP/7._RPP    ! stencil 5
        !c(1,6,6)=  363._RPP/140._RPP  ! stencil 6
        !! 2 => right interface (i+1/2)
        !!  cell  0                  ;    cell  1                  ;    cell  2
        !c(2,0,0)=  363._RPP/140._RPP; c(2,1,0)= -617._RPP/140._RPP; c(2,2,0)=  853._RPP/140._RPP  ! stencil 0
        !c(2,0,1)=    1._RPP/7._RPP  ; c(2,1,1)=  223._RPP/140._RPP; c(2,2,1)= -197._RPP/140._RPP  ! stencil 1
        !c(2,0,2)=   -1._RPP/42._RPP ; c(2,1,2)=   13._RPP/42._RPP ; c(2,2,2)=  153._RPP/140._RPP  ! stencil 2
        !c(2,0,3)=    1._RPP/105._RPP; c(2,1,3)=  -19._RPP/210._RPP; c(2,2,3)=  107._RPP/210._RPP  ! stencil 3
        !c(2,0,4)=   -1._RPP/140._RPP; c(2,1,4)=    5._RPP/84._RPP ; c(2,2,4)= -101._RPP/420._RPP  ! stencil 4
        !c(2,0,5)=    1._RPP/105._RPP; c(2,1,5)=  -31._RPP/420._RPP; c(2,2,5)=  109._RPP/420._RPP  ! stencil 5
        !c(2,0,6)=   -1._RPP/42._RPP ; c(2,1,6)=   37._RPP/210._RPP; c(2,2,6)= -241._RPP/420._RPP  ! stencil 6
        !!  cell  3                  ;    cell  4                  ;    cell  5
        !c(2,3,0)=-2341._RPP/420._RPP; c(2,4,0)=  667._RPP/210._RPP; c(2,5,0)=  -43._RPP/42._RPP   ! stencil 0
        !c(2,3,1)=  153._RPP/140._RPP; c(2,4,1)= -241._RPP/420._RPP; c(2,5,1)=   37._RPP/210._RPP  ! stencil 1
        !c(2,3,2)= -241._RPP/420._RPP; c(2,4,2)=  109._RPP/420._RPP; c(2,5,2)=  -31._RPP/420._RPP  ! stencil 2
        !c(2,3,3)=  319._RPP/420._RPP; c(2,4,3)= -101._RPP/420._RPP; c(2,5,3)=    5._RPP/84._RPP   ! stencil 3
        !c(2,3,4)=  319._RPP/420._RPP; c(2,4,4)=  107._RPP/210._RPP; c(2,5,4)=  -19._RPP/210._RPP  ! stencil 4
        !c(2,3,5)= -241._RPP/420._RPP; c(2,4,5)=  153._RPP/140._RPP; c(2,5,5)=   13._RPP/42._RPP   ! stencil 5
        !c(2,3,6)=  153._RPP/140._RPP; c(2,4,6)= -197._RPP/140._RPP; c(2,5,6)=  223._RPP/140._RPP  ! stencil 6
        !!  cell  6
        !c(2,6,0)=    1._RPP/7._RPP    ! stencil 0
        !c(2,6,1)=   -1._RPP/42._RPP   ! stencil 1
        !c(2,6,2)=    1._RPP/105._RPP  ! stencil 2
        !c(2,6,3)=   -1._RPP/140._RPP  ! stencil 3
        !c(2,6,4)=    1._RPP/105._RPP  ! stencil 4
        !c(2,6,5)=   -1._RPP/42._RPP   ! stencil 5
        !c(2,6,6)=    1._RPP/7._RPP    ! stencil 6
        ! 1 => left interface (i-1/2)
        c(1,0,6)=   60._RPP/ 420._RPP; c(1,0,5)=- 10._RPP/ 420._RPP; c(1,0,4)=   4._RPP/ 420._RPP; c(1,0,3)=-  3._RPP/ 420._RPP
        c(1,1,6)=- 430._RPP/ 420._RPP; c(1,1,5)=  74._RPP/ 420._RPP; c(1,1,4)=- 31._RPP/ 420._RPP; c(1,1,3)=  25._RPP/ 420._RPP
        c(1,2,6)= 1334._RPP/ 420._RPP; c(1,2,5)=-241._RPP/ 420._RPP; c(1,2,4)= 109._RPP/ 420._RPP; c(1,2,3)=-101._RPP/ 420._RPP
        c(1,3,6)=-2341._RPP/ 420._RPP; c(1,3,5)= 459._RPP/ 420._RPP; c(1,3,4)=-241._RPP/ 420._RPP; c(1,3,3)= 319._RPP/ 420._RPP
        c(1,4,6)= 2559._RPP/ 420._RPP; c(1,4,5)=-591._RPP/ 420._RPP; c(1,4,4)= 459._RPP/ 420._RPP; c(1,4,3)= 214._RPP/ 420._RPP
        c(1,5,6)=-1851._RPP/ 420._RPP; c(1,5,5)= 669._RPP/ 420._RPP; c(1,5,4)= 130._RPP/ 420._RPP; c(1,5,3)=- 38._RPP/ 420._RPP
        c(1,6,6)= 1089._RPP/ 420._RPP; c(1,6,5)=  60._RPP/ 420._RPP; c(1,6,4)=- 10._RPP/ 420._RPP; c(1,6,3)=   4._RPP/ 420._RPP

        c(1,0,2)=    4._RPP/ 420._RPP; c(1,0,1)=- 10._RPP/ 420._RPP; c(1,0,0)=  60._RPP/ 420._RPP
        c(1,1,2)=-  38._RPP/ 420._RPP; c(1,1,1)= 130._RPP/ 420._RPP; c(1,1,0)= 669._RPP/ 420._RPP
        c(1,2,2)=  214._RPP/ 420._RPP; c(1,2,1)= 459._RPP/ 420._RPP; c(1,2,0)=-591._RPP/ 420._RPP
        c(1,3,2)=  319._RPP/ 420._RPP; c(1,3,1)=-241._RPP/ 420._RPP; c(1,3,0)= 459._RPP/ 420._RPP
        c(1,4,2)=- 101._RPP/ 420._RPP; c(1,4,1)= 109._RPP/ 420._RPP; c(1,4,0)=-241._RPP/ 420._RPP
        c(1,5,2)=   25._RPP/ 420._RPP; c(1,5,1)=- 31._RPP/ 420._RPP; c(1,5,0)=  74._RPP/ 420._RPP
        c(1,6,2)=-   3._RPP/ 420._RPP; c(1,6,1)=   4._RPP/ 420._RPP; c(1,6,0)=- 10._RPP/ 420._RPP
        ! 2 => right interface (i+1/2)
        c(2,6,0)=   60._RPP/ 420._RPP; c(2,6,1)=- 10._RPP/ 420._RPP; c(2,6,2)=   4._RPP/ 420._RPP; c(2,6,3)=-  3._RPP/ 420._RPP
        c(2,5,0)=- 430._RPP/ 420._RPP; c(2,5,1)=  74._RPP/ 420._RPP; c(2,5,2)=- 31._RPP/ 420._RPP; c(2,5,3)=  25._RPP/ 420._RPP
        c(2,4,0)= 1334._RPP/ 420._RPP; c(2,4,1)=-241._RPP/ 420._RPP; c(2,4,2)= 109._RPP/ 420._RPP; c(2,4,3)=-101._RPP/ 420._RPP
        c(2,3,0)=-2341._RPP/ 420._RPP; c(2,3,1)= 459._RPP/ 420._RPP; c(2,3,2)=-241._RPP/ 420._RPP; c(2,3,3)= 319._RPP/ 420._RPP
        c(2,2,0)= 2559._RPP/ 420._RPP; c(2,2,1)=-591._RPP/ 420._RPP; c(2,2,2)= 459._RPP/ 420._RPP; c(2,2,3)= 214._RPP/ 420._RPP
        c(2,1,0)=-1851._RPP/ 420._RPP; c(2,1,1)= 669._RPP/ 420._RPP; c(2,1,2)= 130._RPP/ 420._RPP; c(2,1,3)=- 38._RPP/ 420._RPP
        c(2,0,0)= 1089._RPP/ 420._RPP; c(2,0,1)=  60._RPP/ 420._RPP; c(2,0,2)=- 10._RPP/ 420._RPP; c(2,0,3)=   4._RPP/ 420._RPP

        c(2,6,4)=    4._RPP/ 420._RPP; c(2,6,5)=- 10._RPP/ 420._RPP; c(2,6,6)=  60._RPP/ 420._RPP
        c(2,5,4)=-  38._RPP/ 420._RPP; c(2,5,5)= 130._RPP/ 420._RPP; c(2,5,6)= 669._RPP/ 420._RPP
        c(2,4,4)=  214._RPP/ 420._RPP; c(2,4,5)= 459._RPP/ 420._RPP; c(2,4,6)=-591._RPP/ 420._RPP
        c(2,3,4)=  319._RPP/ 420._RPP; c(2,3,5)=-241._RPP/ 420._RPP; c(2,3,6)= 459._RPP/ 420._RPP
        c(2,2,4)=- 101._RPP/ 420._RPP; c(2,2,5)= 109._RPP/ 420._RPP; c(2,2,6)=-241._RPP/ 420._RPP
        c(2,1,4)=   25._RPP/ 420._RPP; c(2,1,5)=- 31._RPP/ 420._RPP; c(2,1,6)=  74._RPP/ 420._RPP
        c(2,0,4)=-   3._RPP/ 420._RPP; c(2,0,5)=   4._RPP/ 420._RPP; c(2,0,6)=- 10._RPP/ 420._RPP
      case(8) ! 15th order
        !! 1 => left interface (i-1/2)
        !!  cell  0                  ;    cell  1                  ;    cell  2
        !c(1,0,0)=    7._RPP/56._RPP ; c(1,1,0)= 1443._RPP/840._RPP; c(1,2,0)=-1497._RPP/840._RPP  ! stencil 0
        !c(1,0,1)=   -1._RPP/56._RPP ; c(1,1,1)=   15._RPP/56._RPP ; c(1,2,1)= 1023._RPP/840._RPP  ! stencil 1
        !c(1,0,2)=    1._RPP/168._RPP; c(1,1,2)=  -11._RPP/168._RPP; c(1,2,2)=   73._RPP/168._RPP  ! stencil 2
        !c(1,0,3)=   -3._RPP/840._RPP; c(1,1,3)=   29._RPP/840._RPP; c(1,2,3)= -139._RPP/840._RPP  ! stencil 3
        !c(1,0,4)=    3._RPP/840._RPP; c(1,1,4)=  -27._RPP/840._RPP; c(1,2,4)=  113._RPP/840._RPP  ! stencil 4
        !c(1,0,5)=   -1._RPP/168._RPP; c(1,1,5)=   43._RPP/840._RPP; c(1,2,5)= -167._RPP/840._RPP  ! stencil 5
        !c(1,0,6)=    1._RPP/56._RPP ; c(1,1,6)=  -25._RPP/168._RPP; c(1,2,6)=  463._RPP/840._RPP  ! stencil 6
        !c(1,0,7)=   -7._RPP/56._RPP ; c(1,1,7)=   57._RPP/56._RPP ; c(1,2,7)= -613._RPP/168._RPP  ! stencil 7
        !!  cell  3                  ;    cell  4                  ;    cell  5
        !c(1,3,0)= 1443._RPP/840._RPP; c(1,4,0)=-1007._RPP/840._RPP; c(1,5,0)=  463._RPP/840._RPP  ! stencil 0
        !c(1,3,1)= -657._RPP/840._RPP; c(1,4,1)=  393._RPP/840._RPP; c(1,5,1)= -167._RPP/840._RPP  ! stencil 1
        !c(1,3,2)=  743._RPP/840._RPP; c(1,4,2)= -307._RPP/840._RPP; c(1,5,2)=  113._RPP/840._RPP  ! stencil 2
        !c(1,3,3)=  533._RPP/840._RPP; c(1,4,3)=  533._RPP/840._RPP; c(1,5,3)= -139._RPP/840._RPP  ! stencil 3
        !c(1,3,4)= -307._RPP/840._RPP; c(1,4,4)=  743._RPP/840._RPP; c(1,5,4)=   73._RPP/168._RPP  ! stencil 4
        !c(1,3,5)=  393._RPP/840._RPP; c(1,4,5)= -657._RPP/840._RPP; c(1,5,5)= 1023._RPP/840._RPP  ! stencil 5
        !c(1,3,6)=-1007._RPP/840._RPP; c(1,4,6)= 1443._RPP/840._RPP; c(1,5,6)=-1497._RPP/840._RPP  ! stencil 6
        !c(1,3,7)= 6343._RPP/840._RPP; c(1,4,7)=-8357._RPP/840._RPP; c(1,5,7)= 7323._RPP/840._RPP  ! stencil 7
        !!  cell  6                  ;    cell  7
        !c(1,6,0)=  -25._RPP/168._RPP; c(1,7,0)=    1._RPP/56._RPP   ! stencil 0
        !c(1,6,1)=   43._RPP/840._RPP; c(1,7,1)=   -1._RPP/168._RPP  ! stencil 1
        !c(1,6,2)=  -27._RPP/840._RPP; c(1,7,2)=    3._RPP/840._RPP  ! stencil 2
        !c(1,6,3)=   29._RPP/840._RPP; c(1,7,3)=   -3._RPP/840._RPP  ! stencil 3
        !c(1,6,4)=  -11._RPP/168._RPP; c(1,7,4)=    1._RPP/168._RPP  ! stencil 4
        !c(1,6,5)=   15._RPP/56._RPP ; c(1,7,5)=   -1._RPP/56._RPP   ! stencil 5
        !c(1,6,6)= 1443._RPP/840._RPP; c(1,7,6)=    7._RPP/56._RPP   ! stencil 6
        !c(1,6,7)=-4437._RPP/840._RPP; c(1,7,7)= 2283._RPP/840._RPP  ! stencil 7
        !! 2 => right interface (i+1/2)
        !!  cell  0                  ;    cell  1                  ;    cell  2
        !c(2,0,0)= 2283._RPP/840._RPP; c(2,1,0)=-4437._RPP/840._RPP; c(2,2,0)= 7323._RPP/840._RPP  ! stencil 0
        !c(2,0,1)=    7._RPP/56._RPP ; c(2,1,1)= 1443._RPP/840._RPP; c(2,2,1)=-1497._RPP/840._RPP  ! stencil 1
        !c(2,0,2)=   -1._RPP/56._RPP ; c(2,1,2)=   15._RPP/56._RPP ; c(2,2,2)= 1023._RPP/840._RPP  ! stencil 2
        !c(2,0,3)=    1._RPP/168._RPP; c(2,1,3)=  -11._RPP/168._RPP; c(2,2,3)=   73._RPP/168._RPP  ! stencil 3
        !c(2,0,4)=   -3._RPP/840._RPP; c(2,1,4)=   29._RPP/840._RPP; c(2,2,4)= -139._RPP/840._RPP  ! stencil 4
        !c(2,0,5)=    3._RPP/840._RPP; c(2,1,5)=  -27._RPP/840._RPP; c(2,2,5)=  113._RPP/840._RPP  ! stencil 5
        !c(2,0,6)=   -1._RPP/168._RPP; c(2,1,6)=   43._RPP/840._RPP; c(2,2,6)= -167._RPP/840._RPP  ! stencil 6
        !c(2,0,7)=    1._RPP/56._RPP ; c(2,1,7)=  -25._RPP/168._RPP; c(2,2,7)=  463._RPP/840._RPP  ! stencil 7
        !!  cell  3                  ;    cell  4                  ;    cell  5
        !c(2,3,0)=-8357._RPP/840._RPP; c(2,4,0)= 6343._RPP/840._RPP; c(2,5,0)= -613._RPP/168._RPP  ! stencil 0
        !c(2,3,1)= 1443._RPP/840._RPP; c(2,4,1)=-1007._RPP/840._RPP; c(2,5,1)=  463._RPP/840._RPP  ! stencil 1
        !c(2,3,2)= -657._RPP/840._RPP; c(2,4,2)=  393._RPP/840._RPP; c(2,5,2)= -167._RPP/840._RPP  ! stencil 2
        !c(2,3,3)=  743._RPP/840._RPP; c(2,4,3)= -307._RPP/840._RPP; c(2,5,3)=  113._RPP/840._RPP  ! stencil 3
        !c(2,3,4)=  533._RPP/840._RPP; c(2,4,4)=  533._RPP/840._RPP; c(2,5,4)= -139._RPP/840._RPP  ! stencil 4
        !c(2,3,5)= -307._RPP/840._RPP; c(2,4,5)=  743._RPP/840._RPP; c(2,5,5)=   73._RPP/168._RPP  ! stencil 5
        !c(2,3,6)=  393._RPP/840._RPP; c(2,4,6)= -657._RPP/840._RPP; c(2,5,6)= 1023._RPP/840._RPP  ! stencil 6
        !c(2,3,7)=-1007._RPP/840._RPP; c(2,4,7)= 1443._RPP/840._RPP; c(2,5,7)=-1497._RPP/840._RPP  ! stencil 7
        !!  cell  6                  ;    cell  7
        !c(2,6,0)=   57._RPP/56._RPP ; c(2,7,0)=   -7._RPP/56._RPP   ! stencil 0
        !c(2,6,1)=  -25._RPP/168._RPP; c(2,7,1)=    1._RPP/56._RPP   ! stencil 1
        !c(2,6,2)=   43._RPP/840._RPP; c(2,7,2)=   -1._RPP/168._RPP  ! stencil 2
        !c(2,6,3)=  -27._RPP/840._RPP; c(2,7,3)=    3._RPP/840._RPP  ! stencil 3
        !c(2,6,4)=   29._RPP/840._RPP; c(2,7,4)=   -3._RPP/840._RPP  ! stencil 4
        !c(2,6,5)=  -11._RPP/168._RPP; c(2,7,5)=    1._RPP/168._RPP  ! stencil 5
        !c(2,6,6)=   15._RPP/56._RPP ; c(2,7,6)=   -1._RPP/56._RPP   ! stencil 6
        !c(2,6,7)= 1443._RPP/840._RPP; c(2,7,7)=    7._RPP/56._RPP   ! stencil 7
        ! 1 => left interface (i-1/2)
        c(1,0,7)=- 105._RPP/ 840._RPP; c(1,0,6)=   15._RPP/ 840._RPP; c(1,0,5)=-   5._RPP/ 840._RPP; c(1,0,4)=    3._RPP/ 840._RPP
        c(1,1,7)=  855._RPP/ 840._RPP; c(1,1,6)=- 125._RPP/ 840._RPP; c(1,1,5)=   43._RPP/ 840._RPP; c(1,1,4)=-  27._RPP/ 840._RPP
        c(1,2,7)=-3065._RPP/ 840._RPP; c(1,2,6)=  463._RPP/ 840._RPP; c(1,2,5)=- 167._RPP/ 840._RPP; c(1,2,4)=  113._RPP/ 840._RPP
        c(1,3,7)= 6343._RPP/ 840._RPP; c(1,3,6)=-1007._RPP/ 840._RPP; c(1,3,5)=  393._RPP/ 840._RPP; c(1,3,4)=- 307._RPP/ 840._RPP
        c(1,4,7)=-8357._RPP/ 840._RPP; c(1,4,6)= 1443._RPP/ 840._RPP; c(1,4,5)=- 657._RPP/ 840._RPP; c(1,4,4)=  743._RPP/ 840._RPP
        c(1,5,7)= 7323._RPP/ 840._RPP; c(1,5,6)=-1497._RPP/ 840._RPP; c(1,5,5)= 1023._RPP/ 840._RPP; c(1,5,4)=  365._RPP/ 840._RPP
        c(1,6,7)=-4437._RPP/ 840._RPP; c(1,6,6)= 1443._RPP/ 840._RPP; c(1,6,5)=  225._RPP/ 840._RPP; c(1,6,4)=-  55._RPP/ 840._RPP
        c(1,7,7)= 2283._RPP/ 840._RPP; c(1,7,6)=  105._RPP/ 840._RPP; c(1,7,5)=-  15._RPP/ 840._RPP; c(1,7,4)=    5._RPP/ 840._RPP

        c(1,0,3)=-   3._RPP/ 840._RPP; c(1,0,2)=    5._RPP/ 840._RPP; c(1,0,1)=-  15._RPP/ 840._RPP; c(1,0,0)=  105._RPP/ 840._RPP
        c(1,1,3)=   29._RPP/ 840._RPP; c(1,1,2)=-  55._RPP/ 840._RPP; c(1,1,1)=  225._RPP/ 840._RPP; c(1,1,0)= 1443._RPP/ 840._RPP
        c(1,2,3)=- 139._RPP/ 840._RPP; c(1,2,2)=  365._RPP/ 840._RPP; c(1,2,1)= 1023._RPP/ 840._RPP; c(1,2,0)=-1497._RPP/ 840._RPP
        c(1,3,3)=  533._RPP/ 840._RPP; c(1,3,2)=  743._RPP/ 840._RPP; c(1,3,1)=- 657._RPP/ 840._RPP; c(1,3,0)= 1443._RPP/ 840._RPP
        c(1,4,3)=  533._RPP/ 840._RPP; c(1,4,2)=- 307._RPP/ 840._RPP; c(1,4,1)=  393._RPP/ 840._RPP; c(1,4,0)=-1007._RPP/ 840._RPP
        c(1,5,3)=- 139._RPP/ 840._RPP; c(1,5,2)=  113._RPP/ 840._RPP; c(1,5,1)=- 167._RPP/ 840._RPP; c(1,5,0)=  463._RPP/ 840._RPP
        c(1,6,3)=   29._RPP/ 840._RPP; c(1,6,2)=-  27._RPP/ 840._RPP; c(1,6,1)=   43._RPP/ 840._RPP; c(1,6,0)=- 125._RPP/ 840._RPP
        c(1,7,3)=-   3._RPP/ 840._RPP; c(1,7,2)=    3._RPP/ 840._RPP; c(1,7,1)=-   5._RPP/ 840._RPP; c(1,7,0)=   15._RPP/ 840._RPP
        ! 2 => right interface (i+1/2)
        c(2,7,0)=- 105._RPP/ 840._RPP; c(2,7,1)=   15._RPP/ 840._RPP; c(2,7,2)=-   5._RPP/ 840._RPP; c(2,7,3)=    3._RPP/ 840._RPP
        c(2,6,0)=  855._RPP/ 840._RPP; c(2,6,1)=- 125._RPP/ 840._RPP; c(2,6,2)=   43._RPP/ 840._RPP; c(2,6,3)=-  27._RPP/ 840._RPP
        c(2,5,0)=-3065._RPP/ 840._RPP; c(2,5,1)=  463._RPP/ 840._RPP; c(2,5,2)=- 167._RPP/ 840._RPP; c(2,5,3)=  113._RPP/ 840._RPP
        c(2,4,0)= 6343._RPP/ 840._RPP; c(2,4,1)=-1007._RPP/ 840._RPP; c(2,4,2)=  393._RPP/ 840._RPP; c(2,4,3)=- 307._RPP/ 840._RPP
        c(2,3,0)=-8357._RPP/ 840._RPP; c(2,3,1)= 1443._RPP/ 840._RPP; c(2,3,2)=- 657._RPP/ 840._RPP; c(2,3,3)=  743._RPP/ 840._RPP
        c(2,2,0)= 7323._RPP/ 840._RPP; c(2,2,1)=-1497._RPP/ 840._RPP; c(2,2,2)= 1023._RPP/ 840._RPP; c(2,2,3)=  365._RPP/ 840._RPP
        c(2,1,0)=-4437._RPP/ 840._RPP; c(2,1,1)= 1443._RPP/ 840._RPP; c(2,1,2)=  225._RPP/ 840._RPP; c(2,1,3)=-  55._RPP/ 840._RPP
        c(2,0,0)= 2283._RPP/ 840._RPP; c(2,0,1)=  105._RPP/ 840._RPP; c(2,0,2)=-  15._RPP/ 840._RPP; c(2,0,3)=    5._RPP/ 840._RPP

        c(2,7,4)=-   3._RPP/ 840._RPP; c(2,7,5)=    5._RPP/ 840._RPP; c(2,7,6)=-  15._RPP/ 840._RPP; c(2,7,7)=  105._RPP/ 840._RPP
        c(2,6,4)=   29._RPP/ 840._RPP; c(2,6,5)=-  55._RPP/ 840._RPP; c(2,6,6)=  225._RPP/ 840._RPP; c(2,6,7)= 1443._RPP/ 840._RPP
        c(2,5,4)=- 139._RPP/ 840._RPP; c(2,5,5)=  365._RPP/ 840._RPP; c(2,5,6)= 1023._RPP/ 840._RPP; c(2,5,7)=-1497._RPP/ 840._RPP
        c(2,4,4)=  533._RPP/ 840._RPP; c(2,4,5)=  743._RPP/ 840._RPP; c(2,4,6)=- 657._RPP/ 840._RPP; c(2,4,7)= 1443._RPP/ 840._RPP
        c(2,3,4)=  533._RPP/ 840._RPP; c(2,3,5)=- 307._RPP/ 840._RPP; c(2,3,6)=  393._RPP/ 840._RPP; c(2,3,7)=-1007._RPP/ 840._RPP
        c(2,2,4)=- 139._RPP/ 840._RPP; c(2,2,5)=  113._RPP/ 840._RPP; c(2,2,6)=- 167._RPP/ 840._RPP; c(2,2,7)=  463._RPP/ 840._RPP
        c(2,1,4)=   29._RPP/ 840._RPP; c(2,1,5)=-  27._RPP/ 840._RPP; c(2,1,6)=   43._RPP/ 840._RPP; c(2,1,7)=- 125._RPP/ 840._RPP
        c(2,0,4)=-   3._RPP/ 840._RPP; c(2,0,5)=    3._RPP/ 840._RPP; c(2,0,6)=-   5._RPP/ 840._RPP; c(2,0,7)=   15._RPP/ 840._RPP
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        !  cell  0                    ;     cell  1                   ;     cell  2
        !c(1,0,0)=     1._RPP/9._RPP   ; c(1,1,0)=  4609._RPP/2520._RPP; c(1,2,0)= -5471._RPP/2520._RPP  ! stencil 0
        !c(1,0,1)=    -7._RPP/504._RPP ; c(1,1,1)=   119._RPP/508._RPP ; c(1,2,1)=  3349._RPP/2520._RPP  ! stencil 1
        !c(1,0,2)=     1._RPP/252._RPP ; c(1,1,2)=   -25._RPP/508._RPP ; c(1,2,2)=   191._RPP/504._RPP   ! stencil 2
        !c(1,0,3)=    -1._RPP/504._RPP ; c(1,1,3)=    11._RPP/508._RPP ; c(1,2,3)=   -61._RPP/504._RPP   ! stencil 3
        !c(1,0,4)=     1._RPP/630._RPP ; c(1,1,4)=   -41._RPP/2520._RPP; c(1,2,4)=   199._RPP/2520._RPP  ! stencil 4
        !c(1,0,5)=    -1._RPP/504._RPP ; c(1,1,5)=    49._RPP/2520._RPP; c(1,2,5)=  -221._RPP/2520._RPP  ! stencil 5
        !c(1,0,6)=     1._RPP/252._RPP ; c(1,1,6)=   -19._RPP/508._RPP ; c(1,2,6)=   409._RPP/2520._RPP  ! stencil 6
        !c(1,0,7)=    -7._RPP/504._RPP ; c(1,1,7)=    65._RPP/508._RPP ; c(1,2,7)=  -271._RPP/504._RPP   ! stencil 7
        !c(1,0,8)=     1._RPP/9._RPP   ; c(1,1,8)=  -511._RPP/508._RPP ; c(1,2,8)=  2081._RPP/504._RPP   ! stencil 8
        !!  cell  3                    ;    cell  4                    ;     cell  5
        !c(1,3,0)=  6289._RPP/2520._RPP; c(1,4,0)= -5471._RPP/2520._RPP; c(1,5,0)=  3349._RPP/2520._RPP  ! stencil 0
        !c(1,3,1)= -2531._RPP/2520._RPP; c(1,4,1)=  1879._RPP/2520._RPP; c(1,5,1)= -1061._RPP/2520._RPP  ! stencil 1
        !c(1,3,2)=  2509._RPP/2520._RPP; c(1,4,2)= -1271._RPP/2520._RPP; c(1,5,2)=   619._RPP/2520._RPP  ! stencil 2
        !c(1,3,3)=   275._RPP/504._RPP ; c(1,4,3)=  1879._RPP/2520._RPP; c(1,5,3)=  -641._RPP/2520._RPP  ! stencil 3
        !c(1,3,4)=  -641._RPP/2520._RPP; c(1,4,4)=  1879._RPP/2520._RPP; c(1,5,4)=   275._RPP/504._RPP   ! stencil 4
        !c(1,3,5)=   619._RPP/2520._RPP; c(1,4,5)= -1271._RPP/2520._RPP; c(1,5,5)=  2509._RPP/2520._RPP  ! stencil 5
        !c(1,3,6)= -1061._RPP/2520._RPP; c(1,4,6)=  1879._RPP/2520._RPP; c(1,5,6)= -2531._RPP/2520._RPP  ! stencil 6
        !c(1,3,7)=  3349._RPP/2520._RPP; c(1,4,7)= -5471._RPP/2520._RPP; c(1,5,7)=  6289._RPP/2520._RPP  ! stencil 7
        !c(1,3,8)= -4975._RPP/504._RPP ; c(1,4,8)= 38629._RPP/2520._RPP; c(1,5,8)=-40751._RPP/2520._RPP  ! stencil 8
        !!   cell  6                   ;     cell  7                   ;     cell  8
        !c(1,6,0)=  -271._RPP/504._RPP ; c(1,7,0)=    65._RPP/504._RPP ; c(1,8,0)=    -1._RPP/72._RPP    ! stencil 0
        !c(1,6,1)=   409._RPP/2520._RPP; c(1,7,1)=    19._RPP/504._RPP ; c(1,8,1)=     1._RPP/252._RPP   ! stencil 1
        !c(1,6,2)=  -221._RPP/2520._RPP; c(1,7,2)=    49._RPP/2520._RPP; c(1,8,2)=    -1._RPP/504._RPP   ! stencil 2
        !c(1,6,3)=   199._RPP/2520._RPP; c(1,7,3)=   -41._RPP/2520._RPP; c(1,8,3)=     1._RPP/630._RPP   ! stencil 3
        !c(1,6,4)=   -61._RPP/504._RPP ; c(1,7,4)=    11._RPP/504._RPP ; c(1,8,4)=    -1._RPP/504._RPP   ! stencil 4
        !c(1,6,5)=   955._RPP/2520._RPP; c(1,7,5)=   -25._RPP/504._RPP ; c(1,8,5)=     1._RPP/252._RPP   ! stencil 5
        !c(1,6,6)=  3349._RPP/2520._RPP; c(1,7,6)=   119._RPP/504._RPP ; c(1,8,6)=     -1._RPP/72._RPP   ! stencil 6
        !c(1,6,7)= -5471._RPP/2520._RPP; c(1,7,7)=  4609._RPP/2520._RPP; c(1,8,7)=      7._RPP/63._RPP   ! stencil 7
        !c(1,6,8)= 29809._RPP/2520._RPP; c(1,7,8)=-15551._RPP/2520._RPP; c(1,8,8)=   7129._RPP/2520._RPP ! stencil 8
        !! 2 => right interface (i+1/2)
        !!  cell  0                    ;     cell  1                   ;     cell  2
        !c(2,0,0)=  7129._RPP/2520._RPP; c(2,1,0)=-15551._RPP/2520._RPP; c(2,2,0)= 29809._RPP/2520._RPP  ! stencil 0
        !c(2,0,1)=     7._RPP/63._RPP  ; c(2,1,1)=  4609._RPP/2520._RPP; c(2,2,1)= -5471._RPP/2520._RPP  ! stencil 1
        !c(2,0,2)=    -1._RPP/72._RPP  ; c(2,1,2)=   119._RPP/504._RPP ; c(2,2,2)=  3349._RPP/2520._RPP  ! stencil 2
        !c(2,0,3)=    1._RPP/252._RPP  ; c(2,1,3)=   -25._RPP/504._RPP ; c(2,2,3)=   955._RPP/2520._RPP  ! stencil 3
        !c(2,0,4)=   -1._RPP/504._RPP  ; c(2,1,4)=    11._RPP/504._RPP ; c(2,2,4)=   -61._RPP/504._RPP   ! stencil 4
        !c(2,0,5)=    1._RPP/630._RPP  ; c(2,1,5)=   -41._RPP/2520._RPP; c(2,2,5)=   199._RPP/2520._RPP  ! stencil 5
        !c(2,0,6)=   -1._RPP/504._RPP  ; c(2,1,6)=    49._RPP/2520._RPP; c(2,2,6)=  -221._RPP/2520._RPP  ! stencil 6
        !c(2,0,7)=    1._RPP/252._RPP  ; c(2,1,7)=    19._RPP/504._RPP ; c(2,2,7)=   409._RPP/2520._RPP  ! stencil 7
        !c(2,0,8)=   -1._RPP/72._RPP   ; c(2,1,8)=    65._RPP/504._RPP ; c(2,2,8)=  -271._RPP/504._RPP   ! stencil 8
        !!   cell  3                   ; !  cell  4                    ;     cell  5
        !c(2,3,0)=-40751._RPP/2520._RPP; c(2,4,0)= 38629._RPP/2520._RPP; c(2,5,0)= -4975._RPP/504._RPP   ! stencil 0
        !c(2,3,1)=  6289._RPP/2520._RPP; c(2,4,1)= -5471._RPP/2520._RPP; c(2,5,1)=  3349._RPP/2520._RPP  ! stencil 1
        !c(2,3,2)= -2531._RPP/2520._RPP; c(2,4,2)=  1879._RPP/2520._RPP; c(2,5,2)= -1061._RPP/2520._RPP  ! stencil 2
        !c(2,3,3)=  2509._RPP/2520._RPP; c(2,4,3)= -1271._RPP/2520._RPP; c(2,5,3)=   619._RPP/2520._RPP  ! stencil 3
        !c(2,3,4)=   275._RPP/504._RPP ; c(2,4,4)=  1879._RPP/2520._RPP; c(2,5,4)=  -641._RPP/2520._RPP  ! stencil 4
        !c(2,3,5)=  -641._RPP/2520._RPP; c(2,4,5)=  1879._RPP/2520._RPP; c(2,5,5)=   275._RPP/504._RPP   ! stencil 5
        !c(2,3,6)=   619._RPP/2520._RPP; c(2,4,6)= -1271._RPP/2520._RPP; c(2,5,6)=  2509._RPP/2520._RPP  ! stencil 6
        !c(2,3,7)= -1061._RPP/2520._RPP; c(2,4,7)=  1879._RPP/2520._RPP; c(2,5,7)= -2531._RPP/2520._RPP  ! stencil 7
        !c(2,3,8)=  3349._RPP/2520._RPP; c(2,4,8)= -5471._RPP/2520._RPP; c(2,5,8)=  6289._RPP/2520._RPP  ! stencil 8
        !!   cell  6                   ;     cell  7                   ;     cell  8
        !c(2,6,0)=  2081._RPP/504._RPP ; c(2,7,0)=  -511._RPP/508._RPP ; c(2,8,0)=     1._RPP/9._RPP     ! stencil 0
        !c(2,6,1)=  -271._RPP/504._RPP ; c(2,7,1)=    65._RPP/508._RPP ; c(2,8,1)=    -7._RPP/504._RPP   ! stencil 1
        !c(2,6,2)=   409._RPP/2520._RPP; c(2,7,2)=   -19._RPP/508._RPP ; c(2,8,2)=     1._RPP/252._RPP   ! stencil 2
        !c(2,6,3)=  -221._RPP/2520._RPP; c(2,7,3)=    49._RPP/2520._RPP; c(2,8,3)=    -1._RPP/504._RPP   ! stencil 3
        !c(2,6,4)=   199._RPP/2520._RPP; c(2,7,4)=   -41._RPP/2520._RPP; c(2,8,4)=     1._RPP/630._RPP   ! stencil 4
        !c(2,6,5)=   -61._RPP/504._RPP ; c(2,7,5)=    11._RPP/508._RPP ; c(2,8,5)=    -1._RPP/504._RPP   ! stencil 5
        !c(2,6,6)=   191._RPP/504._RPP ; c(2,7,6)=   -25._RPP/508._RPP ; c(2,8,6)=     1._RPP/252._RPP   ! stencil 6
        !c(2,6,7)=  3349._RPP/2520._RPP; c(2,7,7)=   119._RPP/508._RPP ; c(2,8,7)=    -7._RPP/504._RPP   ! stencil 7
        !c(2,6,8)= -5471._RPP/2520._RPP; c(2,7,8)=  4609._RPP/2520._RPP; c(2,8,8)=     1._RPP/9._RPP     ! stencil 8
        ! 1 => left interface (i-1/2)
        c(1,0,8)=   280._RPP/ 2520._RPP; c(1,0,7)=-   35._RPP/ 2520._RPP; c(1,0,6)=    10._RPP/ 2520._RPP
        c(1,1,8)=- 2555._RPP/ 2520._RPP; c(1,1,7)=   325._RPP/ 2520._RPP; c(1,1,6)=-   95._RPP/ 2520._RPP
        c(1,2,8)= 10405._RPP/ 2520._RPP; c(1,2,7)=- 1355._RPP/ 2520._RPP; c(1,2,6)=   409._RPP/ 2520._RPP
        c(1,3,8)=-24875._RPP/ 2520._RPP; c(1,3,7)=  3349._RPP/ 2520._RPP; c(1,3,6)=- 1061._RPP/ 2520._RPP
        c(1,4,8)= 38629._RPP/ 2520._RPP; c(1,4,7)=- 5471._RPP/ 2520._RPP; c(1,4,6)=  1879._RPP/ 2520._RPP
        c(1,5,8)=-40751._RPP/ 2520._RPP; c(1,5,7)=  6289._RPP/ 2520._RPP; c(1,5,6)=- 2531._RPP/ 2520._RPP
        c(1,6,8)= 29809._RPP/ 2520._RPP; c(1,6,7)=- 5471._RPP/ 2520._RPP; c(1,6,6)=  3349._RPP/ 2520._RPP
        c(1,7,8)=-15551._RPP/ 2520._RPP; c(1,7,7)=  4609._RPP/ 2520._RPP; c(1,7,6)=   595._RPP/ 2520._RPP
        c(1,8,8)=  7129._RPP/ 2520._RPP; c(1,8,7)=   280._RPP/ 2520._RPP; c(1,8,6)=-   35._RPP/ 2520._RPP

        c(1,0,5)=-    5._RPP/ 2520._RPP; c(1,0,4)=     4._RPP/ 2520._RPP; c(1,0,3)=-    5._RPP/ 2520._RPP
        c(1,1,5)=    49._RPP/ 2520._RPP; c(1,1,4)=-   41._RPP/ 2520._RPP; c(1,1,3)=    55._RPP/ 2520._RPP
        c(1,2,5)=-  221._RPP/ 2520._RPP; c(1,2,4)=   199._RPP/ 2520._RPP; c(1,2,3)=-  305._RPP/ 2520._RPP
        c(1,3,5)=   619._RPP/ 2520._RPP; c(1,3,4)=-  641._RPP/ 2520._RPP; c(1,3,3)=  1375._RPP/ 2520._RPP
        c(1,4,5)=- 1271._RPP/ 2520._RPP; c(1,4,4)=  1879._RPP/ 2520._RPP; c(1,4,3)=  1879._RPP/ 2520._RPP
        c(1,5,5)=  2509._RPP/ 2520._RPP; c(1,5,4)=  1375._RPP/ 2520._RPP; c(1,5,3)=-  641._RPP/ 2520._RPP
        c(1,6,5)=   955._RPP/ 2520._RPP; c(1,6,4)=-  305._RPP/ 2520._RPP; c(1,6,3)=   199._RPP/ 2520._RPP
        c(1,7,5)=-  125._RPP/ 2520._RPP; c(1,7,4)=    55._RPP/ 2520._RPP; c(1,7,3)=-   41._RPP/ 2520._RPP
        c(1,8,5)=    10._RPP/ 2520._RPP; c(1,8,4)=-    5._RPP/ 2520._RPP; c(1,8,3)=     4._RPP/ 2520._RPP

        c(1,0,2)=    10._RPP/ 2520._RPP; c(1,0,1)=-   35._RPP/ 2520._RPP; c(1,0,0)=   280._RPP/ 2520._RPP
        c(1,1,2)=-  125._RPP/ 2520._RPP; c(1,1,1)=   595._RPP/ 2520._RPP; c(1,1,0)=  4609._RPP/ 2520._RPP
        c(1,2,2)=   955._RPP/ 2520._RPP; c(1,2,1)=  3349._RPP/ 2520._RPP; c(1,2,0)=- 5471._RPP/ 2520._RPP
        c(1,3,2)=  2509._RPP/ 2520._RPP; c(1,3,1)=- 2531._RPP/ 2520._RPP; c(1,3,0)=  6289._RPP/ 2520._RPP
        c(1,4,2)=- 1271._RPP/ 2520._RPP; c(1,4,1)=  1879._RPP/ 2520._RPP; c(1,4,0)=- 5471._RPP/ 2520._RPP
        c(1,5,2)=   619._RPP/ 2520._RPP; c(1,5,1)=- 1061._RPP/ 2520._RPP; c(1,5,0)=  3349._RPP/ 2520._RPP
        c(1,6,2)=-  221._RPP/ 2520._RPP; c(1,6,1)=   409._RPP/ 2520._RPP; c(1,6,0)=- 1355._RPP/ 2520._RPP
        c(1,7,2)=    49._RPP/ 2520._RPP; c(1,7,1)=-   95._RPP/ 2520._RPP; c(1,7,0)=   325._RPP/ 2520._RPP
        c(1,8,2)=-    5._RPP/ 2520._RPP; c(1,8,1)=    10._RPP/ 2520._RPP; c(1,8,0)=-   35._RPP/ 2520._RPP

        ! 2 => right interface (i+1/2)
        c(2,8,0)=   280._RPP/ 2520._RPP; c(2,8,1)=-   35._RPP/ 2520._RPP; c(2,8,2)=    10._RPP/ 2520._RPP
        c(2,7,0)=- 2555._RPP/ 2520._RPP; c(2,7,1)=   325._RPP/ 2520._RPP; c(2,7,2)=-   95._RPP/ 2520._RPP
        c(2,6,0)= 10405._RPP/ 2520._RPP; c(2,6,1)=- 1355._RPP/ 2520._RPP; c(2,6,2)=   409._RPP/ 2520._RPP
        c(2,5,0)=-24875._RPP/ 2520._RPP; c(2,5,1)=  3349._RPP/ 2520._RPP; c(2,5,2)=- 1061._RPP/ 2520._RPP
        c(2,4,0)= 38629._RPP/ 2520._RPP; c(2,4,1)=- 5471._RPP/ 2520._RPP; c(2,4,2)=  1879._RPP/ 2520._RPP
        c(2,3,0)=-40751._RPP/ 2520._RPP; c(2,3,1)=  6289._RPP/ 2520._RPP; c(2,3,2)=- 2531._RPP/ 2520._RPP
        c(2,2,0)= 29809._RPP/ 2520._RPP; c(2,2,1)=- 5471._RPP/ 2520._RPP; c(2,2,2)=  3349._RPP/ 2520._RPP
        c(2,1,0)=-15551._RPP/ 2520._RPP; c(2,1,1)=  4609._RPP/ 2520._RPP; c(2,1,2)=   595._RPP/ 2520._RPP
        c(2,0,0)=  7129._RPP/ 2520._RPP; c(2,0,1)=   280._RPP/ 2520._RPP; c(2,0,2)=-   35._RPP/ 2520._RPP

        c(2,8,3)=-    5._RPP/ 2520._RPP; c(2,8,4)=     4._RPP/ 2520._RPP; c(2,8,5)=-    5._RPP/ 2520._RPP
        c(2,7,3)=    49._RPP/ 2520._RPP; c(2,7,4)=-   41._RPP/ 2520._RPP; c(2,7,5)=    55._RPP/ 2520._RPP
        c(2,6,3)=-  221._RPP/ 2520._RPP; c(2,6,4)=   199._RPP/ 2520._RPP; c(2,6,5)=-  305._RPP/ 2520._RPP
        c(2,5,3)=   619._RPP/ 2520._RPP; c(2,5,4)=-  641._RPP/ 2520._RPP; c(2,5,5)=  1375._RPP/ 2520._RPP
        c(2,4,3)=- 1271._RPP/ 2520._RPP; c(2,4,4)=  1879._RPP/ 2520._RPP; c(2,4,5)=  1879._RPP/ 2520._RPP
        c(2,3,3)=  2509._RPP/ 2520._RPP; c(2,3,4)=  1375._RPP/ 2520._RPP; c(2,3,5)=-  641._RPP/ 2520._RPP
        c(2,2,3)=   955._RPP/ 2520._RPP; c(2,2,4)=-  305._RPP/ 2520._RPP; c(2,2,5)=   199._RPP/ 2520._RPP
        c(2,1,3)=-  125._RPP/ 2520._RPP; c(2,1,4)=    55._RPP/ 2520._RPP; c(2,1,5)=-   41._RPP/ 2520._RPP
        c(2,0,3)=    10._RPP/ 2520._RPP; c(2,0,4)=-    5._RPP/ 2520._RPP; c(2,0,5)=     4._RPP/ 2520._RPP

        c(2,8,6)=    10._RPP/ 2520._RPP; c(2,8,7)=-   35._RPP/ 2520._RPP; c(2,8,8)=   280._RPP/ 2520._RPP
        c(2,7,6)=-  125._RPP/ 2520._RPP; c(2,7,7)=   595._RPP/ 2520._RPP; c(2,7,8)=  4609._RPP/ 2520._RPP
        c(2,6,6)=   955._RPP/ 2520._RPP; c(2,6,7)=  3349._RPP/ 2520._RPP; c(2,6,8)=- 5471._RPP/ 2520._RPP
        c(2,5,6)=  2509._RPP/ 2520._RPP; c(2,5,7)=- 2531._RPP/ 2520._RPP; c(2,5,8)=  6289._RPP/ 2520._RPP
        c(2,4,6)=- 1271._RPP/ 2520._RPP; c(2,4,7)=  1879._RPP/ 2520._RPP; c(2,4,8)=- 5471._RPP/ 2520._RPP
        c(2,3,6)=   619._RPP/ 2520._RPP; c(2,3,7)=- 1061._RPP/ 2520._RPP; c(2,3,8)=  3349._RPP/ 2520._RPP
        c(2,2,6)=-  221._RPP/ 2520._RPP; c(2,2,7)=   409._RPP/ 2520._RPP; c(2,2,8)=- 1355._RPP/ 2520._RPP
        c(2,1,6)=    49._RPP/ 2520._RPP; c(2,1,7)=-   95._RPP/ 2520._RPP; c(2,1,8)=   325._RPP/ 2520._RPP
        c(2,0,6)=-    5._RPP/ 2520._RPP; c(2,0,7)=    10._RPP/ 2520._RPP; c(2,0,8)=-   35._RPP/ 2520._RPP
      endselect
  endassociate
  endsubroutine create

  pure subroutine compute_int(self, stencil, values)
  !< Compute interpolations (interpolate).
  class(interpolations_rec_js), intent(in)  :: self               !< Interpolations.
  real(RPP),                    intent(in)  :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(RPP),                    intent(out) :: values(0:)         !< Interpolations values.
  ! empty procedure
  endsubroutine compute_int

  pure subroutine compute_rec(self, stencil, values)
  !< Compute interpolations (reconstruct).
  class(interpolations_rec_js), intent(in)  :: self                  !< Interpolations.
  real(RPP),                    intent(in)  :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(RPP),                    intent(out) :: values(1:, 0:)        !< Interpolations values.
  integer(I_P)                              :: s1                    !< Counter.
  integer(I_P)                              :: s2                    !< Counter.
  integer(I_P)                              :: f                     !< Counter.

  values = 0._RPP
  do s1=0, self%S - 1 ! stencils loop
    do s2=0, self%S - 1 ! values loop
      do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        values(f, s1) = values(f, s1) + self%coef(f, s2, s1) * stencil(f, -s2 + s1)
      enddo
    enddo
  enddo
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(interpolations_rec_js), intent(in)           :: self             !< Interpolations.
  character(len=*),             intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                      :: string           !< String-description.
  character(len=:), allocatable                      :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                        :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu beta interpolations object for reconstruction:'//NL
  string = prefix_//string//'  - S   = '//trim(str(self%S))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy interpolations.
  class(interpolations_rec_js), intent(inout) :: self !< Interpolations.

  call self%destroy_
  if (allocated(self%coef)) deallocate(self%coef)
  endsubroutine destroy
endmodule wenoof_interpolations_rec_js
