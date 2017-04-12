!< Jiang-Shu (Lagrange) interpolations object for derivative reconstruction.
module wenoof_interpolations_rec_js
!< Jiang-Shu (Lagrange) interpolations object for derivative reconstruction.
!<
!< @note The provided interpolations implement the Lagrange interpolations defined in *Efficient Implementation
!< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P, str
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
  real(R_P), allocatable :: coef(:,:,:) !< Polynomial coefficients [1:2,0:S-1,0:S-1].
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
        c(1,0,0) =  0.5_R_P; c(1,1,0) =  0.5_R_P ! stencil 0
        c(1,0,1) = -0.5_R_P; c(1,1,1) =  1.5_R_P ! stencil 1
        ! 2 => right interface (i+1/2)
        !  cell  0           ;    cell  1
        c(2,0,0) =  1.5_R_P; c(2,1,0) = -0.5_R_P ! stencil 0
        c(2,0,1) =  0.5_R_P; c(2,1,1) =  0.5_R_P ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        !  cell  0                 ;    cell  1                 ;    cell  2
        c(1,0,0) =  1._R_P/3._R_P; c(1,1,0) =  5._R_P/6._R_P; c(1,2,0) = -1._R_P/6._R_P ! stencil 0
        c(1,0,1) = -1._R_P/6._R_P; c(1,1,1) =  5._R_P/6._R_P; c(1,2,1) =  1._R_P/3._R_P ! stencil 1
        c(1,0,2) =  1._R_P/3._R_P; c(1,1,2) = -7._R_P/6._R_P; c(1,2,2) = 11._R_P/6._R_P ! stencil 2
        ! 2 => right interface (i+1/2)
        !  cell  0                 ;    cell  1                 ;    cell  2
        c(2,0,0) = 11._R_P/6._R_P; c(2,1,0) = -7._R_P/6._R_P; c(2,2,0) =  1._R_P/3._R_P ! stencil 0
        c(2,0,1) =  1._R_P/3._R_P; c(2,1,1) =  5._R_P/6._R_P; c(2,2,1) = -1._R_P/6._R_P ! stencil 1
        c(2,0,2) = -1._R_P/6._R_P; c(2,1,2) =  5._R_P/6._R_P; c(2,2,2) =  1._R_P/3._R_P ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        !  cell  0               ;    cell  1               ;    cell  2               ;    cell  3
        c(1,0,0)=  1._R_P/4._R_P ; c(1,1,0)= 13._R_P/12._R_P; c(1,2,0)= -5._R_P/12._R_P; c(1,3,0)=  1._R_P/12._R_P ! stencil 0
        c(1,0,1)= -1._R_P/12._R_P; c(1,1,1)=  7._R_P/12._R_P; c(1,2,1)=  7._R_P/12._R_P; c(1,3,1)= -1._R_P/12._R_P ! stencil 1
        c(1,0,2)=  1._R_P/12._R_P; c(1,1,2)= -5._R_P/12._R_P; c(1,2,2)= 13._R_P/12._R_P; c(1,3,2)=  1._R_P/4._R_P  ! stencil 2
        c(1,0,3)= -1._R_P/4._R_P ; c(1,1,3)= 13._R_P/12._R_P; c(1,2,3)=-23._R_P/12._R_P; c(1,3,3)= 25._R_P/12._R_P ! stencil 3
        ! 2 => right interface (i+1/2)
        !  cell  0               ;    cell  1               ;   cell  2                ;    cell  3
        c(2,0,0)= 25._R_P/12._R_P; c(2,1,0)=-23._R_P/12._R_P; c(2,2,0)= 13._R_P/12._R_P; c(2,3,0)= -1._R_P/4._R_P  ! stencil 0
        c(2,0,1)=  1._R_P/4._R_P ; c(2,1,1)= 13._R_P/12._R_P; c(2,2,1)= -5._R_P/12._R_P; c(2,3,1)=  1._R_P/12._R_P ! stencil 1
        c(2,0,2)= -1._R_P/12._R_P; c(2,1,2)=  7._R_P/12._R_P; c(2,2,2)=  7._R_P/12._R_P; c(2,3,2)= -1._R_P/12._R_P ! stencil 2
        c(2,0,3)=  1._R_P/12._R_P; c(2,1,3)= -5._R_P/12._R_P; c(2,2,3)= 13._R_P/12._R_P; c(2,3,3)=  1._R_P/4._R_P  ! stencil 3
      case(5) ! 9th order
        ! 1 => left interface (i-1/2)
        !  cell  0                ;    cell  1                ;    cell  2                ;    cell  3
        c(1,0,0)=   1._R_P/5._R_P ; c(1,1,0)=  77._R_P/60._R_P; c(1,2,0)= -43._R_P/60._R_P; c(1,3,0)=  17._R_P/60._R_P  ! stencil 0
        c(1,0,1)=  -1._R_P/20._R_P; c(1,1,1)=   9._R_P/20._R_P; c(1,2,1)=  47._R_P/60._R_P; c(1,3,1)= -13._R_P/60._R_P  ! stencil 1
        c(1,0,2)=   1._R_P/30._R_P; c(1,1,2)= -13._R_P/60._R_P; c(1,2,2)=  47._R_P/60._R_P; c(1,3,2)=   9._R_P/20._R_P  ! stencil 2
        c(1,0,3)=  -1._R_P/20._R_P; c(1,1,3)=  17._R_P/60._R_P; c(1,2,3)= -43._R_P/60._R_P; c(1,3,3)=  77._R_P/60._R_P  ! stencil 3
        c(1,0,4)=   1._R_P/5._R_P ; c(1,1,4)= -21._R_P/20._R_P; c(1,2,4)= 137._R_P/60._R_P; c(1,3,4)=-163._R_P/60._R_P  ! stencil 4
        !  cell  4
        c(1,4,0)=  -1._R_P/20._R_P  ! stencil 0
        c(1,4,1)=   1._R_P/30._R_P  ! stencil 1
        c(1,4,2)=  -1._R_P/20._R_P  ! stencil 2
        c(1,4,3)=   1._R_P/5._R_P   ! stencil 3
        c(1,4,4)= 137._R_P/60._R_P  ! stencil 4
        ! 2 => right interface (i+1/2)
        !  cell  0               ;    cell  1               ;   cell  2                ;    cell  3
        c(2,0,0)= 137._R_P/60._R_P; c(2,1,0)=-163._R_P/60._R_P; c(2,2,0)= 137._R_P/60._R_P; c(2,3,0)= -21._R_P/20._R_P  ! stencil 0
        c(2,0,1)=   1._R_P/5._R_P ; c(2,1,1)=  77._R_P/60._R_P; c(2,2,1)= -43._R_P/60._R_P; c(2,3,1)=  17._R_P/60._R_P  ! stencil 1
        c(2,0,2)=  -1._R_P/20._R_P; c(2,1,2)=   9._R_P/20._R_P; c(2,2,2)=  47._R_P/60._R_P; c(2,3,2)= -13._R_P/60._R_P  ! stencil 2
        c(2,0,3)=   1._R_P/30._R_P; c(2,1,3)= -13._R_P/60._R_P; c(2,2,3)=  47._R_P/60._R_P; c(2,3,3)=   9._R_P/20._R_P  ! stencil 3
        c(2,0,4)=  -1._R_P/20._R_P; c(2,1,4)=  17._R_P/60._R_P; c(2,2,4)= -43._R_P/60._R_P; c(2,3,4)=  77._R_P/60._R_P  ! stencil 4
        !  cell  4
        c(2,4,0)=   1._R_P/5._R_P  ! stencil 0
        c(2,4,1)=  -1._R_P/20._R_P ! stencil 1
        c(2,4,2)=   1._R_P/30._R_P ! stencil 2
        c(2,4,3)=  -1._R_P/20._R_P ! stencil 3
        c(2,4,4)=   1._R_P/5._R_P  ! stencil 4
      case(6) ! 11th order
        ! 1 => left interface (i-1/2)
        !  cell  0                ;    cell  1                ;    cell  2                ;    cell  3
        c(1,0,0)=   1._R_P/6._R_P ; c(1,1,0)=  29._R_P/20._R_P; c(1,2,0)= -21._R_P/20._R_P; c(1,3,0)=  37._R_P/60._R_P  ! stencil 0
        c(1,0,1)=  -1._R_P/30._R_P; c(1,1,1)=  11._R_P/30._R_P; c(1,2,1)=  19._R_P/20._R_P; c(1,3,1)= -23._R_P/60._R_P  ! stencil 1
        c(1,0,2)=   1._R_P/60._R_P; c(1,1,2)=  -2._R_P/15._R_P; c(1,2,2)=  37._R_P/60._R_P; c(1,3,2)=  37._R_P/60._R_P  ! stencil 2
        c(1,0,3)=  -1._R_P/60._R_P; c(1,1,3)=   7._R_P/60._R_P; c(1,2,3)= -23._R_P/60._R_P; c(1,3,3)=  19._R_P/20._R_P  ! stencil 3
        c(1,0,4)=   1._R_P/30._R_P; c(1,1,4)= -13._R_P/60._R_P; c(1,2,4)=  37._R_P/60._R_P; c(1,3,4)= -21._R_P/20._R_P  ! stencil 4
        c(1,0,5)=  -1._R_P/6._R_P ; c(1,1,5)=  31._R_P/30._R_P; c(1,2,5)=-163._R_P/60._R_P; c(1,3,5)=  79._R_P/20._R_P  ! stencil 5
        !  cell  4                ;    cell  5
        c(1,4,0)= -13._R_P/60._R_P; c(1,5,0)=   1._R_P/30._R_P  ! stencil 0
        c(1,4,1)=   7._R_P/60._R_P; c(1,5,1)=  -1._R_P/60._R_P  ! stencil 1
        c(1,4,2)=  -2._R_P/15._R_P; c(1,5,2)=   1._R_P/60._R_P  ! stencil 2
        c(1,4,3)=  11._R_P/30._R_P; c(1,5,3)=  -1._R_P/30._R_P  ! stencil 3
        c(1,4,4)=  29._R_P/20._R_P; c(1,5,4)=   1._R_P/6._R_P   ! stencil 4
        c(1,4,5)= -71._R_P/20._R_P; c(1,5,5)=  49._R_P/20._R_P  ! stencil 5
        ! 2 => right interface (i+1/2)
        !  cell  0                ;    cell  1                ;   cell  2                 ;    cell  3
        c(2,0,0)=  49._R_P/20._R_P; c(2,1,0)= -71._R_P/20._R_P; c(2,2,0)=  79._R_P/20._R_P; c(2,3,0)=-163._R_P/60._R_P  ! stencil 0
        c(2,0,1)=   1._R_P/6._R_P ; c(2,1,1)=  29._R_P/20._R_P; c(2,2,1)= -21._R_P/20._R_P; c(2,3,1)=  37._R_P/60._R_P  ! stencil 1
        c(2,0,2)=  -1._R_P/30._R_P; c(2,1,2)=  11._R_P/30._R_P; c(2,2,2)=  19._R_P/20._R_P; c(2,3,2)= -23._R_P/60._R_P  ! stencil 2
        c(2,0,3)=   1._R_P/60._R_P; c(2,1,3)=  -2._R_P/15._R_P; c(2,2,3)=  37._R_P/60._R_P; c(2,3,3)=  37._R_P/60._R_P  ! stencil 3
        c(2,0,4)=  -1._R_P/60._R_P; c(2,1,4)=   7._R_P/60._R_P; c(2,2,4)= -23._R_P/60._R_P; c(2,3,4)=  19._R_P/20._R_P  ! stencil 4
        c(2,0,5)=   1._R_P/30._R_P; c(2,1,5)= -13._R_P/60._R_P; c(2,2,5)=  37._R_P/60._R_P; c(2,3,5)= -21._R_P/20._R_P  ! stencil 5
        !  cell  4                ;    cell  5
        c(2,4,0)=  31._R_P/30._R_P; c(2,5,0)=  -1._R_P/6._R_P   ! stencil 0
        c(2,4,1)= -13._R_P/60._R_P; c(2,5,1)=   1._R_P/30._R_P  ! stencil 1
        c(2,4,2)=   7._R_P/60._R_P; c(2,5,2)=  -1._R_P/60._R_P  ! stencil 2
        c(2,4,3)=  -2._R_P/15._R_P; c(2,5,3)=   1._R_P/60._R_P  ! stencil 3
        c(2,4,4)=  11._R_P/30._R_P; c(2,5,4)=  -1._R_P/30._R_P  ! stencil 4
        c(2,4,5)=  29._R_P/20._R_P; c(2,5,5)=   1._R_P/6._R_P   ! stencil 5
      case(7) ! 13th order
        ! 1 => left interface (i-1/2)
        !  stencil  6                ;    stencil  5               ;    stencil  4
        c(1,0,6)=   60._R_P/ 420._R_P; c(1,0,5)= -10._R_P/ 420._R_P; c(1,0,4)=   4._R_P/ 420._R_P  ! cell 0
        c(1,1,6)= -430._R_P/ 420._R_P; c(1,1,5)=  74._R_P/ 420._R_P; c(1,1,4)= -31._R_P/ 420._R_P  ! cell 1
        c(1,2,6)= 1334._R_P/ 420._R_P; c(1,2,5)=-241._R_P/ 420._R_P; c(1,2,4)= 109._R_P/ 420._R_P  ! cell 2
        c(1,3,6)=-2341._R_P/ 420._R_P; c(1,3,5)= 459._R_P/ 420._R_P; c(1,3,4)=-241._R_P/ 420._R_P  ! cell 3
        c(1,4,6)= 2559._R_P/ 420._R_P; c(1,4,5)=-591._R_P/ 420._R_P; c(1,4,4)= 459._R_P/ 420._R_P  ! cell 4
        c(1,5,6)=-1851._R_P/ 420._R_P; c(1,5,5)= 669._R_P/ 420._R_P; c(1,5,4)= 130._R_P/ 420._R_P  ! cell 5
        c(1,6,6)= 1089._R_P/ 420._R_P; c(1,6,5)=  60._R_P/ 420._R_P; c(1,6,4)=- 10._R_P/ 420._R_P  ! cell 6
        !  stencil  3                ;    stencil  2               ;    stencil  1
        c(1,0,3)=   -3._R_P/ 420._R_P; c(1,0,2)=   4._R_P/ 420._R_P; c(1,0,1)=- 10._R_P/ 420._R_P  ! cell 0
        c(1,1,3)=   25._R_P/ 420._R_P; c(1,1,2)= -38._R_P/ 420._R_P; c(1,1,1)= 130._R_P/ 420._R_P  ! cell 1
        c(1,2,3)= -101._R_P/ 420._R_P; c(1,2,2)= 214._R_P/ 420._R_P; c(1,2,1)= 459._R_P/ 420._R_P  ! cell 2
        c(1,3,3)=  319._R_P/ 420._R_P; c(1,3,2)= 319._R_P/ 420._R_P; c(1,3,1)=-241._R_P/ 420._R_P  ! cell 3
        c(1,4,3)=  214._R_P/ 420._R_P; c(1,4,2)=-101._R_P/ 420._R_P; c(1,4,1)= 109._R_P/ 420._R_P  ! cell 4
        c(1,5,3)=  -38._R_P/ 420._R_P; c(1,5,2)=  25._R_P/ 420._R_P; c(1,5,1)=- 31._R_P/ 420._R_P  ! cell 5
        c(1,6,3)=    4._R_P/ 420._R_P; c(1,6,2)=-  3._R_P/ 420._R_P; c(1,6,1)=   4._R_P/ 420._R_P  ! cell 6
        !  stencil  0
        c(1,0,0)=   60._R_P/ 420._R_P  ! cell 0
        c(1,1,0)=  669._R_P/ 420._R_P  ! cell 1
        c(1,2,0)= -591._R_P/ 420._R_P  ! cell 2
        c(1,3,0)=  459._R_P/ 420._R_P  ! cell 3
        c(1,4,0)= -241._R_P/ 420._R_P  ! cell 4
        c(1,5,0)=   74._R_P/ 420._R_P  ! cell 5
        c(1,6,0)= - 10._R_P/ 420._R_P  ! cell 6
        ! 2 => right interface (i+1/2)
        !  stencil  0                ;    stencil  1               ;    stencil  2
        c(2,6,0)=   60._R_P/ 420._R_P; c(2,6,1)= -10._R_P/ 420._R_P; c(2,6,2)=   4._R_P/ 420._R_P  ! cell 6
        c(2,5,0)= -430._R_P/ 420._R_P; c(2,5,1)=  74._R_P/ 420._R_P; c(2,5,2)= -31._R_P/ 420._R_P  ! cell 5
        c(2,4,0)= 1334._R_P/ 420._R_P; c(2,4,1)=-241._R_P/ 420._R_P; c(2,4,2)= 109._R_P/ 420._R_P  ! cell 4
        c(2,3,0)=-2341._R_P/ 420._R_P; c(2,3,1)= 459._R_P/ 420._R_P; c(2,3,2)=-241._R_P/ 420._R_P  ! cell 3
        c(2,2,0)= 2559._R_P/ 420._R_P; c(2,2,1)=-591._R_P/ 420._R_P; c(2,2,2)= 459._R_P/ 420._R_P  ! cell 2
        c(2,1,0)=-1851._R_P/ 420._R_P; c(2,1,1)= 669._R_P/ 420._R_P; c(2,1,2)= 130._R_P/ 420._R_P  ! cell 1
        c(2,0,0)= 1089._R_P/ 420._R_P; c(2,0,1)=  60._R_P/ 420._R_P; c(2,0,2)= -10._R_P/ 420._R_P  ! cell 0
        !  stencil  3                ;    stencil  4               ;    stencil  5
        c(2,6,3)=   -3._R_P/ 420._R_P; c(2,6,4)=   4._R_P/ 420._R_P; c(2,6,5)= -10._R_P/ 420._R_P  ! cell 6
        c(2,5,3)=   25._R_P/ 420._R_P; c(2,5,4)= -38._R_P/ 420._R_P; c(2,5,5)= 130._R_P/ 420._R_P  ! cell 5
        c(2,4,3)= -101._R_P/ 420._R_P; c(2,4,4)= 214._R_P/ 420._R_P; c(2,4,5)= 459._R_P/ 420._R_P  ! cell 4
        c(2,3,3)=  319._R_P/ 420._R_P; c(2,3,4)= 319._R_P/ 420._R_P; c(2,3,5)=-241._R_P/ 420._R_P  ! cell 3
        c(2,2,3)=  214._R_P/ 420._R_P; c(2,2,4)=-101._R_P/ 420._R_P; c(2,2,5)= 109._R_P/ 420._R_P  ! cell 2
        c(2,1,3)=  -38._R_P/ 420._R_P; c(2,1,4)=  25._R_P/ 420._R_P; c(2,1,5)= -31._R_P/ 420._R_P  ! cell 1
        c(2,0,3)=    4._R_P/ 420._R_P; c(2,0,4)=  -3._R_P/ 420._R_P; c(2,0,5)=   4._R_P/ 420._R_P  ! cell 0
        !  stencil  6
        c(2,6,6)=  60._R_P/ 420._R_P  ! cell 6
        c(2,5,6)= 669._R_P/ 420._R_P  ! cell 5
        c(2,4,6)=-591._R_P/ 420._R_P  ! cell 4
        c(2,3,6)= 459._R_P/ 420._R_P  ! cell 3
        c(2,2,6)=-241._R_P/ 420._R_P  ! cell 2
        c(2,1,6)=  74._R_P/ 420._R_P  ! cell 1
        c(2,0,6)=- 10._R_P/ 420._R_P  ! cell 0
      case(8) ! 15th order
        ! 1 => left interface (i-1/2)
        !  stencil  7                ;    stencil  6                ;    stencil  5
        c(1,0,7)= -105._R_P/ 840._R_P; c(1,0,6)=   15._R_P/ 840._R_P; c(1,0,5)=-   5._R_P/ 840._R_P  ! cell 0
        c(1,1,7)=  855._R_P/ 840._R_P; c(1,1,6)= -125._R_P/ 840._R_P; c(1,1,5)=   43._R_P/ 840._R_P  ! cell 1
        c(1,2,7)=-3065._R_P/ 840._R_P; c(1,2,6)=  463._R_P/ 840._R_P; c(1,2,5)= -167._R_P/ 840._R_P  ! cell 2
        c(1,3,7)= 6343._R_P/ 840._R_P; c(1,3,6)=-1007._R_P/ 840._R_P; c(1,3,5)=  393._R_P/ 840._R_P  ! cell 3
        c(1,4,7)=-8357._R_P/ 840._R_P; c(1,4,6)= 1443._R_P/ 840._R_P; c(1,4,5)= -657._R_P/ 840._R_P  ! cell 4
        c(1,5,7)= 7323._R_P/ 840._R_P; c(1,5,6)=-1497._R_P/ 840._R_P; c(1,5,5)= 1023._R_P/ 840._R_P  ! cell 5
        c(1,6,7)=-4437._R_P/ 840._R_P; c(1,6,6)= 1443._R_P/ 840._R_P; c(1,6,5)=  225._R_P/ 840._R_P  ! cell 6
        c(1,7,7)= 2283._R_P/ 840._R_P; c(1,7,6)=  105._R_P/ 840._R_P; c(1,7,5)=-  15._R_P/ 840._R_P  ! cell 7
        !  stencil  4                ;    stencil  3                ;    stencil  2
        c(1,0,4)=    3._R_P/ 840._R_P; c(1,0,3)=   -3._R_P/ 840._R_P; c(1,0,2)=    5._R_P/ 840._R_P  ! cell 0
        c(1,1,4)=  -27._R_P/ 840._R_P; c(1,1,3)=   29._R_P/ 840._R_P; c(1,1,2)=  -55._R_P/ 840._R_P  ! cell 1
        c(1,2,4)=  113._R_P/ 840._R_P; c(1,2,3)= -139._R_P/ 840._R_P; c(1,2,2)=  365._R_P/ 840._R_P  ! cell 2
        c(1,3,4)= -307._R_P/ 840._R_P; c(1,3,3)=  533._R_P/ 840._R_P; c(1,3,2)=  743._R_P/ 840._R_P  ! cell 3
        c(1,4,4)=  743._R_P/ 840._R_P; c(1,4,3)=  533._R_P/ 840._R_P; c(1,4,2)= -307._R_P/ 840._R_P  ! cell 4
        c(1,5,4)=  365._R_P/ 840._R_P; c(1,5,3)= -139._R_P/ 840._R_P; c(1,5,2)=  113._R_P/ 840._R_P  ! cell 5
        c(1,6,4)=  -55._R_P/ 840._R_P; c(1,6,3)=   29._R_P/ 840._R_P; c(1,6,2)=  -27._R_P/ 840._R_P  ! cell 6
        c(1,7,4)=    5._R_P/ 840._R_P; c(1,7,3)=   -3._R_P/ 840._R_P; c(1,7,2)=    3._R_P/ 840._R_P  ! cell 7
        !  stencil  1                ;    stencil  0
        c(1,0,1)=  -15._R_P/ 840._R_P; c(1,0,0)=  105._R_P/ 840._R_P  ! cell 0
        c(1,1,1)=  225._R_P/ 840._R_P; c(1,1,0)= 1443._R_P/ 840._R_P  ! cell 1
        c(1,2,1)= 1023._R_P/ 840._R_P; c(1,2,0)=-1497._R_P/ 840._R_P  ! cell 2
        c(1,3,1)= -657._R_P/ 840._R_P; c(1,3,0)= 1443._R_P/ 840._R_P  ! cell 3
        c(1,4,1)=  393._R_P/ 840._R_P; c(1,4,0)=-1007._R_P/ 840._R_P  ! cell 4
        c(1,5,1)= -167._R_P/ 840._R_P; c(1,5,0)=  463._R_P/ 840._R_P  ! cell 5
        c(1,6,1)=   43._R_P/ 840._R_P; c(1,6,0)= -125._R_P/ 840._R_P  ! cell 6
        c(1,7,1)=   -5._R_P/ 840._R_P; c(1,7,0)=   15._R_P/ 840._R_P  ! cell 7
        ! 2 => right interface (i+1/2)
        !  stencil  0                ;    stencil  1                ;    stencil  2
        c(2,7,0)= -105._R_P/ 840._R_P; c(2,7,1)=   15._R_P/ 840._R_P; c(2,7,2)=-   5._R_P/ 840._R_P  ! cell 7
        c(2,6,0)=  855._R_P/ 840._R_P; c(2,6,1)= -125._R_P/ 840._R_P; c(2,6,2)=   43._R_P/ 840._R_P  ! cell 6
        c(2,5,0)=-3065._R_P/ 840._R_P; c(2,5,1)=  463._R_P/ 840._R_P; c(2,5,2)= -167._R_P/ 840._R_P  ! cell 5
        c(2,4,0)= 6343._R_P/ 840._R_P; c(2,4,1)=-1007._R_P/ 840._R_P; c(2,4,2)=  393._R_P/ 840._R_P  ! cell 4
        c(2,3,0)=-8357._R_P/ 840._R_P; c(2,3,1)= 1443._R_P/ 840._R_P; c(2,3,2)= -657._R_P/ 840._R_P  ! cell 3
        c(2,2,0)= 7323._R_P/ 840._R_P; c(2,2,1)=-1497._R_P/ 840._R_P; c(2,2,2)= 1023._R_P/ 840._R_P  ! cell 2
        c(2,1,0)=-4437._R_P/ 840._R_P; c(2,1,1)= 1443._R_P/ 840._R_P; c(2,1,2)=  225._R_P/ 840._R_P  ! cell 1
        c(2,0,0)= 2283._R_P/ 840._R_P; c(2,0,1)=  105._R_P/ 840._R_P; c(2,0,2)=  -15._R_P/ 840._R_P  ! cell 0
        !  stencil  3                ;    stencil  4                ;    stencil  5
        c(2,7,3)=    3._R_P/ 840._R_P; c(2,7,4)=   -3._R_P/ 840._R_P; c(2,7,5)=    5._R_P/ 840._R_P  ! cell 7
        c(2,6,3)=  -27._R_P/ 840._R_P; c(2,6,4)=   29._R_P/ 840._R_P; c(2,6,5)=  -55._R_P/ 840._R_P  ! cell 6
        c(2,5,3)=  113._R_P/ 840._R_P; c(2,5,4)= -139._R_P/ 840._R_P; c(2,5,5)=  365._R_P/ 840._R_P  ! cell 5
        c(2,4,3)= -307._R_P/ 840._R_P; c(2,4,4)=  533._R_P/ 840._R_P; c(2,4,5)=  743._R_P/ 840._R_P  ! cell 4
        c(2,3,3)=  743._R_P/ 840._R_P; c(2,3,4)=  533._R_P/ 840._R_P; c(2,3,5)= -307._R_P/ 840._R_P  ! cell 3
        c(2,2,3)=  365._R_P/ 840._R_P; c(2,2,4)= -139._R_P/ 840._R_P; c(2,2,5)=  113._R_P/ 840._R_P  ! cell 2
        c(2,1,3)=  -55._R_P/ 840._R_P; c(2,1,4)=   29._R_P/ 840._R_P; c(2,1,5)=  -27._R_P/ 840._R_P  ! cell 1
        c(2,0,3)=    5._R_P/ 840._R_P; c(2,0,4)=   -3._R_P/ 840._R_P; c(2,0,5)=    3._R_P/ 840._R_P  ! cell 0
        !  stencil  6                ;    stencil  7
        c(2,7,6)=  -15._R_P/ 840._R_P; c(2,7,7)=  105._R_P/ 840._R_P  ! cell 7
        c(2,6,6)=  225._R_P/ 840._R_P; c(2,6,7)= 1443._R_P/ 840._R_P  ! cell 6
        c(2,5,6)= 1023._R_P/ 840._R_P; c(2,5,7)=-1497._R_P/ 840._R_P  ! cell 5
        c(2,4,6)= -657._R_P/ 840._R_P; c(2,4,7)= 1443._R_P/ 840._R_P  ! cell 4
        c(2,3,6)=  393._R_P/ 840._R_P; c(2,3,7)=-1007._R_P/ 840._R_P  ! cell 3
        c(2,2,6)= -167._R_P/ 840._R_P; c(2,2,7)=  463._R_P/ 840._R_P  ! cell 2
        c(2,1,6)=   43._R_P/ 840._R_P; c(2,1,7)= -125._R_P/ 840._R_P  ! cell 1
        c(2,0,6)=   -5._R_P/ 840._R_P; c(2,0,7)=   15._R_P/ 840._R_P  ! cell 0
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        !  stencil  8                  ;    stencil  7                  ;    stencil  6
        c(1,0,8)=   280._R_P/ 2520._R_P; c(1,0,7)=   -35._R_P/ 2520._R_P; c(1,0,6)=    10._R_P/ 2520._R_P  ! cell 0
        c(1,1,8)= -2555._R_P/ 2520._R_P; c(1,1,7)=   325._R_P/ 2520._R_P; c(1,1,6)=   -95._R_P/ 2520._R_P  ! cell 1
        c(1,2,8)= 10405._R_P/ 2520._R_P; c(1,2,7)= -1355._R_P/ 2520._R_P; c(1,2,6)=   409._R_P/ 2520._R_P  ! cell 2
        c(1,3,8)=-24875._R_P/ 2520._R_P; c(1,3,7)=  3349._R_P/ 2520._R_P; c(1,3,6)= -1061._R_P/ 2520._R_P  ! cell 3
        c(1,4,8)= 38629._R_P/ 2520._R_P; c(1,4,7)= -5471._R_P/ 2520._R_P; c(1,4,6)=  1879._R_P/ 2520._R_P  ! cell 4
        c(1,5,8)=-40751._R_P/ 2520._R_P; c(1,5,7)=  6289._R_P/ 2520._R_P; c(1,5,6)= -2531._R_P/ 2520._R_P  ! cell 5
        c(1,6,8)= 29809._R_P/ 2520._R_P; c(1,6,7)= -5471._R_P/ 2520._R_P; c(1,6,6)=  3349._R_P/ 2520._R_P  ! cell 6
        c(1,7,8)=-15551._R_P/ 2520._R_P; c(1,7,7)=  4609._R_P/ 2520._R_P; c(1,7,6)=   595._R_P/ 2520._R_P  ! cell 7
        c(1,8,8)=  7129._R_P/ 2520._R_P; c(1,8,7)=   280._R_P/ 2520._R_P; c(1,8,6)=-   35._R_P/ 2520._R_P  ! cell 8
        !  stencil  5                  ;    stencil  4                  ;    stencil  3
        c(1,0,5)=    -5._R_P/ 2520._R_P; c(1,0,4)=     4._R_P/ 2520._R_P; c(1,0,3)=    -5._R_P/ 2520._R_P  ! cell 0
        c(1,1,5)=    49._R_P/ 2520._R_P; c(1,1,4)=   -41._R_P/ 2520._R_P; c(1,1,3)=    55._R_P/ 2520._R_P  ! cell 1
        c(1,2,5)=  -221._R_P/ 2520._R_P; c(1,2,4)=   199._R_P/ 2520._R_P; c(1,2,3)=  -305._R_P/ 2520._R_P  ! cell 2
        c(1,3,5)=   619._R_P/ 2520._R_P; c(1,3,4)=  -641._R_P/ 2520._R_P; c(1,3,3)=  1375._R_P/ 2520._R_P  ! cell 3
        c(1,4,5)= -1271._R_P/ 2520._R_P; c(1,4,4)=  1879._R_P/ 2520._R_P; c(1,4,3)=  1879._R_P/ 2520._R_P  ! cell 4
        c(1,5,5)=  2509._R_P/ 2520._R_P; c(1,5,4)=  1375._R_P/ 2520._R_P; c(1,5,3)=  -641._R_P/ 2520._R_P  ! cell 5
        c(1,6,5)=   955._R_P/ 2520._R_P; c(1,6,4)=  -305._R_P/ 2520._R_P; c(1,6,3)=   199._R_P/ 2520._R_P  ! cell 6
        c(1,7,5)=  -125._R_P/ 2520._R_P; c(1,7,4)=    55._R_P/ 2520._R_P; c(1,7,3)=   -41._R_P/ 2520._R_P  ! cell 7
        c(1,8,5)=    10._R_P/ 2520._R_P; c(1,8,4)=    -5._R_P/ 2520._R_P; c(1,8,3)=     4._R_P/ 2520._R_P  ! cell 8
        !  stencil  2                  ;    stencil  1                  ;    stencil  0
        c(1,0,2)=    10._R_P/ 2520._R_P; c(1,0,1)=   -35._R_P/ 2520._R_P; c(1,0,0)=   280._R_P/ 2520._R_P  ! cell 0
        c(1,1,2)=  -125._R_P/ 2520._R_P; c(1,1,1)=   595._R_P/ 2520._R_P; c(1,1,0)=  4609._R_P/ 2520._R_P  ! cell 1
        c(1,2,2)=   955._R_P/ 2520._R_P; c(1,2,1)=  3349._R_P/ 2520._R_P; c(1,2,0)= -5471._R_P/ 2520._R_P  ! cell 2
        c(1,3,2)=  2509._R_P/ 2520._R_P; c(1,3,1)= -2531._R_P/ 2520._R_P; c(1,3,0)=  6289._R_P/ 2520._R_P  ! cell 3
        c(1,4,2)= -1271._R_P/ 2520._R_P; c(1,4,1)=  1879._R_P/ 2520._R_P; c(1,4,0)= -5471._R_P/ 2520._R_P  ! cell 4
        c(1,5,2)=   619._R_P/ 2520._R_P; c(1,5,1)= -1061._R_P/ 2520._R_P; c(1,5,0)=  3349._R_P/ 2520._R_P  ! cell 5
        c(1,6,2)=  -221._R_P/ 2520._R_P; c(1,6,1)=   409._R_P/ 2520._R_P; c(1,6,0)= -1355._R_P/ 2520._R_P  ! cell 6
        c(1,7,2)=    49._R_P/ 2520._R_P; c(1,7,1)=   -95._R_P/ 2520._R_P; c(1,7,0)=   325._R_P/ 2520._R_P  ! cell 7
        c(1,8,2)=    -5._R_P/ 2520._R_P; c(1,8,1)=    10._R_P/ 2520._R_P; c(1,8,0)=   -35._R_P/ 2520._R_P  ! cell 8
        ! 2 => right interface (i+1/2)
        !  stencil  0                  ;    stencil  1                  ;    stencil  2
        c(2,8,0)=   280._R_P/ 2520._R_P; c(2,8,1)=   -35._R_P/ 2520._R_P; c(2,8,2)=    10._R_P/ 2520._R_P  ! cell 8
        c(2,7,0)= -2555._R_P/ 2520._R_P; c(2,7,1)=   325._R_P/ 2520._R_P; c(2,7,2)=   -95._R_P/ 2520._R_P  ! cell 7
        c(2,6,0)= 10405._R_P/ 2520._R_P; c(2,6,1)= -1355._R_P/ 2520._R_P; c(2,6,2)=   409._R_P/ 2520._R_P  ! cell 6
        c(2,5,0)=-24875._R_P/ 2520._R_P; c(2,5,1)=  3349._R_P/ 2520._R_P; c(2,5,2)= -1061._R_P/ 2520._R_P  ! cell 5
        c(2,4,0)= 38629._R_P/ 2520._R_P; c(2,4,1)= -5471._R_P/ 2520._R_P; c(2,4,2)=  1879._R_P/ 2520._R_P  ! cell 4
        c(2,3,0)=-40751._R_P/ 2520._R_P; c(2,3,1)=  6289._R_P/ 2520._R_P; c(2,3,2)= -2531._R_P/ 2520._R_P  ! cell 3
        c(2,2,0)= 29809._R_P/ 2520._R_P; c(2,2,1)= -5471._R_P/ 2520._R_P; c(2,2,2)=  3349._R_P/ 2520._R_P  ! cell 2
        c(2,1,0)=-15551._R_P/ 2520._R_P; c(2,1,1)=  4609._R_P/ 2520._R_P; c(2,1,2)=   595._R_P/ 2520._R_P  ! cell 1
        c(2,0,0)=  7129._R_P/ 2520._R_P; c(2,0,1)=   280._R_P/ 2520._R_P; c(2,0,2)=   -35._R_P/ 2520._R_P  ! cell 0
        !  stencil  3                  ;    stencil  4                  ;    stencil  5
        c(2,8,3)=    -5._R_P/ 2520._R_P; c(2,8,4)=     4._R_P/ 2520._R_P; c(2,8,5)=    -5._R_P/ 2520._R_P  ! cell 8
        c(2,7,3)=    49._R_P/ 2520._R_P; c(2,7,4)=   -41._R_P/ 2520._R_P; c(2,7,5)=    55._R_P/ 2520._R_P  ! cell 7
        c(2,6,3)=  -221._R_P/ 2520._R_P; c(2,6,4)=   199._R_P/ 2520._R_P; c(2,6,5)=  -305._R_P/ 2520._R_P  ! cell 6
        c(2,5,3)=   619._R_P/ 2520._R_P; c(2,5,4)=  -641._R_P/ 2520._R_P; c(2,5,5)=  1375._R_P/ 2520._R_P  ! cell 5
        c(2,4,3)= -1271._R_P/ 2520._R_P; c(2,4,4)=  1879._R_P/ 2520._R_P; c(2,4,5)=  1879._R_P/ 2520._R_P  ! cell 4
        c(2,3,3)=  2509._R_P/ 2520._R_P; c(2,3,4)=  1375._R_P/ 2520._R_P; c(2,3,5)=  -641._R_P/ 2520._R_P  ! cell 3
        c(2,2,3)=   955._R_P/ 2520._R_P; c(2,2,4)=  -305._R_P/ 2520._R_P; c(2,2,5)=   199._R_P/ 2520._R_P  ! cell 2
        c(2,1,3)=  -125._R_P/ 2520._R_P; c(2,1,4)=    55._R_P/ 2520._R_P; c(2,1,5)=   -41._R_P/ 2520._R_P  ! cell 1
        c(2,0,3)=    10._R_P/ 2520._R_P; c(2,0,4)=    -5._R_P/ 2520._R_P; c(2,0,5)=     4._R_P/ 2520._R_P  ! cell 0
        !  stencil  6                  ;    stencil  7                  ;    stencil  8
        c(2,8,6)=    10._R_P/ 2520._R_P; c(2,8,7)=   -35._R_P/ 2520._R_P; c(2,8,8)=   280._R_P/ 2520._R_P  ! cell 8
        c(2,7,6)=  -125._R_P/ 2520._R_P; c(2,7,7)=   595._R_P/ 2520._R_P; c(2,7,8)=  4609._R_P/ 2520._R_P  ! cell 7
        c(2,6,6)=   955._R_P/ 2520._R_P; c(2,6,7)=  3349._R_P/ 2520._R_P; c(2,6,8)= -5471._R_P/ 2520._R_P  ! cell 6
        c(2,5,6)=  2509._R_P/ 2520._R_P; c(2,5,7)= -2531._R_P/ 2520._R_P; c(2,5,8)=  6289._R_P/ 2520._R_P  ! cell 5
        c(2,4,6)= -1271._R_P/ 2520._R_P; c(2,4,7)=  1879._R_P/ 2520._R_P; c(2,4,8)= -5471._R_P/ 2520._R_P  ! cell 4
        c(2,3,6)=   619._R_P/ 2520._R_P; c(2,3,7)= -1061._R_P/ 2520._R_P; c(2,3,8)=  3349._R_P/ 2520._R_P  ! cell 3
        c(2,2,6)=  -221._R_P/ 2520._R_P; c(2,2,7)=   409._R_P/ 2520._R_P; c(2,2,8)= -1355._R_P/ 2520._R_P  ! cell 2
        c(2,1,6)=    49._R_P/ 2520._R_P; c(2,1,7)=   -95._R_P/ 2520._R_P; c(2,1,8)=   325._R_P/ 2520._R_P  ! cell 1
        c(2,0,6)=    -5._R_P/ 2520._R_P; c(2,0,7)=    10._R_P/ 2520._R_P; c(2,0,8)=   -35._R_P/ 2520._R_P  ! cell 0
      endselect
  endassociate
  endsubroutine create

  pure subroutine compute_int(self, stencil, values)
  !< Compute interpolations (interpolate).
  class(interpolations_rec_js), intent(in)  :: self               !< Interpolations.
  real(R_P),                    intent(in)  :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(R_P),                    intent(out) :: values(0:)         !< Interpolations values.
  ! empty procedure
  endsubroutine compute_int

  pure subroutine compute_rec(self, stencil, values)
  !< Compute interpolations (reconstruct).
  class(interpolations_rec_js), intent(in)  :: self                  !< Interpolations.
  real(R_P),                    intent(in)  :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(R_P),                    intent(out) :: values(1:, 0:)        !< Interpolations values.
  integer(I_P)                              :: s1                    !< Counter.
  integer(I_P)                              :: s2                    !< Counter.
  integer(I_P)                              :: f                     !< Counter.

  values = 0._R_P
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
  string = string//prefix_//'  - S   = '//trim(str(self%S))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy interpolations.
  class(interpolations_rec_js), intent(inout) :: self !< Interpolations.

  call self%destroy_
  if (allocated(self%coef)) deallocate(self%coef)
  endsubroutine destroy
endmodule wenoof_interpolations_rec_js
