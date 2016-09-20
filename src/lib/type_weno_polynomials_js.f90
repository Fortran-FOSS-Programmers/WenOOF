module type_weno_polynomials_js
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Lagrange polynomials for Jiang-Shu WENO schemes.
!<
!< @note The provided polynomials implement the Lagrange polynomials defined in *Efficient Implementation
!< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_polynomials_js
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_polynomials) :: weno_polynomials_js
!< Lagrange polynomials for Jiang-Shu WENO schemes object.
!<
!< @note The provided polynomials implement the Lagrange polynomials defined in *Efficient Implementation
!< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130
  private
  real(R_P), allocatable :: coef(:,:)   !< Polynomial coefficients [1:2,0:S-1,0:S-1].
  contains
    procedure, pass(self), public :: destroy
    procedure, pass(self), public :: create
    procedure, pass(self), public :: description
    procedure, pass(self), public :: compute
endtype weno_polynomials_js
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! deferred public methods
  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy Jiang-Shu polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_polynomials_js), intent(inout) :: self   !< WENO smoothenss indicators.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%coef)) deallocate(self%coef)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure subroutine create(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO polynomials coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_polynomials_js), intent(inout) :: self       !< WENO smoothness indicators.
  integer(I_P),               intent(in)    :: S          !< Number of stencils used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  associate(c => self%coef)
    allocate(c(1:2, 0:S - 1, 0:S - 1))
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
        c(1,4,0)= -13._R_P/60._R_P; c(1,5,0)=   1._R_P/6._R_P   ! stencil 0
        c(1,4,1)=   7._R_P/60._R_P; c(1,5,1)=  -1._R_P/60._R_P  ! stencil 1
        c(1,4,2)=  -2._R_P/15._R_P; c(1,5,2)=   1._R_P/60._R_P  ! stencil 2
        c(1,4,3)=  11._R_P/30._R_P; c(1,5,3)=  -1._R_P/30._R_P  ! stencil 3
        c(1,4,4)=  29._R_P/20._R_P; c(1,5,4)=   1._R_P/6._R_P   ! stencil 4
        c(1,4,5)= -71._R_P/20._R_P; c(1,5,5)=  49._R_P/20._R_P  ! stencil 5
        ! 2 => right interface (i+1/2)
        !  cell  0                ;    cell  1                ;   cell  2                 ;    cell  3
        c(2,0,0)=  49._R_P/20._R_P; c(2,1,0)= -71._R_P/20._R_P; c(2,2,0)=  79._R_P/20._R_P; c(2,3,0)=-163._R_P/60._R_P  ! stencil 0
        c(1,0,1)=   1._R_P/6._R_P ; c(1,1,1)=  29._R_P/20._R_P; c(1,2,1)= -21._R_P/20._R_P; c(1,3,1)=  37._R_P/60._R_P  ! stencil 1
        c(1,0,2)=  -1._R_P/30._R_P; c(1,1,2)=  11._R_P/30._R_P; c(1,2,2)=  19._R_P/20._R_P; c(1,3,2)= -23._R_P/60._R_P  ! stencil 2
        c(1,0,3)=   1._R_P/60._R_P; c(1,1,3)=  -2._R_P/15._R_P; c(1,2,3)=  37._R_P/60._R_P; c(1,3,3)=  37._R_P/60._R_P  ! stencil 3
        c(1,0,4)=  -1._R_P/60._R_P; c(1,1,4)=   7._R_P/60._R_P; c(1,2,4)= -23._R_P/60._R_P; c(1,3,4)=  19._R_P/20._R_P  ! stencil 4
        c(1,0,5)=   1._R_P/30._R_P; c(1,1,5)= -13._R_P/60._R_P; c(1,2,5)=  37._R_P/60._R_P; c(1,3,5)= -21._R_P/20._R_P  ! stencil 5
        !  cell  4                ;    cell  5
        c(2,4,0)=  31._R_P/5._R_P ; c(2,5,0)=  -1._R_P/6._R_P   ! stencil 0
        c(1,4,1)= -13._R_P/60._R_P; c(1,5,1)=   1._R_P/6._R_P   ! stencil 1
        c(1,4,2)=   7._R_P/60._R_P; c(1,5,2)=  -1._R_P/60._R_P  ! stencil 2
        c(1,4,3)=  -2._R_P/15._R_P; c(1,5,3)=   1._R_P/60._R_P  ! stencil 3
        c(1,4,4)=  11._R_P/30._R_P; c(1,5,4)=  -1._R_P/30._R_P  ! stencil 4
        c(1,4,5)=  29._R_P/20._R_P; c(1,5,5)=   1._R_P/6._R_P   ! stencil 5
      case(7) ! 13th order
        ! 1 => left interface (i-1/2)
        !  cell  0                  ;    cell  1                  ;    cell  2                  ;    cell  3
        c(1,0,0)=    1._R_P/7._R_P  ; c(1,1,0)=  223._R_P/140._R_P; c(1,2,0)= -197._R_P/140._R_P; c(1,3,0)=  153._R_P/140._R_P  ! stencil 0
        c(1,0,1)=   -1._R_P/42._R_P ; c(1,1,1)=   13._R_P/42._R_P ; c(1,2,1)=  153._R_P/140._R_P; c(1,3,1)= -241._R_P/420._R_P  ! stencil 1
        c(1,0,2)=    1._R_P/105._R_P; c(1,1,2)=  -19._R_P/210._R_P; c(1,2,2)=  107._R_P/210._R_P; c(1,3,2)=  319._R_P/420._R_P  ! stencil 2
        c(1,0,3)=   -1._R_P/140._R_P; c(1,1,3)=    5._R_P/84._R_P ; c(1,2,3)= -101._R_P/420._R_P; c(1,3,3)=  319._R_P/420._R_P  ! stencil 3
        c(1,0,4)=    1._R_P/105._R_P; c(1,1,4)=  -31._R_P/420._R_P; c(1,2,4)=  109._R_P/420._R_P; c(1,3,4)= -241._R_P/420._R_P  ! stencil 4
        c(1,0,5)=   -1._R_P/42._R_P ; c(1,1,5)=  -37._R_P/210._R_P; c(1,2,5)= -241._R_P/420._R_P; c(1,3,5)=  153._R_P/140._R_P  ! stencil 5
        c(1,0,6)=   -1._R_P/7._R_P  ; c(1,1,6)=  -43._R_P/42._R_P ; c(1,2,6)=  667._R_P/210._R_P; c(1,3,6)=-2341._R_P/420._R_P  ! stencil 6
        !  cell  4                  ;    cell  5                  ;    cell  6
        c(1,4,0)= -241._R_P/420._R_P; c(1,5,0)=   37._R_P/210._R_P; c(1,6,0)=   -1._R_P/42._R_P   ! stencil 0
        c(1,4,1)=  109._R_P/420._R_P; c(1,5,1)=  -31._R_P/420._R_P; c(1,6,1)=    1._R_P/105._R_P  ! stencil 1
        c(1,4,2)= -101._R_P/420._R_P; c(1,5,2)=    5._R_P/84._R_P ; c(1,6,2)=   -1._R_P/140._R_P  ! stencil 2
        c(1,4,3)=  107._R_P/210._R_P; c(1,5,3)=  -19._R_P/210._R_P; c(1,6,3)=    1._R_P/105._R_P  ! stencil 3
        c(1,4,4)=  153._R_P/140._R_P; c(1,5,4)=   13._R_P/42._R_P ; c(1,6,4)=   -1._R_P/42._R_P   ! stencil 4
        c(1,4,5)= -197._R_P/140._R_P; c(1,5,5)=  223._R_P/140._R_P; c(1,6,5)=    1._R_P/7._R_P    ! stencil 5
        c(1,4,6)=  853._R_P/140._R_P; c(1,5,6)= -617._R_P/140._R_P; c(1,6,6)=  363._R_P/140._R_P  ! stencil 6
        ! 2 => right interface (i+1/2)
        !  cell  0                  ;    cell  1                  ;    cell  2                  ;    cell  3
        c(1,0,6)=  363._R_P/140._R_P; c(1,1,6)= -617._R_P/140._R_P; c(1,2,6)=  853._R_P/140._R_P; c(1,3,6)=-2341._R_P/420._R_P  ! stencil 0
        c(1,0,0)=    1._R_P/7._R_P  ; c(1,1,0)=  223._R_P/140._R_P; c(1,2,0)= -197._R_P/140._R_P; c(1,3,0)=  153._R_P/140._R_P  ! stencil 1
        c(1,0,1)=   -1._R_P/42._R_P ; c(1,1,1)=   13._R_P/42._R_P ; c(1,2,1)=  153._R_P/140._R_P; c(1,3,1)= -241._R_P/420._R_P  ! stencil 2
        c(1,0,2)=    1._R_P/105._R_P; c(1,1,2)=  -19._R_P/210._R_P; c(1,2,2)=  107._R_P/210._R_P; c(1,3,2)=  319._R_P/420._R_P  ! stencil 3
        c(1,0,3)=   -1._R_P/140._R_P; c(1,1,3)=    5._R_P/84._R_P ; c(1,2,3)= -101._R_P/420._R_P; c(1,3,3)=  319._R_P/420._R_P  ! stencil 4
        c(1,0,4)=    1._R_P/105._R_P; c(1,1,4)=  -31._R_P/420._R_P; c(1,2,4)=  109._R_P/420._R_P; c(1,3,4)= -241._R_P/420._R_P  ! stencil 5
        c(1,0,5)=   -1._R_P/42._R_P ; c(1,1,5)=  -37._R_P/210._R_P; c(1,2,5)= -241._R_P/420._R_P; c(1,3,5)=  153._R_P/140._R_P  ! stencil 6
        !  cell  4                  ;    cell  5                  ;    cell  6
        c(1,4,6)=  667._R_P/210._R_P; c(1,5,6)=  -43._R_P/42._R_P ; c(1,6,6)=    1._R_P/7._R_P    ! stencil 0
        c(1,4,0)= -241._R_P/420._R_P; c(1,5,0)=   37._R_P/210._R_P; c(1,6,0)=   -1._R_P/42._R_P   ! stencil 1
        c(1,4,1)=  109._R_P/420._R_P; c(1,5,1)=  -31._R_P/420._R_P; c(1,6,1)=    1._R_P/105._R_P  ! stencil 2
        c(1,4,2)= -101._R_P/420._R_P; c(1,5,2)=    5._R_P/84._R_P ; c(1,6,2)=   -1._R_P/140._R_P  ! stencil 3
        c(1,4,3)=  107._R_P/210._R_P; c(1,5,3)=  -19._R_P/210._R_P; c(1,6,3)=    1._R_P/105._R_P  ! stencil 4
        c(1,4,4)=  153._R_P/140._R_P; c(1,5,4)=   13._R_P/42._R_P ; c(1,6,4)=   -1._R_P/42._R_P   ! stencil 5
        c(1,4,5)= -197._R_P/140._R_P; c(1,5,5)=  223._R_P/140._R_P; c(1,6,5)=    1._R_P/7._R_P    ! stencil 6
      endselect
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create

  pure subroutine description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_polynomials_js),    intent(in)  :: self   !< WENO polynomial.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  character(len=1), parameter                :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO polynomial'//nl
  string = string//'  Based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    poly(poly_coef,v)'//nl
  string = string//'  where:'//nl
  string = string//'    poly_coef: real(R_P), intent(IN), the polynomial coefficient of the value'//nl
  string = string//'    v: real(R_P), intent(IN), the selected value from the stencil'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure function compute(self, poly_coef, v) result(poly)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the partial value of the interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_polynomials_js), intent(in) :: self       !< WENO alpha coefficient.
  real(R_P),                  intent(in) :: poly_coef  !< Coefficient of the smoothness indicator.
  real(R_P),                  intent(IN) :: v          !< Selected value from the stencil used for the interpolation.
  real(R_P),                             :: poly       !< Partial value of the interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  poly = poly_coef * v
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute

!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_polynomials_js
