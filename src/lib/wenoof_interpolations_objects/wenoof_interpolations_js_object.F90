!< Jiang-Shu (Lagrange) polynomials object for derivative reconstruction.
module wenoof_polynomials_js
!< Jiang-Shu (Lagrange) polynomials object for derivative reconstruction.
!<
!< @note The provided polynomials implement the Lagrange polynomials defined in *Efficient Implementation
!< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P
use wenoof_base_object
use wenoof_interpolations_object

implicit none
private
public :: interpolations_js_object
public :: interpolations_js_object_constructor
public :: create_interpolations_js_object_constructor

type, extends(interpolations_object_constructor) :: interpolations_js_object_constructor
  !< Jiang-Shu (Lagrange) polynomials object constructor.
  logical :: face_left=.true.  !< Activate left-face interpolations computation.
  logical :: face_right=.true. !< Activate right-face interpolations computation.
endtype interpolations_js_object_constructor

type, extends(interpolations_object) :: interpolations_js_object
  !< Jiang-Shu (Lagrange) polynomials object for derivative reconstruction.
  !<
  !< @note The provided polynomials implement the Lagrange polynomials defined in *Efficient Implementation
  !< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  private
  integer(I_P)           :: f1=1_I_P    !< Lower bound of faces index.
  integer(I_P)           :: f2=2_I_P    !< Upper bound of faces index.
  integer(I_P)           :: ff=0_I_P    !< Offset (step) of faces index.
  real(R_P), allocatable :: coef(:,:,:) !< Polynomial coefficients [1:2,0:S-1,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute interpolations.
    procedure, pass(self) :: description !< Return interpolations string-description.
    ! overridden public methods
    procedure, pass(self) :: create  !< Create interpolations.
    procedure, pass(self) :: destroy !< Destroy interpolations.
endtype interpolations_js_object

contains
  ! public non TBP
  subroutine create_interpolations_js_object_constructor(S, face_left, face_right, constructor)
  !< Create interpolations constructor.
  integer(I_P),                                          intent(in)  :: S           !< Stencils dimension.
  logical,                                               intent(in)  :: face_left   !< Activate left-face interpolations.
  logical,                                               intent(in)  :: face_right  !< Activate right-face interpolations.
  class(interpolations_object_constructor), allocatable, intent(out) :: constructor !< Interpolations constructor.

  allocate(interpolations_js_object_constructor :: constructor)
  constructor%S = S
  constructor%face_left = face_left
  constructor%face_right = face_right
  endsubroutine create_interpolations_js_object_constructor

  ! deferred public methods
  pure subroutine compute(self, stencil)
  !< Compute polynomials.
  class(interpolations_js_object), intent(inout) :: self                  !< Interpolations.
  real(R_P),                       intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P)                                   :: s1                    !< Counter.
  integer(I_P)                                   :: s2                    !< Counter.
  integer(I_P)                                   :: f                     !< Counter.

  self%values = 0._R_P
  do s1=0, self%S - 1 ! stencils loop
    do s2=0, self%S - 1 ! values loop
      do f=self%f1, self%f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        self%values(f, s1) = self%values(f, s1) + self%coef(f, s2, s1) * stencil(f + self%ff, -s2 + s1)
      enddo
    enddo
  enddo
  endsubroutine compute

  pure function description(self) result(string)
  !< Return polynomials string-description.
  class(interpolations_js_object), intent(in) :: self             !< Interpolations.
  character(len=:), allocatable               :: string           !< String-description.
  character(len=1), parameter                 :: nl=new_line('a') !< New line character.

  string = 'WENO polynomial'//nl
  string = string//'  Based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and'//nl
  string = string//'  on the work by Gerolimos, Sénéchal and  Vallet  "Very-High-Order WENO Schemes", see '// &
           'JCP, 2009, vol. 228, pp. 8481-8524, doi:10.1016/j.jcp.2009.07.039'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    compute(S, stencil, f1, f2, ff)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(in), the number of the stencils used'//nl
  string = string//'    stencil: real(R_P), intent(IN), the stencil used for the interpolation [1:2, 1-S:-1+S]'//nl
  string = string//'    f1, f2: integer(I_P), intent(in), the faces to be computed (1 => left interface, 2 => right interface)'//nl
  string = string//'    ff: integer(I_P), intent(in), the parameter for the stencil value choice'
  endfunction description

  ! overridden public methods
  pure subroutine create(self, constructor)
  !< Create interpolations.
  class(interpolations_js_object), intent(inout) :: self        !< Interpolations.
  class(base_object_constructor),  intent(in)    :: constructor !< Interpolations constructor.
  integer(I_P)                                   :: S           !< Stencils dimension.

  call self%destroy
  call self%interpolations_object%create(constructor=constructor)
  select type(constructor)
  class is(polynomials_js_constructor)
    S = constructor%S
    if (constructor%face_left.and.constructor%face_right) then
      self%f1 = 1_I_P; self%f2 = 2_I_P; self%ff = 0_I_P
    elseif (constructor%face_left) then
      self%f1 = 1_I_P; self%f2 = 1_I_P; self%ff = 0_I_P
    elseif (constructor%face_right) then
      self%f1 = 2_I_P; self%f2 = 2_I_P; self%ff = -1_I_P
    endif
    allocate(self%coef(1:2, 0:S - 1, 0:S - 1))
  class default
    ! @TODO add error handling
  endselect

  associate(c => self%coef)
    select case(S)
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
        !  cell  0                  ;    cell  1                  ;    cell  2
        c(1,0,0)=    1._R_P/7._R_P  ; c(1,1,0)=  223._R_P/140._R_P; c(1,2,0)= -197._R_P/140._R_P  ! stencil 0
        c(1,0,1)=   -1._R_P/42._R_P ; c(1,1,1)=   13._R_P/42._R_P ; c(1,2,1)=  153._R_P/140._R_P  ! stencil 1
        c(1,0,2)=    1._R_P/105._R_P; c(1,1,2)=  -19._R_P/210._R_P; c(1,2,2)=  107._R_P/210._R_P  ! stencil 2
        c(1,0,3)=   -1._R_P/140._R_P; c(1,1,3)=    5._R_P/84._R_P ; c(1,2,3)= -101._R_P/420._R_P  ! stencil 3
        c(1,0,4)=    1._R_P/105._R_P; c(1,1,4)=  -31._R_P/420._R_P; c(1,2,4)=  109._R_P/420._R_P  ! stencil 4
        c(1,0,5)=   -1._R_P/42._R_P ; c(1,1,5)=  -37._R_P/210._R_P; c(1,2,5)= -241._R_P/420._R_P  ! stencil 5
        c(1,0,6)=   -1._R_P/7._R_P  ; c(1,1,6)=  -43._R_P/42._R_P ; c(1,2,6)=  667._R_P/210._R_P  ! stencil 6
        !  cell  3                       cell  4                  ;    cell  5
        c(1,3,0)=  153._R_P/140._R_P; c(1,4,0)= -241._R_P/420._R_P; c(1,5,0)=   37._R_P/210._R_P  ! stencil 0
        c(1,3,1)= -241._R_P/420._R_P; c(1,4,1)=  109._R_P/420._R_P; c(1,5,1)=  -31._R_P/420._R_P  ! stencil 1
        c(1,3,2)=  319._R_P/420._R_P; c(1,4,2)= -101._R_P/420._R_P; c(1,5,2)=    5._R_P/84._R_P   ! stencil 2
        c(1,3,3)=  319._R_P/420._R_P; c(1,4,3)=  107._R_P/210._R_P; c(1,5,3)=  -19._R_P/210._R_P  ! stencil 3
        c(1,3,4)= -241._R_P/420._R_P; c(1,4,4)=  153._R_P/140._R_P; c(1,5,4)=   13._R_P/42._R_P   ! stencil 4
        c(1,3,5)=  153._R_P/140._R_P; c(1,4,5)= -197._R_P/140._R_P; c(1,5,5)=  223._R_P/140._R_P  ! stencil 5
        c(1,3,6)=-2341._R_P/420._R_P; c(1,4,6)=  853._R_P/140._R_P; c(1,5,6)= -617._R_P/140._R_P  ! stencil 6
        !  cell  6
        c(1,6,0)=   -1._R_P/42._R_P   ! stencil 0
        c(1,6,1)=    1._R_P/105._R_P  ! stencil 1
        c(1,6,2)=   -1._R_P/140._R_P  ! stencil 2
        c(1,6,3)=    1._R_P/105._R_P  ! stencil 3
        c(1,6,4)=   -1._R_P/42._R_P   ! stencil 4
        c(1,6,5)=    1._R_P/7._R_P    ! stencil 5
        c(1,6,6)=  363._R_P/140._R_P  ! stencil 6
        ! 2 => right interface (i+1/2)
        !  cell  0                  ;    cell  1                  ;    cell  2
        c(2,0,0)=  363._R_P/140._R_P; c(2,1,0)= -617._R_P/140._R_P; c(2,2,0)=  853._R_P/140._R_P  ! stencil 0
        c(2,0,1)=    1._R_P/7._R_P  ; c(2,1,1)=  223._R_P/140._R_P; c(2,2,1)= -197._R_P/140._R_P  ! stencil 1
        c(2,0,2)=   -1._R_P/42._R_P ; c(2,1,2)=   13._R_P/42._R_P ; c(2,2,2)=  153._R_P/140._R_P  ! stencil 2
        c(2,0,3)=    1._R_P/105._R_P; c(2,1,3)=  -19._R_P/210._R_P; c(2,2,3)=  107._R_P/210._R_P  ! stencil 3
        c(2,0,4)=   -1._R_P/140._R_P; c(2,1,4)=    5._R_P/84._R_P ; c(2,2,4)= -101._R_P/420._R_P  ! stencil 4
        c(2,0,5)=    1._R_P/105._R_P; c(2,1,5)=  -31._R_P/420._R_P; c(2,2,5)=  109._R_P/420._R_P  ! stencil 5
        c(2,0,6)=   -1._R_P/42._R_P ; c(2,1,6)=   37._R_P/210._R_P; c(2,2,6)= -241._R_P/420._R_P  ! stencil 6
        !  cell  3                  ;    cell  4                  ;    cell  5
        c(2,3,0)=-2341._R_P/420._R_P; c(2,4,0)=  667._R_P/210._R_P; c(2,5,0)=  -43._R_P/42._R_P   ! stencil 0
        c(2,3,1)=  153._R_P/140._R_P; c(2,4,1)= -241._R_P/420._R_P; c(2,5,1)=   37._R_P/210._R_P  ! stencil 1
        c(2,3,2)= -241._R_P/420._R_P; c(2,4,2)=  109._R_P/420._R_P; c(2,5,2)=  -31._R_P/420._R_P  ! stencil 2
        c(2,3,3)=  319._R_P/420._R_P; c(2,4,3)= -101._R_P/420._R_P; c(2,5,3)=    5._R_P/84._R_P   ! stencil 3
        c(2,3,4)=  319._R_P/420._R_P; c(2,4,4)=  107._R_P/210._R_P; c(2,5,4)=  -19._R_P/210._R_P  ! stencil 4
        c(2,3,5)= -241._R_P/420._R_P; c(2,4,5)=  153._R_P/140._R_P; c(2,5,5)=   13._R_P/42._R_P   ! stencil 5
        c(2,3,6)=  153._R_P/140._R_P; c(2,4,6)= -197._R_P/140._R_P; c(2,5,6)=  223._R_P/140._R_P  ! stencil 6
        !  cell  6
        c(2,6,0)=    1._R_P/7._R_P    ! stencil 0
        c(2,6,1)=   -1._R_P/42._R_P   ! stencil 1
        c(2,6,2)=    1._R_P/105._R_P  ! stencil 2
        c(2,6,3)=   -1._R_P/140._R_P  ! stencil 3
        c(2,6,4)=    1._R_P/105._R_P  ! stencil 4
        c(2,6,5)=   -1._R_P/42._R_P   ! stencil 5
        c(2,6,6)=    1._R_P/7._R_P    ! stencil 6
      case(8) ! 15th order
        ! 1 => left interface (i-1/2)
        !  cell  0                  ;    cell  1                  ;    cell  2
        c(1,0,0)=    7._R_P/56._R_P ; c(1,1,0)= 1443._R_P/840._R_P; c(1,2,0)=-1497._R_P/840._R_P  ! stencil 0
        c(1,0,1)=   -1._R_P/56._R_P ; c(1,1,1)=   15._R_P/56._R_P ; c(1,2,1)= 1023._R_P/840._R_P  ! stencil 1
        c(1,0,2)=    1._R_P/168._R_P; c(1,1,2)=  -11._R_P/168._R_P; c(1,2,2)=   73._R_P/168._R_P  ! stencil 2
        c(1,0,3)=   -3._R_P/840._R_P; c(1,1,3)=   29._R_P/840._R_P; c(1,2,3)= -139._R_P/840._R_P  ! stencil 3
        c(1,0,4)=    3._R_P/840._R_P; c(1,1,4)=  -27._R_P/840._R_P; c(1,2,4)=  113._R_P/840._R_P  ! stencil 4
        c(1,0,5)=   -1._R_P/168._R_P; c(1,1,5)=   43._R_P/840._R_P; c(1,2,5)= -167._R_P/840._R_P  ! stencil 5
        c(1,0,6)=    1._R_P/56._R_P ; c(1,1,6)=  -25._R_P/168._R_P; c(1,2,6)=  463._R_P/840._R_P  ! stencil 6
        c(1,0,7)=   -7._R_P/56._R_P ; c(1,1,7)=   57._R_P/56._R_P ; c(1,2,7)= -613._R_P/168._R_P  ! stencil 7
        !  cell  3                  ;    cell  4                  ;    cell  5
        c(1,3,0)= 1443._R_P/840._R_P; c(1,4,0)=-1007._R_P/840._R_P; c(1,5,0)=  463._R_P/840._R_P  ! stencil 0
        c(1,3,1)= -657._R_P/840._R_P; c(1,4,1)=  393._R_P/840._R_P; c(1,5,1)= -167._R_P/840._R_P  ! stencil 1
        c(1,3,2)=  743._R_P/840._R_P; c(1,4,2)= -307._R_P/840._R_P; c(1,5,2)=  113._R_P/840._R_P  ! stencil 2
        c(1,3,3)=  533._R_P/840._R_P; c(1,4,3)=  533._R_P/840._R_P; c(1,5,3)= -139._R_P/840._R_P  ! stencil 3
        c(1,3,4)= -307._R_P/840._R_P; c(1,4,4)=  743._R_P/840._R_P; c(1,5,4)=   73._R_P/168._R_P  ! stencil 4
        c(1,3,5)=  393._R_P/840._R_P; c(1,4,5)= -657._R_P/840._R_P; c(1,5,5)= 1023._R_P/840._R_P  ! stencil 5
        c(1,3,6)=-1007._R_P/840._R_P; c(1,4,6)= 1443._R_P/840._R_P; c(1,5,6)=-1497._R_P/840._R_P  ! stencil 6
        c(1,3,7)= 6343._R_P/840._R_P; c(1,4,7)=-8357._R_P/840._R_P; c(1,5,7)= 7323._R_P/840._R_P  ! stencil 7
        !  cell  6                  ;    cell  7
        c(1,6,0)=  -25._R_P/168._R_P; c(1,7,0)=    1._R_P/56._R_P   ! stencil 0
        c(1,6,1)=   43._R_P/840._R_P; c(1,7,1)=   -1._R_P/168._R_P  ! stencil 1
        c(1,6,2)=  -27._R_P/840._R_P; c(1,7,2)=    3._R_P/840._R_P  ! stencil 2
        c(1,6,3)=   29._R_P/840._R_P; c(1,7,3)=   -3._R_P/840._R_P  ! stencil 3
        c(1,6,4)=  -11._R_P/168._R_P; c(1,7,4)=    1._R_P/168._R_P  ! stencil 4
        c(1,6,5)=   15._R_P/56._R_P ; c(1,7,5)=   -1._R_P/56._R_P   ! stencil 5
        c(1,6,6)= 1443._R_P/840._R_P; c(1,7,6)=    7._R_P/56._R_P   ! stencil 6
        c(1,6,7)=-4437._R_P/840._R_P; c(1,7,7)= 2283._R_P/840._R_P  ! stencil 7
        ! 2 => right interface (i+1/2)
        !  cell  0                  ;    cell  1                  ;    cell  2
        c(2,0,0)= 2283._R_P/840._R_P; c(2,1,0)=-4437._R_P/840._R_P; c(2,2,0)= 7323._R_P/840._R_P  ! stencil 0
        c(2,0,1)=    7._R_P/56._R_P ; c(2,1,1)= 1443._R_P/840._R_P; c(2,2,1)=-1497._R_P/840._R_P  ! stencil 1
        c(2,0,2)=   -1._R_P/56._R_P ; c(2,1,2)=   15._R_P/56._R_P ; c(2,2,2)= 1023._R_P/840._R_P  ! stencil 2
        c(2,0,3)=    1._R_P/168._R_P; c(2,1,3)=  -11._R_P/168._R_P; c(2,2,3)=   73._R_P/168._R_P  ! stencil 3
        c(2,0,4)=   -3._R_P/840._R_P; c(2,1,4)=   29._R_P/840._R_P; c(2,2,4)= -139._R_P/840._R_P  ! stencil 4
        c(2,0,5)=    3._R_P/840._R_P; c(2,1,5)=  -27._R_P/840._R_P; c(2,2,5)=  113._R_P/840._R_P  ! stencil 5
        c(2,0,6)=   -1._R_P/168._R_P; c(2,1,6)=   43._R_P/840._R_P; c(2,2,6)= -167._R_P/840._R_P  ! stencil 6
        c(2,0,7)=    1._R_P/56._R_P ; c(2,1,7)=  -25._R_P/168._R_P; c(2,2,7)=  463._R_P/840._R_P  ! stencil 7
        !  cell  3                  ;    cell  4                  ;    cell  5
        c(2,3,0)=-8357._R_P/840._R_P; c(2,4,0)= 6343._R_P/840._R_P; c(2,5,0)= -613._R_P/168._R_P  ! stencil 0
        c(2,3,1)= 1443._R_P/840._R_P; c(2,4,1)=-1007._R_P/840._R_P; c(2,5,1)=  463._R_P/840._R_P  ! stencil 1
        c(2,3,2)= -657._R_P/840._R_P; c(2,4,2)=  393._R_P/840._R_P; c(2,5,2)= -167._R_P/840._R_P  ! stencil 2
        c(2,3,3)=  743._R_P/840._R_P; c(2,4,3)= -307._R_P/840._R_P; c(2,5,3)=  113._R_P/840._R_P  ! stencil 3
        c(2,3,4)=  533._R_P/840._R_P; c(2,4,4)=  533._R_P/840._R_P; c(2,5,4)= -139._R_P/840._R_P  ! stencil 4
        c(2,3,5)= -307._R_P/840._R_P; c(2,4,5)=  743._R_P/840._R_P; c(2,5,5)=   73._R_P/168._R_P  ! stencil 5
        c(2,3,6)=  393._R_P/840._R_P; c(2,4,6)= -657._R_P/840._R_P; c(2,5,6)= 1023._R_P/840._R_P  ! stencil 6
        c(2,3,7)=-1007._R_P/840._R_P; c(2,4,7)= 1443._R_P/840._R_P; c(2,5,7)=-1497._R_P/840._R_P  ! stencil 7
        !  cell  6                  ;    cell  7
        c(2,6,0)=   57._R_P/56._R_P ; c(2,7,0)=   -7._R_P/56._R_P   ! stencil 0
        c(2,6,1)=  -25._R_P/168._R_P; c(2,7,1)=    1._R_P/56._R_P   ! stencil 1
        c(2,6,2)=   43._R_P/840._R_P; c(2,7,2)=   -1._R_P/168._R_P  ! stencil 2
        c(2,6,3)=  -27._R_P/840._R_P; c(2,7,3)=    3._R_P/840._R_P  ! stencil 3
        c(2,6,4)=   29._R_P/840._R_P; c(2,7,4)=   -3._R_P/840._R_P  ! stencil 4
        c(2,6,5)=  -11._R_P/168._R_P; c(2,7,5)=    1._R_P/168._R_P  ! stencil 5
        c(2,6,6)=   15._R_P/56._R_P ; c(2,7,6)=   -1._R_P/56._R_P   ! stencil 6
        c(2,6,7)= 1443._R_P/840._R_P; c(2,7,7)=    7._R_P/56._R_P   ! stencil 7
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        !  cell  0                    ;     cell  1                   ;     cell  2
        c(1,0,0)=     1._R_P/9._R_P   ; c(1,1,0)=  4609._R_P/2520._R_P; c(1,2,0)= -5471._R_P/2520._R_P  ! stencil 0
        c(1,0,1)=    -7._R_P/504._R_P ; c(1,1,1)=   119._R_P/508._R_P ; c(1,2,1)=  3349._R_P/2520._R_P  ! stencil 1
        c(1,0,2)=     1._R_P/252._R_P ; c(1,1,2)=   -25._R_P/508._R_P ; c(1,2,2)=   191._R_P/504._R_P   ! stencil 2
        c(1,0,3)=    -1._R_P/504._R_P ; c(1,1,3)=    11._R_P/508._R_P ; c(1,2,3)=   -61._R_P/504._R_P   ! stencil 3
        c(1,0,4)=     1._R_P/630._R_P ; c(1,1,4)=   -41._R_P/2520._R_P; c(1,2,4)=   199._R_P/2520._R_P  ! stencil 4
        c(1,0,5)=    -1._R_P/504._R_P ; c(1,1,5)=    49._R_P/2520._R_P; c(1,2,5)=  -221._R_P/2520._R_P  ! stencil 5
        c(1,0,6)=     1._R_P/252._R_P ; c(1,1,6)=   -19._R_P/508._R_P ; c(1,2,6)=   409._R_P/2520._R_P  ! stencil 6
        c(1,0,7)=    -7._R_P/504._R_P ; c(1,1,7)=    65._R_P/508._R_P ; c(1,2,7)=  -271._R_P/504._R_P   ! stencil 7
        c(1,0,8)=     1._R_P/9._R_P   ; c(1,1,8)=  -511._R_P/508._R_P ; c(1,2,8)=  2081._R_P/504._R_P   ! stencil 8
        !  cell  3                    ;    cell  4                    ;     cell  5
        c(1,3,0)=  6289._R_P/2520._R_P; c(1,4,0)= -5471._R_P/2520._R_P; c(1,5,0)=  3349._R_P/2520._R_P  ! stencil 0
        c(1,3,1)= -2531._R_P/2520._R_P; c(1,4,1)=  1879._R_P/2520._R_P; c(1,5,1)= -1061._R_P/2520._R_P  ! stencil 1
        c(1,3,2)=  2509._R_P/2520._R_P; c(1,4,2)= -1271._R_P/2520._R_P; c(1,5,2)=   619._R_P/2520._R_P  ! stencil 2
        c(1,3,3)=   275._R_P/504._R_P ; c(1,4,3)=  1879._R_P/2520._R_P; c(1,5,3)=  -641._R_P/2520._R_P  ! stencil 3
        c(1,3,4)=  -641._R_P/2520._R_P; c(1,4,4)=  1879._R_P/2520._R_P; c(1,5,4)=   275._R_P/504._R_P   ! stencil 4
        c(1,3,5)=   619._R_P/2520._R_P; c(1,4,5)= -1271._R_P/2520._R_P; c(1,5,5)=  2509._R_P/2520._R_P  ! stencil 5
        c(1,3,6)= -1061._R_P/2520._R_P; c(1,4,6)=  1879._R_P/2520._R_P; c(1,5,6)= -2531._R_P/2520._R_P  ! stencil 6
        c(1,3,7)=  3349._R_P/2520._R_P; c(1,4,7)= -5471._R_P/2520._R_P; c(1,5,7)=  6289._R_P/2520._R_P  ! stencil 7
        c(1,3,8)= -4975._R_P/504._R_P ; c(1,4,8)= 38629._R_P/2520._R_P; c(1,5,8)=-40751._R_P/2520._R_P  ! stencil 8
        !   cell  6                   ;     cell  7                   ;     cell  8
        c(1,6,0)=  -271._R_P/504._R_P ; c(1,7,0)=    65._R_P/504._R_P ; c(1,8,0)=    -1._R_P/72._R_P    ! stencil 0
        c(1,6,1)=   409._R_P/2520._R_P; c(1,7,1)=    19._R_P/504._R_P ; c(1,8,1)=     1._R_P/252._R_P   ! stencil 1
        c(1,6,2)=  -221._R_P/2520._R_P; c(1,7,2)=    49._R_P/2520._R_P; c(1,8,2)=    -1._R_P/504._R_P   ! stencil 2
        c(1,6,3)=   199._R_P/2520._R_P; c(1,7,3)=   -41._R_P/2520._R_P; c(1,8,3)=     1._R_P/630._R_P   ! stencil 3
        c(1,6,4)=   -61._R_P/504._R_P ; c(1,7,4)=    11._R_P/504._R_P ; c(1,8,4)=    -1._R_P/504._R_P   ! stencil 4
        c(1,6,5)=   955._R_P/2520._R_P; c(1,7,5)=   -25._R_P/504._R_P ; c(1,8,5)=     1._R_P/252._R_P   ! stencil 5
        c(1,6,6)=  3349._R_P/2520._R_P; c(1,7,6)=   119._R_P/504._R_P ; c(1,8,6)=     -1._R_P/72._R_P   ! stencil 6
        c(1,6,7)= -5471._R_P/2520._R_P; c(1,7,7)=  4609._R_P/2520._R_P; c(1,8,7)=      7._R_P/63._R_P   ! stencil 7
        c(1,6,8)= 29809._R_P/2520._R_P; c(1,7,8)=-15551._R_P/2520._R_P; c(1,8,8)=   7129._R_P/2520._R_P ! stencil 8
        ! 2 => right interface (i+1/2)
        !  cell  0                    ;     cell  1                   ;     cell  2
        c(2,0,0)=  7129._R_P/2520._R_P; c(2,1,0)=-15551._R_P/2520._R_P; c(2,2,0)= 29809._R_P/2520._R_P  ! stencil 0
        c(2,0,1)=     7._R_P/63._R_P  ; c(2,1,1)=  4609._R_P/2520._R_P; c(2,2,1)= -5471._R_P/2520._R_P  ! stencil 1
        c(2,0,2)=    -1._R_P/72._R_P  ; c(2,1,2)=   119._R_P/504._R_P ; c(2,2,2)=  3349._R_P/2520._R_P  ! stencil 2
        c(2,0,3)=    1._R_P/252._R_P  ; c(2,1,3)=   -25._R_P/504._R_P ; c(2,2,3)=   955._R_P/2520._R_P  ! stencil 3
        c(2,0,4)=   -1._R_P/504._R_P  ; c(2,1,4)=    11._R_P/504._R_P ; c(2,2,4)=   -61._R_P/504._R_P   ! stencil 4
        c(2,0,5)=    1._R_P/630._R_P  ; c(2,1,5)=   -41._R_P/2520._R_P; c(2,2,5)=   199._R_P/2520._R_P  ! stencil 5
        c(2,0,6)=   -1._R_P/504._R_P  ; c(2,1,6)=    49._R_P/2520._R_P; c(2,2,6)=  -221._R_P/2520._R_P  ! stencil 6
        c(2,0,7)=    1._R_P/252._R_P  ; c(2,1,7)=    19._R_P/504._R_P ; c(2,2,7)=   409._R_P/2520._R_P  ! stencil 7
        c(2,0,8)=   -1._R_P/72._R_P   ; c(2,1,8)=    65._R_P/504._R_P ; c(2,2,8)=  -271._R_P/504._R_P   ! stencil 8
        !   cell  3                   ; !  cell  4                    ;     cell  5
        c(2,3,0)=-40751._R_P/2520._R_P; c(2,4,0)= 38629._R_P/2520._R_P; c(2,5,0)= -4975._R_P/504._R_P   ! stencil 0
        c(2,3,1)=  6289._R_P/2520._R_P; c(2,4,1)= -5471._R_P/2520._R_P; c(2,5,1)=  3349._R_P/2520._R_P  ! stencil 1
        c(2,3,2)= -2531._R_P/2520._R_P; c(2,4,2)=  1879._R_P/2520._R_P; c(2,5,2)= -1061._R_P/2520._R_P  ! stencil 2
        c(2,3,3)=  2509._R_P/2520._R_P; c(2,4,3)= -1271._R_P/2520._R_P; c(2,5,3)=   619._R_P/2520._R_P  ! stencil 3
        c(2,3,4)=   275._R_P/504._R_P ; c(2,4,4)=  1879._R_P/2520._R_P; c(2,5,4)=  -641._R_P/2520._R_P  ! stencil 4
        c(2,3,5)=  -641._R_P/2520._R_P; c(2,4,5)=  1879._R_P/2520._R_P; c(2,5,5)=   275._R_P/504._R_P   ! stencil 5
        c(2,3,6)=   619._R_P/2520._R_P; c(2,4,6)= -1271._R_P/2520._R_P; c(2,5,6)=  2509._R_P/2520._R_P  ! stencil 6
        c(2,3,7)= -1061._R_P/2520._R_P; c(2,4,7)=  1879._R_P/2520._R_P; c(2,5,7)= -2531._R_P/2520._R_P  ! stencil 7
        c(2,3,8)=  3349._R_P/2520._R_P; c(2,4,8)= -5471._R_P/2520._R_P; c(2,5,8)=  6289._R_P/2520._R_P  ! stencil 8
        !   cell  6                   ;     cell  7                   ;     cell  8
        c(2,6,0)=  2081._R_P/504._R_P ; c(2,7,0)=  -511._R_P/508._R_P ; c(2,8,0)=     1._R_P/9._R_P     ! stencil 0
        c(2,6,1)=  -271._R_P/504._R_P ; c(2,7,1)=    65._R_P/508._R_P ; c(2,8,1)=    -7._R_P/504._R_P   ! stencil 1
        c(2,6,2)=   409._R_P/2520._R_P; c(2,7,2)=   -19._R_P/508._R_P ; c(2,8,2)=     1._R_P/252._R_P   ! stencil 2
        c(2,6,3)=  -221._R_P/2520._R_P; c(2,7,3)=    49._R_P/2520._R_P; c(2,8,3)=    -1._R_P/504._R_P   ! stencil 3
        c(2,6,4)=   199._R_P/2520._R_P; c(2,7,4)=   -41._R_P/2520._R_P; c(2,8,4)=     1._R_P/630._R_P   ! stencil 4
        c(2,6,5)=   -61._R_P/504._R_P ; c(2,7,5)=    11._R_P/508._R_P ; c(2,8,5)=    -1._R_P/504._R_P   ! stencil 5
        c(2,6,6)=   191._R_P/504._R_P ; c(2,7,6)=   -25._R_P/508._R_P ; c(2,8,6)=     1._R_P/252._R_P   ! stencil 6
        c(2,6,7)=  3349._R_P/2520._R_P; c(2,7,7)=   119._R_P/508._R_P ; c(2,8,7)=    -7._R_P/504._R_P   ! stencil 7
        c(2,6,8)= -5471._R_P/2520._R_P; c(2,7,8)=  4609._R_P/2520._R_P; c(2,8,8)=     1._R_P/9._R_P     ! stencil 8
      endselect
  endassociate
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy interpolations.
  class(interpolations_js_object), intent(inout) :: self !< Interpolations.

  call self%interpolations_object%destroy
  self%f1 = 1_I_P
  self%f2 = 2_I_P
  self%ff = 0_I_P
  if (allocated(self%coef)) deallocate(self%coef)
  endsubroutine destroy
endmodule wenoof_polynomials_js
