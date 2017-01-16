!< Jiang-Shu and Gerolymos-Senechal-Vallet smoothness indicators object.
module wenoof_smoothness_indicators_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet smoothness indicators object.
!<
!< @note The provided WENO optimal weights implements the smoothness indicators defined in *Efficient Implementation
!< of Weighted ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Senechal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039

use penf, only : I_P, R_P
use wenoof_base_object
use wenoof_smoothness_indicators

implicit none
private
public :: smoothness_indicators_js
public :: smoothness_indicators_js_constructor
public :: create_smoothness_indicators_js_constructor

type, extends(smoothness_indicators_constructor) :: smoothness_indicators_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet smoothness indicators object constructor.
endtype smoothness_indicators_js_constructor

type, extends(smoothness_indicators) :: smoothness_indicators_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet smoothness indicators object.
  !<
  !< @note The provided WENO optimal weights implements the optimal weights defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Sénéchal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  private
  real(R_P), allocatable :: coef(:,:,:) !< Smoothness indicators coefficients [1:2,0:S-1,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute smoothness indicators.
    procedure, nopass     :: description !< Return smoothness indicators string-description.
    ! overridden public methods
    procedure, pass(self) :: create  !< Create smoothness indicators.
    procedure, pass(self) :: destroy !< Destroy smoothness indicators.
endtype smoothness_indicators_js

contains
  ! public non TBP
  subroutine create_smoothness_indicators_js_constructor(S, constructor)
  !< Create smoothness indicators constructor.
  integer(I_P),                                          intent(in)  :: S           !< Stencils dimension.
  class(smoothness_indicators_constructor), allocatable, intent(out) :: constructor !< Smoothness indicators constructor.

  allocate(smoothness_indicators_js_constructor :: constructor)
  constructor%S = S
  endsubroutine create_smoothness_indicators_js_constructor

  ! deferred public methods
  pure subroutine compute(self, S, stencil, f1, f2, ff)
  !< Compute smoothness indicators.
  class(smoothness_indicators_js), intent(inout) :: self                !< Smoothness indicator.
  integer(I_P),                    intent(in)    :: S                   !< Number of stencils actually used.
  real(R_P),                       intent(in)    :: stencil(1:, 1 - S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P),                    intent(in)    :: f1, f2, ff          !< Faces to be computed.
  integer(I_P)                                   :: s1, s2, s3, f       !< Counters

  do s1=0, S - 1 ! stencils loop
    do f=f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%si(f, s1) = 0._R_P
      do s2=0, S - 1
        do s3=0, S - 1
          self%si(f, s1) = self%si(f, s1) + self%coef(s3, s2, s1) * stencil(f + ff, s1 - s3) * stencil(f + ff, s1 - s2)
        enddo
      enddo
    enddo
  enddo
  endsubroutine compute

  pure function description() result(string)
  !< Return smoothness indicators string-description.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line character.

  string = 'WENO smoothness indicators'//nl
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
  !< Create smoothness indicators.
  class(smoothness_indicators_js), intent(inout) :: self        !< Smoothness indicators.
  class(base_object_constructor),  intent(in)    :: constructor !< Smoothness indicators constructor.
  integer(I_P)                                   :: S           !< Stencils dimension.

  call self%destroy
  call self%smoothness_indicators%create(constructor=constructor)
  select type(constructor)
  class is(smoothness_indicators_js_constructor)
    S = constructor%S
    allocate(self%coef(0:S - 1, 0:S - 1, 0:S - 1))
  class default
    ! @TODO add error handling
  endselect
  associate(c => self%coef)
    select case(S)
    case(2) ! 3rd order
      ! stencil 0
      !       i*i      ;     (i-1)*i
      c(0,0,0) = 1._R_P; c(1,0,0) = -2._R_P
      !      /         ;     (i-1)*(i-1)
      c(0,1,0) = 0._R_P; c(1,1,0) =  1._R_P
      ! stencil 1
      !  (i+1)*(i+1)   ;     (i+1)*i
      c(0,0,1) = 1._R_P; c(1,0,1) = -2._R_P
      !      /         ;      i*i
      c(0,1,1) = 0._R_P; c(1,1,1) =  1._R_P
    case(3) ! 5th order
      ! stencil 0
      !      i*i                ;       (i-1)*i             ;       (i-2)*i
      c(0,0,0) =  10._R_P/3._R_P; c(1,0,0) = -31._R_P/3._R_P; c(2,0,0) =  11._R_P/3._R_P
      !      /                  ;       (i-1)*(i-1)         ;       (i-2)*(i-1)
      c(0,1,0) =   0._R_P       ; c(1,1,0) =  25._R_P/3._R_P; c(2,1,0) = -19._R_P/3._R_P
      !      /                  ;        /                  ;       (i-2)*(i-2)
      c(0,2,0) =   0._R_P       ; c(1,2,0) =   0._R_P       ; c(2,2,0) =   4._R_P/3._R_P
      ! stencil 1
      !     (i+1)*(i+1)         ;        i*(i+1)            ;       (i-1)*(i+1)
      c(0,0,1) =   4._R_P/3._R_P; c(1,0,1) = -13._R_P/3._R_P; c(2,0,1) =   5._R_P/3._R_P
      !      /                  ;        i*i                ;       (i-1)*i
      c(0,1,1) =   0._R_P       ; c(1,1,1) =  13._R_P/3._R_P; c(2,1,1) = -13._R_P/3._R_P
      !      /                  ;        /                  ;       (i-1)*(i-1)
      c(0,2,1) =   0._R_P       ; c(1,2,1) =   0._R_P       ; c(2,2,1) =   4._R_P/3._R_P
      ! stencil 2
      !     (i+2)*(i+2)         ;       (i+1)*(i+2)         ;        i*(i+2)
      c(0,0,2) =   4._R_P/3._R_P; c(1,0,2) = -19._R_P/3._R_P; c(2,0,2) =  11._R_P/3._R_P
      !      /                  ;       (i+1)*(i+1)         ;        i*(i+1)
      c(0,1,2) =   0._R_P       ; c(1,1,2) =  25._R_P/3._R_P; c(2,1,2) = -31._R_P/3._R_P
      !      /                  ;        /                  ;        i*i
      c(0,2,2) =   0._R_P       ; c(1,2,2) =   0._R_P       ; c(2,2,2) =  10._R_P/3._R_P
    case(4) ! 7th order
      ! stencil 0
      !      i*i           ;       (i-1)*i        ;       (i-2)*i         ;       (i-3)*i
      c(0,0,0) = 2107._R_P ; c(1,0,0) =-9402._R_P ; c(2,0,0) =  7042._R_P ; c(3,0,0) = -1854._R_P
      !      /             ;       (i-1)*(i-1)    ;       (i-2)*(i-1)     ;       (i-3)*(i-1)
      c(0,1,0) =    0._R_P ; c(1,1,0) =11003._R_P ; c(2,1,0) =-17246._R_P ; c(3,1,0) =  4642._R_P
      !      /             ;        /             ;       (i-2)*(i-2)     ;       (i-3)*(i-2)
      c(0,2,0) =    0._R_P ; c(1,2,0) =    0._R_P ; c(2,2,0) =  7043._R_P ; c(3,2,0) = -3882._R_P
      !      /             ;        /             ;        /              ;       (i-3)*(i-3)
      c(0,3,0) =    0._R_P ; c(1,3,0) =    0._R_P ; c(2,3,0) =     0._R_P ; c(3,3,0) =   547._R_P
      ! stencil 1
      !    (i+1)*(i+1)    ;        i*(i+1)      ;       (i-1)*(i+1)    ;       (i-2)*(i+1)
      c(0,0,1) =  547._R_P; c(1,0,1) =-2522._R_P; c(2,0,1) =  1922._R_P; c(3,0,1) =  -494._R_P
      !     /             ;        i*i          ;      (i-1)*i         ;       (i-2)*i
      c(0,1,1) =    0._R_P; c(1,1,1) = 3443._R_P; c(2,1,1) = -5966._R_P; c(3,1,1) =  1602._R_P
      !     /             ;        /            ;      (i-1)*(i-1)     ;       (i-2)*(i-1)
      c(0,2,1) =    0._R_P; c(1,2,1) =    0._R_P; c(2,2,1) =  2843._R_P; c(3,2,1) = -1642._R_P
      !     /             ;        /            ;       /              ;       (i-2)*(i-2)
      c(0,3,1) =    0._R_P; c(1,3,1) =    0._R_P; c(2,3,1) =     0._R_P; c(3,3,1) =   267._R_P
      ! stencil 2
      !     (i+2)*(i+2)   ;     (i+1)*(i+2)     ;       i*(i+2)        ;      (i-1)*(i+2)
      c(0,0,2) =  267._R_P; c(1,0,2) =-1642._R_P; c(2,0,2) =  1602._R_P; c(3,0,2) =  -494._R_P
      !     /             ;     (i+1)*(i+1)     ;       i*(i+1)        ;      (i-1)*(i+1)
      c(0,1,2) =    0._R_P; c(1,1,2) = 2843._R_P; c(2,1,2) = -5966._R_P; c(3,1,2) =  1922._R_P
      !     /             ;       /             ;       i*i            ;      (i-1)*i
      c(0,2,2) =    0._R_P; c(1,2,2) =    0._R_P; c(2,2,2) =  3443._R_P; c(3,2,2) = -2522._R_P
      !     /             ;       /             ;       /              ;      (i-1)*(i-1)
      c(0,3,2) =    0._R_P; c(1,3,2) =    0._R_P; c(2,3,2) =     0._R_P; c(3,3,2) =   547._R_P
      ! stencil 3
      !    (i+3)*(i+3)    ;      (i+2)*(i+3)    ;      (i+1)*(i+3)     ;       i*(i+3)
      c(0,0,3) =  547._R_P; c(1,0,3) =-3882._R_P; c(2,0,3) =  4642._R_P; c(3,0,3) = -1854._R_P
      !     /             ;      (i+2)*(i+2)    ;      (i+1)*(i+2)     ;       i*(i+2)
      c(0,1,3) =    0._R_P; c(1,1,3) = 7043._R_P; c(2,1,3) =-17246._R_P; c(3,1,3) =  7042._R_P
      !     /             ;       /             ;      (i+1)*(i+1)     ;       i*(i+1)
      c(0,2,3) =    0._R_P; c(1,2,3) =    0._R_P; c(2,2,3) = 11003._R_P; c(3,2,3) = -9402._R_P
      !     /             ;       /             ;       /              ;       i*i
      c(0,3,3) =    0._R_P; c(1,3,3) =    0._R_P; c(2,3,3) =     0._R_P; c(3,3,3) =  2107._R_P
    case(5) ! 9th order
      ! stencil 0
      !              i*i                 ;             (i-1)*i                ;             (i-2)*i
      c(0,0,0) =   53959._R_P / 2520._R_P; c(1,0,0) = -649501._R_P / 5040._R_P; c(2,0,0) =  252941._R_P / 1680._R_P
      !          (i-3)*i                 ;             (i-4)*i
      c(3,0,0) = -411487._R_P / 5040._R_P; c(4,0,0) =   86329._R_P / 5040._R_P
      !               /                  ;             (i-1)*(i-1)            ;             (i-2)*(i-1)
      c(0,1,0) =       0._R_P            ; c(1,1,0) = 1020563._R_P / 5040._R_P; c(2,1,0) =  -68391._R_P /  140._R_P
      !          (i-3)*(i-1)             ;             (i-4)*(i-1)
      c(3,1,0) =  679229._R_P / 2520._R_P; c(4,1,0) = -288007._R_P / 5040._R_P
      !               /                  ;                  /                 ;             (i-2)*(i-2)
      c(0,2,0) =       0._R_P            ; c(1,2,0) =       0._R_P            ; c(2,2,0) =  507131._R_P / 1680._R_P
      !          (i-3)*(i-2)             ;             (i-4)*(i-2)
      c(3,2,0) = -142033._R_P /  420._R_P; c(4,2,0) =  121621._R_P / 1680._R_P
      !               /                  ;                  /                 ;                  /
      c(0,3,0) =       0._R_P            ; c(1,3,0) =       0._R_P            ; c(2,3,0) =       0._R_P
      !          (i-3)*(i-3)             ;             (i-4)*(i-3)
      c(3,3,0) =  482963._R_P / 5040._R_P; c(4,3,0) = -208501._R_P / 5040._R_P
      !               /                  ;                  /                 ;                  /
      c(0,4,0) =       0._R_P            ; c(1,4,0) =       0._R_P            ; c(2,4,0) =       0._R_P
      !               /                  ;             (i-4)*(i-4)
      c(3,4,0) =       0._R_P            ; c(4,4,0) =   11329._R_P / 2520._R_P
      ! stencil 1
      !          (i+1)*(i+1)             ;                 i*(i+1)            ;             (i-1)*(i+1)
      c(0,0,1) =   11329._R_P / 2520._R_P; c(1,0,1) = -140251._R_P / 5040._R_P; c(2,0,1) =   55051._R_P / 1680._R_P
      !          (i-2)*(i+1)             ;             (i-3)*(i+1)
      c(3,0,1) =  -88297._R_P / 5040._R_P; c(4,0,1) =   18079._R_P / 5040._R_P
      !               /                  ;                 i*i                ;             (i-1)*i
      c(0,1,1) =       0._R_P            ; c(1,1,1) =  242723._R_P / 5040._R_P; c(2,1,1) =  -25499._R_P /  210._R_P
      !          (i-2)*i                 ;             (i-3)*i
      c(3,1,1) =  168509._R_P / 2520._R_P; c(4,1,1) =  -70237._R_P / 5040._R_P
      !               /                  ;                  /                 ;             (i-1)*(i-1)
      c(0,2,1) =       0._R_P            ; c(1,2,1) =       0._R_P            ; c(2,2,1) =  135431._R_P / 1680._R_P
      !          (i-2)*(i-1)             ;             (i-3)*(i-1)
      c(3,2,1) =   -3229._R_P /   35._R_P; c(4,2,1) =   33071._R_P / 1680._R_P
      !               /                  ;                  /                 ;                  /
      c(0,3,1) =       0._R_P            ; c(1,3,1) =       0._R_P            ; c(2,3,1) =       0._R_P
      !          (i-2)*(i-2)             ;             (i-3)*(i-2)
      c(3,3,1) =  138563._R_P / 5040._R_P; c(4,3,1) =  -60871._R_P / 5040._R_P
      !               /                  ;                  /                 ;                  /
      c(0,4,1) =       0._R_P            ; c(1,4,1) =       0._R_P            ; c(2,4,1) =       0._R_P
      !               /                  ;             (i-3)*(i-3)
      c(3,4,1) =       0._R_P            ; c(4,4,1) =    1727._R_P / 1260._R_P
      ! stencil 2
      !          (i+2)*(i+2)             ;             (i+1)*(i+2)            ;                 i*(i+2)
      c(0,0,2) =    1727._R_P / 1260._R_P; c(1,0,2) =  -51001._R_P / 5040._R_P; c(2,0,2) =    7547._R_P /  560._R_P
      !          (i-1)*(i+2)             ;             (i-2)*(i+2)
      c(3,0,2) =  -38947._R_P / 5040._R_P; c(4,0,2) =    8209._R_P / 5040._R_P
      !               /                  ;             (i+1)*(i+1)            ;                 i*(i+1)
      c(0,1,2) =       0._R_P            ; c(1,1,2) =  104963._R_P / 5040._R_P; c(2,1,2) =  -24923._R_P /  420._R_P
      !          (i-1)*(i+1)             ;             (i-2)*(i+1)
      c(3,1,2) =   89549._R_P / 2520._R_P; c(4,1,2) =  -38947._R_P / 5040._R_P
      !               /                  ;                  /                 ;                 i*i
      c(0,2,2) =       0._R_P            ; c(1,2,2) =       0._R_P            ; c(2,2,2) =   77051._R_P / 1680._R_P
      !          (i-1)*i                 ;             (i-2)*i
      c(3,2,2) =  -24923._R_P /  420._R_P; c(4,2,2) =    7547._R_P /  560._R_P
      !               /                  ;                  /                 ;                  /
      c(0,3,2) =       0._R_P            ; c(1,3,2) =       0._R_P            ; c(2,3,2) =       0._R_P
      !          (i-1)*(i-1)             ;             (i-2)*(i-1)
      c(3,3,2) =  104963._R_P / 5040._R_P; c(4,3,2) =  -51001._R_P / 5040._R_P
      !               /                  ;                  /                 ;                  /
      c(0,4,2) =       0._R_P            ; c(1,4,2) =       0._R_P            ; c(2,4,2) =       0._R_P
      !               /                  ;             (i-2)*(i-2)
      c(3,4,2) =       0._R_P            ; c(4,4,2) =    1727._R_P / 1260._R_P
      ! stencil 3
      !          (i+3)*(i+3)             ;             (i+2)*(i+3)            ;             (i+1)*(i+3)
      c(0,0,3) =    1727._R_P / 1260._R_P; c(1,0,3) =  -60871._R_P / 5040._R_P; c(2,0,3) =   33071._R_P / 1680._R_P
      !              i*(i+3)             ;             (i-1)*(i+3)
      c(3,0,3) =  -70237._R_P / 5040._R_P; c(4,0,3) =   18079._R_P / 5040._R_P
      !               /                  ;             (i+2)*(i+2)            ;             (i+1)*(i+2)
      c(0,1,3) =       0._R_P            ; c(1,1,3) =  138563._R_P / 5040._R_P; c(2,1,3) =   -3229._R_P /   35._R_P
      !              i*(i+2)             ;             (i-1)*(i+2)
      c(3,1,3) =  168509._R_P / 2520._R_P; c(4,1,3) =  -88297._R_P / 5040._R_P
      !               /                  ;                  /                 ;             (i+1)*(i+1)
      c(0,2,3) =       0._R_P            ; c(1,2,3) =       0._R_P            ; c(2,2,3) =  135431._R_P / 1680._R_P
      !              i*(i+1)             ;             (i-1)*(i+1)
      c(3,2,3) =  -25499._R_P /  210._R_P; c(4,2,3) =   55051._R_P / 1680._R_P
      !               /                  ;                  /                 ;                  /
      c(0,3,3) =       0._R_P            ; c(1,3,3) =       0._R_P            ; c(2,3,3) =       0._R_P
      !              i*i                 ;             (i-1)*i
      c(3,3,3) =  242723._R_P / 5040._R_P; c(4,3,3) = -140251._R_P / 5040._R_P
      !               /                  ;                  /                 ;                  /
      c(0,4,3) =       0._R_P            ; c(1,4,3) =       0._R_P            ; c(2,4,3) =       0._R_P
      !               /                  ;             (i-1)*(i-1)
      c(3,4,3) =       0._R_P            ; c(4,4,3) =   11329._R_P / 2520._R_P
      ! stencil 4
      !          (i+4)*(i+4)             ;             (i+3)*(i+4)            ;             (i+2)*(i+4)
      c(0,0,4) =   11329._R_P / 2520._R_P; c(1,0,4) = -208501._R_P / 5040._R_P; c(2,0,4) =  121621._R_P / 1680._R_P
      !          (i+1)*(i+4)             ;                 i*(i+4)
      c(3,0,4) = -288007._R_P / 5040._R_P; c(4,0,4) =   86329._R_P / 5040._R_P
      !               /                  ;             (i+3)*(i+3)            ;             (i+2)*(i+3)
      c(0,1,4) =       0._R_P            ; c(1,1,4) =  482963._R_P / 5040._R_P; c(2,1,4) = -142033._R_P /  420._R_P
      !          (i+1)*(i+3)             ;                 i*(i+3)
      c(3,1,4) =  679229._R_P / 2520._R_P; c(4,1,4) = -411487._R_P / 5040._R_P
      !               /                  ;                  /                 ;             (i+1)*(i+2)
      c(0,2,4) =       0._R_P            ; c(1,2,4) =       0._R_P            ; c(2,2,4) =  507131._R_P / 1680._R_P
      !          (i+1)*(i+2)             ;                 i*(i+2)
      c(3,2,4) =  -68391._R_P /  140._R_P; c(4,2,4) =  252941._R_P / 1680._R_P
      !               /                  ;                  /                 ;                  /
      c(0,3,4) =       0._R_P            ; c(1,3,4) =       0._R_P            ; c(2,3,4) =       0._R_P
      !          (i+1)*(i+1)             ;                 i*(i+1)
      c(3,3,4) = 1020563._R_P / 5040._R_P; c(4,3,4) = -649501._R_P / 5040._R_P
      !               /                  ;                  /                 ;                  /
      c(0,4,4) =       0._R_P            ; c(1,4,4) =       0._R_P            ; c(2,4,4) =       0._R_P
      !               /                  ;                 i*i
      c(3,4,4) =       0._R_P            ; c(4,4,4) =   53959._R_P / 2520._R_P
    case(6) ! 11th order
      ! stencil 0
      !                 i*i                  ;                (i-1)*i                 ;                 (i-2)*i
      c(0,0,0) =   6150211._R_P / 120960._R_P; c(1,0,0) =  -2966279._R_P /   7560._R_P; c(2,0,0) =   4762921._R_P /   7560._R_P
      !             (i-3)*i                  ;                (i-4)*i                 ;                 (i-5)*i
      c(3,0,0) = -15848531._R_P /  30240._R_P; c(4,0,0) =   2706017._R_P /  12096._R_P; c(5,0,0) =   -235637._R_P /   6048._R_P

      !                  /                   ;                (i-1)*(i-1)             ;                 (i-2)*(i-1)
      c(0,1,0) =         0._R_P              ; c(1,1,0) =  31617079._R_P /  40320._R_P; c(2,1,0) = -25980937._R_P /  10080._R_P
      !             (i-3)*(i-1)              ;                (i-4)*(i-1)             ;                 (i-5)*(i-1)
      c(3,1,0) =  32862709._R_P /  15120._R_P; c(4,1,0) =  -1048211._R_P /   1120._R_P; c(5,1,0) =    661145._R_P /   4032._R_P

      !                  /                   ;                     /                  ;                 (i-2)*(i-2)
      c(0,2,0) =         0._R_P              ; c(1,2,0) =         0._R_P              ; c(2,2,0) =  21703781._R_P /  10080._R_P
      !              (i-3)*(i-2)             ;                (i-4)*(i-2)             ;                 (i-5)*(i-2)
      c(3,2,0) =  -6937561._R_P /   1890._R_P; c(4,2,0) =   2674951._R_P /   1680._R_P; c(5,2,0) =   -314063._R_P /   1120._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,3,0) =         0._R_P              ; c(1,3,0) =         0._R_P              ; c(2,3,0) =         0._R_P
      !             (i-3)*(i-3)              ;                (i-4)*(i-3)             ;                 (i-5)*(i-3)
      c(3,3,0) =  47689393._R_P /  30240._R_P; c(4,3,0) = -41615261._R_P /  30240._R_P; c(5,3,0) =   1840141._R_P /   7560._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,4,0) =         0._R_P              ; c(1,4,0) =         0._R_P              ; c(2,4,0) =         0._R_P
      !                  /                   ;                (i-4)*(i-4)             ;                 (i-5)*(i-4)
      c(3,4,0) =         0._R_P              ; c(4,4,0) =  12160229._R_P /  40320._R_P; c(5,4,0) =   -539591._R_P /   5040._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,5,0) =         0._R_P              ; c(1,5,0) =         0._R_P              ; c(2,5,0) =         0._R_P
      !                  /                   ;                     /                  ;                 (i-5)*(i-5)
      c(3,5,0) =         0._R_P              ; c(4,5,0) =         0._R_P              ; c(5,5,0) =    384187._R_P /  40320._R_P
      ! stencil 1
      !             (i+1)*(i+1)              ;                    i*(i+1)             ;                 (i-1)*(i+1)
      c(0,0,1) =    384187._R_P /  40320._R_P; c(1,0,1) =  -1139749._R_P /  15120._R_P; c(2,0,1) =     61427._R_P /   504._R_P;
      !             (i-2)*(i+1)              ;                (i-3)*(i+1)             ;                 (i-4)*(i+1)
      c(3,0,1) =  -1015303._R_P /  10080._R_P; c(4,0,1) =   2567287._R_P /  60480._R_P; c(5,0,1) =    -73379._R_P / 10080._R_P;

      !                  /                   ;                    i*i                 ;                 (i-1)*i
      c(0,1,1) =         0._R_P              ; c(1,1,1) =  19365967._R_P / 120960._R_P; c(2,1,1) = -16306061._R_P /  30240._R_P
      !             (i-2)*i                  ;                (i-3)*i                 ;                 (i-4)*i
      c(3,1,1) =   6881719._R_P /  15120._R_P; c(4,1,1) =  -5877617._R_P /  30240._R_P; c(5,1,1) =   2033509._R_P /  60480._R_P

      !                  /                   ;                     /                  ;                 (i-1)*(i-1)
      c(0,2,1) =         0._R_P              ; c(1,2,1) =         0._R_P              ; c(2,2,1) =   4721851._R_P /  10080._R_P
      !             (i-2)*(i-1)              ;                (i-3)*(i-1)             ;                 (i-4)*(i-1)
      c(3,2,1) =   -169859._R_P /    210._R_P; c(4,2,1) =   5300629._R_P /  15120._R_P; c(5,2,1) =    -68601._R_P /   1120._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,3,1) =         0._R_P              ; c(1,3,1) =         0._R_P              ; c(2,3,1) =         0._R_P
      !             (i-2)*(i-2)              ;                (i-3)*(i-2)             ;                 (i-4)*(i-2)
      c(3,3,1) =   1197047._R_P /   3360._R_P; c(4,3,1) =  -9478331._R_P /  30240._R_P; c(5,3,1) =    139471._R_P /   2520._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,4,1) =         0._R_P              ; c(1,4,1) =         0._R_P              ; c(2,4,1) =         0._R_P
      !                  /                   ;                (i-3)*(i-3)             ;                 (i-4)*(i-3)
      c(3,4,1) =         0._R_P              ; c(4,4,1) =   8449957._R_P / 120960._R_P; c(5,4,1) =   -188483._R_P /   7560._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,5,1) =         0._R_P              ; c(1,5,1) =         0._R_P              ; c(2,5,1) =         0._R_P
      !                  /                   ;                     /                  ;                 (i-4)*(i-4)
      c(3,5,1) =         0._R_P              ; c(4,5,1) =         0._R_P              ; c(5,5,1) =     90593._R_P /  40320._R_P
      ! stencil 2
      !             (i+2)*(i+2)              ;                (i+1)*(i+2)             ;                     i*(i+2)
      c(0,0,2) =     90593._R_P /  40320._R_P; c(1,0,2) =     -1240._R_P /     63._R_P; c(2,0,2) =    255397._R_P /   7560._R_P
      !             (i-1)*(i+2)              ;                (i-2)*(i+2)             ;                 (i-3)*(i+2)
      c(3,0,2) =   -288521._R_P /  10080._R_P; c(4,0,2) =    243127._R_P /  20160._R_P; c(5,0,2) =    -12281._R_P /   6048._R_P

      !                  /                   ;                (i+1)*(i+1)             ;                     i*(i+1)
      c(0,1,2) =         0._R_P              ; c(1,1,2) =   1884439._R_P /  40320._R_P; c(2,1,2) =  -5106971._R_P /  30240._R_P
      !             (i-1)*(i+1)              ;                (i-2)*(i+1)             ;                (i-3)*(i+1)
      c(3,1,2) =    248681._R_P /   1680._R_P; c(4,1,2) =   -643999._R_P /  10080._R_P; c(5,1,2) =    662503._R_P /  60480._R_P

      !                  /                   ;                     /                  ;                    i*i
      c(0,2,2) =         0._R_P              ; c(1,2,2) =         0._R_P              ; c(2,2,2) =   4877743._R_P /  30240._R_P
      !             (i-1)*i                  ;                (i-2)*i                 ;                (i-3)*i
      c(3,2,2) =   -559651._R_P /   1890._R_P; c(4,2,2) =   1991239._R_P /  15120._R_P; c(5,2,2) =   -139633._R_P /   6048._R_P

      !                  /                   ;                     /                  ;                     /
      c(0,3,2) =         0._R_P              ; c(1,3,2) =         0._R_P              ; c(2,3,2) =         0._R_P
      !             (i-1)*(i-1)              ;                (i-2)*(i-1)             ;                (i-3)*(i-1)
      c(3,3,2) =    159219._R_P /   1120._R_P; c(4,3,2) =  -1323367._R_P /  10080._R_P; c(5,3,2) =    178999._R_P /   7560._R_P

      !                  /                   ;                     /                  ;                     /
      c(0,4,2) =         0._R_P              ; c(1,4,2) =         0._R_P              ; c(2,4,2) =         0._R_P
      !                  /                   ;                (i-2)*(i-2)             ;                (i-3)*(i-2)
      c(3,4,2) =         0._R_P              ; c(4,4,2) =    141661._R_P /   4480._R_P; c(5,4,2) =   -178747._R_P /  15120._R_P

      !                  /                   ;                     /                  ;                     /
      c(0,5,2) =         0._R_P              ; c(1,5,2) =         0._R_P              ; c(2,5,2) =         0._R_P
      !                  /                   ;                     /                  ;                (i-3)*(i-3)
      c(3,5,2) =         0._R_P              ; c(4,5,2) =         0._R_P              ; c(5,5,2) =    139633._R_P / 120960._R_P
      ! stencil 3
      !             (i+3)*(i+3)              ;                (i+2)*(i+3)             ;                 (i+1)*(i+3)
      c(0,0,3) =    139633._R_P / 120960._R_P; c(1,0,3) =   -178747._R_P /  15120._R_P; c(2,0,3) =    178999._R_P /   7560._R_P
      !                 i*(i+3)              ;                (i-1)*(i+3)             ;                 (i-2)*(i+3)
      c(3,0,3) =   -139633._R_P /   6048._R_P; c(4,0,3) =    662503._R_P /  60480._R_P; c(5,0,3) =    -12281._R_P /   6048._R_P

      !                  /                   ;                (i+2)*(i+2)             ;                 (i+1)*(i+2)
      c(0,1,3) =         0._R_P              ; c(1,1,3) =    141661._R_P /   4480._R_P; c(2,1,3) =  -1323367._R_P /  10080._R_P
      !                 i*(i+2)              ;                (i-1)*(i+2)             ;                 (i-2)*(i+2)
      c(3,1,3) =   1991239._R_P /  15120._R_P; c(4,1,3) =   -643999._R_P /  10080._R_P; c(5,1,3) =    243127._R_P /  20160._R_P

      !                  /                   ;                     /                  ;                    i*(i+1)
      c(0,2,3) =         0._R_P              ; c(1,2,3) =         0._R_P              ; c(2,2,3) =    159219._R_P /   1120._R_P
      !                 i*(i+1)              ;                (i-1)*(i+1)             ;                (i-2)*(i+1)
      c(3,2,3) =   -559651._R_P /   1890._R_P; c(4,2,3) =    248681._R_P /   1680._R_P; c(5,2,3) =   -288521._R_P /  10080._R_P

      !                  /                   ;                     /                  ;                     /
      c(0,3,3) =         0._R_P              ; c(1,3,3) =         0._R_P              ; c(2,3,3) =         0._R_P
      !                 i*i                  ;                (i-1)*i                 ;                (i-2)*i
      c(3,3,3) =   4877743._R_P /  30240._R_P; c(4,3,3) =  -5106971._R_P /  30240._R_P; c(5,3,3) =    255397._R_P /   7560._R_P

      !                  /                   ;                     /                  ;                     /
      c(0,4,3) =         0._R_P              ; c(1,4,3) =         0._R_P              ; c(2,4,3) =         0._R_P
      !                  /                   ;                (i-1)*(i-1)             ;                (i-2)*(i-1)
      c(3,4,3) =         0._R_P              ; c(4,4,3) =   1884439._R_P /  40320._R_P; c(5,4,3) =     -1240._R_P /     63._R_P

      !                  /                   ;                     /                  ;                     /
      c(0,5,3) =         0._R_P              ; c(1,5,3) =         0._R_P              ; c(2,5,3) =         0._R_P
      !                  /                   ;                     /                  ;                (i-2)*(i-2)
      c(3,5,3) =         0._R_P              ; c(4,5,3) =         0._R_P              ; c(5,5,3) =     90593._R_P /  40320._R_P
      ! stencil 4
      !             (i+4)*(i+4)              ;                (i+3)*(i+4)             ;                 (i+2)*(i+4)
      c(0,0,4) =     90593._R_P /  40320._R_P; c(1,0,4) =   -188483._R_P /   7560._R_P; c(2,0,4) =    139471._R_P /   2520._R_P
      !             (i+1)*(i+4)              ;                    i*(i+4)             ;                 (i-1)*(i+4)
      c(3,0,4) =    -68601._R_P /   1120._R_P; c(4,0,4) =   2033509._R_P /  60480._R_P; c(5,0,4) =    -73379._R_P /  10080._R_P

      !                  /                   ;                (i+3)*(i+3)             ;                 (i+2)*(i+3)
      c(0,1,4) =         0._R_P              ; c(1,1,4) =   8449957._R_P / 120960._R_P; c(2,1,4) =  -9478331._R_P /  30240._R_P
      !             (i+1)*(i+3)              ;                    i*(i+3)             ;                 (i-1)*(i+3)
      c(3,1,4) =   5300629._R_P /  15120._R_P; c(4,1,4) =  -5877617._R_P /  30240._R_P; c(5,1,4) =   2567287._R_P /  60480._R_P

      !                  /                   ;                     /                  ;                 (i+2)*(i+2)
      c(0,2,4) =         0._R_P              ; c(1,2,4) =         0._R_P              ; c(2,2,4) =   1197047._R_P /   3360._R_P
      !             (i+1)*(i+2)              ;                    i*(i+2)             ;                 (i-1)*(i+2)
      c(3,2,4) =   -169859._R_P /    210._R_P; c(4,2,4) =   6881719._R_P /  15120._R_P; c(5,2,4) =  -1015303._R_P /  10080._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,3,4) =         0._R_P              ; c(1,3,4) =         0._R_P              ; c(2,3,4) =         0._R_P
      !             (i+1)*(i+1)              ;                    i*(i+1)             ;                 (i-1)*(i+1)
      c(3,3,4) =   4721851._R_P /  10080._R_P; c(4,3,4) = -16306061._R_P /  30240._R_P; c(5,3,4) =     61427._R_P /    504._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,4,4) =         0._R_P              ; c(1,4,4) =         0._R_P              ; c(2,4,4) =         0._R_P
      !                  /                   ;                    i*i                 ;                 (i-1)*i
      c(3,4,4) =         0._R_P              ; c(4,4,4) =  19365967._R_P / 120960._R_P; c(5,4,4) =  -1139749._R_P /  15120._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,5,4) =         0._R_P              ; c(1,5,4) =         0._R_P              ; c(2,5,4) =         0._R_P
      !                  /                   ;                     /                  ;                 (i-1)*(i-1)
      c(3,5,4) =         0._R_P              ; c(4,5,4) =         0._R_P              ; c(5,5,4) =    384187._R_P /  40320._R_P
      ! stencil 5
      !             (i+5)*(i+5)              ;                (i+4)*(i+5)             ;                 (i+3)*(i+5)
      c(0,0,5) =    384187._R_P /  40320._R_P; c(1,0,5) =   -539591._R_P /   5040._R_P; c(2,0,5) =   1840141._R_P /   7560._R_P
      !             (i+2)*(i+5)              ;                (i+1)*(i+5)             ;                     i*(i+5)
      c(3,0,5) =   -314063._R_P /   1120._R_P; c(4,0,5) =    661145._R_P /   4032._R_P; c(5,0,5) =   -235637._R_P /   6048._R_P

      !                  /                   ;                (i+4)*(i+3)             ;                 (i+3)*(i+3)
      c(0,1,5) =         0._R_P              ; c(1,1,5) =  12160229._R_P /  40320._R_P; c(2,1,5) = -41615261._R_P /  30240._R_P
      !             (i+2)*(i+3)              ;                (i+1)*(i+3)             ;                     i*(i+3)
      c(3,1,5) =   2674951._R_P /   1680._R_P; c(4,1,5) =  -1048211._R_P /   1120._R_P; c(5,1,5) =   2706017._R_P /  12096._R_P

      !                  /                   ;                     /                  ;                 (i+3)*(i+2)
      c(0,2,5) =         0._R_P              ; c(1,2,5) =         0._R_P              ; c(2,2,5) =  47689393._R_P /  30240._R_P
      !             (i+2)*(i+2)              ;                (i+1)*(i+2)             ;                     i*(i+2)
      c(3,2,5) =  -6937561._R_P /   1890._R_P; c(4,2,5) =  32862709._R_P /  15120._R_P; c(5,2,5) = -15848531._R_P /  30240._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,3,5) =         0._R_P              ; c(1,3,5) =         0._R_P              ; c(2,3,5) =         0._R_P
      !             (i+2)*(i+1)              ;                (i+1)*(i+1)             ;                     i*(i+1)
      c(3,3,5) =  21703781._R_P /  10080._R_P; c(4,3,5) = -25980937._R_P /  10080._R_P; c(5,3,5) =   4762921._R_P /   7560._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,4,5) =         0._R_P              ; c(1,4,5) =         0._R_P              ; c(2,4,5) =         0._R_P
      !                  /                   ;                (i+1)*i                 ;                     i*i
      c(3,4,5) =         0._R_P              ; c(4,4,5) =  31617079._R_P /  40320._R_P; c(5,4,5) =  -2966279._R_P /   7560._R_P

      !                  /                   ;                     /                  ;                      /
      c(0,5,5) =         0._R_P              ; c(1,5,5) =         0._R_P              ; c(2,5,5) =         0._R_P
      !                  /                   ;                     /                  ;                     i*(i-1)
      c(3,5,5) =         0._R_P              ; c(4,5,5) =         0._R_P              ; c(5,5,5) =   6150211._R_P / 120960._R_P
    case(7) ! 13th order
      ! stencil 0
      !                   i*i                   ;                   (i-1)*i
      c(0,0,0) =     897207163._R_P/7484400._R_P; c(1,0,0) = -22763092357._R_P/19958400._R_P
      !                  (i-2)*i                ;                   (i-3)*i
      c(2,0,0) =  46808583631._R_P/19958400._R_P; c(3,0,0) = -39645439643._R_P/14968800._R_P
      !                  (i-4)*i                ;                   (i-5)*i
      c(4,0,0) =    8579309749._R_P/4989600._R_P; c(5,0,0) =   -2416885043._R_P/3991680._R_P
      !                  (i-6)*i
      c(6,0,0) =   5391528799._R_P/59875200._R_P

      !                    /                    ;                   (i-1)*(i-1)
      c(0,1,0) =                          0._R_P; c(1,1,0) =    6182612731._R_P/2217600._R_P
      !                  (i-2)*(i-1)            ;                   (i-3)*(i-1)
      c(2,1,0) =    -8623431623._R_P/739200._R_P; c(3,1,0) =   66440049371._R_P/4989600._R_P
      !                  (i-4)*(i-1)            ;                   (i-5)*(i-1)
      c(4,1,0) =  -19308505679._R_P/2217600._R_P; c(5,1,0) =    3417057367._R_P/1108800._R_P
      !                  (i-6)*(i-1)
      c(6,1,0) =  -9181961959._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,2,0) =                          0._R_P; c(1,2,0) =                          0._R_P
      !                  (i-2)*(i-2)            ;                   (i-3)*(i-2)
      c(2,2,0) =     1369404749._R_P/110880._R_P; c(3,2,0) =   -28364892607._R_P/997920._R_P
      !                  (i-4)*(i-2)            ;                   (i-5)*(i-2)
      c(4,2,0) =     8290771913._R_P/443520._R_P; c(5,2,0) =  -14734178999._R_P/2217600._R_P
      !                  (i-6)*(i-2)
      c(6,2,0) =    4964771899._R_P/4989600._R_P

      !                    /                    ;                     /
      c(0,3,0) =                          0._R_P; c(1,3,0) =                          0._R_P
      !                    /                    ;                   (i-3)*(i-3)
      c(2,3,0) =                          0._R_P; c(3,3,0) =   49256859919._R_P/2993760._R_P
      !                  (i-4)*(i-3)            ;                   (i-5)*(i-3)
      c(4,3,0) =   -21693002767._R_P/997920._R_P; c(5,3,0) =   38683385051._R_P/4989600._R_P
      !                  (i-6)*(i-3)
      c(6,3,0) = -17425032203._R_P/14968800._R_P

      !                    /                    ;                     /
      c(0,4,0) =                          0._R_P; c(1,4,0) =                          0._R_P
      !                    /                    ;                     /
      c(2,4,0) =                          0._R_P; c(3,4,0) =                          0._R_P
      !                  (i-4)*(i-4)            ;                   (i-5)*(i-4)
      c(4,4,0) =       199730921._R_P/27720._R_P; c(5,4,0) =    -3809437823._R_P/739200._R_P
      !                  (i-6)*(i-4)
      c(6,4,0) =  15476926351._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,5,0) =                          0._R_P; c(1,5,0) =                          0._R_P
      !                    /                    ;                     /
      c(2,5,0) =                          0._R_P; c(3,5,0) =                          0._R_P
      !                    /                    ;                   (i-5)*(i-5)
      c(4,5,0) =                          0._R_P; c(5,5,0) =    2047941883._R_P/2217600._R_P
      !                  (i-6)*(i-5)
      c(6,5,0) = -5556669277._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,6,0) =                          0._R_P; c(1,6,0) =                          0._R_P
      !                    /                    ;                     /
      c(2,6,0) =                          0._R_P; c(3,6,0) =                          0._R_P
      !                    /                    ;                     /
      c(4,6,0) =                          0._R_P; c(5,6,0) =                          0._R_P
      !                  (i-6)*(i-6)
      c(6,6,0) =      62911297._R_P/2993760._R_P
      ! stencil 1
      !                  (i+1)*(i+1)            ;                    i*(i+1)
      c(0,0,1) =      62911297._R_P/2993760._R_P; c(1,0,1) =  -4074544787._R_P/19958400._R_P
      !                  (i-1)*(i+1)            ;                   (i-2)*(i+1)
      c(2,0,1) =    2811067067._R_P/6652800._R_P; c(3,0,1) =  -7124638253._R_P/14968800._R_P
      !                  (i-3)*(i+1)            ;                   (i-4)*(i+1)
      c(4,0,1) =    1531307249._R_P/4989600._R_P; c(5,0,1) =    -712745603._R_P/6652800._R_P
      !                  (i-5)*(i+1)
      c(6,0,1) =    945166329._R_P/59875200._R_P

      !                    /                    ;                       i*(i-1)
      c(0,1,1) =                          0._R_P; c(1,1,1) =       12742497._R_P/246400._R_P
      !                  (i-1)*(i-1)            ;                   (i-2)*(i-1)
      c(2,1,1) =  -14684933057._R_P/6652800._R_P; c(3,1,1) =   12601009501._R_P/4989600._R_P
      !                  (i-3)*(i-1)            ;                   (i-4)*(i-1)
      c(4,1,1) =     -405382961._R_P/246400._R_P; c(5,1,1) =    1924032511._R_P/3326400._R_P
      !                  (i-5)*(i-1)
      c(6,1,1) =    -341910757._R_P/3991680._R_P

      !                    /                    ;                    /
      c(0,2,1) =                          0._R_P; c(1,2,1) =                          0._R_P
      !                  (i-1)*(i-2)            ;                   (i-3)*(i-2)
      c(2,2,1) =      796358777._R_P/332640._R_P; c(3,2,1) =     -616410313._R_P/110880._R_P
      !                  (i-4)*(i-2)            ;                   (i-5)*(i-2)
      c(4,2,1) =    4868089189._R_P/1330560._R_P; c(5,2,1) =   -8619440987._R_P/6652800._R_P
      !                  (i-6)*(i-2)
      c(6,2,1) =     320782183._R_P/1663200._R_P

      !                    /                    ;                     /
      c(0,3,1) =                          0._R_P; c(1,3,1) =                          0._R_P
      !                    /                    ;                   (i-2)*(i-3)
      c(2,3,1) =                          0._R_P; c(3,3,1) =    9780057169._R_P/2993760._R_P
      !                  (i-3)*(i-3)            ;                   (i-4)*(i-3)
      c(4,3,1) =    -4330640057._R_P/997920._R_P; c(5,3,1) =      857838469._R_P/554400._R_P
      !                  (i-5)*(i-3)
      c(6,3,1) =  -3465607493._R_P/14968800._R_P

      !                    /                    ;                     /
      c(0,4,1) =                          0._R_P; c(1,4,1) =                          0._R_P
      !                    /                    ;                     /
      c(2,4,1) =                          0._R_P; c(3,4,1) =                          0._R_P
      !                  (i-3)*(i-4)            ;                   (i-4)*(i-4)
      c(4,4,1) =        53678683._R_P/36960._R_P; c(5,4,1) =   -6932480657._R_P/6652800._R_P
      !                  (i-5)*(i-4)
      c(6,4,1) =  3126718481._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,5,1) =                          0._R_P; c(1,5,1) =                          0._R_P
      !                    /                    ;                     /
      c(2,5,1) =                          0._R_P; c(3,5,1) =                          0._R_P
      !                    /                    ;                   (i-5)*(i-5)
      c(4,5,1) =                          0._R_P; c(5,5,1) =    1250007643._R_P/6652800._R_P
      !                  (i-6)*(i-5)
      c(6,5,1) =    -377474689._R_P/6652800._R_P

      !                    /                    ;                     /
      c(0,6,1) =                          0._R_P; c(1,6,1) =                          0._R_P
      !                    /                    ;                     /
      c(2,6,1) =                          0._R_P; c(3,6,1) =                          0._R_P
      !                    /                    ;                     /
      c(4,6,1) =                          0._R_P; c(5,6,1) =                          0._R_P
      !                  (i-6)*(i-5)
      c(6,6,1) =     64361771._R_P/14968800._R_P
      ! stencil 2
      !                  (i+2)*i                ;                   (i+1)*i
      c(0,0,2) =     64361771._R_P/14968800._R_P; c(1,0,2) =   -295455983._R_P/6652800._R_P
      !                      i*i                ;                   (i-1)*i
      c(2,0,2) =   1894705391._R_P/19958400._R_P; c(3,0,2) = -1618284323._R_P/14968800._R_P
      !                  (i-2)*i                ;                   (i-3)*i
      c(4,0,2) =     115524053._R_P/1663200._R_P; c(5,0,2) =    -95508139._R_P/3991680._R_P
      !                  (i-4)*i
      c(6,0,2) =       8279479._R_P/2395008._R_P

      !                    /                    ;                   (i+1)*(i-1)
      c(0,1,2) =                          0._R_P; c(1,1,2) =     806338417._R_P/6652800._R_P
      !                      i*(i-1)            ;                   (i-1)*(i-1)
      c(2,1,2) =   -3573798407._R_P/6652800._R_P; c(3,1,2) =    1042531337._R_P/1663200._R_P
      !                  (i-2)*(i-1)            ;                   (i-3)*(i-1)
      c(4,1,2) =   -2725575317._R_P/6652800._R_P; c(5,1,2) =     475321093._R_P/3326400._R_P
      !                  (i-4)*(i-1)
      c(6,1,2) =      -15401629._R_P/739200._R_P

      !                    /                    ;                    /
      c(0,2,2) =                          0._R_P; c(1,2,2) =                          0._R_P
          !                  i*(i-1)            ;                   (i-1)*(i-2)
      c(2,2,2) =        34187317._R_P/55440._R_P; c(3,2,2) =    -1476618887._R_P/997920._R_P
      !                  (i-2)*(i-2)            ;                   (i-3)*(i-2)
      c(4,2,2) =    1312114459._R_P/1330560._R_P; c(5,2,2) =    -773749439._R_P/2217600._R_P
      !                  (i-4)*(i-2)
      c(6,2,2) =     256556849._R_P/4989600._R_P

      !                    /                    ;                     /
      c(0,3,2) =                          0._R_P; c(1,3,2) =                          0._R_P
      !                    /                    ;                   (i-1)*(i-3)
      c(2,3,2) =                          0._R_P; c(3,3,2) =    2726585359._R_P/2993760._R_P
      !                  (i-2)*(i-3)            ;                   (i-3)*(i-3)
      c(4,3,2) =     -412424029._R_P/332640._R_P; c(5,3,2) =    2224538011._R_P/4989600._R_P
      !                  (i-4)*(i-3)
      c(6,3,2) =   -995600723._R_P/14968800._R_P

      !                    /                    ;                     /
      c(0,4,2) =                          0._R_P; c(1,4,2) =                          0._R_P
      !                    /                    ;                     /
      c(2,4,2) =                          0._R_P; c(3,4,2) =                          0._R_P
      !                  (i-2)*(i-4)            ;                   (i-3)*(i-4)
      c(4,4,2) =      143270957._R_P/332640._R_P; c(5,4,2) =   -2096571887._R_P/6652800._R_P
      !                  (i-4)*(i-4)
      c(6,4,2) =     105706999._R_P/2217600._R_P

      !                    /                    ;                     /
      c(0,5,2) =                          0._R_P; c(1,5,2) =                          0._R_P
      !                    /                    ;                     /
      c(2,5,2) =                          0._R_P; c(3,5,2) =                          0._R_P
      !                    /                    ;                   (i-3)*(i-5)
      c(4,5,2) =                          0._R_P; c(5,5,2) =     130013563._R_P/2217600._R_P
      !                  (i-4)*(i-5)
      c(6,5,2) =   -359321429._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,6,2) =                          0._R_P; c(1,6,2) =                          0._R_P
      !                    /                    ;                     /
      c(2,6,2) =                          0._R_P; c(3,6,2) =                          0._R_P
      !                    /                    ;                     /
      c(4,6,2) =                          0._R_P; c(5,6,2) =                          0._R_P
      !                  (i-4)*(i-5)
      c(6,6,2) =       2627203._R_P/1871100._R_P
      ! stencil 3
      !                  (i+3)*i                ;                   (i+2)*i
      c(0,0,3) =       2627203._R_P/1871100._R_P; c(1,0,3) =  -323333323._R_P/19958400._R_P
      !                  (i+1)*i                ;                       i*i
      c(2,0,3) =    761142961._R_P/19958400._R_P; c(3,0,3) =  -701563133._R_P/14968800._R_P
      !                  (i-1)*i                ;                   (i-2)*i
      c(4,0,3) =     158544319._R_P/4989600._R_P; c(5,0,3) =  -225623953._R_P/19958400._R_P
      !                  (i-3)*i
      c(6,0,3) =     99022657._R_P/59875200._R_P

      !                    /                    ;                   (i+1)*(i-1)
      c(0,1,3) =                          0._R_P; c(1,1,3) =     108444169._R_P/2217600._R_P
      !                      i*(i-1)            ;                   (i-1)*(i-1)
      c(2,1,3) =     -176498513._R_P/739200._R_P; c(3,1,3) =    1506944981._R_P/4989600._R_P
      !                  (i-2)*(i-1)            ;                   (i-3)*(i-1)
      c(4,1,3) =    -464678369._R_P/2217600._R_P; c(5,1,3) =      84263749._R_P/1108800._R_P
      !                  (i-4)*(i-1)
      c(6,1,3) =   -225623953._R_P/19958400._R_P

      !                    /                    ;                    /
      c(0,2,3) =                          0._R_P; c(1,2,3) =                          0._R_P
          !                  i*(i-1)            ;                   (i-1)*(i-2)
      c(2,2,3) =        16790707._R_P/55440._R_P; c(3,2,3) =     -790531177._R_P/997920._R_P
      !                  (i-2)*(i-2)            ;                   (i-3)*(i-2)
      c(4,2,3) =      250523543._R_P/443520._R_P; c(5,2,3) =    -464678369._R_P/2217600._R_P
      !                  (i-4)*(i-2)
      c(6,2,3) =     158544319._R_P/4989600._R_P

      !                    /                    ;                     /
      c(0,3,3) =                          0._R_P; c(1,3,3) =                          0._R_P
      !                    /                    ;                   (i-1)*(i-3)
      c(2,3,3) =                          0._R_P; c(3,3,3) =    1607739169._R_P/2993760._R_P
      !                  (i-2)*(i-3)            ;                   (i-3)*(i-3)
      c(4,3,3) =     -790531177._R_P/997920._R_P; c(5,3,3) =    1506944981._R_P/4989600._R_P
      !                  (i-4)*(i-3)
      c(6,3,3) =   -701563133._R_P/14968800._R_P

      !                    /                    ;                     /
      c(0,4,3) =                          0._R_P; c(1,4,3) =                          0._R_P
      !                    /                    ;                     /
      c(2,4,3) =                          0._R_P; c(3,4,3) =                          0._R_P
      !                  (i-2)*(i-4)            ;                   (i-3)*(i-4)
      c(4,4,3) =        16790707._R_P/55440._R_P; c(5,4,3) =     -176498513._R_P/739200._R_P
      !                  (i-4)*(i-4)
      c(6,4,3) =    761142961._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,5,3) =                          0._R_P; c(1,5,3) =                          0._R_P
      !                    /                    ;                     /
      c(2,5,3) =                          0._R_P; c(3,5,3) =                          0._R_P
      !                    /                    ;                   (i-3)*(i-5)
      c(4,5,3) =                          0._R_P; c(5,5,3) =     108444169._R_P/2217600._R_P
      !                  (i-4)*(i-5)
      c(6,5,3) =   -323333323._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,6,3) =                          0._R_P; c(1,6,3) =                          0._R_P
      !                    /                    ;                     /
      c(2,6,3) =                          0._R_P; c(3,6,3) =                          0._R_P
      !                    /                    ;                     /
      c(4,6,3) =                          0._R_P; c(5,6,3) =                          0._R_P
      !                  (i-4)*(i-5)
      c(6,6,3) =       2627203._R_P/1871100._R_P
      ! stencil 4
      !                  (i+3)*i                ;                   (i+2)*i
      c(0,0,4) =       2627203._R_P/1871100._R_P; c(1,0,4) =  -359321429._R_P/19958400._R_P
      !                  (i+1)*i                ;                       i*i
      c(2,0,4) =     105706999._R_P/2217600._R_P; c(3,0,4) =  -995600723._R_P/14968800._R_P
      !                  (i-1)*i                ;                   (i-2)*i
      c(4,0,4) =     256556849._R_P/4989600._R_P; c(5,0,4) =     -15401629._R_P/739200._R_P
      !                  (i-3)*i
      c(6,0,4) =       8279479._R_P/2395008._R_P

      !                    /                    ;                   (i+1)*(i-1)
      c(0,1,4) =                          0._R_P; c(1,1,4) =     130013563._R_P/2217600._R_P
      !                      i*(i-1)            ;                   (i-1)*(i-1)
      c(2,1,4) =   -2096571887._R_P/6652800._R_P; c(3,1,4) =    2224538011._R_P/4989600._R_P
      !                  (i-2)*(i-1)            ;                   (i-3)*(i-1)
      c(4,1,4) =    -773749439._R_P/2217600._R_P; c(5,1,4) =     475321093._R_P/3326400._R_P
      !                  (i-4)*(i-1)
      c(6,1,4) =     -95508139._R_P/3991680._R_P

      !                    /                    ;                    /
      c(0,2,4) =                          0._R_P; c(1,2,4) =                          0._R_P
          !                  i*(i-1)            ;                   (i-1)*(i-2)
      c(2,2,4) =      143270957._R_P/332640._R_P; c(3,2,4) =     -412424029._R_P/332640._R_P
      !                  (i-2)*(i-2)            ;                   (i-3)*(i-2)
      c(4,2,4) =    1312114459._R_P/1330560._R_P; c(5,2,4) =   -2725575317._R_P/6652800._R_P
      !                  (i-4)*(i-2)
      c(6,2,4) =     115524053._R_P/1663200._R_P

      !                    /                    ;                     /
      c(0,3,4) =                          0._R_P; c(1,3,4) =                          0._R_P
      !                    /                    ;                   (i-1)*(i-3)
      c(2,3,4) =                          0._R_P; c(3,3,4) =    2726585359._R_P/2993760._R_P
      !                  (i-2)*(i-3)            ;                   (i-3)*(i-3)
      c(4,3,4) =    -1476618887._R_P/997920._R_P; c(5,3,4) =    1042531337._R_P/1663200._R_P
      !                  (i-4)*(i-3)
      c(6,3,4) =  -1618284323._R_P/14968800._R_P

      !                    /                    ;                     /
      c(0,4,4) =                          0._R_P; c(1,4,4) =                          0._R_P
      !                    /                    ;                     /
      c(2,4,4) =                          0._R_P; c(3,4,4) =                          0._R_P
      !                  (i-2)*(i-4)            ;                   (i-3)*(i-4)
      c(4,4,4) =        34187317._R_P/55440._R_P; c(5,4,4) =   -3573798407._R_P/6652800._R_P
      !                  (i-4)*(i-4)
      c(6,4,4) =   1894705391._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,5,4) =                          0._R_P; c(1,5,4) =                          0._R_P
      !                    /                    ;                     /
      c(2,5,4) =                          0._R_P; c(3,5,4) =                          0._R_P
      !                    /                    ;                   (i-3)*(i-5)
      c(4,5,4) =                          0._R_P; c(5,5,4) =     806338417._R_P/6652800._R_P
      !                  (i-4)*(i-5)
      c(6,5,4) =    -295455983._R_P/6652800._R_P

      !                    /                    ;                     /
      c(0,6,4) =                          0._R_P; c(1,6,4) =                          0._R_P
      !                    /                    ;                     /
      c(2,6,4) =                          0._R_P; c(3,6,4) =                          0._R_P
      !                    /                    ;                     /
      c(4,6,4) =                          0._R_P; c(5,6,4) =                          0._R_P
      !                  (i-4)*(i-5)
      c(6,6,4) =     64361771._R_P/14968800._R_P
      ! stencil 5
      !                  (i+3)*i                ;                   (i+2)*i
      c(0,0,5) =     64361771._R_P/14968800._R_P; c(1,0,5) =    -377474689._R_P/6652800._R_P
      !                  (i+1)*i                ;                       i*i
      c(2,0,5) =   3126718481._R_P/19958400._R_P; c(3,0,5) =  -3465607493._R_P/14968800._R_P
      !                  (i-1)*i                ;                   (i-2)*i
      c(4,0,5) =     320782183._R_P/1663200._R_P; c(5,0,5) =    -341910757._R_P/3991680._R_P
      !                  (i-3)*i
      c(6,0,5) =    945155329._R_P/59875200._R_P

      !                    /                    ;                   (i+1)*(i-1)
      c(0,1,5) =                          0._R_P; c(1,1,5) =    1250007643._R_P/6652800._R_P
      !                      i*(i-1)            ;                   (i-1)*(i-1)
      c(2,1,5) =   -6932480657._R_P/6652800._R_P; c(3,1,5) =      857838469._R_P/554400._R_P
      !                  (i-2)*(i-1)            ;                   (i-3)*(i-1)
      c(4,1,5) =   -8619440987._R_P/6652800._R_P; c(5,1,5) =    1924032511._R_P/3326400._R_P
      !                  (i-4)*(i-1)
      c(6,1,5) =    -712745603._R_P/6652800._R_P

      !                    /                    ;                    /
      c(0,2,5) =                          0._R_P; c(1,2,5) =                          0._R_P
          !                  i*(i-1)            ;                   (i-1)*(i-2)
      c(2,2,5) =        53678683._R_P/36960._R_P; c(3,2,5) =    -4330640057._R_P/997920._R_P
      !                  (i-2)*(i-2)            ;                   (i-3)*(i-2)
      c(4,2,5) =    4868089189._R_P/1330560._R_P; c(5,2,5) =     -405382961._R_P/246400._R_P
      !                  (i-4)*(i-2)
      c(6,2,5) =    1531307249._R_P/4989600._R_P

      !                    /                    ;                     /
      c(0,3,5) =                          0._R_P; c(1,3,5) =                          0._R_P
      !                    /                    ;                   (i-1)*(i-3)
      c(2,3,5) =                          0._R_P; c(3,3,5) =    9780057169._R_P/2993760._R_P
      !                  (i-2)*(i-3)            ;                   (i-3)*(i-3)
      c(4,3,5) =     -616410313._R_P/110880._R_P; c(5,3,5) =   12601009501._R_P/4989600._R_P
      !                  (i-4)*(i-3)
      c(6,3,5) =  -7124638253._R_P/14968800._R_P

      !                    /                    ;                     /
      c(0,4,5) =                          0._R_P; c(1,4,5) =                          0._R_P
      !                    /                    ;                     /
      c(2,4,5) =                          0._R_P; c(3,4,5) =                          0._R_P
      !                  (i-2)*(i-4)            ;                   (i-3)*(i-4)
      c(4,4,5) =      796358777._R_P/332640._R_P; c(5,4,5) =  -14684933057._R_P/6652800._R_P
      !                  (i-4)*(i-4)
      c(6,4,5) =    2811067067._R_P/6652800._R_P

      !                    /                    ;                     /
      c(0,5,5) =                          0._R_P; c(1,5,5) =                          0._R_P
      !                    /                    ;                     /
      c(2,5,5) =                          0._R_P; c(3,5,5) =                          0._R_P
      !                    /                    ;                   (i-3)*(i-5)
      c(4,5,5) =                          0._R_P; c(5,5,5) =      127942497._R_P/246400._R_P
      !                  (i-4)*(i-5)
      c(6,5,5) =  -4074544787._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,6,5) =                          0._R_P; c(1,6,5) =                          0._R_P
      !                    /                    ;                     /
      c(2,6,5) =                          0._R_P; c(3,6,5) =                          0._R_P
      !                    /                    ;                     /
      c(4,6,5) =                          0._R_P; c(5,6,5) =                          0._R_P
      !                  (i-4)*(i-5)
      c(6,6,5) =      62911297._R_P/2993760._R_P
      ! stencil 6
      !                  (i+3)*i                ;                   (i+2)*i
      c(0,0,6) =      62911297._R_P/2993760._R_P; c(1,0,6) =  -5556669277._R_P/19958400._R_P
      !                  (i+1)*i                ;                       i*i
      c(2,0,6) =  15476926351._R_P/19958400._R_P; c(3,0,6) = -17425032203._R_P/14968800._R_P
      !                  (i-1)*i                ;                   (i-2)*i
      c(4,0,6) =    4964771899._R_P/4989600._R_P; c(5,0,6) =  -9181961959._R_P/19958400._R_P
      !                  (i-3)*i
      c(6,0,6) =   5391528799._R_P/59875200._R_P

      !                    /                    ;                   (i+1)*(i-1)
      c(0,1,6) =                          0._R_P; c(1,1,6) =    2047941883._R_P/2217600._R_P
      !                      i*(i-1)            ;                   (i-1)*(i-1)
      c(2,1,6) =    -3809437823._R_P/739200._R_P; c(3,1,6) =   38683385051._R_P/4989600._R_P
      !                  (i-2)*(i-1)            ;                   (i-3)*(i-1)
      c(4,1,6) =  -14734178999._R_P/2217600._R_P; c(5,1,6) =    3417057367._R_P/1108800._R_P
      !                  (i-4)*(i-1)
      c(6,1,6) =   -2416885043._R_P/3991680._R_P

      !                    /                    ;                    /
      c(0,2,6) =                          0._R_P; c(1,2,6) =                          0._R_P
          !                  i*(i-1)            ;                   (i-1)*(i-2)
      c(2,2,6) =       199730921._R_P/27720._R_P; c(3,2,6) =   -21693002767._R_P/997920._R_P
      !                  (i-2)*(i-2)            ;                   (i-3)*(i-2)
      c(4,2,6) =     8290771913._R_P/443520._R_P; c(5,2,6) =   -19308505679._R_P/2217600._R_P
      !                  (i-4)*(i-2)
      c(6,2,6) =    8579309749._R_P/4989600._R_P

      !                    /                    ;                     /
      c(0,3,6) =                          0._R_P; c(1,3,6) =                          0._R_P
      !                    /                    ;                   (i-1)*(i-3)
      c(2,3,6) =                          0._R_P; c(3,3,6) =   49256859919._R_P/2993760._R_P
      !                  (i-2)*(i-3)            ;                   (i-3)*(i-3)
      c(4,3,6) =   -28364892607._R_P/997920._R_P; c(5,3,6) =   66440049371._R_P/4989600._R_P
      !                  (i-4)*(i-3)
      c(6,3,6) = -39645439643._R_P/14968800._R_P

      !                    /                    ;                     /
      c(0,4,6) =                          0._R_P; c(1,4,6) =                          0._R_P
      !                    /                    ;                     /
      c(2,4,6) =                          0._R_P; c(3,4,6) =                          0._R_P
      !                  (i-2)*(i-4)            ;                   (i-3)*(i-4)
      c(4,4,6) =     1369404749._R_P/110880._R_P; c(5,4,6) =    -8623431623._R_P/739200._R_P
      !                  (i-4)*(i-4)
      c(6,4,6) =  46808583631._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,5,6) =                          0._R_P; c(1,5,6) =                          0._R_P
      !                    /                    ;                     /
      c(2,5,6) =                          0._R_P; c(3,5,6) =                          0._R_P
      !                    /                    ;                   (i-3)*(i-5)
      c(4,5,6) =                          0._R_P; c(5,5,6) =    6182612731._R_P/2217600._R_P
      !                  (i-4)*(i-5)
      c(6,5,6) = -22763092357._R_P/19958400._R_P

      !                    /                    ;                     /
      c(0,6,6) =                          0._R_P; c(1,6,6) =                          0._R_P
      !                    /                    ;                     /
      c(2,6,6) =                          0._R_P; c(3,6,6) =                          0._R_P
      !                    /                    ;                     /
      c(4,6,6) =                          0._R_P; c(5,6,6) =                          0._R_P
      !                  (i-4)*(i-5)
      c(6,6,6) =     897207163._R_P/7484400._R_P
    case(8) ! 15th order
      ! stencil 0
      !                    /                              ;                      /
      c(0,0,0) =     5870785406797._R_P/  20756736000._R_P; c(1,0,0) =    -3130718954431._R_P/    972972000._R_P
      !                    /                              ;                      /
      c(2,0,0) =    36019630238453._R_P/   4447872000._R_P; c(3,0,0) =    -9030771744409._R_P/    778377600._R_P
      !                    /                              ;                      /
      c(4,0,0) =     7028987165449._R_P/    691891200._R_P; c(5,0,0) =    -5269260407953._R_P/    972972000._R_P
      !                    /                              ;                      /
      c(6,0,0) =    50528822994577._R_P/  31135104000._R_P; c(7,0,0) =     -819100494587._R_P/   3891888000._R_P

      !                    /                              ;                      /
      c(0,1,0) =                                    0._R_P; c(1,1,0) =   581791881407369._R_P/  62270208000._R_P
      !                    /                              ;                      /
      c(2,1,0) =  -185432400549349._R_P/   3891888000._R_P; c(3,1,0) =   428668917728281._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(4,1,0) =     -393303816739._R_P/      6486480._R_P; c(5,1,0) =  1010731494899387._R_P/  31135104000._R_P
      !                    /                              ;                      /
      c(6,1,0) =   -12661520644021._R_P/   1297296000._R_P; c(7,1,0) =    39509061792127._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,2,0) =                                    0._R_P; c(1,2,0) =                                    0._R_P
      !                    /                              ;                      /
      c(2,2,0) =  1272280750118197._R_P/  20756736000._R_P; c(3,2,0) =    -6306477584539._R_P/     35380800._R_P
      !                    /                              ;                      /
      c(4,2,0) =   982150494698309._R_P/   6227020800._R_P; c(5,2,0) =  -109928049802589._R_P/   1297296000._R_P
      !                    /                              ;                      /
      c(6,2,0) =   795325997722517._R_P/  31135104000._R_P; c(7,2,0) =    -6476591199161._R_P/   1945944000._R_P

      !                    /                              ;                      /
      c(0,3,0) =                                    0._R_P; c(1,3,0) =                                    0._R_P
      !                    /                              ;                      /
      c(2,3,0) =                                    0._R_P; c(3,3,0) =     5896382977423._R_P/     45287424._R_P
      !                    /                              ;                      /
      c(4,3,0) =   -35999233471051._R_P/    155675520._R_P; c(5,3,0) =   775760249154827._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(6,3,0) =    -4882688924777._R_P/    129729600._R_P; c(7,3,0) =    10196716797013._R_P/   2075673600._R_P

      !                    /                              ;                      /
      c(0,4,0) =                                    0._R_P; c(1,4,0) =                                    0._R_P
      !                    /                              ;                      /
      c(2,4,0) =                                    0._R_P; c(3,4,0) =                                    0._R_P
      !                    /                              ;                      /
      c(4,4,0) =    23315424178373._R_P/    226437120._R_P; c(5,4,0) =     -983492927359._R_P/      8845200._R_P
      !                    /                              ;                      /
      c(6,4,0) =    41910140004779._R_P/   1245404160._R_P; c(7,4,0) =    -3423798156193._R_P/    778377600._R_P

      !                    /                              ;                      /
      c(0,5,0) =                                    0._R_P; c(1,5,0) =                                    0._R_P
      !                    /                              ;                      /
      c(2,5,0) =                                    0._R_P; c(3,5,0) =                                    0._R_P
      !                    /                              ;                      /
      c(4,5,0) =                                    0._R_P; c(5,5,0) =   624177436330267._R_P/  20756736000._R_P
      !                    /                              ;                      /
      c(6,5,0) =   -70944310593109._R_P/   3891888000._R_P; c(7,5,0) =    10610581100123._R_P/   4447872000._R_P

      !                    /                              ;                      /
      c(0,6,0) =                                    0._R_P; c(1,6,0) =                                    0._R_P
      !                    /                              ;                      /
      c(2,6,0) =                                    0._R_P; c(3,6,0) =                                    0._R_P
      !                    /                              ;                      /
      c(4,6,0) =                                    0._R_P; c(5,6,0) =                                    0._R_P
      !                    /                              ;                      /
      c(6,6,0) =   172229708657639._R_P/  62270208000._R_P; c(7,6,0) =    -1410106709147._R_P/   1945944000._R_P

      !                    /                              ;                      /
      c(0,7,0) =                                    0._R_P; c(1,7,0) =                                    0._R_P
      !                    /                              ;                      /
      c(2,7,0) =                                    0._R_P; c(3,7,0) =                                    0._R_P
      !                    /                              ;                      /
      c(4,7,0) =                                    0._R_P; c(5,7,0) =                                    0._R_P
      !                    /                              ;                      /
      c(6,7,0) =   172229708657639._R_P/  62270208000._R_P; c(7,7,0) =      986005096387._R_P/  20756736000._R_P

      ! stencil 1
      !                    /                              ;                      /
      c(0,0,1) =      986005096387._R_P/  20756736000._R_P; c(1,0,1) =    -1069457397287._R_P/   1945944000._R_P
      !                    /                              ;                      /
      c(2,0,1) =    43315366304381._R_P/  31135104000._R_P; c(3,0,1) =    -1550584925161._R_P/    778377600._R_P
      !                    /                              ;                      /
      c(4,0,1) =      721470910481._R_P/    415134720._R_P; c(5,0,1) =    -1793558121581._R_P/   1945944000._R_P
      !                    /                              ;                      /
      c(6,0,1) =     1221480056521._R_P/   4447872000._R_P; c(7,0,1) =     -137801870867._R_P/   3891888000._R_P

      !                    /                              ;                      /
      c(0,1,1) =                                    0._R_P; c(1,1,1) =   102080471419559._R_P/  62270208000._R_P
      !                    /                              ;                      /
      c(2,1,1) =   -32903428273669._R_P/   3891888000._R_P; c(3,1,1) =    76273513229143._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(4,1,1) =    -1397571412901._R_P/    129729600._R_P; c(5,1,1) =   178922840432597._R_P/  31135104000._R_P
      !                    /                              ;                      /
      c(6,1,1) =    -2230862726341._R_P/   1297296000._R_P; c(7,1,1) =     6925711076497._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,2,1) =                                    0._R_P; c(1,2,1) =                                    0._R_P
      !                    /                              ;                      /
      c(2,2,1) =   229456135916827._R_P/  20756736000._R_P; c(3,2,1) =      -11450077957._R_P/       353808._R_P
      !                    /                              ;                      /
      c(4,2,1) =   178559835040523._R_P/   6227020800._R_P; c(5,2,1) =   -19952704102349._R_P/   1297296000._R_P
      !                    /                              ;                      /
      c(6,2,1) =   143887855797947._R_P/  31135104000._R_P; c(7,2,1) =     -583488131053._R_P/    972972000._R_P

      !                    /                              ;                      /
      c(0,3,1) =                                    0._R_P; c(1,3,1) =                                    0._R_P
      !                    /                              ;                      /
      c(2,3,1) =                                    0._R_P; c(3,3,1) =     5407733702789._R_P/    226437120._R_P
      !                    /                              ;                      /
      c(4,3,1) =    -6630479776771._R_P/    155675520._R_P; c(5,3,1) =   142950967195973._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(6,3,1) =     -224563041869._R_P/     32432400._R_P; c(7,3,1) =       41566759079._R_P/     46126080._R_P

      !                    /                              ;                      /
      c(0,4,1) =                                    0._R_P; c(1,4,1) =                                    0._R_P
      !                    /                              ;                      /
      c(2,4,1) =                                    0._R_P; c(3,4,1) =                                    0._R_P
      !                    /                              ;                      /
      c(4,4,1) =     4322531771339._R_P/    226437120._R_P; c(5,4,1) =      -29244985495._R_P/      1415232._R_P
      !                    /                              ;                      /
      c(6,4,1) =    38941083744793._R_P/   6227020800._R_P; c(7,4,1) =      -90744192823._R_P/    111196800._R_P

      !                    /                              ;                      /
      c(0,5,1) =                                    0._R_P; c(1,5,1) =                                    0._R_P
      !                    /                              ;                      /
      c(2,5,1) =                                    0._R_P; c(3,5,1) =                                    0._R_P
      !                    /                              ;                      /
      c(4,5,1) =                                    0._R_P; c(5,5,1) =   116487285372277._R_P/  20756736000._R_P
      !                    /                              ;                      /
      c(6,5,1) =   -13257668940469._R_P/   3891888000._R_P; c(7,5,1) =    13873328286131._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,6,1) =                                    0._R_P; c(1,6,1) =                                    0._R_P
      !                    /                              ;                      /
      c(2,6,1) =                                    0._R_P; c(3,6,1) =                                    0._R_P
      !                    /                              ;                      /
      c(4,6,1) =                                    0._R_P; c(5,6,1) =                                    0._R_P
      !                    /                              ;                      /
      c(6,6,1) =    32268504444809._R_P/  62270208000._R_P; c(7,6,1) =     -132173819131._R_P/    972972000._R_P

      !                    /                              ;                      /
      c(0,7,1) =                                    0._R_P; c(1,7,1) =                                    0._R_P
      !                    /                              ;                      /
      c(2,7,1) =                                    0._R_P; c(3,7,1) =                                    0._R_P
      !                    /                              ;                      /
      c(4,7,1) =                                    0._R_P; c(5,7,1) =                                    0._R_P
      !                    /                              ;                      /
      c(6,7,1) =                                    0._R_P; c(7,7,1) =       26446172491._R_P/   2965248000._R_P

      ! stencil 2
      !                    /                              ;                      /
      c(0,0,2) =       26446172491._R_P/   2965248000._R_P; c(1,0,2) =     -104391937861._R_P/    972972000._R_P
      !                    /                              ;                      /
      c(2,0,2) =     8624638348211._R_P/  31135104000._R_P; c(3,0,2) =     -310726966393._R_P/    778377600._R_P
      !                    /                              ;                      /
      c(4,0,2) =      721220745563._R_P/   2075673600._R_P; c(5,0,2) =      -25412164549._R_P/    138996000._R_P
      !                    /                              ;                      /
      c(6,0,2) =     1677021138577._R_P/  31135104000._R_P; c(7,0,2) =      -26674345787._R_P/   3891888000._R_P

      !                    /                              ;                      /
      c(0,1,2) =                                    0._R_P; c(1,1,2) =    20863031646089._R_P/  62270208000._R_P
      !                    /                              ;                      /
      c(2,1,2) =    -6905100758509._R_P/   3891888000._R_P; c(3,1,2) =     3240510296069._R_P/   1245404160._R_P
      !                    /                              ;                      /
      c(4,1,2) =      -37187936869._R_P/     16216200._R_P; c(5,1,2) =    37913679009467._R_P/  31135104000._R_P
      !                    /                              ;                      /
      c(6,1,2) =     -468561665821._R_P/   1297296000._R_P; c(7,1,2) =     1438198790527._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,2,2) =                                    0._R_P; c(1,2,2) =                                    0._R_P
      !                    /                              ;                      /
      c(2,2,2) =    49883478342517._R_P/  20756736000._R_P; c(3,2,2) =     -253865691211._R_P/     35380800._R_P
      !                    /                              ;                      /
      c(4,2,2) =    39896100785477._R_P/   6227020800._R_P; c(5,2,2) =    -4456767285989._R_P/   1297296000._R_P
      !                    /                              ;                      /
      c(6,2,2) =    31959522170837._R_P/  31135104000._R_P; c(7,2,2) =     -256879392281._R_P/   1945944000._R_P

      !                    /                              ;                      /
      c(0,3,2) =                                    0._R_P; c(1,3,2) =                                    0._R_P
      !                    /                              ;                      /
      c(2,3,2) =                                    0._R_P; c(3,3,2) =     1231949387723._R_P/    226437120._R_P
      !                    /                              ;                      /
      c(4,3,2) =    -1532094364651._R_P/    155675520._R_P; c(5,3,2) =    33191727291659._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(6,3,2) =      -41643930661._R_P/     25945920._R_P; c(7,3,2) =      431000077397._R_P/   2075673600._R_P

      !                    /                              ;                      /
      c(0,4,2) =                                    0._R_P; c(1,4,2) =                                    0._R_P
      !                    /                              ;                      /
      c(2,4,2) =                                    0._R_P; c(3,4,2) =                                    0._R_P
      !                    /                              ;                      /
      c(4,4,2) =      203912134273._R_P/     45287424._R_P; c(5,4,2) =       -5445142127._R_P/      1105650._R_P
      !                    /                              ;                      /
      c(6,4,2) =     9306913817431._R_P/   6227020800._R_P; c(7,4,2) =     -151441370209._R_P/    778377600._R_P

      !                    /                              ;                      /
      c(0,5,2) =                                    0._R_P; c(1,5,2) =                                    0._R_P
      !                    /                              ;                      /
      c(2,5,2) =                                    0._R_P; c(3,5,2) =                                    0._R_P
      !                    /                              ;                      /
      c(4,5,2) =                                    0._R_P; c(5,5,2) =    28199161918747._R_P/  20756736000._R_P
      !                    /                              ;                      /
      c(6,5,2) =    -3233549114749._R_P/   3891888000._R_P; c(7,5,2) =     3388533713021._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,6,2) =                                    0._R_P; c(1,6,2) =                                    0._R_P
      !                    /                              ;                      /
      c(2,6,2) =                                    0._R_P; c(3,6,2) =                                    0._R_P
      !                    /                              ;                      /
      c(4,6,2) =                                    0._R_P; c(5,6,2) =                                    0._R_P
      !                    /                              ;                      /
      c(6,6,2) =     7965255985319._R_P/  62270208000._R_P; c(7,6,2) =      -65611168187._R_P/   1945944000._R_P

      !                    /                              ;                      /
      c(0,7,2) =                                    0._R_P; c(1,7,2) =                                    0._R_P
      !                    /                              ;                      /
      c(2,7,2) =                                    0._R_P; c(3,7,2) =                                    0._R_P
      !                    /                              ;                      /
      c(4,7,2) =                                    0._R_P; c(5,7,2) =                                    0._R_P
      !                    /                              ;                      /
      c(6,7,2) =                                    0._R_P; c(7,7,2) =       46388292547._R_P/  20756736000._R_P

      ! stencil 3
      !                    /                              ;                      /
      c(0,0,3) =       46388292547._R_P/  20756736000._R_P; c(1,0,3) =      -56245265927._R_P/   1945944000._R_P
      !                    /                              ;                      /
      c(2,0,3) =     2458417783421._R_P/  31135104000._R_P; c(3,0,3) =      -18415814357._R_P/    155675520._R_P
      !                    /                              ;                      /
      c(4,0,3) =       72812006087._R_P/    691891200._R_P; c(5,0,3) =     -108473646221._R_P/   1945944000._R_P
      !                    /                              ;                      /
      c(6,0,3) =      508082860927._R_P/  31135104000._R_P; c(7,0,3) =       -7942541267._R_P/   3891888000._R_P

      !                    /                              ;                      /
      c(0,1,3) =                                    0._R_P; c(1,1,3) =     6047605530599._R_P/  62270208000._R_P
      !                    /                              ;                      /
      c(2,1,3) =    -2129103852829._R_P/   3891888000._R_P; c(3,1,3) =     5227966881367._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(4,1,3) =      -98765696693._R_P/    129729600._R_P; c(5,1,3) =    12752830987157._R_P/  31135104000._R_P
      !                    /                              ;                      /
      c(6,1,3) =     -157580595421._R_P/   1297296000._R_P; c(7,1,3) =      478185649297._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,2,3) =                                    0._R_P; c(1,2,3) =                                    0._R_P
      !                    /                              ;                      /
      c(2,2,3) =    16476387815707._R_P/  20756736000._R_P; c(3,2,3) =       -5527715497._R_P/      2211300._R_P
      !                    /                              ;                      /
      c(4,2,3) =    14416393946891._R_P/   6227020800._R_P; c(5,2,3) =    -1644079167749._R_P/   1297296000._R_P
      !                    /                              ;                      /
      c(6,2,3) =    11870432980667._R_P/  31135104000._R_P; c(7,2,3) =      -47469340603._R_P/    972972000._R_P

      !                    /                              ;                      /
      c(0,3,3) =                                    0._R_P; c(1,3,3) =                                    0._R_P
      !                    /                              ;                      /
      c(2,3,3) =                                    0._R_P; c(3,3,3) =      457249528517._R_P/    226437120._R_P
      !                    /                              ;                      /
      c(4,3,3) =     -595915721251._R_P/    155675520._R_P; c(5,3,3) =      532071643661._R_P/    249080832._R_P
      !                    /                              ;                      /
      c(6,3,3) =      -10590149653._R_P/     16216200._R_P; c(7,3,3) =       25116366157._R_P/    296524800._R_P

      !                    /                              ;                      /
      c(0,4,3) =                                    0._R_P; c(1,4,3) =                                    0._R_P
      !                    /                              ;                      /
      c(2,4,3) =                                    0._R_P; c(3,4,3) =                                    0._R_P
      !                    /                              ;                      /
      c(4,4,3) =      420341161931._R_P/    226437120._R_P; c(5,4,3) =      -74851467823._R_P/     35380800._R_P
      !                    /                              ;                      /
      c(6,4,3) =     4100880843289._R_P/   6227020800._R_P; c(7,4,3) =      -67513265377._R_P/    778377600._R_P

      !                    /                              ;                      /
      c(0,5,3) =                                    0._R_P; c(1,5,3) =                                    0._R_P
      !                    /                              ;                      /
      c(2,5,3) =                                    0._R_P; c(3,5,3) =                                    0._R_P
      !                    /                              ;                      /
      c(4,5,3) =                                    0._R_P; c(5,5,3) =    12780967457077._R_P/  20756736000._R_P
      !                    /                              ;                      /
      c(6,5,3) =    -1521688484269._R_P/   3891888000._R_P; c(7,5,3) =     1631589107891._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,6,3) =                                    0._R_P; c(1,6,3) =                                    0._R_P
      !                    /                              ;                      /
      c(2,6,3) =                                    0._R_P; c(3,6,3) =                                    0._R_P
      !                    /                              ;                      /
      c(4,6,3) =                                    0._R_P; c(5,6,3) =                                    0._R_P
      !                    /                              ;                      /
      c(6,6,3) =     3944861897609._R_P/  62270208000._R_P; c(7,6,3) =       -2407377043._R_P/    138996000._R_P

      !                    /                              ;                      /
      c(0,7,3) =                                    0._R_P; c(1,7,3) =                                    0._R_P
      !                    /                              ;                      /
      c(2,7,3) =                                    0._R_P; c(3,7,3) =                                    0._R_P
      !                    /                              ;                      /
      c(4,7,3) =                                    0._R_P; c(5,7,3) =                                    0._R_P
      !                    /                              ;                      /
      c(6,7,3) =                                    0._R_P; c(7,7,3) =       25116366157._R_P/  20756736000._R_P

      ! stencil 4
      !                    /                              ;                      /
      c(0,0,4) =       25116366157._R_P/  20756736000._R_P; c(1,0,4) =       -2407377043._R_P/    138996000._R_P
      !                    /                              ;                      /
      c(2,0,4) =     1631589107891._R_P/  31135104000._R_P; c(3,0,4) =      -67513265377._R_P/    778377600._R_P
      !                    /                              ;                      /
      c(4,0,4) =       25116366157._R_P/    296524800._R_P; c(5,0,4) =      -47469340603._R_P/    972972000._R_P
      !                    /                              ;                      /
      c(6,0,4) =      478185649297._R_P/  31135104000._R_P; c(7,0,4) =       -7942541267._R_P/   3891888000._R_P

      !                    /                              ;                      /
      c(0,1,4) =                                    0._R_P; c(1,1,4) =     3944861897609._R_P/  62270208000._R_P
      !                    /                              ;                      /
      c(2,1,4) =    -1521688484269._R_P/   3891888000._R_P; c(3,1,4) =     4100880843289._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(4,1,4) =      -10590149653._R_P/     16216200._R_P; c(5,1,4) =    11870432980667._R_P/  31135104000._R_P
      !                    /                              ;                      /
      c(6,1,4) =     -157580595421._R_P/   1297296000._R_P; c(7,1,4) =      508082860927._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,2,4) =                                    0._R_P; c(1,2,4) =                                    0._R_P
      !                    /                              ;                      /
      c(2,2,4) =    12780967457077._R_P/  20756736000._R_P; c(3,2,4) =      -74851467823._R_P/     35380800._R_P
      !                    /                              ;                      /
      c(4,2,4) =      532071643661._R_P/    249080832._R_P; c(5,2,4) =    -1644079167749._R_P/   1297296000._R_P
      !                    /                              ;                      /
      c(6,2,4) =    12752830987157._R_P/  31135104000._R_P; c(7,2,4) =     -108473646221._R_P/   1945944000._R_P

      !                    /                              ;                      /
      c(0,3,4) =                                    0._R_P; c(1,3,4) =                                    0._R_P
      !                    /                              ;                      /
      c(2,3,4) =                                    0._R_P; c(3,3,4) =      420341161931._R_P/    226437120._R_P
      !                    /                              ;                      /
      c(4,3,4) =     -595915721251._R_P/    155675520._R_P; c(5,3,4) =    14416393946891._R_P/   6227020800._R_P
      !                    /                              ;                      /
      c(6,3,4) =      -98765696693._R_P/    129729600._R_P; c(7,3,4) =       72812006087._R_P/    691891200._R_P

      !                    /                              ;                      /
      c(0,4,4) =                                    0._R_P; c(1,4,4) =                                    0._R_P
      !                    /                              ;                      /
      c(2,4,4) =                                    0._R_P; c(3,4,4) =                                    0._R_P
      !                    /                              ;                      /
      c(4,4,4) =      457249528517._R_P/    226437120._R_P; c(5,4,4) =       -5527715497._R_P/      2211300._R_P
      !                    /                              ;                      /
      c(6,4,4) =     5227966881367._R_P/   6227020800._R_P; c(7,4,4) =      -18415814357._R_P/    155675520._R_P

      !                    /                              ;                      /
      c(0,5,4) =                                    0._R_P; c(1,5,4) =                                    0._R_P
      !                    /                              ;                      /
      c(2,5,4) =                                    0._R_P; c(3,5,4) =                                    0._R_P
      !                    /                              ;                      /
      c(4,5,4) =                                    0._R_P; c(5,5,4) =    16476387815707._R_P/  20756736000._R_P
      !                    /                              ;                      /
      c(6,5,4) =    -2129103852829._R_P/   3891888000._R_P; c(7,5,4) =     2458417783421._R_P/  31135104000._R_P

      !                    /                              ;                      /
      c(0,6,4) =                                    0._R_P; c(1,6,4) =                                    0._R_P
      !                    /                              ;                      /
      c(2,6,4) =                                    0._R_P; c(3,6,4) =                                    0._R_P
      !                    /                              ;                      /
      c(4,6,4) =                                    0._R_P; c(5,6,4) =                                    0._R_P
      !                    /                              ;                      /
      c(6,6,4) =     6047605530599._R_P/  62270208000._R_P; c(7,6,4) =      -56245265927._R_P/   1945944000._R_P

      !                    /                              ;                      /
      c(0,7,4) =                                    0._R_P; c(1,7,4) =                                    0._R_P
      !                    /                              ;                      /
      c(2,7,4) =                                    0._R_P; c(3,7,4) =                                    0._R_P
      !                    /                              ;                      /
      c(4,7,4) =                                    0._R_P; c(5,7,4) =                                    0._R_P
      !                    /                              ;                      /
      c(6,7,4) =                                    0._R_P; c(7,7,4) =       46388292547._R_P/  20756736000._R_P
      ! stencil 5
      !                    /                              ;                     /
      c(0,0,5) =       46388292547._R_P/  20756736000._R_P; c(1,0,5) =      -65611168187._R_P/   1945944000._R_P
      !                    /                              ;                     /
      c(2,0,5) =     3388533713021._R_P/  31135104000._R_P; c(3,0,5) =     -151441370209._R_P/    778377600._R_P
      !                    /                              ;                     /
      c(4,0,5) =      431000077397._R_P/   2075673600._R_P; c(5,0,5) =     -256879392281._R_P/   1945944000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,0,5) =     1438198790527._R_P/  31135104000._R_P; c(7,0,5) =      -26674345787._R_P/   3891888000._R_P

      !                    /                              ;                     /
      c(0,1,5) =                                    0._R_P; c(1,1,5) =     7965255985319._R_P/  62270208000._R_P
      !                    /                              ;                     /
      c(2,1,5) =    -3233549114749._R_P/   3891888000._R_P; c(3,1,5) =     9306913817431._R_P/   6227020800._R_P
      !                    /                              ;                     /
      c(4,1,5) =      -41643930661._R_P/     25945920._R_P; c(5,1,5) =    31959522170837._R_P/  31135104000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,1,5) =     -468561665821._R_P/   1297296000._R_P; c(7,1,5) =     1677021138577._R_P/  31135104000._R_P

      !                    /                              ;                     /
      c(0,2,5) =                                    0._R_P; c(1,2,5) =                                    0._R_P
      !                    /                              ;                     /
      c(2,2,5) =    28199161918747._R_P/  20756736000._R_P; c(3,2,5) =       -5445142127._R_P/      1105650._R_P
      !                    /                              ;                     /
      c(4,2,5) =    33191727291659._R_P/   6227020800._R_P; c(5,2,5) =    -4456767285989._R_P/   1297296000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,2,5) =    37913679009467._R_P/  31135104000._R_P; c(7,2,5) =      -25412164549._R_P/    138996000._R_P

      !                    /                              ;                     /
      c(0,3,5) =                                    0._R_P; c(1,3,5) =                                    0._R_P
      !                    /                              ;                     /
      c(2,3,5) =                                    0._R_P; c(3,3,5) =      203912134273._R_P/     45287424._R_P
      !                    /                              ;                     /
      c(4,3,5) =    -1532094364651._R_P/    155675520._R_P; c(5,3,5) =    39896100785477._R_P/   6227020800._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,3,5) =      -37187936869._R_P/     16216200._R_P; c(7,3,5) =      721220745563._R_P/   2075673600._R_P

      !                    /                              ;                     /
      c(0,4,5) =                                    0._R_P; c(1,4,5) =                                    0._R_P
      !                    /                              ;                     /
      c(2,4,5) =                                    0._R_P; c(3,4,5) =                                    0._R_P
      !                    /                              ;                     /
      c(4,4,5) =     1231949387723._R_P/    226437120._R_P; c(5,4,5) =     -253865691211._R_P/     35380800._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,4,5) =     3240510296069._R_P/   1245404160._R_P; c(7,4,5) =     -310726966393._R_P/    778377600._R_P

      !                    /                              ;                     /
      c(0,5,5) =                                    0._R_P; c(1,5,5) =                                    0._R_P
      !                    /                              ;                     /
      c(2,5,5) =                                    0._R_P; c(3,5,5) =                                    0._R_P
      !                    /                              ;                     /
      c(4,5,5) =                                    0._R_P; c(5,5,5) =    49883478342517._R_P/  20756736000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,5,5) =    -6905100758509._R_P/   3891888000._R_P; c(7,5,5) =     8624638348211._R_P/  31135104000._R_P

      !                    /                              ;                     /
      c(0,6,5) =                                    0._R_P; c(1,6,5) =                                    0._R_P
      !                    /                              ;                     /
      c(2,6,5) =                                    0._R_P; c(3,6,5) =                                    0._R_P
      !                    /                              ;                     /
      c(4,6,5) =                                    0._R_P; c(5,6,5) =                                    0._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,6,5) =    20863031646089._R_P/  62270208000._R_P; c(7,6,5) =     -104391937861._R_P/    972972000._R_P

      !                    /                              ;                     /
      c(0,7,5) =                                    0._R_P; c(1,7,5) =                                    0._R_P
      !                    /                              ;                     /
      c(2,7,5) =                                    0._R_P; c(3,7,5) =                                    0._R_P
      !                    /                              ;                     /
      c(4,7,5) =                                    0._R_P; c(5,7,5) =                                    0._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,7,5) =                                    0._R_P; c(7,7,5) =       26446172491._R_P/   2965248000._R_P

      ! stencil 6
      !                    /                              ;                     /
      c(0,0,6) =       26446172491._R_P/   2965248000._R_P; c(1,0,6) =     -132173819131._R_P/    972972000._R_P
      !                    /                              ;                     /
      c(2,0,6) =    13873328286131._R_P/  31135104000._R_P; c(3,0,6) =      -90744192823._R_P/    111196800._R_P
      !                    /                              ;                     /
      c(4,0,6) =       41566759079._R_P/     46126080._R_P; c(5,0,6) =     -583488131053._R_P/    972972000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,0,6) =     6925711076497._R_P/  31135104000._R_P; c(7,0,6) =     -137801870867._R_P/   3891888000._R_P

      !                    /                              ;                     /
      c(0,1,6) =                                    0._R_P; c(1,1,6) =    32268504444809._R_P/  62270208000._R_P
      !                    /                              ;                     /
      c(2,1,6) =   -13257668940469._R_P/   3891888000._R_P; c(3,1,6) =    38941083744793._R_P/   6227020800._R_P
      !                    /                              ;                     /
      c(4,1,6) =     -224563041869._R_P/     32432400._R_P; c(5,1,6) =   143887855797947._R_P/  31135104000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,1,6) =    -2230862726341._R_P/   1297296000._R_P; c(7,1,6) =     1221480056521._R_P/   4447872000._R_P

      !                    /                              ;                     /
      c(0,2,6) =                                    0._R_P; c(1,2,6) =                                    0._R_P
      !                    /                              ;                     /
      c(2,2,6) =   116487285372277._R_P/  20756736000._R_P; c(3,2,6) =      -29244985495._R_P/      1415232._R_P
      !                    /                              ;                     /
      c(4,2,6) =   142950967195973._R_P/   6227020800._R_P; c(5,2,6) =   -19952704102349._R_P/   1297296000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,2,6) =   178922840432597._R_P/  31135104000._R_P; c(7,2,6) =    -1793558121581._R_P/   1945944000._R_P

      !                    /                              ;                     /
      c(0,3,6) =                                    0._R_P; c(1,3,6) =                                    0._R_P
      !                    /                              ;                     /
      c(2,3,6) =                                    0._R_P; c(3,3,6) =     4322531771339._R_P/    226437120._R_P
      !                    /                              ;                     /
      c(4,3,6) =    -6630479776771._R_P/    155675520._R_P; c(5,3,6) =   178559835040523._R_P/   6227020800._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,3,6) =    -1397571412901._R_P/    129729600._R_P; c(7,3,6) =      721470910481._R_P/    415134720._R_P

      !                    /                              ;                     /
      c(0,4,6) =                                    0._R_P; c(1,4,6) =                                    0._R_P
      !                    /                              ;                     /
      c(2,4,6) =                                    0._R_P; c(3,4,6) =                                    0._R_P
      !                    /                              ;                     /
      c(4,4,6) =     5407733702789._R_P/    226437120._R_P; c(5,4,6) =      -11450077957._R_P/       353808._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,4,6) =    76273513229143._R_P/   6227020800._R_P; c(7,4,6) =    -1550584925161._R_P/    778377600._R_P

      !                    /                              ;                     /
      c(0,5,6) =                                    0._R_P; c(1,5,6) =                                    0._R_P
      !                    /                              ;                     /
      c(2,5,6) =                                    0._R_P; c(3,5,6) =                                    0._R_P
      !                    /                              ;                     /
      c(4,5,6) =                                    0._R_P; c(5,5,6) =   229456135916827._R_P/  20756736000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,5,6) =   -32903428273669._R_P/   3891888000._R_P; c(7,5,6) =    43315366304381._R_P/  31135104000._R_P

      !                    /                              ;                     /
      c(0,6,6) =                                    0._R_P; c(1,6,6) =                                    0._R_P
      !                    /                              ;                     /
      c(2,6,6) =                                    0._R_P; c(3,6,6) =                                    0._R_P
      !                    /                              ;                     /
      c(4,6,6) =                                    0._R_P; c(5,6,6) =                                    0._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,6,6) =   102080471419559._R_P/  62270208000._R_P; c(7,6,6) =    -1069457397287._R_P/   1945944000._R_P

      !                    /                              ;                     /
      c(0,7,6) =                                    0._R_P; c(1,7,6) =                                    0._R_P
      !                    /                              ;                     /
      c(2,7,6) =                                    0._R_P; c(3,7,6) =                                    0._R_P
      !                    /                              ;                     /
      c(4,7,6) =                                    0._R_P; c(5,7,6) =                                    0._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,7,6) =                                    0._R_P; c(7,7,6) =      986005096387._R_P/  20756736000._R_P

      ! stencil 7
      !                    /                              ;                     /
      c(0,0,7) =      986005096387._R_P/  20756736000._R_P; c(1,0,7) =    -1410106709147._R_P/   1945944000._R_P
      !                    /                              ;                     /
      c(2,0,7) =    10610581100123._R_P/   4447872000._R_P; c(3,0,7) =    -3423798156193._R_P/    778377600._R_P
      !                    /                              ;                     /
      c(4,0,7) =    10196716797013._R_P/   2075673600._R_P; c(5,0,7) =    -6476591199161._R_P/   1945944000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,0,7) =    39509061792127._R_P/  31135104000._R_P; c(7,0,7) =     -819100494587._R_P/   3891888000._R_P

      !                    /                              ;                     /
      c(0,1,7) =                                    0._R_P; c(1,1,7) =   172229708657639._R_P/  62270208000._R_P
      !                    /                              ;                     /
      c(2,1,7) =   -70944310593109._R_P/   3891888000._R_P; c(3,1,7) =    41910140004779._R_P/   1245404160._R_P
      !                    /                              ;                     /
      c(4,1,7) =    -4882688924777._R_P/    129729600._R_P; c(5,1,7) =   795325997722517._R_P/  31135104000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,1,7) =   -12661520644021._R_P/   1297296000._R_P; c(7,1,7) =    50528822994577._R_P/  31135104000._R_P

      !                    /                              ;                     /
      c(0,2,7) =                                    0._R_P; c(1,2,7) =                                    0._R_P
      !                    /                              ;                     /
      c(2,2,7) =   624177436330267._R_P/  20756736000._R_P; c(3,2,7) =     -983492927359._R_P/      8845200._R_P
      !                    /                              ;                     /
      c(4,2,7) =   775760249154827._R_P/   6227020800._R_P; c(5,2,7) =  -109928049802589._R_P/   1297296000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,2,7) =  1010731494899387._R_P/  31135104000._R_P; c(7,2,7) =    -5269260407953._R_P/    972972000._R_P

      !                    /                              ;                     /
      c(0,3,7) =                                    0._R_P; c(1,3,7) =                                    0._R_P
      !                    /                              ;                     /
      c(2,3,7) =                                    0._R_P; c(3,3,7) =    23315424178373._R_P/    226437120._R_P
      !                    /                              ;                     /
      c(4,3,7) =   -35999233471051._R_P/    155675520._R_P; c(5,3,7) =   982150494698309._R_P/   6227020800._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,3,7) =     -393303816739._R_P/      6486480._R_P; c(7,3,7) =     7028987165449._R_P/    691891200._R_P

      !                    /                              ;                     /
      c(0,4,7) =                                    0._R_P; c(1,4,7) =                                    0._R_P
      !                    /                              ;                     /
      c(2,4,7) =                                    0._R_P; c(3,4,7) =                                    0._R_P
      !                    /                              ;                     /
      c(4,4,7) =     5896382977423._R_P/     45287424._R_P; c(5,4,7) =    -6306477584539._R_P/     35380800._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,4,7) =   428668917728281._R_P/   6227020800._R_P; c(7,4,7) =    -9030771744409._R_P/    778377600._R_P

      !                    /                              ;                     /
      c(0,5,7) =                                    0._R_P; c(1,5,7) =                                    0._R_P
      !                    /                              ;                     /
      c(2,5,7) =                                    0._R_P; c(3,5,7) =                                    0._R_P
      !                    /                              ;                     /
      c(4,5,7) =                                    0._R_P; c(5,5,7) =  1272280750118197._R_P/  20756736000._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,5,7) =  -185432400549349._R_P/   3891888000._R_P; c(7,5,7) =    36019630238453._R_P/   4447872000._R_P

      !                    /                              ;                     /
      c(0,6,7) =                                    0._R_P; c(1,6,7) =                                    0._R_P
      !                    /                              ;                     /
      c(2,6,7) =                                    0._R_P; c(3,6,7) =                                    0._R_P
      !                    /                              ;                     /
      c(4,6,7) =                                    0._R_P; c(5,6,7) =                                    0._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,6,7) =   581791881407369._R_P/  62270208000._R_P; c(7,6,7) =    -3130718954431._R_P/    972972000._R_P

      !                    /                              ;                     /
      c(0,7,7) =                                    0._R_P; c(1,7,7) =                                    0._R_P
      !                    /                              ;                     /
      c(2,7,7) =                                    0._R_P; c(3,7,7) =                                    0._R_P
      !                    /                              ;                     /
      c(4,7,7) =                                    0._R_P; c(5,7,7) =                                    0._R_P
      !                  (i-4)*(i-5)                      ;
      c(6,7,7) =                                    0._R_P; c(7,7,7) =     5870785406797._R_P/  20756736000._R_P
    case(9) ! 17th order
      ! stencil 0
      !                    /                               ;                     /
      c(0,0,0) =    109471139332699._R_P/ 163459296000._R_P;c(1,0,0) =   -894628364420801._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(2,0,0) =  34709567828765989._R_P/1307674368000._R_P;c(3,0,0) = -12083632055537503._R_P/ 261534873600._R_P
      !                    /                               ;                     /
      c(4,0,0) =    534237095117903._R_P/  10461394944._R_P;c(5,0,0) = -47841342141981299._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,0) =  21644628077515483._R_P/1307674368000._R_P;c(7,0,0) =  -5644399400246309._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,0,0) =    129739906408601._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,1,0) =                                     0._R_P;c(1,1,0) =   5602753233305651._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(2,1,0) = -59111412950734301._R_P/ 326918592000._R_P;c(3,1,0) = 207178084258860569._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(4,1,0) = -46020384090357023._R_P/ 130767436800._R_P;c(5,1,0) = 165445178916726479._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,0) = -37531036453047161._R_P/ 326918592000._R_P;c(7,1,0) =    137189721025309._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,0) =  -4517524574525093._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,2,0) =                                     0._R_P;c(1,2,0) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,0) =   4660712172178939._R_P/  16982784000._R_P;c(3,2,0) = -45148728224254817._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,0) =    228786920178433._R_P/    212284800._R_P;c(5,2,0) = -36294580012168613._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(6,2,0) =  33008527082236991._R_P/  93405312000._R_P;c(7,2,0) = -30250052825497529._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,0) =  13952443929995611._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,0) =                                     0._R_P;c(1,3,0) =                                     0._R_P
      !                    /                               ;                     /
      c(2,3,0) =                                     0._R_P;c(3,3,0) = 159646773711558347._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,0) =  -7140074733899851._R_P/   3736212480._R_P;c(5,3,0) =  25802513458691833._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,3,0) = -29387187771747941._R_P/  46702656000._R_P;c(7,3,0) = 107887390486248143._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,0) = -24911758529750003._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,0) =                                     0._R_P;c(1,4,0) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,0) =                                     0._R_P;c(3,4,0) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,0) =   8001879703767347._R_P/   7472424960._R_P;c(5,4,0) =  -5794119024433483._R_P/   3736212480._R_P
      !                    /                               ;                     /
      c(6,4,0) =    150205347326833._R_P/    212284800._R_P;c(7,4,0) = -24293471434588703._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,0) =   1123058785015051._R_P/  52306974720._R_P

      !                    /                               ;                     /
      c(0,5,0) =                                     0._R_P;c(1,5,0) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,0) =                                     0._R_P;c(3,5,0) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,0) =                                     0._R_P;c(5,5,0) = 105045730109557451._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,0) = -23993743892557601._R_P/  46702656000._R_P;c(7,5,0) =  88287149743355417._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,5,0) =   -816990037454483._R_P/  52306974720._R_P

      !                    /                               ;                     /
      c(0,6,0) =                                     0._R_P;c(1,6,0) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,0) =                                     0._R_P;c(3,6,0) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,0) =                                     0._R_P;c(5,6,0) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,0) =   1994952741927931._R_P/  16982784000._R_P;c(7,6,0) = -20204125377340061._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,0) =   9355064903078053._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,7,0) =                                     0._R_P;c(1,7,0) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,0) =                                     0._R_P;c(3,7,0) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,0) =                                     0._R_P;c(5,7,0) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,0) =                                     0._R_P;c(7,7,0) =  10637354815456613._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,7,0) =   -189555672759617._R_P/ 100590336000._R_P

      !                    /                               ;                     /
      c(0,8,0) =                                     0._R_P;c(1,8,0) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,0) =                                     0._R_P;c(3,8,0) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,0) =                                     0._R_P;c(5,8,0) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,0) =                                     0._R_P;c(7,8,0) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,0) =     17848737251203._R_P/ 163459296000._R_P

      ! stencil 1
      !                    /                               ;                     /
      c(0,0,1) =     17848737251203._R_P/ 163459296000._R_P;c(1,0,1) =   -147809125548479._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(2,0,1) =   1152669616433567._R_P/ 261534873600._R_P;c(3,0,1) = -10036258935621221._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(4,0,1) =   2214259153735049._R_P/ 261534873600._R_P;c(5,0,1) =  -7906584673048973._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,1) =   3563951929254757._R_P/1307674368000._R_P;c(7,0,1) =     -7406462028919._R_P/  10461394944._R_P
      !                    /                               ;                     /
      c(8,0,1) =    105994418298211._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,1,1) =                                     0._R_P;c(1,1,1) =   6603455065054091._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(2,1,1) = -10036779580858187._R_P/ 326918592000._R_P;c(3,1,1) =  35272568778872279._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(4,1,1) =  -7832368115834609._R_P/ 130767436800._R_P;c(5,1,1) =  28101378954880001._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,1) =  -6356537203415423._R_P/ 326918592000._R_P;c(7,1,1) =     23159841631123._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,1) =   -760053376543163._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,2,1) =                                     0._R_P;c(1,2,1) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,1) =    800572672346869._R_P/  16982784000._R_P;c(3,2,1) =  -7795675329471191._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,1) =     39564077889589._R_P/    212284800._R_P;c(5,2,1) =  -1254519948165511._R_P/   9340531200._R_P
      !                    /                               ;                     /
      c(6,2,1) =   5694325930465457._R_P/  93405312000._R_P;c(7,2,1) =  -5205585064855199._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,1) =   2394338101248133._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,1) =                                     0._R_P;c(1,3,1) =                                     0._R_P
      !                    /                               ;                     /
      c(2,3,1) =                                     0._R_P;c(3,3,1) =  27770723927721989._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,1) =  -6230647138120121._R_P/  18681062400._R_P;c(5,3,1) =  22533757546843859._R_P/  93405312000._R_P
      !                    /                               ;                     /
      c(6,3,1) =  -5129104009946051._R_P/  46702656000._R_P;c(7,3,1) =  18799624487562689._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,1) =  -4331747069079341._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,1) =                                     0._R_P;c(1,4,1) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,1) =                                     0._R_P;c(3,4,1) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,1) =   1403304354475421._R_P/   7472424960._R_P;c(5,4,1) =  -5091060727437401._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,4,1) =     26403598814209._R_P/    212284800._R_P;c(7,4,1) =  -4266972749341649._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,1) =    984850182064169._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,5,1) =                                     0._R_P;c(1,5,1) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,1) =                                     0._R_P;c(3,5,1) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,1) =                                     0._R_P;c(5,5,1) =  18518028023237957._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,1) =  -4234862610936119._R_P/  46702656000._R_P;c(7,5,1) =   3116380997521963._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,5,1) =  -3601784423075141._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,6,1) =                                     0._R_P;c(1,6,1) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,1) =                                     0._R_P;c(3,6,1) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,1) =                                     0._R_P;c(5,6,1) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,1) =    352812369719413._R_P/  16982784000._R_P;c(7,6,1) =  -3575411646556907._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,1) =   1655072196501883._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,7,1) =                                     0._R_P;c(1,7,1) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,1) =                                     0._R_P;c(3,7,1) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,1) =                                     0._R_P;c(5,7,1) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,1) =                                     0._R_P;c(7,7,1) =    269247491159069._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(8,7,1) =    -33593572337951._R_P/ 100590336000._R_P

      !                    /                               ;                     /
      c(0,8,1) =                                     0._R_P;c(1,8,1) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,1) =                                     0._R_P;c(3,8,1) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,1) =                                     0._R_P;c(5,8,1) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,1) =                                     0._R_P;c(7,8,1) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,1) =      3165355170121._R_P/ 163459296000._R_P

      ! stencil 2
      !                    /                               ;                     /
      c(0,0,2) =      3165355170121._R_P/ 163459296000._R_P;c(1,0,2) =     -3844139848343._R_P/  14370048000._R_P
      !                    /                               ;                     /
      c(2,0,2) =   1063191201446533._R_P/1307674368000._R_P;c(3,0,2) =  -1859899247394491._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(4,0,2) =    409921790776919._R_P/ 261534873600._R_P;c(5,0,2) =  -1457105112643091._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,2) =    652452925567483._R_P/1307674368000._R_P;c(7,0,2) =   -168172381487813._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,0,2) =     19094704104061._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,1,2) =                                     0._R_P;c(1,1,2) =   1239990283564133._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(2,1,2) =  -1918610096603357._R_P/ 326918592000._R_P;c(3,1,2) =   1359891017166853._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(4,1,2) =  -1512744281500799._R_P/ 130767436800._R_P;c(5,1,2) =   5414972538444239._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,2) =  -1218782466526649._R_P/ 326918592000._R_P;c(7,1,2) =      4411553510173._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,2) =     -5748413034701._R_P/  52306974720._R_P

      !                    /                               ;                     /
      c(0,2,2) =                                     0._R_P;c(1,2,2) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,2) =    156622544328763._R_P/  16982784000._R_P;c(3,2,2) =  -1544964557143169._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,2) =      7883820528109._R_P/    212284800._R_P;c(5,2,2) =  -1250454991752101._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(6,2,2) =   1131898542897407._R_P/  93405312000._R_P;c(7,2,2) =  -1029608247917273._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,2) =    470643665358907._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,2) =                                     0._R_P;c(1,3,2) =                                     0._R_P
      !                    /                               ;                     /
      c(2,3,2) =                                     0._R_P;c(3,3,2) =   5599666272693707._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,2) =  -1267992294203351._R_P/  18681062400._R_P;c(5,3,2) =   4601782036044509._R_P/  93405312000._R_P
      !                    /                               ;                     /
      c(6,3,2) =   -209388842757121._R_P/   9340531200._R_P;c(7,3,2) =   3825435713279951._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,2) =   -877252492928723._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,2) =                                     0._R_P;c(1,4,2) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,2) =                                     0._R_P;c(3,4,2) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,2) =    289259235638771._R_P/   7472424960._R_P;c(5,4,2) =  -1056291616534871._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,4,2) =      5489435141989._R_P/    212284800._R_P;c(7,4,2) =   -886173785909759._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,2) =    203891614104599._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,5,2) =                                     0._R_P;c(1,5,2) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,2) =                                     0._R_P;c(3,5,2) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,2) =                                     0._R_P;c(5,5,2) =   3878296682785739._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,2) =   -890937252684641._R_P/  46702656000._R_P;c(7,5,2) =   3281427995720729._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,5,2) =   -757402017640571._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,6,2) =                                     0._R_P;c(1,6,2) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,2) =                                     0._R_P;c(3,6,2) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,2) =                                     0._R_P;c(5,6,2) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,2) =     74730821653819._R_P/  16982784000._R_P;c(7,6,2) =   -759598480120637._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,2) =     70341062456897._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,7,2) =                                     0._R_P;c(1,7,2) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,2) =                                     0._R_P;c(3,7,2) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,2) =                                     0._R_P;c(5,7,2) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,2) =                                     0._R_P;c(7,7,2) =    402355798141541._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,7,2) =     -1026441378647._R_P/  14370048000._R_P

      !                    /                               ;                     /
      c(0,8,2) =                                     0._R_P;c(1,8,2) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,2) =                                     0._R_P;c(3,8,2) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,2) =                                     0._R_P;c(5,8,2) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,2) =                                     0._R_P;c(7,8,2) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,2) =       679328101453._R_P/ 163459296000._R_P

      ! stencil 3
      !                    /                               ;                     /
      c(0,0,3) =       679328101453._R_P/ 163459296000._R_P;c(1,0,3) =     -6056041731167._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(2,0,3) =    247582660569403._R_P/1307674368000._R_P;c(3,0,3) =    -17694932119757._R_P/  52306974720._R_P
      !                    /                               ;                     /
      c(4,0,3) =     19690918384021._R_P/  52306974720._R_P;c(5,0,3) =   -350067382006253._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,3) =    155614950712261._R_P/1307674368000._R_P;c(7,0,3) =    -39587674152443._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,0,3) =       883416230471._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,1,3) =                                     0._R_P;c(1,1,3) =    293675114165963._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(2,1,3) =   -472662830894411._R_P/ 326918592000._R_P;c(3,1,3) =   1720297891825367._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(4,1,3) =   -388442316668753._R_P/ 130767436800._R_P;c(5,1,3) =   1397141337414593._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,3) =   -313421131078079._R_P/ 326918592000._R_P;c(7,1,3) =      1123540717459._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,3) =    -36073774922459._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,2,3) =                                     0._R_P;c(1,2,3) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,3) =     40385614392181._R_P/  16982784000._R_P;c(3,2,3) =   -411721854332951._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,3) =      2145005788633._R_P/    212284800._R_P;c(5,2,3) =   -343655982425891._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(6,2,3) =    311458280689841._R_P/  93405312000._R_P;c(7,2,3) =   -281678601090911._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,3) =    127326292586533._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,3) =                                     0._R_P;c(1,3,3) =                                     0._R_P
      !                    /                               ;                     /
      c(3,3,3) =                                     0._R_P;c(3,3,3) =   1553225813426501._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,3) =    -72310955346373._R_P/   3736212480._R_P;c(5,3,3) =    266698467235063._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,3,3) =   -305368847812163._R_P/  46702656000._R_P;c(7,3,3) =   1114386138224129._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,3) =   -253674820236749._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,3) =                                     0._R_P;c(1,4,3) =                                     0._R_P
      !                    /                               ;                     /
      c(3,4,3) =                                     0._R_P;c(3,4,3) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,3) =     85394018909597._R_P/   7472424960._R_P;c(5,4,3) =    -63811818908581._R_P/   3736212480._R_P
      !                    /                               ;                     /
      c(6,4,3) =      1679094624733._R_P/    212284800._R_P;c(7,4,3) =   -272139518377073._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,3) =      2497209723185._R_P/  10461394944._R_P

      !                    /                               ;                     /
      c(0,5,3) =                                     0._R_P;c(1,5,3) =                                     0._R_P
      !                    /                               ;                     /
      c(3,5,3) =                                     0._R_P;c(3,5,3) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,3) =                                     0._R_P;c(5,5,3) =   1206964694318597._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,3) =   -282622107973367._R_P/  46702656000._R_P;c(7,5,3) =   1051238439516119._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,5,3) =    -48633489917473._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,6,3) =                                     0._R_P;c(1,6,3) =                                     0._R_P
      !                    /                               ;                     /
      c(3,6,3) =                                     0._R_P;c(3,6,3) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,3) =                                     0._R_P;c(5,6,3) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,3) =     24324934655989._R_P/  16982784000._R_P;c(7,6,3) =   -251283767228651._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,3) =    117272649474139._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,7,3) =                                     0._R_P;c(1,7,3) =                                     0._R_P
      !                    /                               ;                     /
      c(3,7,3) =                                     0._R_P;c(3,7,3) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,3) =                                     0._R_P;c(5,7,3) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,3) =                                     0._R_P;c(7,7,3) =    136155780967307._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,7,3) =     -2466233185151._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(0,8,3) =                                     0._R_P;c(1,8,3) =                                     0._R_P
      !                    /                               ;                     /
      c(3,8,3) =                                     0._R_P;c(3,8,3) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,3) =                                     0._R_P;c(5,8,3) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,3) =                                     0._R_P;c(7,8,3) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,3) =       238114846399._R_P/ 163459296000._R_P

      ! stencil 4
      !                    /                               ;                     /
      c(0,0,4) =       238114846399._R_P/ 163459296000._R_P;c(1,0,4) =     -2297804363777._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(2,0,4) =     20216075320673._R_P/ 261534873600._R_P;c(3,0,4) =   -192700060973723._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(4,0,4) =     45272942020727._R_P/ 261534873600._R_P;c(5,0,4) =   -167888314942259._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,4) =     76858903972891._R_P/1307674368000._R_P;c(7,0,4) =     -3976300410337._R_P/ 261534873600._R_P
      !                    /                               ;                     /
      c(8,0,4) =      2227506474493._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,0,4) =                                     0._R_P;c(1,1,4) =    119979314906981._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(2,1,4) =   -207359252612669._R_P/ 326918592000._R_P;c(3,1,4) =    161084839253509._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(4,1,4) =   -192310346872991._R_P/ 130767436800._R_P;c(5,1,4) =    723357784442063._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,4) =   -167690675241113._R_P/ 326918592000._R_P;c(7,1,4) =       613753663261._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,4) =     -3976300410337._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,2,4) =                                     0._R_P;c(1,2,4) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,4) =     19010310966523._R_P/  16982784000._R_P;c(3,2,4) =   -207059158040897._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,4) =      1143576251161._R_P/    212284800._R_P;c(5,2,4) =    -38450763316993._R_P/   9340531200._R_P
      !                    /                               ;                     /
      c(6,2,4) =    180786151740479._R_P/  93405312000._R_P;c(7,2,4) =   -167690675241113._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,4) =     76858903972891._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,4) =                                     0._R_P;c(1,3,4) =                                     0._R_P
      !                    /                               ;                     /
      c(2,3,4) =                                     0._R_P;c(3,3,4) =    836484368637131._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,4) =   -207139067201783._R_P/  18681062400._R_P;c(5,3,4) =    805195803373277._R_P/  93405312000._R_P
      !                    /                               ;                     /
      c(6,3,4) =    -38450763316993._R_P/   9340531200._R_P;c(7,3,4) =    723357784442063._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,4) =   -167888314942259._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,4) =                                     0._R_P;c(1,4,4) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,4) =                                     0._R_P;c(3,4,4) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,4) =     52297392889139._R_P/   7472424960._R_P;c(5,4,4) =   -207139067201783._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,4,4) =      1143576251161._R_P/    212284800._R_P;c(7,4,4) =   -192310346872991._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,4) =     45272942020727._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,5,4) =                                     0._R_P;c(1,5,4) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,4) =                                     0._R_P;c(3,5,4) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,4) =                                     0._R_P;c(5,5,4) =    836484368637131._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,4) =   -207059158040897._R_P/  46702656000._R_P;c(7,5,4) =    161084839253509._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,5,4) =   -192700060973723._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,6,4) =                                     0._R_P;c(1,6,4) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,4) =                                     0._R_P;c(3,6,4) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,4) =                                     0._R_P;c(5,6,4) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,4) =     19010310966523._R_P/  16982784000._R_P;c(7,6,4) =   -207359252612669._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,4) =     20216075320673._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,7,4) =                                     0._R_P;c(1,7,4) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,4) =                                     0._R_P;c(3,7,4) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,4) =                                     0._R_P;c(5,7,4) =                                     0._R_P
      !                    /                               ;                     /
      c(0,7,4) =                                     0._R_P;c(7,7,4) =    119979314906981._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,7,4) =     -2297804363777._R_P/ 100590336000._R_P

      !                    /                               ;                     /
      c(0,8,4) =                                     0._R_P;c(1,8,4) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,4) =                                     0._R_P;c(3,8,4) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,4) =                                     0._R_P;c(5,8,4) =                                     0._R_P
      !                    /                               ;                     /
      c(0,8,4) =                                     0._R_P;c(7,8,4) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,4) =       238114846399._R_P/ 163459296000._R_P

      ! stencil 5
      !                    /                               ;                     /
      c(0,0,5) =       238114846399._R_P/ 163459296000._R_P;c(1,0,5) =     -2466233185151._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(2,0,5) =    117272649474139._R_P/1307674368000._R_P;c(3,0,5) =    -48633489917473._R_P/ 261534873600._R_P
      !                    /                               ;                     /
      c(4,0,5) =      2497209723185._R_P/  10461394944._R_P;c(5,0,5) =   -253674820236749._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,5) =    127326292586533._R_P/1307674368000._R_P;c(7,0,5) =    -36073774922459._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,0,5) =       883416230471._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,1,5) =                                     0._R_P;c(1,1,5) =    136155780967307._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(2,1,5) =   -251283767228651._R_P/ 326918592000._R_P;c(3,1,5) =   1051238439516119._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(4,1,5) =   -272139518377073._R_P/ 130767436800._R_P;c(5,1,5) =   1114386138224129._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,5) =   -281678601090911._R_P/ 326918592000._R_P;c(7,1,5) =      1123540717459._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,5) =    -39587674152443._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,2,5) =                                     0._R_P;c(1,2,5) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,5) =     24324934655989._R_P/  16982784000._R_P;c(3,2,5) =   -282622107973367._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,5) =      1679094624733._R_P/    212284800._R_P;c(5,2,5) =   -305368847812163._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(6,2,5) =    311458280689841._R_P/  93405312000._R_P;c(7,2,5) =   -313421131078079._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,5) =    155614950712261._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,5) =                                     0._R_P;c(1,3,5) =                                     0._R_P
      !                    /                               ;                     /
      c(2,3,5) =                                     0._R_P;c(3,3,5) =   1206964694318597._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,5) =    -63811818908581._R_P/   3736212480._R_P;c(5,3,5) =    266698467235063._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,3,5) =   -343655982425891._R_P/  46702656000._R_P;c(7,3,5) =   1397141337414593._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,5) =   -350067382006253._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,5) =                                     0._R_P;c(1,4,5) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,5) =                                     0._R_P;c(3,4,5) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,5) =     85394018909597._R_P/   7472424960._R_P;c(5,4,5) =    -72310955346373._R_P/   3736212480._R_P
      !                    /                               ;                     /
      c(6,4,5) =      2145005788633._R_P/    212284800._R_P;c(7,4,5) =   -388442316668753._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,5) =     19690918384021._R_P/  52306974720._R_P

      !                    /                               ;                     /
      c(0,5,5) =                                     0._R_P;c(1,5,5) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,5) =                                     0._R_P;c(3,5,5) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,5) =                                     0._R_P;c(5,5,5) =   1553225813426501._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,5) =   -411721854332951._R_P/  46702656000._R_P;c(7,5,5) =   1720297891825367._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,5,5) =    -17694932119757._R_P/  52306974720._R_P

      !                    /                               ;                     /
      c(0,6,5) =                                     0._R_P;c(1,6,5) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,5) =                                     0._R_P;c(3,6,5) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,5) =                                     0._R_P;c(5,6,5) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,5) =     40385614392181._R_P/  16982784000._R_P;c(7,6,5) =   -472662830894411._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,5) =    247582660569403._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,7,5) =                                     0._R_P;c(1,7,5) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,5) =                                     0._R_P;c(3,7,5) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,5) =                                     0._R_P;c(5,7,5) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,5) =                                     0._R_P;c(7,7,5) =    293675114165963._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,7,5) =     -6056041731167._R_P/ 100590336000._R_P

      !                    /                               ;                     /
      c(0,8,5) =                                     0._R_P;c(1,8,5) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,5) =                                     0._R_P;c(3,8,5) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,5) =                                     0._R_P;c(5,8,5) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,5) =                                     0._R_P;c(7,8,5) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,5) =       679328101453._R_P/ 163459296000._R_P

      ! stencil 6
      !                    /                               ;                     /
      c(0,0,6) =       679328101453._R_P/ 163459296000._R_P;c(1,0,6) =     -1026441378647._R_P/  14370048000._R_P
      !                    /                               ;                     /
      c(2,0,6) =     70341062456897._R_P/ 261534873600._R_P;c(3,0,6) =   -757402017640571._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(4,0,6) =    203891614104599._R_P/ 261534873600._R_P;c(5,0,6) =   -877252492928723._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,6) =    470643665358907._R_P/1307674368000._R_P;c(7,0,6) =     -5748413034701._R_P/  52306974720._R_P
      !                    /                               ;                     /
      c(8,0,6) =     19094704104061._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,1,6) =                                     0._R_P;c(1,1,6) =    402355798141541._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(2,1,6) =   -759598480120637._R_P/ 326918592000._R_P;c(3,1,6) =   3281427995720729._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(4,1,6) =   -886173785909759._R_P/ 130767436800._R_P;c(5,1,6) =   3825435713279951._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,6) =  -1029608247917273._R_P/ 326918592000._R_P;c(7,1,6) =      4411553510173._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,6) =   -168172381487813._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,2,6) =                                     0._R_P;c(1,2,6) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,6) =     74730821653819._R_P/  16982784000._R_P;c(3,2,6) =   -890937252684641._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,6) =      5489435141989._R_P/    212284800._R_P;c(5,2,6) =   -209388842757121._R_P/   9340531200._R_P
      !                    /                               ;                     /
      c(6,2,6) =   1131898542897407._R_P/  93405312000._R_P;c(7,2,6) =  -1218782466526649._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,6) =    652452925567483._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,6) =                                     0._R_P;c(1,3,6) =                                     0._R_P
      !                    /                               ;                     /
      c(2,3,6) =                                     0._R_P;c(3,3,6) =   3878296682785739._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,6) =  -1056291616534871._R_P/  18681062400._R_P;c(5,3,6) =   4601782036044509._R_P/  93405312000._R_P
      !                    /                               ;                     /
      c(6,3,6) =  -1250454991752101._R_P/  46702656000._R_P;c(7,3,6) =   5414972538444239._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,6) =  -1457105112643091._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,6) =                                     0._R_P;c(1,4,6) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,6) =                                     0._R_P;c(3,4,6) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,6) =    289259235638771._R_P/   7472424960._R_P;c(5,4,6) =  -1267992294203351._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,4,6) =      7883820528109._R_P/    212284800._R_P;c(7,4,6) =  -1512744281500799._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,6) =    409921790776919._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,5,6) =                                     0._R_P;c(1,5,6) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,6) =                                     0._R_P;c(3,5,6) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,6) =                                     0._R_P;c(5,5,6) =   5599666272693707._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,6) =  -1544964557143169._R_P/  46702656000._R_P;c(7,5,6) =   1359891017166853._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,5,6) =  -1859899247394491._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,6,6) =                                     0._R_P;c(1,6,6) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,6) =                                     0._R_P;c(3,6,6) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,6) =                                     0._R_P;c(5,6,6) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,6) =    156622544328763._R_P/  16982784000._R_P;c(7,6,6) =  -1918610096603357._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,6) =   1063191201446533._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,7,6) =                                     0._R_P;c(1,7,6) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,6) =                                     0._R_P;c(3,7,6) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,6) =                                     0._R_P;c(5,7,6) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,6) =                                     0._R_P;c(7,7,6) =   1239990283564133._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,7,6) =     -3844139848343._R_P/  14370048000._R_P

      !                    /                               ;                     /
      c(0,8,6) =                                     0._R_P;c(1,8,6) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,6) =                                     0._R_P;c(3,8,6) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,6) =                                     0._R_P;c(5,8,6) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,6) =                                     0._R_P;c(7,8,6) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,6) =      3165355170121._R_P/ 163459296000._R_P

      ! stencil 7
      !                    /                               ;                     /
      c(0,0,7) =      3165355170121._R_P/ 163459296000._R_P;c(1,0,7) =    -33593572337951._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(2,0,7) =   1655072196501883._R_P/1307674368000._R_P;c(3,0,7) =  -3601784423075141._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(4,0,7) =    984850182064169._R_P/ 261534873600._R_P;c(5,0,7) =  -4331747069079341._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,7) =   2394338101248133._R_P/1307674368000._R_P;c(7,0,7) =   -760053376543163._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,0,7) =    105994418298211._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,1,7) =                                     0._R_P;c(1,1,7) =    269247491159069._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(2,1,7) =  -3575411646556907._R_P/ 326918592000._R_P;c(3,1,7) =   3116380997521963._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(4,1,7) =  -4266972749341649._R_P/ 130767436800._R_P;c(5,1,7) =  18799624487562689._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,7) =  -5205585064855199._R_P/ 326918592000._R_P;c(7,1,7) =     23159841631123._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,7) =     -7406462028919._R_P/  10461394944._R_P

      !                    /                               ;                     /
      c(0,2,7) =                                     0._R_P;c(1,2,7) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,7) =    352812369719413._R_P/  16982784000._R_P;c(3,2,7) =  -4234862610936119._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,7) =     26403598814209._R_P/    212284800._R_P;c(5,2,7) =  -5129104009946051._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(6,2,7) =   5694325930465457._R_P/  93405312000._R_P;c(7,2,7) =  -6356537203415423._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,7) =   3563951929254757._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,2,7) =                                     0._R_P;c(1,3,7) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,7) =                                     0._R_P;c(3,3,7) =  18518028023237957._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,7) =  -5091060727437401._R_P/  18681062400._R_P;c(5,3,7) =  22533757546843859._R_P/  93405312000._R_P
      !                    /                               ;                     /
      c(6,3,7) =  -1254519948165511._R_P/   9340531200._R_P;c(7,3,7) =  28101378954880001._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,7) =  -7906584673048973._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,7) =                                     0._R_P;c(1,4,7) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,7) =                                     0._R_P;c(3,4,7) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,7) =   1403304354475421._R_P/   7472424960._R_P;c(5,4,7) =  -6230647138120121._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,4,7) =     39564077889589._R_P/    212284800._R_P;c(7,4,7) =  -7832368115834609._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,7) =   2214259153735049._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,5,7) =                                     0._R_P;c(1,5,7) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,7) =                                     0._R_P;c(3,5,7) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,7) =                                     0._R_P;c(5,5,7) =  27770723927721989._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,7) =  -7795675329471191._R_P/  46702656000._R_P;c(7,5,7) =  35272568778872279._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,5,7) = -10036258935621221._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,6,7) =                                     0._R_P;c(1,6,7) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,7) =                                     0._R_P;c(3,6,7) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,7) =                                     0._R_P;c(5,6,7) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,7) =    800572672346869._R_P/  16982784000._R_P;c(7,6,7) = -10036779580858187._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,7) =   1152669616433567._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,7,7) =                                     0._R_P;c(1,7,7) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,7) =                                     0._R_P;c(3,7,7) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,7) =                                     0._R_P;c(5,7,7) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,7) =                                     0._R_P;c(7,7,7) =   6603455065054091._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,7,7) =   -147809125548479._R_P/ 100590336000._R_P

      !                    /                               ;                     /
      c(0,8,7) =                                     0._R_P;c(1,8,7) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,7) =                                     0._R_P;c(3,8,7) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,7) =                                     0._R_P;c(5,8,7) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,7) =                                     0._R_P;c(7,8,7) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,7) =     17848737251203._R_P/ 163459296000._R_P

      ! stencil 8
      !                    /                               ;                     /
      c(0,0,8) =     17848737251203._R_P/ 163459296000._R_P;c(1,0,8) =   -189555672759617._R_P/ 100590336000._R_P
      !                    /                               ;                     /
      c(2,0,8) =   9355064903078053._R_P/1307674368000._R_P;c(3,0,8) =   -816990037454483._R_P/  52306974720._R_P
      !                    /                               ;                     /
      c(4,0,8) =   1123058785015051._R_P/  52306974720._R_P;c(5,0,8) = -24911758529750003._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(6,0,8) =  13952443929995611._R_P/1307674368000._R_P;c(7,0,8) =  -4517524574525093._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(8,0,8) =    129739906408601._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,1,8) =                                     0._R_P;c(1,1,8) =  10637354815456613._R_P/1307674368000._R_P
      !                    /                               ;                     /
      c(2,1,8) = -20204125377340061._R_P/ 326918592000._R_P;c(3,1,8) =  88287149743355417._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(4,1,8) = -24293471434588703._R_P/ 130767436800._R_P;c(5,1,8) = 107887390486248143._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(6,1,8) = -30250052825497529._R_P/ 326918592000._R_P;c(7,1,8) =    137189721025309._R_P/   4572288000._R_P
      !                    /                               ;                     /
      c(8,1,8) =  -5644399400246309._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,2,8) =                                     0._R_P;c(1,2,8) =                                     0._R_P
      !                    /                               ;                     /
      c(2,2,8) =   1994952741927931._R_P/  16982784000._R_P;c(3,2,8) = -23993743892557601._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(4,2,8) =    150205347326833._R_P/    212284800._R_P;c(5,2,8) = -29387187771747941._R_P/  46702656000._R_P
      !                    /                               ;                     /
      c(6,2,8) =  33008527082236991._R_P/  93405312000._R_P;c(7,2,8) = -37531036453047161._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,2,8) =  21644628077515483._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,3,8) =                                     0._R_P;c(1,3,8) =                                     0._R_P
      !                    /                               ;                     /
      c(2,3,8) =                                     0._R_P;c(3,3,8) = 105045730109557451._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(4,3,8) =  -5794119024433483._R_P/   3736212480._R_P;c(5,3,8) =  25802513458691833._R_P/  18681062400._R_P
      !                    /                               ;                     /
      c(6,3,8) = -36294580012168613._R_P/  46702656000._R_P;c(7,3,8) = 165445178916726479._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,3,8) = -47841342141981299._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,4,8) =                                     0._R_P;c(1,4,8) =                                     0._R_P
      !                    /                               ;                     /
      c(2,4,8) =                                     0._R_P;c(3,4,8) =                                     0._R_P
      !                    /                               ;                     /
      c(4,4,8) =   8001879703767347._R_P/   7472424960._R_P;c(5,4,8) =  -7140074733899851._R_P/   3736212480._R_P
      !                    /                               ;                     /
      c(6,4,8) =    228786920178433._R_P/    212284800._R_P;c(7,4,8) = -46020384090357023._R_P/ 130767436800._R_P
      !                    /                               ;                     /
      c(8,4,8) =    534237095117903._R_P/  10461394944._R_P

      !                    /                               ;                     /
      c(0,5,8) =                                     0._R_P;c(1,5,8) =                                     0._R_P
      !                    /                               ;                     /
      c(2,5,8) =                                     0._R_P;c(3,5,8) =                                     0._R_P
      !                    /                               ;                     /
      c(4,5,8) =                                     0._R_P;c(5,5,8) = 159646773711558347._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(6,5,8) = -45148728224254817._R_P/  46702656000._R_P;c(7,5,8) = 207178084258860569._R_P/ 653837184000._R_P
      !                    /                               ;                     /
      c(8,5,8) = -12083632055537503._R_P/ 261534873600._R_P

      !                    /                               ;                     /
      c(0,6,8) =                                     0._R_P;c(1,6,8) =                                     0._R_P
      !                    /                               ;                     /
      c(2,6,8) =                                     0._R_P;c(3,6,8) =                                     0._R_P
      !                    /                               ;                     /
      c(4,6,8) =                                     0._R_P;c(5,6,8) =                                     0._R_P
      !                    /                               ;                     /
      c(6,6,8) =   4660712172178939._R_P/  16982784000._R_P;c(7,6,8) = -59111412950734301._R_P/ 326918592000._R_P
      !                    /                               ;                     /
      c(8,6,8) =  34709567828765989._R_P/1307674368000._R_P

      !                    /                               ;                     /
      c(0,7,8) =                                     0._R_P;c(1,7,8) =                                     0._R_P
      !                    /                               ;                     /
      c(2,7,8) =                                     0._R_P;c(3,7,8) =                                     0._R_P
      !                    /                               ;                     /
      c(4,7,8) =                                     0._R_P;c(5,7,8) =                                     0._R_P
      !                    /                               ;                     /
      c(6,7,8) =                                     0._R_P;c(7,7,8) =   5602753233305651._R_P/ 186810624000._R_P
      !                    /                               ;                     /
      c(8,7,8) =   -894628364420801._R_P/ 100590336000._R_P

      !                    /                               ;                     /
      c(0,8,8) =                                     0._R_P;c(1,8,8) =                                     0._R_P
      !                    /                               ;                     /
      c(2,8,8) =                                     0._R_P;c(3,8,8) =                                     0._R_P
      !                    /                               ;                     /
      c(4,8,8) =                                     0._R_P;c(5,8,8) =                                     0._R_P
      !                    /                               ;                     /
      c(6,8,8) =                                     0._R_P;c(7,8,8) =                                     0._R_P
      !                    /                               ;                     /
      c(8,8,8) =    109471139332699._R_P/ 163459296000._R_P
    endselect
  endassociate
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy smoothness indicators.
  class(smoothness_indicators_js), intent(inout) :: self !< Smoothenss indicators.

  call self%smoothness_indicators%destroy
  if (allocated(self%coef)) deallocate(self%coef)
  endsubroutine destroy
endmodule wenoof_smoothness_indicators_js
