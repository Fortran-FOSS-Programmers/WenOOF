!< Jiang-Shu (Lagrange) interpolations object for function interpolation.
module wenoof_interpolations_int_js
!< Jiang-Shu (Lagrange) interpolations object for function interpolation.
!<
!< @note The provided interpolations implement the Lagrange interpolations defined in *High Order Weighted Essentially
!< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
!< doi:10.1137/070679065.

use penf, only : I_P, R_P, str
use wenoof_base_object, only : base_object, base_object_constructor
use wenoof_interpolations_object, only : interpolations_object, interpolations_object_constructor

implicit none
private
public :: interpolations_int_js
public :: interpolations_int_js_constructor

type, extends(interpolations_object_constructor) :: interpolations_int_js_constructor
  !< Jiang-Shu (Lagrange) interpolations object for function interpolation constructor.
  real(R_P)              :: x_target   !< Coordinate of the interpolation point.
  contains
    ! public deferred methods
    procedure, pass(lhs) :: constr_assign_constr !< `=` operator.
endtype interpolations_int_js_constructor

type, extends(interpolations_object) :: interpolations_int_js
  !< Jiang-Shu (Lagrange) interpolations object for function interpolation.
  !<
  !< @note The provided interpolations implement the Lagrange interpolations defined in *High Order Weighted Essentially
  !< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
  !< doi:10.1137/070679065.
  real(R_P), allocatable :: coef(:,:,:) !< Polynomial coefficients [0:S-1,0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: create               !< Create interpolations.
    procedure, pass(self) :: compute_int          !< Compute interpolations (interpolate).
    procedure, pass(self) :: compute_rec          !< Compute interpolations (reconstruct).
    procedure, pass(self) :: description          !< Return object string-description.
    procedure, pass(self) :: destroy              !< Destroy interpolations.
    procedure, pass(lhs)  :: object_assign_object !< `=` operator.
endtype interpolations_int_js

contains
  ! constructor

  ! deferred public methods
  subroutine constr_assign_constr(lhs, rhs)
  !< `=` operator.
  class(interpolations_int_js_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor),           intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(interpolations_int_js_constructor)
     lhs%x_target = rhs%x_target
  endselect
  endsubroutine constr_assign_constr

  ! public deferred methods
  subroutine create(self, constructor)
  !< Create interpolations.
  class(interpolations_int_js),   intent(inout) :: self        !< Interpolations.
  class(base_object_constructor), intent(in)    :: constructor !< Interpolations constructor.
  integer(I_P)                                  :: i           !< Counter.

  call self%destroy
  call self%create_(constructor=constructor)
  if (self%ror) then
    allocate(self%coef(0:self%S - 1, 0:self%S - 1, 2:self%S))
    select type(constructor)
    type is(interpolations_int_js_constructor)
      do i=2, self%S
        call assign_interp_coeff(c=self%coef(:,:,i), x_target=constructor%x_target, S=i)
      enddo
    endselect
  else
    allocate(self%coef(0:self%S - 1, 0:self%S - 1, 1))
    select type(constructor)
    type is(interpolations_int_js_constructor)
      call assign_interp_coeff(c=self%coef(:,:,1), x_target=constructor%x_target, S=self%S)
    endselect
  endif
  endsubroutine create

  pure subroutine compute_int(self, ord, stencil, values)
  !< Compute interpolations (interpolation).
  class(interpolations_int_js), intent(in)  :: self            !< Interpolations.
  integer(I_P),                 intent(in)  :: ord             !< Interpolation order.
  real(R_P),                    intent(in)  :: stencil(1-ord:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(R_P),                    intent(out) :: values(0:)      !< Interpolations values.
  integer(I_P)                              :: s1              !< Counter.
  integer(I_P)                              :: s2              !< Counter.

  values = 0._R_P
  !do s1=0, self%S - 1 ! stencils loop
  !  do s2=0, self%S - 1 ! values loop
  do s1=0, ord-1  ! stencils loop
    do s2=0, ord-1  ! values loop
      values(s1) = values(s1) + self%coef(s2, s1, ord) * stencil(-s2 + s1)
    enddo
  enddo
  endsubroutine compute_int

  pure subroutine compute_rec(self, ord, stencil, values)
  !< Compute interpolations (reconstruct).
  class(interpolations_int_js), intent(in)  :: self               !< Interpolations.
  integer(I_P),                 intent(in)  :: ord                !< Interpolation order.
  real(R_P),                    intent(in)  :: stencil(1:,1-ord:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(R_P),                    intent(out) :: values(1:, 0:)     !< Interpolations values.
  ! empty procedure
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(interpolations_int_js), intent(in)           :: self             !< Interpolations.
  character(len=*),             intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                      :: string           !< String-description.
  character(len=:), allocatable                      :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                        :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu beta interpolations object for interpolation:'//NL
  string = string//prefix_//'  - S   = '//trim(str(self%S))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy interpolations.
  class(interpolations_int_js), intent(inout) :: self !< Interpolations.

  call self%destroy_
  if (allocated(self%coef)) deallocate(self%coef)
  endsubroutine destroy

  pure subroutine object_assign_object(lhs, rhs)
  !< `=` operator.
  class(interpolations_int_js), intent(inout) :: lhs !< Left hand side.
  class(base_object),           intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(interpolations_int_js)
     if (allocated(rhs%coef)) then
        lhs%coef = rhs%coef
     else
        if (allocated(lhs%coef)) deallocate(lhs%coef)
     endif
  endselect
  endsubroutine object_assign_object

  ! private non TBP
  pure subroutine assign_interp_coeff(c, x_target, S)
  !< Assign the value of beta coefficients.
  real(R_P),    intent(inout) :: c(0:,0:)    !< Coefficient values.
  real(R_P),    intent(in)    :: x_target    !< Interpolation coordinate.
  integer(I_P), intent(in)    :: S           !< Counter.
  real(R_P)                   :: prod, c_sum !< Temporary variables.
  integer(I_P)                :: i, j, k     !< Counters.

  if(x_target==-0.5_R_P) then
    ! left interface (i-1/2)
    select case(S)
      case(2) ! 3rd order
        !  cell  1      ;    cell  0
        c(1,0)=  0.5_R_P; c(0,0)=  0.5_R_P ! stencil 0
        c(1,1)=  1.5_R_P; c(0,1)= -0.5_R_P ! stencil 1
      case(3) ! 5th order
        !  cell  2            ;    cell  1            ;    cell  0
        c(2,0)= -1._R_P/8._R_P; c(1,0)=  3._R_P/4._R_P; c(0,0)=  3._R_P/8._R_P ! stencil 0
        c(2,1)=  3._R_P/8._R_P; c(1,1)=  3._R_P/4._R_P; c(0,1)= -1._R_P/8._R_P ! stencil 1
        c(2,2)= 15._R_P/8._R_P; c(1,2)= -5._R_P/4._R_P; c(0,2)=  3._R_P/8._R_P ! stencil 2
      case(4) ! 7th order
        !  cell  3             ;    cell  2             ;    cell  1             ;    cell  0
        c(3,0)=  1._R_P/16._R_P; c(2,0)= -5._R_P/16._R_P; c(1,0)= 15._R_P/16._R_P; c(0,0)=  5._R_P/16._R_P ! stencil 0
        c(3,1)= -1._R_P/16._R_P; c(2,1)=  9._R_P/16._R_P; c(1,1)=  9._R_P/16._R_P; c(0,1)= -1._R_P/16._R_P ! stencil 1
        c(3,2)=  5._R_P/16._R_P; c(2,2)= 15._R_P/16._R_P; c(1,2)= -5._R_P/16._R_P; c(0,2)=  1._R_P/16._R_P ! stencil 2
        c(3,3)= 35._R_P/16._R_P; c(2,3)=-35._R_P/16._R_P; c(1,3)= 21._R_P/16._R_P; c(0,3)= -5._R_P/16._R_P ! stencil 3
      case(5) ! 9th order
        !  cell  4               ;    cell  3               ;    cell  2               ;    cell  1
        c(4,0)=  -5._R_P/128._R_P; c(3,0)=   7._R_P/32._R_P ; c(2,0)= -35._R_P/64._R_P ; c(1,0)=  35._R_P/32._R_P  ! stencil 0
        c(4,1)=   3._R_P/128._R_P; c(3,1)=  -5._R_P/32._R_P ; c(2,1)=  45._R_P/64._R_P ; c(1,1)=  15._R_P/32._R_P  ! stencil 1
        c(4,2)=  -5._R_P/128._R_P; c(3,2)=  15._R_P/32._R_P ; c(2,2)=  45._R_P/64._R_P ; c(1,2)=  -5._R_P/32._R_P  ! stencil 2
        c(4,3)=  35._R_P/128._R_P; c(3,3)=  35._R_P/32._R_P ; c(2,3)= -35._R_P/64._R_P ; c(1,3)=   7._R_P/32._R_P  ! stencil 3
        c(4,4)= 315._R_P/128._R_P; c(3,4)=-105._R_P/32._R_P ; c(2,4)= 189._R_P/64._R_P ; c(1,4)= -45._R_P/32._R_P  ! stencil 4
        !  cell  0
        c(0,0)=  35._R_P/128._R_P  ! stencil 0
        c(0,1)=  -5._R_P/128._R_P  ! stencil 1
        c(0,2)=   3._R_P/128._R_P  ! stencil 2
        c(0,3)=  -5._R_P/128._R_P  ! stencil 3
        c(0,4)=  35._R_P/128._R_P  ! stencil 4
      case(6) ! 11th order
        !  cell  5                ;    cell  4                ;    cell  3
        c(5,0)=    7._R_P/256._R_P; c(4,0)=  -45._R_P/256._R_P; c(3,0)=   63._R_P/128._R_P  ! stencil 0
        c(5,1)=   -3._R_P/256._R_P; c(4,1)=   21._R_P/256._R_P; c(3,1)=  -35._R_P/128._R_P  ! stencil 1
        c(5,2)=    3._R_P/256._R_P; c(4,2)=  -25._R_P/256._R_P; c(3,2)=   75._R_P/128._R_P  ! stencil 2
        c(5,3)=   -7._R_P/256._R_P; c(4,3)=  105._R_P/256._R_P; c(3,3)=  105._R_P/128._R_P  ! stencil 3
        c(5,4)=   63._R_P/256._R_P; c(4,4)=  315._R_P/256._R_P; c(3,4)= -105._R_P/128._R_P  ! stencil 4
        c(5,5)=  693._R_P/256._R_P; c(4,5)=-1155._R_P/256._R_P; c(3,5)=  693._R_P/128._R_P  ! stencil 5
        !   cell  2               ;    cell  1                ;    cell  0
        c(2,0)= -105._R_P/128._R_P; c(1,0)=  315._R_P/256._R_P; c(0,0)=   63._R_P/256._R_P  ! stencil 0
        c(2,1)=  105._R_P/128._R_P; c(1,1)=  105._R_P/256._R_P; c(0,1)=   -7._R_P/256._R_P  ! stencil 1
        c(2,2)=   75._R_P/128._R_P; c(1,2)=  -25._R_P/256._R_P; c(0,2)=    3._R_P/256._R_P  ! stencil 2
        c(2,3)=  -35._R_P/128._R_P; c(1,3)=   21._R_P/256._R_P; c(0,3)=   -3._R_P/256._R_P  ! stencil 3
        c(2,4)=   63._R_P/128._R_P; c(1,4)=  -45._R_P/256._R_P; c(0,4)=    7._R_P/256._R_P  ! stencil 4
        c(2,5)= -495._R_P/128._R_P; c(1,5)=  385._R_P/256._R_P; c(0,5)=  -63._R_P/256._R_P  ! stencil 5
      case(7) ! 13th order
        !  cell  6                ;    cell  5                ;    cell  4
        c(6,0)=  -21._R_P/1024._R_P; c(5,0)=   77._R_P/512._R_P ; c(4,0)= -495._R_P/1024._R_P  ! stencil 0
        c(6,1)=    7._R_P/1024._R_P; c(5,1)=  -27._R_P/512._R_P ; c(4,1)=  189._R_P/1024._R_P  ! stencil 1
        c(6,2)=   -5._R_P/1024._R_P; c(5,2)=   21._R_P/512._R_P ; c(4,2)= -175._R_P/1024._R_P  ! stencil 2
        c(6,3)=    7._R_P/1024._R_P; c(5,3)=  -35._R_P/512._R_P ; c(4,3)=  525._R_P/1024._R_P  ! stencil 3
        c(6,4)=  -21._R_P/1024._R_P; c(5,4)=  189._R_P/512._R_P ; c(4,4)=  945._R_P/1024._R_P  ! stencil 4
        c(6,5)=  231._R_P/1024._R_P; c(5,5)=  693._R_P/512._R_P ; c(4,5)=-1155._R_P/1024._R_P  ! stencil 5
        c(6,6)= 3003._R_P/1024._R_P; c(5,6)=-3003._R_P/512._R_P ; c(4,6)= 9009._R_P/1024._R_P  ! stencil 6
        !  cell  3                 ;   cell  2                ;    cell  1
        c(3,0)=  231._R_P/256._R_P ; c(2,0)=-1155._R_P/1024._R_P; c(1,0)=  693._R_P/512._R_P   ! stencil 0
        c(3,1)= -105._R_P/256._R_P ; c(2,1)=  945._R_P/1024._R_P; c(1,1)=  189._R_P/512._R_P   ! stencil 1
        c(3,2)=  175._R_P/256._R_P ; c(2,2)=  525._R_P/1024._R_P; c(1,2)=  -35._R_P/512._R_P   ! stencil 2
        c(3,3)=  175._R_P/256._R_P ; c(2,3)= -175._R_P/1024._R_P; c(1,3)=   21._R_P/512._R_P   ! stencil 3
        c(3,4)= -105._R_P/256._R_P ; c(2,4)=  189._R_P/1024._R_P; c(1,4)=  -27._R_P/512._R_P   ! stencil 4
        c(3,5)=  231._R_P/256._R_P ; c(2,5)= -495._R_P/1024._R_P; c(1,5)=   77._R_P/512._R_P   ! stencil 5
        c(3,6)=-2145._R_P/256._R_P ; c(2,6)= 5005._R_P/1024._R_P; c(1,6)= -819._R_P/512._R_P   ! stencil 6
        !  cell  0
        c(0,0)=  231._R_P/1024._R_P  ! stencil 0
        c(0,1)=  -21._R_P/1024._R_P  ! stencil 1
        c(0,2)=    7._R_P/1024._R_P  ! stencil 2
        c(0,3)=   -5._R_P/1024._R_P  ! stencil 3
        c(0,4)=    7._R_P/1024._R_P  ! stencil 4
        c(0,5)=  -21._R_P/1024._R_P  ! stencil 5
        c(0,6)=  231._R_P/1024._R_P  ! stencil 6
      case(8) ! 15th order
        !  cell  7                ;    cell  6                ;    cell  5
        c(7,0)=    33._R_P/2048._R_P; c(6,0)=  -273._R_P/2048._R_P; c(5,0)=  1001._R_P/2048._R_P  ! stencil 0
        c(7,1)=    -9._R_P/2048._R_P; c(6,1)=    77._R_P/2048._R_P; c(5,1)=  -297._R_P/2048._R_P  ! stencil 1
        c(7,2)=     5._R_P/2048._R_P; c(6,2)=   -45._R_P/2048._R_P; c(5,2)=   189._R_P/2048._R_P  ! stencil 2
        c(7,3)=    -5._R_P/2048._R_P; c(6,3)=    49._R_P/2048._R_P; c(5,3)=  -245._R_P/2048._R_P  ! stencil 3
        c(7,4)=     9._R_P/2048._R_P; c(6,4)=  -105._R_P/2048._R_P; c(5,4)=   945._R_P/2048._R_P  ! stencil 4
        c(7,5)=   -33._R_P/2048._R_P; c(6,5)=   693._R_P/2048._R_P; c(5,5)=  2079._R_P/2048._R_P  ! stencil 5
        c(7,6)=   429._R_P/2048._R_P; c(6,6)=  3003._R_P/2048._R_P; c(5,6)= -3003._R_P/2048._R_P  ! stencil 6
        c(7,7)=  6435._R_P/2048._R_P; c(6,7)=-15015._R_P/2048._R_P; c(5,7)= 27027._R_P/2048._R_P  ! stencil 7
        !  cell  4                ;    cell  3                ;    cell  2
        c(4,0)= -2145._R_P/2048._R_P; c(3,0)=  3003._R_P/2048._R_P; c(2,0)= -3003._R_P/2048._R_P  ! stencil 0
        c(4,1)=   693._R_P/2048._R_P; c(3,1)= -1155._R_P/2048._R_P; c(2,1)=  2079._R_P/2048._R_P  ! stencil 1
        c(4,2)=  -525._R_P/2048._R_P; c(3,2)=  1575._R_P/2048._R_P; c(2,2)=   945._R_P/2048._R_P  ! stencil 2
        c(4,3)=  1225._R_P/2048._R_P; c(3,3)=  1225._R_P/2048._R_P; c(2,3)=  -245._R_P/2048._R_P  ! stencil 3
        c(4,4)=  1575._R_P/2048._R_P; c(3,4)=  -525._R_P/2048._R_P; c(2,4)=   189._R_P/2048._R_P  ! stencil 4
        c(4,5)= -1155._R_P/2048._R_P; c(3,5)=   693._R_P/2048._R_P; c(2,5)=  -297._R_P/2048._R_P  ! stencil 5
        c(4,6)=  3003._R_P/2048._R_P; c(3,6)= -2145._R_P/2048._R_P; c(2,6)=  1001._R_P/2048._R_P  ! stencil 6
        c(4,7)=-32175._R_P/2048._R_P; c(3,7)= 25025._R_P/2048._R_P; c(2,7)=-12285._R_P/2048._R_P  ! stencil 7
        !  cell  1                ;    cell  0
        c(1,0)=  3003._R_P/2048._R_P; c(0,0)=   429._R_P/2048._R_P  ! stencil 0
        c(1,1)=   693._R_P/2048._R_P; c(0,1)=   -33._R_P/2048._R_P  ! stencil 1
        c(1,2)=  -105._R_P/2048._R_P; c(0,2)=     9._R_P/2048._R_P  ! stencil 2
        c(1,3)=    49._R_P/2048._R_P; c(0,3)=    -5._R_P/2048._R_P  ! stencil 3
        c(1,4)=   -45._R_P/2048._R_P; c(0,4)=     5._R_P/2048._R_P  ! stencil 4
        c(1,5)=    77._R_P/2048._R_P; c(0,5)=    -9._R_P/2048._R_P  ! stencil 5
        c(1,6)=  -273._R_P/2048._R_P; c(0,6)=    33._R_P/2048._R_P  ! stencil 6
        c(1,7)=  3465._R_P/2048._R_P; c(0,7)=  -429._R_P/2048._R_P  ! stencil 7
      case(9) ! 17th order
        !  cell  8                    ;     cell  7                   ;     cell  6
        c(8,0)=   -429._R_P/32768._R_P; c(7,0)=    495._R_P/4096._R_P ; c(6,0)=  -4095._R_P/8192._R_P   ! stencil 0
        c(8,1)=     99._R_P/32768._R_P; c(7,1)=   -117._R_P/4096._R_P ; c(6,1)=   1001._R_P/8192._R_P   ! stencil 1
        c(8,2)=    -45._R_P/32768._R_P; c(7,2)=     55._R_P/4096._R_P ; c(6,2)=   -495._R_P/8192._R_P   ! stencil 2
        c(8,3)=     35._R_P/32768._R_P; c(7,3)=    -45._R_P/4096._R_P ; c(6,3)=    441._R_P/8192._R_P   ! stencil 3
        c(8,4)=    -45._R_P/32768._R_P; c(7,4)=     63._R_P/4096._R_P ; c(6,4)=   -735._R_P/8192._R_P   ! stencil 4
        c(8,5)=     99._R_P/32768._R_P; c(7,5)=   -165._R_P/4096._R_P ; c(6,5)=   3465._R_P/8192._R_P   ! stencil 5
        c(8,6)=   -429._R_P/32768._R_P; c(7,6)=   1287._R_P/4096._R_P ; c(6,6)=   9009._R_P/8192._R_P   ! stencil 6
        c(8,7)=   6435._R_P/32768._R_P; c(7,7)=   6435._R_P/4096._R_P ; c(6,7)= -15015._R_P/8192._R_P   ! stencil 7
        c(8,8)= 109395._R_P/32768._R_P; c(7,8)= -36465._R_P/4096._R_P ; c(6,8)= 153153._R_P/8192._R_P   ! stencil 8
        !  cell  5                    ;    cell  4                    ;     cell  3
        c(5,0)=   5005._R_P/4096._R_P ; c(4,0)= -32175._R_P/16384._R_P; c(3,0)=   9009._R_P/4096._R_P   ! stencil 0
        c(5,1)=  -1287._R_P/4096._R_P ; c(4,1)=   9009._R_P/16384._R_P; c(3,1)=  -3003._R_P/4096._R_P   ! stencil 1
        c(5,2)=    693._R_P/4096._R_P ; c(4,2)=  -5775._R_P/16384._R_P; c(3,2)=   3465._R_P/4096._R_P   ! stencil 2
        c(5,3)=   -735._R_P/4096._R_P ; c(4,3)=  11025._R_P/16384._R_P; c(3,3)=   2205._R_P/4096._R_P   ! stencil 3
        c(5,4)=   2205._R_P/4096._R_P ; c(4,4)=  11025._R_P/16384._R_P; c(3,4)=   -735._R_P/4096._R_P   ! stencil 4
        c(5,5)=   3465._R_P/4096._R_P ; c(4,5)=  -5775._R_P/16384._R_P; c(3,5)=    693._R_P/4096._R_P   ! stencil 5
        c(5,6)=  -3003._R_P/4096._R_P ; c(4,6)=   9009._R_P/16384._R_P; c(3,6)=  -1287._R_P/4096._R_P   ! stencil 6
        c(5,7)=   9009._R_P/4096._R_P ; c(4,7)= -32175._R_P/16384._R_P; c(3,7)=   5005._R_P/4096._R_P   ! stencil 7
        c(5,8)=-109395._R_P/4096._R_P ; c(4,8)= 425425._R_P/16384._R_P; c(3,8)= -69615._R_P/4096._R_P   ! stencil 8
        !   cell  2                   ;     cell  1                   ;     cell  0
        c(2,0)= -15015._R_P/8192._R_P ; c(1,0)=   6435._R_P/4096._R_P ; c(0,0)=   6435._R_P/32768._R_P  ! stencil 0
        c(2,1)=   9009._R_P/8192._R_P ; c(1,1)=   1287._R_P/4096._R_P ; c(0,1)=   -429._R_P/32768._R_P  ! stencil 1
        c(2,2)=   3465._R_P/8192._R_P ; c(1,2)=   -165._R_P/4096._R_P ; c(0,2)=     99._R_P/32768._R_P  ! stencil 2
        c(2,3)=   -735._R_P/8192._R_P ; c(1,3)=     63._R_P/4096._R_P ; c(0,3)=    -45._R_P/32768._R_P  ! stencil 3
        c(2,4)=    441._R_P/8192._R_P ; c(1,4)=    -45._R_P/4096._R_P ; c(0,4)=     35._R_P/32768._R_P  ! stencil 4
        c(2,5)=   -495._R_P/8192._R_P ; c(1,5)=     55._R_P/4096._R_P ; c(0,5)=    -45._R_P/32768._R_P  ! stencil 5
        c(2,6)=   1001._R_P/8192._R_P ; c(1,6)=   -117._R_P/4096._R_P ; c(0,6)=     99._R_P/32768._R_P  ! stencil 6
        c(2,7)=  -4095._R_P/8192._R_P ; c(1,7)=    495._R_P/4096._R_P ; c(0,7)=   -429._R_P/32768._R_P  ! stencil 7
        c(2,8)=  58905._R_P/8192._R_P ; c(1,8)=  -7293._R_P/4096._R_P ; c(0,8)=   6435._R_P/32768._R_P  ! stencil 8
    endselect
  elseif(x_target==0.5_R_P) then
    ! right interface (i+1/2)
    select case(S)
      case(2) ! 3rd order
        !  cell  1      ;    cell  0
        c(1,0)= -0.5_R_P; c(0,0)=  1.5_R_P ! stencil 0
        c(1,1)=  0.5_R_P; c(0,1)=  0.5_R_P ! stencil 1
      case(3) ! 5th order
        !  cell  2            ;    cell  1            ;    cell  0
        c(2,0)=  3._R_P/8._R_P; c(1,0)= -5._R_P/4._R_P; c(0,0)= 15._R_P/8._R_P ! stencil 0
        c(2,1)= -1._R_P/8._R_P; c(1,1)=  3._R_P/4._R_P; c(0,1)=  3._R_P/8._R_P ! stencil 1
        c(2,2)=  3._R_P/8._R_P; c(1,2)=  3._R_P/4._R_P; c(0,2)= -1._R_P/8._R_P ! stencil 2
      case(4) ! 7th order
        !  cell  3             ;    cell  2             ;   cell  1              ;    cell  0
        c(3,0)= -5._R_P/16._R_P; c(2,0)= 21._R_P/16._R_P; c(1,0)=-35._R_P/16._R_P; c(0,0)= 35._R_P/16._R_P ! stencil 0
        c(3,1)=  1._R_P/16._R_P; c(2,1)= -5._R_P/16._R_P; c(1,1)= 15._R_P/16._R_P; c(0,1)=  5._R_P/16._R_P ! stencil 1
        c(3,2)= -1._R_P/16._R_P; c(2,2)=  9._R_P/16._R_P; c(1,2)=  9._R_P/16._R_P; c(0,2)= -1._R_P/16._R_P ! stencil 2
        c(3,3)=  5._R_P/16._R_P; c(2,3)= 15._R_P/16._R_P; c(1,3)= -5._R_P/16._R_P; c(0,3)=  1._R_P/16._R_P ! stencil 3
      case(5) ! 9th order
        !  cell  4               ;    cell  3               ;   cell  2                ;    cell  1
        c(4,0)=  35._R_P/128._R_P; c(3,0)= -45._R_P/32._R_P ; c(2,0)= 189._R_P/64._R_P ; c(1,0)=-105._R_P/32._R_P  ! stencil 0
        c(4,1)=  -5._R_P/128._R_P; c(3,1)=   7._R_P/32._R_P ; c(2,1)= -35._R_P/64._R_P ; c(1,1)=  35._R_P/32._R_P  ! stencil 1
        c(4,2)=   3._R_P/128._R_P; c(3,2)=  -5._R_P/32._R_P ; c(2,2)=  45._R_P/64._R_P ; c(1,2)=  15._R_P/32._R_P  ! stencil 2
        c(4,3)=  -5._R_P/128._R_P; c(3,3)=  15._R_P/32._R_P ; c(2,3)=  45._R_P/64._R_P ; c(1,3)=  -5._R_P/32._R_P  ! stencil 3
        c(4,4)=  35._R_P/128._R_P; c(3,4)=  35._R_P/32._R_P ; c(2,4)= -35._R_P/64._R_P ; c(1,4)=   7._R_P/32._R_P  ! stencil 4
        !  cell  0
        c(0,0)= 315._R_P/128._R_P ! stencil 0
        c(0,1)=  35._R_P/128._R_P ! stencil 1
        c(0,2)=  -5._R_P/128._R_P ! stencil 2
        c(0,3)=   3._R_P/128._R_P ! stencil 3
        c(0,4)=  -5._R_P/128._R_P ! stencil 4
      case(6) ! 11th order
        !  cell  5                ;    cell  4                ;   cell  3
        c(5,0)=  -63._R_P/256._R_P; c(4,0)=  385._R_P/256._R_P; c(3,0)= -495._R_P/128._R_P  ! stencil 0
        c(5,1)=    7._R_P/256._R_P; c(4,1)=  -45._R_P/256._R_P; c(3,1)=   63._R_P/128._R_P  ! stencil 1
        c(5,2)=   -3._R_P/256._R_P; c(4,2)=   21._R_P/256._R_P; c(3,2)=  -35._R_P/128._R_P  ! stencil 2
        c(5,3)=    3._R_P/256._R_P; c(4,3)=  -25._R_P/256._R_P; c(3,3)=   75._R_P/128._R_P  ! stencil 3
        c(5,4)=   -7._R_P/256._R_P; c(4,4)=  105._R_P/256._R_P; c(3,4)=  105._R_P/128._R_P  ! stencil 4
        c(5,5)=   63._R_P/256._R_P; c(4,5)=  315._R_P/256._R_P; c(3,5)= -105._R_P/128._R_P  ! stencil 5
        !  cell  2                ;    cell  1                ;   cell  0
        c(2,0)=  693._R_P/128._R_P; c(1,0)=-1155._R_P/256._R_P; c(0,0)=  693._R_P/256._R_P  ! stencil 0
        c(2,1)= -105._R_P/128._R_P; c(1,1)=  315._R_P/256._R_P; c(0,1)=   63._R_P/256._R_P  ! stencil 1
        c(2,2)=  105._R_P/128._R_P; c(1,2)=  105._R_P/256._R_P; c(0,2)=   -7._R_P/256._R_P  ! stencil 2
        c(2,3)=   75._R_P/128._R_P; c(1,3)=  -25._R_P/256._R_P; c(0,3)=    3._R_P/256._R_P  ! stencil 3
        c(2,4)=  -35._R_P/128._R_P; c(1,4)=   21._R_P/256._R_P; c(0,4)=   -3._R_P/256._R_P  ! stencil 4
        c(2,5)=   63._R_P/128._R_P; c(1,5)=  -45._R_P/256._R_P; c(0,5)=    7._R_P/256._R_P  ! stencil 5
      case(7) ! 13th order
        !  cell  6                ;    cell  5                ;    cell  4
        c(6,0)=  231._R_P/1024._R_P; c(5,0)= -819._R_P/512._R_P ; c(4,0)= 5005._R_P/1024._R_P  ! stencil 0
        c(6,1)=  -21._R_P/1024._R_P; c(5,1)=   77._R_P/512._R_P ; c(4,1)= -495._R_P/1024._R_P  ! stencil 1
        c(6,2)=    7._R_P/1024._R_P; c(5,2)=  -27._R_P/512._R_P ; c(4,2)=  189._R_P/1024._R_P  ! stencil 2
        c(6,3)=   -5._R_P/1024._R_P; c(5,3)=   21._R_P/512._R_P ; c(4,3)= -175._R_P/1024._R_P  ! stencil 3
        c(6,4)=    7._R_P/1024._R_P; c(5,4)=  -35._R_P/512._R_P ; c(4,4)=  525._R_P/1024._R_P  ! stencil 4
        c(6,5)=  -21._R_P/1024._R_P; c(5,5)=  189._R_P/512._R_P ; c(4,5)=  945._R_P/1024._R_P  ! stencil 5
        c(6,6)=  231._R_P/1024._R_P; c(5,6)=  693._R_P/512._R_P ; c(4,6)=-1155._R_P/1024._R_P  ! stencil 6
        !  cell  3                ;    cell  2                ;    cell  1
        c(3,0)=-2145._R_P/256._R_P ; c(2,0)= 9009._R_P/1024._R_P; c(1,0)=-3003._R_P/512._R_P   ! stencil 0
        c(3,1)=  231._R_P/256._R_P ; c(2,1)=-1155._R_P/1024._R_P; c(1,1)=  693._R_P/512._R_P   ! stencil 1
        c(3,2)= -105._R_P/256._R_P ; c(2,2)=  945._R_P/1024._R_P; c(1,2)=  189._R_P/512._R_P   ! stencil 2
        c(3,3)=  175._R_P/256._R_P ; c(2,3)=  525._R_P/1024._R_P; c(1,3)=  -35._R_P/512._R_P   ! stencil 3
        c(3,4)=  175._R_P/256._R_P ; c(2,4)= -175._R_P/1024._R_P; c(1,4)=   21._R_P/512._R_P   ! stencil 4
        c(3,5)= -105._R_P/256._R_P ; c(2,5)=  189._R_P/1024._R_P; c(1,5)=  -27._R_P/512._R_P   ! stencil 5
        c(3,6)=  231._R_P/256._R_P ; c(2,6)= -495._R_P/1024._R_P; c(1,6)=   77._R_P/512._R_P   ! stencil 6
        !  cell  0
        c(0,0)= 3003._R_P/1024._R_P  ! stencil 0
        c(0,1)=  231._R_P/1024._R_P  ! stencil 1
        c(0,2)=  -21._R_P/1024._R_P  ! stencil 2
        c(0,3)=    7._R_P/1024._R_P  ! stencil 3
        c(0,4)=   -5._R_P/1024._R_P  ! stencil 4
        c(0,5)=    7._R_P/1024._R_P  ! stencil 5
        c(0,6)=  -21._R_P/1024._R_P  ! stencil 6
      case(8) ! 15th order
        !  cell  7                ;    cell  6                ;    cell  5
        c(7,0)=  -429._R_P/2048._R_P; c(6,0)=  3465._R_P/2048._R_P; c(5,0)=-12285._R_P/2048._R_P  ! stencil 0
        c(7,1)=    33._R_P/2048._R_P; c(6,1)=  -273._R_P/2048._R_P; c(5,1)=  1001._R_P/2048._R_P  ! stencil 1
        c(7,2)=    -9._R_P/2048._R_P; c(6,2)=    77._R_P/2048._R_P; c(5,2)=  -297._R_P/2048._R_P  ! stencil 2
        c(7,3)=     5._R_P/2048._R_P; c(6,3)=   -45._R_P/2048._R_P; c(5,3)=   189._R_P/2048._R_P  ! stencil 3
        c(7,4)=    -5._R_P/2048._R_P; c(6,4)=    49._R_P/2048._R_P; c(5,4)=  -245._R_P/2048._R_P  ! stencil 4
        c(7,5)=     9._R_P/2048._R_P; c(6,5)=  -105._R_P/2048._R_P; c(5,5)=   945._R_P/2048._R_P  ! stencil 5
        c(7,6)=   -33._R_P/2048._R_P; c(6,6)=   693._R_P/2048._R_P; c(5,6)=  2079._R_P/2048._R_P  ! stencil 6
        c(7,7)=   429._R_P/2048._R_P; c(6,7)=  3003._R_P/2048._R_P; c(5,7)= -3003._R_P/2048._R_P  ! stencil 7
        !  cell  4                ;    cell  3                ;    cell  2
        c(4,0)= 25025._R_P/2048._R_P; c(3,0)=-32175._R_P/2048._R_P; c(2,0)= 27027._R_P/2048._R_P  ! stencil 0
        c(4,1)= -2145._R_P/2048._R_P; c(3,1)=  3003._R_P/2048._R_P; c(2,1)= -3003._R_P/2048._R_P  ! stencil 1
        c(4,2)=   693._R_P/2048._R_P; c(3,2)= -1155._R_P/2048._R_P; c(2,2)=  2079._R_P/2048._R_P  ! stencil 2
        c(4,3)=  -525._R_P/2048._R_P; c(3,3)=  1575._R_P/2048._R_P; c(2,3)=   945._R_P/2048._R_P  ! stencil 3
        c(4,4)=  1225._R_P/2048._R_P; c(3,4)=  1225._R_P/2048._R_P; c(2,4)=  -245._R_P/2048._R_P  ! stencil 4
        c(4,5)=  1575._R_P/2048._R_P; c(3,5)=  -525._R_P/2048._R_P; c(2,5)=   189._R_P/2048._R_P  ! stencil 5
        c(4,6)= -1155._R_P/2048._R_P; c(3,6)=   693._R_P/2048._R_P; c(2,6)=  -297._R_P/2048._R_P  ! stencil 6
        c(4,7)=  3003._R_P/2048._R_P; c(3,7)= -2145._R_P/2048._R_P; c(2,7)=  1001._R_P/2048._R_P  ! stencil 7
        !  cell  1                  ;    cell  0
        c(1,0)=-15015._R_P/2048._R_P; c(0,0)=  6435._R_P/2048._R_P  ! stencil 0
        c(1,1)=  3003._R_P/2048._R_P; c(0,1)=   429._R_P/2048._R_P  ! stencil 1
        c(1,2)=   693._R_P/2048._R_P; c(0,2)=   -33._R_P/2048._R_P  ! stencil 2
        c(1,3)=  -105._R_P/2048._R_P; c(0,3)=     9._R_P/2048._R_P  ! stencil 3
        c(1,4)=    49._R_P/2048._R_P; c(0,4)=    -5._R_P/2048._R_P  ! stencil 4
        c(1,5)=   -45._R_P/2048._R_P; c(0,5)=     5._R_P/2048._R_P  ! stencil 5
        c(1,6)=    77._R_P/2048._R_P; c(0,6)=    -9._R_P/2048._R_P  ! stencil 6
        c(1,7)=  -273._R_P/2048._R_P; c(0,7)=    33._R_P/2048._R_P  ! stencil 7
      case(9) ! 17th order
        !  cell  8                  ;     cell  7                 ;     cell  6
        c(8,0)=   6435._R_P/32768._R_P; c(7,0)=  -7293._R_P/ 4096._R_P; c(6,0)=  58905._R_P/ 8192._R_P  ! stencil 0
        c(8,1)=   -429._R_P/32768._R_P; c(7,1)=    495._R_P/ 4096._R_P; c(6,1)=  -4095._R_P/ 8192._R_P  ! stencil 1
        c(8,2)=     99._R_P/32768._R_P; c(7,2)=   -117._R_P/ 4096._R_P; c(6,2)=   1001._R_P/ 8192._R_P  ! stencil 2
        c(8,3)=    -45._R_P/32768._R_P; c(7,3)=     55._R_P/ 4096._R_P; c(6,3)=   -495._R_P/ 8192._R_P  ! stencil 3
        c(8,4)=     35._R_P/32768._R_P; c(7,4)=    -45._R_P/ 4096._R_P; c(6,4)=    441._R_P/ 8192._R_P  ! stencil 4
        c(8,5)=    -45._R_P/32768._R_P; c(7,5)=     63._R_P/ 4096._R_P; c(6,5)=   -735._R_P/ 8192._R_P  ! stencil 5
        c(8,6)=     99._R_P/32768._R_P; c(7,6)=   -165._R_P/ 4096._R_P; c(6,6)=   3465._R_P/ 8192._R_P  ! stencil 6
        c(8,7)=   -429._R_P/32768._R_P; c(7,7)=   1287._R_P/ 4096._R_P; c(6,7)=   9009._R_P/ 8192._R_P  ! stencil 7
        c(8,8)=   6435._R_P/32768._R_P; c(7,8)=   6435._R_P/ 4096._R_P; c(6,8)= -15015._R_P/ 8192._R_P  ! stencil 8
        !   cell  5                 ; !  cell  4                  ;     cell  3
        c(5,0)= -69615._R_P/ 4096._R_P; c(4,0)= 425425._R_P/16384._R_P; c(3,0)=-109395._R_P/ 4096._R_P  ! stencil 0
        c(5,1)=   5005._R_P/ 4096._R_P; c(4,1)= -32175._R_P/16384._R_P; c(3,1)=   9009._R_P/ 4096._R_P  ! stencil 1
        c(5,2)=  -1287._R_P/ 4096._R_P; c(4,2)=   9009._R_P/16384._R_P; c(3,2)=  -3003._R_P/ 4096._R_P  ! stencil 2
        c(5,3)=    693._R_P/ 4096._R_P; c(4,3)=  -5775._R_P/16384._R_P; c(3,3)=   3465._R_P/ 4096._R_P  ! stencil 3
        c(5,4)=   -735._R_P/ 4096._R_P; c(4,4)=  11025._R_P/16384._R_P; c(3,4)=   2205._R_P/ 4096._R_P  ! stencil 4
        c(5,5)=   2205._R_P/ 4096._R_P; c(4,5)=  11025._R_P/16384._R_P; c(3,5)=   -735._R_P/ 4096._R_P  ! stencil 5
        c(5,6)=   3465._R_P/ 4096._R_P; c(4,6)=  -5775._R_P/16384._R_P; c(3,6)=    693._R_P/ 4096._R_P  ! stencil 6
        c(5,7)=  -3003._R_P/ 4096._R_P; c(4,7)=   9009._R_P/16384._R_P; c(3,7)=  -1287._R_P/ 4096._R_P  ! stencil 7
        c(5,8)=   9009._R_P/ 4096._R_P; c(4,8)= -32175._R_P/16384._R_P; c(3,8)=   5005._R_P/ 4096._R_P  ! stencil 8
        !   cell  2                 ;     cell  1                 ;     cell  0
        c(2,0)= 153153._R_P/ 8192._R_P; c(1,0)= -36465._R_P/ 4096._R_P; c(0,0)= 109395._R_P/32768._R_P  ! stencil 0
        c(2,1)= -15015._R_P/ 8192._R_P; c(1,1)=   6435._R_P/ 4096._R_P; c(0,1)=   6435._R_P/32768._R_P  ! stencil 1
        c(2,2)=   9009._R_P/ 8192._R_P; c(1,2)=   1287._R_P/ 4096._R_P; c(0,2)=   -429._R_P/32768._R_P  ! stencil 2
        c(2,3)=   3465._R_P/ 8192._R_P; c(1,3)=   -165._R_P/ 4096._R_P; c(0,3)=     99._R_P/32768._R_P  ! stencil 3
        c(2,4)=   -735._R_P/ 8192._R_P; c(1,4)=     63._R_P/ 4096._R_P; c(0,4)=    -45._R_P/32768._R_P  ! stencil 4
        c(2,5)=    441._R_P/ 8192._R_P; c(1,5)=    -45._R_P/ 4096._R_P; c(0,5)=     35._R_P/32768._R_P  ! stencil 5
        c(2,6)=   -495._R_P/ 8192._R_P; c(1,6)=     55._R_P/ 4096._R_P; c(0,6)=    -45._R_P/32768._R_P  ! stencil 6
        c(2,7)=   1001._R_P/ 8192._R_P; c(1,7)=   -117._R_P/ 4096._R_P; c(0,7)=     99._R_P/32768._R_P  ! stencil 7
        c(2,8)=  -4095._R_P/ 8192._R_P; c(1,8)=    495._R_P/ 4096._R_P; c(0,8)=   -429._R_P/32768._R_P  ! stencil 8
    endselect
  else
    ! internal point
    do k=0,S-1  !stencils loop
      c_sum = 0._R_P
      do j=0,S-2  !values loop
        prod = 1._R_P
        do i=0,S-1
          if (i==j) cycle
          prod = prod * ((x_target - (-S+k+i+1)) / ((-S+k+j+1) - (-S+k+i+1)))
        enddo
        c(S-1-j,k) = prod
        c_sum = c_sum + prod
      enddo
      c(0,k) = 1._R_P - c_sum
    enddo
  endif
  endsubroutine assign_interp_coeff

endmodule wenoof_interpolations_int_js
