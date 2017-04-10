!< Jiang-Shu (Lagrange) interpolations object for function interpolation.
module wenoof_interpolations_int_js
!< Jiang-Shu (Lagrange) interpolations object for function interpolation.
!<
!< @note The provided interpolations implement the Lagrange interpolations defined in *High Order Weighted Essentially
!< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
!< doi:10.1137/070679065.

#ifdef r16p
use penf, only : I_P, RPP=>R16P, str
#else
use penf, only : I_P, RPP=>R8P, str
#endif
use wenoof_base_object, only : base_object_constructor
use wenoof_interpolations_object, only : interpolations_object, interpolations_object_constructor

implicit none
private
public :: interpolations_int_js
public :: interpolations_int_js_constructor

type, extends(interpolations_object_constructor) :: interpolations_int_js_constructor
  !< Jiang-Shu (Lagrange) interpolations object for function interpolation constructor.
endtype interpolations_int_js_constructor

type, extends(interpolations_object) :: interpolations_int_js
  !< Jiang-Shu (Lagrange) interpolations object for function interpolation.
  !<
  !< @note The provided interpolations implement the Lagrange interpolations defined in *High Order Weighted Essentially
  !< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
  !< doi:10.1137/070679065.
  real(RPP), allocatable :: coef(:,:) !< Polynomial coefficients [0:S-1,0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: create      !< Create interpolations.
    procedure, pass(self) :: compute_int !< Compute interpolations (interpolate).
    procedure, pass(self) :: compute_rec !< Compute interpolations (reconstruct).
    procedure, pass(self) :: description !< Return object string-description.
    procedure, pass(self) :: destroy     !< Destroy interpolations.
endtype interpolations_int_js

contains
  ! public deferred methods
  subroutine create(self, constructor)
  !< Create interpolations.
  class(interpolations_int_js),   intent(inout) :: self        !< Interpolations.
  class(base_object_constructor), intent(in)    :: constructor !< Interpolations constructor.
  real(RPP), allocatable                        :: f(:,:)      !< Temporary variable.
  real(RPP)                                     :: prod        !< Temporary variable.
  real(RPP)                                     :: c_sum       !< Temporary variable.
  integer(I_P)                                  :: i, j, k     !< Counters.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%coef(0:self%S - 1, 0:self%S - 1))
  select type(constructor)
  type is(interpolations_int_js_constructor)
    associate(S => self%S, c => self%coef, stencil => constructor%stencil, x_target => constructor%x_target)
      if(x_target==-0.5_RPP) then
        ! left interface (i-1/2)
        select case(S)
          case(2) ! 3rd order
            !  cell  1      ;    cell  0
            c(1,0)=  0.5_RPP; c(0,0)=  0.5_RPP ! stencil 0
            c(1,1)=  1.5_RPP; c(0,1)= -0.5_RPP ! stencil 1
          case(3) ! 5th order
            !  cell  2            ;    cell  1            ;    cell  0
            c(2,0)= -1._RPP/8._RPP; c(1,0)=  3._RPP/4._RPP; c(0,0)=  3._RPP/8._RPP ! stencil 0
            c(2,1)=  3._RPP/8._RPP; c(1,1)=  3._RPP/4._RPP; c(0,1)= -1._RPP/8._RPP ! stencil 1
            c(2,2)= 15._RPP/8._RPP; c(1,2)= -5._RPP/4._RPP; c(0,2)=  3._RPP/8._RPP ! stencil 2
          case(4) ! 7th order
            !  cell  3             ;    cell  2             ;    cell  1             ;    cell  0
            c(3,0)=  1._RPP/16._RPP; c(2,0)= -5._RPP/16._RPP; c(1,0)= 15._RPP/16._RPP; c(0,0)=  5._RPP/16._RPP ! stencil 0
            c(3,1)= -1._RPP/16._RPP; c(2,1)=  9._RPP/16._RPP; c(1,1)=  9._RPP/16._RPP; c(0,1)= -1._RPP/16._RPP ! stencil 1
            c(3,2)=  5._RPP/16._RPP; c(2,2)= 15._RPP/16._RPP; c(1,2)= -5._RPP/16._RPP; c(0,2)=  1._RPP/16._RPP ! stencil 2
            c(3,3)= 35._RPP/16._RPP; c(2,3)=-35._RPP/16._RPP; c(1,3)= 21._RPP/16._RPP; c(0,3)= -5._RPP/16._RPP ! stencil 3
          case(5) ! 9th order
            !  cell  4               ;    cell  3               ;    cell  2               ;    cell  1
            c(4,0)=  -5._RPP/128._RPP; c(3,0)=   7._RPP/32._RPP ; c(2,0)= -35._RPP/64._RPP ; c(1,0)=  35._RPP/32._RPP  ! stencil 0
            c(4,1)=   3._RPP/128._RPP; c(3,1)=  -5._RPP/32._RPP ; c(2,1)=  45._RPP/64._RPP ; c(1,1)=  15._RPP/32._RPP  ! stencil 1
            c(4,2)=  -5._RPP/128._RPP; c(3,2)=  15._RPP/32._RPP ; c(2,2)=  45._RPP/64._RPP ; c(1,2)=  -5._RPP/32._RPP  ! stencil 2
            c(4,3)=  35._RPP/128._RPP; c(3,3)=  35._RPP/32._RPP ; c(2,3)= -35._RPP/64._RPP ; c(1,3)=   7._RPP/32._RPP  ! stencil 3
            c(4,4)= 315._RPP/128._RPP; c(3,4)=-105._RPP/32._RPP ; c(2,4)= 189._RPP/64._RPP ; c(1,4)= -45._RPP/32._RPP  ! stencil 4
            !  cell  0
            c(0,0)=  35._RPP/128._RPP  ! stencil 0
            c(0,1)=  -5._RPP/128._RPP  ! stencil 1
            c(0,2)=   3._RPP/128._RPP  ! stencil 2
            c(0,3)=  -5._RPP/128._RPP  ! stencil 3
            c(0,4)=  35._RPP/128._RPP  ! stencil 4
          case(6) ! 11th order
            !  cell  5                ;    cell  4                ;    cell  3
            c(5,0)=    7._RPP/256._RPP; c(4,0)=  -45._RPP/256._RPP; c(3,0)=   63._RPP/128._RPP  ! stencil 0
            c(5,1)=   -3._RPP/256._RPP; c(4,1)=   21._RPP/256._RPP; c(3,1)=  -35._RPP/128._RPP  ! stencil 1
            c(5,2)=    3._RPP/256._RPP; c(4,2)=  -25._RPP/256._RPP; c(3,2)=   75._RPP/128._RPP  ! stencil 2
            c(5,3)=   -7._RPP/256._RPP; c(4,3)=  105._RPP/256._RPP; c(3,3)=  105._RPP/128._RPP  ! stencil 3
            c(5,4)=   63._RPP/256._RPP; c(4,4)=  315._RPP/256._RPP; c(3,4)= -105._RPP/128._RPP  ! stencil 4
            c(5,5)=  693._RPP/256._RPP; c(4,5)=-1155._RPP/256._RPP; c(3,5)=  693._RPP/128._RPP  ! stencil 5
            !   cell  2               ;    cell  1                ;    cell  0
            c(2,0)= -105._RPP/128._RPP; c(1,0)=  315._RPP/256._RPP; c(0,0)=   63._RPP/256._RPP  ! stencil 0
            c(2,1)=  105._RPP/128._RPP; c(1,1)=  105._RPP/256._RPP; c(0,1)=   -7._RPP/256._RPP  ! stencil 1
            c(2,2)=   75._RPP/128._RPP; c(1,2)=  -25._RPP/256._RPP; c(0,2)=    3._RPP/256._RPP  ! stencil 2
            c(2,3)=  -35._RPP/128._RPP; c(1,3)=   21._RPP/256._RPP; c(0,3)=   -3._RPP/256._RPP  ! stencil 3
            c(2,4)=   63._RPP/128._RPP; c(1,4)=  -45._RPP/256._RPP; c(0,4)=    7._RPP/256._RPP  ! stencil 4
            c(2,5)= -495._RPP/128._RPP; c(1,5)=  385._RPP/256._RPP; c(0,5)=  -63._RPP/256._RPP  ! stencil 5
          case(7) ! 13th order
            !  cell  6                ;    cell  5                ;    cell  4
            c(6,0)=  -21._RPP/1024._RPP; c(5,0)=   77._RPP/512._RPP ; c(4,0)= -495._RPP/1024._RPP  ! stencil 0
            c(6,1)=    7._RPP/1024._RPP; c(5,1)=  -27._RPP/512._RPP ; c(4,1)=  189._RPP/1024._RPP  ! stencil 1
            c(6,2)=   -5._RPP/1024._RPP; c(5,2)=   21._RPP/512._RPP ; c(4,2)= -175._RPP/1024._RPP  ! stencil 2
            c(6,3)=    7._RPP/1024._RPP; c(5,3)=  -35._RPP/512._RPP ; c(4,3)=  525._RPP/1024._RPP  ! stencil 3
            c(6,4)=  -21._RPP/1024._RPP; c(5,4)=  189._RPP/512._RPP ; c(4,4)=  945._RPP/1024._RPP  ! stencil 4
            c(6,5)=  231._RPP/1024._RPP; c(5,5)=  693._RPP/512._RPP ; c(4,5)=-1155._RPP/1024._RPP  ! stencil 5
            c(6,6)= 3003._RPP/1024._RPP; c(5,6)=-3003._RPP/512._RPP ; c(4,6)= 9009._RPP/1024._RPP  ! stencil 6
            !  cell  3                 ;   cell  2                ;    cell  1
            c(3,0)=  231._RPP/256._RPP ; c(2,0)=-1155._RPP/1024._RPP; c(1,0)=  693._RPP/512._RPP   ! stencil 0
            c(3,1)= -105._RPP/256._RPP ; c(2,1)=  945._RPP/1024._RPP; c(1,1)=  189._RPP/512._RPP   ! stencil 1
            c(3,2)=  175._RPP/256._RPP ; c(2,2)=  525._RPP/1024._RPP; c(1,2)=  -35._RPP/512._RPP   ! stencil 2
            c(3,3)=  175._RPP/256._RPP ; c(2,3)= -175._RPP/1024._RPP; c(1,3)=   21._RPP/512._RPP   ! stencil 3
            c(3,4)= -105._RPP/256._RPP ; c(2,4)=  189._RPP/1024._RPP; c(1,4)=  -27._RPP/512._RPP   ! stencil 4
            c(3,5)=  231._RPP/256._RPP ; c(2,5)= -495._RPP/1024._RPP; c(1,5)=   77._RPP/512._RPP   ! stencil 5
            c(3,6)=-2145._RPP/256._RPP ; c(2,6)= 5005._RPP/1024._RPP; c(1,6)= -819._RPP/512._RPP   ! stencil 6
            !  cell  0
            c(0,0)=  231._RPP/1024._RPP  ! stencil 0
            c(0,1)=  -21._RPP/1024._RPP  ! stencil 1
            c(0,2)=    7._RPP/1024._RPP  ! stencil 2
            c(0,3)=   -5._RPP/1024._RPP  ! stencil 3
            c(0,4)=    7._RPP/1024._RPP  ! stencil 4
            c(0,5)=  -21._RPP/1024._RPP  ! stencil 5
            c(0,6)=  231._RPP/1024._RPP  ! stencil 6
          case(8) ! 15th order
            !  cell  7                ;    cell  6                ;    cell  5
            c(7,0)=    33._RPP/2048._RPP; c(6,0)=  -273._RPP/2048._RPP; c(5,0)=  1001._RPP/2048._RPP  ! stencil 0
            c(7,1)=    -9._RPP/2048._RPP; c(6,1)=    77._RPP/2048._RPP; c(5,1)=  -297._RPP/2048._RPP  ! stencil 1
            c(7,2)=     5._RPP/2048._RPP; c(6,2)=   -45._RPP/2048._RPP; c(5,2)=   189._RPP/2048._RPP  ! stencil 2
            c(7,3)=    -5._RPP/2048._RPP; c(6,3)=    49._RPP/2048._RPP; c(5,3)=  -245._RPP/2048._RPP  ! stencil 3
            c(7,4)=     9._RPP/2048._RPP; c(6,4)=  -105._RPP/2048._RPP; c(5,4)=   945._RPP/2048._RPP  ! stencil 4
            c(7,5)=   -33._RPP/2048._RPP; c(6,5)=   693._RPP/2048._RPP; c(5,5)=  2079._RPP/2048._RPP  ! stencil 5
            c(7,6)=   429._RPP/2048._RPP; c(6,6)=  3003._RPP/2048._RPP; c(5,6)= -3003._RPP/2048._RPP  ! stencil 6
            c(7,7)=  6435._RPP/2048._RPP; c(6,7)=-15015._RPP/2048._RPP; c(5,7)= 27027._RPP/2048._RPP  ! stencil 7
            !  cell  4                ;    cell  3                ;    cell  2
            c(4,0)= -2145._RPP/2048._RPP; c(3,0)=  3003._RPP/2048._RPP; c(2,0)= -3003._RPP/2048._RPP  ! stencil 0
            c(4,1)=   693._RPP/2048._RPP; c(3,1)= -1155._RPP/2048._RPP; c(2,1)=  2079._RPP/2048._RPP  ! stencil 1
            c(4,2)=  -525._RPP/2048._RPP; c(3,2)=  1575._RPP/2048._RPP; c(2,2)=   945._RPP/2048._RPP  ! stencil 2
            c(4,3)=  1225._RPP/2048._RPP; c(3,3)=  1225._RPP/2048._RPP; c(2,3)=  -245._RPP/2048._RPP  ! stencil 3
            c(4,4)=  1575._RPP/2048._RPP; c(3,4)=  -525._RPP/2048._RPP; c(2,4)=   189._RPP/2048._RPP  ! stencil 4
            c(4,5)= -1155._RPP/2048._RPP; c(3,5)=   693._RPP/2048._RPP; c(2,5)=  -297._RPP/2048._RPP  ! stencil 5
            c(4,6)=  3003._RPP/2048._RPP; c(3,6)= -2145._RPP/2048._RPP; c(2,6)=  1001._RPP/2048._RPP  ! stencil 6
            c(4,7)=-32175._RPP/2048._RPP; c(3,7)= 25025._RPP/2048._RPP; c(2,7)=-12285._RPP/2048._RPP  ! stencil 7
            !  cell  1                ;    cell  0
            c(1,0)=  3003._RPP/2048._RPP; c(0,0)=   429._RPP/2048._RPP  ! stencil 0
            c(1,1)=   693._RPP/2048._RPP; c(0,1)=   -33._RPP/2048._RPP  ! stencil 1
            c(1,2)=  -105._RPP/2048._RPP; c(0,2)=     9._RPP/2048._RPP  ! stencil 2
            c(1,3)=    49._RPP/2048._RPP; c(0,3)=    -5._RPP/2048._RPP  ! stencil 3
            c(1,4)=   -45._RPP/2048._RPP; c(0,4)=     5._RPP/2048._RPP  ! stencil 4
            c(1,5)=    77._RPP/2048._RPP; c(0,5)=    -9._RPP/2048._RPP  ! stencil 5
            c(1,6)=  -273._RPP/2048._RPP; c(0,6)=    33._RPP/2048._RPP  ! stencil 6
            c(1,7)=  3465._RPP/2048._RPP; c(0,7)=  -429._RPP/2048._RPP  ! stencil 7
          case(9) ! 17th order
            !  cell  8                    ;     cell  7                   ;     cell  6
            c(8,0)=   -429._RPP/32768._RPP; c(7,0)=    495._RPP/4096._RPP ; c(6,0)=  -4095._RPP/8192._RPP   ! stencil 0
            c(8,1)=     99._RPP/32768._RPP; c(7,1)=   -117._RPP/4096._RPP ; c(6,1)=   1001._RPP/8192._RPP   ! stencil 1
            c(8,2)=    -45._RPP/32768._RPP; c(7,2)=     55._RPP/4096._RPP ; c(6,2)=   -495._RPP/8192._RPP   ! stencil 2
            c(8,3)=     35._RPP/32768._RPP; c(7,3)=    -45._RPP/4096._RPP ; c(6,3)=    441._RPP/8192._RPP   ! stencil 3
            c(8,4)=    -45._RPP/32768._RPP; c(7,4)=     63._RPP/4096._RPP ; c(6,4)=   -735._RPP/8192._RPP   ! stencil 4
            c(8,5)=     99._RPP/32768._RPP; c(7,5)=   -165._RPP/4096._RPP ; c(6,5)=   3465._RPP/8192._RPP   ! stencil 5
            c(8,6)=   -429._RPP/32768._RPP; c(7,6)=   1287._RPP/4096._RPP ; c(6,6)=   9009._RPP/8192._RPP   ! stencil 6
            c(8,7)=   6435._RPP/32768._RPP; c(7,7)=   6435._RPP/4096._RPP ; c(6,7)= -15015._RPP/8192._RPP   ! stencil 7
            c(8,8)= 109395._RPP/32768._RPP; c(7,8)= -36465._RPP/4096._RPP ; c(6,8)= 153153._RPP/8192._RPP   ! stencil 8
            !  cell  5                    ;    cell  4                    ;     cell  3
            c(5,0)=   5005._RPP/4096._RPP ; c(4,0)= -32175._RPP/16384._RPP; c(3,0)=   9009._RPP/4096._RPP   ! stencil 0
            c(5,1)=  -1287._RPP/4096._RPP ; c(4,1)=   9009._RPP/16384._RPP; c(3,1)=  -3003._RPP/4096._RPP   ! stencil 1
            c(5,2)=    693._RPP/4096._RPP ; c(4,2)=  -5775._RPP/16384._RPP; c(3,2)=   3465._RPP/4096._RPP   ! stencil 2
            c(5,3)=   -735._RPP/4096._RPP ; c(4,3)=  11025._RPP/16384._RPP; c(3,3)=   2205._RPP/4096._RPP   ! stencil 3
            c(5,4)=   2205._RPP/4096._RPP ; c(4,4)=  11025._RPP/16384._RPP; c(3,4)=   -735._RPP/4096._RPP   ! stencil 4
            c(5,5)=   3465._RPP/4096._RPP ; c(4,5)=  -5775._RPP/16384._RPP; c(3,5)=    693._RPP/4096._RPP   ! stencil 5
            c(5,6)=  -3003._RPP/4096._RPP ; c(4,6)=   9009._RPP/16384._RPP; c(3,6)=  -1287._RPP/4096._RPP   ! stencil 6
            c(5,7)=   9009._RPP/4096._RPP ; c(4,7)= -32175._RPP/16384._RPP; c(3,7)=   5005._RPP/4096._RPP   ! stencil 7
            c(5,8)=-109395._RPP/4096._RPP ; c(4,8)= 425425._RPP/16384._RPP; c(3,8)= -69615._RPP/4096._RPP   ! stencil 8
            !   cell  2                   ;     cell  1                   ;     cell  0
            c(2,0)= -15015._RPP/8192._RPP ; c(1,0)=   6435._RPP/4096._RPP ; c(0,0)=   6435._RPP/32768._RPP  ! stencil 0
            c(2,1)=   9009._RPP/8192._RPP ; c(1,1)=   1287._RPP/4096._RPP ; c(0,1)=   -429._RPP/32768._RPP  ! stencil 1
            c(2,2)=   3465._RPP/8192._RPP ; c(1,2)=   -165._RPP/4096._RPP ; c(0,2)=     99._RPP/32768._RPP  ! stencil 2
            c(2,3)=   -735._RPP/8192._RPP ; c(1,3)=     63._RPP/4096._RPP ; c(0,3)=    -45._RPP/32768._RPP  ! stencil 3
            c(2,4)=    441._RPP/8192._RPP ; c(1,4)=    -45._RPP/4096._RPP ; c(0,4)=     35._RPP/32768._RPP  ! stencil 4
            c(2,5)=   -495._RPP/8192._RPP ; c(1,5)=     55._RPP/4096._RPP ; c(0,5)=    -45._RPP/32768._RPP  ! stencil 5
            c(2,6)=   1001._RPP/8192._RPP ; c(1,6)=   -117._RPP/4096._RPP ; c(0,6)=     99._RPP/32768._RPP  ! stencil 6
            c(2,7)=  -4095._RPP/8192._RPP ; c(1,7)=    495._RPP/4096._RPP ; c(0,7)=   -429._RPP/32768._RPP  ! stencil 7
            c(2,8)=  58905._RPP/8192._RPP ; c(1,8)=  -7293._RPP/4096._RPP ; c(0,8)=   6435._RPP/32768._RPP  ! stencil 8
        endselect
      elseif(x_target==0.5_RPP) then
        ! right interface (i+1/2)
        select case(self%S)
          case(2) ! 3rd order
            !  cell  1      ;    cell  0
            c(1,0)= -0.5_RPP; c(0,0)=  1.5_RPP ! stencil 0
            c(1,1)=  0.5_RPP; c(0,1)=  0.5_RPP ! stencil 1
          case(3) ! 5th order
            !  cell  2            ;    cell  1            ;    cell  0
            c(2,0)=  3._RPP/8._RPP; c(1,0)= -5._RPP/4._RPP; c(0,0)= 15._RPP/8._RPP ! stencil 0
            c(2,1)= -1._RPP/8._RPP; c(1,1)=  3._RPP/4._RPP; c(0,1)=  3._RPP/8._RPP ! stencil 1
            c(2,2)=  3._RPP/8._RPP; c(1,2)=  3._RPP/4._RPP; c(0,2)= -1._RPP/8._RPP ! stencil 2
          case(4) ! 7th order
            !  cell  3             ;    cell  2             ;   cell  1              ;    cell  0
            c(3,0)= -5._RPP/16._RPP; c(2,0)= 21._RPP/16._RPP; c(1,0)=-35._RPP/16._RPP; c(0,0)= 35._RPP/16._RPP ! stencil 0
            c(3,1)=  1._RPP/16._RPP; c(2,1)= -5._RPP/16._RPP; c(1,1)= 15._RPP/16._RPP; c(0,1)=  5._RPP/16._RPP ! stencil 1
            c(3,2)= -1._RPP/16._RPP; c(2,2)=  9._RPP/16._RPP; c(1,2)=  9._RPP/16._RPP; c(0,2)= -1._RPP/16._RPP ! stencil 2
            c(3,3)=  5._RPP/16._RPP; c(2,3)= 15._RPP/16._RPP; c(1,3)= -5._RPP/16._RPP; c(0,3)=  1._RPP/16._RPP ! stencil 3
          case(5) ! 9th order
            !  cell  4               ;    cell  3               ;   cell  2                ;    cell  1
            c(4,0)=  35._RPP/128._RPP; c(3,0)= -45._RPP/32._RPP ; c(2,0)= 189._RPP/64._RPP ; c(1,0)=-105._RPP/32._RPP  ! stencil 0
            c(4,1)=  -5._RPP/128._RPP; c(3,1)=   7._RPP/32._RPP ; c(2,1)= -35._RPP/64._RPP ; c(1,1)=  35._RPP/32._RPP  ! stencil 1
            c(4,2)=   3._RPP/128._RPP; c(3,2)=  -5._RPP/32._RPP ; c(2,2)=  45._RPP/64._RPP ; c(1,2)=  15._RPP/32._RPP  ! stencil 2
            c(4,3)=  -5._RPP/128._RPP; c(3,3)=  15._RPP/32._RPP ; c(2,3)=  45._RPP/64._RPP ; c(1,3)=  -5._RPP/32._RPP  ! stencil 3
            c(4,4)=  35._RPP/128._RPP; c(3,4)=  35._RPP/32._RPP ; c(2,4)= -35._RPP/64._RPP ; c(1,4)=   7._RPP/32._RPP  ! stencil 4
            !  cell  0
            c(0,0)= 315._RPP/128._RPP ! stencil 0
            c(0,1)=  35._RPP/128._RPP ! stencil 1
            c(0,2)=  -5._RPP/128._RPP ! stencil 2
            c(0,3)=   3._RPP/128._RPP ! stencil 3
            c(0,4)=  -5._RPP/128._RPP ! stencil 4
          case(6) ! 11th order
            !  cell  5                ;    cell  4                ;   cell  3
            c(5,0)=  -63._RPP/256._RPP; c(4,0)=  385._RPP/256._RPP; c(3,0)= -495._RPP/128._RPP  ! stencil 0
            c(5,1)=    7._RPP/256._RPP; c(4,1)=  -45._RPP/256._RPP; c(3,1)=   63._RPP/128._RPP  ! stencil 1
            c(5,2)=   -3._RPP/256._RPP; c(4,2)=   21._RPP/256._RPP; c(3,2)=  -35._RPP/128._RPP  ! stencil 2
            c(5,3)=    3._RPP/256._RPP; c(4,3)=  -25._RPP/256._RPP; c(3,3)=   75._RPP/128._RPP  ! stencil 3
            c(5,4)=   -7._RPP/256._RPP; c(4,4)=  105._RPP/256._RPP; c(3,4)=  105._RPP/128._RPP  ! stencil 4
            c(5,5)=   63._RPP/256._RPP; c(4,5)=  315._RPP/256._RPP; c(3,5)= -105._RPP/128._RPP  ! stencil 5
            !  cell  2                ;    cell  1                ;   cell  0
            c(2,0)=  693._RPP/128._RPP; c(1,0)=-1155._RPP/256._RPP; c(0,0)=  693._RPP/256._RPP  ! stencil 0
            c(2,1)= -105._RPP/128._RPP; c(1,1)=  315._RPP/256._RPP; c(0,1)=   63._RPP/256._RPP  ! stencil 1
            c(2,2)=  105._RPP/128._RPP; c(1,2)=  105._RPP/256._RPP; c(0,2)=   -7._RPP/256._RPP  ! stencil 2
            c(2,3)=   75._RPP/128._RPP; c(1,3)=  -25._RPP/256._RPP; c(0,3)=    3._RPP/256._RPP  ! stencil 3
            c(2,4)=  -35._RPP/128._RPP; c(1,4)=   21._RPP/256._RPP; c(0,4)=   -3._RPP/256._RPP  ! stencil 4
            c(2,5)=   63._RPP/128._RPP; c(1,5)=  -45._RPP/256._RPP; c(0,5)=    7._RPP/256._RPP  ! stencil 5
          case(7) ! 13th order
            !  cell  6                ;    cell  5                ;    cell  4
            c(6,0)=  231._RPP/1024._RPP; c(5,0)= -819._RPP/512._RPP ; c(4,0)= 5005._RPP/1024._RPP  ! stencil 0
            c(6,1)=  -21._RPP/1024._RPP; c(5,1)=   77._RPP/512._RPP ; c(4,1)= -495._RPP/1024._RPP  ! stencil 1
            c(6,2)=    7._RPP/1024._RPP; c(5,2)=  -27._RPP/512._RPP ; c(4,2)=  189._RPP/1024._RPP  ! stencil 2
            c(6,3)=   -5._RPP/1024._RPP; c(5,3)=   21._RPP/512._RPP ; c(4,3)= -175._RPP/1024._RPP  ! stencil 3
            c(6,4)=    7._RPP/1024._RPP; c(5,4)=  -35._RPP/512._RPP ; c(4,4)=  525._RPP/1024._RPP  ! stencil 4
            c(6,5)=  -21._RPP/1024._RPP; c(5,5)=  189._RPP/512._RPP ; c(4,5)=  945._RPP/1024._RPP  ! stencil 5
            c(6,6)=  231._RPP/1024._RPP; c(5,6)=  693._RPP/512._RPP ; c(4,6)=-1155._RPP/1024._RPP  ! stencil 6
            !  cell  3                ;    cell  2                ;    cell  1
            c(3,0)=-2145._RPP/256._RPP ; c(2,0)= 9009._RPP/1024._RPP; c(1,0)=-3003._RPP/512._RPP   ! stencil 0
            c(3,1)=  231._RPP/256._RPP ; c(2,1)=-1155._RPP/1024._RPP; c(1,1)=  693._RPP/512._RPP   ! stencil 1
            c(3,2)= -105._RPP/256._RPP ; c(2,2)=  945._RPP/1024._RPP; c(1,2)=  189._RPP/512._RPP   ! stencil 2
            c(3,3)=  175._RPP/256._RPP ; c(2,3)=  525._RPP/1024._RPP; c(1,3)=  -35._RPP/512._RPP   ! stencil 3
            c(3,4)=  175._RPP/256._RPP ; c(2,4)= -175._RPP/1024._RPP; c(1,4)=   21._RPP/512._RPP   ! stencil 4
            c(3,5)= -105._RPP/256._RPP ; c(2,5)=  189._RPP/1024._RPP; c(1,5)=  -27._RPP/512._RPP   ! stencil 5
            c(3,6)=  231._RPP/256._RPP ; c(2,6)= -495._RPP/1024._RPP; c(1,6)=   77._RPP/512._RPP   ! stencil 6
            !  cell  0
            c(0,0)= 3003._RPP/1024._RPP  ! stencil 0
            c(0,1)=  231._RPP/1024._RPP  ! stencil 1
            c(0,2)=  -21._RPP/1024._RPP  ! stencil 2
            c(0,3)=    7._RPP/1024._RPP  ! stencil 3
            c(0,4)=   -5._RPP/1024._RPP  ! stencil 4
            c(0,5)=    7._RPP/1024._RPP  ! stencil 5
            c(0,6)=  -21._RPP/1024._RPP  ! stencil 6
          case(8) ! 15th order
            !  cell  7                ;    cell  6                ;    cell  5
            c(7,0)=  -429._RPP/2048._RPP; c(6,0)=  3465._RPP/2048._RPP; c(5,0)=-12285._RPP/2048._RPP  ! stencil 0
            c(7,1)=    33._RPP/2048._RPP; c(6,1)=  -273._RPP/2048._RPP; c(5,1)=  1001._RPP/2048._RPP  ! stencil 1
            c(7,2)=    -9._RPP/2048._RPP; c(6,2)=    77._RPP/2048._RPP; c(5,2)=  -297._RPP/2048._RPP  ! stencil 2
            c(7,3)=     5._RPP/2048._RPP; c(6,3)=   -45._RPP/2048._RPP; c(5,3)=   189._RPP/2048._RPP  ! stencil 3
            c(7,4)=    -5._RPP/2048._RPP; c(6,4)=    49._RPP/2048._RPP; c(5,4)=  -245._RPP/2048._RPP  ! stencil 4
            c(7,5)=     9._RPP/2048._RPP; c(6,5)=  -105._RPP/2048._RPP; c(5,5)=   945._RPP/2048._RPP  ! stencil 5
            c(7,6)=   -33._RPP/2048._RPP; c(6,6)=   693._RPP/2048._RPP; c(5,6)=  2079._RPP/2048._RPP  ! stencil 6
            c(7,7)=   429._RPP/2048._RPP; c(6,7)=  3003._RPP/2048._RPP; c(5,7)= -3003._RPP/2048._RPP  ! stencil 7
            !  cell  4                ;    cell  3                ;    cell  2
            c(4,0)= 25025._RPP/2048._RPP; c(3,0)=-32175._RPP/2048._RPP; c(2,0)= 27027._RPP/2048._RPP  ! stencil 0
            c(4,1)= -2145._RPP/2048._RPP; c(3,1)=  3003._RPP/2048._RPP; c(2,1)= -3003._RPP/2048._RPP  ! stencil 1
            c(4,2)=   693._RPP/2048._RPP; c(3,2)= -1155._RPP/2048._RPP; c(2,2)=  2079._RPP/2048._RPP  ! stencil 2
            c(4,3)=  -525._RPP/2048._RPP; c(3,3)=  1575._RPP/2048._RPP; c(2,3)=   945._RPP/2048._RPP  ! stencil 3
            c(4,4)=  1225._RPP/2048._RPP; c(3,4)=  1225._RPP/2048._RPP; c(2,4)=  -245._RPP/2048._RPP  ! stencil 4
            c(4,5)=  1575._RPP/2048._RPP; c(3,5)=  -525._RPP/2048._RPP; c(2,5)=   189._RPP/2048._RPP  ! stencil 5
            c(4,6)= -1155._RPP/2048._RPP; c(3,6)=   693._RPP/2048._RPP; c(2,6)=  -297._RPP/2048._RPP  ! stencil 6
            c(4,7)=  3003._RPP/2048._RPP; c(3,7)= -2145._RPP/2048._RPP; c(2,7)=  1001._RPP/2048._RPP  ! stencil 7
            !  cell  1                  ;    cell  0
            c(1,0)=-15015._RPP/2048._RPP; c(0,0)=  6435._RPP/2048._RPP  ! stencil 0
            c(1,1)=  3003._RPP/2048._RPP; c(0,1)=   429._RPP/2048._RPP  ! stencil 1
            c(1,2)=   693._RPP/2048._RPP; c(0,2)=   -33._RPP/2048._RPP  ! stencil 2
            c(1,3)=  -105._RPP/2048._RPP; c(0,3)=     9._RPP/2048._RPP  ! stencil 3
            c(1,4)=    49._RPP/2048._RPP; c(0,4)=    -5._RPP/2048._RPP  ! stencil 4
            c(1,5)=   -45._RPP/2048._RPP; c(0,5)=     5._RPP/2048._RPP  ! stencil 5
            c(1,6)=    77._RPP/2048._RPP; c(0,6)=    -9._RPP/2048._RPP  ! stencil 6
            c(1,7)=  -273._RPP/2048._RPP; c(0,7)=    33._RPP/2048._RPP  ! stencil 7
          case(9) ! 17th order
            !  cell  8                  ;     cell  7                 ;     cell  6
            c(8,0)=   6435._RPP/32768._RPP; c(7,0)=  -7293._RPP/ 4096._RPP; c(6,0)=  58905._RPP/ 8192._RPP  ! stencil 0
            c(8,1)=   -429._RPP/32768._RPP; c(7,1)=    495._RPP/ 4096._RPP; c(6,1)=  -4095._RPP/ 8192._RPP  ! stencil 1
            c(8,2)=     99._RPP/32768._RPP; c(7,2)=   -117._RPP/ 4096._RPP; c(6,2)=   1001._RPP/ 8192._RPP  ! stencil 2
            c(8,3)=    -45._RPP/32768._RPP; c(7,3)=     55._RPP/ 4096._RPP; c(6,3)=   -495._RPP/ 8192._RPP  ! stencil 3
            c(8,4)=     35._RPP/32768._RPP; c(7,4)=    -45._RPP/ 4096._RPP; c(6,4)=    441._RPP/ 8192._RPP  ! stencil 4
            c(8,5)=    -45._RPP/32768._RPP; c(7,5)=     63._RPP/ 4096._RPP; c(6,5)=   -735._RPP/ 8192._RPP  ! stencil 5
            c(8,6)=     99._RPP/32768._RPP; c(7,6)=   -165._RPP/ 4096._RPP; c(6,6)=   3465._RPP/ 8192._RPP  ! stencil 6
            c(8,7)=   -429._RPP/32768._RPP; c(7,7)=   1287._RPP/ 4096._RPP; c(6,7)=   9009._RPP/ 8192._RPP  ! stencil 7
            c(8,8)=   6435._RPP/32768._RPP; c(7,8)=   6435._RPP/ 4096._RPP; c(6,8)= -15015._RPP/ 8192._RPP  ! stencil 8
            !   cell  5                 ; !  cell  4                  ;     cell  3
            c(5,0)= -69615._RPP/ 4096._RPP; c(4,0)= 425425._RPP/16384._RPP; c(3,0)=-109395._RPP/ 4096._RPP  ! stencil 0
            c(5,1)=   5005._RPP/ 4096._RPP; c(4,1)= -32175._RPP/16384._RPP; c(3,1)=   9009._RPP/ 4096._RPP  ! stencil 1
            c(5,2)=  -1287._RPP/ 4096._RPP; c(4,2)=   9009._RPP/16384._RPP; c(3,2)=  -3003._RPP/ 4096._RPP  ! stencil 2
            c(5,3)=    693._RPP/ 4096._RPP; c(4,3)=  -5775._RPP/16384._RPP; c(3,3)=   3465._RPP/ 4096._RPP  ! stencil 3
            c(5,4)=   -735._RPP/ 4096._RPP; c(4,4)=  11025._RPP/16384._RPP; c(3,4)=   2205._RPP/ 4096._RPP  ! stencil 4
            c(5,5)=   2205._RPP/ 4096._RPP; c(4,5)=  11025._RPP/16384._RPP; c(3,5)=   -735._RPP/ 4096._RPP  ! stencil 5
            c(5,6)=   3465._RPP/ 4096._RPP; c(4,6)=  -5775._RPP/16384._RPP; c(3,6)=    693._RPP/ 4096._RPP  ! stencil 6
            c(5,7)=  -3003._RPP/ 4096._RPP; c(4,7)=   9009._RPP/16384._RPP; c(3,7)=  -1287._RPP/ 4096._RPP  ! stencil 7
            c(5,8)=   9009._RPP/ 4096._RPP; c(4,8)= -32175._RPP/16384._RPP; c(3,8)=   5005._RPP/ 4096._RPP  ! stencil 8
            !   cell  2                 ;     cell  1                 ;     cell  0
            c(2,0)= 153153._RPP/ 8192._RPP; c(1,0)= -36465._RPP/ 4096._RPP; c(0,0)= 109395._RPP/32768._RPP  ! stencil 0
            c(2,1)= -15015._RPP/ 8192._RPP; c(1,1)=   6435._RPP/ 4096._RPP; c(0,1)=   6435._RPP/32768._RPP  ! stencil 1
            c(2,2)=   9009._RPP/ 8192._RPP; c(1,2)=   1287._RPP/ 4096._RPP; c(0,2)=   -429._RPP/32768._RPP  ! stencil 2
            c(2,3)=   3465._RPP/ 8192._RPP; c(1,3)=   -165._RPP/ 4096._RPP; c(0,3)=     99._RPP/32768._RPP  ! stencil 3
            c(2,4)=   -735._RPP/ 8192._RPP; c(1,4)=     63._RPP/ 4096._RPP; c(0,4)=    -45._RPP/32768._RPP  ! stencil 4
            c(2,5)=    441._RPP/ 8192._RPP; c(1,5)=    -45._RPP/ 4096._RPP; c(0,5)=     35._RPP/32768._RPP  ! stencil 5
            c(2,6)=   -495._RPP/ 8192._RPP; c(1,6)=     55._RPP/ 4096._RPP; c(0,6)=    -45._RPP/32768._RPP  ! stencil 6
            c(2,7)=   1001._RPP/ 8192._RPP; c(1,7)=   -117._RPP/ 4096._RPP; c(0,7)=     99._RPP/32768._RPP  ! stencil 7
            c(2,8)=  -4095._RPP/ 8192._RPP; c(1,8)=    495._RPP/ 4096._RPP; c(0,8)=   -429._RPP/32768._RPP  ! stencil 8
        endselect
      else
        ! internal point
        allocate(f(0:S-1, 0:S-1))
        do k=0,S-1  !stencils loop
          c_sum = 0._RPP
          do j=0,S-2  !values loop
            prod = 1._RPP
            do i=0,S-1
              if (i==j) cycle
              prod = prod * ((x_target - stencil(-S+k+i+1)) / (stencil(-S+k+j+1) - stencil(-S+k+i+1)))
            enddo
            f(j,k) = prod
            c_sum = c_sum + prod
          enddo
          f(S-1,k) = 1._RPP - c_sum
        enddo
        do k=0,S-1
          do j=0,S-1
            c(j,k) = f(S-1-j,k)
          enddo
        enddo
      endif
    endassociate
  endselect
  endsubroutine create

  pure subroutine compute_int(self, stencil, values)
  !< Compute interpolations (interpolation).
  class(interpolations_int_js), intent(in)  :: self               !< Interpolations.
  real(RPP),                    intent(in)  :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(RPP),                    intent(out) :: values(0:)         !< Interpolations values.
  integer(I_P)                              :: s1                 !< Counter.
  integer(I_P)                              :: s2                 !< Counter.

  values = 0._RPP
  do s1=0, self%S - 1 ! stencils loop
    do s2=0, self%S - 1 ! values loop
      values(s1) = values(s1) + self%coef(s2, s1) * stencil(-s2 + s1)
    enddo
  enddo
  endsubroutine compute_int

  pure subroutine compute_rec(self, stencil, values)
  !< Compute interpolations (reconstruct).
  class(interpolations_int_js), intent(in)  :: self                  !< Interpolations.
  real(RPP),                    intent(in)  :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(RPP),                    intent(out) :: values(1:, 0:)        !< Interpolations values.
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
  string = prefix_//string//'  - S   = '//trim(str(self%S))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy interpolations.
  class(interpolations_int_js), intent(inout) :: self !< Interpolations.

  call self%destroy_
  if (allocated(self%coef)) deallocate(self%coef)
  endsubroutine destroy
endmodule wenoof_interpolations_int_js
