!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for interpolation.
module wenoof_kappa_int_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for interpolation.
!<
!< @note The provided WENO kappa implements the linear weights defined in *High Order Weighted Essentially
!< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
!< doi:10.1137/070679065.

use penf, only : I_P, R_P, str
use wenoof_base_object, only : base_object, base_object_constructor
use wenoof_interpolations_factory, only : interpolations_factory
use wenoof_interpolations_object, only : interpolations_object, interpolations_object_constructor
use wenoof_interpolations_int_js, only : interpolations_int_js
use wenoof_kappa_object, only : kappa_object, kappa_object_constructor

implicit none
private
public :: kappa_int_js
public :: kappa_int_js_constructor

type, extends(kappa_object_constructor) :: kappa_int_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal kappa object constructor.
  class(interpolations_object_constructor), allocatable :: interpolations_constructor !< Interpolations coefficients constructor.
  real(R_P), allocatable                                :: stencil(:)                 !< Stencil used for interpolation, [1-S:S-1].
  real(R_P)                                             :: x_target                   !< Coordinate of the interpolation point.
  contains
    ! public deferred methods
    procedure, pass(lhs) :: constr_assign_constr !< `=` operator.
endtype kappa_int_js_constructor

type, extends(kappa_object):: kappa_int_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet kappa object.
  !<
  !< @note The provided WENO kappa implements the linear weights defined in *High Order Weighted Essentially
  !< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
  !< doi:10.1137/070679065.
  class(interpolations_object), allocatable :: interpolations !< Interpolations object for getting coefficients.
  real(R_P),                    allocatable :: values(:)      !< Kappa coefficients values [0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: create               !< Create kappa.
    procedure, pass(self) :: compute_int          !< Compute kappa (interpolate).
    procedure, pass(self) :: compute_rec          !< Compute kappa (reconstruct).
    procedure, pass(self) :: description          !< Return object string-description.
    procedure, pass(self) :: destroy              !< Destroy kappa.
    procedure, pass(lhs)  :: object_assign_object !< `=` operator.
endtype kappa_int_js

contains
  ! constructor

  ! deferred public methods
  subroutine constr_assign_constr(lhs, rhs)
  !< `=` operator.
  class(kappa_int_js_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(kappa_int_js_constructor)
     if (allocated(rhs%interpolations_constructor)) then
        if (.not.allocated(lhs%interpolations_constructor)) &
           allocate(lhs%interpolations_constructor, mold=rhs%interpolations_constructor)
           lhs%interpolations_constructor = rhs%interpolations_constructor
     else
        if (allocated(lhs%interpolations_constructor)) deallocate(lhs%interpolations_constructor)
     endif
     if (allocated(rhs%stencil)) then
           lhs%stencil = rhs%stencil
     else
        if (allocated(lhs%stencil)) deallocate(lhs%stencil)
     endif
     lhs%x_target = rhs%x_target
  endselect
  endsubroutine constr_assign_constr

  ! deferred public methods
  subroutine create(self, constructor)
  !< Create kappa.
  !<
  !< @note The kappa coefficients are also computed, they being constants.
  class(kappa_int_js),            intent(inout) :: self        !< Kappa.
  class(base_object_constructor), intent(in)    :: constructor !< Kappa constructor.
  type(interpolations_factory)                  :: i_factory   !< Interpolations factory.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(0:self%S - 1))
  select type(constructor)
  type is(kappa_int_js_constructor)
    call i_factory%create(constructor=constructor%interpolations_constructor, object=self%interpolations)
    call self%compute_int(stencil=constructor%stencil, x_target=constructor%x_target, values=self%values)
  endselect
  endsubroutine create

  pure subroutine compute_int(self, stencil, x_target, values)
  !< Compute kappa.
  class(kappa_int_js), intent(in)  :: self                        !< Kappa.
  real(R_P),           intent(in)  :: stencil(1-self%S:)          !< Stencil used for interpolation, [1-S:S-1].
  real(R_P),           intent(in)  :: x_target                    !< Coordinate of the interpolation point.
  real(R_P),           intent(out) :: values(0:)                  !< Kappa values.
  real(R_P)                        :: coeff(0:2*self%S-2)         !< Interpolation coefficients on the whole stencil.
  real(R_P)                        :: coef(0:self%S-1,0:self%S-1) !< Temporary variable.
  real(R_P)                        :: prod                        !< Temporary variable.
  real(R_P)                        :: coeff_t                     !< Temporary variable.
  real(R_P)                        :: val_sum                     !< Temporary variable.
  integer(I_P)                     :: i, j, k                     !< Counters.

  associate(S=>self%S, interp=>self%interpolations)
    if (x_target == -0.5_R_P) then ! left interface (i-1/2)
      select case(S)
        case(2) ! 3rd order
          values(0) = 3._R_P/4._R_P ! stencil 0
          values(1) = 1._R_P/4._R_P ! stencil 1
        case(3) ! 5th order
          values(0) = 5._R_P/16._R_P ! stencil 0
          values(1) = 5._R_P/8._R_P  ! stencil 1
          values(2) = 1._R_P/16._R_P ! stencil 2
        case(4) ! 7th order
          values(0) =  7._R_P/64._R_P ! stencil 0
          values(1) = 35._R_P/64._R_P ! stencil 1
          values(2) = 21._R_P/64._R_P ! stencil 2
          values(3) =  1._R_P/64._R_P ! stencil 3
        case(5) ! 9th order
          values(0) =  9._R_P/256._R_P ! stencil 0
          values(1) = 21._R_P/64._R_P  ! stencil 1
          values(2) = 63._R_P/128._R_P ! stencil 2
          values(3) =  9._R_P/64._R_P  ! stencil 3
          values(4) =  1._R_P/256._R_P ! stencil 4
        case(6) ! 11th order
          values(0) =  11._R_P/1024._R_P  ! stencil 0
          values(1) = 165._R_P/1024._R_P  ! stencil 1
          values(2) = 231._R_P/512._R_P   ! stencil 2
          values(3) = 165._R_P/512._R_P   ! stencil 3
          values(4) =  55._R_P/1024._R_P  ! stencil 4
          values(5) =   1._R_P/1024._R_P  ! stencil 5
        case(7) ! 13th order
          values(0) =   13._R_P/4096._R_P  ! stencil 0
          values(1) =  143._R_P/2048._R_P  ! stencil 1
          values(2) = 1287._R_P/4096._R_P  ! stencil 2
          values(3) =  429._R_P/1024._R_P  ! stencil 3
          values(4) =  179._R_P/1024._R_P  ! stencil 4
          values(5) =   39._R_P/2048._R_P  ! stencil 5
          values(6) =    1._R_P/4096._R_P  ! stencil 6
        case(8) ! 15th order
          values(0) =   15._R_P/16384._R_P ! stencil 0
          values(1) =  455._R_P/16384._R_P ! stencil 1
          values(2) = 3003._R_P/16384._R_P ! stencil 2
          values(3) = 6435._R_P/16384._R_P ! stencil 3
          values(4) = 5005._R_P/16384._R_P ! stencil 4
          values(5) = 1365._R_P/16384._R_P ! stencil 5
          values(6) =  105._R_P/16384._R_P ! stencil 6
          values(7) =    1._R_P/16384._R_P ! stencil 7
        case(9) ! 17th order
          values(0) =    17._R_P/65536._R_P  ! stencil 0
          values(1) =    85._R_P/8192._R_P   ! stencil 1
          values(2) =  1547._R_P/16384._R_P  ! stencil 2
          values(3) =  2431._R_P/8192._R_P   ! stencil 3
          values(4) = 12155._R_P/32768._R_P  ! stencil 4
          values(5) =  1547._R_P/8192._R_P   ! stencil 5
          values(6) =   595._R_P/16384._R_P  ! stencil 6
          values(7) =    17._R_P/8192._R_P   ! stencil 7
          values(8) =     1._R_P/65536._R_P  ! stencil 8
      endselect
    elseif(x_target == 0.5_R_P) then ! right interface (i+1/2)
      select case(S)
        case(2) ! 3rd order
          values(0) = 1._R_P/4._R_P ! stencil 0
          values(1) = 3._R_P/4._R_P ! stencil 1
        case(3) ! 5th order
          values(0) = 1._R_P/16._R_P ! stencil 0
          values(1) = 5._R_P/8._R_P  ! stencil 1
          values(2) = 5._R_P/16._R_P ! stencil 2
        case(4) ! 7th order
          values(0) =  1._R_P/64._R_P ! stencil 0
          values(1) = 21._R_P/64._R_P ! stencil 1
          values(2) = 35._R_P/64._R_P ! stencil 2
          values(3) =  7._R_P/64._R_P ! stencil 3
        case(5) ! 9th order
          values(0) =  1._R_P/256._R_P ! stencil 0
          values(1) =  9._R_P/64._R_P  ! stencil 1
          values(2) = 63._R_P/128._R_P  ! stencil 2
          values(3) = 21._R_P/64._R_P  ! stencil 3
          values(4) =  9._R_P/256._R_P ! stencil 4
        case(6) ! 11th order
          values(0) =   1._R_P/1024._R_P  ! stencil 0
          values(1) =  55._R_P/1024._R_P  ! stencil 1
          values(2) = 165._R_P/512._R_P   ! stencil 2
          values(3) = 231._R_P/512._R_P   ! stencil 3
          values(4) = 165._R_P/1024._R_P  ! stencil 4
          values(5) =  11._R_P/1024._R_P  ! stencil 5
        case(7) ! 13th order
          values(0) =    1._R_P/4096._R_P  ! stencil 0
          values(1) =   39._R_P/2048._R_P  ! stencil 1
          values(2) =  179._R_P/1024._R_P  ! stencil 2
          values(3) =  429._R_P/1024._R_P  ! stencil 3
          values(4) = 1287._R_P/4096._R_P  ! stencil 4
          values(5) =  143._R_P/2048._R_P  ! stencil 5
          values(6) =   13._R_P/4096._R_P  ! stencil 6
        case(8) ! 15th order
          values(0) =    1._R_P/16384._R_P ! stencil 0
          values(1) =  105._R_P/16384._R_P ! stencil 1
          values(2) = 1365._R_P/16384._R_P ! stencil 2
          values(3) = 5005._R_P/16384._R_P ! stencil 3
          values(4) = 6435._R_P/16384._R_P ! stencil 4
          values(5) = 3003._R_P/16384._R_P ! stencil 5
          values(6) =  455._R_P/16384._R_P ! stencil 6
          values(7) =   15._R_P/16384._R_P ! stencil 7
        case(9) ! 17th order
          values(0) =     1._R_P/65536._R_P  ! stencil 0
          values(1) =    17._R_P/8192._R_P   ! stencil 1
          values(2) =   595._R_P/16384._R_P  ! stencil 2
          values(3) =  1547._R_P/8192._R_P   ! stencil 3
          values(4) = 12155._R_P/32768._R_P  ! stencil 4
          values(5) =  2431._R_P/8192._R_P   ! stencil 5
          values(6) =  1547._R_P/16384._R_P  ! stencil 6
          values(7) =    85._R_P/8192._R_P   ! stencil 7
          values(8) =    17._R_P/65536._R_P  ! stencil 8
      endselect
    elseif(x_target == 0._R_P) then
      values = 1._R_P / S
    else
      ! internal point
      val_sum = 0._R_P
      do j=0, 2 * S - 3 ! values loop
        prod = 1._R_P
        do i=0, 2 * S - 2
          if (i==j) cycle
          prod = prod * ((x_target - stencil(-S+i+1)) / (stencil(-S+j+1) - stencil(-S+i+1)))
        enddo
        coeff(j) = prod
        val_sum = val_sum + coeff(j)
      enddo
      coeff(2*S-2) = 1._R_P - val_sum
      select type(interp)
        type is(interpolations_int_js)
          val_sum = 0._R_P
          do k=0,S-1
            do j=0,S-1
              coef(j, k) = interp%coef(S-1-j,S-1-k)
            enddo
          enddo
          do j = 0,S-2
            coeff_t = 0._R_P
            k = j
            do i = 0,j-1
              coeff_t = coeff_t + values(i) * coef(k,i)
              k = k - 1
            enddo
            values(j) = (coeff(j) - coeff_t) / coef(0,j)
            val_sum = val_sum + values(j)
          enddo
          values(S-1) = 1._R_P - val_sum
      endselect
    endif
  endassociate
  endsubroutine compute_int

  pure subroutine compute_rec(self, values)
  !< Compute kappa (reconstruct).
  class(kappa_int_js), intent(in)  :: self          !< Kappa.
  real(R_P),           intent(out) :: values(1:,0:) !< Kappa values.
  ! empty procedure
  endsubroutine compute_rec

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(kappa_int_js), intent(in)           :: self             !< Kappa coefficient.
  character(len=*),    intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable             :: string           !< String-description.
  character(len=:), allocatable             :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter               :: NL=new_line('a') !< New line char.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu kappa coefficients object for interpolation:'//NL
  string = string//prefix_//'  - S   = '//trim(str(self%S))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy kappa.
  class(kappa_int_js), intent(inout) :: self !< Kappa.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy

  pure subroutine object_assign_object(lhs, rhs)
  !< `=` operator.
  class(kappa_int_js), intent(inout) :: lhs !< Left hand side.
  class(base_object),  intent(in)    :: rhs !< Right hand side.

  call lhs%assign_(rhs=rhs)
  select type(rhs)
  type is(kappa_int_js)
     if (allocated(rhs%interpolations)) then
        if (.not.allocated(lhs%interpolations)) allocate(lhs%interpolations, mold=rhs%interpolations)
        lhs%interpolations = rhs%interpolations
     else
        if (allocated(lhs%interpolations)) deallocate(lhs%interpolations)
     endif
     if (allocated(rhs%values)) then
        lhs%values = rhs%values
     else
        if (allocated(lhs%values)) deallocate(lhs%values)
     endif
  endselect
  endsubroutine object_assign_object
endmodule wenoof_kappa_int_js
