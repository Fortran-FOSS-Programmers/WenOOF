!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for interpolation.
module wenoof_kappa_int_js
!< Jiang-Shu and Gerolymos-Senechal-Vallet kappa coefficients for interpolation.
!<
!< @note The provided WENO kappa implements the linear weights defined in *High Order Weighted Essentially
!< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
!< doi:10.1137/070679065.

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
use wenoof_base_object
use wenoof_kappa_object

implicit none
private
public :: kappa_int_js
public :: kappa_int_js_constructor

type, extends(kappa_object_constructor) :: kappa_int_js_constructor
  !< Jiang-Shu and Gerolymos-Senechal-Vallet optimal kappa object constructor.
endtype kappa_int_js_constructor

type, extends(kappa_object):: kappa_int_js
  !< Jiang-Shu and Gerolymos-Senechal-Vallet kappa object.
  !<
  !< @note The provided WENO kappa implements the linear weights defined in *High Order Weighted Essentially
  !< Nonoscillatory Schemes for Convection Dominated Problems*, Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126,
  !< doi:10.1137/070679065.
  contains
    ! public deferred methods
    procedure, pass(self) :: create      !< Create kappa.
    procedure, pass(self) :: compute     !< Compute kappa.
    procedure, pass(self) :: description !< Return kappa string-description.
    procedure, pass(self) :: destroy     !< Destroy kappa.
endtype kappa_int_js

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create kappa.
  !<
  !< @note The kappa coefficients are also computed, they being constants.
  class(kappa_int_js),            intent(inout) :: self        !< Kappa.
  class(base_object_constructor), intent(in)    :: constructor !< Kappa constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._RPP
  call self%compute
  endsubroutine create

  pure subroutine compute(self)
  !< Compute kappa.
  class(kappa_int_js), intent(inout) :: self !< Kappa.

  associate(val => self%values)
    select case(self%S)
      case(2) ! 3rd order
        ! 1 => left interface (i-1/2)
        val(1, 0) = 3._RPP/4._RPP ! stencil 0
        val(1, 1) = 1._RPP/4._RPP ! stencil 1
        ! 2 => right interface (i+1/2)
        val(2, 0) = 1._RPP/4._RPP ! stencil 0
        val(2, 1) = 3._RPP/4._RPP ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        val(1, 0) = 5._RPP/16._RPP ! stencil 0
        val(1, 1) = 5._RPP/8._RPP  ! stencil 1
        val(1, 2) = 1._RPP/16._RPP ! stencil 2
        ! 2 => right interface (i+1/2)
        val(2, 0) = 1._RPP/16._RPP ! stencil 0
        val(2, 1) = 5._RPP/8._RPP  ! stencil 1
        val(2, 2) = 5._RPP/16._RPP ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =  7._RPP/64._RPP ! stencil 0
        val(1, 1) = 35._RPP/64._RPP ! stencil 1
        val(1, 2) = 21._RPP/64._RPP ! stencil 2
        val(1, 3) =  1._RPP/64._RPP ! stencil 3
        ! 2 => right interface (i+1/2)
        val(2, 0) =  1._RPP/64._RPP ! stencil 0
        val(2, 1) = 21._RPP/64._RPP ! stencil 1
        val(2, 2) = 35._RPP/64._RPP ! stencil 2
        val(2, 3) =  7._RPP/64._RPP ! stencil 3
      case(5) ! 9th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =  9._RPP/256._RPP ! stencil 0
        val(1, 1) = 21._RPP/64._RPP  ! stencil 1
        val(1, 2) = 63._RPP/128._RPP ! stencil 2
        val(1, 3) =  9._RPP/64._RPP  ! stencil 3
        val(1, 4) =  1._RPP/256._RPP ! stencil 4
        ! 2 => right interface (i+1/2)
        val(2, 0) =  1._RPP/256._RPP ! stencil 0
        val(2, 1) =  9._RPP/64._RPP  ! stencil 1
        val(2, 2) = 63._RPP/128._RPP  ! stencil 2
        val(2, 3) = 21._RPP/64._RPP  ! stencil 3
        val(2, 4) =  9._RPP/256._RPP ! stencil 4
      case(6) ! 11th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =  11._RPP/1024._RPP  ! stencil 0
        val(1, 1) = 165._RPP/1024._RPP  ! stencil 1
        val(1, 2) = 231._RPP/512._RPP   ! stencil 2
        val(1, 3) = 165._RPP/512._RPP   ! stencil 3
        val(1, 4) =  55._RPP/1024._RPP  ! stencil 4
        val(1, 5) =   1._RPP/1024._RPP  ! stencil 5
        ! 2 => right interface (i+1/2)
        val(2, 0) =   1._RPP/1024._RPP  ! stencil 0
        val(2, 1) =  55._RPP/1024._RPP  ! stencil 1
        val(2, 2) = 165._RPP/512._RPP   ! stencil 2
        val(2, 3) = 231._RPP/512._RPP   ! stencil 3
        val(2, 4) = 165._RPP/1024._RPP  ! stencil 4
        val(2, 5) =  11._RPP/1024._RPP  ! stencil 5
      case(7) ! 13th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   13._RPP/4096._RPP  ! stencil 0
        val(1, 1) =  143._RPP/2048._RPP  ! stencil 1
        val(1, 2) = 1287._RPP/4096._RPP  ! stencil 2
        val(1, 3) =  429._RPP/1024._RPP  ! stencil 3
        val(1, 4) =  179._RPP/1024._RPP  ! stencil 4
        val(1, 5) =   39._RPP/2048._RPP  ! stencil 5
        val(1, 6) =    1._RPP/4096._RPP  ! stencil 6
        ! 2 => right interface (i+1/2)
        val(2, 0) =    1._RPP/4096._RPP  ! stencil 0
        val(2, 1) =   39._RPP/2048._RPP  ! stencil 1
        val(2, 2) =  179._RPP/1024._RPP  ! stencil 2
        val(2, 3) =  429._RPP/1024._RPP  ! stencil 3
        val(2, 4) = 1287._RPP/4096._RPP  ! stencil 4
        val(2, 5) =  143._RPP/2048._RPP  ! stencil 5
        val(2, 6) =   13._RPP/4096._RPP  ! stencil 6
      case(8) ! 15th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =   15._RPP/16384._RPP ! stencil 0
        val(1, 1) =  455._RPP/16384._RPP ! stencil 1
        val(1, 2) = 3003._RPP/16384._RPP ! stencil 2
        val(1, 3) = 6435._RPP/16384._RPP ! stencil 3
        val(1, 4) = 5005._RPP/16384._RPP ! stencil 4
        val(1, 5) = 1365._RPP/16384._RPP ! stencil 5
        val(1, 6) =  105._RPP/16384._RPP ! stencil 6
        val(1, 7) =    1._RPP/16384._RPP ! stencil 7
        ! 2 => right interface (i+1/2)
        val(2, 0) =    1._RPP/16384._RPP ! stencil 0
        val(2, 1) =  105._RPP/16384._RPP ! stencil 1
        val(2, 2) = 1365._RPP/16384._RPP ! stencil 2
        val(2, 3) = 5005._RPP/16384._RPP ! stencil 3
        val(2, 4) = 6435._RPP/16384._RPP ! stencil 4
        val(2, 5) = 3003._RPP/16384._RPP ! stencil 5
        val(2, 6) =  455._RPP/16384._RPP ! stencil 6
        val(2, 7) =   15._RPP/16384._RPP ! stencil 7
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        val(1, 0) =    17._RPP/65536._RPP  ! stencil 0
        val(1, 1) =    85._RPP/8192._RPP   ! stencil 1
        val(1, 2) =  1547._RPP/16384._RPP  ! stencil 2
        val(1, 3) =  2431._RPP/8192._RPP   ! stencil 3
        val(1, 4) = 12155._RPP/32768._RPP  ! stencil 4
        val(1, 5) =  1547._RPP/8192._RPP   ! stencil 5
        val(1, 6) =   595._RPP/16384._RPP  ! stencil 6
        val(1, 7) =    17._RPP/8192._RPP   ! stencil 7
        val(1, 8) =     1._RPP/65536._RPP  ! stencil 8
        ! 2 => right interface (i+1/2)
        val(2, 0) =     1._RPP/65536._RPP  ! stencil 0
        val(2, 1) =    17._RPP/8192._RPP   ! stencil 1
        val(2, 2) =   595._RPP/16384._RPP  ! stencil 2
        val(2, 3) =  1547._RPP/8192._RPP   ! stencil 3
        val(2, 4) = 12155._RPP/32768._RPP  ! stencil 4
        val(2, 5) =  2431._RPP/8192._RPP   ! stencil 5
        val(2, 6) =  1547._RPP/16384._RPP  ! stencil 6
        val(2, 7) =    85._RPP/8192._RPP   ! stencil 7
        val(2, 8) =    17._RPP/65536._RPP  ! stencil 8
    endselect
  endassociate
  endsubroutine compute

  pure function description(self) result(string)
  !< Return string-description of kappa.
  class(kappa_int_js), intent(in) :: self   !< Kappa.
  character(len=:), allocatable   :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'kappa_int_js%description to be implemented, do not use!'
#endif
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy kappa.
  class(kappa_int_js), intent(inout) :: self !< Kappa.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_kappa_int_js
