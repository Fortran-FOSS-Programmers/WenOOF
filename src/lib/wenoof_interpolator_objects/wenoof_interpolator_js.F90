!< Jiang-Shu (upwind) interpolator object.
module wenoof_interpolator_js
!< Jiang-Shu (upwind) interpolator object.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P, str
use wenoof_interpolator_object

implicit none
private
public :: interpolator_js
public :: interpolator_js_constructor

type, extends(interpolator_constructor) :: interpolator_js_constructor
  !< Jiang-Shu (upwind) interpolator object constructor.
endtype interpolator_js_constructor

type, extends(interpolator) :: interpolator_js
  !< Jiang-Shu (upwind) interpolator object.
  contains
    ! public deferred methods
    procedure, pass(self) :: description          !< Return interpolator string-description.
    procedure, pass(self) :: interpolate_standard !< Interpolate values (without providing debug values).
    procedure, pass(self) :: interpolate_debug    !< Interpolate values (providing also debug values).
endtype interpolator_js

contains
  ! public deferred methods
  pure function description(self) result(string)
  !< Return interpolator string-descripition.
  class(interpolator_js), intent(in) :: self             !< Interpolator.
  character(len=:), allocatable      :: string           !< String-description.
  character(len=1), parameter        :: nl=new_line('a') !< New line character.
  character(len=:), allocatable      :: dummy_string     !< Dummy string.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator_js to be implemented, do not use!'
#endif
  endfunction description

  pure subroutine interpolate_standard(self, stencil, interpolation)
  !< Interpolate values (without providing debug values).
  class(interpolator_js), intent(inout) :: self                     !< Interpolator.
  real(R_P),              intent(in)    :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),              intent(out)   :: interpolation(1:)        !< Result of the interpolation, [1:2].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator_js to be implemented, do not use!'
#endif
  endsubroutine interpolate_standard

  pure subroutine interpolate_debug(self, stencil, interpolation, si, weights)
  !< Interpolate values (providing also debug values).
  class(interpolator_js), intent(inout) :: self                     !< Interpolator.
  real(R_P),              intent(in)    :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),              intent(out)   :: interpolation(1:)        !< Result of the interpolation, [1:2].
  real(R_P),              intent(out)   :: si(1:, 0:)               !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(R_P),              intent(out)   :: weights(1:, 0:)          !< Weights of the stencils, [1:2, 0:S-1].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator_js to be implemented, do not use!'
#endif
  endsubroutine interpolate_debug
endmodule wenoof_interpolator_js
