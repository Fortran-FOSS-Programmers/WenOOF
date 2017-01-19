!< Jiang-Shu (upwind) reconstructor object.
module wenoof_reconstructor_js
!< Jiang-Shu (upwind) reconstructor object.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P, str
use wenoof_interpolator_object

implicit none
private
public :: reconstructor_js
public :: reconstructor_js_constructor

type, extends(interpolator_object_constructor) :: reconstructor_js_constructor
  !< Jiang-Shu (upwind) reconstructor object constructor.
endtype reconstructor_js_constructor

type, extends(interpolator) :: reconstructor_js
  !< Jiang-Shu (upwind) reconstructor object.
  !<
  !< @note Provide the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  !<
  !< @note The supported accuracy formal order are: 3rd, 5th, 7th, 9th, 11th, 13th, 15th, 17th  corresponding to use 2, 3, 4, 5, 6,
  !< 7, 8, 9 stencils composed of 2, 3, 4, 5, 6, 7, 8, 9 values, respectively.
  contains
    ! public deferred methods
    procedure, pass(self) :: description          !< Return reconstructor string-description.
    procedure, pass(self) :: interpolate_standard !< Interpolate values (without providing debug values).
    procedure, pass(self) :: interpolate_debug    !< Interpolate values (providing also debug values).
endtype reconstructor_js

contains
  ! public deferred methods
  pure function description(self) result(string)
  !< Return reconstructor string-descripition.
  class(reconstructor_js), intent(in) :: self             !< Reconstructor.
  character(len=:), allocatable       :: string           !< String-description.
  character(len=1), parameter         :: nl=new_line('a') !< New line character.
  character(len=:), allocatable       :: dummy_string     !< Dummy string.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'reconstructor_js%description to be implemented, do not use!'
#endif
  endfunction description

  pure subroutine interpolate_standard(self, stencil, interpolation)
  !< Interpolate values (without providing debug values).
  class(reconstructor_js), intent(inout) :: self                     !< Reconstructor.
  real(R_P),               intent(in)    :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),               intent(out)   :: interpolation(1:)        !< Result of the interpolation, [1:2].
  integer(I_P)                           :: f, s                     !< Counters.

  call self%interpolations%compute(stencil=stencil)
  call self%weights%compute(stencil=stencil)
  interpolation = 0._R_P
  do s=0, self%S - 1 ! stencils loop
    do f=self%f1, self%f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      interpolation(f + self%ff) = interpolation(f + self%ff) + &
                                   self%weights%values(f + self%ff, s) * self%interpolations%values(f, s)
    enddo
  enddo
  endsubroutine interpolate_standard

  pure subroutine interpolate_debug(self, stencil, interpolation, si, weights)
  !< Interpolate values (providing also debug values).
  class(reconstructor_js), intent(inout) :: self                     !< Reconstructor.
  real(R_P),               intent(in)    :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),               intent(out)   :: interpolation(1:)        !< Result of the interpolation, [1:2].
  real(R_P),               intent(out)   :: si(1:, 0:)               !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(R_P),               intent(out)   :: weights(1:, 0:)          !< Weights of the stencils, [1:2, 0:S-1].

  call self%interpolate_standard(stencil=stencil, interpolation=interpolation)
  si = 0._R_P ! @TODO implement beta extraction
  weights = self%weights%values
  endsubroutine interpolate_debug
endmodule wenoof_reconstructor_js
