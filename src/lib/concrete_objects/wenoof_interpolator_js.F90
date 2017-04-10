!< Jiang-Shu (upwind) interpolator object.
module wenoof_interpolator_js
!< Jiang-Shu (upwind) interpolator object.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
#ifdef r16p
use penf, only: I_P, RPP=>R16P, str
#else
use penf, only: I_P, RPP=>R8P, str
#endif
use wenoof_base_object, only : base_object_constructor
use wenoof_interpolations_factory, only : interpolations_factory
use wenoof_interpolations_object, only : interpolations_object
use wenoof_interpolator_object, only : interpolator_object, interpolator_object_constructor
use wenoof_weights_factory, only : weights_factory
use wenoof_weights_object, only : weights_object

implicit none
private
public :: interpolator_js
public :: interpolator_js_constructor

type, extends(interpolator_object_constructor) :: interpolator_js_constructor
  !< Jiang-Shu (upwind) interpolator object constructor.
endtype interpolator_js_constructor

type, extends(interpolator_object) :: interpolator_js
  !< Jiang-Shu (upwind) interpolator object.
  !<
  !< @note Provide the *High Order Weighted Essentially Nonoscillatory Schemes for Convection Dominated Problems*,
  !< Chi-Wang Shu, SIAM Review, 2009, vol. 51, pp. 82--126, doi:10.1137/070679065.
  !<
  !< @note The supported accuracy formal order are: 3rd, 5th, 7th, 9th, 11th, 13th, 15th, 17th  corresponding to use 2, 3, 4, 5, 6,
  !< 7, 8, 9 stencils composed of 2, 3, 4, 5, 6, 7, 8, 9 values, respectively.
  contains
    ! public deferred methods
    procedure, pass(self) :: create                   !< Create interpolator.
    procedure, pass(self) :: description              !< Return object string-description.
    procedure, pass(self) :: destroy                  !< Destroy interpolator.
    procedure, pass(self) :: interpolate_int_debug    !< Interpolate values (providing also debug values, interpolate).
    procedure, pass(self) :: interpolate_int_standard !< Interpolate values (without providing debug values, interpolate).
    procedure, pass(self) :: interpolate_rec_debug    !< Interpolate values (providing also debug values, reconstruct).
    procedure, pass(self) :: interpolate_rec_standard !< Interpolate values (without providing debug values, reconstruct).
endtype interpolator_js

contains
  ! public deferred methods
  subroutine create(self, constructor)
  !< Create interpolator.
  class(interpolator_js),         intent(inout) :: self        !< Interpolator.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.
  type(interpolations_factory)                  :: i_factory   !< Inteprolations factory.
  type(weights_factory)                         :: w_factory   !< Weights factory.

  call self%destroy
  call self%create_(constructor=constructor)
  select type(constructor)
  class is(interpolator_object_constructor)
    call i_factory%create(constructor=constructor%interpolations_constructor, object=self%interpolations)
    call w_factory%create(constructor=constructor%weights_constructor, object=self%weights)
  endselect
  endsubroutine create

  pure function description(self, prefix) result(string)
  !< Return object string-descripition.
  class(interpolator_js), intent(in)           :: self             !< Interpolator.
  character(len=*),       intent(in), optional :: prefix           !< Prefixing string.
  character(len=:), allocatable                :: string           !< String-description.
  character(len=:), allocatable                :: prefix_          !< Prefixing string, local variable.
  character(len=1), parameter                  :: NL=new_line('a') !< New line character.

  prefix_ = '' ; if (present(prefix)) prefix_ = prefix
  string = prefix_//'Jiang-Shu interpolator:'//NL
  string = prefix_//string//'  - S   = '//trim(str(self%S))//NL
  string = prefix_//string//self%weights%description(prefix=prefix_//'  ')
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy interpolator.
  class(interpolator_js), intent(inout) :: self !< Interpolator.

  call self%destroy_
  if (allocated(self%interpolations)) deallocate(self%interpolations)
  if (allocated(self%weights)) deallocate(self%weights)
  endsubroutine destroy

  pure subroutine interpolate_int_debug(self, stencil, interpolation, si, weights)
  !< Interpolate values (providing also debug values, interpolate).
  class(interpolator_js), intent(in)  :: self                         !< Interpolator.
  real(RPP),              intent(in)  :: stencil(1 - self%S:)         !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(RPP),              intent(out) :: interpolation                !< Result of the interpolation.
  real(RPP),              intent(out) :: si(0:)                       !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(RPP),              intent(out) :: weights(0:)                  !< Weights of the stencils, [1:2, 0:S-1].
  real(RPP)                           :: interpolations(0:self%S - 1) !< Stencils interpolations.
  integer(I_P)                        :: s                            !< Counters.

  call self%interpolations%compute(stencil=stencil, values=interpolations)
  call self%weights%compute(stencil=stencil, values=weights)
  ! call self%weights%smoothness_indicators_of_rank_1(si=si)
  interpolation = 0._RPP
  do s=0, self%S - 1 ! stencils loop
    interpolation = interpolation + weights(s) * interpolations(s)
  enddo
  endsubroutine interpolate_int_debug

  pure subroutine interpolate_int_standard(self, stencil, interpolation)
  !< Interpolate values (without providing debug values, interpolate).
  class(interpolator_js), intent(in)  :: self                         !< Interpolator.
  real(RPP),              intent(in)  :: stencil(1 - self%S:)         !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(RPP),              intent(out) :: interpolation                !< Result of the interpolation.
  real(RPP)                           :: interpolations(0:self%S - 1) !< Stencils interpolations.
  real(RPP)                           :: weights(0:self%S - 1)        !< Weights of stencils interpolations.
  integer(I_P)                        :: s                            !< Counters.

  call self%interpolations%compute(stencil=stencil, values=interpolations)
  call self%weights%compute(stencil=stencil, values=weights)
  interpolation = 0._RPP
  do s=0, self%S - 1 ! stencils loop
    interpolation = interpolation + weights(s) * interpolations(s)
  enddo
  endsubroutine interpolate_int_standard

  pure subroutine interpolate_rec_debug(self, stencil, interpolation, si, weights)
  !< Interpolate values (providing also debug values, reconstruct).
  class(interpolator_js), intent(in)  :: self                     !< Reconstructor.
  real(RPP),              intent(in)  :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(RPP),              intent(out) :: interpolation(1:)        !< Result of the interpolation, [1:2].
  real(RPP),              intent(out) :: si(1:, 0:)               !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(RPP),              intent(out) :: weights(1:, 0:)          !< Weights of the stencils, [1:2, 0:S-1].
  ! empty procedure
  endsubroutine interpolate_rec_debug

  pure subroutine interpolate_rec_standard(self, stencil, interpolation)
  !< Interpolate values (without providing debug values, reconstruct).
  class(interpolator_js), intent(in)  :: self                     !< Reconstructor.
  real(RPP),              intent(in)  :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(RPP),              intent(out) :: interpolation(1:)        !< Result of the interpolation, [1:2].
  ! empty procedure
  endsubroutine interpolate_rec_standard
endmodule wenoof_interpolator_js
