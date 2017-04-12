!< Abstract interpolator object.
module wenoof_interpolator_object
!< Abstract interpolator object.

use penf, only : R_P
use wenoof_base_object, only : base_object, base_object_constructor
use wenoof_interpolations_object, only : interpolations_object, interpolations_object_constructor
use wenoof_weights_object, only : weights_object, weights_object_constructor

implicit none
private
public :: interpolator_object
public :: interpolator_object_constructor

type, extends(base_object_constructor), abstract :: interpolator_object_constructor
  !< Abstract interpolator object constructor.
  !<
  !< @note Every concrete WENO interpolator implementations must define their own constructor type.
  class(interpolations_object_constructor), allocatable :: interpolations_constructor !< Stencil interpolations constructor.
  class(weights_object_constructor),        allocatable :: weights_constructor        !< Weights of interpolations constructor.
endtype interpolator_object_constructor

type, extends(base_object), abstract :: interpolator_object
  !< Abstract interpolator object.
  !<
  !< @note Do not implement any actual interpolator: provide the interface for the different interpolators implemented.
  class(interpolations_object), allocatable :: interpolations !< Stencil interpolations.
  class(weights_object),        allocatable :: weights        !< Weights of interpolations.
  contains
    ! public methods
    generic :: interpolate => interpolate_int_debug,    interpolate_rec_debug, &
                              interpolate_int_standard, interpolate_rec_standard !< Interpolate.
    ! public deferred methods
    procedure(interpolate_int_debug_interface),    pass(self), deferred :: interpolate_int_debug    !< Interpolate (int) debug.
    procedure(interpolate_int_standard_interface), pass(self), deferred :: interpolate_int_standard !< Interpolate (int) standard.
    procedure(interpolate_rec_debug_interface),    pass(self), deferred :: interpolate_rec_debug    !< Interpolate (rec) debug.
    procedure(interpolate_rec_standard_interface), pass(self), deferred :: interpolate_rec_standard !< Interpolate (rec) standard.
endtype interpolator_object

abstract interface
  !< Abstract interfaces of [[interpolator_object]].
  pure subroutine interpolate_int_debug_interface(self, stencil, interpolation, si, weights)
  !< Interpolate (interpolate) values (providing also debug values).
  import :: interpolator_object, R_P
  class(interpolator_object), intent(in)  :: self                 !< Interpolator.
  real(R_P),                  intent(in)  :: stencil(1 - self%S:) !< Stencil of the interpolation [1-S:-1+S].
  real(R_P),                  intent(out) :: interpolation        !< Result of the interpolation.
  real(R_P),                  intent(out) :: si(0:)               !< Computed values of smoothness indicators [0:S-1].
  real(R_P),                  intent(out) :: weights(0:)          !< Weights of the stencils, [0:S-1].
  endsubroutine interpolate_int_debug_interface

  pure subroutine interpolate_int_standard_interface(self, stencil, interpolation)
  !< Interpolate (interpolate) values (without providing debug values).
  import :: interpolator_object, R_P
  class(interpolator_object), intent(in)  :: self                 !< Interpolator.
  real(R_P),                  intent(in)  :: stencil(1 - self%S:) !< Stencil of the interpolation [1-S:-1+S].
  real(R_P),                  intent(out) :: interpolation        !< Result of the interpolation.
  endsubroutine interpolate_int_standard_interface

  pure subroutine interpolate_rec_debug_interface(self, stencil, interpolation, si, weights)
  !< Interpolate (reconstruct) values (providing also debug values).
  import :: interpolator_object, R_P
  class(interpolator_object), intent(in)  :: self                     !< Interpolator.
  real(R_P),                  intent(in)  :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),                  intent(out) :: interpolation(1:)        !< Result of the interpolation, [1:2].
  real(R_P),                  intent(out) :: si(1:, 0:)               !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(R_P),                  intent(out) :: weights(1:, 0:)          !< Weights of the stencils, [1:2, 0:S-1].
  endsubroutine interpolate_rec_debug_interface

  pure subroutine interpolate_rec_standard_interface(self, stencil, interpolation)
  !< Interpolate (reconstruct) values (without providing debug values).
  import :: interpolator_object, R_P
  class(interpolator_object), intent(in)  :: self                     !< Interpolator.
  real(R_P),                  intent(in)  :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),                  intent(out) :: interpolation(1:)        !< Result of the interpolation, [1:2].
  endsubroutine interpolate_rec_standard_interface
endinterface
endmodule wenoof_interpolator_object
