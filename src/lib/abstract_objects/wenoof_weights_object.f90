!< Abstract weights object.
module wenoof_weights_object
!< Abstract weights object.

use penf, only : I_P, R_P
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: weights_object
public :: weights_object_constructor

type, extends(base_object_constructor), abstract :: weights_object_constructor
  !< Abstract weights object constructor.
endtype weights_object_constructor

type, extends(base_object), abstract :: weights_object
  !< Weights of stencil interpolations object.
  contains
    ! public methods
    generic :: compute               => compute_int, compute_rec                             !< Compute weights.
    generic :: smoothness_indicators => smoothness_indicators_int, smoothness_indicators_rec !< Return IS.
    ! deferred public methods
    procedure(compute_int_interface),               pass(self), deferred :: compute_int               !< Compute weights (interp.).
    procedure(compute_rec_interface),               pass(self), deferred :: compute_rec               !< Compute weights (recons.).
    procedure(smoothness_indicators_int_interface), pass(self), deferred :: smoothness_indicators_int !< Return IS.
    procedure(smoothness_indicators_rec_interface), pass(self), deferred :: smoothness_indicators_rec !< Return IS.
endtype weights_object

abstract interface
  !< Abstract interfaces of [[weights_object]].
  pure subroutine compute_int_interface(self, ord, stencil, values)
  !< Compute weights (interpolate).
  import :: weights_object, I_P, R_P
  class(weights_object), intent(in)  :: self            !< Weights.
  integer(I_P),          intent(in)  :: ord             !< Interpolation order.
  real(R_P),             intent(in)  :: stencil(1-ord:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(R_P),             intent(out) :: values(0:)      !< Weights values of stencil interpolations.
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, ord, stencil, values)
  !< Compute beta (reconstruct).
  import :: weights_object, I_P, R_P
  class(weights_object), intent(in)  :: self               !< Weights.
  integer(I_P),          intent(in)  :: ord                !< Reconstruction order.
  real(R_P),             intent(in)  :: stencil(1:,1-ord:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(R_P),             intent(out) :: values(1:,0:)      !< Weights values of stencil interpolations.
  endsubroutine compute_rec_interface

  pure subroutine smoothness_indicators_int_interface(self, si)
  !< Return smoothness indicators (interpolate).
  import :: weights_object, R_P
  class(weights_object),  intent(in)  :: self  !< Weights.
  real(R_P),              intent(out) :: si(:) !< Smoothness indicators.
  endsubroutine smoothness_indicators_int_interface

  pure subroutine smoothness_indicators_rec_interface(self, si)
  !< Return smoothness indicators (reconstruct).
  import :: weights_object, R_P
  class(weights_object),  intent(in)  :: self      !< Weights.
  real(R_P),              intent(out) :: si(1:,0:) !< Smoothness indicators.
  endsubroutine smoothness_indicators_rec_interface
endinterface
endmodule wenoof_weights_object
