!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.
module wenoof_beta_object
!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.

use penf, only : I_P, R_P
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: beta_object
public :: beta_object_constructor

type, extends(base_object_constructor), abstract :: beta_object_constructor
  !< Abstract Beta coefficients object constructor.
endtype beta_object_constructor

type, extends(base_object), abstract :: beta_object
  !< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.
  contains
    ! public methods
    generic :: compute => compute_int, compute_rec !< Compute beta.
    ! deferred public methods
    procedure(compute_int_interface), pass(self), deferred :: compute_int !< Compute beta (interpolate).
    procedure(compute_rec_interface), pass(self), deferred :: compute_rec !< Compute beta (reconstruct).
endtype beta_object

abstract interface
  !< Abstract interfaces of [[beta_object]].
  pure subroutine compute_int_interface(self, ord, stencil, values)
  !< Compute beta (interpolate).
  import :: beta_object, I_P, R_P
  class(beta_object), intent(in)  :: self               !< Beta.
  integer(I_P),       intent(in)  :: ord                !< Order of interpolation.
  real(R_P),          intent(in)  :: stencil(1-ord:)    !< Stencil used for the interpolation, [1-S:-1+S].
  real(R_P),          intent(out) :: values(0:)         !< Beta values [0:S-1].
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, ord, stencil, values)
  !< Compute beta (reconstruct).
  import :: beta_object, I_P, R_P
  class(beta_object), intent(in)  :: self                  !< Beta.
  integer(I_P),       intent(in)  :: ord                   !< Order of reconstruction.
  real(R_P),          intent(in)  :: stencil(1:,1-ord:)    !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(R_P),          intent(out) :: values(1:,0:)         !< Beta values [1:2,0:S-1].
  endsubroutine compute_rec_interface
endinterface
endmodule wenoof_beta_object
