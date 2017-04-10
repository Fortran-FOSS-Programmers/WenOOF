!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.
module wenoof_beta_object
!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
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
  pure subroutine compute_int_interface(self, stencil, values)
  !< Compute beta (interpolate).
  import :: beta_object, RPP
  class(beta_object), intent(in)  :: self               !< Beta.
  real(RPP),          intent(in)  :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(RPP),          intent(out) :: values(0:)         !< Beta values [0:S-1].
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, stencil, values)
  !< Compute beta (reconstruct).
  import :: beta_object, RPP
  class(beta_object), intent(in)  :: self                  !< Beta.
  real(RPP),          intent(in)  :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(RPP),          intent(out) :: values(1:,0:)         !< Beta values [1:2,0:S-1].
  endsubroutine compute_rec_interface
endinterface
endmodule wenoof_beta_object
