!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.
module wenoof_kappa_object
!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.

use penf, only : I_P, R_P
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: kappa_object
public :: kappa_object_constructor

type, extends(base_object_constructor), abstract :: kappa_object_constructor
  !< Abstract kappa object constructor.
endtype kappa_object_constructor

type, extends(base_object), abstract :: kappa_object
  !< Kappa (optimal, linear weights of stencil interpolations) object.
  contains
    ! public methods
    generic :: compute => compute_rec, compute_int !< Compute kappa.
    ! deferred public methods
    procedure(compute_int_interface), pass(self), deferred :: compute_int !< Compute kappa (interpolate).
    procedure(compute_rec_interface), pass(self), deferred :: compute_rec !< Compute kappa (reconstruct).
endtype kappa_object

abstract interface
  !< Abstract interfaces of [[kappa_object]].
  pure subroutine compute_int_interface(self, x_target, values)
  !< Compute kappa (interpolate).
  import :: kappa_object, R_P
  class(kappa_object), intent(in)  :: self               !< Kappa.
  real(R_P),           intent(in)  :: x_target           !< Coordinate of the interpolation point.
  real(R_P),           intent(out) :: values(0:)         !< Kappa values.
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, values)
  !< Compute kappa (reconstruct).
  import :: kappa_object, R_P
  class(kappa_object), intent(in)  :: self          !< Kappa.
  real(R_P),           intent(out) :: values(1:,0:) !< Kappa values.
  endsubroutine compute_rec_interface
endinterface
endmodule wenoof_kappa_object
