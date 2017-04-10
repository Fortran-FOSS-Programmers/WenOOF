!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.
module wenoof_kappa_object
!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
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
  pure subroutine compute_int_interface(self, stencil, x_target, values)
  !< Compute kappa (interpolate).
  import :: kappa_object, I_P, RPP
  class(kappa_object), intent(in)  :: self               !< Kappa.
  real(RPP),           intent(in)  :: stencil(1-self%S:) !< Stencil used for interpolation, [1-S:S-1].
  real(RPP),           intent(in)  :: x_target           !< Coordinate of the interpolation point.
  real(RPP),           intent(out) :: values(0:)         !< Kappa values.
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, values)
  !< Compute kappa (reconstruct).
  import :: kappa_object, RPP
  class(kappa_object), intent(in)  :: self          !< Kappa.
  real(RPP),           intent(out) :: values(1:,0:) !< Kappa values.
  endsubroutine compute_rec_interface
endinterface
endmodule wenoof_kappa_object
