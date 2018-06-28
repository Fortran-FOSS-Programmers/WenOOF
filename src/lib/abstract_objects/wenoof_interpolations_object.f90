!< Abstract interpolations object.
module wenoof_interpolations_object
!< Abstract interpolations object.

use penf, only : I_P, R_P
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: interpolations_object
public :: interpolations_object_constructor

type, extends(base_object_constructor), abstract :: interpolations_object_constructor
  !< Abstract interpolations object constructor.
endtype interpolations_object_constructor

type, extends(base_object), abstract :: interpolations_object
  !< Abstract interpolations object.
  contains
    ! public methods
    generic :: compute => compute_int, compute_rec !< Compute interpolations.
    ! deferred public methods
    procedure(compute_int_interface), pass(self), deferred :: compute_int !< Compute interpolations (interpolate).
    procedure(compute_rec_interface), pass(self), deferred :: compute_rec !< Compute interpolations (reconstruct).
endtype interpolations_object

abstract interface
  !< Abstract interfaces of [[interpolations_object]].
  pure subroutine compute_int_interface(self, ord, stencil, values)
  !< Compute interpolations (interpolate).
  import :: interpolations_object, I_P, R_P
  class(interpolations_object), intent(in)  :: self            !< Interpolations.
  integer(I_P),                 intent(in)  :: ord             !< Interpolation order.
  real(R_P),                    intent(in)  :: stencil(1-ord:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(R_P),                    intent(out) :: values(0:)      !< Interpolations values.
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, ord, stencil, values)
  !< Compute interpolations (reconstruct).
  import :: interpolations_object, I_P, R_P
  class(interpolations_object), intent(in)  :: self               !< Interpolations.
  integer(I_P),                 intent(in)  :: ord                !< Interpolation order.
  real(R_P),                    intent(in)  :: stencil(1:,1-ord:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(R_P),                    intent(out) :: values(1:, 0:)     !< Interpolations values.
  endsubroutine compute_rec_interface
endinterface
endmodule wenoof_interpolations_object
