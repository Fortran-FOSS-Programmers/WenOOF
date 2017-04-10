!< Abstract interpolations object.
module wenoof_interpolations_object
!< Abstract interpolations object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: interpolations_object
public :: interpolations_object_constructor

type, extends(base_object_constructor), abstract :: interpolations_object_constructor
  !< Abstract interpolations object constructor.
  real(RPP), allocatable :: stencil(:) !< Stencil used for interpolation, [1-S:S-1].
  real(RPP)              :: x_target   !< Coordinate of the interpolation point.
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
  pure subroutine compute_int_interface(self, stencil, values)
  !< Compute interpolations (interpolate).
  import :: interpolations_object, RPP
  class(interpolations_object), intent(in)  :: self               !< Interpolations.
  real(RPP),                    intent(in)  :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  real(RPP),                    intent(out) :: values(0:)         !< Interpolations values.
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, stencil, values)
  !< Compute interpolations (reconstruct).
  import :: interpolations_object, RPP
  class(interpolations_object), intent(in)  :: self                  !< Interpolations.
  real(RPP),                    intent(in)  :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(RPP),                    intent(out) :: values(1:, 0:)        !< Interpolations values.
  endsubroutine compute_rec_interface
endinterface
endmodule wenoof_interpolations_object
