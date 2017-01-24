!< Abstract interpolations object.
module wenoof_interpolations_object
!< Abstract interpolations object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: interpolations_object
public :: interpolations_object_constructor

type, extends(base_object_constructor) :: interpolations_object_constructor
  !< Abstract interpolations object constructor.
endtype interpolations_object_constructor

type, extends(base_object), abstract :: interpolations_object
  !< Abstract interpolations object.
  real(R_P), allocatable :: values(:,:) !< Stencil interpolations values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute beta.
endtype interpolations_object

abstract interface
  !< Abstract interfaces of [[interpolations_object]].
  pure subroutine compute_interface(self, stencil)
  !< Compute interpolations.
  import :: interpolations_object, R_P
  class(interpolations_object), intent(inout) :: self                  !< Interpolations.
  real(R_P),                    intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_interface
endinterface

endmodule wenoof_interpolations_object
