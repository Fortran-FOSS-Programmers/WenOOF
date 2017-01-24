!< Abstract weights object.
module wenoof_weights_object
!< Abstract weights object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: weights_object
public :: weights_object_constructor

type, extends(base_object_constructor) :: weights_object_constructor
  !< Abstract weights object constructor.
endtype weights_object_constructor

type, extends(base_object), abstract :: weights_object
  !< Weights of stencil interpolations object.
  real(R_P), allocatable :: values(:,:) !< Weights values of stencil interpolations [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute weights.
endtype weights_object

abstract interface
  !< Abstract interfaces of [[weights_object]].
  pure subroutine compute_interface(self, stencil)
  !< Compute beta.
  import :: weights_object, R_P
  class(weights_object), intent(inout) :: self                  !< Weights.
  real(R_P),             intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_interface
endinterface

endmodule wenoof_weights_object
