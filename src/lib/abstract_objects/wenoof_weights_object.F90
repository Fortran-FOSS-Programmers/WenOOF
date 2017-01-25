!< Abstract weights object.
module wenoof_weights_object
!< Abstract weights object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
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
  real(RPP), allocatable :: values(:,:) !< Weights values of stencil interpolations [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure(compute_interface),               pass(self), deferred :: compute               !< Compute weights.
    procedure(smoothness_indicators_interface), pass(self), deferred :: smoothness_indicators !< Return smoothness indicators.
endtype weights_object

abstract interface
  !< Abstract interfaces of [[weights_object]].
  pure subroutine compute_interface(self, stencil)
  !< Compute beta.
  import :: weights_object, RPP
  class(weights_object), intent(inout) :: self                  !< Weights.
  real(RPP),             intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_interface

  pure function smoothness_indicators_interface(self) result(si)
  !< Return smoothness indicators.
  import :: weights_object, RPP
  class(weights_object), intent(in) :: self    !< Weights.
  real(RPP), allocatable            :: si(:,:) !< Smoothness indicators.
  endfunction smoothness_indicators_interface
endinterface

endmodule wenoof_weights_object
