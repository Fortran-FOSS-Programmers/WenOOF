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
  contains
    ! public methods
    generic :: compute               => compute_with_stencil_of_rank_1, compute_with_stencil_of_rank_2
    generic :: smoothness_indicators => smoothness_indicators_of_rank_1, smoothness_indicators_of_rank_2
    ! deferred public methods
    procedure(compute_with_stencil_of_rank_1_interface),  pass(self), deferred :: compute_with_stencil_of_rank_1  !< Compute beta.
    procedure(compute_with_stencil_of_rank_2_interface),  pass(self), deferred :: compute_with_stencil_of_rank_2  !< Compute beta.
    procedure(smoothness_indicators_of_rank_1_interface), pass(self), deferred :: smoothness_indicators_of_rank_1 !< Return IS.
    procedure(smoothness_indicators_of_rank_2_interface), pass(self), deferred :: smoothness_indicators_of_rank_2 !< Return IS.
endtype weights_object

abstract interface
  !< Abstract interfaces of [[weights_object]].
  pure subroutine compute_with_stencil_of_rank_1_interface(self, stencil)
  !< Compute beta.
  import :: weights_object, RPP
  class(weights_object), intent(inout) :: self               !< Weights.
  real(RPP),             intent(in)    :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  endsubroutine compute_with_stencil_of_rank_1_interface

  pure subroutine compute_with_stencil_of_rank_2_interface(self, stencil)
  !< Compute beta.
  import :: weights_object, RPP
  class(weights_object), intent(inout) :: self                  !< Weights.
  real(RPP),             intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_with_stencil_of_rank_2_interface

  pure subroutine smoothness_indicators_of_rank_1_interface(self, si)
  !< Return smoothness indicators.
  import :: weights_object, RPP
  class(weights_object),  intent(in)  :: self  !< Weights.
  real(RPP), allocatable, intent(out) :: si(:) !< Smoothness indicators.
  endsubroutine smoothness_indicators_of_rank_1_interface

  pure subroutine smoothness_indicators_of_rank_2_interface(self, si)
  !< Return smoothness indicators.
  import :: weights_object, RPP
  class(weights_object),  intent(in)  :: self    !< Weights.
  real(RPP), allocatable, intent(out) :: si(:,:) !< Smoothness indicators.
  endsubroutine smoothness_indicators_of_rank_2_interface
endinterface

endmodule wenoof_weights_object
