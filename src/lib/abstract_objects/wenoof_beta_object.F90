!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.
module wenoof_beta_object
!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
use wenoof_base_object

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
    generic :: compute => compute_with_stencil_of_rank_1, compute_with_stencil_of_rank_2
    ! deferred public methods
    procedure(compute_with_stencil_of_rank_1_interface), pass(self), deferred :: compute_with_stencil_of_rank_1!< Compute beta.
    procedure(compute_with_stencil_of_rank_2_interface), pass(self), deferred :: compute_with_stencil_of_rank_2!< Compute beta.
endtype beta_object

abstract interface
  !< Abstract interfaces of [[beta_object]].
  pure subroutine compute_with_stencil_of_rank_1_interface(self, stencil)
  !< Compute beta.
  import :: beta_object, RPP
  class(beta_object), intent(inout) :: self               !< Beta.
  real(RPP),          intent(in)    :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  endsubroutine compute_with_stencil_of_rank_1_interface

  pure subroutine compute_with_stencil_of_rank_2_interface(self, stencil)
  !< Compute beta.
  import :: beta_object, RPP
  class(beta_object), intent(inout) :: self                  !< Beta.
  real(RPP),          intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_with_stencil_of_rank_2_interface
endinterface

endmodule wenoof_beta_object
