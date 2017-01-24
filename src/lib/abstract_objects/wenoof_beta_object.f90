!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.
module wenoof_beta_object
!< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: beta_object
public :: beta_object_constructor

type, extends(base_object_constructor) :: beta_object_constructor
  !< Abstract Beta coefficients object constructor.
endtype beta_object_constructor

type, extends(base_object), abstract :: beta_object
  !< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.
  real(R_P), allocatable :: values(:,:) !< Beta values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute beta.
endtype beta_object

abstract interface
  !< Abstract interfaces of [[beta_object]].
  pure subroutine compute_interface(self, stencil)
  !< Compute beta.
  import :: beta_object, R_P
  class(beta_object), intent(inout) :: self                  !< Beta.
  real(R_P),          intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_interface
endinterface

endmodule wenoof_beta_object
