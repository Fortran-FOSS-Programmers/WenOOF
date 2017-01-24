!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.
module wenoof_kappa_object
!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: kappa_object
public :: kappa_object_constructor

type, extends(base_object_constructor) :: kappa_object_constructor
  !< Abstract kappa object constructor.
endtype kappa_object_constructor

type, extends(base_object), abstract :: kappa_object
  !< Kappa (optimal, linear weights of stencil interpolations) object.
  real(R_P), allocatable :: values(:,:) !< Kappa coefficients values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute kappa.
endtype kappa_object

abstract interface
  !< Abstract interfaces of [[kappa_object]].
  pure subroutine compute_interface(self)
  !< Compute kappa.
  import :: kappa_object
  class(kappa_object), intent(inout) :: self !< Kappa.
  endsubroutine compute_interface
endinterface

endmodule wenoof_kappa_object
