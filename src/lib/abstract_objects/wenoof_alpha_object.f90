!< Abstract alpha (non linear weights) object.
module wenoof_alpha_object
!< Abstract alpha (non linear weights) object.

use penf, only : I_P, R_P
use wenoof_base_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_object
public :: alpha_object_constructor

type, extends(base_object_constructor) :: alpha_object_constructor
  !< Abstract alpha (non linear weights) object constructor.
  contains
endtype alpha_object_constructor

type, extends(base_object), abstract :: alpha_object
  !< Abstract alpha (non linear weights) object.
  real(R_P), allocatable :: values(:,:)   !< Alpha coefficients [1:2,0:S-1].
  real(R_P), allocatable :: values_sum(:) !< Sum of alpha coefficients [1:2].
  contains
    ! public deferred methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute alpha.
endtype alpha_object

abstract interface
  !< Abstract interfaces of [[alpha_object]].
  pure subroutine compute_interface(self, beta, kappa)
  !< Compute alpha.
  import :: alpha_object, beta_object, kappa_object
  class(alpha_object), intent(inout) :: self  !< Alpha.
  class(beta_object),  intent(in)    :: beta  !< Beta.
  class(kappa_object), intent(in)    :: kappa !< Kappa.
  endsubroutine compute_interface
endinterface

endmodule wenoof_alpha_object
