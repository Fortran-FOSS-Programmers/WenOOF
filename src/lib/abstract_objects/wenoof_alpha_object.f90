!< Abstract alpha (non linear weights) object.
module wenoof_alpha_object
!< Abstract alpha (non linear weights) object.

use penf, only : R_P
use wenoof_base_object, only : base_object, base_object_constructor

implicit none
private
public :: alpha_object
public :: alpha_object_constructor

type, extends(base_object_constructor), abstract :: alpha_object_constructor
  !< Abstract alpha (non linear weights) object constructor.
  contains
endtype alpha_object_constructor

type, extends(base_object), abstract :: alpha_object
  !< Abstract alpha (non linear weights) object.
  contains
    ! public methods
    generic :: compute => compute_int, compute_rec !< Compute alpha.
    ! public deferred methods
    procedure(compute_int_interface), pass(self), deferred :: compute_int !< Compute alpha (interpolate).
    procedure(compute_rec_interface), pass(self), deferred :: compute_rec !< Compute alpha (reconstruct).
endtype alpha_object

abstract interface
  !< Abstract interfaces of [[alpha_object]].
  pure subroutine compute_int_interface(self, beta, kappa, values)
  !< Compute alpha (interpolate).
  import :: alpha_object, R_P
  class(alpha_object), intent(in)  :: self       !< Alpha.
  real(R_P),           intent(in)  :: beta(0:)   !< Beta [0:S-1].
  real(R_P),           intent(in)  :: kappa(0:)  !< Kappa [0:S-1].
  real(R_P),           intent(out) :: values(0:) !< Alpha values [0:S-1].
  endsubroutine compute_int_interface

  pure subroutine compute_rec_interface(self, beta, kappa, values)
  !< Compute alpha (reconstruct).
  import :: alpha_object, R_P
  class(alpha_object), intent(in)  :: self          !< Alpha.
  real(R_P),           intent(in)  :: beta(1:,0:)   !< Beta [1:2,0:S-1].
  real(R_P),           intent(in)  :: kappa(1:,0:)  !< Kappa [1:2,0:S-1].
  real(R_P),           intent(out) :: values(1:,0:) !< Alpha values [1:2,0:S-1].
  endsubroutine compute_rec_interface
endinterface
endmodule wenoof_alpha_object
