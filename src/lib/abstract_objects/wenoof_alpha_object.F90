!< Abstract alpha (non linear weights) object.
module wenoof_alpha_object
!< Abstract alpha (non linear weights) object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
use wenoof_base_object
use wenoof_beta_object
use wenoof_kappa_object

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
    ! public methods
    generic :: compute => compute_alpha_int, compute_alpha_rec
    ! deferred public methods
    procedure(compute_alpha_int_interface), pass(self), deferred :: compute_alpha_int!< Compute beta.
    procedure(compute_alpha_rec_interface), pass(self), deferred :: compute_alpha_rec!< Compute beta.
endtype alpha_object

abstract interface
  !< Abstract interfaces of [[alpha_object]].
  pure subroutine compute_alpha_int_interface(self, beta, kappa)
  !< Compute alpha.
  import :: alpha_object, beta_object, kappa_object
  class(alpha_object), intent(inout) :: self  !< Alpha.
  class(beta_object),  intent(in)    :: beta  !< Beta.
  class(kappa_object), intent(in)    :: kappa !< Kappa.
  endsubroutine compute_alpha_int_interface

  pure subroutine compute_alpha_rec_interface(self, beta, kappa)
  !< Compute alpha.
  import :: alpha_object, beta_object, kappa_object
  class(alpha_object), intent(inout) :: self  !< Alpha.
  class(beta_object),  intent(in)    :: beta  !< Beta.
  class(kappa_object), intent(in)    :: kappa !< Kappa.
  endsubroutine compute_alpha_rec_interface
endinterface

endmodule wenoof_alpha_object
