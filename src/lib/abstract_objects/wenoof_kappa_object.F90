!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.
module wenoof_kappa_object
!< Abstract Kappa (optimal, linear weights of stencil interpolations) object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
use wenoof_base_object

implicit none
private
public :: kappa_object
public :: kappa_object_constructor

type, extends(base_object_constructor), abstract :: kappa_object_constructor
  !< Abstract kappa object constructor.
endtype kappa_object_constructor

type, extends(base_object), abstract :: kappa_object
  !< Kappa (optimal, linear weights of stencil interpolations) object.
  contains
    ! public methods
    generic :: compute => compute_kappa_int, compute_kappa_rec
    ! deferred public methods
    procedure(compute_kappa_int_interface), pass(self), deferred :: compute_kappa_int!< Compute beta.
    procedure(compute_kappa_rec_interface), pass(self), deferred :: compute_kappa_rec!< Compute beta.
endtype kappa_object

abstract interface
  !< Abstract interfaces of [[kappa_object]].
  pure subroutine compute_kappa_int_interface(self)
  !< Compute kappa.
  import :: kappa_object
  class(kappa_object), intent(inout) :: self !< Kappa.
  endsubroutine compute_kappa_int_interface

  pure subroutine compute_kappa_rec_interface(self)
  !< Compute kappa.
  import :: kappa_object
  class(kappa_object), intent(inout) :: self !< Kappa.
  endsubroutine compute_kappa_rec_interface
endinterface

endmodule wenoof_kappa_object
