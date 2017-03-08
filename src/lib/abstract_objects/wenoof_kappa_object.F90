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
  real(RPP), allocatable :: stencil(:)                 !< Stencil used for interpolation, [1-S:S-1].
  real(RPP)              :: x_target                   !< Coordinate of the interpolation point.
endtype kappa_object_constructor

type, extends(base_object), abstract :: kappa_object
  !< Kappa (optimal, linear weights of stencil interpolations) object.
  real(RPP), allocatable :: values_rank_1(:)   !< Kappa coefficients values [0:S-1].
  real(RPP), allocatable :: values_rank_2(:,:) !< Kappa coefficients values [1:2,0:S-1].
  contains
    ! public methods
    generic :: compute => compute_kappa_rec, compute_kappa_int
    ! deferred public methods
    procedure(compute_kappa_rec_interface), pass(self), deferred :: compute_kappa_rec!< Compute interp.
    procedure(compute_kappa_int_interface), pass(self), deferred :: compute_kappa_int!< Compute interp.
endtype kappa_object

abstract interface
  !< Abstract interfaces of [[kappa_object]].
  pure subroutine compute_kappa_rec_interface(self)
  !< Compute kappa.
  import :: kappa_object
  class(kappa_object), intent(inout) :: self        !< Kappa.
  endsubroutine compute_kappa_rec_interface

  pure subroutine compute_kappa_int_interface(self, stencil, x_target)
  !< Compute kappa.
  import :: kappa_object, RPP
  class(kappa_object), intent(inout) :: self        !< Kappa.
  real(RPP),           intent(in)    :: stencil(:)  !< Stencil used for interpolation, [1-S:S-1].
  real(RPP),           intent(in)    :: x_target    !< Coordinate of the interpolation point.
  endsubroutine compute_kappa_int_interface
endinterface

endmodule wenoof_kappa_object
