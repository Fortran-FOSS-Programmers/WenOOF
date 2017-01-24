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
    ! public overridable methods
    procedure, pass(self) :: create  !< Create kappa.
    procedure, pass(self) :: destroy !< Destroy kappa.
endtype kappa_object

abstract interface
  !< Abstract interfaces of [[kappa_object]].
  ! subroutine create_interface(self, constructor)
  ! !< Create kappa.
  ! import :: kappa_object, base_object_constructor
  ! class(kappa_object),            intent(inout) :: self        !< Kappa.
  ! class(base_object_constructor), intent(in)    :: constructor !< Kappa constructor.
  ! endsubroutine create_interface

  pure subroutine compute_interface(self)
  !< Compute kappa.
  import :: kappa_object
  class(kappa_object), intent(inout) :: self !< Kappa.
  endsubroutine compute_interface

  ! pure function description_interface(self) result(string)
  ! !< Return beta string-description.
  ! import :: kappa_object
  ! class(kappa_object), intent(in) :: self   !< Kappa.
  ! character(len=:), allocatable   :: string !< String-description.
  ! endfunction description_interface

  ! elemental subroutine destroy_interface(self)
  ! !< Destroy kappa.
  ! import :: kappa_object
  ! class(kappa_object), intent(inout) :: self !< Kappa.
  ! endsubroutine destroy_interface
endinterface

contains
  ! public overridable methods
  subroutine create(self, constructor)
  !< Create kappa.
  class(kappa_object),            intent(inout) :: self        !< Kappa.
  class(base_object_constructor), intent(in)    :: constructor !< Kappa constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy kappa.
  class(kappa_object), intent(inout) :: self !< Kappa.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_kappa_object
