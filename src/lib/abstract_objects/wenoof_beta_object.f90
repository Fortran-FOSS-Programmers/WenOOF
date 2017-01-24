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
    ! public overridable methods
    procedure, pass(self) :: create  !< Create beta.
    procedure, pass(self) :: destroy !< Destroy beta.
endtype beta_object

abstract interface
  !< Abstract interfaces of [[beta_object]].
  ! subroutine create_interface(self, constructor)
  ! !< Create beta.
  ! import :: beta_object, base_object_constructor
  ! class(beta_object),              intent(inout) :: self        !< Beta.
  ! class(base_object_constructor),  intent(in)    :: constructor !< Beta constructor.
  ! endsubroutine create_interface

  pure subroutine compute_interface(self, stencil)
  !< Compute beta.
  import :: beta_object, R_P
  class(beta_object), intent(inout) :: self                  !< Beta.
  real(R_P),          intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_interface

  ! pure function description_interface(self) result(string)
  ! !< Return beta string-description.
  ! import :: beta_object
  ! class(beta_object), intent(in) :: self   !< Beta.
  ! character(len=:), allocatable  :: string !< String-description.
  ! endfunction description_interface

  ! elemental subroutine destroy_interface(self)
  ! !< Destroy beta.
  ! import :: beta_object
  ! class(beta_object), intent(inout) :: self !< Beta.
  ! endsubroutine destroy_interface
endinterface

contains
  ! public overridable methods
  subroutine create(self, constructor)
  !< Create beta.
  class(beta_object),             intent(inout) :: self        !< Beta.
  class(base_object_constructor), intent(in)    :: constructor !< Beta constructor.

  call self%destroy_
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy beta.
  class(beta_object), intent(inout) :: self !< Beta.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_beta_object
