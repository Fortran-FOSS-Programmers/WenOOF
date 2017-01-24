!< Abstract interpolations object.
module wenoof_interpolations_object
!< Abstract interpolations object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: interpolations_object
public :: interpolations_object_constructor

type, extends(base_object_constructor) :: interpolations_object_constructor
  !< Abstract interpolations object constructor.
endtype interpolations_object_constructor

type, extends(base_object), abstract :: interpolations_object
  !< Abstract interpolations object.
  real(R_P), allocatable :: values(:,:) !< Stencil interpolations values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute beta.
    ! public overridable methods
    procedure, pass(self) :: create  !< Create interpolations.
    procedure, pass(self) :: destroy !< Destroy interpolations.
endtype interpolations_object

abstract interface
  !< Abstract interfaces of [[interpolations_object]].
  ! subroutine create_interface(self, constructor)
  ! !< Create interpolations.
  ! import :: interpolations_object, base_object_constructor
  ! class(interpolations_object),    intent(inout) :: self        !< Interpolations.
  ! class(base_object_constructor),  intent(in)    :: constructor !< Interpolations constructor.
  ! endsubroutine create_interface

  pure subroutine compute_interface(self, stencil)
  !< Compute interpolations.
  import :: interpolations_object, R_P
  class(interpolations_object), intent(inout) :: self                  !< Interpolations.
  real(R_P),                    intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_interface

  ! pure function description_interface(self) result(string)
  ! !< Return interpolations string-description.
  ! import :: interpolations_object
  ! class(interpolations_object), intent(in) :: self   !< Interpolations.
  ! character(len=:), allocatable            :: string !< String-description.
  ! endfunction description_interface

  ! elemental subroutine destroy_interface(self)
  ! !< Destroy interpolations.
  ! import :: interpolations_object
  ! class(interpolations_object), intent(inout) :: self !< Interpolations.
  ! endsubroutine destroy_interface
endinterface

contains
  ! public overridable methods
  subroutine create(self, constructor)
  !< Create interpolations.
  class(interpolations_object),   intent(inout) :: self        !< Interpolations.
  class(base_object_constructor), intent(in)    :: constructor !< Interpolations constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy interpolations.
  class(interpolations_object), intent(inout) :: self !< Interpolations.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_interpolations_object
