!< Abstract weights object.
module wenoof_weights_object
!< Abstract weights object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: weights_object
public :: weights_object_constructor

type, extends(base_object_constructor) :: weights_object_constructor
  !< Abstract weights object constructor.
endtype weights_object_constructor

type, extends(base_object), abstract :: weights_object
  !< Weights of stencil interpolations object.
  real(R_P), allocatable :: values(:,:) !< Weights values of stencil interpolations [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute weights.
    ! public overridable methods
    procedure, pass(self) :: create  !< Create weights.
    procedure, pass(self) :: destroy !< Destroy weights.
endtype weights_object

abstract interface
  !< Abstract interfaces of [[weights_object]].
  ! subroutine create_interface(self, constructor)
  ! !< Create weights.
  ! import :: weights_object, base_object_constructor
  ! class(weights_object),          intent(inout) :: self        !< Weights.
  ! class(base_object_constructor), intent(in)    :: constructor !< Weights constructor.
  ! endsubroutine create_interface

  pure subroutine compute_interface(self, stencil)
  !< Compute beta.
  import :: weights_object, R_P
  class(weights_object), intent(inout) :: self                  !< Weights.
  real(R_P),             intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  endsubroutine compute_interface

  ! pure function description_interface(self) result(string)
  ! !< Return weights string-description.
  ! import :: weights_object
  ! class(weights_object), intent(in) :: self   !< Weights.
  ! character(len=:), allocatable     :: string !< String-description.
  ! endfunction description_interface

  ! elemental subroutine destroy_interface(self)
  ! !< Destroy weights.
  ! import :: weights_object
  ! class(beta_object), intent(inout) :: self !< Weights.
  ! endsubroutine destroy_interface
endinterface

contains
  ! public overridable methods
  subroutine create(self, constructor)
  !< Create weights.
  class(weights_object),          intent(inout) :: self        !< Weights.
  class(base_object_constructor), intent(in)    :: constructor !< Weights constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy weights.
  class(weights_object), intent(inout) :: self !< Weights.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_weights_object
