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
  integer(I_P) :: S=0_I_P !< Stencils dimension.
endtype weights_object_constructor

type, extends(base_object) :: weights_object
  !< Weights of stencil interpolations object.
  integer(I_P)           :: S=0_I_P     !< Stencils dimension.
  real(R_P), allocatable :: values(:,:) !< Weights values of stencil interpolations [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute weights.
    procedure, pass(self) :: description !< Return weights string-description.
    ! public methods
    procedure, pass(self) :: create  !< Createte weights.
    procedure, pass(self) :: destroy !< Destroy weights.
endtype weights_object

contains
  ! deferred public methods
  pure subroutine compute(self, stencil)
  !< Compute weights.
  class(weights_object), intent(inout) :: self                  !< Weights.
  real(R_P),             intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'optimal_weights%compute to be implemented by your concrete optimal weights object'
#endif
  endsubroutine compute

  pure function description(self) result(string)
  !< Return weights string-description.
  class(weights_object), intent(in) :: self   !< Weights.
  character(len=:), allocatable     :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'optimal_weights%description to be implemented by your concrete optimal weights object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create weights.
  class(weights_object),          intent(inout) :: self        !< Weights.
  class(base_object_constructor), intent(in)    :: constructor !< Weights constructor.

  call self%destroy
  select type(constructor)
  class is(optimal_weights_constructor)
    self%S = constructor%S
    allocate(self%values(1:2, 0:constructor%S - 1))
  class default
    ! @TODO add error handling
  endselect
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy weights.
  class(weights_object), intent(inout) :: self !< Weights.

  self%S = 0_I_P
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_weights_object
