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
  integer(I_P) :: S=0_I_P !< Stencils dimension.
endtype interpolations_object_constructor

type, extends(base_object) :: interpolations_object
  !< Abstract interpolations object.
  integer(I_P)           :: S=0_I_P     !< Stencils dimension.
  real(R_P), allocatable :: values(:,:) !< Stencil interpolations values [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute interpolations.
    procedure, nopass     :: description !< Return interpolations string-description.
    ! public methods
    procedure, pass(self) :: create  !< Createte interpolations.
    procedure, pass(self) :: destroy !< Destroy interpolations.
endtype interpolations_object

contains
  ! deferred public methods
  pure subroutine compute(self, stencil)
  !< Compute interpolations.
  class(interpolation_objects), intent(inout) :: self                  !< Interpolations.
  real(R_P),                    intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolations%compute to be implemented by your concrete interpolations object'
#endif
  endsubroutine compute

  pure function description(self) result(string)
  !< Return interpolations string-description.
  class(interpolation_objects), intent(in) :: self   !< Interpolations.
  character(len=:), allocatable            :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolations%description to be implemented by your concrete interpolations object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create interpolations.
  class(interpolations_object),   intent(inout) :: self        !< Interpolations.
  class(base_object_constructor), intent(in)    :: constructor !< Interpolations constructor.

  call self%destroy
  select type(constructor)
  class is(interpolations_constructor)
    self%S = constructor%S
    allocate(self%poly(1:2, 0:constructor%S - 1))
  class default
    ! @TODO add error handling
  endselect
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy interpolations.
  class(interpolation_objects), intent(inout) :: self !< Interpolations.

  self%S = 0_I_P
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_interpolations_object
