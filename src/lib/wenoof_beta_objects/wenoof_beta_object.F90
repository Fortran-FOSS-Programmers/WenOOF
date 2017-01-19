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

type, extends(base_object) :: beta_object
  !< Abstract Beta coefficients (smoothness indicators of stencil interpolations) object.
  real(R_P), allocatable :: values(:,:) !< Beta values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: compute     !< Compute beta.
    procedure, pass(self) :: description !< Return beta string-description.
    ! public methods
    procedure, pass(self) :: create  !< Createte beta.
    procedure, pass(self) :: destroy !< Destroy beta.
endtype beta_object

contains
  ! public deferred methods
  pure subroutine compute(self, stencil)
  !< Compute beta.
  class(beta_object), intent(inout) :: self                  !< Beta.
  real(R_P),          intent(in)    :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'beta_object%compute to be implemented by your concrete beta object'
#endif
  endsubroutine compute

  pure function description(self) result(string)
  !< Return beta string-description.
  class(beta_object), intent(in) :: self   !< Beta.
  character(len=:), allocatable  :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'beta_object%description to be implemented by your concrete beta object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create beta.
  class(beta_object),             intent(inout) :: self        !< beta.
  class(base_object_constructor), intent(in)    :: constructor !< beta constructor.

  call self%destroy
  call self%base_object%create(constructor=constructor)
  allocate(self%values(1:2, 0:self%S-1))
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy beta.
  class(beta_object), intent(inout) :: self !< beta.

  call self%base_object%destroy
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_beta_object
