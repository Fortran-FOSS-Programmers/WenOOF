!< Abstract Kappa coefficients (optimal, linear weights of stencil interpolations) object.
module wenoof_kappa_object
!< Abstract Kappa coefficients (optimal, linear weights of stencil interpolations) object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: kappa_object
public :: kappa_object_constructor

type, extends(base_object_constructor) :: kappa_object_constructor
  !< Abstract kappa object constructor.
endtype kappa_object_constructor

type, extends(base_object) :: kappa_object
  !< Kappa coefficients (optimal, linear weights of stencil interpolations) object.
  real(R_P), allocatable :: values(:,:) !< Kappa coefficients values [1:2,0:S-1].
  contains
    ! public deferred methods
    procedure, pass(self) :: compute     !< Compute kappa.
    procedure, pass(self) :: description !< Return kappa string-description.
    ! public methods
    procedure, pass(self) :: create  !< Createte kappa.
    procedure, pass(self) :: destroy !< Destroy kappa.
endtype kappa_object

contains
  ! public deferred methods
  pure subroutine compute(self)
  !< Compute kappa.
  class(kappa_object), intent(inout) :: self !< Kappa.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'kappa%compute to be implemented by your concrete kappa object'
#endif
  endsubroutine compute

  pure function description(self) result(string)
  !< Return kappa string-description.
  class(kappa_object), intent(in) :: self   !< Kappa.
  character(len=:), allocatable   :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'kappa%description to be implemented by your concrete kappa object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create kappa.
  class(kappa_object),            intent(inout) :: self        !< Kappa.
  class(base_object_constructor), intent(in)    :: constructor !< Kappa constructor.

  call self%destroy
  call self%base_object%create(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  self%values = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy kappa.
  class(kappa_object), intent(inout) :: self !< Kappa.

  call self%base_object%destroy
  if (allocated(self%values)) deallocate(self%values)
  endsubroutine destroy
endmodule wenoof_kappa_object
