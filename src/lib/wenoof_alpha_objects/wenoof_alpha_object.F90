!< Abstract alpha coefficients (non linear weights) object.
module wenoof_alpha_object
!< Abstract alpha coefficients (non linear weights) object.

use penf, only : I_P, R_P
use wenoof_base_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_object
public :: alpha_object_constructor

type, extends(base_object_constructor) :: alpha_object_constructor
  !< Abstract alpha coefficients (non linear weights) object constructor.
endtype alpha_object_constructor

type, extends(base_object) :: alpha_object
  !< Abstract alpha coefficients (non linear weights) object.
  real(R_P), allocatable :: values(:,:)   !< Alpha coefficients [1:2,0:S-1].
  real(R_P), allocatable :: values_sum(:) !< Sum of alpha coefficients [1:2].
  contains
    ! public deferred methods
    procedure, pass(self) :: compute     !< Compute alpha coefficients.
    procedure, nopass     :: description !< Return alpha coefficients string-description.
    ! public methods
    procedure, pass(self) :: create  !< Create alpha coefficients.
    procedure, pass(self) :: destroy !< Destroy alpha coefficients.
endtype alpha_object

contains
  ! public deferred methods
  pure subroutine compute(self, beta, kappa)
  !< Compute alpha coefficients.
  class(alpha_object), intent(inout) :: self  !< Alpha coefficients.
  class(beta_object),  intent(in)    :: beta  !< Beta coefficients.
  class(kappa_object), intent(in)    :: kappa !< Kappa coefficients.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'alpha_object%compute to be implemented by your concrete alpha coefficients object'
#endif
  endsubroutine compute

  pure function description(self) result(string)
  !< Return alpha coefficients string-description.
  class(alpha_object), intent(in) :: self   !< Alpha coefficients.
  character(len=:), allocatable   :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'alpha_object%description to be implemented by your concrete alpha coefficients object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create alpha coefficients.
  class(alpha_object),            intent(inout) :: self        !< Alpha coefficients.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha coefficients constructor.

  call self%destroy
  call self%base_object%create(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  allocate(self%values_sum(1:2))
  self%values = 0._R_P
  self%values_sum = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy alpha coefficients.
  class(alpha_object), intent(inout) :: self !< Alpha coefficients.

  call self%base_object%destroy
  if (allocated(self%values)) deallocate(self%values)
  if (allocated(self%values_sum)) deallocate(self%values_sum)
  endsubroutine destroy
endmodule wenoof_alpha_object
