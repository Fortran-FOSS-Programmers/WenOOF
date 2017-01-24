!< Abstract alpha (non linear weights) object.
module wenoof_alpha_object
!< Abstract alpha (non linear weights) object.

use penf, only : I_P, R_P
use wenoof_base_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_object
public :: alpha_object_constructor

type, extends(base_object_constructor) :: alpha_object_constructor
  !< Abstract alpha (non linear weights) object constructor.
  contains
endtype alpha_object_constructor

type, extends(base_object), abstract :: alpha_object
  !< Abstract alpha (non linear weights) object.
  real(R_P), allocatable :: values(:,:)   !< Alpha coefficients [1:2,0:S-1].
  real(R_P), allocatable :: values_sum(:) !< Sum of alpha coefficients [1:2].
  contains
    ! public deferred methods
    procedure(compute_interface), pass(self), deferred :: compute !< Compute alpha.
    ! public overridable methods
    procedure, pass(self) :: create  !< Create alpha.
    procedure, pass(self) :: destroy !< Destroy alpha.
endtype alpha_object

abstract interface
  !< Abstract interfaces of [[alpha_object]].
  ! subroutine create_interface(self, constructor)
  ! !< Create alpha.
  ! import :: alpha_object, base_object_constructor
  ! class(alpha_object),            intent(inout) :: self        !< Alpha.
  ! class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.
  ! endsubroutine create_interface

  pure subroutine compute_interface(self, beta, kappa)
  !< Compute alpha.
  import :: alpha_object, beta_object, kappa_object
  class(alpha_object), intent(inout) :: self  !< Alpha.
  class(beta_object),  intent(in)    :: beta  !< Beta.
  class(kappa_object), intent(in)    :: kappa !< Kappa.
  endsubroutine compute_interface

  ! pure function description_interface(self) result(string)
  ! !< Return alpha string-description.
  ! import :: alpha_object
  ! class(alpha_object), intent(in) :: self   !< Alpha.
  ! character(len=:), allocatable   :: string !< String-description.
  ! endfunction description_interface

  ! elemental subroutine destroy_interface(self)
  ! !< Destroy alpha.
  ! import :: alpha_object
  ! class(alpha_object), intent(inout) :: self !< Alpha.
  ! endsubroutine destroy_interface
endinterface

contains
  ! public overridable methods
  subroutine create(self, constructor)
  !< Creat alpha.
  class(alpha_object),            intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  allocate(self%values_sum(1:2))
  self%values = 0._R_P
  self%values_sum = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy alpha.
  class(alpha_object), intent(inout) :: self !< Alpha.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  if (allocated(self%values_sum)) deallocate(self%values_sum)
  endsubroutine destroy
endmodule wenoof_alpha_object
