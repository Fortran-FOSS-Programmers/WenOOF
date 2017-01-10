!< Abstract alpha coefficients object.
module wenoof_alpha_coefficients
!< Abstract alpha coefficients object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: alpha_coefficients
public :: alpha_coefficients_constructor

type, extends(base_object_constructor) :: alpha_coefficients_constructor
  !< Abstract alpha coefficients object constructor.
  integer(I_P) :: S = 0 !< Stencils dimension.
endtype alpha_coefficients_constructor

type, extends(base_object) :: alpha_coefficients
  !< Abstract alpha coefficients object.
  !<
  !< @note Do not implement any real alpha coefficient, but provide the interface for the different alpha coefficient implemented.
  real(R_P), allocatable :: alpha_coef(:,:) !< Alpha coefficients [1:2,0:S-1]
  real(R_P), allocatable :: alpha_tot(:)    !< Sum of alpha coefficients
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute alpha coefficients.
    procedure, nopass     :: description !< Return alpha coefficients string-description.
    ! public methods
    procedure, pass(self) :: create  !< Create alpha coefficients.
    procedure, pass(self) :: destroy !< Destroy alpha coefficients.
endtype alpha_coefficients

contains
  ! deferred public methods
  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute alpha coefficients.
  class(alpha_coefficients), intent(inout) :: self                   !< Alpha coefficients.
  integer(I_P),              intent(in)    :: S                      !< Number of stencils used.
  real(R_P),                 intent(in)    :: weight_opt(1:2, 0:S-1) !< Optimal weight of the stencil.
  real(R_P),                 intent(in)    :: IS(1:2, 0:S-1)         !< Smoothness indicators of the stencils.
  real(R_P),                 intent(in)    :: eps                    !< Parameter for avoiding divided by zero.
  integer(I_P),              intent(in)    :: f1, f2                 !< Faces to be computed.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'alpha_coefficients%compute to be implemented by your concrete alpha coefficients object'
#endif
  endsubroutine compute

  pure subroutine description(string)
  !< Return alpha coefficients string-description.
  character(len=:), allocatable  :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'alpha_coefficients%description to be implemented by your concrete alpha coefficients object'
#endif
  endsubroutine description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create alpha coefficients.
  class(alpha_coefficients),             intent(inout) :: self        !< Alpha coefficients.
  class(alpha_coefficients_constructor), intent(in)    :: constructor !< Alpha coefficients constructor.

  call self%destroy
  allocate(self%alpha_coef(1:2, 0:constructor%S - 1))
  allocate(self%alpha_tot(1:2))
  self%alpha_coef(:,:) = 0._R_P
  self%alpha_tot(:) = 0._R_P
  endsubroutine create

  pure subroutine destroy(self)
  !< Destroy alpha coefficients.
  class(alpha_coefficients), intent(inout) :: self !< Alpha coefficients.

  if (allocated(self%alpha_coef)) deallocate(self%alpha_coef)
  if (allocated(self%alpha_tot)) deallocate(self%alpha_tot)
  endsubroutine destroy
endmodule wenoof_alpha_coefficients_abstract
