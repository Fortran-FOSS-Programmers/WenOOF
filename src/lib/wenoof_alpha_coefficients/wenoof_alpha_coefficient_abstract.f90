!< Abstract WENO alpha coefficient object.
module wenoof_alpha_coefficient_abstract
!< Abstract WENO alpha coefficient object.
!<
!< @TODO Add equations.

use penf, only : I_P, R_P

implicit none
private
public :: alpha_coefficient

type :: alpha_coefficient
  !< WENO weights.
  !<
  !< @note Do not implement any real alpha coefficient, but provide the interface for the different alpha coefficient implemented.
  real(R_P), allocatable :: alpha_coef(:,:) !< Alpha coefficients [1:2,0:S-1]
  real(R_P), allocatable :: alpha_tot(:)    !< Sum of alpha coefficients
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute coefficients.
    procedure, nopass     :: description !< Return string-description of coefficients.
    ! public methods
    procedure, pass(self) :: alloc   !< Allocate coefficients.
    procedure, pass(self) :: destroy !< Destroy coefficients.
endtype alpha_coefficient

contains
  ! deferred public methods
  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute coefficients.
  class(alpha_coefficient), intent(inout) :: self                   !< WENO alpha coefficient.
  integer(I_P),             intent(in)    :: S                      !< Number of stencils used.
  real(R_P),                intent(in)    :: weight_opt(1:2, 0:S-1) !< Optimal weight of the stencil.
  real(R_P),                intent(in)    :: IS(1:2, 0:S-1)         !< Smoothness indicators of the stencils.
  real(R_P),                intent(in)    :: eps                    !< Parameter for avoiding divided by zero.
  integer(I_P),             intent(in)    :: f1, f2                 !< Faces to be computed.

  error stop 'alpha_coefficient%compute to be implemented by your concrete alpha coefficient object'
  endsubroutine compute

  pure subroutine description(string)
  !< Return string-description of coefficients.
  !<
  !< @TODO make a function.
  character(len=:), allocatable, intent(out) :: string !< String returned.

  error stop 'alpha_coefficient%description to be implemented by your concrete alpha coefficient object'
  endsubroutine description

  ! public methods
  pure subroutine alloc(self, S)
  !< Alloc coefficients.
  class(alpha_coefficient), intent(inout) :: self !< WENO alpha coefficients.
  integer(I_P),             intent(in)    :: S    !< Number of stencils used.

  call self%destroy
  allocate(self%alpha_coef(1:2, 0:S - 1))
  allocate(self%alpha_tot(1:2))
  self%alpha_coef(:,:) = 0._R_P
  self%alpha_tot(:) = 0._R_P
  endsubroutine alloc

  pure subroutine destroy(self)
  !< Destroy coefficients.
  class(alpha_coefficient), intent(inout) :: self !< WENO alpha coefficients.

  if (allocated(self%alpha_coef)) deallocate(self%alpha_coef)
  if (allocated(self%alpha_tot)) deallocate(self%alpha_tot)
  endsubroutine destroy
endmodule wenoof_alpha_coefficient_abstract
