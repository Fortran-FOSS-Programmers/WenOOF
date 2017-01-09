!< Abstract WENO alpha coefficient object.
module wenoof_alpha_coefficient_abstract
!< Abstract WENO alpha coefficient object.
!<
!< @TODO Add equations.

use penf, only : I_P, R_P

implicit none
private
save
public :: alpha_coefficient

type, abstract :: alpha_coefficient
  !< WENO weights.
  !<
  !< @note Do not implement any real alpha coefficient, but provide the interface for the different alpha coefficient implemented.
  real(R_P), allocatable :: alpha_coef(:,:) !< Alpha coefficients [1:2,0:S-1]
  real(R_P), allocatable :: alpha_tot(:)    !< Sum of alpha coefficients
  contains
    procedure(constructor_interface), pass(self), deferred :: create      !< Create coefficients.
    procedure(compute_interface),     pass(self), deferred :: compute     !< Compute coefficients.
    procedure(description_interface), nopass,     deferred :: description !< Return string-description of coefficients.
    procedure(destructor_interface),  pass(self), deferred :: destroy     !< Destroy coefficients.
endtype alpha_coefficient

abstract interface
  !< Create coefficients.
  pure subroutine constructor_interface(self, S)
  !< Create coefficients.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  import :: alpha_coefficient, I_P
  class(alpha_coefficient), intent(inout) :: self !< WENO alpha coefficients.
  integer(I_P),             intent(in)    :: S    !< Number of stencils used.
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Compute coefficients.
  pure subroutine compute_interface(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute coefficients.
  import :: alpha_coefficient, I_P, R_P
  class(alpha_coefficient), intent(inout) :: self                   !< WENO alpha coefficient.
  integer(I_P),             intent(in)    :: S                      !< Number of stencils used.
  real(R_P),                intent(in)    :: weight_opt(1:2, 0:S-1) !< Optimal weight of the stencil.
  real(R_P),                intent(in)    :: IS(1:2, 0:S-1)         !< Smoothness indicators of the stencils.
  real(R_P),                intent(in)    :: eps                    !< Parameter for avoiding divided by zero.
  integer(I_P),             intent(in)    :: f1, f2                 !< Faces to be computed.
  endsubroutine compute_interface
endinterface

abstract interface
  !< Return string-description of coefficients.
  pure subroutine description_interface(string)
  !< Return string-description of coefficients.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  endsubroutine description_interface
endinterface

abstract interface
  !< Destroy coefficients.
  pure subroutine destructor_interface(self)
  !< Destroy coefficients.
  import :: alpha_coefficient
  class(alpha_coefficient), intent(inout) :: self   !< WENO alpha coefficients.
  endsubroutine destructor_interface
endinterface
endmodule wenoof_alpha_coefficient_abstract
