module type_weno_alpha_coefficient
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO alpha coefficient object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_alpha_coefficient
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_alpha_coefficient
  !< WENO weights.
  !<
  !< @note Do not implement any real alpha coefficient provide the interface for the different alpha coefficient implemented.
  real(R_P), allocatable :: alpha(:,:)   !< Alpha coefficients [1:2,0:S-1]
  real(R_P), allocatable :: alpha_tot(:) !< Sum of alpha coefficients
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
    procedure(description_interface), nopass,     deferred, public :: description
    procedure(compute_interface),     pass(self), deferred, public :: compute
endtype weno_alpha_coefficient

abstract interface
  !< Destroy WENO alpha coefficients.
  pure subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO alpha coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_alpha_coefficient
  class(weno_alpha_coefficient), intent(inout) :: self   !< WENO alpha coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  !< Create WENO alpha coefficients.
  pure subroutine constructor_interface(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO alpha coefficients.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_alpha_coefficient, I_P
  class(weno_alpha_coefficient), intent(inout) :: self        !< WENO alpha coefficients.
  integer(I_P),                  intent(in)    :: S           !< Number of stencils used.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Return a string describing WENO alpha coefficient.
  pure subroutine description_interface(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  pure subroutine compute_interface(self, S, weight_opt, IS, eps, f1, f2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_alpha_coefficient, I_P, R_P
  class(weno_alpha_coefficient), intent(inout) :: self                         !< WENO alpha coefficient.
  integer(I_P),                  intent(in)    :: S                            !< Number of stencils used.
  real(R_P),                     intent(in)    :: weight_opt(1: 2, 0: S - 1)   !< Optimal weight of the stencil.
  real(R_P),                     intent(in)    :: IS(1: 2, 0: S - 1)           !< Smoothness indicators of the stencils.
  real(R_P),                     intent(in)    :: eps                          !< Parameter for avoiding divided by zero.
  integer(I_P),                  intent(in)    :: f1, f2                       !< Faces to be computed.
  integer(I_P)                                 :: f, s1                        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_alpha_coefficient
