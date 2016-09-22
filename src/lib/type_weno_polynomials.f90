module type_weno_polynomials
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO polynomials object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_polynomials
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_polynomials
  !< WENO polynomials.
  !<
  !< @note Do not implement any real polynomial: provide the interface for the different polynomials implemented.
  private
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
    procedure(description_interface), pass(self), deferred, public :: description
    procedure(compute_interface),     pass(self), deferred, public :: compute
endtype weno_polynomials

abstract interface
  !< Destroy WENO polynomials.
  pure subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO polynomials.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_polynomials
  class(weno_polynomials), intent(inout) :: self   !< WENO polynomials.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  !< Create WENO polynomials.
  pure subroutine constructor_interface(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO polynomials.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_polynomials, I_P
  class(weno_polynomials),             intent(inout) :: self          !< WENO polynomials.
  integer(I_P),                        intent(in)    :: S             !< Stencil dimension.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Return a string describing WENO polynomials.
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO polynomials.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_polynomials
  class(weno_polynomials),       intent(in)  :: self   !< WENO polynomials.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  !< Compute the partial value of the WENO interpolating polynomial.
  pure function compute_interface(self, poly_coef, v) result(poly)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the partial value of the interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_polynomials, I_P, R_P
  class(weno_polynomials), intent(in) :: self        !< WENO polynomial.
  real(R_P),               intent(in) :: poly_coef   !< Polynomila coefficient for the value v.
  real(R_P),               intent(in) :: v           !< Single value of the interpolation stencil.
  real(R_P)                           :: poly        !< Partial value of the interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_polynomials
