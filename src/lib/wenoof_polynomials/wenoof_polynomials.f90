module weno_polynomials
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
  real(R_P), allocatable :: poly(:,:)     !< Polynomial reconstructions [1:2,0:S-1].
  real(R_P), allocatable :: coef(:,:,:)   !< Polynomial coefficients [1:2,0:S-1,0:S-1].
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
    procedure(description_interface), nopass,     deferred, public :: description
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
  pure subroutine description_interface(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO polynomials.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  !< Compute the partial value of the WENO interpolating polynomial.
  pure subroutine compute_interface(self, S, stencil, f1, f2, ff)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the partial value of the interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_polynomials, I_P, R_P
  class(weno_polynomials), intent(inout) :: self                    !< WENO polynomial.
  integer(I_P),            intent(in)    :: S                       !< Number of stencils actually used.
  real(R_P),               intent(in)    :: stencil(1:, 1 - S:)     !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P),            intent(in)    :: f1, f2, ff              !< Faces to be computed.
  integer(I_P)                           :: s1, s2, f               !< Counters
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule weno_polynomials
