!< Abstract WENO polynomials object.
module wenoof_polynomials_abstract
!< Abstract WENO polynomials object.

use penf, only : I_P, R_P

implicit none
private
public :: polynomials

type, abstract :: polynomials
  !< Abstract WENO polynomials object.
  !<
  !< @note Do not implement any real polynomial: provide the interface for the different polynomials implemented.
  real(R_P), allocatable :: poly(:,:)   !< Polynomial reconstructions [1:2,0:S-1].
  real(R_P), allocatable :: coef(:,:,:) !< Polynomial coefficients [1:2,0:S-1,0:S-1].
  contains
    procedure(compute_interface),     pass(self), deferred :: compute     !< Compute polynomials.
    procedure(create_interface),      pass(self), deferred :: create      !< Create polynomials.
    procedure(description_interface), nopass,     deferred :: description !< Return string-description of polynomials.
    procedure(destroy_interface),     pass(self), deferred :: destroy     !< Destroy polynomials.
endtype polynomials

abstract interface
  !< Compute polynomials.
  pure subroutine compute_interface(self, S, stencil, f1, f2, ff)
  !< Compute polynomials.
  import :: polynomials, I_P, R_P
  class(polynomials), intent(inout) :: self                !< WENO polynomial.
  integer(I_P),       intent(in)    :: S                   !< Number of stencils actually used.
  real(R_P),          intent(in)    :: stencil(1:, 1 - S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P),       intent(in)    :: f1, f2, ff          !< Faces to be computed.
  integer(I_P)                      :: s1, s2, f           !< Counters
  endsubroutine compute_interface
endinterface

abstract interface
  !< Create polynomials.
  pure subroutine create_interface(self, S)
  !< Create polynomials.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  import :: polynomials, I_P
  class(polynomials), intent(inout) :: self !< WENO polynomials.
  integer(I_P),       intent(in)    :: S    !< Stencil dimension.
  endsubroutine create_interface
endinterface

abstract interface
  !< Return string-description of polynomials.
  pure subroutine description_interface(string)
  !< Return string-description of polynomials.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  endsubroutine description_interface
endinterface

abstract interface
  !< Destroy polynomials.
  pure subroutine destroy_interface(self)
  !< Destroy polynomials.
  import :: polynomials
  class(polynomials), intent(inout) :: self !< WENO polynomials.
  endsubroutine destroy_interface
endinterface
endmodule wenoof_polynomials_abstract
