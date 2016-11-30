module wenoof_smoothness_indicators_abstract
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO smoothness indicators object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: IS
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: IS
  !< WENO smoothness indicators.
  !<
  !< @note Do not implement any real smoothness indicators: provide the interface for the different smoothness_indicators implemented.
  real(R_P), allocatable :: si(:,:)       !< Smoothness indicators [1:2,0:S-1].
  real(R_P), allocatable :: coef(:,:,:)   !< Smoothness indicators coefficients [1:2,0:S-1,0:S-1].
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
    procedure(description_interface), nopass,     deferred, public :: description
    procedure(compute_interface),     pass(self), deferred, public :: compute
endtype IS

abstract interface
  !< Destroy WENO polynomial coefficients.
  pure subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: IS
  class(IS), intent(inout) :: self   !< WENO smoothenss indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  !< Create WENO polynomial coefficients.
  pure subroutine constructor_interface(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO smoothness indicators coefficients.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: IS, I_P
  class(IS),    intent(inout) :: self        !< WENO smoothness indicators.
  integer(I_P), intent(in)    :: S           !< Number of stencils used.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Return a string describing WENO smoothness_indicators.
  pure subroutine description_interface(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO smoothness_indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface

abstract interface
  !< Compute the smoothness indicators of the WENO interpolating polynomial.
  pure subroutine compute_interface(self, S, stencil, f1, f2, ff)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the partial value of the smoothness indicator of a single WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: IS, I_P, R_P
  class(IS),    intent(inout) :: self                    !< WENO smoothness indicator.
  integer(I_P), intent(in)    :: S                       !< Number of stencils actually used.
  real(R_P),    intent(in)    :: stencil(1:, 1 - S:)     !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P), intent(in)    :: f1, f2, ff              !< Faces to be computed.
  integer(I_P)                :: s1, s2, s3, f           !< Counters
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule wenoof_smoothness_indicators_abstract
