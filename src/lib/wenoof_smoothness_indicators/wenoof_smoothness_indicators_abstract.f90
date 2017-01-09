module wenoof_smoothness_indicators_abstract
!< Abstract WENO smoothness indicators object.

use penf, only : I_P, R_P

implicit none
private
save
public :: IS

type, abstract :: IS
  !< WENO smoothness indicators.
  !<
  !< @note Do not implement any real smoothness indicators: provide the interface for the different
  !< smoothness_indicators implemented.
  real(R_P), allocatable :: si(:,:)     !< Smoothness indicators [1:2,0:S-1].
  real(R_P), allocatable :: coef(:,:,:) !< Smoothness indicators coefficients [1:2,0:S-1,0:S-1].
  contains
    procedure(compute_interface),     pass(self), deferred :: compute     !< Compute IS.
    procedure(create_interface),      pass(self), deferred :: create      !< Create IS.
    procedure(description_interface), nopass,     deferred :: description !< Return string-description of IS.
    procedure(destroy_interface),     pass(self), deferred :: destroy     !< Destroy IS.
endtype IS

abstract interface
  !< Compute IS.
  pure subroutine compute_interface(self, S, stencil, f1, f2, ff)
  !< Compute IS.
  import :: IS, I_P, R_P
  class(IS),    intent(inout) :: self                !< WENO smoothness indicator.
  integer(I_P), intent(in)    :: S                   !< Number of stencils actually used.
  real(R_P),    intent(in)    :: stencil(1:, 1 - S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P), intent(in)    :: f1, f2, ff          !< Faces to be computed.
  integer(I_P)                :: s1, s2, s3, f       !< Counters
  endsubroutine compute_interface
endinterface

abstract interface
  !< Create IS.
  pure subroutine create_interface(self, S)
  !< Create IS.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  import :: IS, I_P
  class(IS),    intent(inout) :: self !< WENO smoothness indicators.
  integer(I_P), intent(in)    :: S    !< Number of stencils used.
  endsubroutine create_interface
endinterface

abstract interface
  !< Return string-description of IS.
  pure subroutine description_interface(string)
  !< Return string-description of IS.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  endsubroutine description_interface
endinterface

abstract interface
  !< Destroy IS.
  pure subroutine destroy_interface(self)
  !< Destroy IS.
  import :: IS
  class(IS), intent(inout) :: self !< WENO smoothenss indicators.
  endsubroutine destroy_interface
endinterface
endmodule wenoof_smoothness_indicators_abstract
