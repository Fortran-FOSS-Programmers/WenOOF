module type_weno_IS
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
public :: weno_IS
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_IS_constructor
  !< Abstract type used for create new concrete WENO smoothness indicators.
  !<
  !< @note Every concrete WENO smoothness indicators implementation must define its own constructor type.
  private
endtype weno_IS_constructor

type, abstract :: weno_IS
  !< WENO smoothness indicators.
  !<
  !< @note Do not implement any real smoothness indicators: provide the interface for the different smoothness_indicators implemented.
  private
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destructor
    procedure(constructor_interface), pass(self), deferred, public :: constructor
    procedure(description_interface), pass(self), deferred, public :: description
endtype weno_IS

abstract interface
  !< Destroy WENO polynomial coefficients.
  pure subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS
  class(weno_IS), intent(inout) :: self   !< WENO smoothenss indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  !< Create WENO polynomial coefficients.
  pure subroutine constructor_interface(self,constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO polynomial coefficients.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS_constructor
  import :: weno_IS
  class(weno_IS),             intent(inout)  :: self          !< WENO smoothness indicators.
  class(weno_IS_constructor), intent(inout)  :: constructor   !< WENO smoothness indicators constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Return a string describing WENO smoothness_indicators.
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO smoothness_indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS
  class(weno_IS),                intent(in)  :: self   !< WENO smoothness indicator.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_smoothness_indicators
