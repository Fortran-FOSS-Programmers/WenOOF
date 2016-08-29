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
    procedure(IS_abstract_destructor),  pass(self), deferred, public :: IS_destructor
    procedure(IS_abstract_constructor), pass(self), deferred, public :: IS_constructor
    procedure(IS_abstract_description), pass(self), deferred, public :: IS_description
endtype weno_IS

abstract interface

  pure subroutine IS_abstract_destructor(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS
  class(weno_IS), intent(INOUT) :: self   !< WENO smoothenss indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IS_abstract_destructor

  pure subroutine IS_abstract_constructor(self,constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO polynomial coefficients.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS_constructor
  import :: weno_IS
  class(weno_IS),             intent(INOUT)  :: self          !< WENO smoothness indicators.
  class(weno_IS_constructor), intent(INOUT)  :: constructor   !< WENO smoothness indicators constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IS_abstract_constructor

  pure subroutine IS_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO smoothness_indicators.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_IS
  class(weno_IS),                intent(IN)  :: self   !< WENO smoothness indicator.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine IS_abstract_description

endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_smoothness_indicators
