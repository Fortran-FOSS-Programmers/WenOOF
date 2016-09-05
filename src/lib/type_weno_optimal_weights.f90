module type_weno_optimal_weights
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO optimal weights object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_optimal_weights
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_optimal_weights_constructor
  !< Abstract type used for create new concrete WENO optimal weights.
  !<
  !< @note Every concrete WENO optimale weights implementation must define its own constructor type.
  private
endtype weno_optimal_weights_constructor

type, abstract :: weno_optimal_weights
  !< WENO optimal weights.
  !<
  !< @note Do not implement any real optimal weight: provide the interface for the different optimal weights implemented.
  private
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destructor
    procedure(constructor_interface), pass(self), deferred, public :: constructor
    procedure(description_interface), pass(self), deferred, public :: description
endtype weno_optimal_weights

abstract interface
  !< Destroy WENO optimal weights.
  pure subroutine destructor_interface(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_optimal_weights
  class(weno_optimal_weights), intent(inout)  :: self   !< WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destructor_interface
endinterface

abstract interface
  !< Create WENO optimal weights.
  pure subroutine constructor_interface(self, constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO optimal weights.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_optimal_weights
  class(weno_optimal_weights),             intent(inout)  :: self          !< WENO optimal weights.
  class(weno_optimal_weights_constructor), intent(inout)  :: constructor   !< WENO optimal weights constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Return a string describing WENO optimal weights.
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_optimal_weights
  class(weno_optimal_weights),   intent(in)  :: self   !< WENO optimal weights.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_optimal_weights
