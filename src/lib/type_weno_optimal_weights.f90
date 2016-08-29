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
    procedure(optimal_weights_abstract_destructor),  pass(self), deferred, public :: optimal_weights_destructor
    procedure(optimal_weights_abstract_constructor), pass(self), deferred, public :: optimal_weights_constructor
    procedure(optimal_weights_abstract_description), pass(self), deferred, public :: optimal_weights_description
endtype weno_optimal_weights

abstract interface

  pure subroutine optimal_weights_abstract_destructor(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_optimal_weights
  class(weno_optimal_weights),   intent(INOUT)  :: self   !< WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine optimal_weights_abstract_destructor

  pure subroutine optimal_weights_abstract_constructor(self, constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO optimal weights.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_optimal_weights
  class(weno_optimal_weights),             intent(INOUT)  :: self          !< WENO optimal weights.
  class(weno_optimal_weights_constructor), intent(INOUT)  :: constructor   !< WENO optimal weights constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine optimal_weights_abstract_constructor

  pure subroutine optimal_weights_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_optimal_weights
  class(weno_optimal_weights),   intent(IN)  :: self   !< WENO optimal weights.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine optimal_weights_abstract_description

endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_optimal_weights
