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
type, abstract :: weno_optimal_weights
  !< WENO optimal weights.
  !<
  !< @note Do not implement any real optimal weight: provide the interface for the different optimal weights implemented.
  real(R_P), allocatable :: opt(:,:)   !< Optimal weights                    [1:2,0:S-1].
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destroy
    procedure(constructor_interface), pass(self), deferred, public :: create
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
  pure subroutine constructor_interface(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO optimal weights.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_optimal_weights, I_P
  class(weno_optimal_weights), intent(inout) :: self       !< WENO optimal weights.
  integer(I_P),                intent(in)    :: S          !< Number of stencils used.
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
