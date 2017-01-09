!< Abstract WENO optimal weights object.
module wenoof_optimal_weights_abstract
!< Abstract WENO optimal weights object.

use penf, only : I_P, R_P

implicit none
private
public :: optimal_weights

type, abstract :: optimal_weights
  !< WENO optimal weights.
  !<
  !< @note Do not implement any real optimal weight: provide the interface for the different optimal weights implemented.
  real(R_P), allocatable :: opt(:,:) !< Optimal weights [1:2,0:S-1].
  contains
    procedure(create_interface),      pass(self), deferred :: create      !< Create weights.
    procedure(description_interface), nopass,     deferred :: description !< Return string-description of weights.
    procedure(destroy_interface),     pass(self), deferred :: destroy     !< Destroy weights.
endtype optimal_weights

abstract interface
  !< Create weights.
  pure subroutine create_interface(self, S)
  !< Create weights.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  import :: optimal_weights, I_P
  class(optimal_weights), intent(inout) :: self !< WENO optimal weights.
  integer(I_P),           intent(in)    :: S    !< Number of stencils used.
  endsubroutine create_interface
endinterface

abstract interface
  !< Return string-description of weights.
  pure subroutine description_interface(string)
  !< Return string-description of weights.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  endsubroutine description_interface
endinterface

abstract interface
  !< Destroy weights.
  pure subroutine destroy_interface(self)
  !< Destroy weights.
  import :: optimal_weights
  class(optimal_weights), intent(inout) :: self !< WENO optimal weights.
  endsubroutine destroy_interface
endinterface
endmodule wenoof_optimal_weights_abstract
