module type_weno_polynomials
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
type, abstract :: weno_polynomials_constructor
  !< Abstract type used for create new concrete WENO polynomials.
  !<
  !< @note Every concrete WENO polynomials implementation must define its own constructor type.
  private
endtype weno_polynomials_constructor

type, abstract :: weno_polynomials
  !< WENO polynomials.
  !<
  !< @note Do not implement any real polynomial: provide the interface for the different polynomials implemented.
  private
  contains
    procedure(destructor_interface),  pass(self), deferred, public :: destructor
    procedure(constructor_interface), pass(self), deferred, public :: constructor
    procedure(description_interface), pass(self), deferred, public :: description
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
  pure subroutine constructor_interface(self,constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO polynomials.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_polynomials_constructor
  import :: weno_polynomials
  class(weno_polynomials),             intent(inout) :: self          !< WENO polynomials.
  class(weno_polynomials_constructor), intent(inout) :: constructor   !< WENO polynomials constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine constructor_interface
endinterface

abstract interface
  !< Return a string describing WENO polynomials.
  pure subroutine description_interface(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO polynomials.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_polynomials
  class(weno_polynomials),       intent(in)  :: self   !< WENO polynomials.
  character(len=:), allocatable, intent(out) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description_interface
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_polynomials
