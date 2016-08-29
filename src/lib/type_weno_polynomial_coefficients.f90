module type_weno_poly_coefficients
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO polynomial coefficients object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_poly_coefficients
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_poly_coefficients_constructor
  !< Abstract type used for create new concrete WENO polynomial coefficients.
  !<
  !< @note Every concrete WENO polynomial coefficients implementation must define its own constructor type.
  private
endtype weno_poly_coefficients_constructor

type, abstract :: weno_poly_coefficients
  !< WENO polynomial coefficients.
  !<
  !< @note Do not implement any real polynomial coefficient: provide the interface for the different polynomial coefficients implemented.
  private
  contains
    procedure(poly_coefficients_abstract_destructor),  pass(self), deferred, public :: poly_coefficients_destructor
    procedure(poly_coefficients_abstract_constructor), pass(self), deferred, public :: poly_coefficients_constructor
    procedure(poly_coefficients_abstract_description), pass(self), deferred, public :: poly_coefficients_description
endtype weno_poly_coefficients

abstract interface

  pure subroutine poly_coefficients_abstract_destructor(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_poly_coefficients
  class(weno_poly_coefficients), intent(INOUT)  :: self   !< WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine poly_coefficients_abstract_destructor

  pure subroutine poly_coefficients_abstract_constructor(self,constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO polynomial coefficients.
  !
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_poly_coefficients_constructor
  import :: weno_poly_coefficients
  class(weno_poly_coefficients),             intent(INOUT)  :: self          !< WENO polynomial coefficients.
  class(weno_poly_coefficients_constructor), intent(INOUT)  :: constructor   !< WENO polynomial coefficients constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine poly_coefficients_abstract_constructor

  pure subroutine poly_coefficients_abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO polynomial coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_poly_coefficients
  class(weno_poly_coefficients), intent(IN)  :: self   !< WENO polynomial coefficients.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine poly_coefficients_abstract_description

endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_poly_coefficients
