module type_weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO interpolator object,
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_interpolator, weno_constructor
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract :: weno_interpolator
  !< WENO interpolator factory object.
  !<
  !< @note Do not implement any real interpolator: provide the interface for the different interpolators implemented.
  !<
  !< @note A concrete implementation must define its own **interpolate** method: it is not provide as abstract deferred type
  !< because it could have very diffirent signature for each different interpolator.
  private
  contains
    procedure(abstract_destructor),  pass(self), deferred :: destroy
    procedure(abstract_constructor), pass(self), deferred :: create
    procedure(abstract_description), pass(self), deferred :: description
endtype weno_interpolator

type, abstract :: weno_constructor
  !< Abstract type used for create new concrete WENO interpolators.
  !<
  !< @note Every concrete WENO interpolator implementations must define their own constructor type.
  private
endtype weno_constructor

abstract interface
  elemental subroutine abstract_destructor(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  class(weno_interpolator), intent(INOUT) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_destructor

  subroutine abstract_constructor(self, constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create a WENO interpolator.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  import :: weno_constructor
  class(weno_interpolator), intent(INOUT) :: self        !< WENO interpolator.
  class(weno_constructor),  intent(IN)    :: constructor !< WENO constructor.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_constructor

  pure subroutine abstract_description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing a WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------
  import :: weno_interpolator
  class(weno_interpolator),      intent(IN)  :: self   !< WENO interpolator.
  character(len=:), allocatable, intent(OUT) :: string !< String returned.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_description
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_interpolator
