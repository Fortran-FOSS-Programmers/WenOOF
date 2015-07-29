module type_weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------
!< Abstract WENO interpolator object,
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, abstract, public :: weno_interpolator
  private
  contains
    procedure(abstract_interpolate), deferred :: interpolate
endtype weno_interpolator

abstract interface
  pure subroutine abstract_interpolate(self, stencil, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stecil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  import weno_interpolator
  class(weno_interpolator), intent(INOUT) :: self          !< weno interpolator.
  real,                     intent(IN)    :: stencil(1:)   !< Stencil used for the interpolation.
  real,                     intent(OUT)   :: interpolation !< result of the interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine abstract_interpolate
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_interpolator
