!< WenOOF, WENO interpolation Object Oriented Fortran library
module wenoof
!< WenOOF, WENO interpolation Object Oriented Fortran library

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
use wenoof_interpolator_object
use wenoof_objects_factory

implicit none
private
public :: interpolator_object
public :: wenoof_create

interface wenoof_create
  module procedure wenoof_create_reconstructor
  module procedure wenoof_create_interpolator
end interface wenoof_create

contains
  subroutine wenoof_create_reconstructor(interpolator_type, S, interpolator, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                            intent(in)           :: interpolator_type          !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                          !< Stencil dimension.
  class(interpolator_object), allocatable, intent(out)          :: interpolator               !< The concrete WENO interpolator.
  real(RPP),                               intent(in), optional :: eps                        !< Small epsilon to avoid zero-div.
  type(objects_factory)                                         :: factory                    !< The factory.

  call factory%create(interpolator_type=interpolator_type, &
                      S=S,                                 &
                      interpolator=interpolator,           &
                      eps=eps)
  endsubroutine wenoof_create_reconstructor

  subroutine wenoof_create_interpolator(interpolator_type, S, x_target, interpolator, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                            intent(in)           :: interpolator_type  !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                  !< Stencil dimension.
  real(RPP),                               intent(in)           :: x_target           !< Coordinate of the interpolation point.
  class(interpolator_object), allocatable, intent(out)          :: interpolator       !< The concrete WENO interpolator.
  real(RPP),                               intent(in), optional :: eps                !< Small epsilon to avoid zero-div.
  real(RPP),                  allocatable                       :: stencil(:)         !< Stencil used for interpolation, [1-S:-1+S].
  integer(I_P)                                                  :: i                  !< Counter.
  type(objects_factory)                                         :: factory            !< The factory.

  if ((x_target < -0.5_RPP).or.(x_target > 0.5_RPP)) then
    error stop 'error: x_target must be between -0.5 and 0.5, that represent left and right cell interfaces'
  endif
  allocate(stencil(1-S:S-1))
  do i=1,2*S-1
    stencil(-S+i) = 1.0_RPP - S + i
  enddo
  call factory%create(interpolator_type=interpolator_type, &
                      S=S,                                 &
                      interpolator=interpolator,           &
                      stencil=stencil,                     &
                      x_target=x_target,                   &
                      eps=eps)
  endsubroutine wenoof_create_interpolator
endmodule wenoof
