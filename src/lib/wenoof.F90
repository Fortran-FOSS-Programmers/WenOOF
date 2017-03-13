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

  call factory%create(interpolator_type='reconstructor-'//interpolator_type, &
                      S=S,                                                   &
                      interpolator=interpolator,                             &
                      eps=eps)
  endsubroutine wenoof_create_reconstructor

  subroutine wenoof_create_interpolator(interpolator_type, S, interpolator, stencil, x_target, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                            intent(in)           :: interpolator_type  !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                  !< Stencil dimension.
  class(interpolator_object), allocatable, intent(out)          :: interpolator       !< The concrete WENO interpolator.
  real(RPP),                               intent(in)           :: stencil(1-S:)      !< Stencil used for interpolation, [1-S:-1+S].
  real(RPP),                               intent(in)           :: x_target           !< Coordinate of the interpolation point.
  real(RPP),                               intent(in), optional :: eps                !< Small epsilon to avoid zero-div.
  type(objects_factory)                                         :: factory            !< The factory.

  call factory%create(interpolator_type='interpolator-'//interpolator_type, &
                      S=S,                                                  &
                      interpolator=interpolator,                            &
                      stencil=stencil,                                      &
                      x_target=x_target,                                    &
                      eps=eps)
  endsubroutine wenoof_create_interpolator
endmodule wenoof
