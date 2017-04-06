!< WenOOF, WENO interpolation Object Oriented Fortran library
module wenoof
!< WenOOF, WENO interpolation Object Oriented Fortran library

use, intrinsic :: iso_c_binding
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

contains
  subroutine wenoof_create(interpolator_type, S, interpolator, face_left, face_right, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                            intent(in)           :: interpolator_type          !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                          !< Stencil dimension.
  class(interpolator_object), allocatable, intent(out)          :: interpolator               !< The concrete WENO interpolator.
  logical,                                 intent(in), optional :: face_left                  !< Activate left-face interpolations.
  logical,                                 intent(in), optional :: face_right                 !< Activate right-face interpolations.
  real(RPP),                               intent(in), optional :: eps                        !< Small epsilon to avoid zero-div.
  type(objects_factory)                                         :: factory                    !< The factory.

  call factory%create(interpolator_type=interpolator_type, &
                      S=S,                                 &
                      interpolator=interpolator,           &
                      face_left=face_left,                 &
                      face_right=face_right,               &
                      eps=eps)
  endsubroutine wenoof_create

  subroutine wenoof_reconstructor(S, stencil, interpolation) bind(c, name='wenoof_reconstructor')
  integer(C_INT32_T), intent(in)          :: S                      !< Stencils number/dimension.
  real(C_DOUBLE),     intent(in)          :: stencil(1:2, 1-S:-1+S) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(C_DOUBLE),     intent(out)         :: interpolation(1:2)     !< Result of the interpolation, [1:2].
  class(interpolator_object), allocatable :: reconstructor          !< The WENO reconstructor.

  call wenoof_create(interpolator_type='reconstructor-JS', S=S, interpolator=reconstructor)
  call reconstructor%interpolate(stencil=stencil, interpolation=interpolation)
  endsubroutine wenoof_reconstructor
endmodule wenoof
