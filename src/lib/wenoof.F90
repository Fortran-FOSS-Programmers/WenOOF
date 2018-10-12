!< WenOOF, WENO interpolation Object Oriented Fortran library
module wenoof
!< WenOOF, WENO interpolation Object Oriented Fortran library

!use, intrinsic :: iso_c_binding
use penf, only : I_P, R_P
use wenoof_interpolator_object, only : interpolator_object
use wenoof_objects_factory, only : objects_factory

implicit none
private
public :: interpolator_object
public :: wenoof_create
!#ifndef r16p
!public :: wenoof_initialize_c_wrap
!public :: wenoof_interpolate_c_wrap
!public :: wenoof_reconstruct_c_wrap
!#endif

interface wenoof_create
  module procedure wenoof_create_reconstructor
  module procedure wenoof_create_interpolator
end interface wenoof_create

!#ifndef r16p
!class(interpolator_object), allocatable :: interpolator_c_wrap !< The WENO interpolator/reconstructor for C/Python wrappers.
!#endif

contains
  subroutine wenoof_create_reconstructor(interpolator_type, S, interpolator, ror, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                            intent(in)           :: interpolator_type   !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                   !< Stencil dimension.
  class(interpolator_object), allocatable, intent(out)          :: interpolator        !< The concrete WENO interpolator.
  logical,                                 intent(in), optional :: ror                 !< Activate or not ROR strategy.
  real(R_P),                               intent(in), optional :: eps                 !< Small epsilon to avoid zero-div.
  type(objects_factory)                                         :: factory             !< The factory.

  call factory%create(interpolator_type=interpolator_type, &
                      S=S,                                 &
                      interpolator=interpolator,           &
                      ror=ror,                             &
                      eps=eps)
  endsubroutine wenoof_create_reconstructor

  subroutine wenoof_create_interpolator(interpolator_type, S, x_target, interpolator, ror, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                            intent(in)           :: interpolator_type  !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                  !< Stencil dimension.
  real(R_P),                               intent(in)           :: x_target           !< Coordinate of the interpolation point.
  class(interpolator_object), allocatable, intent(out)          :: interpolator       !< The concrete WENO interpolator.
  logical,                                 intent(in), optional :: ror                !< Activate or not ROR strategy.
  real(R_P),                               intent(in), optional :: eps                !< Small epsilon to avoid zero-div.
  integer(I_P)                                                  :: i                  !< Counter.
  type(objects_factory)                                         :: factory            !< The factory.

  if ((x_target < -0.5_R_P).or.(x_target > 0.5_R_P)) then
    error stop 'error: x_target must be between -0.5 and 0.5, that represent left and right cell interfaces'
  endif
  call factory%create(interpolator_type=interpolator_type, &
                      S=S,                                 &
                      x_target=x_target,                   &
                      interpolator=interpolator,           &
                      ror=ror,                             &
                      eps=eps)
  endsubroutine wenoof_create_interpolator

!#ifndef r16p
!  ! procedure for C/Python wrappers
!  subroutine wenoof_initialize_c_wrap(interpolator_type, S, ror, eps, x_target, verbose) bind(c, name='wenoof_initialize_c_wrap')
!  !< Intialize the WENO interpolator for C/Python wrappers.
!  !<
!  !< @note For the arguments the following conventions hold:
!  !<
!  !<+ `interpolator_type` must be null-terminated or new-line-terminated;
!  !<+ `ror` must be logical value, otherwise default value is taken;
!  !<+ `eps` must be positive, otherwise default value is taken;
!  !<+ `x_target` must be in [-0.5,0.5], otherwise a reconstruction is initialized instead of an interpolation;
!  !<+ `verbose` must be equal to 1, otherwise quite mode is activated.
!  character(kind=C_CHAR), intent(in)        :: interpolator_type(*) !< Type of the interpolator.
!  integer(C_INT32_T),     intent(in), value :: S                    !< Stencils number/dimension.
!  logical(C_BOOL),        intent(in), value :: ror                  !< ROR strategy switch.
!  real(C_DOUBLE),         intent(in), value :: eps                  !< Small epsilon to avoid zero-div.
!  real(C_DOUBLE),         intent(in), value :: x_target             !< Coordinate of the interpolation point.
!  integer(C_INT32_T),     intent(in), value :: verbose              !< Activate verbose output.
!  character(len=:), allocatable             :: interpolator_type_   !< Type of the interpolator, local variable.
!  integer(I_P)                              :: c                    !< Counter.
!
!  interpolator_type_ = ''
!  c = 1
!  do
!    if (interpolator_type(c) == new_line('a').or.interpolator_type(c) == char(0)) exit
!    interpolator_type_ = interpolator_type_//interpolator_type(c)
!    c = c + 1
!  enddo
!  if (-0.5_C_DOUBLE <= x_target .and. x_target <= 0.5_C_DOUBLE) then
!    if (verbose==1) print '(A,F23.15)', 'interpolate at ', x_target
!    if (eps > 0._C_DOUBLE) then
!       if (ror .eqv. .false.) then
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, x_target=x_target, interpolator=interpolator_c_wrap, &
!                            ror=ror, eps=eps)
!       else
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, x_target=x_target, interpolator=interpolator_c_wrap, &
!                            eps=eps)
!      endif
!    else
!       if (ror .eqv. .false.) then
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, x_target=x_target, interpolator=interpolator_c_wrap, &
!                            ror=ror)
!       else
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, x_target=x_target, interpolator=interpolator_c_wrap)
!      endif
!    endif
!  else
!    if (eps > 0._C_DOUBLE) then
!       if (ror .eqv. .false.) then
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, interpolator=interpolator_c_wrap, ror=ror, eps=eps)
!       else
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, interpolator=interpolator_c_wrap, eps=eps)
!      endif
!    else
!       if (ror .eqv. .false.) then
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, interpolator=interpolator_c_wrap, ror=ror)
!       else
!         call wenoof_create(interpolator_type=interpolator_type_, S=S, interpolator=interpolator_c_wrap)
!      endif
!    endif
!  endif
!  if (verbose==1) print '(A)', interpolator_c_wrap%description()
!  endsubroutine wenoof_initialize_c_wrap
!
!  subroutine wenoof_interpolate_c_wrap(S, ord, stencil, interpolation) bind(c, name='wenoof_interpolate_c_wrap')
!  !< Interpolate over stencils values by means of WenOOF interpolator.
!  integer(C_INT32_T), intent(in), value  :: S                 !< Stencils number/dimension.
!  integer(C_INT32_T), intent(in), value  :: ord               !< Interpolation order.
!  real(C_DOUBLE),     intent(in)         :: stencil(1-S:-1+S) !< Stencil of the interpolation.
!  real(C_DOUBLE),     intent(out)        :: interpolation     !< Result of the interpolation.
!
!  call interpolator_c_wrap%interpolate(ord=ord, stencil=stencil, interpolation=interpolation)
!  endsubroutine wenoof_interpolate_c_wrap
!
!  subroutine wenoof_reconstruct_c_wrap(S, ord, stencil, interpolation) bind(c, name='wenoof_reconstruct_c_wrap')
!  !< Reconstruct over stencils values by means of WenOOF reconstructor.
!  integer(C_INT32_T), intent(in), value   :: S                      !< Stencils number/dimension.
!  integer(C_INT32_T), intent(in), value   :: ord                    !< Reconstruction order.
!  real(C_DOUBLE),     intent(in)          :: stencil(1:2, 1-S:-1+S) !< Stencil of the interpolation.
!  real(C_DOUBLE),     intent(out)         :: interpolation(1:2)     !< Result of the interpolation.
!
!  call interpolator_c_wrap%interpolate(ord=ord, stencil=stencil, interpolation=interpolation)
!  endsubroutine wenoof_reconstruct_c_wrap
!#endif
endmodule wenoof
