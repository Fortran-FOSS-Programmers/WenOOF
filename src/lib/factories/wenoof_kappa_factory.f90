!< Wenoof kappa factory.
module wenoof_kappa_factory
!< Wenoof kappa factory.

use, intrinsic :: iso_c_binding, only : C_BOOL
use penf, only : I_P, R_P
use wenoof_interpolations_object, only : interpolations_object_constructor
use wenoof_kappa_object, only : kappa_object, kappa_object_constructor
use wenoof_kappa_rec_js, only : kappa_rec_js, kappa_rec_js_constructor
use wenoof_kappa_int_js, only : kappa_int_js, kappa_int_js_constructor

implicit none
private
public :: kappa_factory

type :: kappa_factory
  !< Factory, create an instance of concrete extension of [[kappa_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create                                          !< Create a concrete instance of [[kappa_object]].
    procedure, nopass :: create_constructor_int                          !< Create [[kappa_object_constructor]] (interpolate).
    procedure, nopass :: create_constructor_rec                          !< Create [[kappa_object_constructor]] (reconstruct).
    generic           :: create_constructor => create_constructor_int, & !< Create a concrete instance
                                               create_constructor_rec    !< of [[kappa_object_constructor]].
endtype kappa_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[kappa_object]] given its constructor.
  class(kappa_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(kappa_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(kappa_rec_js_constructor)
    allocate(kappa_rec_js :: object)
  type is(kappa_int_js_constructor)
    allocate(kappa_int_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor_rec(interpolator_type, S, constructor, ror)
  !< Create an instance of concrete extension of [[kappa_object_constructor]].
  character(*),                                 intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)           :: S                 !< Stencils dimension.
  class(kappa_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  logical(kind=C_BOOL),                         intent(in), optional :: ror               !< Activate or not ROR strategy.

  allocate(kappa_rec_js_constructor :: constructor)
  call constructor%create(S=S, ror=ror)
  endsubroutine create_constructor_rec

  subroutine create_constructor_int(interpolator_type, S, x_target, interpolations_constructor, constructor, ror)
  !< Create an instance of concrete extension of [[kappa_object_constructor]].
  character(*),                                 intent(in)           :: interpolator_type          !< Type of the interpolator.
  integer(I_P),                                 intent(in)           :: S                          !< Stencils dimension.
  real(R_P),                                    intent(in)           :: x_target                   !< Coordinate of the interp point.
  class(interpolations_object_constructor),     intent(in)           :: interpolations_constructor !< interpolations constructor.
  class(kappa_object_constructor), allocatable, intent(out)          :: constructor                !< Constructor.
  logical(kind=C_BOOL),                         intent(in), optional :: ror                        !< Activate or not ROR strategy.

  allocate(kappa_int_js_constructor :: constructor)
  select type(constructor)
  type is(kappa_int_js_constructor)
    constructor%x_target = x_target
    call constructor%create(S=S, ror=ror)
    allocate(constructor%interpolations_constructor, mold=interpolations_constructor)
    constructor%interpolations_constructor = interpolations_constructor
  endselect
  endsubroutine create_constructor_int
endmodule wenoof_kappa_factory
