!< Wenoof kappa factory.
module wenoof_kappa_factory
!< Wenoof kappa factory.

use penf, only: I_P
use wenoof_kappa_object
use wenoof_kappa_rec_js
use wenoof_kappa_int_js

implicit none
private
public :: kappa_factory

type :: kappa_factory
  !< Factory, create an instance of concrete extension of [[kappa_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create                                          !< Create a concrete instance of [[kappa_object]].
    procedure, nopass :: create_constructor_rec
    procedure, nopass :: create_constructor_int
    generic           :: create_constructor => create_constructor_rec, & !< Create a concrete instance
                                               create_constructor_int    !< of [[kappa_object_constructor]].
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

  subroutine create_constructor_rec(interpolator_type, S, constructor)
  !< Create an instance of concrete extension of [[kappa_object_constructor]].
  character(*),                                 intent(in)  :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)  :: S                 !< Stencils dimension.
  class(kappa_object_constructor), allocatable, intent(out) :: constructor       !< Constructor.

  allocate(kappa_rec_js_constructor :: constructor)
  call constructor%create(S=S)
  endsubroutine create_constructor_rec

  subroutine create_constructor_int(interpolator_type, S, stencil, x_target, constructor)
  !< Create an instance of concrete extension of [[kappa_object_constructor]].
  character(*),                                 intent(in)  :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)  :: S                 !< Stencils dimension.
  real(RPP),                                    intent(in)  :: stencil(1-S:)     !< Stencil used for inter, [1-S:-1+S].
  real(RPP),                                    intent(in)  :: x_target          !< Coordinate of the interp point.
  class(kappa_object_constructor), allocatable, intent(out) :: constructor       !< Constructor.

  allocate(kappa_int_js_constructor :: constructor)
  allocate(stencil  :: constructor%stencil)
  constructor%x_target = x_target
  call constructor%create(S=S)
  endsubroutine create_constructor_int
endmodule wenoof_kappa_factory
