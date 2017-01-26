!< Wenoof kappa factory.
module wenoof_kappa_factory
!< Wenoof kappa factory.

use penf, only: I_P
use wenoof_kappa_object
use wenoof_kappa_rec_js

implicit none
private
public :: kappa_factory

type :: kappa_factory
  !< Factory, create an instance of concrete extension of [[kappa_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create             !< Create a concrete instance of [[kappa_object]].
    procedure, nopass :: create_constructor !< Create a concrete instance of [[kappa_object_constructor]].
endtype kappa_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[kappa_object]] given its constructor.
  class(kappa_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(kappa_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(kappa_rec_js_constructor)
    allocate(kappa_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create

  subroutine create_constructor(interpolator_type, S, constructor)
  !< Create an instance of concrete extension of [[kappa_object_constructor]].
  character(*),                                 intent(in)  :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)  :: S                 !< Stencils dimension.
  class(kappa_object_constructor), allocatable, intent(out) :: constructor       !< Constructor.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    ! @TODO implement this
    error stop 'interpolator-JS to be implemented'
  case('reconstructor-JS')
    allocate(kappa_rec_js_constructor :: constructor)
  case('reconstructor-M-JS')
    allocate(kappa_rec_js_constructor :: constructor)
  case('reconstructor-M-Z')
    allocate(kappa_rec_js_constructor :: constructor)
  case('reconstructor-Z')
    allocate(kappa_rec_js_constructor :: constructor)
  endselect
  call constructor%create(S=S)
  endsubroutine create_constructor
endmodule wenoof_kappa_factory
