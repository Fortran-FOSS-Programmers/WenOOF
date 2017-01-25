!< Wenoof kappa factory.
module wenoof_kappa_factory
!< Wenoof kappa factory.

! use wenoof_base_object
use wenoof_kappa_object
use wenoof_kappa_rec_js

implicit none
private
public :: kappa_factory

type :: kappa_factory
  !< Factory, create an instance of concrete extension of [[kappa_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create !< Create a concrete instance of [[kappa_object]].
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
endmodule wenoof_kappa_factory
