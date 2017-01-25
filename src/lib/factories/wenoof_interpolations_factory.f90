!< Wenoof interpolations factory.
module wenoof_interpolations_factory
!< Wenoof interpolations factory.

use wenoof_interpolations_object
use wenoof_interpolations_rec_js

implicit none
private
public :: interpolations_factory

type :: interpolations_factory
  !< Factory, create an instance of concrete extension of [[interpolations_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create !< Create a concrete instance of [[interpolations_object]].
endtype interpolations_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[interpolations_object]] given its constructor.
  class(interpolations_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolations_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(interpolations_rec_js_constructor)
    allocate(interpolations_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create
endmodule wenoof_interpolations_factory
