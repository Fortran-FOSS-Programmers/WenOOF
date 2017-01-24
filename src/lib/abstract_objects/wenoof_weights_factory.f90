!< Weights factory.
module wenoof_weights_factory
!< Weights factory.

use wenoof_base_object
use wenoof_weights_object
use wenoof_weights_js

implicit none
private
public :: weights_factory

type :: weights_factory
  !< Weights factory, create an instance of concrete extension of [[weights_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create !< Create a concrete instance of [[weights_object]].
endtype weights_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[weights_object]] given its constructor.
  class(weights_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(weights_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(weights_js_constructor)
    allocate(weights_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create
endmodule wenoof_weights_factory
