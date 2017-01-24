!< Wenoof interpolator factory.
module wenoof_interpolator_factory
!< Wenoof interpolator factory.

! use wenoof_base_object
use wenoof_interpolator_object
use wenoof_interpolator_js
use wenoof_reconstructor_js

implicit none
private
public :: interpolator_factory

type :: interpolator_factory
  !< Factory, create an instance of concrete extension of [[interpolator_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create !< Create a concrete instance of [[interpolator_object]].
endtype interpolator_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[interpolator_object]] given its constructor.
  class(interpolator_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolator_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(interpolator_js_constructor)
    allocate(interpolator_js :: object)
  type is(reconstructor_js_constructor)
    allocate(reconstructor_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create
endmodule wenoof_interpolator_factory
