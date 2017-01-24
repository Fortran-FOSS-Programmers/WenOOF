!< Wenoof beta factory.
module wenoof_beta_factory
!< Wenoof beta factory.

! use wenoof_base_object
use wenoof_beta_object
use wenoof_beta_rec_js

implicit none
private
public :: beta_factory

type :: beta_factory
  !< Factory, create an instance of concrete extension of [[beta_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create !< Create a concrete instance of [[beta_object]].
endtype beta_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[beta_object]] given its constructor.
  class(beta_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(beta_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(beta_rec_js_constructor)
    allocate(beta_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create
endmodule wenoof_beta_factory
