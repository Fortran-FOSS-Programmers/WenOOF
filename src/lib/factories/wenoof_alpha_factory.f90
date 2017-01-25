!< Wenoof alpha factory.
module wenoof_alpha_factory
!< Wenoof alpha factory.

use wenoof_alpha_object
use wenoof_alpha_rec_js
use wenoof_alpha_rec_m
use wenoof_alpha_rec_z
! use wenoof_base_object

implicit none
private
public :: alpha_factory

type :: alpha_factory
  !< Factory, create an instance of concrete extension of [[alpha_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create !< Create a concrete instance of [[alpha_object]].
endtype alpha_factory

contains
  subroutine create(constructor, object)
  !< Create an instance of concrete extension of [[alpha_object]] given its constructor.
  class(alpha_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(alpha_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(alpha_rec_js_constructor)
    allocate(alpha_rec_js :: object)
  type is(alpha_rec_m_constructor)
    allocate(alpha_rec_m :: object)
  type is(alpha_rec_z_constructor)
    allocate(alpha_rec_z :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create
endmodule wenoof_alpha_factory
