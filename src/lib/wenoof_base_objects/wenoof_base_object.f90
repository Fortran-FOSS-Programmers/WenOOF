!< Abstract base object, the ancestor of all.
module wenoof_base_object
!< Abstract base object, the ancestor of all.
!<
!< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...

implicit none
private
public :: base_object
public :: base_object_constructor

type :: base_object_constructor
  !< Abstract base object constructor.
endtype base_object_constructor

type, abstract :: base_object
  !< Abstract base object, the ancestor of all.
  !<
  !< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...
  contains
    ! deferred public methods
    procedure(create_interface),      pass(self), deferred :: create      !< Create object.
    procedure(description_interface), pass(self), deferred :: description !< Return object string-description.
    procedure(destroy_interface),     pass(self), deferred :: destroy     !< Destroy object.
endtype base_object

abstract interface
  !< Abstract interface of [base_object] methods.

  subroutine create_interface(self, constructor)
  !< Create object.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  import :: base_object
  import :: base_object_constructor
  class(base_object),              intent(inout) :: self        !< Object.
  class(base_object_constructor),  intent(in)    :: constructor !< Object constructor.
  endsubroutine create_interface

  pure function description_interface(self) result(string)
  !< Return object string-description.
  import :: base_object
  class(base_object), intent(in) :: self   !< Object.
  character(len=:), allocatable  :: string !< String-description.
  endfunction description_interface

  elemental subroutine destroy_interface(self)
  !< Destroy object
  import :: base_object
  class(base_object), intent(inout) :: self !< Object.
  endsubroutine destroy_interface
endinterface
endmodule wenoof_base_object
