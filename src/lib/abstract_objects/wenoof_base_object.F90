!< Abstract base object, the ancestor of all.
module wenoof_base_object
!< Abstract base object, the ancestor of all.
!<
!< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif

implicit none
private
public :: base_object
public :: base_object_constructor

real(RPP), parameter :: EPS_DEF=10._RPP**(-6) !< Small epsilon to avoid division by zero, default value.

type, abstract :: base_object_constructor
  !< Abstract base object constructor.
  integer(I_P) :: S=0_I_P           !< Stencils dimension.
  real(RPP)    :: eps=EPS_DEF       !< Small epsilon to avoid division by zero.
  contains
    procedure, pass(self) :: create => create_base_object_constructor
endtype base_object_constructor

type, abstract :: base_object
  !< Abstract base object, the ancestor of all.
  !<
  !< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...
  integer(I_P) :: S=0_I_P     !< Stencils dimension.
  real(RPP)    :: eps=EPS_DEF !< Small epsilon to avoid division by zero.
  contains
    ! public deferred methods
    procedure(create_interface),      pass(self), deferred :: create      !< Create object.
    procedure(description_interface), pass(self), deferred :: description !< Return object string-description.
    procedure(destroy_interface),     pass(self), deferred :: destroy     !< Destroy object.
    ! public non overridable methods
    procedure, pass(self), non_overridable :: create_  !< Create object.
    procedure, pass(self), non_overridable :: destroy_ !< Destroy object.
endtype base_object

abstract interface
  !< Abstract interfaces of [[base_object]].
  subroutine create_interface(self, constructor)
  !< Create object.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  import :: base_object, base_object_constructor
  class(base_object),              intent(inout) :: self        !< Object.
  class(base_object_constructor),  intent(in)    :: constructor !< Object constructor.
  endsubroutine create_interface

  pure function description_interface(self, prefix) result(string)
  !< Return object string-description.
  import :: base_object
  class(base_object), intent(in)           :: self   !< Object.
  character(len=*),   intent(in), optional :: prefix !< Prefixing string.
  character(len=:), allocatable            :: string !< String-description.
  endfunction description_interface

  elemental subroutine destroy_interface(self)
  !< Destroy object.
  import :: base_object
  class(base_object), intent(inout) :: self !< Object.
  endsubroutine destroy_interface
endinterface

contains
  ! base object constructor

  ! public methods
  subroutine create_base_object_constructor(self, S, eps)
  !< Create alpha constructor.
  class(base_object_constructor), intent(inout)        :: self       !< Constructor.
  integer(I_P),                   intent(in)           :: S          !< Stencils dimension.
  real(RPP),                      intent(in), optional :: eps        !< Small epsilon to avoid division by zero.

  self%S = S
  if (present(eps)) self%eps = eps
  endsubroutine create_base_object_constructor

  ! base object

  ! public non overridable methods
  subroutine create_(self, constructor)
  !< Create object.
  class(base_object),             intent(inout) :: self        !< Object.
  class(base_object_constructor), intent(in)    :: constructor !< Object constructor.

  call self%destroy_
  select type(constructor)
  class is(base_object_constructor)
    self%S = constructor%S
    self%eps = constructor%eps
  endselect
  endsubroutine create_

  elemental subroutine destroy_(self)
  !< Destroy object.
  class(base_object), intent(inout) :: self !< Object.

  self%S = 0_I_P
  self%eps = EPS_DEF
  endsubroutine destroy_
endmodule wenoof_base_object
