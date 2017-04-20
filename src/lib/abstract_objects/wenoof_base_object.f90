!< Abstract base object, the ancestor of all.
module wenoof_base_object
!< Abstract base object, the ancestor of all.
!<
!< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...

use penf, only : I_P, R_P

implicit none
private
public :: base_object
public :: base_object_constructor

real(R_P), parameter :: EPS_DEF=10._R_P**(-6) !< Small epsilon to avoid division by zero, default value.

type, abstract :: base_object_constructor
  !< Abstract base object constructor.
  integer(I_P) :: S=0_I_P     !< Stencils dimension.
  real(R_P)    :: eps=EPS_DEF !< Small epsilon to avoid division by zero.
  contains
    ! public methods
    procedure, pass(self) :: create => create_base_object_constructor
    ! public operators
    generic :: assignment(=) => constr_assign_constr !< `=` overloading.
    ! public deferred methods
    procedure(constr_assign_constr_interface), pass(lhs),  deferred :: constr_assign_constr !< `=` operator.
    ! public non overridable methods
    procedure, pass(lhs),  non_overridable :: assign_ => assign_constr_ !< Assign object.
endtype base_object_constructor

type, abstract :: base_object
  !< Abstract base object, the ancestor of all.
  !<
  !< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...
  integer(I_P) :: S=0_I_P     !< Stencils dimension.
  real(R_P)    :: eps=EPS_DEF !< Small epsilon to avoid division by zero.
  contains
    ! public operators
    generic :: assignment(=) => object_assign_object !< `=` overloading.
    ! public deferred methods
    procedure(create_interface),               pass(self), deferred :: create               !< Create object.
    procedure(description_interface),          pass(self), deferred :: description          !< Return object string-description.
    procedure(destroy_interface),              pass(self), deferred :: destroy              !< Destroy object.
    procedure(object_assign_object_interface), pass(lhs),  deferred :: object_assign_object !< `=` operator.
    ! public non overridable methods
    procedure, pass(lhs),  non_overridable :: assign_  !< Assign object.
    procedure, pass(self), non_overridable :: create_  !< Create object.
    procedure, pass(self), non_overridable :: destroy_ !< Destroy object.
endtype base_object

abstract interface
  !< Abstract interfaces of [[base_object_constructor]].
  subroutine constr_assign_constr_interface(lhs, rhs)
  !< `=` operator.
  import :: base_object_constructor
  class(base_object_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor), intent(in)    :: rhs !< Right hand side.
  endsubroutine constr_assign_constr_interface
endinterface

abstract interface
  !< Abstract interfaces of [[base_object]].
  subroutine object_assign_object_interface(lhs, rhs)
  !< `=` operator.
  import :: base_object
  class(base_object), intent(inout) :: lhs !< Left hand side.
  class(base_object), intent(in)    :: rhs !< Right hand side.
  endsubroutine object_assign_object_interface

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
  real(R_P),                      intent(in), optional :: eps        !< Small epsilon to avoid division by zero.

  self%S = S
  if (present(eps)) self%eps = eps
  endsubroutine create_base_object_constructor

  ! public non overridable methods
  subroutine assign_constr_(lhs, rhs)
  !< Assign object constructor.
  class(base_object_constructor), intent(inout) :: lhs !< Left hand side.
  class(base_object_constructor), intent(in)    :: rhs !< Right hand side.

  lhs%S = rhs%S
  lhs%eps = rhs%eps
  endsubroutine assign_constr_

  ! base object

  ! public non overridable methods
  subroutine assign_(lhs, rhs)
  !< Assign object.
  class(base_object), intent(inout) :: lhs !< Left hand side.
  class(base_object), intent(in)    :: rhs !< Right hand side.

  lhs%S = rhs%S
  lhs%eps = rhs%eps
  endsubroutine assign_

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
