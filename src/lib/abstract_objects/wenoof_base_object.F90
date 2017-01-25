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

type :: base_object_constructor
  !< Abstract base object constructor.
  integer(I_P) :: S=0_I_P           !< Stencils dimension.
  logical      :: face_left=.true.  !< Activate left-face interpolation computation.
  logical      :: face_right=.true. !< Activate right-face interpolation computation.
  real(RPP)    :: eps=EPS_DEF       !< Small epsilon to avoid division by zero.
endtype base_object_constructor

type, abstract :: base_object
  !< Abstract base object, the ancestor of all.
  !<
  !< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...
  integer(I_P) :: S=0_I_P     !< Stencils dimension.
  integer(I_P) :: f1=1_I_P    !< Lower bound of faces index.
  integer(I_P) :: f2=2_I_P    !< Upper bound of faces index.
  integer(I_P) :: ff=0_I_P    !< Offset (step) of faces index.
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

  pure function description_interface(self) result(string)
  !< Return object string-description.
  import :: base_object
  class(base_object), intent(in) :: self   !< Object.
  character(len=:), allocatable  :: string !< String-description.
  endfunction description_interface

  elemental subroutine destroy_interface(self)
  !< Destroy object.
  import :: base_object
  class(base_object), intent(inout) :: self !< Object.
  endsubroutine destroy_interface
endinterface

contains
  ! public non overridable methods
  subroutine create_(self, constructor)
  !< Create object.
  class(base_object),             intent(inout) :: self        !< Object.
  class(base_object_constructor), intent(in)    :: constructor !< Object constructor.

  call self%destroy_
  select type(constructor)
  class is(base_object_constructor)
    self%S = constructor%S
    if (constructor%face_left.and.constructor%face_right) then
      self%f1 = 1_I_P; self%f2 = 2_I_P; self%ff = 0_I_P
    elseif (constructor%face_left) then
      self%f1 = 1_I_P; self%f2 = 1_I_P; self%ff = 0_I_P
    elseif (constructor%face_right) then
      self%f1 = 2_I_P; self%f2 = 2_I_P; self%ff = -1_I_P
    endif
    self%eps = constructor%eps
  endselect
  endsubroutine create_

  elemental subroutine destroy_(self)
  !< Destroy object.
  class(base_object), intent(inout) :: self !< Object.

  self%S = 0_I_P
  self%f1 = 1_I_P
  self%f2 = 2_I_P
  self%ff = 0_I_P
  self%eps = EPS_DEF
  endsubroutine destroy_
endmodule wenoof_base_object
