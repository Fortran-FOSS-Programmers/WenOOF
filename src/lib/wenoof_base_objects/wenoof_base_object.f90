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
  integer(I_P) :: S=0_I_P           !< Stencils dimension.
  logical      :: face_left=.true.  !< Activate left-face interpolation computation.
  logical      :: face_right=.true. !< Activate right-face interpolation computation.
  contains
    ! public methods
    procedure, pass(self) :: create  => create_base_object_constructor  !< Create object constructor.
    procedure, pass(self) :: destroy => destroy_base_object_constructor !< Destroy object constructor.
endtype base_object_constructor

type :: base_object
  !< Abstract base object, the ancestor of all.
  !<
  !< Define a minimal, base, object that is used as ancestor of all objects, e.g. smoothness indicator, optimal weights, etc...
  integer(I_P) :: S=0_I_P  !< Stencils dimension.
  integer(I_P) :: f1=1_I_P !< Lower bound of faces index.
  integer(I_P) :: f2=2_I_P !< Upper bound of faces index.
  integer(I_P) :: ff=0_I_P !< Offset (step) of faces index.
  contains
    ! public deferred methods
    procedure, pass(self), deferred :: description !< Return object string-description.
    ! public methods
    procedure, pass(self), deferred :: create  !< Create object.
    procedure, pass(self), deferred :: destroy !< Destroy object.
endtype base_object

contains
  ! constructor methods

  ! public methods
  subroutine create_base_object_constructor(self, S, face_left, face_right)
  !< Create object constructor.
  class(base_object_constructor), intent(inout)        :: self       !< Object constructor.
  integer(I_P),                   intent(in)           :: S          !< Stencils dimension.
  logical,                        intent(in), optional :: face_left  !< Activate left-face interpolations.
  logical,                        intent(in), optional :: face_right !< Activate right-face interpolations.

  constructor%S = S
  constructor%face_left = .true.  ; if (present(face_left)) constructor%face_left = face_left
  constructor%face_right = .true. ; if (present(face_right)) constructor%face_right = face_right
  endsubroutine create_base_object_constructor

  pure subroutine destroy_base_object_constructor(self)
  !< Destroy object constructor.
  class(base_object_constructor), intent(inout) :: self !< Object constructor.

  self%S = 0_I_P
  constructor%face_left = .true.
  constructor%face_right = .true.
  endsubroutine destroy_base_object_constructor

  ! base_object methods

  ! public deferred methods
  pure function description(self) result(string)
  !< Return object string-description.
  class(base_object), intent(in) :: self   !< Object.
  character(len=:), allocatable  :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'base_object%description to be implemented by a concrete extension of base_object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create object.
  !<
  !< @note Before call this method a concrete constructor must be instantiated.
  class(base_object),              intent(inout) :: self        !< Object.
  class(base_object_constructor),  intent(in)    :: constructor !< Object constructor.

  call self%destroy
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
  class default
    ! @TODO add error handling
  endselect
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy object
  class(base_object), intent(inout) :: self !< Object.

  self%S = 0_I_P
  self%f1 = 1_I_P
  self%f2 = 2_I_P
  self%ff = 0_I_P
  endsubroutine destroy
endmodule wenoof_base_object
