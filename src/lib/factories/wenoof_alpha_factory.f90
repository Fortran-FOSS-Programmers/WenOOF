!< Wenoof alpha factory.
module wenoof_alpha_factory
!< Wenoof alpha factory.

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
use wenoof_alpha_object
use wenoof_alpha_rec_js
use wenoof_alpha_rec_m
use wenoof_alpha_rec_z

implicit none
private
public :: alpha_factory

type :: alpha_factory
  !< Factory, create an instance of concrete extension of [[alpha_object]] given its constructor.
  contains
    ! public methods
    procedure, nopass :: create             !< Create a concrete instance of [[alpha_object]].
    procedure, nopass :: create_constructor !< Create a concrete instance of [[alpha_object_constructor]].
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

  subroutine create_constructor(interpolator_type, S, constructor, face_left, face_right, eps)
  !< Create an instance of concrete extension of [[alpha_object_constructor]].
  character(*),                                 intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)           :: S                 !< Stencils dimension.
  class(alpha_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  logical,                                      intent(in), optional :: face_left         !< Activate left-face interpolations.
  logical,                                      intent(in), optional :: face_right        !< Activate right-face interpolations.
  real(RPP),                                    intent(in), optional :: eps               !< Small epsilon to avoid zero/division.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    allocate(alpha_rec_js_constructor :: constructor)
  case('interpolator-M-JS')
    allocate(alpha_rec_m_constructor :: constructor)
    select type(constructor)
    type is(alpha_rec_m_constructor)
      constructor%base_type = 'JS'
    endselect
  case('interpolator-M-Z')
    allocate(alpha_rec_m_constructor :: constructor)
    select type(constructor)
    type is(alpha_rec_m_constructor)
      constructor%base_type = 'Z'
    endselect
  case('interpolator-Z')
    allocate(alpha_rec_z_constructor :: constructor)
  case('reconstructor-JS')
    allocate(alpha_rec_js_constructor :: constructor)
  case('reconstructor-M-JS')
    allocate(alpha_rec_m_constructor :: constructor)
    select type(constructor)
    type is(alpha_rec_m_constructor)
      constructor%base_type = 'JS'
    endselect
  case('reconstructor-M-Z')
    allocate(alpha_rec_m_constructor :: constructor)
    select type(constructor)
    type is(alpha_rec_m_constructor)
      constructor%base_type = 'Z'
    endselect
  case('reconstructor-Z')
    allocate(alpha_rec_z_constructor :: constructor)
  endselect
  call constructor%create(S=S, face_left=face_left, face_right=face_right, eps=eps)
  endsubroutine create_constructor
endmodule wenoof_alpha_factory
