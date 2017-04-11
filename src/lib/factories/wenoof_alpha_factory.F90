!< Wenoof alpha factory.
module wenoof_alpha_factory
!< Wenoof alpha factory.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
use wenoof_alpha_object, only : alpha_object, alpha_object_constructor
use wenoof_alpha_int_js, only : alpha_int_js, alpha_int_js_constructor
use wenoof_alpha_int_m, only : alpha_int_m, alpha_int_m_constructor
use wenoof_alpha_int_z, only : alpha_int_z, alpha_int_z_constructor
use wenoof_alpha_rec_js, only : alpha_rec_js, alpha_rec_js_constructor
use wenoof_alpha_rec_m, only : alpha_rec_m, alpha_rec_m_constructor
use wenoof_alpha_rec_z, only : alpha_rec_z, alpha_rec_z_constructor

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
  type is(alpha_int_js_constructor)
    allocate(alpha_int_js :: object)
  type is(alpha_int_m_constructor)
    allocate(alpha_int_m :: object)
  type is(alpha_int_z_constructor)
    allocate(alpha_int_z :: object)
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

  subroutine create_constructor(interpolator_type, S, constructor, eps)
  !< Create an instance of concrete extension of [[alpha_object_constructor]].
  character(*),                                 intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)           :: S                 !< Stencils dimension.
  class(alpha_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  real(RPP),                                    intent(in), optional :: eps               !< Small epsilon to avoid zero/division.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    allocate(alpha_int_js_constructor :: constructor)
  case('interpolator-M-JS')
    allocate(alpha_int_m_constructor :: constructor)
    select type(constructor)
    type is(alpha_int_m_constructor)
      constructor%base_type = 'JS'
    endselect
  case('interpolator-M-Z')
    allocate(alpha_int_m_constructor :: constructor)
    select type(constructor)
    type is(alpha_int_m_constructor)
      constructor%base_type = 'Z'
    endselect
  case('interpolator-Z')
    allocate(alpha_int_z_constructor :: constructor)
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
  case default
    write(stderr, '(A)') 'error: interpolator type "'//trim(adjustl(interpolator_type))//'" is unknown!'
    stop
  endselect
  call constructor%create(S=S, eps=eps)
  endsubroutine create_constructor
endmodule wenoof_alpha_factory
