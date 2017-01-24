!< Wenoof objects factory.
module wenoof_objects_factory
!< Wenoof factory.

use wenoof_alpha_object
use wenoof_alpha_rec_js
use wenoof_alpha_rec_m
use wenoof_alpha_rec_z
use wenoof_base_object
use wenoof_beta_object
use wenoof_beta_rec_js
use wenoof_kappa_object
use wenoof_kappa_rec_js
use wenoof_interpolations_object
use wenoof_interpolations_rec_js
use wenoof_weights_object
use wenoof_weights_js

implicit none
private
public :: objects_factory

type :: objects_factory
  !< Factory, create an instance of concrete extension of [[base_object]] given its constructor.
  contains
    ! public methods
    generic :: create => create_alpha_object,          &
                         create_beta_object,           &
                         create_kappa_object,          &
                         create_interpolations_object, &
                         create_weights_object !< Create a concrete instance of [[alpha_object]], [[beta_object]],
                                               !< [[kappa_object]], [[interpolations_object]] or [[weights_object]].
    procedure, nopass :: create_base_object    !< Create a concrete instance of [[base_object]].
    ! private methods
    procedure, nopass, private :: create_alpha_object          !< Create a concrete instance of [[alpha_object]].
    procedure, nopass, private :: create_beta_object           !< Create a concrete instance of [[beta_object]].
    procedure, nopass, private :: create_kappa_object          !< Create a concrete instance of [[kappa_object]].
    procedure, nopass, private :: create_interpolations_object !< Create a concrete instance of [[interpolations_object]].
    procedure, nopass, private :: create_weights_object        !< Create a concrete instance of [[weights_object]].
endtype objects_factory

contains
  subroutine create_base_object(constructor, object)
  !< Create an instance of concrete extension of [[base_object]] given its constructor.
  class(base_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(base_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(alpha_rec_js_constructor)
    allocate(alpha_rec_js :: object)
  type is(alpha_rec_m_constructor)
    allocate(alpha_rec_m :: object)
  type is(alpha_rec_z_constructor)
    allocate(alpha_rec_z :: object)
  type is(beta_rec_js_constructor)
    allocate(beta_rec_js :: object)
  type is(kappa_rec_js_constructor)
    allocate(kappa_rec_js :: object)
  type is(interpolations_rec_js_constructor)
    allocate(interpolations_rec_js :: object)
  type is(weights_js_constructor)
    allocate(weights_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_base_object

  subroutine create_alpha_object(constructor, object)
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
  endsubroutine create_alpha_object

  subroutine create_beta_object(constructor, object)
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
  endsubroutine create_beta_object

  subroutine create_kappa_object(constructor, object)
  !< Create an instance of concrete extension of [[kappa_object]] given its constructor.
  class(kappa_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(kappa_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(kappa_rec_js_constructor)
    allocate(kappa_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_kappa_object

  subroutine create_interpolations_object(constructor, object)
  !< Create an instance of concrete extension of [[interpolations_object]] given its constructor.
  class(interpolations_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolations_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(interpolations_rec_js_constructor)
    allocate(interpolations_rec_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_interpolations_object

  subroutine create_weights_object(constructor, object)
  !< Create an instance of concrete extension of [[weights_object]] given its constructor.
  class(weights_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(weights_object), allocatable, intent(out) :: object      !< Object.

  select type(constructor)
  type is(weights_js_constructor)
    allocate(weights_js :: object)
  class default
    error stop 'error: WenOOF object factory do NOT know the constructor given'
  endselect
  call object%create(constructor=constructor)
  endsubroutine create_weights_object
endmodule wenoof_objects_factory
