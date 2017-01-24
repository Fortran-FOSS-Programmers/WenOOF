!< Wenoof objects factory.
module wenoof_objects_factory
!< Wenoof factory.

use wenoof_alpha_factory
use wenoof_alpha_object
use wenoof_beta_factory
use wenoof_beta_object
use wenoof_kappa_factory
use wenoof_kappa_object
use wenoof_interpolations_factory
use wenoof_interpolations_object
use wenoof_interpolator_factory
use wenoof_interpolator_object
use wenoof_weights_factory
use wenoof_weights_object

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
                         create_interpolator_object,   &
                         create_weights_object !< Create a concrete instance of [[alpha_object]], [[beta_object]],
                                               !< [[kappa_object]], [[interpolations_object]], [[interpolator_object]] or
                                               !< [[weights_object]].
    ! private methods
    procedure, nopass, private :: create_alpha_object          !< Create a concrete instance of [[alpha_object]].
    procedure, nopass, private :: create_beta_object           !< Create a concrete instance of [[beta_object]].
    procedure, nopass, private :: create_kappa_object          !< Create a concrete instance of [[kappa_object]].
    procedure, nopass, private :: create_interpolations_object !< Create a concrete instance of [[interpolations_object]].
    procedure, nopass, private :: create_interpolator_object   !< Create a concrete instance of [[interpolator_object]].
    procedure, nopass, private :: create_weights_object        !< Create a concrete instance of [[weights_object]].
endtype objects_factory

contains
  subroutine create_alpha_object(constructor, object)
  !< Create an instance of concrete extension of [[alpha_object]] given its constructor.
  class(alpha_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(alpha_object), allocatable, intent(out) :: object      !< Object.
  type(alpha_factory)                           :: factory     !< The factory.

  call factory%create(constructor=constructor, object=object)
  endsubroutine create_alpha_object

  subroutine create_beta_object(constructor, object)
  !< Create an instance of concrete extension of [[beta_object]] given its constructor.
  class(beta_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(beta_object), allocatable, intent(out) :: object      !< Object.
  type(beta_factory)                           :: factory     !< The factory.

  call factory%create(constructor=constructor, object=object)
  endsubroutine create_beta_object

  subroutine create_kappa_object(constructor, object)
  !< Create an instance of concrete extension of [[kappa_object]] given its constructor.
  class(kappa_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(kappa_object), allocatable, intent(out) :: object      !< Object.
  type(kappa_factory)                           :: factory     !< The factory.

  call factory%create(constructor=constructor, object=object)
  endsubroutine create_kappa_object

  subroutine create_interpolations_object(constructor, object)
  !< Create an instance of concrete extension of [[interpolations_object]] given its constructor.
  class(interpolations_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolations_object), allocatable, intent(out) :: object      !< Object.
  type(interpolations_factory)                           :: factory     !< The factory.

  call factory%create(constructor=constructor, object=object)
  endsubroutine create_interpolations_object

  subroutine create_interpolator_object(constructor, object)
  !< Create an instance of concrete extension of [[interpolator_object]] given its constructor.
  class(interpolator_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(interpolator_object), allocatable, intent(out) :: object      !< Object.
  type(interpolator_factory)                           :: factory     !< The factory.

  call factory%create(constructor=constructor, object=object)
  endsubroutine create_interpolator_object

  subroutine create_weights_object(constructor, object)
  !< Create an instance of concrete extension of [[weights_object]] given its constructor.
  class(weights_object_constructor),  intent(in)  :: constructor !< Constructor.
  class(weights_object), allocatable, intent(out) :: object      !< Object.
  type(weights_factory)                           :: factory     !< The factory.

  call factory%create(constructor=constructor, object=object)
  endsubroutine create_weights_object
endmodule wenoof_objects_factory
