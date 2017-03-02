!< Wenoof objects factory.
module wenoof_objects_factory
!< Wenoof factory.

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
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
                         create_reconstructor,         &
                         create_interpolator,          &
                         create_interpolator_object,   &
                         create_weights_object !< Create a concrete instance of [[alpha_object]], [[beta_object]],
                                               !< [[kappa_object]], [[interpolations_object]], [[interpolator_object]] or
                                               !< [[weights_object]].
    generic :: create_constructor => create_alpha_object_constructor,              &
                                     create_beta_object_constructor,               &
                                     create_kappa_rec_object_constructor,          &
                                     create_kappa_int_object_constructor,          &
                                     create_interpolations_rec_object_constructor, &
                                     create_interpolations_int_object_constructor, &
                                     create_interpolator_object_constructor,       &
                                     create_weights_object_constructor !< Create a concrete instance of
                                                                       !< [[alpha_object_constructor]], [[beta_object_constructor]],
                                                                       !< [[kappa_object_constructor]],
                                                                       !< [[interpolations_object_constructor]],
                                                                       !< [[interpolator_object_constructor]] or
                                                                       !< [[weights_object_constructor]].
    ! private methods
    procedure, nopass,     private :: create_alpha_object                      !< Create [[alpha_object]] instance
    procedure, nopass,     private :: create_beta_object                       !< Create [[beta_object]] instance.
    procedure, nopass,     private :: create_kappa_object                      !< Create [[kappa_object]] instance.
    procedure, nopass,     private :: create_interpolations_object             !< Create [[interpolations_object]] instance.
    procedure, pass(self), private :: create_interpolator                      !< Create [[interpolator_object]] instance.
    procedure, nopass,     private :: create_interpolator_object               !< Create [[interpolator_object]] instance.
    procedure, nopass,     private :: create_weights_object                    !< Create [[weights_object]] instance.
    procedure, nopass,     private :: create_alpha_object_constructor          !< Create [[alpha_object_constructor]] instance.
    procedure, nopass,     private :: create_beta_object_constructor           !< Create [[beta_object_constructor]] instance.
    procedure, nopass,     private :: create_kappa_object_constructor          !< Create [[kappa_object_constructor]] instance.
    procedure, nopass,     private :: create_interpolations_object_constructor !< Create [[interpolations_object_constructor]] inst.
    procedure, nopass,     private :: create_interpolator_object_constructor   !< Create [[interpolator_object_constructor]] inst.
    procedure, nopass,     private :: create_weights_object_constructor        !< Create [[weights_object_constructor]] instance.
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

  subroutine create_reconstructor(self, interpolator_type, S, interpolator, eps)
  !< Create an instance of concrete extension of [[interpolator_object]] given user options.
  class(objects_factory),                  intent(in)           :: self                       !< The factory.
  character(*),                            intent(in)           :: interpolator_type          !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                          !< Stencils dimension.
  class(interpolator_object), allocatable, intent(out)          :: interpolator               !< Interpolator.
  real(RPP),                               intent(in), optional :: eps                        !< Small epsilon to avoid zero/div.
  class(alpha_object_constructor),          allocatable         :: alpha_constructor          !< Alpha constructor.
  class(beta_object_constructor),           allocatable         :: beta_constructor           !< Beta constructor.
  class(interpolations_object_constructor), allocatable         :: interpolations_constructor !< Interpolations constructor.
  class(kappa_object_constructor),          allocatable         :: kappa_constructor          !< Kappa constructor.
  class(weights_object_constructor),        allocatable         :: weights_constructor        !< Weights constructor.
  class(interpolator_object_constructor),   allocatable         :: interpolator_constructor   !< Interpolator constructor.

  call self%create_constructor(interpolator_type=interpolator_type, &
                               S=S,                                 &
                               constructor=alpha_constructor,       &
                               eps=eps)

  call self%create_constructor(interpolator_type=interpolator_type, &
                               S=S,                                 &
                               constructor=beta_constructor)

  call self%create_constructor(interpolator_type=interpolator_type, S=S, constructor=kappa_constructor)

  call self%create_constructor(interpolator_type=interpolator_type, &
                               S=S,                                 &
                               alpha_constructor=alpha_constructor, &
                               beta_constructor=beta_constructor,   &
                               kappa_constructor=kappa_constructor, &
                               constructor=weights_constructor)

  call self%create_constructor(interpolator_type=interpolator_type,    &
                               S=S,                                    &
                               constructor=interpolations_constructor)

  call self%create_constructor(interpolator_type=interpolator_type,                   &
                               S=S,                                                   &
                               interpolations_constructor=interpolations_constructor, &
                               weights_constructor=weights_constructor,               &
                               constructor=interpolator_constructor)

  call self%create_interpolator_object(constructor=interpolator_constructor, object=interpolator)
  endsubroutine create_reconstructor

  subroutine create_interpolator(self, interpolator_type, S, interpolator, stencil, xtarget, eps)
  !< Create an instance of concrete extension of [[interpolator_object]] given user options.
  class(objects_factory),                  intent(in)           :: self                       !< The factory.
  character(*),                            intent(in)           :: interpolator_type          !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                          !< Stencils dimension.
  class(interpolator_object), allocatable, intent(out)          :: interpolator               !< Interpolator.
  real(RPP),                               intent(in)           :: stencil(1-S:)              !< Stencil used for inter, [1-S:-1+S].
  real(RPP),                               intent(in)           :: x_target                   !< Coordinate of the interp point.
  real(RPP),                               intent(in), optional :: eps                        !< Small epsilon to avoid zero/div.
  class(alpha_object_constructor),          allocatable         :: alpha_constructor          !< Alpha constructor.
  class(beta_object_constructor),           allocatable         :: beta_constructor           !< Beta constructor.
  class(interpolations_object_constructor), allocatable         :: interpolations_constructor !< Interpolations constructor.
  class(kappa_object_constructor),          allocatable         :: kappa_constructor          !< Kappa constructor.
  class(weights_object_constructor),        allocatable         :: weights_constructor        !< Weights constructor.
  class(interpolator_object_constructor),   allocatable         :: interpolator_constructor   !< Interpolator constructor.

  call self%create_constructor(interpolator_type=interpolator_type, &
                               S=S,                                 &
                               constructor=alpha_constructor,       &
                               eps=eps)

  call self%create_constructor(interpolator_type=interpolator_type, &
                               S=S,                                 &
                               constructor=beta_constructor)

  call self%create_constructor(interpolator_type=interpolator_type,    &
                               S=S,                                    &
                               stencil=stencil,                        &
                               x_target=x_target,                      &
                               constructor=interpolations_constructor)

  call self%create_constructor(interpolator_type=interpolator_type, &
                               S=S,                                 &
                               stencil=stencil,                     &
                               x_target=x_target,                   &
                               constructor=kappa_constructor)

  call self%create_constructor(interpolator_type=interpolator_type, &
                               S=S,                                 &
                               alpha_constructor=alpha_constructor, &
                               beta_constructor=beta_constructor,   &
                               kappa_constructor=kappa_constructor, &
                               constructor=weights_constructor)

  call self%create_constructor(interpolator_type=interpolator_type,                   &
                               S=S,                                                   &
                               interpolations_constructor=interpolations_constructor, &
                               weights_constructor=weights_constructor,               &
                               constructor=interpolator_constructor)

  call self%create_interpolator_object(constructor=interpolator_constructor, object=interpolator)
  endsubroutine create_interpolator

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

  subroutine create_alpha_object_constructor(interpolator_type, S, constructor, eps)
  !< Create an instance of concrete extension of [[alpha_object_constructor]].
  character(*),                                 intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)           :: S                 !< Stencils dimension.
  class(alpha_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  real(RPP),                                    intent(in), optional :: eps               !< Small epsilon to avoid zero/division.
  type(alpha_factory)                                                :: factory           !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type, &
                                  S=S,                                 &
                                  constructor=constructor,             &
                                  eps=eps)
  endsubroutine create_alpha_object_constructor

  subroutine create_beta_object_constructor(interpolator_type, S, constructor)
  !< Create an instance of concrete extension of [[beta_object_constructor]].
  character(*),                                intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                intent(in)           :: S                 !< Stencils dimension.
  class(beta_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  type(beta_factory)                                                :: factory           !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type, &
                                  S=S,                                 &
                                  constructor=constructor)
  endsubroutine create_beta_object_constructor

  subroutine create_kappa_int_object_constructor(interpolator_type, S, stencil, x_target, constructor)
  !< Create an instance of concrete extension of [[kappa_object_constructor]].
  character(*),                                 intent(in)  :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)  :: S                 !< Stencils dimension.
  real(RPP),                                    intent(in)  :: stencil(1-S:)     !< Stencil used for inter, [1-S:-1+S].
  real(RPP),                                    intent(in)  :: x_target          !< Coordinate of the interp point.
  class(kappa_object_constructor), allocatable, intent(out) :: constructor       !< Constructor.
  type(kappa_factory)                                       :: factory           !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type, &
                                  S=S,                                 &
                                  stencil=stencil,                     &
                                  x_target=x_target,                   &
                                  constructor=constructor)
  endsubroutine create_kappa_int_object_constructor

  subroutine create_kappa_rec_object_constructor(interpolator_type, S, constructor)
  !< Create an instance of concrete extension of [[kappa_object_constructor]].
  character(*),                                 intent(in)  :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                 intent(in)  :: S                 !< Stencils dimension.
  class(kappa_object_constructor), allocatable, intent(out) :: constructor       !< Constructor.
  type(kappa_factory)                                       :: factory           !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type, &
                                  S=S,                                 &
                                  constructor=constructor)
  endsubroutine create_kappa_rec_object_constructor

  subroutine create_interpolations_rec_object_constructor(interpolator_type, S, constructor)
  !< Create an instance of concrete extension of [[interpolations_object_constructor]].
  character(*),                                          intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                          intent(in)           :: S                 !< Stencils dimension.
  class(interpolations_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  type(interpolations_factory)                                                :: factory           !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type, &
                                  S=S,                                 &
                                  constructor=constructor)
  endsubroutine create_interpolations_rec_object_constructor

  subroutine create_interpolations_int_object_constructor(interpolator_type, S, stencil, x_target, constructor)
  !< Create an instance of concrete extension of [[interpolations_object_constructor]].
  character(*),                                          intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                          intent(in)           :: S                 !< Stencils dimension.
  real(RPP),                                             intent(in)           :: stencil(1-S:)     !< Stencil used for inter, [1-S:-1+S].
  real(RPP),                                             intent(in)           :: x_target          !< Coordinate of the interp point.
  class(interpolations_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  type(interpolations_factory)                                                :: factory           !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type, &
                                  S=S,                                 &
                                  stencil=stencil,                     &
                                  x_target=x_target,                   &
                                  constructor=constructor)
  endsubroutine create_interpolations_int_object_constructor

  subroutine create_interpolator_object_constructor(interpolator_type, S, interpolations_constructor, weights_constructor, &
                                                    constructor)
  !< Create an instance of concrete extension of [[interpolator_object_constructor]].
  character(*),                                        intent(in)           :: interpolator_type          !< Type of interpolator.
  integer(I_P),                                        intent(in)           :: S                          !< Stencils dimension.
  class(interpolations_object_constructor),            intent(in)           :: interpolations_constructor !< Interpolations const.
  class(weights_object_constructor),                   intent(in)           :: weights_constructor        !< Weights constructor.
  class(interpolator_object_constructor), allocatable, intent(out)          :: constructor                !< Constructor.
  type(interpolator_factory)                                                :: factory                    !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type,                   &
                                  S=S,                                                   &
                                  interpolations_constructor=interpolations_constructor, &
                                  weights_constructor=weights_constructor,               &
                                  constructor=constructor)
  endsubroutine create_interpolator_object_constructor

  subroutine create_weights_object_constructor(interpolator_type, S, alpha_constructor, beta_constructor, kappa_constructor, &
                                               constructor)
  !< Create an instance of concrete extension of [[weights_object_constructor]].
  character(*),                                   intent(in)           :: interpolator_type !< Type of the interpolator.
  integer(I_P),                                   intent(in)           :: S                 !< Stencils dimension.
  class(alpha_object_constructor),                intent(in)           :: alpha_constructor !< Alpha constructor.
  class(beta_object_constructor),                 intent(in)           :: beta_constructor  !< Beta constructor.
  class(kappa_object_constructor),                intent(in)           :: kappa_constructor !< kappa constructor.
  class(weights_object_constructor), allocatable, intent(out)          :: constructor       !< Constructor.
  type(weights_factory)                                                :: factory           !< The factory.

  call factory%create_constructor(interpolator_type=interpolator_type, &
                                  S=S,                                 &
                                  alpha_constructor=alpha_constructor, &
                                  beta_constructor=beta_constructor,   &
                                  kappa_constructor=kappa_constructor, &
                                  constructor=constructor)
  endsubroutine create_weights_object_constructor
endmodule wenoof_objects_factory
