!< Abstract interpolator object.
module wenoof_interpolator_object
!< Abstract interpolator object.

use penf, only : I_P, R_P
use wenoof_base_object
use wenoof_interpolations_obejct
use wenoof_objects_factory
use wenoof_weights_object

implicit none
private
public :: interpolator_object
public :: interpolator_object_constructor

real(R_P), parameter :: EPS_DEF=10._R_P**(-6) !< Small epsilon to avoid division be zero, default value.

type, extends(base_object_constructor) :: interpolator_object_constructor
  !< Abstract interpolator object constructor.
  !<
  !< @note Every concrete WENO interpolator implementations must define their own constructor type.
  real(R_P)                                             :: eps=EPS_DEF                !< Small epsilon to avoid division by zero.
  class(interpolations_object_constructor), allocatable :: interpolations_constructor !< Stencil interpolations constructor.
  class(weights_object_constructor),        allocatable :: weights_constructor        !< Weights of interpolations constructor.
endtype interpolator_object_constructor

type, extends(base_object) :: interpolator_object
  !< Abstract interpolator object.
  !<
  !< @note Do not implement any actual interpolator: provide the interface for the different interpolators implemented.
  real(R_P)                                 :: eps=EPS_DEF    !< Small epsilon to avoid division by zero.
  class(interpolations_object), allocatable :: interpolations !< Stencil interpolations.
  class(weights_object),        allocatable :: weights        !< Weights of interpolations.
  contains
    ! public deferred methods
    procedure, pass(self) :: description          !< Return interpolator string-description.
    procedure, pass(self) :: interpolate_standard !< Interpolate values (without providing debug values).
    procedure, pass(self) :: interpolate_debug    !< Interpolate values (providing also debug values).
    ! public methods
    generic               :: interpolate => interpolate_standard, interpolate_debug !< Interpolate values.
    procedure, pass(self) :: create                                                 !< Create interpolator.
    procedure, pass(self) :: destroy                                                !< Destroy interpolator.
endtype interpolator_object

contains
  ! constructor methods

  ! public methods
  subroutine create_interpolator_constructor(self, S, interpolations_constructor, weights_constructor, eps)
  !< Create interpolator constructor.
  class(interpolator_object_constructor),   intent(inout)        :: self                       !< Interpolator constructor.
  integer(I_P),                             intent(in)           :: S                          !< Stencils dimension.
  class(interpolations_constructor_object), intent(in)           :: interpolations_constructor !< Interpolations constructor.
  class(weights_constructor_object),        intent(in)           :: weights_constructor        !< Weights constructor.
  real(R_P),                                intent(in), optional :: eps                        !< Small epsilon to avoid / by zero.

  call self%destroy
  self%S = S
  self%eps = EPS_DEF ; if present(eps) self%eps = eps
  allocate(self%interpolations_constructor, source=interpolations_constructor)
  allocate(self%weights_constructor       , source=weights_constructor       )
  endsubroutine create_interpolator_constructor

  pure subroutine destroy_interpolator_constructor(self)
  !< Destroy interpolator constructor.
  class(interpolator_object_constructor), intent(inout) :: self !< Interpolator constructor.

  call self%base_object_constructor%destroy
  self%eps = EPS_DEF
  if (allocated(self%interpolations_constructor)) deallocate(self%interpolations_constructor)
  if (allocated(self%weights_constructor       )) deallocate(self%weights_constructor       )
  endsubroutine destroy_interpolator_constructor

  ! interpolator methods

  ! deferred public methods
  pure function description(self) result(string)
  !< Return interpolator string-description.
  class(interpolator_object), intent(in) :: self   !< The interpolator.
  character(len=:), allocatable          :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator%description to be implemented by your concrete interpolator object'
#endif
  endfunction description

  pure subroutine interpolate_standard(self, stencil, interpolation)
  !< Interpolate values (without providing debug values).
  class(interpolator_object), intent(inout) :: self                     !< Interpolator.
  real(R_P),                  intent(in)    :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),                  intent(out)   :: interpolation(1:)        !< Result of the interpolation, [1:2].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator%interpolate_standard to be implemented by your concrete interpolator object'
#endif
  endsubroutine interpolate_standard

  pure subroutine interpolate_debug(self, stencil, interpolation, si, weights)
  !< Interpolate values (providing also debug values).
  class(interpolator_object), intent(inout) :: self                     !< Interpolator.
  real(R_P),                  intent(in)    :: stencil(1:, 1 - self%S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  real(R_P),                  intent(out)   :: interpolation(1:)        !< Result of the interpolation, [1:2].
  real(R_P),                  intent(out)   :: si(1:, 0:)               !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(R_P),                  intent(out)   :: weights(1:, 0:)          !< Weights of the stencils, [1:2, 0:S-1].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator%interpolate_debug to be implemented by your concrete interpolator object'
#endif
  endsubroutine interpolate_debug

  ! public methods
  subroutine create(self, constructor)
  !< Create interpolator.
  class(interpolator_object),     intent(inout) :: self        !< Interpolator.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.
  type(objects_factory)                         :: factory     !< Objects factory.

  call self%destroy
  call self%base_object%create(constructor=constructor)
  select type(constructor)
  class is(interpolator_constructor)
    self%eps = constructors%eps
    call factory%create(constructor=constructor%interpolations_constructor, object=self%interpolations)
    call factory%create(constructor=constructor%weights_constructor, object=self%weights)
  class default
    ! @TODO add error handling
  endselect
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy interpolator
  class(interpolator_object), intent(inout) :: self !< Interpolator.

  call self%base_object%destroy
  self%S = EPS_DEF
  if (allocated(self%interpolations)) deallocate(self%interpolations)
  if (allocated(self%weights)) deallocate(self%weights)
  endsubroutine destroy
endmodule wenoof_interpolator_object
