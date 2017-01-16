!< Abstract interpolator object.
module wenoof_interpolator
!< Abstract interpolator object.

use penf, only : I_P, R_P
use wenoof_alpha_coefficients
use wenoof_base_object
use wenoof_objects_factory
use wenoof_optimal_weights
use wenoof_smoothness_indicators
use wenoof_polynomials

implicit none
private
public :: interpolator
public :: interpolator_constructor

type, extends(base_object_constructor) :: interpolator_constructor
  !< Abstract interpolator object constructor.
  !<
  !< @note Every concrete WENO interpolator implementations must define their own constructor type.
  class(smoothness_indicators_constructor), allocatable :: is      !< Smoothness indicators constructor.
  class(alpha_coefficients_constructor),    allocatable :: alpha   !< Alpha coefficients constructor.
  class(optimal_weights_constructor),       allocatable :: weights !< Optimal weights constructor.
  class(polynomials_constructor),           allocatable :: polynom !< Polynomilas constructor.
  contains
    ! public methods
    procedure, pass(self) :: create  => create_interpolator_constructor  !< Create interpolator constructor.
    procedure, pass(self) :: destroy => destroy_interpolator_constructor !< Destroy interpolator constructor.
endtype interpolator_constructor

type, extends(base_object) :: interpolator
  !< Abstract interpolator object.
  !<
  !< @note Do not implement any actual interpolator: provide the interface for the different interpolators implemented.
  class(smoothness_indicators), allocatable :: is      !< Smoothness indicators.
  class(alpha_coefficients),    allocatable :: alpha   !< Alpha coefficients.
  class(optimal_weights),       allocatable :: weights !< Optimal weights.
  class(polynomials),           allocatable :: polynom !< Polynomilas.
  contains
    ! public deferred methods
    procedure, nopass     :: description          !< Return interpolator string-description.
    procedure, pass(self) :: interpolate_standard !< Interpolate values (without providing debug values).
    procedure, pass(self) :: interpolate_debug    !< Interpolate values (providing also debug values).
    ! public methods
    generic               :: interpolate => interpolate_standard, interpolate_debug !< Interpolate values.
    procedure, pass(self) :: create                                                 !< Create interpolator.
    procedure, pass(self) :: destroy                                                !< Destroy interpolator.
endtype interpolator

contains
  ! constructor methods

  ! public methods
  subroutine create_interpolator_constructor(self, is, alpha, weights, polynom)
  !< Create interpolator constructor.
  class(interpolator_constructor),          intent(inout) :: self    !< Interpolator constructor.
  class(smoothness_indicators_constructor), intent(in)    :: is      !< Smoothness indicators constructor.
  class(alpha_coefficients_constructor),    intent(in)    :: alpha   !< Alpha coefficients constructor.
  class(optimal_weights_constructor),       intent(in)    :: weights !< Optimal weights constructor.
  class(polynomials_constructor),           intent(in)    :: polynom !< Polynomilas constructor.

  call self%destroy
  allocate(self%is,      source=is     )
  allocate(self%alpha,   source=alpha  )
  allocate(self%weights, source=weights)
  allocate(self%polynom, source=polynom)
  endsubroutine create_interpolator_constructor

  pure subroutine destroy_interpolator_constructor(self)
  !< Destroy interpolator constructor.
  class(interpolator_constructor), intent(inout) :: self !< Interpolator.

  if (allocated(self%is))      deallocate(self%is)
  if (allocated(self%alpha))   deallocate(self%alpha)
  if (allocated(self%weights)) deallocate(self%weights)
  if (allocated(self%polynom)) deallocate(self%polynom)
  endsubroutine destroy_interpolator_constructor

  ! interpolator methods

  ! deferred public methods
  pure function description() result(string)
  !< Return interpolator string-description.
  character(len=:), allocatable  :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator%description to be implemented by your concrete interpolator object'
#endif
  endfunction description

  pure subroutine interpolate_standard(self, S, stencil, location, interpolation)
  !< Interpolate values (without providing debug values).
  class(interpolator), intent(inout) :: self                  !< Interpolator.
  integer(I_P),        intent(in)    :: S                     !< Number of stencils actually used.
  real(R_P),           intent(in)    :: stencil(1:, 1 - S:)   !< Stencil of the interpolation [1:2, 1-S:-1+S].
  character(*),        intent(in)    :: location              !< Location of interpolation: left, right, both.
  real(R_P),           intent(out)   :: interpolation(1:)     !< Result of the interpolation, [1:2].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator%interpolate to be implemented by your concrete interpolator object'
#endif
  endsubroutine interpolate_standard

  pure subroutine interpolate_debug(self, S, stencil, location, interpolation, si, weights)
  !< Interpolate values (providing also debug values).
  class(interpolator), intent(inout) :: self                !< Interpolator.
  integer(I_P),        intent(in)    :: S                   !< Number of stencils actually used.
  real(R_P),           intent(in)    :: stencil(1:, 1 - S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  character(*),        intent(in)    :: location            !< Location of interpolation: left, right, both.
  real(R_P),           intent(out)   :: interpolation(1:)   !< Result of the interpolation, [1:2].
  real(R_P),           intent(out)   :: si(1:, 0:)          !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(R_P),           intent(out)   :: weights(1:, 0:)     !< Weights of the stencils, [1:2, 0:S-1].

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator%interpolate to be implemented by your concrete interpolator object'
#endif
  endsubroutine interpolate_debug

  ! public methods
  subroutine create(self, constructor)
  !< Create interpolator.
  class(interpolator),            intent(inout) :: self        !< Interpolator.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.
  type(objects_factory)                         :: factory     !< Objects factory.

  call self%destroy
  select type(constructor)
  class is(interpolator_constructor)
    call factory%create(constructor=constructor%is,      object=self%is)
    call factory%create(constructor=constructor%alpha,   object=self%alpha)
    call factory%create(constructor=constructor%weights, object=self%weights)
    call factory%create(constructor=constructor%polynom, object=self%polynom)
  class default
    ! @TODO add error handling
  endselect
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy interpolator
  class(interpolator), intent(inout) :: self !< Interpolator.

  if (allocated(self%is)) deallocate(self%is)
  if (allocated(self%alpha)) deallocate(self%alpha)
  if (allocated(self%weights)) deallocate(self%weights)
  if (allocated(self%polynom)) deallocate(self%polynom)
  endsubroutine destroy
endmodule wenoof_interpolator
