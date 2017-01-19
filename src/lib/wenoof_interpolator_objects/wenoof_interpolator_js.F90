!< Jiang-Shu (upwind) interpolator object.
module wenoof_interpolator_js
!< Jiang-Shu (upwind) interpolator object.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P, str
use wenoof_alpha_coefficients
use wenoof_alpha_coefficients_m
use wenoof_alpha_coefficients_z
use wenoof_alpha_coefficients_js
use wenoof_base_object
use wenoof_interpolator
use wenoof_optimal_weights
use wenoof_optimal_weights_js
use wenoof_polynomials
use wenoof_polynomials_js
use wenoof_smoothness_indicators
use wenoof_smoothness_indicators_js

implicit none
private
public :: interpolator_js
public :: interpolator_js_constructor
public :: create_interpolator_js_constructor

type, extends(interpolator_constructor) :: interpolator_js_constructor
  !< Jiang-Shu (upwind) interpolator object constructor.
  !<
  !< @note The constructed WENO interpolator implements the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  private
  integer(I_P) :: S = 0               !< Stencils dimension.
  real(R_P)    :: eps = 10._R_P**(-6) !< Parameter for avoiding division by zero.
endtype interpolator_js_constructor

type, extends(interpolator) :: interpolator_js
  !< Jiang-Shu (upwind) interpolator object.
  !<
  !< @note The WENO interpolator implemented is the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  !<
  !< @note The supported accuracy formal order are: 3rd, 5th, 7th, 9th, 11th, 13th, 15th, 17th  corresponding to use 2, 3, 4, 5, 6,
  !< 7, 8, 9 stencils composed of 2, 3, 4, 5, 6, 7, 8, 9 values, respectively.
  private
  integer(I_P) :: S = 0_I_P    !< Stencil dimension.
  real(R_P)    :: eps = 0._R_P !< Parameter for avoiding divisiion by zero.
  contains
    ! public deferred methods
    procedure, nopass     :: description          !< Return interpolator string-description.
    procedure, pass(self) :: interpolate_standard !< Interpolate values (without providing debug values).
    procedure, pass(self) :: interpolate_debug    !< Interpolate values (providing also debug values).
    ! public methods
    procedure, pass(self) :: create  !< Create interpolator.
    procedure, pass(self) :: destroy !< Destroy interpolator.
endtype interpolator_js

contains
  ! public non TBP
  subroutine create_interpolator_js_constructor(is, alpha, weights, polynom, S, constructor, eps)
  !< Return an instance of [interpolator_js_constructor].
  class(smoothness_indicators_constructor),     intent(in)           :: is          !< Smoothness indicators constructor.
  class(alpha_coefficients_constructor),        intent(in)           :: alpha       !< Alpha coefficients constructor.
  class(optimal_weights_constructor),           intent(in)           :: weights     !< Optimal weights constructor.
  class(polynomials_constructor),               intent(in)           :: polynom     !< Polynomilas constructor.
  integer(I_P),                                 intent(in)           :: S           !< Stencil dimension.
  class(interpolator_constructor), allocatable, intent(out)          :: constructor !< Interpolator constructor.
  real(R_P),                                    intent(in), optional :: eps         !< Parameter for avoiding division by zero.

  allocate(interpolator_js_constructor :: constructor)
  call constructor%create(is=is, alpha=alpha, weights=weights, polynom=polynom)
  select type(constructor)
  class is(interpolator_js_constructor)
    constructor%S = S
    if (present(eps)) constructor%eps = eps
  endselect
  endsubroutine create_interpolator_js_constructor

  ! public deferred methods
  pure function description() result(string)
  !< Return interpolator string-descripition.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line character.
  character(len=:), allocatable :: dummy_string     !< Dummy string.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator_js to be implemented'
#endif
  endfunction description

  pure subroutine interpolate_standard(self, S, stencil, location, interpolation)
  !< Interpolate values (without providing debug values).
  class(interpolator_js), intent(inout) :: self                  !< Interpolator.
  integer(I_P),           intent(in)    :: S                     !< Number of stencils actually used.
  real(R_P),              intent(in)    :: stencil(1:, 1 - S:)   !< Stencil of the interpolation [1:2, 1-S:-1+S].
  character(*),           intent(in)    :: location              !< Location of interpolation: left, right, both.
  real(R_P),              intent(out)   :: interpolation(1:)     !< Result of the interpolation, [1:2].
  real(R_P)                             :: weights(1:2, 0:S - 1) !< Weights of the stencils, [1:2, 0:S-1].
  integer(I_P)                          :: f1, f2, ff            !< Faces to be computed.
  integer(I_P)                          :: f, k                  !< Counters.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator_js to be implemented'
#endif
  endsubroutine interpolate_standard

  pure subroutine interpolate_debug(self, S, stencil, location, interpolation, si, weights)
  !< Interpolate values (providing also debug values).
  class(interpolator_js), intent(inout) :: self                !< Interpolator.
  integer(I_P),           intent(in)    :: S                   !< Number of stencils actually used.
  real(R_P),              intent(in)    :: stencil(1:, 1 - S:) !< Stencil of the interpolation [1:2, 1-S:-1+S].
  character(*),           intent(in)    :: location            !< Location of interpolation: left, right, both.
  real(R_P),              intent(out)   :: interpolation(1:)   !< Result of the interpolation, [1:2].
  real(R_P),              intent(out)   :: si(1:, 0:)          !< Computed values of smoothness indicators [1:2, 0:S-1].
  real(R_P),              intent(out)   :: weights(1:, 0:)     !< Weights of the stencils, [1:2, 0:S-1].
  integer(I_P)                          :: f1, f2, ff          !< Faces to be computed.
  integer(I_P)                          :: f, k                !< Counters.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'interpolator_js to be implemented'
#endif
  endsubroutine interpolate_debug

  ! overridden methods
  subroutine create(self, constructor)
  !< Create interpolator.
  class(interpolator_js),         intent(inout) :: self        !< Interpolator.
  class(base_object_constructor), intent(in)    :: constructor !< Constructor.

  call self%destroy
  call self%interpolator%create(constructor=constructor)
  select type(constructor)
  type is(interpolator_js_constructor)
    self%S = constructor%S
    self%eps = constructor%eps
  endselect
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destoy interpolator.
  class(interpolator_js), intent(inout) :: self !< Interpolator.

  call self%interpolator%destroy
  self%S = 0_I_P
  self%eps = 0._R_P
  endsubroutine destroy
endmodule wenoof_interpolator_js
