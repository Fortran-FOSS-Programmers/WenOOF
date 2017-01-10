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

type, extends(interpolator_constructor) :: interpolator_js_constructor
  !< Jiang-Shu (upwind) interpolator object constructor.
  !<
  !< @note The constructed WENO interpolator implements the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  private
  integer(I_P) :: S = 0               !< Stencils dimension.
  real(R_P)    :: eps = 10._R_P**(-6) !< Parameter for avoiding division by zero.
endtype interpolator_js_constructor

interface interpolator_js_constructor
  procedure interpolator_js_constructor_
endinterface

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
    procedure, nopass     :: description !< Return interpolator string-description.
    procedure, pass(self) :: interpolate !< Interpolate values.
    ! public methods
    procedure, pass(self) :: create  !< Create interpolator.
    procedure, pass(self) :: destroy !< Destroy interpolator.
endtype interpolator_js

contains
  ! function-constructor
  function interpolator_js_constructor_(is, alpha, weights, polynom, S, eps) result(constructor)
  !< Return an instance of [interpolator_js_constructor].
  class(base_object_constructor), intent(in)           :: is          !< Smoothness indicators constructor.
  class(base_object_constructor), intent(in)           :: alpha       !< Alpha coefficients constructor.
  class(base_object_constructor), intent(in)           :: weights     !< Optimal weights constructor.
  class(base_object_constructor), intent(in)           :: polynom     !< Polynomilas constructor.
  integer(I_P),                   intent(in)           :: S           !< Stencil dimension.
  real(R_P),                      intent(in), optional :: eps         !< Parameter for avoiding division by zero.
  class(interpolator_constructor), allocatable         :: constructor !< Interpolator constructor.

  allocate(interpolator_js_constructor :: constructor)
  call constructor%create(is=is, alpha=alpha, weights=weights, polynom=polynom)
  select type(constructor)
  class is(interpolator_js_constructor)
    constructor%S = S
    if (present(eps)) constructor%eps = eps
  endselect
  endfunction interpolator_js_constructor_

  ! public deferred methods
  pure function description() result(string)
  !< Return interpolator string-descripition.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line character.
  character(len=:), allocatable :: dummy_string     !< Dummy string.

  string = 'Jiang-Shu WENO upwind-biased interpolator'//nl
  string = string//'  Based on the scheme proposed by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  ! string = string//'  Provide a formal order of accuracy equals to: '//trim(str(2*self%S - 1, .true.))//nl
  ! string = string//'  Use '//trim(str(self%S, .true.))//' stencils composed by '//trim(str(self%S, .true.))//' values'//nl
  ! string = string//'  The eps value used for avoiding division by zero is '//trim(str(self%eps, .true.))//nl
  string = string//'  The "interpolate" method has the following public API'//nl
  string = string//'    interpolate(S, stencil, location, interpolation)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(in), the number of stencils actually used'//nl
  string = string//'    stencil(1:, 1-S:-1+S): real(R_P), intent(in), the stencils used'//nl
  string = string//'    location: character(*), intent(in), the location of interpolation {left, right, both}'//nl
  string = string//'    interpolation(1:, 1-S:-1+S): realR_P, intent(out), the interpolated values'//nl
  string = string//'  The alpha coefficient are evaluated by the following method'//nl
  ! string = string//self%alpha%description()//nl
  string = string//'  The smoothness indicators are evaluated by the following method'//nl
  ! string = string//self%IS%description()//nl
  string = string//'  The polynomials are evaluated by the following method'//nl
  ! string = string//self%polynom%description()//nl
  string = string//'  The optimal weights are evaluated by the following method'//nl
  ! string = string//self%weights%description()
  endfunction description

  pure subroutine interpolate(self, S, stencil, location, interpolation)
  !< Interpolate values.
  class(interpolator_js), intent(inout) :: self                  !< Interpolator.
  integer(I_P),           intent(in)    :: S                     !< Number of stencils actually used.
  real(R_P),              intent(in)    :: stencil(1:, 1 - S:)   !< Stencil of the interpolation [1:2, 1-S:-1+S].
  character(*),           intent(in)    :: location              !< Location of interpolation: left, right, both.
  real(R_P),              intent(out)   :: interpolation(1:)     !< Result of the interpolation, [1:2].
  real(R_P)                             :: weights(1:2, 0:S - 1) !< Weights of the stencils, [1:2, 0:S-1].
  integer(I_P)                          :: f1, f2, ff            !< Faces to be computed.
  integer(I_P)                          :: f, k                  !< Counters.

  select case(location)
  case('both', 'b')
    f1=1_I_P; f2=2_I_P; ff=0_I_P
  case('left', 'l')
    f1=1_I_P; f2=1_I_P; ff=0_I_P
  case('right', 'r')
    f1=2_I_P; f2=2_I_P; ff=-1_I_P
  endselect

  call self%polynom%compute(S=S, stencil=stencil, f1=f1, f2=f2, ff = ff)

  call self%is%compute(S=S, stencil=stencil, f1=f1, f2=f2, ff = ff)

  call self%alpha%compute(S=S, weight_opt=self%weights%opt, IS = self%IS%si, eps = self%eps, f1=f1, f2=f2)

  ! computing the weights
  do k = 0, S - 1 ! stencils loop
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      weights(f, k) = self%alpha%alpha_coef(f, k) / self%alpha%alpha_tot(f)
    enddo
  enddo

  ! computing the convolution
  interpolation = 0.
  do k = 0, S - 1 ! stencils loop
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      interpolation(f + ff) = interpolation(f + ff) + weights(f, k) * self%polynom%poly(f, k)
    enddo
  enddo
  endsubroutine interpolate

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
