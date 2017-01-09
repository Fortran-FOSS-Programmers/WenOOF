!< Jiang-Shu WENO upwind interpolator object.
module wenoof_interpolator_js
!< Jiang-Shu WENO upwind interpolator object.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P, str
use wenoof_interpolator_abstract
use wenoof_smoothness_indicators_abstract
use wenoof_alpha_coefficient_abstract
use wenoof_optimal_weights_abstract
use wenoof_polynomials_abstract
use wenoof_alpha_coefficient_m
use wenoof_alpha_coefficient_z
use wenoof_alpha_coefficient_js
use wenoof_optimal_weights_js
use wenoof_smoothness_indicators_js
use wenoof_polynomials_js

implicit none
private
public :: wenoof_constructor_upwind
public :: wenoof_interpolator_upwind

type, extends(wenoof_constructor) :: wenoof_constructor_upwind
  !< Jiang-Shu WENO upwind interpolator object constructor.
  !<
  !< @note The constructed WENO interpolator implements the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  integer(I_P) :: S = 0               !< Stencils dimension.
  real(R_P)    :: eps = 10._R_P**(-6) !< Parameter for avoiding divided by zero when computing smoothness indicators.
endtype wenoof_constructor_upwind

interface wenoof_constructor_upwind
  procedure wenoof_constructor_upwind_init
endinterface

type, extends(wenoof_interpolator) :: wenoof_interpolator_upwind
  !< Jiang-Shu WENO upwind interpolator object.
  !<
  !< @note The WENO interpolator implemented is the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  !<
  !< @note The supported accuracy formal order are: 3rd, 5th, 7th, 9th, 11th, 13th, 15th, 17th  corresponding to use 2, 3, 4, 5, 6,
  !< 7, 8, 9 stencils composed of 2, 3, 4, 5, 6, 7, 8, 9 values, respectively.
  private
  integer(I_P) :: S = 0_I_P    !< Stencil dimension.
  real(R_P)    :: eps = 0._R_P !< Parameter for avoiding divided by zero when computing IS.
  contains
    ! public methods
    procedure, pass(self) :: create      !< Create interpolator.
    procedure, pass(self) :: description !< Return string-description of interpolator.
    procedure, pass(self) :: destroy     !< Destroy interpolator.
    procedure, pass(self) :: init_error  !< Return error code for wrong intialization.
    procedure, pass(self) :: interpolate !< Interpolate values.
endtype wenoof_interpolator_upwind
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! constructor_upwind
  elemental function wenoof_constructor_upwind_init(S, eps) result(constructor)
  !< Create (initialize) the WENO interpolator.
  !<
  !< @note For this class of interpolators it is sufficient to provide the maximum number of stencils used (that is also the
  !< dimension, i.e. number of values, of each stencil). During the actual interpolation phase the client code can specify, for each
  !< intepolation a different number of stencil bounded by this maximum value. This is useful for coupling the interpolator with
  !< algorithm like the Recursive Order Reduction (ROR) strategy.
  integer(I_P), intent(IN)           :: S           !< Maximum stencils dimension.
  real(R_P),    intent(IN), optional :: eps         !< Parameter for avoiding divided by zero when computing smoothness indicators.
  type(wenoof_constructor_upwind)    :: constructor !<WENO constructor.

  constructor%S = S
  if (present(eps)) constructor%eps = eps
  endfunction  wenoof_constructor_upwind_init

  ! interpolator_upwind
  ! public methods
  subroutine create(self, constructor, IS_type, alpha_type, weights_opt_type, polynomial_type)
  !< Create interpolator.
  class(wenoof_interpolator_upwind), intent(inout) :: self             !< WENO interpolator.
  class(wenoof_constructor),         intent(in)    :: constructor      !< WENO constructor.
  class(IS),                         intent(in)    :: IS_type          !< The concrete WENO smoothness indicator.
  class(alpha_coefficient),          intent(in)    :: alpha_type       !< The concrete WENO alpha coefficient.
  class(optimal_weights),            intent(in)    :: weights_opt_type !< The concrete WENO optimal weights.
  class(polynomials),                intent(in)    :: polynomial_type  !< The concrete WENO polynomial.

  select type(constructor)
  type is(wenoof_constructor_upwind)
    call self%destroy
    allocate(self%IS,      source=IS_type)
    allocate(self%alpha,   source=alpha_type)
    allocate(self%weights, source=weights_opt_type)
    allocate(self%polynom, source=polynomial_type)
    self%S = constructor%S
    self%eps = constructor%eps
    call self%IS%create(S=self%S)
    select type(alpha_type)
    type is(alpha_coefficient_js)
      call self%alpha%alloc(S=self%S)
    type is(alpha_coefficient_z)
      call self%alpha%alloc(S=self%S)
    type is(alpha_coefficient_m)
      call self%alpha%alloc(S=self%S)
    endselect
    call self%weights%create(S=self%S)
    call self%polynom%create(S=self%S)
  endselect
  endsubroutine create

  pure subroutine description(self, string)
  !< Return string-description of interpolator.
  !<
  !< @TODO make a function.
  class(wenoof_interpolator_upwind), intent(in)  :: self             !< WENO interpolator.
  character(len=:), allocatable,     intent(out) :: string           !< String returned.
  character(len=:), allocatable                  :: dummy_string     !< Dummy string.
  character(len=1), parameter                    :: nl=new_line('a') !< New line character.

  string = 'Jiang-Shu WENO upwind-biased interpolator'//nl
  string = string//'  Based on the scheme proposed by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  string = string//'  Provide a formal order of accuracy equals to: '//trim(str(2*self%S - 1, .true.))//nl
  string = string//'  Use '//trim(str(self%S, .true.))//' stencils composed by '//trim(str(self%S, .true.))//' values'//nl
  string = string//'  The eps value used for avoiding division by zero is '//trim(str(self%eps, .true.))//nl
  string = string//'  The "interpolate" method has the following public API'//nl
  string = string//'    interpolate(S, stencil, location, interpolation)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(in), the number of stencils actually used'//nl
  string = string//'    stencil(1:, 1-S:-1+S): real(R_P), intent(in), the stencils used'//nl
  string = string//'    location: character(*), intent(in), the location of interpolation {left, right, both}'//nl
  string = string//'    interpolation(1:, 1-S:-1+S): realR_P, intent(out), the interpolated values'//nl
  string = string//'  The alpha coefficient are evaluated by the following method'//nl
  call self%alpha%description(string=dummy_string)
  string = string//dummy_string//nl
  string = string//'  The smoothness indicators are evaluated by the following method'//nl
  call self%IS%description(string=dummy_string)
  string = string//dummy_string//nl
  string = string//'  The polynomials are evaluated by the following method'//nl
  call self%polynom%description(string=dummy_string)
  string = string//dummy_string//nl
  string = string//'  The optimal weights are evaluated by the following method'//nl
  call self%weights%description(string=dummy_string)
  string = string//dummy_string
  endsubroutine description

  elemental subroutine destroy(self)
  !< Destoy interpolator.
  class(wenoof_interpolator_upwind), intent(inout) :: self !< WENO interpolator.

  self%S = 0_I_P
  self%eps = 0._R_P
  !< Destroy WENO smoothness indicators object.
  if (allocated(self%IS)) then
    call self%IS%destroy
    deallocate (self%IS)
  endif
  !< Destroy WENO alpha object.
  if (allocated(self%alpha)) deallocate(self%alpha)
  !< Destroy WENO optimal weights object.
  if (allocated(self%weights)) then
    call self%weights%destroy
    deallocate(self%weights)
  endif
  !< Destroy WENO polynomials object.
  if (allocated(self%polynom)) then
    call self%polynom%destroy
    deallocate(self%polynom)
  endif
  endsubroutine destroy

  subroutine init_error(self, error_code)
  !< Return a string describing the WENO interpolator upwind.
  !<
  !< @TODO refactor.
  class(wenoof_interpolator_upwind), intent(inout) :: self             !< WENO interpolator.
  integer(I_P),                      intent(in)    :: error_code       !< Error code.
  character(len=:), allocatable                    :: string           !< Printed string.
  character(len=1), parameter                      :: nl=new_line('a') !< New line character.

  select case(error_code)
    case(0_I_P)
      string = 'The alpha base type is not present'//nl
      string = string//'  Please choose a valid alpha type for the base calculus of alphas before the WENO M procedure'//nl
      string = string//'  The available alpha base types are:'//nl
      string = string//'  "JS" for Jiang and Shu alpha evaluation and "Z" for WENO Z alpha evaluation'
    case(1_I_P)
      string = 'The alpha base type selected is invalid'//nl
      string = string//'  Please choose a valid alpha type for the base calculus of alphs before the WENO M procedure'//nl
      string = string//'  The available alpha base types are:'//nl
      string = string//'  "JS" for Jiang and Shu alpha evaluation and "Z" for WENO Z alpha evaluation'
    case(2_I_P)
      string = 'The alpha type selected is invalid'//nl
      string = string//'  Please choose a valid alpha type for the calculus of alpha coefficients'//nl
      string = string//'  The available alpha types are:'//nl
      string = string//'  "JS" for Jiang and Shu alpha evaluation, "Z" for WENO Z alpha evaluation'//nl
      string = string//'  and "M" for WENO M alpha evaluation'
    case(3_I_P)
      string = 'The smoothness indicators type selected is invalid'//nl
      string = string//'  Please choose a valid smoothness indicators type for the calculus of smoothness indicators'//nl
      string = string//'  The available alpha types are:'//nl
      string = string//'  "JS" for Jiang and Shu smoothness indicators'
    case(4_I_P)
      string = 'The optimal weights type selected is invalid'//nl
      string = string//'  Please choose a valid optimal weights type'//nl
      string = string//'  The available optimal weights are:'//nl
      string = string//'  "JS" for Jiang and Shu optimal weights'
    case(5_I_P)
      string = 'The polynomials type selected is invalid'//nl
      string = string//'  Please choose a valid polynomials type'//nl
      string = string//'  The available polynomials are:'//nl
      string = string//'  "JS" for Jiang and Shu polynomials'
  endselect
  write(stderr, '(A)') string
  call self%destroy
  endsubroutine init_error

  pure subroutine interpolate(self, S, stencil, location, interpolation)
  !< Interpolate values.
  class(wenoof_interpolator_upwind), intent(inout) :: self                  !< WENO interpolator.
  integer(I_P),                      intent(in)    :: S                     !< Number of stencils actually used.
  real(R_P),                         intent(in)    :: stencil(1:, 1 - S:)   !< Stencil of the interpolation [1:2, 1-S:-1+S].
  character(*),                      intent(in)    :: location              !< Location of interpolation: left, right, both.
  real(R_P),                         intent(out)   :: interpolation(1:)     !< Result of the interpolation, [1:2].
  real(R_P)                                        :: weights(1:2, 0:S - 1) !< Weights of the stencils, [1:2, 0:S-1].
  integer(I_P)                                     :: f1, f2, ff            !< Faces to be computed.
  integer(I_P)                                     :: f, k                  !< Counters.

  select case(location)
  case('both', 'b')
    f1=1_I_P; f2=2_I_P; ff=0_I_P
  case('left', 'l')
    f1=1_I_P; f2=1_I_P; ff=0_I_P
  case('right', 'r')
    f1=2_I_P; f2=2_I_P; ff=-1_I_P
  endselect

  call self%polynom%compute(S=S, stencil=stencil, f1=f1, f2=f2, ff = ff)

  call self%IS%compute(S=S, stencil=stencil, f1=f1, f2=f2, ff = ff)

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
endmodule wenoof_interpolator_js
