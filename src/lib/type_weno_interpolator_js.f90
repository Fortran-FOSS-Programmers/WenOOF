module type_weno_interpolator_js
!-----------------------------------------------------------------------------------------------------------------------------------
!< Concrete WENO Jiang-Shu upwind interpolator object.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P, str

use type_weno_interpolator
use type_weno_smoothness_indicators
use type_weno_alpha_coefficient
use type_weno_optimal_weights
use type_weno_polynomials

use type_weno_alpha_coefficient_m
use type_weno_alpha_coefficient_z
use type_weno_alpha_coefficient_js
use type_weno_optimal_weights_js
use type_weno_smoothness_indicators_js
use type_weno_polynomials_js
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_interpolator, weno_constructor
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_constructor) :: weno_constructor_upwind
  !< Upwind biased WENO interpolator constructor,
  !<
  !< @note The constructed WENO interpolator implements the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  integer(I_P) :: S = 0               !< Stencils dimension.
  real(R_P)    :: eps = 10._R_P**(-6) !< Parameter for avoiding divided by zero when computing smoothness indicators.
endtype weno_constructor_upwind

interface weno_constructor_upwind
  procedure weno_constructor_upwind_init
endinterface

type, extends(weno_interpolator) :: weno_interpolator_upwind
  !< Upwind biased WENO interpolator object,
  !<
  !< @note The WENO interpolator implemented is the *Efficient Implementation of Weighted ENO Schemes*,
  !< Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  !<
  !< @note The supported accuracy formal order are: 3rd, 5th, 7th, 9th, 11th, 13th corresponding to use 2, 3, 4, 5, 6, 7 stencils
  !< composed of 2, 3, 4, 5, 6, 7 values, respectively.
  private
  integer(I_P)                               :: S = 0_I_P    !< Stencil dimension.
  real(R_P)                                  :: eps = 0._R_P !< Parameter for avoiding divided by zero when computing IS.
  contains
    ! public methods
    procedure, pass(self), public :: destroy
    procedure, pass(self), public :: create
    procedure, pass(self), public :: description
    procedure, pass(self), public :: init_error
    procedure, pass(self), public :: interpolate
    ! private methods
    final :: finalize
endtype weno_interpolator_upwind

interface associate_WENO_IS
  module procedure associate_WENO_IS_js
end interface

interface associate_WENO_alpha
  module procedure associate_WENO_alpha_js, &
                   associate_WENO_alpha_z,  &
                   associate_WENO_alpha_m
end interface

interface associate_WENO_weights
  module procedure associate_WENO_weights_js
end interface

interface associate_WENO_polynomials
  module procedure associate_WENO_polynomials_js
end interface
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! weno_constructor_upwind
  elemental function  weno_constructor_upwind_init(S, eps) result(constructor)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create (initialize) the WENO interpolator.
  !<
  !< @note For this class of interpolators it is sufficient to provide the maximum number of stencils used (that is also the
  !< dimension, i.e. number of values, of each stencil). During the actual interpolation phase the client code can specify, for each
  !< intepolation a different number of stencil bounded by this maximum value. This is useful for coupling the interpolator with
  !< algorithm like the Recursive Order Reduction (ROR) strategy.
  !---------------------------------------------------------------------------------------------------------------------------------
  integer(I_P), intent(IN)           :: S           !< Maximum stencils dimension.
  real(R_P),    intent(IN), optional :: eps         !< Parameter for avoiding divided by zero when computing smoothness indicators.
  type(weno_constructor_upwind)      :: constructor !<WENO constructor.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  constructor%S = S
  if (present(eps)) constructor%eps = eps
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction  weno_constructor_upwind_init

  ! weno_interpolator_upwind
  ! public methods
  function associate_WENO_IS_js(IS_input) result(IS_pointer)
    !< Check the type of the smoothness indicator passed as input and return a Jiang-Shu smoothness indicator associated to the smoothness indicator.
    class(weno_IS), intent(in), target  :: IS_input   !< Input smoothness indicator.
    class(weno_IS_js),          pointer :: IS_pointer !< Jiang Shu smoothness indicator.

    select type(IS_input)
      type is(weno_IS_js)
        IS_pointer => IS_input
      class default
        write(stderr, '(A)')'error: wrong smoothness indicator type chosen'
        stop
    end select
  end function associate_WENO_IS_js

  function associate_WENO_alpha_js(alpha_input) result(alpha_pointer)
    !< Check the type of the alpha coefficient passed as input and return a Jiang-Shu alpha coefficient associated to the alpha coefficient.
    class(weno_alpha_coefficient), intent(in), target  :: alpha_input   !< Input alpha coefficient.
    class(weno_alpha_coefficient_js),          pointer :: alpha_pointer !< Jiang Shu alpha coefficients.

    select type(alpha_input)
      type is(weno_alpha_coefficient_js)
        alpha_pointer => alpha_input
      class default
        write(stderr, '(A)')'error: wrong alpha coefficient type chosen'
        stop
    end select
  end function associate_WENO_alpha_js

  function associate_WENO_alpha_z(alpha_input) result(alpha_pointer)
    !< Check the type of the alpha coefficient passed as input and return a WENO Z alpha coefficient associated to the alpha coefficient.
    class(weno_alpha_coefficient), intent(in), target  :: alpha_input   !< Input alpha coefficient.
    class(weno_alpha_coefficient_z),           pointer :: alpha_pointer !< WENO Z alpha coefficients.

    select type(alpha_input)
      type is(weno_alpha_coefficient_z)
        alpha_pointer => alpha_input
      class default
        write(stderr, '(A)')'error: wrong alpha coefficient type chosen'
        stop
    end select
  end function associate_WENO_alpha_z

  function associate_WENO_alpha_m(alpha_input) result(alpha_pointer)
    !< Check the type of the alpha coefficient passed as input and return a WENO M alpha coefficient associated to the alpha coefficient.
    class(weno_alpha_coefficient), intent(in), target  :: alpha_input   !< Input alpha coefficient.
    class(weno_alpha_coefficient_m),           pointer :: alpha_pointer !< WENO M alpha coefficients.

    select type(alpha_input)
      type is(weno_alpha_coefficient_m)
        alpha_pointer => alpha_input
      class default
        write(stderr, '(A)')'error: wrong alpha coefficient type chosen'
        stop
    end select
  end function associate_WENO_alpha_m

  function associate_WENO_weights_js(weights_input) result(weights_pointer)
    !< Check the type of the optimal weights passed as input and return Jiang-Shu optimal weights associated to the optimal weights.
    class(weno_optimal_weights), intent(in), target  :: weights_input   !< Input optimal weights.
    class(weno_optimal_weights_js),          pointer :: weights_pointer !< Jiang Shu optimal weights.

    select type(weights_input)
      type is(weno_optimal_weights_js)
        weights_pointer => weights_input
      class default
        write(stderr, '(A)')'error: wrong optimal weights type chosen'
        stop
    end select
  end function associate_WENO_weights_js

  function associate_WENO_polynomials_js(polyn_input) result(polyn_pointer)
    !< Check the type of the polynomials passed as input and return Jiang-Shu polynomials associated to the polynomials.
    class(weno_polynomials), intent(in), target  :: polyn_input   !< Input optimal weights.
    class(weno_polynomials_js),          pointer :: polyn_pointer !< Jiang Shu optimal weights.

    select type(polyn_input)
      type is(weno_polynomials_js)
        polyn_pointer => polyn_input
      class default
        write(stderr, '(A)')'error: wrong polynomials type chosen'
        stop
    end select
  end function associate_WENO_polynomials_js

  elemental subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destoy the WENO interpolator upwind.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(inout) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%S = 0_I_P
  self%eps = 0._R_P
  !< Destroy WENO smoothness indicators object.
  if (associated(self%IS)) then
    call self%IS%destroy
    deallocate (self%IS)
    nullify (self%IS)
  endif
  !< Destroy WENO alpha object.
  if (associated(self%alpha)) deallocate(self%alpha)
  !< Destroy WENO optimal weights object.
  if (associated(self%weights)) then
    call self%weights%destroy
    deallocate(self%weights)
    nullify(self%weights)
  endif
  !< Destroy WENO polynomials object.
  if (associated(self%polynom)) then
    call self%polynom%destroy
    deallocate(self%polynom)
    nullify(self%polynom)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  subroutine create(self, constructor, IS_type, alpha_type, alpha_base_type, weights_opt_type, polynomial_type)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the WENO interpolator upwind.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(inout) :: self        !< WENO interpolator.
  class(weno_constructor),         intent(in)    :: constructor !< WENO constructor.
  character(*),             intent(in)               :: IS_type           !< The concrete WENO smoothness indicator.
  character(*),             intent(in)               :: alpha_type        !< The concrete WENO alpha coefficient.
  character(*),             intent(in), optional     :: alpha_base_type   !< The WENO alpha coefficient base for WENO Mapped.
  character(*),             intent(in)               :: weights_opt_type  !< The concrete WENO optimal weights.
  character(*),             intent(in)               :: polynomial_type   !< The concrete WENO polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select type(constructor)
  type is(weno_constructor_upwind)
    call self%destroy
    self%S = constructor%S
    self%eps = constructor%eps
    !< Create WENO smoothness indicators object.
    select case(IS_type)
    case('JS')
      self%IS => associate_WENO_IS(IS_input=weno_IS_js)
    endselect
    !< Create WENO alpha object.
    select case(alpha_type)
    case('JS')
      self%alpha => associate_WENO_alpha(alpha_input=weno_alpha_coefficient_js)
    case('Z')
      self%alpha => associate_WENO_alpha(alpha_input=weno_alpha_coefficient_z)
    case('M')
      self%alpha => associate_WENO_alpha(alpha_input=weno_alpha_coefficient_m)
    endselect
    !< Create WENO optimal weights object.
    select case(weights_opt_type)
    case('JS')
      self%weights => associate_WENO_weights(weights_input=weno_optimal_weights_js)
    endselect
    !< Create WENO polynomials object.
    select case(polynomial_type)
    case('JS')
      self%polynom => associate_WENO_polynom(polyn_input=weno_polynomials_js)
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create

  pure subroutine description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing the WENO interpolator upwind.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(in)  :: self              !< WENO interpolator.
  character(len=:), allocatable,   intent(out) :: string            !< String returned.
  character(len=:), allocatable                :: dummy_string      !< Dummy string.
  character(len=1), parameter                  :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO upwind-biased interpolator'//nl
  string = string//'  Based on the scheme proposed by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  string = string//'  Provide a formal order of accuracy equals to: '//trim(str(2*self%S - 1, .true.))//nl
  string = string//'  Use '//trim(str(self%S, .true.))//' stencils composed by '//trim(str(self%S, .true.))//' values'//nl
  string = string//'  The eps value used for avoiding division by zero is '//trim(str(self%eps, .true.))//nl
  string = string//'  The "interpolate" method has the following public API'//nl
  string = string//'    interpolate(S, stencil, location, interpolation)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(IN), the number of stencils actually used'//nl
  string = string//'    stencil(1:, 1-S:-1+S): real(R_P), intent(IN), the stencils used'//nl
  string = string//'    location: character(*), intent(IN), the location of interpolation {left, right, both}'//nl
  string = string//'    interpolation(1:, 1-S:-1+S): realR_P, intent(OUT), the interpolated values'//nl
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  subroutine init_error(self, error_code)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing the WENO interpolator upwind.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(inout) :: self              !< WENO interpolator.
  integer(I_P),                    intent(in)    :: error_code        !< Error code.
  character(len=:), allocatable                  :: string            !< Printed string.
  character(len=1), parameter                    :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
  !print *, string
  write(stderr, '(A)') string
  call self%destroy
  stop
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init_error

  pure subroutine interpolate(self, S, stencil, location, interpolation)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Interpolate the stencil input values computing the actual interpolation.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_interpolator_upwind), intent(in)  :: self                      !< WENO interpolator.
  integer,                         intent(in)  :: S                         !< Number of stencils actually used.
  real(R_P),                       intent(in)  :: stencil(1:, 1 - S:)       !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  character(*),                    intent(in)  :: location                  !< Location of interpolated value(s): left, right, both.
  real(R_P),                       intent(out) :: interpolation(1:)         !< Result of the interpolation, [1:2].
  real(R_P)                                    :: polynomials(1:2, 0:S - 1) !< Polynomial reconstructions.
  real(R_P)                                    :: weights(1:2, 0:S - 1)     !< Weights of the stencils, [1:2, 0:S-1 ].
  real(R_P)                                    :: a(1:2, 0:S - 1)           !< Alpha coefficients for the weights.
  real(R_P)                                    :: a_tot(1:2)                !< Sum of the alpha coefficients.
  real(R_P)                                    :: IS(1:2, 0:S - 1)          !< Smoothness indicators of the stencils, [1:2, 0:S-1].
  integer(I_P)                                 :: f1, f2, ff                !< Faces to be computed.
  integer(I_P)                                 :: s1, s2, s3, f, k          !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------
  select case(location)
  case('both', 'b')
    f1=1_I_P; f2=2_I_P; ff=0_I_P
  case('left', 'l')
    f1=1_I_P; f2=1_I_P; ff=0_I_P
  case('right', 'r')
    f1=2_I_P; f2=2_I_P; ff=-1_I_P
  endselect

  ! computing polynomials reconstructions.
  polynomials = 0.
  do s1 = 0, S - 1 ! stencils loop
    do s2 = 0, S - 1 ! values loop
      do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        polynomials(f, s1) = polynomials(f, s1) + &
                             self%polynom%compute(poly_coef=self%polynom%coef(f, s2, s1), v=stencil(f + ff, -s2 + s1))
      enddo
    enddo
  enddo

  ! computing smoothness indicators
  do s1 = 0, S - 1 ! stencils loop
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      IS(f, s1) = 0.
      do s2 = 0, S - 1
        do s3 = 0, S - 1
          IS(f, s1)=IS(f, s1) + &
                    self%IS%compute(smooth_coef=self%IS%coef(s3, s2, s1), v1=stencil(f + ff, s1 - s3), v2=stencil(f + ff, s1 - s2))
        enddo
      enddo
    enddo
  enddo

  ! computing alpha coefficients
  a_tot = 0.
  do s1 = 0, S - 1 ! stencil loops
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      a(f, s1) = self%alpha%compute(S=S, weight_opt=self%weight%w(f, s1), IS = IS(f, s1))
      a_tot(f) = a_tot(f) + a(f, s1)
    enddo
  enddo

  ! computing the weights
  do s1 = 0, S - 1 ! stencils loop
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      weights(f, s1) = a(f, s1) / a_tot(f)
    enddo
  enddo

  ! computing the convultion
  interpolation = 0.
  do k = 0, S - 1 ! stencils loop
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      interpolation(f + ff) = interpolation(f + ff) + weights(f, k) * polynomials(f, k)
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine interpolate

  elemental subroutine finalize(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Finalize object.
  !---------------------------------------------------------------------------------------------------------------------------------
  type(weno_interpolator_upwind), intent(INOUT) :: self !< WENO interpolator.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule type_weno_interpolator_js
