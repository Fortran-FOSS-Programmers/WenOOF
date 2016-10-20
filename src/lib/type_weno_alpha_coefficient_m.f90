module type_weno_alpha_coefficient_m
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Henrick alpha coefficient for WENO schemes.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
!< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers, JCP,
!< 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P
use type_weno_alpha_coefficient
use type_weno_alpha_coefficient_js
use type_weno_alpha_coefficient_z
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_alpha_coefficient_m, associate_WENO_alpha_m, initialize
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_alpha_coefficient_z) :: weno_alpha_coefficient_m
  !< Henrick WENO alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
  !< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
  !< JCP, 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023.
  class(weno_alpha_coefficient), allocatable :: alpha_base !< To be set into [[initialize]] method.
  contains
    ! deferred public methods
    procedure, nopass,     public :: description
    procedure, pass(self), public :: compute
    ! public methods
    procedure, pass(self), public :: initialize
endtype weno_alpha_coefficient_m
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public, non TBP
  function associate_WENO_alpha_m(alpha_input) result(alpha_pointer)
    !< Check the type of alpha coefficient passed as input and return a WENO M alpha coefficient associated to alpha coefficient.
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

  ! deferred public methods
  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy Henrick WENO alpha coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_m), intent(inout) :: self   !< WENO alpha coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%alpha_coef)) deallocate(self%alpha_coef)
  if (allocated(self%alpha_tot)) deallocate(self%alpha_tot)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure subroutine create(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO alpha coefficients.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_m), intent(inout) :: self       !< WENO alpha coefficients.
  integer(I_P),                    intent(in)    :: S          !< Number of stencils used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  allocate(self%alpha_coef(1:2, 0:S - 1))
  allocate(self%alpha_tot(1:2))
  self%alpha_coef(:,:) = 100000._R_P
  self%alpha_tot(:) = 0._R_P
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create

  pure subroutine description(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable,   intent(out) :: string !< String returned.
  character(len=1), parameter                  :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO alpha coefficient'//nl
  string = string//'  Based on the work by Henrick, Aslam and Powers "Mapped weighted essentially non-oscillatory schemes: '// &
                   'Achieving optimal order near critical points", see JCP, 2005, vol. 207, pp. 542--567, '// &
                   'doi:10.1006/jcph.1996.0130'//nl
  string = string//'  The "alpha_m" method has the following public API'//nl
  string = string//'    alpha_m(IS,eps)'//nl
  string = string//'  where:'//nl
  string = string//'    weight_opt: real(R_P), intent(IN), the optimal weight of the actual stencil'//nl
  string = string//'    IS: real(R_P), intent(IN), the smoothness indicator of the actual stencil'//nl
  string = string//'    eps: real(R_P), intent(IN), the coefficient to avoid zero division used'//nl
  string = string//'    alpha: realR_P, intent(OUT), the alpha value'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_m), intent(inout) :: self                         !< WENO alpha coefficient.
  integer(I_P),                    intent(in)    :: S                            !< Number of stencils used.
  real(R_P),                       intent(in)    :: weight_opt(1: 2, 0: S - 1)   !< Optimal weight of the stencil.
  real(R_P),                       intent(in)    :: IS(1: 2, 0: S - 1)           !< Smoothness indicators of the stencils.
  real(R_P),                       intent(in)    :: eps                          !< Parameter for avoiding divided by zero.
  integer(I_P),                    intent(in)    :: f1, f2                       !< Faces to be computed.
  integer(I_P)                                   :: f, s1                        !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%alpha_tot = 0._R_P
  call self%alpha_base%compute(S=S, weight_opt=weight_opt, IS=IS, eps=eps, f1=f1, f2=f2)
  do s1 = 0, S - 1 ! stencil loops
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%alpha_coef(f, s1) = (self%alpha_base%alpha_coef(f, s1) * (weight_opt(f, s1) + weight_opt(f, s1) * weight_opt(f, s1) - &
                          3._R_P * weight_opt(f, s1) * self%alpha_base%alpha_coef(f, s1) + self%alpha_base%alpha_coef(f, s1) * &
                          self%alpha_base%alpha_coef(f, s1))) / &
                          (weight_opt(f, s1) * weight_opt(f, s1) + self%alpha_base%alpha_coef(f, s1) * &
                          (1._R_P - 2._R_P * weight_opt(f, s1)))
      self%alpha_tot(f) = self%alpha_tot(f) + self%alpha_coef(f, s1)
    enddo
  enddo
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute

  ! public methods
  subroutine initialize(self, alpha_base)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the alias for the base alpha coefficient function.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_m), intent(inout) :: self
  character(*),                    intent(in)    :: alpha_base
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(alpha_base)
  case('JS')
    if (allocated(self%alpha_base)) deallocate(self%alpha_base)
    allocate(weno_alpha_coefficient_js :: self%alpha_base)
  case('Z')
    if (allocated(self%alpha_base)) deallocate(self%alpha_base)
    allocate(weno_alpha_coefficient_z :: self%alpha_base)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine initialize
endmodule type_weno_alpha_coefficient_m
