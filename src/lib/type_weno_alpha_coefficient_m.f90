module type_weno_alpha_coefficient_m
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Henrick alpha coefficient for WENO schemes.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
!< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers, JCP,
!< 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
use type_weno_alpha_coefficient
use type_weno_alpha_coefficient_js
use type_weno_alpha_coefficient_z
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_alpha_coefficient_m
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_alpha_coefficient_z) :: weno_alpha_coefficient_m
  !< Henrick WENO alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
  !< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
  !< JCP, 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023.
  class(weno_alpha_coefficient), allocatable :: alpha_base !< To be set into [[initialize]] method.
  private
  contains
    ! deferred public methods
    procedure, pass(self), public :: description
    procedure, pass(self), public :: compute
    ! public methods
    procedure, pass(self), public :: initialize
endtype weno_alpha_coefficient_m
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! deferred public methods
  pure subroutine description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO alpha coefficient.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_m), intent(in)  :: self   !< WENO alpha coefficient.
  character(len=:), allocatable,   intent(out) :: string !< String returned.
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

  pure function compute(self, S, weight_opt, IS, eps) result(alpha)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the alpha coefficient of the WENO interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_m), intent(in) :: self       !< WENO alpha coefficient.
  integer(I_P),                    intent(in) :: S          !< Number of stencils used.
  real(R_P),                       intent(in) :: weight_opt !< Optimal weight of the stencil.
  real(R_P),                       intent(in) :: IS         !< Smoothness indicator of the stencil.
  real(R_P),                       intent(in) :: eps        !< Parameter for avoiding divided by zero.
  real(R_P),                                  :: alpha      !< Alpha coefficient of the stencil.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  associate(alpha_base => self%alpha_base%compute)
    m_alpha = (alpha_base(S,weight_opt,IS,eps) * (weight_opt + weight_opt * weight_opt - 3._R_P * weight_opt *          &
               alpha_base(S,weight_opt,IS,eps) + alpha_base(S,weight_opt,IS,eps) * alpha_base(S,weight_opt,IS,eps) ) )/ &
              (weight_opt * weight_opt + alpha_base(S,weight_opt,IS,eps) * (1._R_P - 2._R_P * weight_opt))
  endassociate
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction compute

  ! public methods
  subroutine initialize(self, alpha_base)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create the alias for the base alpha coefficient function.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_alpha_coefficient_m), intent(inout) :: self
  character(*)                     intent(in)    :: alpha_base
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
