!< Henrick WENO alpha coefficient object.
module wenoof_alpha_coefficient_m
!< Henrick WENO alpha coefficient object.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
!< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers, JCP,
!< 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023

use penf, only : I_P, R_P
use wenoof_alpha_coefficient_abstract
use wenoof_alpha_coefficient_js
use wenoof_alpha_coefficient_z

implicit none
private
public :: alpha_coefficient_m

type, extends(alpha_coefficient_z) :: alpha_coefficient_m
  !< Henrick WENO alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
  !< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
  !< JCP, 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023.
  class(alpha_coefficient), allocatable :: alpha_base !< To be set into [[alpha_coefficient_m:initialize]] method.
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Return string-description of coefficients.
    procedure, nopass     :: description !< Compute coefficients.
    ! public methods
    procedure, pass(self) :: initialize !< Initialize the base alpha coefficient function.
    ! overridden public methods
    procedure, pass(self) :: destroy !< Destroy coefficients.
endtype alpha_coefficient_m

contains
  ! deferred public methods
  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute coefficients.
  class(alpha_coefficient_m), intent(inout) :: self                       !< WENO alpha coefficient.
  integer(I_P),               intent(in)    :: S                          !< Number of stencils used.
  real(R_P),                  intent(in)    :: weight_opt(1: 2, 0: S - 1) !< Optimal weight of the stencil.
  real(R_P),                  intent(in)    :: IS(1: 2, 0: S - 1)         !< Smoothness indicators of the stencils.
  real(R_P),                  intent(in)    :: eps                        !< Parameter for avoiding divided by zero.
  integer(I_P),               intent(in)    :: f1, f2                     !< Faces to be computed.
  integer(I_P)                              :: f, s1                      !< Counters.

  self%alpha_tot = 0._R_P
  call self%alpha_base%compute(S=S, weight_opt=weight_opt, IS=IS, eps=eps, f1=f1, f2=f2)
  do s1=0, S - 1 ! stencil loops
    do f=f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%alpha_coef(f, s1) =                                                                                &
        (self%alpha_base%alpha_coef(f, s1) * (weight_opt(f, s1) + weight_opt(f, s1) * weight_opt(f, s1) -     &
         3._R_P * weight_opt(f, s1) * self%alpha_base%alpha_coef(f, s1) + self%alpha_base%alpha_coef(f, s1) * &
         self%alpha_base%alpha_coef(f, s1))) /                                                                &
         (weight_opt(f, s1) * weight_opt(f, s1) + self%alpha_base%alpha_coef(f, s1) *                         &
         (1._R_P - 2._R_P * weight_opt(f, s1)))
      self%alpha_tot(f) = self%alpha_tot(f) + self%alpha_coef(f, s1)
    enddo
  enddo
  endsubroutine compute

  pure subroutine description(string)
  !< Return string-descripition of coefficients.
  !<
  !< @TODO change to function.
  character(len=:), allocatable, intent(out) :: string           !< String returned.
  character(len=1), parameter                :: nl=new_line('a') !< New line character.

  string = 'WENO alpha coefficient'//nl
  string = string//'  Based on the work by Henrick, Aslam and Powers "Mapped weighted essentially non-oscillatory schemes: '// &
                   'Achieving optimal order near critical points", see JCP, 2005, vol. 207, pp. 542--567, '// &
                   'doi:10.1006/jcph.1996.0130'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    compute(S, weigt_opt, IS, eps, f1, f2)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(in), the number of the stencils used'//nl
  string = string//'    weight_opt: real(R_P), intent(in), the optimal weight of the actual stencil'//nl
  string = string//'    IS: real(R_P), intent(in), the smoothness indicator of the actual stencil'//nl
  string = string//'    eps: real(R_P), intent(in), the coefficient to avoid zero division used'//nl
  string = string//'    f1, f2: integer(I_P), intent(in), the faces to be computed (1 => left interface, 2 => right interface)'
  endsubroutine description

  ! public methods
  subroutine initialize(self, alpha_base)
  !< Initialize the base alpha coefficient function.
  class(alpha_coefficient_m), intent(inout) :: self       !< WENO alpha coefficient.
  character(*),               intent(in)    :: alpha_base !< Base alpha coefficient type.

  select case(alpha_base)
  case('JS')
    if (allocated(self%alpha_base)) deallocate(self%alpha_base)
    allocate(alpha_coefficient_js :: self%alpha_base)
  case('Z')
    if (allocated(self%alpha_base)) deallocate(self%alpha_base)
    allocate(alpha_coefficient_z :: self%alpha_base)
  endselect
  endsubroutine initialize

  ! overridden methods
  pure subroutine destroy(self)
  !< Destroy coefficients.
  class(alpha_coefficient_m), intent(inout) :: self !< WENO alpha coefficients.

  call self%alpha_coefficient%destroy
  if (allocated(self%alpha_base)) deallocate(self%alpha_base)
  endsubroutine destroy
endmodule wenoof_alpha_coefficient_m
