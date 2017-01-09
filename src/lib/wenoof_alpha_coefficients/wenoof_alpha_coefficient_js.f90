!< Jiang-Shu WENO alpha coefficient object.
module wenoof_alpha_coefficient_js
!< Jiang-Shu WENO alpha coefficient object.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.

use penf, only : I_P, R_P
use wenoof_alpha_coefficient_abstract

implicit none
private
public :: alpha_coefficient_js

type, extends(alpha_coefficient) :: alpha_coefficient_js
  !< Jiang-Shu WENO alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted
  !< ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  private
  contains
    ! deferred public methods
    procedure, pass(self) :: create      !< Create coefficients.
    procedure, pass(self) :: compute     !< Compute coefficients.
    procedure, nopass     :: description !< Return string-description of coefficients.
    procedure, pass(self) :: destroy     !< Destroy coefficients.
endtype alpha_coefficient_js
contains
  ! deferred public methods
  pure subroutine destroy(self)
  !< Destroy coefficients.
  class(alpha_coefficient_js), intent(inout) :: self !< WENO alpha coefficients.

  if (allocated(self%alpha_coef)) deallocate(self%alpha_coef)
  if (allocated(self%alpha_tot)) deallocate(self%alpha_tot)
  endsubroutine destroy

  pure subroutine create(self, S)
  !< Create coefficients.
  class(alpha_coefficient_js), intent(inout) :: self !< WENO alpha coefficients.
  integer(I_P),                intent(in)    :: S    !< Number of stencils used.

  call self%destroy
  allocate(self%alpha_coef(1:2, 0:S - 1))
  allocate(self%alpha_tot(1:2))
  self%alpha_coef(:,:) = 100000._R_P
  self%alpha_tot(:) = 0._R_P
  endsubroutine create

  pure subroutine description(string)
  !< Return string-descripition of coefficients.
  !<
  !< @TODO change to function.
  character(len=:), allocatable, intent(out) :: string           !< String returned.
  character(len=1), parameter                :: nl=new_line('a') !< New line character.

  string = 'WENO alpha coefficient'//nl
  string = string//'  Based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    alpha(S, weigt_opt, IS, eps)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(in), the number of the stencils used'//nl
  string = string//'    weight_opt: real(R_P), intent(in), the optimal weight of the actual stencil'//nl
  string = string//'    IS: real(R_P), intent(in), the smoothness indicator of the actual stencil'//nl
  string = string//'    eps: real(R_P), intent(in), the coefficient to avoid zero division used'
  endsubroutine description

  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute coefficients.
  class(alpha_coefficient_js), intent(inout) :: self                       !< WENO alpha coefficient.
  integer(I_P),                intent(in)    :: S                          !< Number of stencils used.
  real(R_P),                   intent(in)    :: weight_opt(1: 2, 0: S - 1) !< Optimal weight of the stencil.
  real(R_P),                   intent(in)    :: IS(1: 2, 0: S - 1)         !< Smoothness indicators of the stencils.
  real(R_P),                   intent(in)    :: eps                        !< Parameter for avoiding divided by zero.
  integer(I_P),                intent(in)    :: f1, f2                     !< Faces to be computed.
  integer(I_P)                               :: f, s1                      !< Counters.

  self%alpha_tot = 0._R_P
  do s1=0, S - 1 ! stencil loops
    do f=f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%alpha_coef(f, s1) = weight_opt(f, s1) * (1._R_P/(eps + IS(f, s1))**S)
      self%alpha_tot(f) = self%alpha_tot(f) + self%alpha_coef(f, s1)
    enddo
  enddo
  endsubroutine compute
endmodule wenoof_alpha_coefficient_js
