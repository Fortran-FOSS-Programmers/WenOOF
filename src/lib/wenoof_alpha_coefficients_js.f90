!< Jiang-Shu alpha coefficients object.
module wenoof_alpha_coefficients_js
!< Jiang-Shu alpha coefficients object.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.

use penf, only : I_P, R_P
use wenoof_alpha_coefficients

implicit none
private
public :: alpha_coefficients_js
public :: alpha_coefficients_js_constructor

type, extends(alpha_coefficients_constructor) :: alpha_coefficients_js_constructor
  !< Jiang-Shu alpha coefficient object constructor.
endtype alpha_coefficients_js_constructor

interface  alpha_coefficients_js_constructor
  procedure alpha_coefficients_js_constructor_
endinterface

type, extends(alpha_coefficients) :: alpha_coefficients_js
  !< Jiang-Shu alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted
  !< ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  private
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute alpha coefficients.
    procedure, nopass     :: description !< Return alpha coefficients string-description.
endtype alpha_coefficients_js

contains
  ! function-constructor
  function alpha_coefficients_js_constructor_(S) result(constructor)
  !< Return an instance of [alpha_coefficients_js_constructor].
  integer(I_P), intent(in)                           :: S           !< Maximum stencils dimension.
  class(alpha_coefficients_constructor), allocatable :: constructor !< Alpha coefficients constructor.

  allocate(alpha_coefficients_js_constructor :: constructor)
  constructor%S = S
  endfunction alpha_coefficients_js_constructor_

  ! deferred public methods
  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute alpha coefficients.
  class(alpha_coefficient_js), intent(inout) :: self                       !< Alpha coefficient.
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

  pure function description() result(string)
  !< Return alpha coefficients string-descripition.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line character.

  string = 'WENO alpha coefficient'//nl
  string = string//'  Based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    alpha(S, weigt_opt, IS, eps, f1, f2)'//nl
  string = string//'  where:'//nl
  string = string//'    S: integer(I_P), intent(in), the number of the stencils used'//nl
  string = string//'    weight_opt: real(R_P), intent(in), the optimal weight of the actual stencil'//nl
  string = string//'    IS: real(R_P), intent(in), the smoothness indicator of the actual stencil'//nl
  string = string//'    eps: real(R_P), intent(in), the coefficient to avoid zero division used'//nl
  string = string//'    f1, f2: integer(I_P), intent(in), the faces to be computed (1 => left interface, 2 => right interface)'
  endfunction description
endmodule wenoof_alpha_coefficient_js
