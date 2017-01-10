!< Henrick alpha coefficients object.
module wenoof_alpha_coefficients_m
!< Henrick alpha coefficients object.
!<
!< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
!< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers, JCP,
!< 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023

use penf, only : I_P, R_P
use wenoof_alpha_coefficients
use wenoof_alpha_coefficients_js
use wenoof_alpha_coefficients_z

implicit none
private
public :: alpha_coefficients_m
public :: alpha_coefficients_m_constructor

type, extends(alpha_coefficients_constructor) :: alpha_coefficients_m_constructor
  !< Henrick alpha coefficients object constructor.
  character(len=:), allocatable :: base_type !< Base alpha coefficient type.
endtype alpha_coefficients_m_constructor

interface  alpha_coefficients_m_constructor
  procedure alpha_coefficients_m_constructor_
endinterface

type, extends(alpha_coefficients) :: alpha_coefficients_m
  !< Henrick alpha coefficients object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Mapped weighted essentially
  !< non-oscillatory schemes: Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
  !< JCP, 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023.
  class(alpha_coefficients), allocatable :: alpha_base !< Base alpha coefficients to be re-mapped.
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute alpha coefficients.
    procedure, nopass     :: description !< Return alpha coefficients string-description.
    ! public methods
    procedure, pass(self) :: create  !< Create alpha coefficients.
    procedure, pass(self) :: destroy !< Destroy alpha coefficients.
endtype alpha_coefficients_m

contains
  ! function-constructor
  elemental function alpha_coefficients_m_constructor_(S, base_type) result(constructor)
  !< Return an instance of [alpha_coefficients_m_constructor].
  integer(I_P), intent(in)               :: S           !< Maximum stencils dimension.
  character(*), intent(in), optional     :: base_type   !< Base alpha coefficients type.
  type(alpha_coefficients_m_constructor) :: constructor !< Alpha coefficients constructor.

  constructor%S = S
  if (present(base_type)) constructor%base_type = base_type
  endfunction alpha_coefficients_m_constructor_

  ! deferred public methods
  pure subroutine compute(self, S, weight_opt, IS, eps, f1, f2)
  !< Compute alpha coefficients.
  class(alpha_coefficients_m), intent(inout) :: self                       !< Alpha coefficients.
  integer(I_P),                intent(in)    :: S                          !< Number of stencils used.
  real(R_P),                   intent(in)    :: weight_opt(1: 2, 0: S - 1) !< Optimal weight of the stencil.
  real(R_P),                   intent(in)    :: IS(1: 2, 0: S - 1)         !< Smoothness indicators of the stencils.
  real(R_P),                   intent(in)    :: eps                        !< Parameter for avoiding divided by zero.
  integer(I_P),                intent(in)    :: f1, f2                     !< Faces to be computed.
  integer(I_P)                               :: f, s1                      !< Counters.

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

  pure function description() result(string)
  !< Return alpha coefficients string-descripition.
  character(len=:), allocatable :: string           !< String-description.
  character(len=1), parameter   :: nl=new_line('a') !< New line character.

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
  endfunction description

  ! overridden methods
  pure subroutine create(self, constructor)
  !< Create alpha coefficients.
  class(alpha_coefficients_m),           intent(inout) :: self        !< Alpha coefficients.
  class(alpha_coefficients_constructor), intent(in)    :: constructor !< Alpha coefficients constructor.

  call self%destroy
  call self%alpha_coefficients%create(constructor=constructor)
  select type(constructor)
  type is(alpha_coefficients_m_constructor)
    if (allocated(constructor%base_type)) then
      select case(constructor%base_type)
      case('JS')
        if (allocated(self%alpha_base)) deallocate(self%alpha_base)
        allocate(alpha_coefficients_js :: self%alpha_base)
        call self%alpha_base%create(constructor=constructor)
      case('Z')
        if (allocated(self%alpha_base)) deallocate(self%alpha_base)
        allocate(alpha_coefficients_z :: self%alpha_base)
        call self%alpha_base%create(constructor=constructor)
      endselect
    endif
  endselect
  endsubroutine create

  pure subroutine destroy(self)
  !< Destroy alpha coefficients.
  class(alpha_coefficients_m), intent(inout) :: self !< Alpha coefficients.

  call self%alpha_coefficients%destroy
  if (allocated(self%alpha_base)) deallocate(self%alpha_base)
  endsubroutine destroy
endmodule wenoof_alpha_coefficients_m
