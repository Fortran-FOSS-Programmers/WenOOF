!< Jiang-Shu alpha coefficients (non linear weights) object.
module wenoof_alpha_rec_js
!< Jiang-Shu alpha coefficients (non linear weights) object.
!<
!< @note The provided alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.

use penf, only : I_P, R_P
use wenoof_alpha_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_rec_js
public :: alpha_rec_js_constructor

type, extends(alpha_object_constructor) :: alpha_rec_js_constructor
  !< Jiang-Shu alpha coefficient object constructor.
endtype alpha_rec_js_constructor

type, extends(alpha_object) :: alpha_rec_js
  !< Jiang-Shu alpha coefficient object.
  !<
  !< @note The provided WENO alpha coefficient implements the alpha coefficients defined in *Efficient Implementation of Weighted
  !< ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  contains
    ! public deferred methods
    procedure, pass(self) :: compute     !< Compute alpha coefficients.
    procedure, nopass     :: description !< Return alpha coefficients string-description.
endtype alpha_rec_js

contains
  ! deferred public methods
  pure subroutine compute(self, beta, kappa)
  !< Compute alpha coefficients.
  class(alpha_rec_js), intent(inout) :: self  !< Alpha coefficient.
  class(beta_object),  intent(in)    :: beta  !< Beta coefficients.
  class(kappa_object), intent(in)    :: kappa !< Kappa coefficients.
  integer(I_P)                       :: f, s1 !< Counters.

  self%values_sum = 0._R_P
  do s1=0, self%S - 1 ! stencil loops
    do f=self%f1, self%f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%values(f, s1) =  kappa%values(f, s1)/(self%eps + beta%values(f, s1)) ** self%S
      self%values_sum(f) = self%values_sum(f) + self%values(f, s1)
    enddo
  enddo
  endsubroutine compute

  pure function description(self) result(string)
  !< Return alpha coefficients string-descripition.
  class(alpha_rec_js), intent(in) :: self             !< Alpha coefficient.
  character(len=:), allocatable   :: string           !< String-description.
  character(len=1), parameter     :: nl=new_line('a') !< New line character.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'alpha_rec_js%description to be implemented, do not use!'
#endif
  endfunction description
endmodule wenoof_alpha_rec_js
