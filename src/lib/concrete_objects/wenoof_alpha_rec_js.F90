!< Jiang-Shu alpha (non linear weights) object.
module wenoof_alpha_rec_js
!< Jiang-Shu alpha (non linear weights) object.
!<
!< @note The provided alpha implements the alpha coefficients defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.

#ifdef r16p
use penf, only: I_P, RPP=>R16P, str
#else
use penf, only: I_P, RPP=>R8P, str
#endif
use wenoof_alpha_object
use wenoof_base_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_rec_js
public :: alpha_rec_js_constructor

type, extends(alpha_object_constructor) :: alpha_rec_js_constructor
  !< Jiang-Shu alpha object constructor.
endtype alpha_rec_js_constructor

type, extends(alpha_object) :: alpha_rec_js
  !< Jiang-Shu alpha object.
  !<
  !< @note The provided WENO alpha implements the alpha coefficients defined in *Efficient Implementation of Weighted
  !< ENO Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130.
  real(RPP), allocatable :: values(:,:)   !< Alpha coefficients [1:2,0:S-1].
  real(RPP), allocatable :: values_sum(:) !< Sum of alpha coefficients [1:2].
  contains
    ! public deferred methods
    procedure, pass(self) :: create                        !< Create alpha.
    procedure, pass(self) :: compute => compute_alpha_rec  !< Compute alpha.
    procedure, pass(self) :: description                   !< Return alpha string-description.
    procedure, pass(self) :: destroy                       !< Destroy alpha.
endtype alpha_rec_js

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create alpha.
  class(alpha_rec_js),            intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  allocate(self%values_sum(1:2))
  self%values = 0._RPP
  self%values_sum = 0._RPP
  endsubroutine create

  pure subroutine compute_alpha_rec(self, beta, kappa)
  !< Compute alpha.
  class(alpha_rec_js), intent(inout) :: self  !< Alpha coefficient.
  class(beta_object),  intent(in)    :: beta  !< Beta coefficients.
  class(kappa_object), intent(in)    :: kappa !< Kappa coefficients.
  integer(I_P)                       :: f, s1 !< Counters.

  self%values_sum = 0._RPP
  do s1=0, self%S - 1 ! stencil loops
    do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%values(f, s1) = kappa%values_rank_2(f, s1)/(self%eps + beta%values(f, s1)) ** self%S
      self%values_sum(f) = self%values_sum(f) + self%values(f, s1)
    enddo
  enddo
  endsubroutine compute_alpha_rec

  pure function description(self) result(string)
  !< Return alpha string-descripition.
  class(alpha_rec_js), intent(in) :: self             !< Alpha coefficient.
  character(len=:), allocatable   :: string           !< String-description.
  character(len=1), parameter     :: nl=new_line('a') !< New line char.

  string = '    Jiang-Shu alpha coefficients for reconstructor:'//nl
  string = string//'      - S   = '//trim(str(self%S))//nl
  string = string//'      - eps = '//trim(str(self%eps))
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy alpha.
  class(alpha_rec_js), intent(inout) :: self !< Alpha.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  if (allocated(self%values_sum)) deallocate(self%values_sum)
  endsubroutine destroy
endmodule wenoof_alpha_rec_js
