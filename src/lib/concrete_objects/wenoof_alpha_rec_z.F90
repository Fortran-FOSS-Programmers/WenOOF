!< Borges alpha (non linear weights) object.
module wenoof_alpha_rec_z
!< Borges alpha (non linear weights) object.
!<
!< @note The provided WENO alpha implements the alpha coefficients defined in *An improved weighted essentially non-oscillatory
!< scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP, 2008,
!< vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.

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
public :: alpha_rec_z
public :: alpha_rec_z_constructor

type, extends(alpha_object_constructor) :: alpha_rec_z_constructor
  !< Borges alpha (non linear weights) object constructor.
endtype alpha_rec_z_constructor

type, extends(alpha_object) :: alpha_rec_z
  !< Borges alpha (non linear weights) object.
  !<
  !< @note The provided alpha implements the alpha coefficients defined in *An improved weighted essentially non-oscillatory
  !< scheme for hyperbolic conservation laws*, Rafael Borges, Monique Carmona, Bruno Costa and Wai Sun Don, JCP,
  !< 2008, vol. 227, pp. 3191-3211, doi: 10.1016/j.jcp.2007.11.038.
  real(RPP), allocatable :: values(:,:)   !< Alpha coefficients [1:2,0:S-1].
  real(RPP), allocatable :: values_sum(:) !< Sum of alpha coefficients [1:2].
  contains
    ! public deferred methods
    procedure, pass(self) :: create                        !< Create alpha.
    procedure, pass(self) :: compute => compute_alpha_rec  !< Compute alpha.
    procedure, pass(self) :: description                   !< Return alpha string-description.
    procedure, pass(self) :: destroy                       !< Destroy alpha.
endtype alpha_rec_z
contains
  ! public deferred methods
  subroutine create(self, constructor)
  !< Create alpha.
  class(alpha_rec_z),             intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values_rank_2(1:2, 0:self%S - 1))
  allocate(self%values_sum_rank_2(1:2))
  associate(val => self%values_rank_2, val_sum => self%values_sum_rank_2)
    val = 0._RPP
    val_sum = 0._RPP
  endassociate
  endsubroutine create

  pure subroutine compute_alpha_rec(self, beta, kappa)
  !< Compute alpha.
  class(alpha_rec_z),  intent(inout) :: self  !< Alpha.
  class(beta_object),  intent(in)    :: beta  !< Beta.
  class(kappa_object), intent(in)    :: kappa !< Kappa.
  integer(I_P)                       :: f, s1 !< Counters.

  associate(val => self%values_rank_2, val_sum => self%values_sum_rank_2)
    val_sum = 0._RPP
    do s1=0, self%S - 1 ! stencil loops
      do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        val(f, s1) = kappa%values_rank_2(f, s1) *                         &
                     ((1._RPP + (tau(S=self%S, beta=beta%values_rank_2) / &
                     (self%eps + beta%values_rank_2(f, s1)))) ** (weno_exp(self%S)))
        val_sum(f) = val_sum(f) + val(f, s1)
      enddo
    enddo
  endassociate
  endsubroutine compute_alpha_rec

  pure function description(self) result(string)
  !< Return alpha string-descripition.
  class(alpha_rec_z), intent(in) :: self             !< Alpha coefficients.
  character(len=:), allocatable  :: string           !< String-description.
  character(len=1), parameter    :: nl=new_line('a') !< New line char.

  string = '    Borges alpha coefficients for reconstructor:'//nl
  string = string//'      - S   = '//trim(str(self%S))//nl
  string = string//'      - eps = '//trim(str(self%eps))

  endfunction description

  elemental subroutine destroy(self)
  !< Destroy alpha.
  class(alpha_rec_z), intent(inout) :: self !< Alpha.

  call self%destroy_
  if (allocated(self%values_rank_2)) deallocate(self%values_rank_2)
  if (allocated(self%values_sum_rank_2)) deallocate(self%values_sum_rank_2)
  endsubroutine destroy

  ! private non TBP
  pure function tau(S, beta) result(w_tau)
  !< Compute the tau coefficient used in the WENO-Z alpha coefficients.
  integer(I_P), intent(in) :: S           !< Number of stencils used.
  real(RPP),    intent(in) :: beta(0:S-1) !< Smoothness indicators.
  real(RPP)                :: w_tau       !< Tau coefficient.

  w_tau = abs(beta(0) -                           &
              (1_I_P - weno_odd(S)) * beta(1) -   &
              (1_I_P - weno_odd(S)) * beta(S-2) + &
              (1_I_P - 2_I_P * weno_odd(S)) * beta(S-1))
  endfunction tau

  pure function weno_exp(S) result(w_exp)
  !< Compute the exponent used in the alpha function.
  integer(I_P), intent(in) :: S     !< Number of stencils used.
  integer(I_P)             :: w_exp !< Exponent used in the alpha function.

  w_exp = int(S, I_P)
  endfunction weno_exp

  pure function weno_odd(S) result(w_odd)
  !< Compute the distinguisher between odd and even number of stencils.
  integer(I_P), intent(in) :: S     !< Number of stencils used.
  integer(I_P)             :: w_odd !< Distinguishing between odd and even number of stencils.

  w_odd = int(mod(S, 2_I_P), I_P)
  endfunction weno_odd
endmodule wenoof_alpha_rec_z
