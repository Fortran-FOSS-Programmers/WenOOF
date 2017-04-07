!< Henrick alpha (non linear weights) object.
module wenoof_alpha_int_m
!< Henrick alpha (non linear weights) object.
!<
!< @note The provided alpha implements the alpha coefficients defined in *Mapped weighted essentially non-oscillatory schemes:
!< Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers, JCP,
!< 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023

#ifdef r16p
use penf, only: I_P, RPP=>R16P, str
#else
use penf, only: I_P, RPP=>R8P, str
#endif
use wenoof_alpha_object
use wenoof_alpha_rec_js
use wenoof_alpha_rec_z
use wenoof_base_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_int_m
public :: alpha_int_m_constructor

type, extends(alpha_object_constructor) :: alpha_int_m_constructor
  !< Henrick alpha (non linear weights) object constructor.
  character(len=:), allocatable :: base_type !< Base alpha coefficient type.
endtype alpha_int_m_constructor

type, extends(alpha_object) :: alpha_int_m
  !< Henrick alpha (non linear weights) object.
  !<
  !< @note The provided alpha implements the alpha coefficients defined in *Mapped weighted essentially non-oscillatory schemes:
  !< Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
  !< JCP, 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023.
  real(RPP),           allocatable :: values(:)   !< Alpha coefficients [0:S-1].
  real(RPP)                        :: values_sum  !< Sum of alpha coefficients.
  class(alpha_object), allocatable :: alpha_base  !< Base alpha to be re-mapped.
  contains
    ! public deferred methods
    procedure, pass(self) :: create                        !< Create alpha.
    procedure, pass(self) :: compute => compute_alpha_int  !< Compute alpha.
    procedure, pass(self) :: description                   !< Return alpha string-description.
    procedure, pass(self) :: destroy                       !< Destroy alpha.
endtype alpha_int_m

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create alpha.
  class(alpha_int_m),             intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values_rank_1(0:self%S - 1))
  associate(val => self%values_rank_1, val_sum => self%values_sum_rank_1)
    val = 0._RPP
    val_sum = 0._RPP
  endassociate
  select type(constructor)
  type is(alpha_int_m_constructor)
    if (allocated(constructor%base_type)) then
      select case(constructor%base_type)
      case('JS')
        if (allocated(self%alpha_base)) deallocate(self%alpha_base)
        allocate(alpha_rec_js :: self%alpha_base)
        call self%alpha_base%create(constructor=constructor)
      case('Z')
        if (allocated(self%alpha_base)) deallocate(self%alpha_base)
        allocate(alpha_rec_z :: self%alpha_base)
        call self%alpha_base%create(constructor=constructor)
      endselect
    endif
  class default
    ! @TODO add error handling
  endselect
  endsubroutine create

  pure subroutine compute_alpha_int(self, beta, kappa)
  !< Compute alpha.
  class(alpha_int_m),  intent(inout) :: self        !< Alpha.
  class(beta_object),  intent(in)    :: beta        !< Beta.
  class(kappa_object), intent(in)    :: kappa       !< Kappa.
  real(RPP)                          :: kappa_base  !< Kappa evaluated from the base alphas.
  integer(I_P)                       :: s1          !< Counter.

  associate(val => self%values_rank_1, val_sum => self%values_sum_rank_1)
    val_sum = 0._RPP
    call self%alpha_base%compute(beta=beta, kappa=kappa)
    do s1=0, self%S - 1 ! stencil loops
      kappa_base = self%alpha_base%values_rank_1(s1) / self%alpha_base%values_sum_rank_1
      val(s1) =                                                                                      &
        (kappa_base * (kappa%values_rank_1(s1) + kappa%values_rank_1(s1) * kappa%values_rank_1(s1) - &
         3._RPP * kappa%values_rank_1(s1) * kappa_base + kappa_base * kappa_base)) /                 &
         (kappa%values_rank_1(s1) * kappa%values_rank_1(s1) + kappa_base *                           &
         (1._RPP - 2._RPP * kappa%values_rank_1(s1)))
      val_sum = val_sum + val(s1)
    enddo
  endassociate
  endsubroutine compute_alpha_int

  pure function description(self) result(string)
  !< Return alpha string-descripition.
  class(alpha_int_m), intent(in) :: self             !< Alpha.
  character(len=:), allocatable  :: string           !< String-description.
  character(len=1), parameter    :: nl=new_line('a') !< New line char.

  string = '    Henrick alpha coefficients for reconstructor:'//nl
  string = string//'      - S   = '//trim(str(self%S))//nl
  string = string//'      - eps = '//trim(str(self%eps))//nl
  associate(alpha_base=>self%alpha_base)
    select type(alpha_base)
    type is(alpha_rec_js)
      string = string//'      - base-mapped-alpha type = Jiang-Shu'
    type is(alpha_rec_z)
      string = string//'      - base-mapped-alpha type = Bogeg'
    endselect
  endassociate
  endfunction description

  elemental subroutine destroy(self)
  !< Destroy alpha.
  class(alpha_int_m), intent(inout) :: self !< Alpha.

  call self%destroy_
  if (allocated(self%values_rank_1)) deallocate(self%values_rank_1)
  if (allocated(self%alpha_base)) deallocate(self%alpha_base)
  endsubroutine destroy
endmodule wenoof_alpha_int_m
