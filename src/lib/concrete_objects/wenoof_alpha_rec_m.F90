!< Henrick alpha (non linear weights) object.
module wenoof_alpha_rec_m
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
public :: alpha_rec_m
public :: alpha_rec_m_constructor

type, extends(alpha_object_constructor) :: alpha_rec_m_constructor
  !< Henrick alpha (non linear weights) object constructor.
  character(len=:), allocatable :: base_type !< Base alpha coefficient type.
endtype alpha_rec_m_constructor

type, extends(alpha_object) :: alpha_rec_m
  !< Henrick alpha (non linear weights) object.
  !<
  !< @note The provided alpha implements the alpha coefficients defined in *Mapped weighted essentially non-oscillatory schemes:
  !< Achieving optimal order near critical points*, Andrew K. Henrick, Tariq D. Aslam, Joseph M. Powers,
  !< JCP, 2005, vol. 207, pp. 542-567, doi:10.1016/j.jcp.2005.01.023.
  real(RPP),           allocatable :: values(:,:)   !< Alpha coefficients [1:2,0:S-1].
  real(RPP),           allocatable :: values_sum(:) !< Sum of alpha coefficients [1:2].
  class(alpha_object), allocatable :: alpha_base    !< Base alpha to be re-mapped.
  contains
    ! public deferred methods
    procedure, pass(self) :: create                        !< Create alpha.
    procedure, pass(self) :: compute => compute_alpha_rec  !< Compute alpha.
    procedure, pass(self) :: description                   !< Return alpha string-description.
    procedure, pass(self) :: destroy                       !< Destroy alpha.
endtype alpha_rec_m

contains
  ! deferred public methods
  subroutine create(self, constructor)
  !< Create alpha.
  class(alpha_rec_m),             intent(inout) :: self        !< Alpha.
  class(base_object_constructor), intent(in)    :: constructor !< Alpha constructor.

  call self%destroy
  call self%create_(constructor=constructor)
  allocate(self%values(1:2, 0:self%S - 1))
  allocate(self%values_sum(1:2))
  self%values = 0._RPP
  self%values_sum = 0._RPP
  select type(constructor)
  type is(alpha_rec_m_constructor)
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

  pure subroutine compute_alpha_rec(self, beta, kappa)
  !< Compute alpha.
  class(alpha_rec_m),  intent(inout) :: self        !< Alpha.
  class(beta_object),  intent(in)    :: beta        !< Beta.
  class(kappa_object), intent(in)    :: kappa       !< Kappa.
  real(RPP)                          :: kappa_base  !< Kappa evaluated from the base alphas.
  integer(I_P)                       :: f, s1       !< Counters.

  self%values_sum = 0._RPP
  call self%alpha_base%compute(beta=beta, kappa=kappa)
  do s1=0, self%S - 1 ! stencil loops
    do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      kappa_base = self%alpha_base%values(f, s1) / self%alpha_base%values_sum(f)
      self%values(f, s1) =                                                                                    &
        (kappa_base * (kappa%values_rank_2(f, s1) + kappa%values_rank_2(f, s1) * kappa%values_rank_2(f, s1) - &
         3._RPP * kappa%values_rank_2(f, s1) * kappa_base + kappa_base *                                      &
         kappa_base)) /                                                                                       &
         (kappa%values_rank_2(f, s1) * kappa%values_rank_2(f, s1) + kappa_base *                              &
         (1._RPP - 2._RPP * kappa%values_rank_2(f, s1)))
      self%values_sum(f) = self%values_sum(f) + self%values(f, s1)
    enddo
  enddo
  endsubroutine compute_alpha_rec

  pure function description(self) result(string)
  !< Return alpha string-descripition.
  class(alpha_rec_m), intent(in) :: self             !< Alpha.
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
  class(alpha_rec_m), intent(inout) :: self !< Alpha.

  call self%destroy_
  if (allocated(self%values)) deallocate(self%values)
  if (allocated(self%values_sum)) deallocate(self%values_sum)
  if (allocated(self%alpha_base)) deallocate(self%alpha_base)
  endsubroutine destroy
endmodule wenoof_alpha_rec_m
