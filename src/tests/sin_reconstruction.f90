!< WenOOF test: reconstruction of sin function.

module test_module
!< Auxiliary module defining the test class.

use flap, only : command_line_interface
use penf, only : I_P, R_P, FR_P, str, strz
use pyplot_module, only :  pyplot
use wenoof, only : interpolator, wenoof_create

implicit none
private
public :: test

character(99), parameter :: interpolators(1:4) = ["all ", &
                                                  "JS  ", &
                                                  "JS-Z", &
                                                  "JS-M"] !< List of available interpolators.
real(R_P), parameter     :: pi = 4._R_P * atan(1._R_P)  !< Pi greek.

type :: solution_data
  !< Class to handle solution data.
  real(R_P), allocatable :: x_cell(:)        !< Cell domain [1-S:points_number+S].
  real(R_P), allocatable :: fx_cell(:)       !< Cell refecence values [1-S:points_number+S].
  real(R_P), allocatable :: x_face(:)        !< Face domain [1:points_number].
  real(R_P), allocatable :: fx_face(:)       !< Face reference values [1:points_number].
  real(R_P), allocatable :: interpolation(:) !< Interpolated values [1:points_number].
  real(R_P), allocatable :: si(:,:)          !< Computed smoothness indicators [1:points_number,0:S-1].
  real(R_P)              :: error_L2         !< L2 norm of the numerical error.
endtype solution_data

type :: test
  !< Class to handle test(s).
  !<
  !< Test is driven by the Command Line Interface (CLI) options.
  !<
  !< Test has only 1 public method `execute`: it executes test(s) accordingly to cli options.
  private
  type(command_line_interface)     :: cli                     !< Command line interface handler.
  integer(I_P)                     :: error=0                 !< Error handler.
  character(99)                    :: interpolator_type='JS'  !< Interpolator used.
  character(99)                    :: output_bname='unset'    !< Output files basename.
  character(99)                    :: output_dir=''           !< Output directory.
  integer(I_P)                     :: pn_number               !< Number of different points-number tested.
  integer(I_P), allocatable        :: points_number(:)        !< Points number used to discretize the domain.
  integer(I_P)                     :: S_number                !< Number of different stencils tested.
  integer(I_P), allocatable        :: S(:)                    !< Stencils used.
  real(R_P)                        :: eps                     !< Smal episol to avoid zero-division.
  type(solution_data), allocatable :: solution(:,:)           !< Solution [1:pn_number, 1:S_number].
  real(R_P), allocatable           :: accuracy(:,:)           !< Accuracy (measured) [1:pn_number-1, 1:S_number].
  logical                          :: errors_analysis=.false. !< Flag for activating errors analysis.
  logical                          :: plots=.false.           !< Flag for activating plots saving.
  logical                          :: results=.false.         !< Flag for activating results saving.
  contains
    ! public methods
    procedure, pass(self) :: execute !< Execute selected test(s).
    ! private methods
    procedure, pass(self), private :: allocate_solution_data     !< Allocate solution data.
    procedure, pass(self), private :: analize_errors             !< Analize errors.
    procedure, pass(self), private :: compute_reference_solution !< Compute reference solution.
    procedure, pass(self), private :: deallocate_solution_data   !< Deallocate solution data.
    procedure, pass(self), private :: initialize                 !< Initialize test(s).
    procedure, pass(self), private :: perform                    !< Perform test(s).
    procedure, pass(self), private :: save_results_and_plots     !< Save results and plots.
endtype test

contains
  ! public methods
  subroutine execute(self)
  !< Execute test(s).
  class(test), intent(inout) :: self !< Test.
  integer(I_P)               :: s    !< Counter.

  call self%initialize
  if (trim(adjustl(self%interpolator_type))/='all') then
    call self%perform
  else
    do s=2, size(interpolators, dim=1)
      self%interpolator_type = trim(adjustl(interpolators(s)))
      call self%perform
    enddo
  endif
  endsubroutine execute

  ! private methods
  subroutine allocate_solution_data(self)
  !< Allocate solution data.
  class(test),  intent(inout) :: self !< Test.
  integer(I_P)                :: s    !< Counter.
  integer(I_P)                :: pn   !< Counter.

  call self%deallocate_solution_data
  self%pn_number = size(self%points_number, dim=1)
  self%S_number = size(self%S, dim=1)
  allocate(self%solution(1:self%pn_number, 1:self%S_number))
  if (self%pn_number>1) then
    allocate(self%accuracy(1:self%pn_number, 1:self%S_number))
    self%accuracy = 0._R_P
  endif
  do s=1, self%S_number
    do pn=1, self%pn_number
      allocate(self%solution(pn, s)%x_cell( 1-self%S(s):self%points_number(pn)+self%S(s)              ))
      allocate(self%solution(pn, s)%fx_cell(1-self%S(s):self%points_number(pn)+self%S(s)              ))
      allocate(self%solution(pn, s)%x_face(           1:self%points_number(pn)                        ))
      allocate(self%solution(pn, s)%fx_face(          1:self%points_number(pn)                        ))
      allocate(self%solution(pn, s)%interpolation(    1:self%points_number(pn)                        ))
      allocate(self%solution(pn, s)%si(               1:self%points_number(pn),          0:self%S(s)-1))
      self%solution(pn, s)%x_cell        = 0._R_P
      self%solution(pn, s)%fx_cell       = 0._R_P
      self%solution(pn, s)%x_face        = 0._R_P
      self%solution(pn, s)%fx_face       = 0._R_P
      self%solution(pn, s)%interpolation = 0._R_P
      self%solution(pn, s)%si            = 0._R_P
    enddo
  enddo
  endsubroutine allocate_solution_data

  subroutine compute_reference_solution(self)
  !< Allocate solution data.
  class(test),  intent(inout) :: self !< Test.
  integer(I_P)                :: s    !< Counter.
  integer(I_P)                :: pn   !< Counter.
  integer(I_P)                :: i    !< Counter.

  call self%allocate_solution_data
  do s=1, self%S_number
    do pn=1, self%pn_number
      ! compute the values used for the reconstruction of sin function: cell values
      do i=1 - self%S(s), self%points_number(pn) + self%S(s)
        self%solution(pn, s)%x_cell(i) = i * 2 * pi / self%points_number(pn)
        self%solution(pn, s)%fx_cell(i) = sin(self%solution(pn, s)%x_cell(i))
      enddo
      ! face values to which the reconstruction should tend
      do i = 1, self%points_number(pn)
        self%solution(pn, s)%x_face(i) = self%solution(pn, s)%x_cell(i) + pi / self%points_number(pn)
        self%solution(pn, s)%fx_face(i) = sin(self%solution(pn, s)%x_face(i))
      enddo
    enddo
  enddo
  endsubroutine compute_reference_solution

  subroutine deallocate_solution_data(self)
  !< Deallocate solution data.
  class(test), intent(inout) :: self !< Test.

  if (allocated(self%solution)) deallocate(self%solution)
  if (allocated(self%accuracy)) deallocate(self%accuracy)
  endsubroutine deallocate_solution_data

  subroutine initialize(self)
  !< Initialize test: set Command Line Interface, parse it and check its validity.
  class(test), intent(inout) :: self !< Test.

  call set_cli
  call parse_cli
  contains
    subroutine set_cli()
    !< Set Command Line Interface.

    associate(cli => self%cli)
      call cli%init(progname    = 'sin reconstruction',                                 &
                    authors     = 'Fortran-FOSS-Programmers',                           &
                    license     = 'GNU GPLv3',                                          &
                    description = 'Test WenOOF library on sin function reconstruction', &
                    examples    = ["sin_reconstruction --interpolator JS --results",    &
                                   "sin_reconstruction --interpolator JS-Z -r     ",    &
                                   "sin_reconstruction --interpolator JS-M        ",    &
                                   "sin_reconstruction --interpolator all -p -r   "])
      call cli%add(switch='--interpolator', switch_ab='-i', help='WENO interpolator type', required=.false., def='JS', act='store')
      call cli%add(switch='--points_number', switch_ab='-pn', nargs='+', help='Number of points used to discretize the domain', &
                   required=.false., act='store', def='50')
      call cli%add(switch='--stencils', switch_ab='-s', nargs='+', help='Stencils dimensions (and number)', &
                   required=.false., act='store', def='2', choices='2, 3, 4, 5, 6, 7, 8, 9')
      call cli%add(switch='--eps', help='Small epsilon to avoid zero-division', required=.false., act='store', def='1.e-6')
      call cli%add(switch='--output_dir', help='Output directory', required=.false., act='store', def='./')
      call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--plots', switch_ab='-p', help='Save plots', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--output', help='Output files basename', required=.false., act='store', def='sin_reconstruction')
      call cli%add(switch='--errors_analysis', help='Peform errors analysis', required=.false., act='store_true', def='.false.')
    endassociate
    endsubroutine set_cli

    subroutine parse_cli()
    !< Parse Command Line Interface and check its validity.
    character(len=:), allocatable :: valid_solvers_list !< Pretty printed list of available solvers.

    call self%cli%parse(error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-i', val=self%interpolator_type, error=self%error) ; if (self%error/=0) stop
    call self%cli%get_varying(switch='-pn', val=self%points_number, error=self%error) ; if (self%error/=0) stop
    call self%cli%get_varying(switch='-s', val=self%S, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--eps', val=self%eps, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--output_dir', val=self%output_dir, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-r', val=self%results, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='-p', val=self%plots, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--output', val=self%output_bname, error=self%error) ; if (self%error/=0) stop
    call self%cli%get(switch='--errors_analysis', val=self%errors_analysis, error=self%error) ; if (self%error/=0) stop

    if (.not.is_interpolator_valid()) then
      print "(A)", 'error: the interpolator type "'//trim(adjustl(self%interpolator_type))//'" is unknown!'
      print "(A)", list_interpolators()
      stop
    endif
    endsubroutine parse_cli

    function is_interpolator_valid()
    !< Verify if the selected interpolator is valid.
    logical      :: is_interpolator_valid !< Return true is the selected interpolator is available.
    integer(I_P) :: s                   !< Counter.

    is_interpolator_valid = .false.
    do s=1, size(interpolators, dim=1)
      is_interpolator_valid = (trim(adjustl(self%interpolator_type))==trim(adjustl(interpolators(s))))
      if (is_interpolator_valid) exit
    enddo
    endfunction is_interpolator_valid

    function list_interpolators() result(list)
    !< List available solvers.
    character(len=:), allocatable :: list !< Pretty printed list of available interpolators.
    integer(I_P)                  :: s    !< Counter.

    list = 'Valid interpolator names are:' // new_line('a')
    do s=1, ubound(interpolators, dim=1)
      list = list // '  + ' // trim(adjustl(interpolators(s))) // new_line('a')
    enddo
    endfunction list_interpolators
  endsubroutine initialize

  subroutine perform(self)
  !< Perform the test.
  class(test), intent(inout)       :: self               !< Test.
  real(R_P), allocatable           :: error(:,:)         !< Error (norm L2) with respect the exact solution.
  real(R_P), allocatable           :: order(:,:)         !< Observed order based on subsequent refined solutions.
  class(interpolator), allocatable :: weno_interpolator  !< WENO interpolator.
  integer(I_P)                     :: s                  !< Counter.
  integer(I_P)                     :: pn                 !< Counter.
  integer(I_P)                     :: i                  !< Counter.

  call self%compute_reference_solution
  do s=1, self%S_number
    call wenoof_create(interpolator_type=trim(adjustl(self%interpolator_type)), &
                       S=self%S(s),                                             &
                       eps=self%eps,                                            &
                       wenoof_interpolator=weno_interpolator)
    do pn=1, self%pn_number
      do i=1, self%points_number(pn)
        call weno_interpolator%interpolate(S=self%S(s),                                                                      &
                                           stencil=reshape(source=self%solution(pn, s)%fx_cell(i+1-self%S(s):i-1+self%S(s)), &
                                                           shape=[1,2*self%S(s)-1]),                                         &
                                           location='right',                                                                 &
                                           interpolation=self%solution(pn, s)%interpolation(i:i),                            &
                                           si=self%solution(pn, s)%si(i:i, 0:self%S(s)-1))
      enddo
    enddo
  enddo
  call self%analize_errors
  call self%save_results_and_plots
  endsubroutine perform

  subroutine save_results_and_plots(self)
  !< Save results and plots.
  class(test), intent(inout)    :: self       !< Test.
  type(pyplot)                  :: plt        !< Plot handler.
  character(len=:), allocatable :: buffer     !< Buffer string.
  character(len=:), allocatable :: output_dir !< Output directory.
  character(len=:), allocatable :: file_bname !< File base name.
  integer(I_P)                  :: file_unit  !< File unit.
  integer(I_P)                  :: s          !< Counter.
  integer(I_P)                  :: pn         !< Counter.
  integer(I_P)                  :: i          !< Counter.
  integer(I_P)                  :: ss         !< Counter.

  output_dir = trim(adjustl(self%output_dir))//'/'
  if (self%results.or.self%plots) call execute_command_line('mkdir -p '//output_dir)
  file_bname = output_dir//trim(adjustl(self%output_bname))//'-'//trim(adjustl(self%interpolator_type))

  if (self%results) then
    do s=1, self%S_number
      open(newunit=file_unit, file=file_bname//'-S_'//trim(str(self%S(s), .true.))//'.dat')
      buffer = 'VARIABLES = "x" "sin(x)" "weno_interpolation"'
      do ss=0, self%S(s)-1
        buffer = buffer//' "si-'//trim(str(ss, .true.))//'"'
      enddo
      write(file_unit, "(A)") buffer
      do pn=1, self%pn_number
        write(file_unit, "(A)") 'ZONE T = "'//'S_'//trim(str(self%S(s), .true.))//&
                                              '-Np_'//trim(str(self%points_number(pn), .true.))//'"'
        do i = 1, self%points_number(pn)
          write(file_unit, "("//trim(str(3+self%S(s), .true.))//"("//FR_P//",1X))") &
             self%solution(pn, s)%x_face(i),        &
             self%solution(pn, s)%fx_face(i),       &
             self%solution(pn, s)%interpolation(i), &
            (self%solution(pn, s)%si(i, ss), ss=0, self%S(s)-1)
        enddo
      enddo
      close(file_unit)
    enddo

    if (self%errors_analysis.and.self%pn_number>1) then
      open(newunit=file_unit, file=file_bname//'-accuracy.dat')
      write(file_unit, "(A)") 'VARIABLES = "S" "Np" "error (L2)" "observed order"'
      do s=1, self%S_number
        do pn=1, self%pn_number
          write(file_unit, "(2(I5,1X),2("//FR_P//",1X))") self%S(s),                     &
                                                          self%points_number(pn),        &
                                                          self%solution(pn, s)%error_L2, &
                                                          self%accuracy(pn, s)
        enddo
      enddo
      close(file_unit)
    endif
  endif

  if (self%plots) then
    do s=1, self%S_number
      do pn=1, self%pn_number
        buffer = 'WENO interpolation of $\sin(x)$; '//&
                 'S='//trim(str(self%S(s), .true.))//'Np='//trim(str(self%points_number(pn), .true.))
        call plt%initialize(grid=.true., xlabel='angle (rad)', title=buffer, legend=.true.)
        call plt%add_plot(x=self%solution(pn, s)%x_face(:),  &
                          y=self%solution(pn, s)%fx_face(:), &
                          label='$\sin(x)$',                 &
                          linestyle='k-',                    &
                          linewidth=2,                       &
                          ylim=[-1.1_R_P, 1.1_R_P])
        call plt%add_plot(x=self%solution(pn, s)%x_face(:),        &
                          y=self%solution(pn, s)%interpolation(:), &
                          label='WENO interpolation',              &
                          linestyle='ro',                          &
                          markersize=6,                            &
                          ylim=[-1.1_R_P, 1.1_R_P])
        call plt%savefig(file_bname//&
                         '-S_'//trim(str(self%S(s), .true.))//'-Np_'//trim(str(self%points_number(pn), .true.))//'.png')
      enddo
    enddo
  endif
  endsubroutine save_results_and_plots

  subroutine analize_errors(self)
  !< Analize errors.
  class(test), intent(inout) :: self !< Test.
  integer(I_P)               :: s    !< Counter.
  integer(I_P)               :: pn   !< Counter.
  integer(I_P)               :: i    !< Counter.

  if (self%errors_analysis) then
    do s=1, self%S_number
      do pn=1, self%pn_number
        associate(error_L2=>self%solution(pn, s)%error_L2, &
                  fx_face=>self%solution(pn, s)%fx_face,   &
                  interpolation=>self%solution(pn, s)%interpolation)
          error_L2 = 0._R_P
          do i=1, self%points_number(pn)
            error_L2 = error_L2 + (interpolation(i) - fx_face(i))**2
          enddo
          error_L2 = sqrt(error_L2)
        endassociate
      enddo
    enddo
    if (self%pn_number>1) then
      do s=1, self%S_number
        do pn=2, self%pn_number
          self%accuracy(pn, s) = log(self%solution(pn - 1, s)%error_L2 / self%solution(pn, s)%error_L2) / &
                                 log((1._R_P*self%points_number(pn)) / self%points_number(pn - 1))
        enddo
      enddo
    endif
  endif
  endsubroutine analize_errors
endmodule test_module

program sin_reconstruction
!< WenOOF test: reconstruction of sin function.

use test_module

implicit none
type(test) :: sin_test

call sin_test%execute
endprogram sin_reconstruction
