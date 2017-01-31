!< WenOOF test UI: definition of common User Interface (UI) for WenOOF tests.
module wenoof_test_ui
!< WenOOF test UI: definition of common User Interface (UI) for WenOOF tests.

use flap, only : command_line_interface
#ifdef r16p
use penf, only: I_P, RPP=>R16P, FRPP=>FR16P
#else
use penf, only: I_P, RPP=>R8P, FRPP=>FR8P
#endif

implicit none
private
public :: test_ui

character(99), parameter :: interpolators(1:4) = ["reconstructor-JS  ", &
                                                  "reconstructor-M-JS", &
                                                  "reconstructor-M-Z ", &
                                                  "reconstructor-Z   "] !< List of available interpolators.

type :: test_ui
  !< Class to handle test(s) User Interface (UI).
  type(command_line_interface) :: cli                     !< Command line interface handler.
  integer(I_P)                 :: error=0                 !< Error handler.
  character(99)                :: interpolator_type='JS'  !< Interpolator used.
  character(99)                :: output_bname='unset'    !< Output files basename.
  character(99)                :: output_dir=''           !< Output directory.
  integer(I_P)                 :: pn_number               !< Number of different points-number tested.
  integer(I_P), allocatable    :: points_number(:)        !< Points number used to discretize the domain.
  integer(I_P)                 :: S_number                !< Number of different stencils tested.
  integer(I_P), allocatable    :: S(:)                    !< Stencils used.
  real(RPP)                    :: eps                     !< Smal episol to avoid zero-division.
  logical                      :: errors_analysis=.false. !< Flag for activating errors analysis.
  logical                      :: plots=.false.           !< Flag for activating plots saving.
  logical                      :: results=.false.         !< Flag for activating results saving.
  logical                      :: verbose=.false.         !< Flag for activating verbose output.
  contains
    ! public methods
    procedure, pass(self) :: get               !< Get user options.
    procedure, pass(self) :: loop_interpolator !< Loop over available interpolators.
endtype test_ui

contains
  ! public methods
  subroutine get(self)
  !< Get user options.
  class(test_ui), intent(inout) :: self !< Test UI.

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
      call cli%add(switch='--interpolator', switch_ab='-i', help='WENO interpolator type', required=.false., &
                   def='reconstructor-JS', act='store',                                                      &
                   choices='all,reconstructor-JS,reconstructor-M-JS,reconstructor-M-Z,reconstructor-Z')
      call cli%add(switch='--points_number', switch_ab='-pn', nargs='+', help='Number of points used to discretize the domain', &
                   required=.false., act='store', def='50 100')
      call cli%add(switch='--stencils', switch_ab='-s', nargs='+', help='Stencils dimensions (and number)', &
                   required=.false., act='store', def='2 3 4 5 6 7 8 9', choices='2, 3, 4, 5, 6, 7, 8, 9')
      call cli%add(switch='--eps', help='Small epsilon to avoid zero-division', required=.false., act='store', def='1.e-6')
      call cli%add(switch='--output_dir', help='Output directory', required=.false., act='store', def='./')
      call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--plots', switch_ab='-p', help='Save plots', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--output', help='Output files basename', required=.false., act='store', def='sin_reconstruction')
      call cli%add(switch='--errors_analysis', help='Peform errors analysis', required=.false., act='store_true', def='.false.')
      call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.')
    endassociate
    endsubroutine set_cli

    subroutine parse_cli()
    !< Parse Command Line Interface and check its validity.

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
    call self%cli%get(switch='--verbose', val=self%verbose, error=self%error) ; if (self%error/=0) stop

    self%pn_number = size(self%points_number, dim=1)
    self%S_number = size(self%S, dim=1)
    endsubroutine parse_cli
  endsubroutine get

  function loop_interpolator(self, interpolator) result(again)
  !< Loop over available interpolators.
  class(test_ui), intent(in)  :: self         !< Test UI.
  character(99),  intent(out) :: interpolator !< Interpolator name.
  logical                     :: again        !< Flag continuing the loop.
  integer(I_P), save          :: i = 0        !< Counter.

  again = .false.
  if (i==0) then
    i = 1
    interpolator = trim(adjustl(interpolators(i)))
    again = .true.
  elseif (i<size(interpolators, dim=1)) then
    i = i + 1
    interpolator = trim(adjustl(interpolators(i)))
    again = .true.
  else
    i = 0
    interpolator = ''
    again = .false.
  endif
  endfunction loop_interpolator
endmodule wenoof_test_ui
