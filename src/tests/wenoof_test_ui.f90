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

character(99), parameter :: interpolators(1:4) = ["JS  ", &
                                                  "M-JS", &
                                                  "M-Z ", &
                                                  "Z   "] !< List of available interpolators.

type :: test_ui
  !< Class to handle test(s) User Interface (UI).
  type(command_line_interface) :: cli                                  !< Command line interface handler.
  integer(I_P)                 :: error=0                              !< Error handler.
  character(99)                :: interpolator_type='JS'               !< Interpolator used.
  character(99)                :: output_bname='unset'                 !< Output files basename.
  character(99)                :: output_dir=''                        !< Output directory.
  integer(I_P)                 :: pn_number                            !< Number of different points-number tested.
  integer(I_P), allocatable    :: points_number(:)                     !< Points number used to discretize the domain.
  integer(I_P)                 :: S_number                             !< Number of different stencils tested.
  integer(I_P), allocatable    :: S(:)                                 !< Stencils used.
  real(RPP)                    :: eps                                  !< Small epsilon to avoid zero-division.
  real(RPP)                    :: x_target                             !< Interpolation target coordinate.
  logical                      :: interpolate                          !< Flag for activating interpolation.
  logical                      :: errors_analysis=.false.              !< Flag for activating errors analysis.
  logical                      :: plots=.false.                        !< Flag for activating plots saving.
  logical                      :: results=.false.                      !< Flag for activating results saving.
  logical                      :: verbose=.false.                      !< Flag for activating verbose output.
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
      call cli%init(progname    = 'WenOOF Test',                                        &
                    authors     = 'Fortran-FOSS-Programmers',                           &
                    license     = 'GNU GPLv3',                                          &
                    description = 'Test WenOOF library on function reconstruction',     &
                    examples    = ["$EXECUTABLE interpolate -i JS --results",    &
                                   "$EXECUTABLE interpolate -i JS-Z -r     ",    &
                                   "$EXECUTABLE reconstruct -i JS-M        ",    &
                                   "$EXECUTABLE reconstruct -i all -p -r   "])

      call cli%add_group(group='interpolate', description='perform WENO interpolation')
      call cli%add      (group='interpolate', switch='--x_target', switch_ab='-x',                                         &
                         help='WENO interpolation target point coordinate',                                                &
                         required=.true., def='0', act='store')
      call cli%add      (group='interpolate', switch='--interpolator', switch_ab='-i',                                     &
                         help='WENO interpolator type', required=.false.,                                                  &
                         def='JS', act='store', choices='all,JS,M-JS,M-Z,Z')
      call cli%add      (group='interpolate', switch='--points_number', switch_ab='-pn', nargs='+',                        &
                         help='Number of points used to discretize the domain',                                            &
                         required=.false., act='store', def='50 100')
      call cli%add      (group='interpolate', switch='--stencils', switch_ab='-s', nargs='+',                              &
                         help='Stencils dimensions (and number)', required=.false., act='store',                           &
                         def='2 3 4 5 6 7 8 9', choices='2, 3, 4, 5, 6, 7, 8, 9')
      call cli%add      (group='interpolate', switch='--eps', help='Small epsilon to avoid zero-division',                 &
                         required=.false., act='store', def='1.e-6')
      call cli%add      (group='interpolate', switch='--output_dir', help='Output directory', required=.false.,            &
                         act='store', def='./')
      call cli%add      (group='interpolate', switch='--results', switch_ab='-r', help='Save results', required=.false.,   &
                         act='store_true', def='.false.')
      call cli%add      (group='interpolate', switch='--plots', switch_ab='-p', help='Save plots', required=.false.,       &
                         act='store_true', def='.false.')
      call cli%add      (group='interpolate', switch='--output', help='Output files basename', required=.false.,           &
                         act='store', def='output')
      call cli%add      (group='interpolate', switch='--errors_analysis', help='Peform errors analysis', required=.false., &
                         act='store_true', def='.false.')
      call cli%add      (group='interpolate', switch='--verbose', help='Verbose output', required=.false.,                 &
                         act='store_true', def='.false.')

      call cli%add_group(group='reconstruct', description='perform WENO reconstruction')
      call cli%add      (group='reconstruct', switch='--interpolator', switch_ab='-i',                                     &
                         help='WENO interpolator type', required=.false.,                                                  &
                         def='JS', act='store', choices='all,JS,M-JS,M-Z,Z')
      call cli%add      (group='reconstruct', switch='--points_number', switch_ab='-pn', nargs='+',                        &
                         help='Number of points used to discretize the domain',                                            &
                         required=.false., act='store', def='50 100')
      call cli%add      (group='reconstruct', switch='--stencils', switch_ab='-s', nargs='+',                              &
                         help='Stencils dimensions (and number)', required=.false., act='store',                           &
                         def='2 3 4 5 6 7 8 9', choices='2, 3, 4, 5, 6, 7, 8, 9')
      call cli%add      (group='reconstruct', switch='--eps', help='Small epsilon to avoid zero-division',                 &
                         required=.false., act='store', def='1.e-6')
      call cli%add      (group='reconstruct', switch='--output_dir', help='Output directory', required=.false.,            &
                         act='store', def='./')
      call cli%add      (group='reconstruct', switch='--results', switch_ab='-r', help='Save results', required=.false.,   &
                         act='store_true', def='.false.')
      call cli%add      (group='reconstruct', switch='--plots', switch_ab='-p', help='Save plots', required=.false.,       &
                         act='store_true', def='.false.')
      call cli%add      (group='reconstruct', switch='--output', help='Output files basename', required=.false.,           &
                         act='store', def='output')
      call cli%add      (group='reconstruct', switch='--errors_analysis', help='Peform errors analysis', required=.false., &
                         act='store_true', def='.false.')
      call cli%add      (group='reconstruct', switch='--verbose', help='Verbose output', required=.false.,                 &
                         act='store_true', def='.false.')
    endassociate
    endsubroutine set_cli

    subroutine parse_cli()
    !< Parse Command Line Interface and check its validity.

    call self%cli%parse(error=self%error) ; if (self%error/=0) stop
    if (self%cli%run_command(group='interpolate')) then
      self%interpolate=.true.
      call self%cli%get(group='interpolate', switch='-x', val=self%x_target, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='-i', val=self%interpolator_type, error=self%error) ; if (self%error/=0) stop
      call self%cli%get_varying(group='interpolate', switch='-pn', val=self%points_number, error=self%error)
      if (self%error/=0) stop
      call self%cli%get_varying(group='interpolate', switch='-s', val=self%S, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='--eps', val=self%eps, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='--output_dir', val=self%output_dir, error=self%error)
      if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='-r', val=self%results, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='-p', val=self%plots, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='--output', val=self%output_bname, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='--errors_analysis', val=self%errors_analysis, error=self%error)
      if (self%error/=0) stop
      call self%cli%get(group='interpolate', switch='--verbose', val=self%verbose, error=self%error) ; if (self%error/=0) stop
    elseif (self%cli%run_command(group='reconstruct')) then
      self%interpolate=.false.
      call self%cli%get(group='reconstruct', switch='-i', val=self%interpolator_type, error=self%error) ; if (self%error/=0) stop
      call self%cli%get_varying(group='reconstruct', switch='-pn', val=self%points_number, error=self%error)
      if (self%error/=0) stop
      call self%cli%get_varying(group='reconstruct', switch='-s', val=self%S, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='reconstruct', switch='--eps', val=self%eps, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='reconstruct', switch='--output_dir', val=self%output_dir, error=self%error)
      if (self%error/=0) stop
      call self%cli%get(group='reconstruct', switch='-r', val=self%results, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='reconstruct', switch='-p', val=self%plots, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='reconstruct', switch='--output', val=self%output_bname, error=self%error) ; if (self%error/=0) stop
      call self%cli%get(group='reconstruct', switch='--errors_analysis', val=self%errors_analysis, error=self%error)
      if (self%error/=0) stop
      call self%cli%get(group='reconstruct', switch='--verbose', val=self%verbose, error=self%error) ; if (self%error/=0) stop
    else
#ifndef DEBUG
      ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
      call self%cli%print_usage
      error stop 'error: action not present; choose the correct action between "interpolate" and "reconstruct"'
#endif
    endif

    self%pn_number = size(self%points_number, dim=1)
    self%S_number = size(self%S, dim=1)
    endsubroutine parse_cli
  endsubroutine get

  function loop_interpolator(self, interpolator) result(again)
  !< Loop over available interpolators.
  class(test_ui), intent(in)  :: self          !< Test UI.
  character(99),  intent(out) :: interpolator  !< Interpolator name.
  character(99)               :: approximation !< Approximation type.
  logical                     :: again         !< Flag continuing the loop.
  integer(I_P), save          :: i = 0         !< Counter.

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
