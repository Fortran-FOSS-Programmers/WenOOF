!< WenOOF test: reconstruction of cosine function.
module cos_test_module
!< Auxiliary module defining the test class.

use flap, only : command_line_interface
#ifdef r16p
use penf, only: I_P, RPP=>R16P, FRPP=>FR16P, str, strz
#else
use penf, only: I_P, RPP=>R8P, FRPP=>FR8P, str, strz
#endif
use pyplot_module, only : pyplot
use wenoof, only : interpolator_object, wenoof_create
use wenoof_test_ui, only : test_ui

implicit none
private
public :: test

real(RPP), parameter     :: pi = 4._RPP * atan(1._RPP)  !< Pi greek.

type :: solution_data
  !< Class to handle solution data.
  real(RPP), allocatable :: x_cell(:)           !< Cell domain [1-S:points_number+S].
  real(RPP), allocatable :: fx_cell(:)          !< Cell refecence values [1-S:points_number+S].
  real(RPP), allocatable :: x_face(:,:)         !< Face domain [1:2,1:points_number].
  real(RPP), allocatable :: fx_face(:,:)        !< Face reference values [1:2,1:points_number].
  real(RPP), allocatable :: dfx_cell(:)         !< Cell refecence values of df/dx [1:points_number].
  real(RPP), allocatable :: interpolations(:,:) !< Interpolated values [1:2,1:points_number].
  real(RPP), allocatable :: reconstruction(:)   !< Reconstruction values [1:2,1:points_number].
  real(RPP), allocatable :: si(:,:,:)           !< Computed smoothness indicators [1:2,1:points_number,0:S-1].
  real(RPP), allocatable :: weights(:,:,:)      !< Computed weights [1:2,1:points_number,0:S-1].
  real(RPP)              :: Dx=0._RPP           !< Space step (spatial resolution).
  real(RPP)              :: error_L2=0._RPP     !< L2 norm of the numerical error.
endtype solution_data

type :: test
  !< Class to handle test(s).
  !<
  !< Test is driven by the Command Line Interface (CLI) options.
  !<
  !< Test has only 1 public method `execute`: it executes test(s) accordingly to cli options.
  private
  type(test_ui)                    :: ui            !< Command line interface handler.
  type(solution_data), allocatable :: solution(:,:) !< Solution [1:pn_number, 1:S_number].
  real(RPP), allocatable           :: accuracy(:,:) !< Accuracy (measured) [1:pn_number-1, 1:S_number].
  contains
    ! public methods
    procedure, pass(self) :: execute !< Execute selected test(s).
    ! private methods
    procedure, pass(self), private :: allocate_solution_data     !< Allocate solution data.
    procedure, pass(self), private :: analize_errors             !< Analize errors.
    procedure, pass(self), private :: compute_reference_solution !< Compute reference solution.
    procedure, pass(self), private :: deallocate_solution_data   !< Deallocate solution data.
    procedure, pass(self), private :: perform                    !< Perform test(s).
    procedure, pass(self), private :: save_results_and_plots     !< Save results and plots.
endtype test

contains
  ! public methods
  subroutine execute(self)
  !< Execute test(s).
  class(test), intent(inout) :: self !< Test.
  integer(I_P)               :: s    !< Counter.

  call self%ui%get
  if (trim(adjustl(self%ui%interpolator_type))/='all') then
    call self%perform
  else
    do while(self%ui%loop_interpolator(interpolator=self%ui%interpolator_type))
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
  allocate(self%solution(1:self%ui%pn_number, 1:self%ui%S_number))
  if (self%ui%pn_number>1) then
    allocate(self%accuracy(1:self%ui%pn_number, 1:self%ui%S_number))
    self%accuracy = 0._RPP
  endif
  do s=1, self%ui%S_number
    do pn=1, self%ui%pn_number
      allocate(self%solution(pn, s)%x_cell( 1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
      allocate(self%solution(pn, s)%fx_cell(1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
      allocate(self%solution(pn, s)%x_face(        1:2,  1:self%ui%points_number(pn)                              ))
      allocate(self%solution(pn, s)%fx_face(       1:2,  1:self%ui%points_number(pn)                              ))
      allocate(self%solution(pn, s)%dfx_cell(            1:self%ui%points_number(pn)                              ))
      allocate(self%solution(pn, s)%interpolations(1:2,  1:self%ui%points_number(pn)                              ))
      allocate(self%solution(pn, s)%reconstruction(      1:self%ui%points_number(pn)                              ))
      allocate(self%solution(pn, s)%si(            1:2,  1:self%ui%points_number(pn),             0:self%ui%S(s)-1))
      allocate(self%solution(pn, s)%weights(       1:2,  1:self%ui%points_number(pn),             0:self%ui%S(s)-1))
      self%solution(pn, s)%x_cell         = 0._RPP
      self%solution(pn, s)%fx_cell        = 0._RPP
      self%solution(pn, s)%x_face         = 0._RPP
      self%solution(pn, s)%fx_face        = 0._RPP
      self%solution(pn, s)%dfx_cell       = 0._RPP
      self%solution(pn, s)%interpolations = 0._RPP
      self%solution(pn, s)%reconstruction = 0._RPP
      self%solution(pn, s)%si             = 0._RPP
      self%solution(pn, s)%weights        = 0._RPP
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
  do s=1, self%ui%S_number
    do pn=1, self%ui%pn_number
      self%solution(pn, s)%Dx = 2 * pi / self%ui%points_number(pn)
      ! compute the values used for the interpolation/reconstruction of cos function: cell values
      do i=1 - self%ui%S(s), self%ui%points_number(pn) + self%ui%S(s)
        self%solution(pn, s)%x_cell(i) = i * self%solution(pn, s)%Dx - self%solution(pn, s)%Dx / 2._RPP
        self%solution(pn, s)%fx_cell(i) = sin(self%solution(pn, s)%x_cell(i))
      enddo
      ! values to which the interpolation/reconstruction should tend
      do i = 1, self%ui%points_number(pn)
        self%solution(pn, s)%x_face(1,i) = self%solution(pn, s)%x_cell(i) - self%solution(pn, s)%Dx / 2._RPP
        self%solution(pn, s)%x_face(2,i) = self%solution(pn, s)%x_cell(i) + self%solution(pn, s)%Dx / 2._RPP
        self%solution(pn, s)%fx_face(1,i) = sin(self%solution(pn, s)%x_face(1,i))
        self%solution(pn, s)%fx_face(2,i) = sin(self%solution(pn, s)%x_face(2,i))
        self%solution(pn, s)%dfx_cell(i) = cos(self%solution(pn, s)%x_cell(i))
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

  subroutine perform(self)
  !< Perform the test.
  class(test), intent(inout)              :: self         !< Test.
  real(RPP), allocatable                  :: error(:,:)   !< Error (norm L2) with respect the exact solution.
  real(RPP), allocatable                  :: order(:,:)   !< Observed order based on subsequent refined solutions.
  class(interpolator_object), allocatable :: interpolator !< WENO interpolator.
  real(RPP), allocatable                  :: stencil(:,:) !< Stencils used.
  integer(I_P)                            :: s            !< Counter.
  integer(I_P)                            :: pn           !< Counter.
  integer(I_P)                            :: i            !< Counter.

  call self%compute_reference_solution
  do s=1, self%ui%S_number
    call wenoof_create(interpolator_type='reconstructor-'//trim(adjustl(self%ui%interpolator_type)), &
                       S=self%ui%S(s),                                                               &
                       interpolator=interpolator,                                                    &
                       eps=self%ui%eps)
    if (self%ui%verbose) print '(A)', interpolator%description()
    allocate(stencil(1:2, 1-self%ui%S(s):-1+self%ui%S(s)))
    do pn=1, self%ui%pn_number
      do i=1, self%ui%points_number(pn)
        stencil(1,:) = self%solution(pn, s)%fx_cell(i+1-self%ui%S(s):i-1+self%ui%S(s))
        stencil(2,:) = self%solution(pn, s)%fx_cell(i+1-self%ui%S(s):i-1+self%ui%S(s))
        call interpolator%interpolate(stencil=stencil,                                        &
                                      interpolation=self%solution(pn, s)%interpolations(:,i), &
                                      si=self%solution(pn, s)%si(:, i, 0:self%ui%S(s)-1),     &
                                      weights=self%solution(pn, s)%weights(:, i, 0:self%ui%S(s)-1))
        self%solution(pn, s)%reconstruction(i) = &
          (self%solution(pn, s)%interpolations(2,i) - self%solution(pn, s)%interpolations(1,i))/self%solution(pn, s)%Dx
      enddo
    enddo
    deallocate(stencil)
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
  integer(I_P)                  :: f          !< Counter.

  output_dir = trim(adjustl(self%ui%output_dir))//'/'
  if (self%ui%results.or.self%ui%plots) call execute_command_line('mkdir -p '//output_dir)
  file_bname = output_dir//trim(adjustl(self%ui%output_bname))//'-'//trim(adjustl(self%ui%interpolator_type))

  if (self%ui%results) then
    do s=1, self%ui%S_number
      do pn=1, self%ui%pn_number
        open(newunit=file_unit, file=file_bname//'-S_'//trim(str(self%ui%S(s), .true.))//&
                                     '-Np_'//trim(str(self%ui%points_number(pn), .true.))//'.dat')
        buffer = 'VARIABLES = "x" "sin(x)" "cos(x)" "x_left" "x_right" "sin(x)_left" "sin(x)_right"'
        buffer = buffer//' "reconstruction_left" "reconstruction_right" "cos_reconstruction"'
        do ss=0, self%ui%S(s)-1
          buffer = buffer//' "si-'//trim(str(ss, .true.))//'_left"'//' "si-'//trim(str(ss, .true.))//'_right"'
        enddo
        do ss=0, self%ui%S(s)-1
          buffer = buffer//' "W-'//trim(str(ss, .true.))//'_left"'//' "W-'//trim(str(ss, .true.))//'_right"'
        enddo
        write(file_unit, "(A)") buffer
        write(file_unit, "(A)") 'ZONE T = "'//'S_'//trim(str(self%ui%S(s), .true.))//&
                                              '-Np_'//trim(str(self%ui%points_number(pn), .true.))//'"'
        associate(x_cell         => self%solution(pn, s)%x_cell,         &
                  fx_cell        => self%solution(pn, s)%fx_cell,        &
                  dfx_cell       => self%solution(pn, s)%dfx_cell,       &
                  x_face         => self%solution(pn, s)%x_face,         &
                  fx_face        => self%solution(pn, s)%fx_face,        &
                  interpolations => self%solution(pn, s)%interpolations, &
                  reconstruction => self%solution(pn, s)%reconstruction, &
                  si             => self%solution(pn, s)%si,             &
                  weights        => self%solution(pn, s)%weights,        &
                  Dx             => self%solution(pn, s)%Dx)
          do i = 1, self%ui%points_number(pn)
            write(file_unit, "("//trim(str(10+4*self%ui%S(s), .true.))//"("//FRPP//",1X))") &
               x_cell(i),                                                                   &
               fx_cell(i),                                                                  &
               dfx_cell(i),                                                                 &
              (x_face(f,i), f=1, 2),                                                        &
              (fx_face(f,i), f=1, 2),                                                       &
              (interpolations(f,i), f=1, 2),                                                &
               reconstruction(i),                                                           &
             ((si(f, i, ss), f=1, 2), ss=0, self%ui%S(s)-1),                                &
             ((weights(f, i, ss), f=1, 2), ss=0, self%ui%S(s)-1)
          enddo
        endassociate
        close(file_unit)
      enddo
    enddo

    if (self%ui%errors_analysis.and.self%ui%pn_number>1) then
      open(newunit=file_unit, file=file_bname//'-accuracy.dat')
      write(file_unit, "(A)") 'VARIABLES = "S" "Np" "error (L2)" "observed order" "formal order"'
      do s=1, self%ui%S_number
        do pn=1, self%ui%pn_number
          write(file_unit, "(2(I5,1X),"//FRPP//",1X,F5.2,1X,I3)") self%ui%S(s),                  &
                                                                  self%ui%points_number(pn),     &
                                                                  self%solution(pn, s)%error_L2, &
                                                                  self%accuracy(pn, s),          &
                                                                  2*self%ui%S(s)-1
        enddo
      enddo
      close(file_unit)
    endif
  endif

#ifndef r16p
  ! pyplot fortran does not support 128 bit reals
  if (self%ui%plots) then
    do s=1, self%ui%S_number
      do pn=1, self%ui%pn_number
        buffer = 'WENO reconstruction of $d \sin(x)/Dx=\cos(x)$; '//&
                 'S='//trim(str(self%ui%S(s), .true.))//'Np='//trim(str(self%ui%points_number(pn), .true.))
        call plt%initialize(grid=.true., xlabel='angle (rad)', title=buffer, legend=.true.)
        call plt%add_plot(x=self%solution(pn, s)%x_cell(1:self%ui%points_number(pn)), &
                          y=self%solution(pn, s)%dfx_cell(:),                         &
                          label='$\cos(x)$',                                          &
                          linestyle='k-',                                             &
                          linewidth=2,                                                &
                          ylim=[-1.1_RPP, 1.1_RPP])
        call plt%add_plot(x=self%solution(pn, s)%x_cell(1:self%ui%points_number(pn)),                            &
                          y=self%solution(pn, s)%reconstruction(:),                                              &
                          label='WENO reconstruction',                                                           &
                          linestyle='ro',                                                                        &
                          markersize=6,                                                                          &
                          ylim=[-1.1_RPP, 1.1_RPP])
        call plt%savefig(file_bname//&
                         '-S_'//trim(str(self%ui%S(s), .true.))//'-Np_'//trim(str(self%ui%points_number(pn), .true.))//'.png')
      enddo
    enddo
  endif
#endif
  endsubroutine save_results_and_plots

  subroutine analize_errors(self)
  !< Analize errors.
  class(test), intent(inout) :: self !< Test.
  integer(I_P)               :: s    !< Counter.
  integer(I_P)               :: pn   !< Counter.
  integer(I_P)               :: i    !< Counter.

  if (self%ui%errors_analysis) then
    do s=1, self%ui%S_number
      do pn=1, self%ui%pn_number
        associate(error_L2=>self%solution(pn, s)%error_L2, &
                  Dx=>self%solution(pn, s)%Dx, &
                  dfx_cell=>self%solution(pn, s)%dfx_cell, &
                  reconstruction=>self%solution(pn, s)%reconstruction)
          error_L2 = 0._RPP
          do i=1, self%ui%points_number(pn)
            error_L2 = error_L2 + (reconstruction(i) - dfx_cell(i))**2
          enddo
          error_L2 = sqrt(error_L2)
        endassociate
      enddo
    enddo
    if (self%ui%pn_number>1) then
      do s=1, self%ui%S_number
        do pn=2, self%ui%pn_number
          self%accuracy(pn, s) = log(self%solution(pn - 1, s)%error_L2 / self%solution(pn, s)%error_L2) / &
                                 log(self%solution(pn - 1, s)%Dx / self%solution(pn, s)%Dx)
        enddo
      enddo
    endif
  endif
  endsubroutine analize_errors
endmodule cos_test_module

program cos_reconstruction
!< WenOOF test: reconstruction of cosine function.

use cos_test_module

implicit none
type(test) :: cos_test

call cos_test%execute
endprogram cos_reconstruction
