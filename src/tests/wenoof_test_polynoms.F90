!< WenOOF test: interpolation or reconstruction of polynomial functions.

module wenoof_test_polynoms_module
!< Auxiliary module defining the test class.

use flap, only : command_line_interface
use penf, only : FR_P, I_P, R_P, str, strz
use pyplot_module, only : pyplot
use wenoof, only : interpolator_object, wenoof_create
use wenoof_test_ui, only : test_ui

implicit none
private
public :: test

type :: solution_data
  !< Class to handle solution data.
  real(R_P), allocatable :: x_cell(:)           !< Cell domain [1-S:points_number+S].
  real(R_P), allocatable :: fx_cell(:)          !< Cell refecence values [1-S:points_number+S].
  real(R_P), allocatable :: x_face(:,:)         !< Face domain [1:2,1:points_number].
  real(R_P), allocatable :: fx_face(:,:)        !< Face reference values [1:2,1:points_number].
  real(R_P), allocatable :: x_int(:)            !< Interpolation domain [1-S:points_number+S].
  real(R_P), allocatable :: fx_int(:)           !< Interpolation refecence values [1-S:points_number+S].
  real(R_P), allocatable :: dfx_cell(:)         !< Cell refecence values of df/dx [1:points_number].
  real(R_P), allocatable :: interpolations(:,:) !< Interpolated values [1:2,1:points_number].
  real(R_P), allocatable :: reconstruction(:)   !< Reconstruction values [1:2,1:points_number].
  real(R_P), allocatable :: si_r(:,:,:)         !< Computed smoothness indicators [1:2,1:points_number,0:S-1].
  real(R_P), allocatable :: weights_r(:,:,:)    !< Computed weights [1:2,1:points_number,0:S-1].
  real(R_P), allocatable :: interpolation(:)    !< Interpolated values [1:points_number].
  real(R_P), allocatable :: si_i(:,:)           !< Computed smoothness indicators [1:points_number,0:S-1].
  real(R_P), allocatable :: weights_i(:,:)      !< Computed weights [1:points_number,0:S-1].
  real(R_P)              :: error_L2            !< L2 norm of the numerical error.
  real(R_P)              :: x_target            !< Abscissa of the interpolation [-0.5:0.5].
  real(R_P)              :: Dx=0._R_P           !< Space step (spatial resolution).
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
  real(R_P), allocatable           :: accuracy(:,:) !< Accuracy (measured) [1:pn_number-1, 1:S_number].
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

  call self%ui%get
  if (self%ui%interpolate.and.self%ui%reconstruct) then
    call subexecute
    self%ui%interpolate = .false.
    call subexecute
  else
    call subexecute
  endif
  contains
     subroutine subexecute
     !< Subexecute test(s).

     if (trim(adjustl(self%ui%interpolator_type))/='all') then
       call self%perform
     else
       do while(self%ui%loop_interpolator(interpolator=self%ui%interpolator_type))
         call self%perform
       enddo
     endif
     endsubroutine subexecute
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
    self%accuracy = 0._R_P
  endif
  if (self%ui%interpolate) then
    self%solution%x_target = self%ui%x_target
    do s=1, self%ui%S_number
      do pn=1, self%ui%pn_number
        allocate(self%solution(pn, s)%x_cell( 1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
        allocate(self%solution(pn, s)%fx_cell(1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
        allocate(self%solution(pn, s)%x_int(  1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
        allocate(self%solution(pn, s)%fx_int( 1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
        allocate(self%solution(pn, s)%interpolation(       1:self%ui%points_number(pn)                              ))
        allocate(self%solution(pn, s)%si_i(                1:self%ui%points_number(pn),             0:self%ui%S(s)-1))
        allocate(self%solution(pn, s)%weights_i(           1:self%ui%points_number(pn),             0:self%ui%S(s)-1))
        self%solution(pn, s)%x_cell         = 0._R_P
        self%solution(pn, s)%fx_cell        = 0._R_P
        self%solution(pn, s)%x_int          = 0._R_P
        self%solution(pn, s)%fx_int         = 0._R_P
        self%solution(pn, s)%interpolation  = 0._R_P
        self%solution(pn, s)%si_i           = 0._R_P
        self%solution(pn, s)%weights_i      = 0._R_P
        self%solution(pn, s)%error_L2       = 0._R_P
      enddo
    enddo
  else
    do s=1, self%ui%S_number
      do pn=1, self%ui%pn_number
        allocate(self%solution(pn, s)%x_cell( 1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
        allocate(self%solution(pn, s)%fx_cell(1-self%ui%S(s):self%ui%points_number(pn)+self%ui%S(s)                 ))
        allocate(self%solution(pn, s)%x_face(        1:2,  1:self%ui%points_number(pn)                              ))
        allocate(self%solution(pn, s)%fx_face(       1:2,  1:self%ui%points_number(pn)                              ))
        allocate(self%solution(pn, s)%dfx_cell(            1:self%ui%points_number(pn)                              ))
        allocate(self%solution(pn, s)%interpolations(1:2,  1:self%ui%points_number(pn)                              ))
        allocate(self%solution(pn, s)%reconstruction(      1:self%ui%points_number(pn)                              ))
        allocate(self%solution(pn, s)%si_r(          1:2,  1:self%ui%points_number(pn),             0:self%ui%S(s)-1))
        allocate(self%solution(pn, s)%weights_r(     1:2,  1:self%ui%points_number(pn),             0:self%ui%S(s)-1))
        self%solution(pn, s)%x_cell         = 0._R_P
        self%solution(pn, s)%fx_cell        = 0._R_P
        self%solution(pn, s)%x_face         = 0._R_P
        self%solution(pn, s)%fx_face        = 0._R_P
        self%solution(pn, s)%dfx_cell       = 0._R_P
        self%solution(pn, s)%interpolations = 0._R_P
        self%solution(pn, s)%reconstruction = 0._R_P
        self%solution(pn, s)%si_r           = 0._R_P
        self%solution(pn, s)%weights_r      = 0._R_P
        self%solution(pn, s)%error_L2       = 0._R_P
      enddo
    enddo
  endif
  endsubroutine allocate_solution_data

  subroutine compute_reference_solution(self)
  !< Allocate solution data.
  class(test),  intent(inout) :: self !< Test.
  integer(I_P)                :: s    !< Counter.
  integer(I_P)                :: pn   !< Counter.
  integer(I_P)                :: i    !< Counter.

  call self%allocate_solution_data
  if (self%ui%interpolate) then
    do s=1, self%ui%S_number
      do pn=1, self%ui%pn_number
        self%solution(pn, s)%Dx = 1._R_P / self%ui%points_number(pn)
        ! compute the values used for the interpolation of polynomials function: cell values
        do i=1 - self%ui%S(s), self%ui%points_number(pn) + self%ui%S(s)
          self%solution(pn, s)%x_cell(i) = i * self%solution(pn, s)%Dx - self%solution(pn, s)%Dx / 2._R_P
          self%solution(pn, s)%fx_cell(i) = interface_function(x=self%solution(pn, s)%x_cell(i), o=2*self%ui%S(s)+2)
        enddo
        ! values to which the interpolation should tend
        do i = 1, self%ui%points_number(pn)
          self%solution(pn, s)%x_int(i) = self%solution(pn, s)%x_cell(i) + self%solution(pn, s)%x_target * self%solution(pn, s)%Dx
          self%solution(pn, s)%fx_int(i) = interface_function(self%solution(pn, s)%x_int(i), o=2*self%ui%S(s)+2)
        enddo
      enddo
    enddo
  else
    do s=1, self%ui%S_number
      do pn=1, self%ui%pn_number
        self%solution(pn, s)%Dx = 1._R_P / self%ui%points_number(pn)
        ! compute the values used for the interpolation/reconstruction of polynomials function: cell values
        do i=1 - self%ui%S(s), self%ui%points_number(pn) + self%ui%S(s)
          self%solution(pn, s)%x_cell(i) = i * self%solution(pn, s)%Dx - self%solution(pn, s)%Dx / 2._R_P
          self%solution(pn, s)%fx_cell(i) = interface_function(x=self%solution(pn, s)%x_cell(i), o=2*self%ui%S(s)+2)
        enddo
        ! values to which the interpolation/reconstruction should tend
        do i = 1, self%ui%points_number(pn)
          self%solution(pn, s)%x_face(1,i) = self%solution(pn, s)%x_cell(i) - self%solution(pn, s)%Dx / 2._R_P
          self%solution(pn, s)%x_face(2,i) = self%solution(pn, s)%x_cell(i) + self%solution(pn, s)%Dx / 2._R_P
          self%solution(pn, s)%fx_face(1,i) = interface_function(self%solution(pn, s)%x_face(1,i), o=2*self%ui%S(s)+2)
          self%solution(pn, s)%fx_face(2,i) = interface_function(self%solution(pn, s)%x_face(2,i), o=2*self%ui%S(s)+2)
          self%solution(pn, s)%dfx_cell(i) = dinterface_function_dx(self%solution(pn, s)%x_cell(i), o=2*self%ui%S(s)+2)
        enddo
      enddo
    enddo
  endif
  endsubroutine compute_reference_solution

  subroutine deallocate_solution_data(self)
  !< Deallocate solution data.
  class(test), intent(inout) :: self !< Test.

  if (allocated(self%solution)) deallocate(self%solution)
  if (allocated(self%accuracy)) deallocate(self%accuracy)
  endsubroutine deallocate_solution_data

  subroutine perform(self)
  !< Perform the test.
  class(test), intent(inout)              :: self           !< Test.
  class(interpolator_object), allocatable :: interpolator   !< WENO interpolator.
  real(R_P), allocatable                  :: stencil_i(:)   !< Stencils used for interpolation.
  real(R_P), allocatable                  :: stencil_r(:,:) !< Stencils used for reconstruction.
  integer(I_P)                            :: s              !< Counter.
  integer(I_P)                            :: pn             !< Counter.
  integer(I_P)                            :: i              !< Counter.

  call self%compute_reference_solution
  if (self%ui%interpolate) then
    do s=1, self%ui%S_number
    call wenoof_create(interpolator_type='interpolator-'//trim(adjustl(self%ui%interpolator_type)), &
                       S=self%ui%S(s),                                                              &
                       x_target=0.3_R_P,                                                            &
                       interpolator=interpolator,                                                   &
                       eps=self%ui%eps)
    if (self%ui%verbose) print '(A)', interpolator%description()
    allocate(stencil_i(1-self%ui%S(s):-1+self%ui%S(s)))
    do pn=1, self%ui%pn_number
      do i=1, self%ui%points_number(pn)
        stencil_i(:) = self%solution(pn, s)%fx_cell(i+1-self%ui%S(s):i-1+self%ui%S(s))
        call interpolator%interpolate(stencil=stencil_i,                                   &
                                      interpolation=self%solution(pn, s)%interpolation(i), &
                                      si=self%solution(pn, s)%si_i(i, 0:self%ui%S(s)-1),   &
                                      weights=self%solution(pn, s)%weights_i(i, 0:self%ui%S(s)-1))
      enddo
    enddo
    deallocate(stencil_i)
  enddo
  else
    do s=1, self%ui%S_number
      call wenoof_create(interpolator_type='reconstructor-'//trim(adjustl(self%ui%interpolator_type)), &
                         S=self%ui%S(s),                                                               &
                         interpolator=interpolator,                                                    &
                         eps=self%ui%eps)
      if (self%ui%verbose) print '(A)', interpolator%description()
      allocate(stencil_r(1:2, 1-self%ui%S(s):-1+self%ui%S(s)))
      do pn=1, self%ui%pn_number
        do i=1, self%ui%points_number(pn)
          stencil_r(1,:) = self%solution(pn, s)%fx_cell(i+1-self%ui%S(s):i-1+self%ui%S(s))
          stencil_r(2,:) = self%solution(pn, s)%fx_cell(i+1-self%ui%S(s):i-1+self%ui%S(s))
          call interpolator%interpolate(stencil=stencil_r,                                        &
                                        interpolation=self%solution(pn, s)%interpolations(:,i), &
                                        si=self%solution(pn, s)%si_r(:, i, 0:self%ui%S(s)-1),   &
                                        weights=self%solution(pn, s)%weights_r(:, i, 0:self%ui%S(s)-1))
          self%solution(pn, s)%reconstruction(i) = &
            (self%solution(pn, s)%interpolations(2,i) - self%solution(pn, s)%interpolations(1,i))/self%solution(pn, s)%Dx
        enddo
      enddo
      deallocate(stencil_r)
    enddo
  endif
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
  integer(I_P)                  :: istat      !< Status IO error.

  output_dir = trim(adjustl(self%ui%output_dir))//'/'
  if (self%ui%results.or.self%ui%plots) call execute_command_line('mkdir -p '//output_dir)

  if (self%ui%interpolate) then
    file_bname = output_dir//trim(adjustl(self%ui%output_bname))//'-interpolator-'//trim(adjustl(self%ui%interpolator_type))
      if (self%ui%results) then
        do s=1, self%ui%S_number
          do pn=1, self%ui%pn_number
            open(newunit=file_unit, file=file_bname//'-S_'//trim(str(self%ui%S(s), .true.))//&
                                         '-Np_'//trim(str(self%ui%points_number(pn), .true.))//'.dat')
            buffer = 'VARIABLES = "x" "f(x)" "x_int" "f(x)_int" "interpolation"'
            do ss=0, self%ui%S(s)-1
              buffer = buffer//' "si-'//trim(str(ss, .true.))
            enddo
            do ss=0, self%ui%S(s)-1
              buffer = buffer//' "W-'//trim(str(ss, .true.))
            enddo
            write(file_unit, "(A)") buffer
            write(file_unit, "(A)") 'ZONE T = "'//'S_'//trim(str(self%ui%S(s), .true.))//&
                                                  '-Np_'//trim(str(self%ui%points_number(pn), .true.))//'"'
            associate(x_cell        => self%solution(pn, s)%x_cell,        &
                      fx_cell       => self%solution(pn, s)%fx_cell,       &
                      x_int         => self%solution(pn, s)%x_int,         &
                      fx_int        => self%solution(pn, s)%fx_int,        &
                      interpolation => self%solution(pn, s)%interpolation, &
                      si            => self%solution(pn, s)%si_i,          &
                      weights       => self%solution(pn, s)%weights_i,     &
                      Dx            => self%solution(pn, s)%Dx)
              do i = 1, self%ui%points_number(pn)
                write(file_unit, "("//trim(str(5+2*self%ui%S(s), .true.))//"("//FR_P//",1X))") &
                   x_cell(i),                                                                  &
                   fx_cell(i),                                                                 &
                   x_int(i),                                                                   &
                   fx_int(i),                                                                  &
                   interpolation(i),                                                           &
                  (si(i, ss), ss=0, self%ui%S(s)-1),                                           &
                  (weights(i, ss), ss=0, self%ui%S(s)-1)
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
              write(file_unit, "(2(I5,1X),"//FR_P//",1X,F5.2,1X,I3)") self%ui%S(s),                  &
                                                                      self%ui%points_number(pn),     &
                                                                      self%solution(pn, s)%error_L2, &
                                                                      self%accuracy(pn, s),          &
                                                                      2*self%ui%S(s)-1
            enddo
          enddo
          close(file_unit)
        endif
      endif
  else
    file_bname = output_dir//trim(adjustl(self%ui%output_bname))//'-reconstructor-'//trim(adjustl(self%ui%interpolator_type))
    if (self%ui%results) then
      do s=1, self%ui%S_number
        do pn=1, self%ui%pn_number
          open(newunit=file_unit, file=file_bname//'-S_'//trim(str(self%ui%S(s), .true.))//&
                                       '-Np_'//trim(str(self%ui%points_number(pn), .true.))//'.dat')
          buffer = 'VARIABLES = "x" "f(x)" "df_dx(x)" "x_left" "x_right" "f(x)_left" "f(x)_right"'
          buffer = buffer//' "reconstruction_left" "reconstruction_right" "df_dx_reconstruction"'
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
                    si             => self%solution(pn, s)%si_r,           &
                    weights        => self%solution(pn, s)%weights_r,      &
                    Dx             => self%solution(pn, s)%Dx)
            do i = 1, self%ui%points_number(pn)
              write(file_unit, "("//trim(str(10+4*self%ui%S(s), .true.))//"("//FR_P//",1X))") &
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
            write(file_unit, "(2(I5,1X),"//FR_P//",1X,F5.2,1X,I3)") self%ui%S(s),                  &
                                                                    self%ui%points_number(pn),     &
                                                                    self%solution(pn, s)%error_L2, &
                                                                    self%accuracy(pn, s),          &
                                                                    2*self%ui%S(s)-1
          enddo
        enddo
        close(file_unit)
      endif
    endif
  endif

#ifndef r16p
  ! pyplot fortran does not support 128 bit reals
  if (self%ui%plots) then
    if (self%ui%interpolate) then
      do s=1, self%ui%S_number
        do pn=1, self%ui%pn_number
          buffer = 'WENO interpolation of polynomial function; '//&
                   'S='//trim(str(self%ui%S(s), .true.))//'Np='//trim(str(self%ui%points_number(pn), .true.))
          call plt%initialize(grid=.true., xlabel='x (m)', title=buffer, legend=.true.)
          call plt%add_plot(x=self%solution(pn, s)%x_cell(1:self%ui%points_number(pn)), &
                            y=self%solution(pn, s)%fx_cell(:),                          &
                            label='polynom',                                            &
                            linestyle='k-',                                             &
                            linewidth=2,                                                &
                            ylim=[-1.1_R_P, 1.1_R_P], istat=istat)
          call plt%add_plot(x=self%solution(pn, s)%x_int(1:self%ui%points_number(pn)), &
                            y=self%solution(pn, s)%interpolation(:),                   &
                            label='WENO interpolation',                                &
                            linestyle='ro',                                            &
                            markersize=6,                                              &
                            ylim=[-1.1_R_P, 1.1_R_P], istat=istat)
          call plt%savefig(file_bname//                                                                                          &
                           '-S_'//trim(str(self%ui%S(s), .true.))//'-Np_'//trim(str(self%ui%points_number(pn), .true.))//'.png', &
                           istat=istat)
        enddo
      enddo
    else
      do s=1, self%ui%S_number
        do pn=1, self%ui%pn_number
          buffer = 'WENO reconstruction of $d \p(x)/Dx; '//&
                   'S='//trim(str(self%ui%S(s), .true.))//'Np='//trim(str(self%ui%points_number(pn), .true.))
          call plt%initialize(grid=.true., xlabel='x (m)', title=buffer, legend=.true.)
          call plt%add_plot(x=self%solution(pn, s)%x_cell(1:self%ui%points_number(pn)),  &
                            y=self%solution(pn, s)%dfx_cell(:),                          &
                            label='dp',                                                  &
                            linestyle='k-',                                              &
                            linewidth=2,                                                 &
                            ylim=[-1.1_R_P, 1.1_R_P], istat=istat)
          call plt%add_plot(x=self%solution(pn, s)%x_cell(1:self%ui%points_number(pn)),  &
                            y=self%solution(pn, s)%reconstruction(:),                    &
                            label='WENO reconstruction',                                 &
                            linestyle='ro',                                              &
                            markersize=6,                                                &
                            ylim=[-1.1_R_P, 1.1_R_P], istat=istat)
          call plt%savefig(file_bname//                                                                                          &
                           '-S_'//trim(str(self%ui%S(s), .true.))//'-Np_'//trim(str(self%ui%points_number(pn), .true.))//'.png', &
                           istat=istat)
        enddo
      enddo
    endif
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
    if (self%ui%interpolate) then
      do s=1, self%ui%S_number
        do pn=1, self%ui%pn_number
          associate(error_L2=>self%solution(pn, s)%error_L2, &
                    Dx=>self%solution(pn, s)%Dx, &
                    fx_int=>self%solution(pn, s)%fx_int, &
                    interpolation=>self%solution(pn, s)%interpolation)
            error_L2 = 0._R_P
            do i=1, self%ui%points_number(pn)
              error_L2 = error_L2 + (interpolation(i) - fx_int(i))**2
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
    else
      do s=1, self%ui%S_number
        do pn=1, self%ui%pn_number
          associate(error_L2=>self%solution(pn, s)%error_L2, &
                    Dx=>self%solution(pn, s)%Dx, &
                    dfx_cell=>self%solution(pn, s)%dfx_cell, &
                    reconstruction=>self%solution(pn, s)%reconstruction)
            error_L2 = 0._R_P
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
  endif
  endsubroutine analize_errors

  ! non TBP
  pure function interface_function(x, o) result(y)
  !< Interface function.
  real(R_P),    intent(in) :: x !< X value.
  integer(I_P), intent(in) :: o !< Polynomial order.
  real(R_P)                :: y !< Y value.
  integer(I_P)             :: i !< Counter.

  y = 0._R_P
  do i=1, o
    y = y + i * (x ** i)
  enddo
  endfunction interface_function

  pure function dinterface_function_dx(x, o) result(y)
  !< Derivative of interface function.
  real(R_P),    intent(in) :: x !< X value.
  integer(I_P), intent(in) :: o !< Polynomial order.
  real(R_P)                :: y !< Y value.
  integer(I_P)             :: i !< Counter.

  y = 0._R_P
  do i=1, o
    y = y + i * i * (x ** (i - 1))
  enddo
  endfunction dinterface_function_dx
endmodule wenoof_test_polynoms_module

program wenoof_test_polynoms
!< WenOOF test: interpolation of polynomial functions.

use wenoof_test_polynoms_module

implicit none
type(test) :: polynoms_test

call polynoms_test%execute
endprogram wenoof_test_polynoms
