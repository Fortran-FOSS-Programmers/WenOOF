!< Abstract smoothness indicator object.
module wenoof_smoothness_indicators
!< Abstract smoothness indicator object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: smoothness_indicators
public :: smoothness_indicators_constructor

type, extends(base_object_constructor) :: smoothness_indicators_constructor
  !< Abstract smoothness indicators object constructor.
  integer(I_P) :: S = 0 !< Stencils dimension.
endtype smoothness_indicators_constructor

type, extends(base_object) :: smoothness_indicators
  !< Abstract smoothness indicator object.
  real(R_P), allocatable :: si(:,:) !< Smoothness indicators [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute smoothness indicators.
    procedure, nopass     :: description !< Return smoothness indicators string-description.
    ! public methods
    procedure, pass(self) :: create  !< Createte smoothness indicators.
    procedure, pass(self) :: destroy !< Destroy smoothness indicators.
endtype smoothness_indicators

contains
  ! deferred public methods
  pure subroutine compute(self, S, stencil, f1, f2, ff)
  !< Compute smoothness indicators.
  class(smoothness_indicators), intent(inout) :: self                !< Smoothness indicators.
  integer(I_P),                 intent(in)    :: S                   !< Number of stencils used.
  real(R_P),                    intent(in)    :: stencil(1:, 1 - S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P),                 intent(in)    :: f1, f2, ff          !< Faces to be computed.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'smoothness_indicators%compute to be implemented by your concrete smoothness indicators object'
#endif
  endsubroutine compute

  pure function description() result(string)
  !< Return smoothness indicators string-description.
  character(len=:), allocatable  :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'smoothness_indicators%description to be implemented by your concrete smoothness indicators object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create smoothness indicators.
  class(smoothness_indicators),   intent(inout) :: self        !< Smoothness indicators.
  class(base_object_constructor), intent(in)    :: constructor !< Smoothness indicators constructor.

  call self%destroy
  select type(constructor)
  class is(smoothness_indicators_constructor)
    allocate(self%si(1:2, 0:constructor%S - 1))
  class default
    ! @TODO add error handling
  endselect
  self%si = 0._R_P
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy smoothness indicators.
  class(smoothness_indicators), intent(inout) :: self !< Smoothness indicators.

  if (allocated(self%si)) deallocate(self%si)
  endsubroutine destroy
endmodule wenoof_smoothness_indicators
