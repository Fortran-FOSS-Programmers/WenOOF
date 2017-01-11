!< Abstract optimal weights object.
module wenoof_optimal_weights
!< Abstract optimal weights object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: optimal_weights
public :: optimal_weights_constructor

type, extends(base_object_constructor) :: optimal_weights_constructor
  !< Abstract optimal weights object constructor.
  integer(I_P) :: S = 0 !< Stencils dimension.
endtype optimal_weights_constructor

type, extends(base_object) :: optimal_weights
  !< Optimal weights object.
  real(R_P), allocatable :: opt(:,:) !< Optimal weights [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute weights.
    procedure, nopass     :: description !< Return weights string-description.
    ! public methods
    procedure, pass(self) :: create  !< Createte weights.
    procedure, pass(self) :: destroy !< Destroy weights.
endtype optimal_weights

contains
  ! deferred public methods
  pure subroutine compute(self, S)
  !< Compute weights.
  class(optimal_weights), intent(inout) :: self !< Optimal weights.
  integer(I_P),           intent(in)    :: S    !< Number of stencils used.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'optimal_weights%compute to be implemented by your concrete optimal weights object'
#endif
  endsubroutine compute

  pure function description() result(string)
  !< Return weights string-description.
  character(len=:), allocatable  :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'optimal_weights%description to be implemented by your concrete optimal weights object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create weights.
  !<
  !< @note During creation the weights are also computed.
  class(optimal_weights),         intent(inout) :: self        !< Optimal weights.
  class(base_object_constructor), intent(in)    :: constructor !< Optimal weights constructor.

  call self%destroy
  select type(constructor)
  class is(optimal_weights_constructor)
    allocate(self%opt(1:2, 0:constructor%S - 1))
    call self%compute(S=constructor%S)
  class default
    ! @TODO add error handling
  endselect
  endsubroutine create

  elemental subroutine destroy(self)
  !< Destroy weights.
  class(optimal_weights), intent(inout) :: self !< Optimial weights.

  if (allocated(self%opt)) deallocate(self%opt)
  endsubroutine destroy
endmodule wenoof_optimal_weights
