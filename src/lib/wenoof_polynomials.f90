!< Abstract polynomials object.
module wenoof_polynomials
!< Abstract polynomials object.

use penf, only : I_P, R_P
use wenoof_base_object

implicit none
private
public :: polynomials
public :: polynomials_constructor

type, extends(base_object_constructor) :: polynomials_constructor
  !< Abstract polynomials object constructor.
  integer(I_P) :: S = 0 !< Stencils dimension.
endtype polynomials_constructor

type, extends(base_object) :: polynomials
  !< Abstract polynomials object.
  real(R_P), allocatable :: poly(:,:)   !< Polynomial reconstructions [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self) :: compute     !< Compute polynomials.
    procedure, nopass     :: description !< Return polynomials string-description.
    ! public methods
    procedure, pass(self) :: create  !< Createte polynomials.
    procedure, pass(self) :: destroy !< Destroy polynomials.
endtype polynomials

contains
  ! deferred public methods
  pure subroutine compute(self, S, stencil, f1, f2, ff)
  !< Compute polynomials.
  class(polynomials), intent(inout) :: self                !< Polynomials.
  integer(I_P),       intent(in)    :: S                   !< Number of stencils used.
  real(R_P),          intent(in)    :: stencil(1:, 1 - S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P),       intent(in)    :: f1, f2, ff          !< Faces to be computed.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'polynomials%compute to be implemented by your concrete polynomials object'
#endif
  endsubroutine compute

  pure function description() result(string)
  !< Return polynomials string-description.
  character(len=:), allocatable  :: string !< String-description.

#ifndef DEBUG
  ! error stop in pure procedure is a F2015 feature not yet supported in debug mode
  error stop 'polynomials%description to be implemented by your concrete polynomials object'
#endif
  endfunction description

  ! public methods
  pure subroutine create(self, constructor)
  !< Create polynomials.
  class(polynomials),             intent(inout) :: self        !< Polynomials.
  class(polynomials_constructor), intent(in)    :: constructor !< Polynomials constructor.

  call self%destroy
  allocate(self%poly(1:2, 0:constructor%S - 1))
  self%poly = 0._R_P
  endsubroutine create

  pure subroutine destroy(self)
  !< Destroy polynomials.
  class(polynomials), intent(inout) :: self !< Polynomials.

  if (allocated(self%poly)) deallocate(self%poly)
  endsubroutine destroy
endmodule wenoof_polynomials
