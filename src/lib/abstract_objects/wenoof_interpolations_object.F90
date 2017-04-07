!< Abstract interpolations object.
module wenoof_interpolations_object
!< Abstract interpolations object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
use wenoof_base_object

implicit none
private
public :: interpolations_object
public :: interpolations_object_constructor

type, extends(base_object_constructor), abstract :: interpolations_object_constructor
  !< Abstract interpolations object constructor.
  real(RPP), allocatable :: stencil(:)   !< Stencil used for interpolation, [1-S:S-1].
  real(RPP)              :: x_target     !< Coordinate of the interpolation point.
endtype interpolations_object_constructor

type, extends(base_object), abstract :: interpolations_object
  !< Abstract interpolations object.
  real(RPP), allocatable :: values_rank_1(:)   !< Stencil interpolations values [0:S-1].
  ! real(RPP), allocatable :: values_rank_2(:,:) !< Stencil interpolations values [1:2,0:S-1].
  contains
    ! public methods
    generic :: compute => compute_stencil_rank_1, compute_stencil_rank_2
    ! deferred public methods
    procedure(compute_stencil_rank_1_interface), pass(self), deferred :: compute_stencil_rank_1 !< Compute interpolation values.
    procedure(compute_stencil_rank_2_interface), pass(self), deferred :: compute_stencil_rank_2 !< Compute interpolation values.
endtype interpolations_object

abstract interface
  !< Abstract interfaces of [[interpolations_object]].
  pure subroutine compute_stencil_rank_1_interface(self, stencil)
  !< Compute interpolations.
  import :: interpolations_object, RPP
  class(interpolations_object), intent(inout) :: self               !< Interpolations.
  real(RPP),                    intent(in)    :: stencil(1-self%S:) !< Stencil used for the interpolation, [1-S:-1+S].
  endsubroutine compute_stencil_rank_1_interface

  pure subroutine compute_stencil_rank_2_interface(self, stencil, values)
  !< Compute interpolations.
  import :: interpolations_object, RPP
  class(interpolations_object), intent(in)  :: self                  !< Interpolations.
  real(RPP),                    intent(in)  :: stencil(1:,1-self%S:) !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  real(RPP),                    intent(out) :: values(1:, 0:)        !< Interpolations values.
  endsubroutine compute_stencil_rank_2_interface
endinterface
endmodule wenoof_interpolations_object
