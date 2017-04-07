!< Abstract alpha (non linear weights) object.
module wenoof_alpha_object
!< Abstract alpha (non linear weights) object.

#ifdef r16p
use penf, only: RPP=>R16P
#else
use penf, only: RPP=>R8P
#endif
use wenoof_base_object
use wenoof_beta_object
use wenoof_kappa_object

implicit none
private
public :: alpha_object
public :: alpha_object_constructor

type, extends(base_object_constructor), abstract :: alpha_object_constructor
  !< Abstract alpha (non linear weights) object constructor.
  contains
endtype alpha_object_constructor

type, extends(base_object), abstract :: alpha_object
  !< Abstract alpha (non linear weights) object.
  ! real(RPP), allocatable :: values_rank_1(:)      !< Alpha values [0:S-1].
  ! real(RPP)              :: values_sum_rank_1     !< Sum of alpha coefficients.
  ! real(RPP), allocatable :: values_rank_2(:,:)    !< Alpha values [1:2,0:S-1].
  ! real(RPP), allocatable :: values_sum_rank_2(:)  !< Sum of alpha coefficients [1:2].
  contains
    ! public deferred methods
    procedure(compute_interpolate_interface), pass(self), deferred :: compute_interpolate !< Compute alpha (interpolate).
    procedure(compute_reconstruct_interface), pass(self), deferred :: compute_reconstruct !< Compute alpha (reconstruct).
    generic :: compute => compute_interpolate, compute_reconstruct
endtype alpha_object

abstract interface
  !< Abstract interfaces of [[alpha_object]].
  pure subroutine compute_interpolate_interface(self, beta, kappa, values)
  !< Compute alpha.
  import :: alpha_object, beta_object, kappa_object
  class(alpha_object), intent(in)  :: self       !< Alpha.
  real(RPP),           intent(in)  :: beta(0:)   !< Beta [0:S-1].
  real(RPP),           intent(in)  :: kappa(0:)  !< Kappa [0:S-1].
  real(RPP),           intent(out) :: values(0:) !< Alpha values [0:S-1].
  endsubroutine compute_interpolate_interface

  pure subroutine compute_reconstruct_interface(self, beta, kappa, values)
  !< Compute alpha.
  import :: alpha_object, beta_object, kappa_object
  class(alpha_object), intent(in)  :: self          !< Alpha.
  real(RPP),           intent(in)  :: beta(1:,0:)   !< Beta [1:2,0:S-1].
  real(RPP),           intent(in)  :: kappa(1:,0:)  !< Kappa [1:2,0:S-1].
  real(RPP),           intent(out) :: values(1:,0:) !< Alpha values [1:2,0:S-1].
  endsubroutine compute_reconstruct_interface
endinterface
endmodule wenoof_alpha_object
