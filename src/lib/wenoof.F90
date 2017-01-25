!< WenOOF, WENO interpolation Object Oriented Fortran library
module wenoof
!< WenOOF, WENO interpolation Object Oriented Fortran library

#ifdef r16p
use penf, only: I_P, RPP=>R16P
#else
use penf, only: I_P, RPP=>R8P
#endif
use wenoof_alpha_object
use wenoof_alpha_rec_js
use wenoof_alpha_rec_z
use wenoof_alpha_rec_m
use wenoof_beta_object
use wenoof_beta_rec_js
use wenoof_interpolator_object
! use wenoof_interpolator_js
use wenoof_interpolations_object
use wenoof_interpolations_rec_js
use wenoof_kappa_object
use wenoof_kappa_rec_js
use wenoof_reconstructor_js
use wenoof_weights_object
use wenoof_weights_js

implicit none
private
public :: interpolator_object
public :: wenoof_create

contains
  subroutine wenoof_create(interpolator_type, S, interpolator, face_left, face_right, eps)
  !< WenOOF creator, create a concrete WENO interpolator object.
  character(*),                            intent(in)           :: interpolator_type          !< Type of the interpolator.
  integer(I_P),                            intent(in)           :: S                          !< Stencil dimension.
  class(interpolator_object), allocatable, intent(out)          :: interpolator               !< The concrete WENO interpolator.
  logical,                                 intent(in), optional :: face_left                  !< Activate left-face interpolations.
  logical,                                 intent(in), optional :: face_right                 !< Activate right-face interpolations.
  real(RPP),                               intent(in), optional :: eps                        !< Small epsilon to avoid zero-div.
  class(alpha_object_constructor),          allocatable         :: alpha_constructor          !< Alpha constructor.
  class(beta_object_constructor),           allocatable         :: beta_constructor           !< Beta constructor.
  class(interpolations_object_constructor), allocatable         :: interpolations_constructor !< Interpolations constructor.
  class(kappa_object_constructor),          allocatable         :: kappa_constructor          !< Kappa constructor.
  class(weights_object_constructor),        allocatable         :: weights_constructor        !< Weights constructor.
  class(interpolator_object_constructor),   allocatable         :: interpolator_constructor   !< Interpolator constructor.

  select case(trim(adjustl(interpolator_type)))
  case('interpolator-JS')
    ! @TODO implement this
    error stop 'interpolator-JS to be implemented'
  case('reconstructor-JS')
    call create_alpha_rec_js_constructor(S=S,                           &
                                         constructor=alpha_constructor, &
                                         face_left=face_left,           &
                                         face_right=face_right,         &
                                         eps=eps)
    call create_beta_rec_js_constructor(S=S,                          &
                                        constructor=beta_constructor, &
                                        face_left=face_left,          &
                                        face_right=face_right)
    call create_interpolations_rec_js_constructor(S=S,                                    &
                                                  constructor=interpolations_constructor, &
                                                  face_left=face_left,                    &
                                                  face_right=face_right)
    call create_kappa_rec_js_constructor(S=S, constructor=kappa_constructor)
    call create_weights_js_constructor(S=S,                                 &
                                       alpha_constructor=alpha_constructor, &
                                       beta_constructor=beta_constructor,   &
                                       kappa_constructor=kappa_constructor, &
                                       constructor=weights_constructor,     &
                                       face_left=face_left,                 &
                                       face_right=face_right,               &
                                       eps=eps)
    call create_reconstructor_js_constructor(S=S,                                                   &
                                             interpolations_constructor=interpolations_constructor, &
                                             weights_constructor=weights_constructor,               &
                                             constructor=interpolator_constructor,                  &
                                             face_left=face_left,                                   &
                                             face_right=face_right,                                 &
                                             eps=eps)
    allocate(reconstructor_js :: interpolator)
    call interpolator%create(constructor=interpolator_constructor)
  case('JS-Z')
    ! @TODO add Z support
    error stop 'JS-Z to be implemented'
  case('JS-M')
    ! @TODO add M support
    error stop 'JS-M to be implemented'
  case default
    ! @TODO add error handling
  endselect
  endsubroutine wenoof_create
endmodule wenoof
