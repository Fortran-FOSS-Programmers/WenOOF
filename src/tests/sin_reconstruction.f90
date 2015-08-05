program sin_reconstruction
!-----------------------------------------------------------------------------------------------------------------------------------
!< Reconstruct sin function.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use IR_Precision, only : I_P, R_P, FR_P
use wenoof, only : weno_factory, weno_constructor_upwind, weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(weno_factory)                    :: factory                    !< WENO factory.
class(weno_interpolator), allocatable :: interpolator               !< WENO interpolator.
integer(I_P), parameter               :: S = 3_I_P                  !< Stencils used.
integer(I_P), parameter               :: Nv = 30_I_P                !< Number of discretized values to be interpolated.
real(R_P),    parameter               :: pi = 4._R_P * atan(1._R_P) !< Extent of domain.
real(R_P)                             :: x(1-S:Nv+S)                !< Whole domain.
real(R_P)                             :: fx(1-S:Nv+S)               !< Discretized values to be interpolated.
real(R_P)                             :: xi(1:Nv)                   !< Domain of the interpolation.
real(R_P)                             :: fx_ref(1:Nv)               !< Reference values.
real(R_P)                             :: interpolation(1:1, 1:Nv)   !< Interpolated values.
integer                               :: i, j, f                    !< Counters.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! build the values used for the reconstruction of sin function: nodal values
x = 0.
do i = 1 - S, Nv + S
  x(i) = i * 2 * pi / Nv
  fx(i) = sin(x(i))
enddo
! face values to which the reconstruction should tend
do i = 1, Nv
  xi(i) = x(i) + pi / Nv
  fx_ref(i) = sin(xi(i))
enddo

! prepare the weno interpolator
call factory%create(constructor=weno_constructor_upwind(S=S), interpolator=interpolator)

! interpolate values
interpolation = 0.
do i = 1, Nv ! interpolated values loop
  call interpolator%interpolate(S=S,                                                      &
                                stencil=reshape(source=fx(i+1-S:i-1+S), shape=[1,2*S-1]), &
                                location='right',                                         &
                                interpolation=interpolation(1:1, i))
enddo

! print results
print "(A)", '# x, sin(x), weno_interpolation(x)'
do i = 1, Nv
  print "(3("//FR_P//", 1X))", xi(i), fx_ref(i), interpolation(1, i)
enddo
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram sin_reconstruction

