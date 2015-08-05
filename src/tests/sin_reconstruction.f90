program sin_reconstruction
!-----------------------------------------------------------------------------------------------------------------------------------
!< Reconstruct sin function.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use wenoof, only : weno_factory, weno_constructor_upwind, weno_interpolator
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(weno_factory)                    :: factory                  !< WENO factory.
class(weno_interpolator), allocatable :: interpolator             !< WENO interpolator.
integer, parameter                    :: S = 3                    !< Stencils used.
integer, parameter                    :: Nv = 30                  !< Number of discretized values to be interpolated.
real,    parameter                    :: pi = 4. * atan(1.)       !< Extent of domain.
real                                  :: x(1-S:Nv+S)              !< Whole domain.
real                                  :: fx(1-S:Nv+S)             !< Discretized values to be interpolated.
real                                  :: xi(1:Nv)                 !< Domain of the interpolation.
real                                  :: fx_ref(1:Nv)             !< Reference values.
real                                  :: interpolation(1:1, 1:Nv) !< Interpolated values.
integer                               :: i, j, f                  !< Counters.
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
  print "(3(E13.6E2, 1X))", xi(i), fx_ref(i), interpolation(1, i)
enddo
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram sin_reconstruction

