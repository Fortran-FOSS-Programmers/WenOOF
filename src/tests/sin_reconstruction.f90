program sin_reconstruction
!-----------------------------------------------------------------------------------------------------------------------------------
!< Reconstruct sin function.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use wenoof
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(weno_constructor_upwind)         :: upwind_constructor       !< WENO constructor.
class(weno_interpolator), allocatable :: weno                     !< WENO interpolator.
integer, parameter                    :: S = 3                    !< Stencils used.
integer, parameter                    :: Nv = 30                  !< Number of discretized values to be interpolated.
real,    parameter                    :: pi = 4. * atan(1.)       !< Extent of domain.
real                                  :: x(1-S:Nv+S)              !< Whole domain.
real                                  :: fx(1-S:Nv+S)             !< Discretized values to be interpolated.
real                                  :: stencil(1:2, 1-S:-1+S)   !< Left and right interpolated values.
real                                  :: xi(1:Nv)                 !< Domain of the interpolation.
real                                  :: fx_ref(1:Nv)             !< Reference values.
real                                  :: interpolation(1:2, 1:Nv) !< Left and right interpolated values.
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
allocate(weno_interpolator_upwind :: weno) ! this could be avoid if weno is defined of type(weno_interpolator_upwind)
upwind_constructor = weno_constructor_upwind(S=S)
call weno%create(constructor=upwind_constructor)

! interpolate values
interpolation = 0.
select type(weno)
type is(weno_interpolator_upwind)
  do i = 1, Nv ! interpolated values loop
    ! prepare stencils
    do j = i + 1 - S, i - 1 + S ! stencil values loop
      do f = 1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        stencil(f, j - i) = fx(j)
      enddo
    enddo
    ! interpolate left and right
    do f = 1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      call weno%interpolate(S=S, stencil=stencil, interpolation=interpolation(1:2, i))
    enddo
  enddo
endselect

! print results
print "(A)", '# x, sin(x), weno_interpolation(x)'
do i = 1, Nv
  print "(3(E13.6E2, 1X))", xi(i), fx_ref(i), interpolation(2, i)
enddo
stop
!-----------------------------------------------------------------------------------------------------------------------------------
endprogram sin_reconstruction

