module type_weno_optimal_weights_js
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Jiang-Shu and Gerolymos-Sénéchal-Vallet optimal weights for WENO schemes.
!<
!< @note The provided WENO optimal weights implements the optimal weights defined in *Efficient Implementation of Weighted ENO
!< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
!< *Very-high-order weno schemes*, G. A. Gerolymos, D. Sénéchal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
!< doi:10.1016/j.jcp.2009.07.039
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use penf, only : I_P, R_P
use type_weno_optimal_weights
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_optimal_weights_js
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_optimal_weights):: weno_optimal_weights_js
  !< Jiang-Shu and Gerolymos-Sénéchal-Vallet WENO optimal weights object.
  !<
  !< @note The provided WENO optimal weights implements the optimal weights defined in *Efficient Implementation of Weighted ENO
  !< Schemes*, Guang-Shan Jiang, Chi-Wang Shu, JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and
  !< *Very-high-order weno schemes*, G. A. Gerolymos, D. Sénéchal, I. Vallet, JCP, 2009, vol. 228, pp. 8481-8524,
  !< doi:10.1016/j.jcp.2009.07.039
  private
  real(R_P), allocatable :: opt(:,:)   !< Optimal weights                    [1:2,0:S-1].
  contains
    ! deferred public methods
    procedure, pass(self), public :: destroy
    procedure, pass(self), public :: create
    procedure, pass(self), public :: description
endtype weno_optimal_weights_js
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! deferred public methods
  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy Jiang-Shu and Gerolymos-Sénéchal-Vallet WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_optimal_weights_js), intent(inout)  :: self   !< WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%opt)) deallocate(self%opt)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure subroutine create(self,S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_optimal_weights_js), intent(inout) :: self       !< WENO Jiang-Shu optimal weights.
  integer(I_P),                   intent(in)    :: S          !< Number of stencils used.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  allocate(self%opt(1:2, 0:S - 1))
  associate(opt => self%opt)
    select case(S)
      case(2) ! 3rd order
        ! 1 => left interface (i-1/2)
        opt(1, 0) = 2._R_P/3._R_P ! stencil 0
        opt(1, 1) = 1._R_P/3._R_P ! stencil 1
        ! 2 => right interface (i+1/2)
        opt(2, 0) = 1._R_P/3._R_P ! stencil 0
        opt(2, 1) = 2._R_P/3._R_P ! stencil 1
      case(3) ! 5th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) = 0.3_R_P ! stencil 0
        opt(1, 1) = 0.6_R_P ! stencil 1
        opt(1, 2) = 0.1_R_P ! stencil 2
        ! 2 => right interface (i+1/2)
        opt(2, 0) = 0.1_R_P ! stencil 0
        opt(2, 1) = 0.6_R_P ! stencil 1
        opt(2, 2) = 0.3_R_P ! stencil 2
      case(4) ! 7th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =  4._R_P/35._R_P ! stencil 0
        opt(1, 1) = 18._R_P/35._R_P ! stencil 1
        opt(1, 2) = 12._R_P/35._R_P ! stencil 2
        opt(1, 3) =  1._R_P/35._R_P ! stencil 3
        ! 2 => right interface (i+1/2)
        opt(2, 0) =  1._R_P/35._R_P ! stencil 0
        opt(2, 1) = 12._R_P/35._R_P ! stencil 1
        opt(2, 2) = 18._R_P/35._R_P ! stencil 2
        opt(2, 3) =  4._R_P/35._R_P ! stencil 3
      case(5) ! 9th order
        ! 1 => left interface (i-1/2)
        opt(1, 0) =  5._R_P/126._R_P ! stencil 0
        opt(1, 1) = 20._R_P/63._R_P  ! stencil 1
        opt(1, 2) = 10._R_P/21._R_P  ! stencil 2
        opt(1, 3) = 10._R_P/63._R_P  ! stencil 3
        opt(1, 4) =  1._R_P/126._R_P ! stencil 4
        ! 2 => right interface (i+1/2)
        opt(2, 0) =  1._R_P/126._R_P ! stencil 0
        opt(2, 1) = 10._R_P/63._R_P  ! stencil 1
        opt(2, 2) = 10._R_P/21._R_P  ! stencil 2
        opt(2, 3) = 20._R_P/63._R_P  ! stencil 3
        opt(2, 4) =  5._R_P/126._R_P ! stencil 4
      case(6) ! 11th order
        ! 1 => left interface (i-1/2)
        opt(2, 0) =   1._R_P/77._R_P  ! stencil 0
        opt(2, 1) =  25._R_P/154._R_P ! stencil 1
        opt(2, 2) = 100._R_P/231._R_P ! stencil 2
        opt(2, 3) =  25._R_P/77._R_P  ! stencil 3
        opt(2, 4) =   5._R_P/77._R_P  ! stencil 4
        opt(2, 5) =   1._R_P/462._R_P ! stencil 5
        ! 2 => right interface (i+1/2)
        opt(2, 0) =   1._R_P/462._R_P ! stencil 0
        opt(2, 1) =   5._R_P/77._R_P  ! stencil 1
        opt(2, 2) =  25._R_P/77._R_P  ! stencil 2
        opt(2, 3) = 100._R_P/231._R_P ! stencil 3
        opt(2, 4) =  25._R_P/154._R_P ! stencil 4
        opt(2, 5) =   1._R_P/77._R_P  ! stencil 5
      case(7) ! 13th order
        ! 1 => left interface (i-1/2)
        opt(2, 0) =   7._R_P/1716._R_P ! stencil 0
        opt(2, 1) =  21._R_P/286._R_P  ! stencil 1
        opt(2, 2) = 175._R_P/572._R_P  ! stencil 2
        opt(2, 3) = 175._R_P/429._R_P  ! stencil 3
        opt(2, 4) = 105._R_P/572._R_P  ! stencil 4
        opt(2, 5) =   7._R_P/286._R_P  ! stencil 5
        opt(2, 6) =   1._R_P/1716._R_P ! stencil 6
        ! 2 => right interface (i+1/2)
        opt(2, 0) =   1._R_P/1716._R_P ! stencil 0
        opt(2, 1) =   7._R_P/286._R_P  ! stencil 1
        opt(2, 2) = 105._R_P/572._R_P  ! stencil 2
        opt(2, 3) = 175._R_P/429._R_P  ! stencil 3
        opt(2, 4) = 175._R_P/572._R_P  ! stencil 4
        opt(2, 5) =  21._R_P/286._R_P  ! stencil 5
        opt(2, 6) =   7._R_P/1716._R_P ! stencil 6
      case(8) ! 15th order
        ! 1 => left interface (i-1/2)
        opt(2, 0) =   8._R_P/6435._R_P ! stencil 0
        opt(2, 1) = 196._R_P/6435._R_P ! stencil 1
        opt(2, 2) = 392._R_P/2145._R_P ! stencil 2
        opt(2, 3) = 490._R_P/1287._R_P ! stencil 3
        opt(2, 4) = 392._R_P/1287._R_P ! stencil 4
        opt(2, 5) = 196._R_P/2145._R_P ! stencil 5
        opt(2, 6) =  56._R_P/6435._R_P ! stencil 6
        opt(2, 7) =   1._R_P/6435._R_P ! stencil 7
        ! 2 => right interface (i+1/2)
        opt(2, 0) =   1._R_P/6435._R_P ! stencil 0
        opt(2, 1) =  56._R_P/6435._R_P ! stencil 1
        opt(2, 2) = 196._R_P/2145._R_P ! stencil 2
        opt(2, 3) = 392._R_P/1287._R_P ! stencil 3
        opt(2, 4) = 490._R_P/1287._R_P ! stencil 4
        opt(2, 5) = 392._R_P/2145._R_P ! stencil 5
        opt(2, 6) = 196._R_P/6435._R_P ! stencil 6
        opt(2, 7) =   8._R_P/6435._R_P ! stencil 7
      case(9) ! 17th order
        ! 1 => left interface (i-1/2)
        opt(2, 0) =    1._R_P/24310._R_P ! stencil 0
        opt(2, 1) =    9._R_P/24310._R_P ! stencil 1
        opt(2, 2) =  144._R_P/12155._R_P ! stencil 2
        opt(2, 3) = 1176._R_P/12155._R_P ! stencil 3
        opt(2, 4) = 3528._R_P/12155._R_P ! stencil 4
        opt(2, 5) =  882._R_P/2431._R_P  ! stencil 5
        opt(2, 6) = 2352._R_P/12155._R_P ! stencil 6
        opt(2, 7) =  504._R_P/12155._R_P ! stencil 7
        opt(2, 8) =   36._R_P/12155._R_P ! stencil 8
        ! 2 => right interface (i+1/2)
        opt(2, 0) =    1._R_P/24310._R_P ! stencil 0
        opt(2, 1) =   36._R_P/12155._R_P ! stencil 1
        opt(2, 2) =  504._R_P/12155._R_P ! stencil 2
        opt(2, 3) = 2352._R_P/12155._R_P ! stencil 3
        opt(2, 4) =  882._R_P/2431._R_P  ! stencil 4
        opt(2, 5) = 3528._R_P/12155._R_P ! stencil 5
        opt(2, 6) = 1176._R_P/12155._R_P ! stencil 6
        opt(2, 7) =  144._R_P/12155._R_P ! stencil 7
        opt(2, 8) =    9._R_P/24310._R_P ! stencil 8
    endselect
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create

  pure subroutine description(self, string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing Jiang-Shu and Gerolymos-Sénéchal-Vallet WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_optimal_weights_js), intent(in)  :: self   !< WENO optimal weights.
  character(len=:), allocatable,  intent(out) :: string !< String returned.
  character(len=1), parameter                 :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO optimal weights,'//nl
  string = string//'  based on the work by Jiang and Shu "Efficient Implementation of Weighted ENO Schemes", see '// &
           'JCP, 1996, vol. 126, pp. 202--228, doi:10.1006/jcph.1996.0130 and'//nl
  string = string//'  on the work by Gerolymos, Sénéchal and Vallet "Very-high-order weno schemes", see '// &
           'JCP, 2009, vol. 228, pp. 8481--8524, doi:10.1016/j.jcp.2009.07.039'//nl
  string = string//'    The optimal weights are allocated in a two-dimensional array, in which the first index'//nl
  string = string//'    is the face selected (1 => i-1/2, 2 => i+1/2) and the second index is the number of the stencil '//nl
  string = string//'    (from 0 to S-1)'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_optimal_weights_js
