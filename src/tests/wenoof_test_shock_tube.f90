!< WenOOF test: shock tube tester, 1D Euler equation.

module foreseer_euler_1d
!< Definition of Euler 1D class for WenOOF test.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use flow, only : conservative_object, conservative_compressible, primitive_compressible,         &
                 conservative_to_primitive_compressible, primitive_to_conservative_compressible, &
                 eos_object, eos_compressible
use foreseer, only : riemann_solver_object, riemann_solver_compressible_exact,          &
                     riemann_solver_compressible_hllc, riemann_solver_compressible_llf, &
                     riemann_solver_compressible_pvl, riemann_solver_compressible_roe
use penf, only : I4P, R8P
use foodie, only : integrand
use vecfor, only : ex, vector
use wenoof, only : interpolator_object, wenoof_create

implicit none
private
public :: euler_1d

type, extends(integrand) :: euler_1d
   !< Euler 1D PDEs system field.
   !<
   !< It is a FOODIE integrand class concrete extension.
   !<
   !<### 1D Euler PDEs system
   !< The 1D Euler PDEs system considered is a non linear, hyperbolic (inviscid) system of conservation laws for compressible gas
   !< dynamics, that reads as
   !<$$
   !<\begin{matrix}
   !<U_t = R(U)  \Leftrightarrow U_t = F(U)_x \\
   !<U = \begin{bmatrix}
   !<\rho \\
   !<\rho u \\
   !<\rho E
   !<\end{bmatrix}\;\;\;
   !<F(U) = \begin{bmatrix}
   !<\rho u \\
   !<\rho u^2 + p \\
   !<\rho u H
   !<\end{bmatrix}
   !<\end{matrix}
   !<$$
   !< where \(\rho\) is the density, \(u\) is the velocity, \(p\) the pressure, \(E\) the total internal specific energy and \(H\)
   !< the total specific enthalpy. The PDEs system must completed with the proper initial and boundary conditions. Moreover, an
   !< ideal (thermally and calorically perfect) gas is considered
   !<$$
   !<\begin{matrix}
   !<R = c_p - c_v \\
   !<\gamma = \frac{c_p}{c_v}\\
   !<e = c_v T \\
   !<h = c_p T
   !<\end{matrix}
   !<$$
   !< where *R* is the gas constant, \(c_p\,c_v\) are the specific heats at constant pressure and volume (respectively), *e* is the
   !< internal energy, *h* is the internal enthalpy and *T* is the temperature. The following addition equations of state hold:
   !<$$
   !<\begin{matrix}
   !<T = \frac{p}{\rho R} \\
   !<E = \rho e + \frac{1}{2} \rho u^2 \\
   !<H = \rho h + \frac{1}{2} \rho u^2 \\
   !<a = \sqrt{\frac{\gamma p}{\rho}}
   !<\end{matrix}
   !<$$
   !<
   !<#### Numerical grid organization
   !< The finite volume, Godunov's like approach is employed. The conservative variables (and the primitive ones) are co-located at
   !< the cell center. The cell and (inter)faces numeration is as follow.
   !<```
   !<                cell            (inter)faces
   !<                 |                   |
   !<                 v                   v
   !<     |-------|-------|-.....-|-------|-------|-------|-------|-.....-|-------|-------|-------|-.....-|-------|-------|
   !<     | 1-Ng  | 2-Ng  | ..... |  -1   |   0   |   1   |  2    | ..... |  Ni   | Ni+1  | Ni+1  | ..... |Ni+Ng-1| Ni+Ng |
   !<     |-------|-------|-.....-|-------|-------|-------|-------|-.....-|-------|-------|-------|-.....-|-------|-------|
   !<    0-Ng                             -1      0       1       2      Ni-1     Ni                                    Ni+Ng
   !<```
   !< Where *Ni* are the finite volumes (cells) used for discretizing the domain and *Ng* are the ghost cells used for imposing the
   !< left and right boundary conditions (for a total of *2Ng* cells).
   integer(I4P)                                 :: weno_order=0                          !< WENO reconstruction order.
   integer(I4P)                                 :: Ni=0                                  !< Space dimension.
   integer(I4P)                                 :: Ng=0                                  !< Ghost cells number.
   real(R8P)                                    :: Dx=0._R8P                             !< Space step.
   type(eos_compressible)                       :: eos                                   !< Equation of state.
   type(conservative_compressible), allocatable :: U(:)                                  !< Integrand (state) variables.
   character(:),                    allocatable :: BC_L                                  !< Left boundary condition type.
   character(:),                    allocatable :: BC_R                                  !< Right boundary condition type.
   class(interpolator_object),      allocatable :: interpolator                          !< WENO interpolator.
   procedure(reconstruct_interfaces_), pointer  :: reconstruct_interfaces=>&
                                                   reconstruct_interfaces_characteristic !< Reconstruct interface states.
   class(riemann_solver_object), allocatable    :: riemann_solver                        !< Riemann solver.
   contains
      ! auxiliary methods
      procedure, pass(self) :: initialize       !< Initialize field.
      procedure, pass(self) :: destroy          !< Destroy field.
      procedure, pass(self) :: output           !< Extract Euler field.
      procedure, pass(self) :: dt => compute_dt !< Compute the current time step, by means of CFL condition.
      ! ADT integrand deferred methods
      procedure, pass(self) :: t => dEuler_dt                                       !< Time derivative, residuals function.
      procedure, pass(lhs)  :: local_error => euler_local_error                     !< Operator `||euler-euler||`.
      procedure, pass(lhs)  :: integrand_multiply_integrand => euler_multiply_euler !< Operator `*`.
      procedure, pass(lhs)  :: integrand_multiply_real => euler_multiply_real       !< Operator `euler * real`.
      procedure, pass(rhs)  :: real_multiply_integrand => real_multiply_euler       !< Operator `real * euler`.
      procedure, pass(lhs)  :: add => add_euler                                     !< Operator `+`.
      procedure, pass(lhs)  :: sub => sub_euler                                     !< Operator `-`.
      procedure, pass(lhs)  :: assign_integrand => euler_assign_euler               !< Operator `=`.
      procedure, pass(lhs)  :: assign_real => euler_assign_real                     !< Operator `euler = real`.
      ! private methods
      procedure, pass(self), private :: impose_boundary_conditions            !< Impose boundary conditions.
      procedure, pass(self), private :: reconstruct_interfaces_characteristic !< Reconstruct (charc.) interface states.
      procedure, pass(self), private :: reconstruct_interfaces_conservative   !< Reconstruct (cons.) interface states.
      procedure, pass(self), private :: reconstruct_interfaces_primitive      !< Reconstruct (prim.) interface states.
endtype euler_1d

abstract interface
   !< Abstract interfaces of [[euler_1d]] pointer methods.
   subroutine reconstruct_interfaces_(self, conservative, r_conservative)
   !< Reconstruct interface states.
   import :: conservative_compressible, euler_1d, primitive_compressible
   class(euler_1d),                 intent(in)    :: self                     !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:) !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)   !< Reconstructed conservative variables.
   endsubroutine reconstruct_interfaces_
endinterface

contains
   ! auxiliary methods
   subroutine initialize(self, Ni, Dx, BC_L, BC_R, initial_state, eos, weno_order, weno_variables, riemann_solver_scheme)
   !< Initialize field.
   class(euler_1d),              intent(inout)        :: self                   !< Euler field.
   integer(I4P),                 intent(in)           :: Ni                     !< Space dimension.
   real(R8P),                    intent(in)           :: Dx                     !< Space step.
   character(*),                 intent(in)           :: BC_L                   !< Left boundary condition type.
   character(*),                 intent(in)           :: BC_R                   !< Right boundary condition type.
   type(primitive_compressible), intent(in)           :: initial_state(1:)      !< Initial state of primitive variables.
   type(eos_compressible),       intent(in)           :: eos                    !< Equation of state.
   integer(I4P),                 intent(in), optional :: weno_order             !< WENO reconstruction order.
   character(*),                 intent(in), optional :: weno_variables         !< Variables on which WENO reconstruction is done.
   character(*),                 intent(in), optional :: riemann_solver_scheme  !< Riemann solver scheme.
   character(:), allocatable                          :: weno_variables_        !< WENO Variables, local variable.
   character(:), allocatable                          :: riemann_solver_scheme_ !< Riemann solver scheme, local variable.
   integer(I4P)                                       :: i                      !< Space couner.

   call self%destroy
   self%weno_order = 1 ; if (present(weno_order)) self%weno_order = weno_order
   self%Ni = Ni
   self%Ng = (self%weno_order + 1) / 2
   self%Dx = Dx
   self%eos = eos
   if (allocated(self%U)) deallocate(self%U) ; allocate(self%U(1-self%Ng:self%Ni+self%Ng))
   do i=1, Ni
      self%U(i) = primitive_to_conservative_compressible(primitive=initial_state(i), eos=eos)
   enddo
   self%BC_L = BC_L
   self%BC_R = BC_R

   if (self%weno_order>1) call wenoof_create(interpolator_type='reconstructor-JS', S=self%Ng, interpolator=self%interpolator)
   weno_variables_ = 'characteristic'
   if (present(weno_variables)) weno_variables_ = trim(adjustl(weno_variables))
   select case(weno_variables_)
   case('characteristic')
      self%reconstruct_interfaces => reconstruct_interfaces_characteristic
   case('conservative')
      self%reconstruct_interfaces => reconstruct_interfaces_conservative
   case('primitive')
      self%reconstruct_interfaces => reconstruct_interfaces_primitive
   case default
      write(stderr, '(A)') 'error: WENO reconstruction variables set "'//weno_variables_//'" unknown!'
      stop
   endselect

   riemann_solver_scheme_ = 'llf'
   if (present(riemann_solver_scheme)) riemann_solver_scheme_ = trim(adjustl(riemann_solver_scheme))
   select case(riemann_solver_scheme_)
   case('exact')
      allocate(riemann_solver_compressible_exact :: self%riemann_solver)
   case('hllc')
      allocate(riemann_solver_compressible_hllc :: self%riemann_solver)
   case('llf')
      allocate(riemann_solver_compressible_llf :: self%riemann_solver)
   case('pvl')
      allocate(riemann_solver_compressible_pvl :: self%riemann_solver)
   case('roe')
      allocate(riemann_solver_compressible_roe :: self%riemann_solver)
   case default
      write(stderr, '(A)') 'error: Riemann Solver scheme "'//riemann_solver_scheme_//'" unknown!'
      stop
   endselect
   endsubroutine initialize

   pure subroutine destroy(self)
   !< Destroy field.
   class(euler_1d), intent(inout) :: self !< Euler field.

   self%weno_order = 0
   self%Ni = 0
   self%Ng = 0
   self%Dx = 0._R8P
   if (allocated(self%U)) deallocate(self%U)
   if (allocated(self%BC_L)) deallocate(self%BC_L)
   if (allocated(self%BC_R)) deallocate(self%BC_R)
   if (allocated(self%interpolator)) deallocate(self%interpolator)
   if (allocated(self%riemann_solver)) deallocate(self%riemann_solver)
   endsubroutine destroy

   pure function output(self, is_primitive) result(state)
   !< Output the Euler field state.
   class(euler_1d), intent(in)           :: self          !< Euler field.
   logical,         intent(in), optional :: is_primitive  !< Output in primitive variables.
   real(R8P), allocatable                :: state(:,:)    !< Euler state vector.
   real(R8P), allocatable                :: state_(:)     !< Euler state vector, local variable.
   logical                               :: is_primitive_ !< Output in primitive variables, local variable.
   type(primitive_compressible)          :: primitive     !< Primitive state.
   integer(I4P)                          :: i             !< Counter.

   is_primitive_ = .false. ; if (present(is_primitive)) is_primitive_ = is_primitive
   if (is_primitive_) then
      allocate(state(1:size(primitive%array(), dim=1), 1:self%Ni))
      do i=1, self%Ni
         primitive = conservative_to_primitive_compressible(conservative=self%U(i), eos=self%eos)
         state_ = primitive%array()
      enddo
   else
      allocate(state(1:size(self%U(1)%array(), dim=1), 1:self%Ni))
      do i=1, self%Ni
         state_ = self%U(i)%array()
         state(:, i) = state_
      enddo
   endif
   endfunction output

   pure function compute_dt(self, steps_max, t_max, t, CFL) result(Dt)
   !< Compute the current time step by means of CFL condition.
   class(euler_1d), intent(in) :: self      !< Euler field.
   integer(I4P),    intent(in) :: steps_max !< Maximun number of time steps.
   real(R8P),       intent(in) :: t_max     !< Maximum integration time.
   real(R8P),       intent(in) :: t         !< Time.
   real(R8P),       intent(in) :: CFL       !< CFL value.
   real(R8P)                   :: Dt        !< Time step.
   type(vector)                :: u         !< Velocity vector.
   real(R8P)                   :: a         !< Speed of sound.
   real(R8P)                   :: vmax      !< Maximum propagation speed of signals.
   integer(I4P)                :: i         !< Counter.

   associate(Ni=>self%Ni, Dx=>self%Dx)
      vmax = 0._R8P
      do i=1, Ni
         u = self%U(i)%velocity()
         a = self%eos%speed_of_sound(density=self%U(i)%density, pressure=self%U(i)%pressure(eos=self%eos))
         vmax = max(vmax, u%normL2() + a)
      enddo
      Dt = Dx * CFL / vmax
      if (steps_max <= 0 .and. t_max > 0._R8P) then
         if ((t + Dt) > t_max) Dt = t_max - t
      endif
   endassociate
   endfunction compute_dt

   ! ADT integrand deferred methods
   function dEuler_dt(self, t) result(dState_dt)
   !< Time derivative of Euler field, the residuals function.
   class(euler_1d), intent(in)           :: self                         !< Euler field.
   real(R8P),       intent(in), optional :: t                            !< Time.
   class(integrand), allocatable         :: dState_dt                    !< Euler field time derivative.
   type(conservative_compressible)       :: U(1-self%Ng:self%Ni+self%Ng) !< Conservative variables.
   type(conservative_compressible)       :: UR(1:2,0:self%Ni+1)          !< Reconstructed conservative variables.
   type(conservative_compressible)       :: F(0:self%Ni)                 !< Fluxes of conservative variables.
   integer(I4P)                          :: i                            !< Counter.

   do i=1, self%Ni
      U(i) = self%U(i)
   enddo
   call self%impose_boundary_conditions(U=U)
   call self%reconstruct_interfaces(conservative=U, r_conservative=UR)
   do i=0, self%Ni
      call self%riemann_solver%solve(eos_left=self%eos,  state_left=UR( 2, i  ), &
                                     eos_right=self%eos, state_right=UR(1, i+1), normal=ex, fluxes=F(i))
   enddo
   allocate(euler_1d :: dState_dt)
   select type(dState_dt)
   class is(euler_1d)
      dState_dt = self
      do i=1, self%Ni
          dState_dt%U(i) = (F(i - 1) - F(i)) / self%Dx
      enddo
   endselect
   endfunction dEuler_dt

  function euler_local_error(lhs, rhs) result(error)
  !< Estimate local truncation error between 2 euler approximations.
  !<
  !< The estimation is done by norm L2 of U:
  !<
  !< $$ error = \sqrt{ \sum_i{\sum_i{ \frac{(lhs\%U_i - rhs\%U_i)^2}{lhs\%U_i^2} }} } $$
  class(euler_1d),  intent(in) :: lhs      !< Left hand side.
  class(integrand), intent(in) :: rhs      !< Right hand side.
  real(R8P)                    :: error    !< Error estimation.
  integer(I4P)                 :: i        !< Space counter.

  select type(rhs)
  class is (euler_1d)
     error = 0._R8P
     do i=1, lhs%Ni
        error = error + (lhs%U(i)%density - rhs%U(i)%density) ** 2 / lhs%U(i)%density ** 2
     enddo
     error = sqrt(error)
  endselect
  endfunction euler_local_error

  function euler_multiply_euler(lhs, rhs) result(opr)
  !< Multiply an Euler field by another one.
  class(euler_1d),  intent(in)  :: lhs !< Left hand side.
  class(integrand), intent(in)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate(euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
     select type(rhs)
     class is (euler_1d)
        do i=1, lhs%Ni
           opr%U(i) = lhs%U(i) * rhs%U(i)
        enddo
     endselect
  endselect
  endfunction euler_multiply_euler

  function euler_multiply_real(lhs, rhs) result(opr)
  !< Multiply an Euler field by a real scalar.
  class(euler_1d), intent(in)   :: lhs !< Left hand side.
  real(R8P),       intent(in)   :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate(euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
     do i=1, lhs%Ni
        opr%U(i) = rhs * lhs%U(i)
     enddo
  endselect
  endfunction euler_multiply_real

  function real_multiply_euler(lhs, rhs) result(opr)
  !< Multiply a real scalar by an Euler field.
  real(R8P),       intent(in)   :: lhs !< Left hand side.
  class(euler_1d), intent(in)   :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate(euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = rhs
     do i=1, rhs%Ni
        opr%U(i) = lhs * rhs%U(i)
     enddo
  endselect
  endfunction real_multiply_euler

  function add_euler(lhs, rhs) result(opr)
  !< Add two Euler fields.
  class(euler_1d),  intent(in)  :: lhs !< Left hand side.
  class(integrand), intent(in)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate (euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
     select type(rhs)
     class is (euler_1d)
        do i=1, lhs%Ni
           opr%U(i) = lhs%U(i) + rhs%U(i)
        enddo
     endselect
  endselect
  endfunction add_euler

  function sub_euler(lhs, rhs) result(opr)
  !< Subtract two Euler fields.
  class(euler_1d),  intent(in)  :: lhs !< Left hand side.
  class(integrand), intent(in)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate (euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
     select type(rhs)
     class is (euler_1d)
        do i=1, lhs%Ni
           opr%U(i) = lhs%U(i) - rhs%U(i)
        enddo
     endselect
  endselect
  endfunction sub_euler

  subroutine euler_assign_euler(lhs, rhs)
  !< Assign one Euler field to another.
  class(euler_1d),  intent(inout) :: lhs !< Left hand side.
  class(integrand), intent(in)    :: rhs !< Right hand side.
  integer(I4P)                    :: i   !< Counter.

  select type(rhs)
  class is(euler_1d)
     lhs%weno_order = rhs%weno_order
     lhs%Ni         = rhs%Ni
     lhs%Ng         = rhs%Ng
     lhs%Dx         = rhs%Dx
     lhs%eos        = rhs%eos
     if (allocated(rhs%U)) then
        if (allocated(lhs%U)) deallocate(lhs%U) ; allocate(lhs%U(1-lhs%Ng:lhs%Ni+lhs%Ng))
        do i=1, lhs%Ni
           lhs%U(i) = rhs%U(i)
        enddo
     endif
     if (allocated(rhs%BC_L)) lhs%BC_L = rhs%BC_L
     if (allocated(rhs%BC_R)) lhs%BC_R = rhs%BC_R
     if (allocated(rhs%interpolator)) then
        if (allocated(lhs%interpolator)) deallocate(lhs%interpolator)
        allocate(lhs%interpolator, source=rhs%interpolator)
     endif
     if (associated(rhs%reconstruct_interfaces)) lhs%reconstruct_interfaces => rhs%reconstruct_interfaces
     if (allocated(rhs%riemann_solver)) then
        if (allocated(lhs%riemann_solver)) deallocate(lhs%riemann_solver) ; allocate(lhs%riemann_solver, source=rhs%riemann_solver)
     endif
  endselect
  endsubroutine euler_assign_euler

  subroutine euler_assign_real(lhs, rhs)
  !< Assign one real to an Euler field.
  class(euler_1d), intent(inout) :: lhs !< Left hand side.
  real(R8P),       intent(in)    :: rhs !< Right hand side.
  integer(I4P)                   :: i   !< Counter.

  if (allocated(lhs%U)) then
     do i=1, lhs%Ni
        lhs%U(i)%density = rhs
        lhs%U(i)%momentum = rhs
        lhs%U(i)%energy = rhs
     enddo
  endif
  endsubroutine euler_assign_real

   ! private methods
   pure subroutine impose_boundary_conditions(self, U)
   !< Impose boundary conditions.
   !<
   !< The boundary conditions are imposed on the primitive variables by means of the ghost cells approach.
   class(euler_1d),              intent(in)    :: self          !< Euler field.
   type(conservative_compressible), intent(inout) :: U(1-self%Ng:) !< Conservative variables.
   ! type(primitive_compressible), intent(inout) :: P(1-self%Ng:) !< Primitive variables.
   integer(I4P)                                :: i             !< Space counter.

   select case(trim(adjustl(self%BC_L)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=1-self%Ng, 0
            ! P(i) = P(-i+1)
            U(i) = U(-i+1)
         enddo
      case('REF') ! reflective BC
         do i=1-self%Ng, 0
            ! P(i)%density  =   P(-i+1)%density
            ! P(i)%velocity = - P(-i+1)%velocity
            ! P(i)%pressure =   P(-i+1)%pressure
            U(i)%density  =   U(-i+1)%density
            U(i)%momentum = - U(-i+1)%momentum
            U(i)%energy   =   U(-i+1)%energy
         enddo
   endselect

   select case(trim(adjustl(self%BC_R)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=self%Ni+1, self%Ni+self%Ng
            ! P(i) = P(self%Ni-(i-self%Ni-1))
            U(i) = U(self%Ni-(i-self%Ni-1))
         enddo
      case('REF') ! reflective BC
         do i=self%Ni+1, self%Ni+self%Ng
            ! P(i)%density  =   P(self%Ni-(i-self%Ni-1))%density
            ! P(i)%velocity = - P(self%Ni-(i-self%Ni-1))%velocity
            ! P(i)%pressure =   P(self%Ni-(i-self%Ni-1))%pressure
            U(i)%density  =   U(self%Ni-(i-self%Ni-1))%density
            U(i)%momentum = - U(self%Ni-(i-self%Ni-1))%momentum
            U(i)%energy   =   U(self%Ni-(i-self%Ni-1))%energy
         enddo
   endselect
   endsubroutine impose_boundary_conditions

   subroutine reconstruct_interfaces_characteristic(self, conservative, r_conservative)
   !< Reconstruct interfaces states.
   !<
   !< The reconstruction is done in pseudo characteristic variables.
   class(euler_1d),                 intent(in)    :: self                                 !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:)             !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)               !< Reconstructed conservative vars.
   type(primitive_compressible)                   :: primitive(1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
   type(primitive_compressible)                   :: r_primitive(1:2, 0:self%Ni+1)        !< Reconstructed primitive variables.
   type(primitive_compressible)                   :: Pm(1:2)                              !< Mean of primitive variables.
   real(R8P)                                      :: LPm(1:3, 1:3, 1:2)                   !< Mean left eigenvectors matrix.
   real(R8P)                                      :: RPm(1:3, 1:3, 1:2)                   !< Mean right eigenvectors matrix.
   real(R8P)                                      :: C(1:2, 1-self%Ng:-1+self%Ng, 1:3)    !< Pseudo characteristic variables.
   real(R8P)                                      :: CR(1:2, 1:3)                         !< Pseudo characteristic reconst.
   real(R8P)                                      :: buffer(1:3)                          !< Dummy buffer.
   integer(I4P)                                   :: i                                    !< Counter.
   integer(I4P)                                   :: j                                    !< Counter.
   integer(I4P)                                   :: f                                    !< Counter.
   integer(I4P)                                   :: v                                    !< Counter.
   class(interpolator_object), allocatable        :: interpolator                         !< WENO interpolator.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_conservative(1, i) = conservative(i)
         r_conservative(2, i) = r_conservative(1, i)
      enddo
   case(3, 5, 7, 9, 11, 13, 15, 17) ! 3rd-17th order WENO reconstruction
      call wenoof_create(interpolator_type='reconstructor-JS', S=self%Ng, interpolator=interpolator)
      do i=1-self%Ng, self%Ni+self%Ng
         primitive(i) = conservative_to_primitive_compressible(conservative=conservative(i), eos=self%eos)
      enddo
      do i=0, self%Ni+1
         ! compute pseudo charteristic variables
         do f=1, 2
            Pm(f) = 0.5_R8P * (primitive(i+f-2) + primitive(i+f-1))
         enddo
         do f=1, 2
            LPm(:, :, f) = Pm(f)%left_eigenvectors(eos=self%eos)
            RPm(:, :, f) = Pm(f)%right_eigenvectors(eos=self%eos)
         enddo
         do j=i+1-self%Ng, i-1+self%Ng
            do f=1, 2
               do v=1, 3
                  C(f, j-i, v) = dot_product(LPm(v, :, f), [primitive(j)%density, primitive(j)%velocity%x, primitive(j)%pressure])
               enddo
            enddo
         enddo
         ! compute WENO reconstruction of pseudo charteristic variables
         do v=1, 3
            call interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
         enddo
         ! trasform back reconstructed pseudo charteristic variables to primitive ones
         do f=1, 2
            do v=1, 3
               buffer(v) = dot_product(RPm(v, :, f), CR(f, :))
            enddo
            r_primitive(f, i)%density  = buffer(1)
            r_primitive(f, i)%velocity = buffer(2) * ex
            r_primitive(f, i)%pressure = buffer(3)
         enddo
      enddo
      do i=0, self%Ni+1
         r_conservative(1, i) = primitive_to_conservative_compressible(primitive=r_primitive(1, i), eos=self%eos)
         r_conservative(2, i) = primitive_to_conservative_compressible(primitive=r_primitive(2, i), eos=self%eos)
      enddo
   endselect
   endsubroutine reconstruct_interfaces_characteristic

   subroutine reconstruct_interfaces_conservative(self, conservative, r_conservative)
   !< Reconstruct interfaces states.
   !<
   !< The reconstruction is done in conservative variables.
   class(euler_1d),                 intent(in)    :: self                                 !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:)             !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)               !< Reconstructed conservative vars.
   real(R8P), allocatable                         :: U(:)                                 !< Serialized conservative variables.
   real(R8P)                                      :: C(1:2, 1-self%Ng:-1+self%Ng, 1:3)    !< Pseudo characteristic variables.
   real(R8P)                                      :: CR(1:2, 1:3)                         !< Pseudo characteristic reconst.
   integer(I4P)                                   :: i                                    !< Counter.
   integer(I4P)                                   :: j                                    !< Counter.
   integer(I4P)                                   :: f                                    !< Counter.
   integer(I4P)                                   :: v                                    !< Counter.
   class(interpolator_object), allocatable        :: interpolator                         !< WENO interpolator.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_conservative(1, i) = conservative(i)
         r_conservative(2, i) = r_conservative(1, i)
      enddo
   case(3, 5, 7, 9, 11, 13, 15, 17) ! 3rd-17th order WENO reconstruction
      call wenoof_create(interpolator_type='reconstructor-JS', S=self%Ng, interpolator=interpolator)
      do i=0, self%Ni+1
         do j=i+1-self%Ng, i-1+self%Ng
             U = conservative(j)%array()
            do f=1, 2
               C(f, j-i, 1) = U(1)
               C(f, j-i, 2) = U(2)
               C(f, j-i, 3) = U(5)
            enddo
         enddo
         do v=1, 3
            call interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
         enddo
         do f=1, 2
            r_conservative(f, i)%density  = CR(f, 1)
            r_conservative(f, i)%momentum = CR(f, 2) * ex
            r_conservative(f, i)%energy   = CR(f, 3)
         enddo
      enddo
   endselect
   endsubroutine reconstruct_interfaces_conservative

   subroutine reconstruct_interfaces_primitive(self, conservative, r_conservative)
   !< Reconstruct interfaces states.
   !<
   !< The reconstruction is done in primitive variables.
   class(euler_1d),                 intent(in)    :: self                                 !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:)             !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)               !< Reconstructed conservative vars.
   type(primitive_compressible)                   :: primitive(1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
   type(primitive_compressible)                   :: r_primitive(1:2, 0:self%Ni+1)        !< Reconstructed primitive variables.
   real(R8P), allocatable                         :: P(:)                                 !< Serialized primitive variables.
   real(R8P)                                      :: C(1:2, 1-self%Ng:-1+self%Ng, 1:3)    !< Pseudo characteristic variables.
   real(R8P)                                      :: CR(1:2, 1:3)                         !< Pseudo characteristic reconst.
   integer(I4P)                                   :: i                                    !< Counter.
   integer(I4P)                                   :: j                                    !< Counter.
   integer(I4P)                                   :: f                                    !< Counter.
   integer(I4P)                                   :: v                                    !< Counter.
   class(interpolator_object), allocatable        :: interpolator                         !< WENO interpolator.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_conservative(1, i) = conservative(i)
         r_conservative(2, i) = r_conservative(1, i)
      enddo
   case(3, 5, 7, 9, 11, 13, 15, 17) ! 3rd-17th order WENO reconstruction
      call wenoof_create(interpolator_type='reconstructor-JS', S=self%Ng, interpolator=interpolator)
      do i=1-self%Ng, self%Ni+self%Ng
         primitive(i) = conservative_to_primitive_compressible(conservative=conservative(i), eos=self%eos)
      enddo
      do i=0, self%Ni+1
         do j=i+1-self%Ng, i-1+self%Ng
             P = primitive(j)%array()
            do f=1, 2
               C(f, j-i, 1) = P(1)
               C(f, j-i, 2) = P(2)
               C(f, j-i, 3) = P(5)
            enddo
         enddo
         do v=1, 3
            call interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
         enddo
         do f=1, 2
            r_primitive(f, i)%density  = CR(f, 1)
            r_primitive(f, i)%velocity = CR(f, 2) * ex
            r_primitive(f, i)%pressure = CR(f, 3)
         enddo
      enddo
      do i=0, self%Ni+1
         r_conservative(1, i) = primitive_to_conservative_compressible(primitive=r_primitive(1, i), eos=self%eos)
         r_conservative(2, i) = primitive_to_conservative_compressible(primitive=r_primitive(2, i), eos=self%eos)
      enddo
   endselect
   endsubroutine reconstruct_interfaces_primitive
endmodule foreseer_euler_1d

program foreseer_test_shock_tube
!< WenOOF test: shock tube tester, 1D Euler equation.

use flap, only : command_line_interface
use foodie, only : tvd_runge_kutta_integrator, emd_runge_kutta_integrator
use flow, only : conservative_compressible, primitive_compressible,                              &
                 conservative_to_primitive_compressible, primitive_to_conservative_compressible, &
                 eos_compressible
use foreseer_euler_1d, only : euler_1d
use penf, only : cton, FR8P, I4P, R8P, str
use vecfor, only : ex, vector

implicit none
integer(I4P)                     :: weno_order                !< WENO reconstruction order.
character(len=:), allocatable    :: weno_variables            !< Variables set on which WENO reconstruction is done.
character(len=:), allocatable    :: rk_scheme                 !< Runge-Kutta scheme type: TVD, embedded.
type(emd_runge_kutta_integrator) :: emd_rk_integrator         !< Embedded Runge-Kutta integrator.
type(tvd_runge_kutta_integrator) :: tvd_rk_integrator         !< TVD Runge-Kutta integrator.
integer(I4P)                     :: rk_stages_number          !< Runge-Kutta stages number.
type(euler_1d), allocatable      :: rk_stage(:)               !< Runge-Kutta stages.
real(R8P)                        :: dt                        !< Time step.
real(R8P)                        :: t                         !< Time.
integer(I4P)                     :: step                      !< Time steps counter.
type(euler_1d)                   :: domain                    !< Domain of Euler equations.
real(R8P)                        :: CFL                       !< CFL value.
character(3)                     :: BC_L                      !< Left boundary condition type.
character(3)                     :: BC_R                      !< Right boundary condition type.
integer(I4P)                     :: Ni                        !< Number of grid cells.
real(R8P)                        :: Dx                        !< Space step discretization.
real(R8P), allocatable           :: x(:)                      !< Cell center x-abscissa values.
integer(I4P)                     :: steps_max                 !< Maximum number of time steps.
real(R8P)                        :: t_max                     !< Maximum integration time.
character(99), allocatable       :: riemann_solver_schemes(:) !< Riemann Problem solver scheme(s).
character(99)                    :: s_scheme                  !< Space integration scheme.
character(99)                    :: t_scheme                  !< Time integration scheme.
logical                          :: results                   !< Flag for activating results saving.
logical                          :: time_serie                !< Flag for activating time serie-results saving.
logical                          :: verbose                   !< Flag for activating more verbose output.
integer(I4P)                     :: s                         !< Schemes counter.

call parse_command_line_interface
select case(rk_scheme)
case('tvd')
   do s=1, size(riemann_solver_schemes, dim=1)
      if (verbose) print "(A)", 'Use Riemann Problem solver "'//trim(adjustl(riemann_solver_schemes(s)))//'"'
      call initialize(riemann_solver_scheme=riemann_solver_schemes(s))
      call save_time_serie(filename='euler_1D-'//&
                                    trim(adjustl(s_scheme))//'-'//&
                                    trim(adjustl(t_scheme))//'-'//&
                                    trim(adjustl(riemann_solver_schemes(s)))//'.dat', t=t)
      step = 0
      t = 0._R8P
      tvd_time_loop: do
         step = step + 1
         dt = domain%dt(steps_max=steps_max, t_max=t_max, t=t, CFL=CFL)
         call tvd_rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
         t = t + dt
         call save_time_serie(t=t)
         if (verbose) print "(A)", 'step = '//str(n=step)//', time step = '//str(n=dt)//', time = '//str(n=t)
         if ((t == t_max).or.(step == steps_max)) exit tvd_time_loop
      enddo tvd_time_loop
   enddo
case('emd')
   do s=1, size(riemann_solver_schemes, dim=1)
      if (verbose) print "(A)", 'Use Riemann Problem solver "'//trim(adjustl(riemann_solver_schemes(s)))//'"'
      call initialize(riemann_solver_scheme=riemann_solver_schemes(s))
      call save_time_serie(filename='euler_1D-'//&
                                    trim(adjustl(s_scheme))//'-'//&
                                    trim(adjustl(t_scheme))//'-'//&
                                    trim(adjustl(riemann_solver_schemes(s)))//'.dat', t=t)
      step = 0
      t = 0._R8P
      emd_time_loop: do
         step = step + 1
         dt = domain%dt(steps_max=steps_max, t_max=t_max, t=t, CFL=CFL)
         call emd_rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
         t = t + dt
         call save_time_serie(t=t)
         if (verbose) print "(A)", 'step = '//str(n=step)//', time step = '//str(n=dt)//', time = '//str(n=t)
         if ((t == t_max).or.(step == steps_max)) exit emd_time_loop
      enddo emd_time_loop
   enddo
endselect

contains
   subroutine initialize(riemann_solver_scheme)
   !< Initialize the test.
   character(*), intent(in)                  :: riemann_solver_scheme !< Riemann Problem solver scheme.
   type(primitive_compressible), allocatable :: initial_state(:)      !< Initial state of primitive variables.
   integer(I4P)                              :: i                     !< Space counter.

   if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:rk_stages_number))
   select case(rk_scheme)
   case('tvd')
      call tvd_rk_integrator%init(stages=rk_stages_number)
   case('emd')
      call emd_rk_integrator%init(stages=rk_stages_number, tolerance=1.e-12_R8P)
   endselect
   t = 0._R8P
   if (allocated(x)) deallocate(x) ; allocate(x(1:Ni))
   if (allocated(initial_state)) deallocate(initial_state) ; allocate(initial_state(1:Ni))
   Dx = 1._R8P / Ni
   ! Sod's problem
   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R8P * Dx
      initial_state(i)%density  = 1._R8P
      initial_state(i)%velocity = 0._R8P
      initial_state(i)%pressure = 1._R8P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R8P * Dx
      initial_state(i)%density  = 0.125_R8P
      initial_state(i)%velocity = 0._R8P
      initial_state(i)%pressure = 0.1_R8P
   enddo
   call domain%initialize(Ni=Ni, Dx=Dx,                                         &
                          BC_L=BC_L, BC_R=BC_R,                                 &
                          initial_state=initial_state,                          &
                          eos=eos_compressible(cp=1040.004_R8P, cv=742.86_R8P), &
                          weno_order=weno_order,                                &
                          weno_variables=weno_variables,                        &
                          riemann_solver_scheme=riemann_solver_scheme)
   endsubroutine initialize

   subroutine parse_command_line_interface()
   !< Parse Command Line Interface (CLI).
   type(command_line_interface)  :: cli                   !< Command line interface handler.
   character(99)                 :: riemann_solver_scheme !< Riemann Problem solver scheme.
   integer(I4P)                  :: error                 !< Error handler.
   character(len=:), allocatable :: buffer                !< String buffer.

   call cli%init(description = 'WenOOF test: shock tube tester, 1D Euler equations', &
                 examples    = ["foreseer_test_shock_tube         ",                 &
                                "foreseer_test_shock_tube --tserie"])
   call cli%add(switch='--Ni', help='Number finite volumes used', required=.false., act='store', def='100')
   call cli%add(switch='--steps', help='Number time steps performed', required=.false., act='store', def='60')
   call cli%add(switch='--t-max', help='Maximum integration time', required=.false., act='store', def='0.')
   call cli%add(switch='--riemann', help='Riemann Problem solver', required=.false., act='store', def='all', &
                choices='all,exact,hllc,llf,pvl,roe')
   call cli%add(switch='--s-scheme', help='Space intergation scheme', required=.false., act='store', def='weno-char-1',           &
     choices='weno-char-1,weno-char-3,weno-char-5,weno-char-7,weno-char-9,weno-char-11,weno-char-13,weno-char-15,weno-char-17,'// &
             'weno-cons-1,weno-cons-3,weno-cons-5,weno-cons-7,weno-cons-9,weno-cons-11,weno-cons-13,weno-cons-15,weno-cons-17,'// &
             'weno-prim-1,weno-prim-3,weno-prim-5,weno-prim-7,weno-prim-9,weno-prim-11,weno-prim-13,weno-prim-15,weno-prim-17')
   call cli%add(switch='--t-scheme', help='Time intergation scheme', required=.false., act='store', def='tvd-rk-1', &
                choices='tvd-rk-1,tvd-rk-2,tvd-rk-3,tvd-rk-5,'// &
                        'emd-rk-2,emd-rk-6,emd-rk-7,emd-rk-9,emd-rk-17')
   call cli%add(switch='--cfl', help='CFL value', required=.false., act='store', def='0.7')
   call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-result', required=.false., act='store_true', def='.false.')
   call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.')
   call cli%parse(error=error)
   call cli%get(switch='--Ni',       val=Ni,                    error=error) ; if (error/=0) stop
   call cli%get(switch='--steps',    val=steps_max,             error=error) ; if (error/=0) stop
   call cli%get(switch='--t-max',    val=t_max,                 error=error) ; if (error/=0) stop
   call cli%get(switch='--riemann',  val=riemann_solver_scheme, error=error) ; if (error/=0) stop
   call cli%get(switch='--s-scheme', val=s_scheme,              error=error) ; if (error/=0) stop
   call cli%get(switch='--t-scheme', val=t_scheme,              error=error) ; if (error/=0) stop
   call cli%get(switch='--cfl',      val=CFL,                   error=error) ; if (error/=0) stop
   call cli%get(switch='--tserie',   val=time_serie,            error=error) ; if (error/=0) stop
   call cli%get(switch='--verbose',  val=verbose,               error=error) ; if (error/=0) stop

   if (t_max > 0._R8P) steps_max = 0

   buffer = trim(adjustl(s_scheme))
   select case(buffer(6:9))
   case('char')
      weno_variables = 'characteristic'
   case('cons')
      weno_variables = 'conservative'
   case('prim')
      weno_variables = 'primitive'
   endselect
   weno_order = cton(buffer(11:), knd=1_I4P)

   buffer = trim(adjustl(t_scheme))
   rk_scheme = buffer(1:3)
   rk_stages_number = cton(buffer(8:), knd=1_I4P)

   if (trim(adjustl(riemann_solver_scheme))=='all') then
      riemann_solver_schemes = ['exact', 'hllc ', 'llf  ', 'pvl  ', 'roe  ']
   else
      riemann_solver_schemes = [trim(adjustl(riemann_solver_scheme))]
   endif
   endsubroutine parse_command_line_interface

   subroutine save_time_serie(filename, finish, t)
   !< Save time-serie results.
   character(*), intent(in), optional :: filename  !< Output filename.
   logical,      intent(in), optional :: finish    !< Flag for triggering the file closing.
   real(R8P),    intent(in)           :: t         !< Current integration time.
   integer(I4P), save                 :: tsfile    !< File unit for saving time serie results.
   type(primitive_compressible)       :: primitive !< Primitive variables.
   integer(I4P)                       :: i         !< Counter.

   if (time_serie) then
      if (present(filename)) then
         open(newunit=tsfile, file=filename)
      endif
      write(tsfile, '(A)')'VARIABLES = "x" "rho" "u" "p"'
      write(tsfile, '(A)')'ZONE T="'//str(n=t)//'"'
      do i=1, Ni
         primitive = conservative_to_primitive_compressible(conservative=domain%U(i), eos=domain%eos)
         write(tsfile, '(4'//'('//FR8P//',1X))')x(i), primitive%density, primitive%velocity%x, primitive%pressure
      enddo
      if (present(finish)) then
         if (finish) close(tsfile)
      endif
   endif
   endsubroutine save_time_serie
endprogram foreseer_test_shock_tube
