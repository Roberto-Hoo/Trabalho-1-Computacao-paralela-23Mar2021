
! **************************************************************************
!
!   NUMERICAL SOLUTION OF  u  =  ( F(x,y) u_x )  +  ( G(x,y) u_y )
!                           t                  x                 y
!             where  F(x,y) = G(x,y) = | grad u |^(p-2)
!
! **************************************************************************
!
!  This code INITIALIZE.f90 must be run first (to initialize all variables)
!
!          The code runs from t = t0 = 0  to t = tF = 0
!
!  INPUT DATA: this code does not read from any input data
!
!  OUTPUT FILES: this code generates the following data files on disk:
!
!  (1) u24h_omp10_INPUT.dat : output file with all variables initialized
!               in order to allow the code u24h_omp10.f90 to be executed
!  (2) u24h_omp10_MATLAB.dat: the SINGLE PRECISION version (i.e.: REAL*4)
!   of data file "u24h_omp10_INPUT.dat" to be used by MATLAB for plotting
!
! ***************************************************************************

PROGRAM initialize

IMPLICIT NONE

LOGICAL, PARAMETER :: debug = .true. ! To write information in debug mode, switch to .true.
CHARACTER*1  char1
INTEGER, PARAMETER :: i64 = selected_int_kind(18)
INTEGER (i64)      :: Qt_Atualizacoes_u  ! <-- Total amount of u updates

!  COMPUTATIONAL REGION will be: [-10, 10] x [-10, 10],  with h = 0.1:

INTEGER, PARAMETER ::  M1 = -120, M2 = 120   ! <--- x-grid is: x(M1),...,x(M2)
INTEGER, PARAMETER ::  N1 = -120, N2 = 120   ! <--- y-grid is: y(N1),...,y(N2

INTEGER, PARAMETER ::  Nthreads = 16        ! <--- no. of threads in OMP parallel regions

REAL*8,  PARAMETER ::  tx_fv1 = 0.5D0  ! <--- % position of far field ZONE 1
REAL*8,  PARAMETER ::  tx_fv2 = 0.7D0  ! <--- % position of far field ZONE 2
REAL*8,  PARAMETER ::  tx_fv3 = 0.9D0  ! <--- % position of far field ZONE 3

REAL*8  refv1   ! <--- main position of far field ZONE 1, refv1 = L*tx_fv1
REAL*8  refv2   ! <--- main position of far field ZONE 2, refv2 = L*tx_fv2
REAL*8  refv3   ! <--- main position of far field ZONE 3, refv3 = L*tx_fv3

REAL*8  L   ! <--- L = min( abs(M1), abs(M2), abs(N1), abs(N2) ), lado do menor quadrado


INTEGER, PARAMETER ::  n = 2         ! <--- space dimension
REAL*8,  PARAMETER ::  p = 2.1D0     ! <--- p must be greater than n

REAL*8,  PARAMETER ::  h = 0.1D0     ! <--- spatial meshgrid spacing

REAL*8,  PARAMETER ::  cfl = 0.01D0  ! <--- cfl number (will apply to: "u24h_serial.f90")

REAL*8   x(M1:M2), y(N1:N2)          ! <--- grid point coordinates
REAL*8   u0(M1:M2,N1:N2)             ! <--- initial solution values

REAL*8,  PARAMETER ::  t0 = 0D0      ! <--- initial time
REAL*8,  PARAMETER ::  tF = 0D0      ! <--- new solution in NOT computed here

INTEGER, PARAMETER ::  no_runs = 0   ! <--- again: new solution is NOT computed here
INTEGER, PARAMETER ::  count = 0     ! <--- again: new solution is NOT computed here

REAL*8,  PARAMETER ::  dt_dump = 0.01D0   ! time between consecutive dumpings
! of solution statistics
REAL*8   x1, x2, x3, x4    ! x-coordinates of singularities
REAL*8   y1, y2, y3, y4    ! y-coordinates of singularities
INTEGER  i1, i2, i3, i4    ! i-indices of singularities
INTEGER  j1, j2, j3, j4    ! j-indices of singularities

REAL*8   b1, b2, b3, b4    ! solution values at the singularities
REAL*8   b                 ! initial guess of the far-field values
REAL*8   mass1, mass2, mass3, mass4    ! mass of initial cones relative to height b

INTEGER  i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b    ! i-indices of far-field zones
INTEGER  j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b    ! j-indices of far-field zones
INTEGER  length_i1, length_i2, length_i3                   ! i-length of far-field zones
INTEGER  length_j1, length_j2, length_j3                   ! j-length of far-field zones

REAL*8   ff1_value_u0, ff2_value_u0, ff3_value_u0          ! average values of u0 in far-field zones

REAL*8   variation_mass_0, variation_sup_0

REAL*8   min_u0, max_u0    ! minimum and maximum values of initial solution over the grid
REAL*8   mass_u0           ! mass of initial solution u0 relative to reference height b

! Extra variables needed to compute
! variation_mass_0, variation_sup_0:
! (note that  cfl  has ALREADY been defined above)

REAL*8   u(M1-1:M2+1,N1-1:N2+1)
REAL*8   v(M1-1:M2+1,N1-1:N2+1)
REAL*8   F(M1:M2+1,N1:N2)
REAL*8   G(M1:M2,N1:N2+1)

REAL*8   t
REAL*8   dt
REAL*8   time_next_dump

REAL*8   q
REAL*8   h_to_pm2

REAL*8   x_min, x_max
REAL*8   y_min, y_max

! Local values of some additional mesh parameters:

REAL*8   minimum_distance_to_travel_1
REAL*8   minimum_distance_to_travel_2
REAL*8   minimum_distance_to_travel_3
REAL*8   minimum_distance_to_travel_4

REAL*8   maximum_distance_to_travel_1
REAL*8   maximum_distance_to_travel_2
REAL*8   maximum_distance_to_travel_3
REAL*8   maximum_distance_to_travel_4

REAL*8   t1_min, t2_min, t3_min, t4_min
REAL*8   t1_max, t2_max, t3_max, t4_max

REAL*8   d1, d2, d3, d4

REAL*8   r12, r13, r14, r23, r24, r34
REAL*8   r0

INTEGER  DATE_start (8)
INTEGER  DATE_finish(8)

REAL*8   time_start, time_finish
REAL*8   total_elapsed_time

! Local variables:

INTEGER  i, j

REAL*8,  PARAMETER ::  pi = 3.141592653589793238D0;

REAL*8   expoente, common_factor
REAL*8   C1, C2, C3, C4
REAL*8   cc1, cc2, cc3, cc4
REAL*8   Rx, Ry

REAL*8   expected_b

! Intrinsic functions used:

REAL*8, INTRINSIC :: sum
REAL*8, INTRINSIC :: abs, sqrt
REAL*8, INTRINSIC :: min, max
REAL*8, INTRINSIC :: minval, maxval
REAL*8, INTRINSIC :: gamma
REAL*4, INTRINSIC :: real

INTRINSIC  DATE_AND_TIME
INTRINSIC  CPU_TIME
!  INTRINSIC  DBLE


!******************************************************************
!*                                                                *
!*                  MAIN PROGRAM starts here:                     *
!*                                                                *
!******************************************************************


CALL DATE_AND_TIME(VALUES = DATE_start)

CALL CPU_TIME(time_start)

x_min = M1*h; x_max = M2*h;
y_min = N1*h; y_max = N2*h;

L  = min( abs(x_min), abs(x_max), abs(y_min), abs(y_max))

refv1 = L*tx_fv1 !<--- main position of far field ZONE 1
refv2 = L*tx_fv2 !<--- main position of far field ZONE 2
refv3 = L*tx_fv3 !<--- main position of far field ZONE 3


if ( max(refv1, refv2, refv3) > L ) then
	write(6,*) '*******************************************************'
	write(6,*) '****** ERROR ***** ERROR ***** ERROR ***** ERROR ******'
	write(6,*)
	write(6,*) 'Far-field zones INCONSISTENT with computational region'
	write(6,*)
	write(6,*) 'Please reset the reference values'
	write(6,*) '                                   refv1, refv2, refv3'
	write(6,*) 'which define the far-field zones'
	write(6,*) 'in code: "initialize.f90"'
	write(6,*)
	write(6,*) '*******************************************************'
	write(6,*)
	write(6,*) 'Execution will be ABORTED because of this error'
	write(6,*)
	write(6,*) '*******************************************************'

	STOP

endif

!--------------------------------------------------------------
! Location of singularities and corresponding solution values:
!--------------------------------------------------------------

x1 = +1.0D0;   ! <--- x-coordinate of 1st singularity
y1 = +0.5D0;   ! <--- y-coordinate of 1st singularity

x2 = -1.5D0;   ! <--- x-coordinate of 2nd singularity
y2 = +1.0D0;   ! <--- y-coordinate of 2nd singularity

x3 = +0.5D0;   ! <--- x-coordinate of 3rd singularity
y3 = -1.5D0;   ! <--- y-coordinate of 3rd singularity

x4 = -1.5D0;   ! <--- x-coordinate of 4rd singularity
y4 = -1.5D0;   ! <--- y-coordinate of 4rd singularity


b1 = 5D0;

b2 = 2D0;

b3 = 1D0;

b4 = -2D0;

b = (b1 + b2 + b3 +b4)/4 +0.5;

! expected far-field solution value is: (b1 + b2 + b3 + b4)/4
expected_b = (b1 + b2 + b3 +b4)/4;
!--------------------------------------------------------------


i1 = M1; i2 = M1; i3 = M1; i4 = M1;
j1 = N1; j2 = N1; j3 = N1; j4 = N1;

DO i = M1, M2
	x(i) = i*h;
	IF ( x(i) < x1 + h/100D0 ) i1 = i;    ! <--- i-index of 1st singularity
	IF ( x(i) < x2 + h/100D0 ) i2 = i;    ! <--- i-index of 2nd singularity
	IF ( x(i) < x3 + h/100D0 ) i3 = i;    ! <--- i-index of 3rd singularity
	IF ( x(i) < x4 + h/100D0 ) i4 = i;    ! <--- i-index of 4th singularity
ENDDO

DO j = N1, N2
	y(j) = j*h;
	IF ( y(j) < y1 + h/100D0 ) j1 = j;    ! <--- j-index of 1st singularity
	IF ( y(j) < y2 + h/100D0 ) j2 = j;    ! <--- j-index of 2nd singularity
	IF ( y(j) < y3 + h/100D0 ) j3 = j;    ! <--- j-index of 3rd singularity
	IF ( y(j) < y4 + h/100D0 ) j4 = j;    ! <--- j-index of 4th singularity
ENDDO

!--------------------------------------------------
!
!             INITIAL APPROXIMATION
!                     to
!             STEADY STATE SOLUTION
!
!--------------------------------------------------

r12 = sqrt( (x1 - x2)**2 + (y1 - y2)**2 );
r13 = sqrt( (x1 - x3)**2 + (y1 - y3)**2 );
r14 = sqrt( (x1 - x4)**2 + (y1 - y4)**2 );
r23 = sqrt( (x2 - x3)**2 + (y2 - y3)**2 );
r24 = sqrt( (x2 - x4)**2 + (y2 - y4)**2 );
r34 = sqrt( (x3 - x4)**2 + (y3 - y4)**2 );

r0 = 1D0/2D0 * min(r12, r13, r14, r23, r24, r34);

!--------------------------------------------------

! Computing the mass of each initial cone
! relatively to the base plane z = b:

mass1 = 1D0/3D0 * pi*r0**2 * abs( b1 - b );
mass2 = 1D0/3D0 * pi*r0**2 * abs( b2 - b );
mass3 = 1D0/3D0 * pi*r0**2 * abs( b3 - b );
mass4 = 1D0/3D0 * pi*r0**2 * abs( b4 - b );

! The radius of the support of each initial cone
! for large t is asymptotically given by
! Rj = ccj * t^(1/(n*(p-2)+p))
! where ccj is given below (using Barenblatt solution):

expoente = p/(p + n*(p-2)) * (p-2)/(p-1);
common_factor = ( 1/(n*(p-2)+p)**(n/p) * p/(p-1) / (2*pi) * &
	((p-2)/p)**(n*(p-1)/p) / beta( n*(p-1)/p, (2*p-3)/(p-2) ) )**expoente;

C1 = common_factor * abs(mass1)**expoente; if (C1 < 1.0D-3) C1 = 1.0D-3;
C2 = common_factor * abs(mass2)**expoente; if (C2 < 1.0D-3) C2 = 1.0D-3;
C3 = common_factor * abs(mass3)**expoente; if (C3 < 1.0D-3) C3 = 1.0D-3;
C4 = common_factor * abs(mass4)**expoente; if (C4 < 1.0D-3) C4 = 1.0D-3;

cc1 = 1/(n*(p-2)+p)**(1/p) * ((p-2)/p)**((p-1)/p) / C1**((p-1)/p);
cc2 = 1/(n*(p-2)+p)**(1/p) * ((p-2)/p)**((p-1)/p) / C2**((p-1)/p);
cc3 = 1/(n*(p-2)+p)**(1/p) * ((p-2)/p)**((p-1)/p) / C3**((p-1)/p);
cc4 = 1/(n*(p-2)+p)**(1/p) * ((p-2)/p)**((p-1)/p) / C4**((p-1)/p);

! Estimating time need to reach the boundary of
! computational region [-Rx, Rx] x [-Ry, Ry]:

minimum_distance_to_travel_1 = min(x_max-x1, x1-x_min, y_max-y1, y1-y_min);
minimum_distance_to_travel_2 = min(x_max-x2, x2-x_min, y_max-y2, y2-y_min);
minimum_distance_to_travel_3 = min(x_max-x3, x3-x_min, y_max-y3, y3-y_min);
minimum_distance_to_travel_4 = min(x_max-x4, x4-x_min, y_max-y4, y4-y_min);

maximum_distance_to_travel_1 = max(x_max-x1, x1-x_min, y_max-y1, y1-y_min);
maximum_distance_to_travel_2 = max(x_max-x2, x2-x_min, y_max-y2, y2-y_min);
maximum_distance_to_travel_3 = max(x_max-x3, x3-x_min, y_max-y3, y3-y_min);
maximum_distance_to_travel_4 = max(x_max-x4, x4-x_min, y_max-y4, y4-y_min);

t1_min = ( cc1 * minimum_distance_to_travel_1 )**(n*(p-2)+p);
t2_min = ( cc2 * minimum_distance_to_travel_2 )**(n*(p-2)+p);
t3_min = ( cc3 * minimum_distance_to_travel_3 )**(n*(p-2)+p);
t4_min = ( cc4 * minimum_distance_to_travel_4 )**(n*(p-2)+p);

t1_max = ( cc1 * maximum_distance_to_travel_1 )**(n*(p-2)+p);
t2_max = ( cc2 * maximum_distance_to_travel_2 )**(n*(p-2)+p);
t3_max = ( cc3 * maximum_distance_to_travel_3 )**(n*(p-2)+p);
t4_max = ( cc4 * maximum_distance_to_travel_4 )**(n*(p-2)+p);

write(6,*) '-----------------------------------------------------------'
write(6,*) 'Estimated minimum time to reach the computational boundary:'
write(6,100) t1_min, t2_min, t3_min, t4_min

100 FORMAT( 1X, E15.6 )

write(6,*) 'Estimated time to pass the computational boundary entirely:'
write(6,100) t1_max, t2_max, t3_max, t4_max
write(6,*) '-----------------------------------------------------------'

!--------------------------------------------------
!  Setting the values of initial solution u0:
!--------------------------------------------------

u0 = b;

DO i = M1, M2
	DO j = N1, N2

		d1 = sqrt( (x(i) - x1)**2 + (y(j) - y1)**2 );

		d2 = sqrt( (x(i) - x2)**2 + (y(j) - y2)**2 );

		d3 = sqrt( (x(i) - x3)**2 + (y(j) - y3)**2 );

		d4 = sqrt( (x(i) - x4)**2 + (y(j) - y4)**2 );

		if ( d1 < r0 ) then
			u0(i,j) = b1 - (b1 - b)*d1/r0;
		endif

		if ( d2 < r0 ) then
			u0(i,j) = b2 - (b2 - b)*d2/r0;
		endif

		if ( d3 < r0 ) then
			u0(i,j) = b3 - (b3 - b)*d3/r0;
		endif

		if ( d4 < r0 ) then
			u0(i,j) = b4 - (b4 - b)*d4/r0;
		endif

	ENDDO

ENDDO

u0(i1,j1) = b1;
u0(i2,j2) = b2;
u0(i3,j3) = b3;
u0(i4,j4) = b4;

!--------------------------------------------------
!  Finding the indices of far-field zones:
!--------------------------------------------------

i_ff1a = M1; i_ff1b = M1;
i_ff2a = M1; i_ff2b = M1;
i_ff3a = M1; i_ff3b = M1;

DO i = M1, M2
	if ( x(i) < -refv1 + h/100D0 ) i_ff1a = i;
	if ( x(i) < -refv2 + h/100D0 ) i_ff2a = i;
	if ( x(i) < -refv3 + h/100D0 ) i_ff3a = i;
	if ( x(i) < +refv1 + h/100D0 ) i_ff1b = i;
	if ( x(i) < +refv2 + h/100D0 ) i_ff2b = i;
	if ( x(i) < +refv3 + h/100D0 ) i_ff3b = i;
ENDDO

j_ff1a = N1; j_ff1b = N1;
j_ff2a = N1; j_ff2b = N1;
j_ff3a = N1; j_ff3b = N1;

DO j = N1, N2
	if ( y(j) < -refv1 + h/100D0 ) j_ff1a = j;
	if ( y(j) < -refv2 + h/100D0 ) j_ff2a = j;
	if ( y(j) < -refv3 + h/100D0 ) j_ff3a = j;
	if ( y(j) < +refv1 + h/100D0 ) j_ff1b = j;
	if ( y(j) < +refv2 + h/100D0 ) j_ff2b = j;
	if ( y(j) < +refv3 + h/100D0 ) j_ff3b = j;
ENDDO

length_i1 = i_ff1b - i_ff1a + 1;
length_i2 = i_ff2b - i_ff2a + 1;
length_i3 = i_ff3b - i_ff3a + 1;

length_j1 = j_ff1b - j_ff1a + 1;
length_j2 = j_ff2b - j_ff2a + 1;
length_j3 = j_ff3b - j_ff3a + 1;


!-----------------------------------------------------
!  Computing initial values of u0 on far-field zones:
!-----------------------------------------------------

ff1_value_u0  = ( sum( u0(i_ff1a:i_ff1a+10,j_ff1a:j_ff1b) )/(11*length_j1) &
	+ sum( u0(i_ff1b-10:i_ff1b,j_ff1a:j_ff1b) )/(11*length_j1) &
	+ sum( u0(i_ff1a:i_ff1b,j_ff1a:j_ff1a+10) )/(11*length_i1) &
	+ sum( u0(i_ff1a:i_ff1b,j_ff1b-10:j_ff1b) )/(11*length_i1) )/4D0;

ff2_value_u0  = ( sum( u0(i_ff2a:i_ff2a+10,j_ff2a:j_ff2b) )/(11*length_j2) &
	+ sum( u0(i_ff2b-10:i_ff2b,j_ff2a:j_ff2b) )/(11*length_j2) &
	+ sum( u0(i_ff2a:i_ff2b,j_ff2a:j_ff2a+10) )/(11*length_i2) &
	+ sum( u0(i_ff2a:i_ff2b,j_ff2b-10:j_ff2b) )/(11*length_i2) )/4D0;

ff3_value_u0  = ( sum( u0(i_ff3a:i_ff3a+10,j_ff3a:j_ff3b) )/(11*length_j3) &
	+ sum( u0(i_ff3b-10:i_ff3b,j_ff3a:j_ff3b) )/(11*length_j3) &
	+ sum( u0(i_ff3a:i_ff3b,j_ff3a:j_ff3a+10) )/(11*length_i3) &
	+ sum( u0(i_ff3a:i_ff3b,j_ff3b-10:j_ff3b) )/(11*length_i3) )/4D0;

!--------------------------------------------------

min_u0 = minval( u0 );
max_u0 = maxval( u0 );

mass_u0 = sum( u0 - b )*h**2;


!******************************************************************
!*                                                                *
!*     The following computation determines a starting value      *
!*     for the quantities  variation_mass, variation_sup:         *
!*                                                                *
!******************************************************************

!----------------------------------------------------------
!      Extending solution values to extended grid:
!----------------------------------------------------------

u(M1:M2,N1:N2) = u0;

u(M1-1,N1:N2) = u(M1,N1:N2);
u(M2+1,N1:N2) = u(M2,N1:N2);
u(M1-1:M2+1,N1-1) = u(M1-1:M2+1,N1);
u(M1-1:M2+1,N2+1) = u(M1-1:M2+1,N2);

!----------------------------------------------------------
!----------------------------------------------------------

write(6,*) '-----------------------------------------------------------'
write(6,*) ' este eh o NOVO programa initialize.f90'
write(6,*) '-----------------------------------------------------------'

write(6,*) ' Computing from t = 0 to t = 0.01... '

q = (p - 2D0)/2D0;

h_to_pm2 = h**(p-2);

time_next_dump = t0 + dt_dump;

dt = cfl*h*h;

t = 0;
Qt_Atualizacoes_u = 0;

DO WHILE ( t < time_next_dump - dt/4 )

	t = t + dt;


	! ****************************************************************** !
	!                computation of F values on the grid:                !
	! ****************************************************************** !
	! ******  F(i,j) for M1 <= i <= M2 + 1 and N1 <= j <= N2:

	F(M1:M2+1,N1:N2) = ( (u(M1:M2+1,N1:N2) - u(M1-1:M2,N1:N2))**2 + &
		( (u(M1-1:M2,N1+1:N2+1)+u(M1:M2+1,N1+1:N2+1))-(u(M1-1:M2,N1-1:N2-1)+u(M1:M2+1,N1-1:N2-1)) )**2 /16 )**q / h_to_pm2;

	! ****************************************************************** !
	!                computation of G values on the grid:                !
	! ****************************************************************** !
	! ******  G(i,j) for M1 <= i <= M2 and N1 <= j <= N2 + 1:

	G(M1:M2,N1:N2+1) = ( (u(M1:M2,N1:N2+1) - u(M1:M2,N1-1:N2))**2 + &
		( (u(M1+1:M2+1,N1-1:N2)+u(M1+1:M2+1,N1:N2+1))-(u(M1-1:M2-1,N1-1:N2)+u(M1-1:M2-1,N1:N2+1)) )**2 /16 )**q / h_to_pm2;


	! ****************************************************************** !
	!      computation of v = [ new u values at the new time level ]:    !
	! ****************************************************************** !
	! ******  v(i,j) for M1 <= i <= M2 and N1 <= j <= N2:

	v(M1:M2,N1:N2) = u(M1:M2,N1:N2) + &
		cfl*( F(M1+1:M2+1,N1:N2)*(u(M1+1:M2+1,N1:N2)-u(M1:M2,N1:N2)) - F(M1:M2,N1:N2)*(u(M1:M2,N1:N2)-u(M1-1:M2-1,N1:N2)) ) + &
		cfl*( G(M1:M2,N1+1:N2+1)*(u(M1:M2,N1+1:N2+1)-u(M1:M2,N1:N2)) - G(M1:M2,N1:N2)*(u(M1:M2,N1:N2)-u(M1:M2,N1-1:N2-1)) );

	! ******  extending  v  to the extra grid points:

	v(M1-1,N1:N2) = v(M1,N1:N2);
	v(M2+1,N1:N2) = v(M2,N1:N2);
	v(M1-1:M2+1,N1-1) = v(M1-1:M2+1,N1);
	v(M1-1:M2+1,N2+1) = v(M1-1:M2+1,N2);

	! ****************************************************************** !
	!            computation of new u values completed!                  !
	! ****************************************************************** !

	u = v;    ! <--- updating u (at the new time level)

	! ******  correcting u values at the singular points:

	u(i1,j1) = b1;

	u(i2,j2) = b2;

	u(i3,j3) = b3;

	u(i4,j4) = b4;

	Qt_Atualizacoes_u = Qt_Atualizacoes_u +1;

ENDDO

!******************************************************************
!*                                                                *
!*               Computation of u COMPLETED                       *
!*                                                                *
!******************************************************************

!--- Computing starting values of  variation_mass, variation_sup:

variation_mass_0 = sum( u(M1:M2,N1:N2) - u0 )*h**2;

variation_sup_0 = maxval( abs( u(M1:M2,N1:N2) - u0 ) );

!--- Final computations:

CALL DATE_AND_TIME(VALUES = DATE_finish)

!  time_finish = dble(DATE_finish(7)) + dble(DATE_finish(8))/1000D0;   ! <--- not appropriate
CALL CPU_TIME(time_finish)

total_elapsed_time = time_finish - time_start;   ! <--- unit: seconds

write(6,*) ' Done! '
write(6,*) '-----------------------------------------------------------'
write(6,*) ' Saving data to disk files  "u24h_omp10_INPUT.dat" '
write(6,*) ' and "u24h_omp10_MATLAB.dat"(for visualization)... '

!--------------------------------------------------

!--------------------------------------------------
!  Saving important data to output file
! (to be used by code u24h_omp10.f90):
!--------------------------------------------------

OPEN(10, FILE = 'u24h_omp10_INPUT.dat', ACTION = 'WRITE')

510 FORMAT( 10(1X,I10) )
520 FORMAT( 10(1X,E25.16) )
525 FORMAT( 10(1X,I12) )

write(10,510) M1, M2
write(10,510) N1, N2
write(10,510) Nthreads

write(10,510) n

write(10,520) p
write(10,520) h

write(10,520) cfl

write(10,520) x_min, x_max
write(10,520) y_min, y_max

write(10,520) t0, tF, dt_dump

write(10,520) x1, x2, x3, x4
write(10,520) y1, y2, y3, y4
write(10,510) i1, i2, i3, i4
write(10,510) j1, j2, j3, j4

write(10,520) b1, b2, b3, b4
write(10,520) b

write(10,520) refv1, refv2, refv3

write(10,510) count
write(10,510) no_runs

write(10,510) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
write(10,510) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
write(10,510) length_i1, length_i2, length_i3
write(10,510) length_j1, length_j2, length_j3

write(10,520) ff1_value_u0, ff2_value_u0, ff3_value_u0

write(10,520) min_u0, max_u0
write(10,520) mass_u0

write(10,520) t1_min, t2_min, t3_min, t4_min
write(10,520) t1_max, t2_max, t3_max, t4_max

write(10,520)

! Saving LOG variable values:

write(10,510) (DATE_start (i), i = 1,8)
write(10,510) (DATE_finish(i), i = 1,8)

write(10,520) tF

write(10,520) ff1_value_u0
write(10,520) ff2_value_u0
write(10,520) ff3_value_u0

write(10,520) variation_mass_0
write(10,520) variation_sup_0

write(10,520) min_u0
write(10,520) max_u0

write(10,520) total_elapsed_time

write(10,520)

write(10,520) (x(i), i = M1, M2)
write(10,520)
write(10,520) (y(j), j = N1, N2)

write(10,520)

DO i = M1, M2
	write(10,520) (u0(i,j), j = N1, N2)
ENDDO

! Saving SOLUTION STATISTICS variables:

write(10,510)
write(10,510) count
write(10,510)

write(10,520) ff1_value_u0, ff2_value_u0, ff3_value_u0

write(10,520)
write(10,520) variation_mass_0, variation_sup_0

write(10,520)
write(10,520) min_u0, max_u0

write(10,520)
write(10,520) total_elapsed_time

write(10,520)
write(10,525) Qt_Atualizacoes_u

CLOSE(10)

!--------------------------------------------------
!  Saving data to MATLAB file for plotting:
!--------------------------------------------------

OPEN(20, FILE = 'u24h_omp10_MATLAB.dat', ACTION = 'WRITE')

530 FORMAT( 10(1X,I10) )
535 FORMAT( 10(1X,I12) )
540 FORMAT( 10(1X,E20.8) )

write(20,530) M1, M2
write(20,530) N1, N2
write(20,530) Nthreads

write(20,530) n

write(20,540) real(p)
write(20,540) real(h)

write(20,540) real(cfl)

write(20,540) real(x_min), real(x_max)
write(20,540) real(y_min), real(y_max)

write(20,540) real(t0), real(tF), real(dt_dump)

write(20,540) real(x1), real(x2), real(x3), real(x4)
write(20,540) real(y1), real(y2), real(y3), real(y4)
write(20,530) i1, i2, i3, i4
write(20,530) j1, j2, j3, j4

write(20,540) real(b1), real(b2), real(b3), real(b4)
write(20,540) real(b)

write(20,540) real(refv1), real(refv2), real(refv3)

write(20,530) count
write(20,530) no_runs

write(20,530) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
write(20,530) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
write(20,530) length_i1, length_i2, length_i3
write(20,530) length_j1, length_j2, length_j3

write(20,540) real(ff1_value_u0), real(ff2_value_u0), real(ff3_value_u0)

write(20,540) real(min_u0), real(max_u0)
write(20,540) real(mass_u0)

write(20,540) real(t1_min), real(t2_min), real(t3_min), real(t4_min)
write(20,540) real(t1_max), real(t2_max), real(t3_max), real(t4_max)

write(20,540)

! Saving LOG variable values:

write(20,530) (DATE_start (i), i = 1,8)
write(20,530) (DATE_finish(i), i = 1,8)

write(20,540) real(tF)

write(20,540) real(ff1_value_u0)
write(20,540) real(ff2_value_u0)
write(20,540) real(ff3_value_u0)

write(20,540) real(variation_mass_0)
write(20,540) real(variation_sup_0 )

write(20,540) real(min_u0)
write(20,540) real(max_u0)

write(20,540) real(total_elapsed_time)

write(20,540)

write(20,540) (real(x(i)), i = M1, M2)
write(20,540)
write(20,540) (real(y(j)), j = N1, N2)

write(20,540)

DO i = M1, M2
	write(20,540) (real(u0(i,j)), j = N1, N2)
ENDDO

! Saving SOLUTION STATISTICS variables:

write(20,530)
write(20,530) count
write(20,530)

write(20,540) real(ff1_value_u0), real(ff2_value_u0), real(ff3_value_u0)

write(20,540)
write(20,540) real(variation_mass_0), real(variation_sup_0)

write(20,540)
write(20,540) real(min_u0), real(max_u0)

write(20,540)
write(20,540) real(total_elapsed_time)

write(20,530)
write(20,535) Qt_Atualizacoes_u

CLOSE(20)

WRITE(6,*) ' Done! '
write(6,*) '-----------------------------------------------------------'
write(6,*) '***** SUCCESSFUL EXECUTION ***** SUCCESSFUL EXECUTION *****'
write(6,*) '-----------------------------------------------------------'

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

CONTAINS

FUNCTION beta(x,y)

real*8  x, y, beta

beta = gamma(x)*gamma(y)/gamma(x+y);

END FUNCTION beta

END PROGRAM initialize

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
