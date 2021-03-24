!
! **************************************************************************
!
!   NUMERICAL SOLUTION OF  u  =  ( F(x,y) u_x )  +  ( G(x,y) u_y )
!                           t                  x                 y
!             where  F(x,y) = G(x,y) = | grad u |^(p-2)
!
! **************************************************************************
!
!  Code INITIALIZE.f90 must have already been run !!!!! <---IMPORTANT
!
!          The code runs from t = t0  to t = tF
!
!  INPUT DATA: this code reads data from file "u24h_serial_input.dat" (*)
!
!  OUTPUT FILES: this code generates the following data files on disk:
!
!  (1) u24h_serial_input_previous.dat : a safe copy of the latest input file
!                                    used by this code (i.e., file (*) above)
!  (2) u24h_serial_INPUT.dat : output file with all major quantities computed
!                  by this code PLUS all previous input data contained in (*)
!  (3) u24h_serial_LOG.txt : a log file with running statistics and a summary
!              of major output results computed by this and previous runnings
!  (4) u24h_serial_ERROR.txt : an ERROR log file explaining some instances of
!           program error (inconsistencies or lack of proper input data file)
!
! ***************************************************************************
!

PROGRAM u24h_serial4P   ! <--- to be executed ONLY after program
!       "initialize.f90" has been run
IMPLICIT NONE

LOGICAL, PARAMETER :: debug = .false. ! To write information in debug mode, switch to .true.
CHARACTER*1  char1
INTEGER, PARAMETER :: i64 = selected_int_kind(18) !integer de 64bits, 8bytes
INTEGER (i64)      :: Qt_Atualizacoes_u  ! <--- Total amount of u updates


INTEGER  M1, M2   ! horizontal grid is: x(M1), ..., x(M2)
INTEGER  N1, N2   ! vertical grid is:   y(N1), ..., y(N2)

INTEGER  duracao !<--- TempoAtual-TempoAnterior(em segundos)

INTEGER  n     ! <--- n is the space dimension
REAL*8   p     ! <--- p must be greater than n

REAL*8   h     ! <--- spatial grid spacing

REAL*8   cfl   ! <--- Courant-Friedrichs-Lewy number (defined by code: "initialize.f90")

REAL*8   x_min, x_max   ! <--- computational region is:
REAL*8   y_min, y_max   ! [ x_min, x_max ] x [ y_min, y_max ]

REAL*8   dt    ! <--- timestep length

REAL*8   refv1, refv2, refv3   ! <--- reference values defining far-field zones

REAL*8,  DIMENSION(:),   ALLOCATABLE ::  x, y         ! <--- grid points
REAL*8,  DIMENSION(:,:), ALLOCATABLE ::  u            ! <--- extended solution
REAL*8,  DIMENSION(:,:), ALLOCATABLE ::  v            ! <--- (temporary) extended solution at new time level
REAL*8,  DIMENSION(:,:), ALLOCATABLE ::  u_previous   ! <--- solution at the previous data dumping

REAL*8,  DIMENSION(:,:), ALLOCATABLE ::  F
REAL*8,  DIMENSION(:,:), ALLOCATABLE ::  G

REAL*8   t0, tF
REAL*8   t

REAL*8   dt_dump           ! <--- specified by the program: initialize.f90
REAL*8   time_next_dump    !           (usually defined to be: dt_dump = 0.01)

REAL*8   x1, x2, x3, x4
REAL*8   y1, y2, y3, y4
INTEGER  i1, i2, i3, i4
INTEGER  j1, j2, j3, j4

REAL*8   b1, b2, b3, b4
REAL*8   b

INTEGER  i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
INTEGER  j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
INTEGER  length_i1, length_i2, length_i3
INTEGER  length_j1, length_j2, length_j3

REAL*8   t1_min, t2_min, t3_min, t4_min
REAL*8   t1_max, t2_max, t3_max, t4_max

REAL*8   ff1_value_u0, ff2_value_u0, ff3_value_u0   ! <--- far field values of initial data

REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff1_value_u,  ff2_value_u,  ff3_value_u
REAL*8,  DIMENSION(:), ALLOCATABLE ::  min_u, max_u
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_mass, variation_sup
REAL*8,  DIMENSION(:), ALLOCATABLE ::  time_per_cycle

REAL*8   total_elapsed_time

INTEGER  count

REAL*8   min_u0, max_u0
REAL*8   mass_u0         ! <--- not used

INTEGER  no_runs

! PREVIOUS data (already computed) to be retained:

REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff1_value_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff2_value_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff3_value_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  min_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  max_u_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_mass_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_sup_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  time_per_cycle_previous

INTEGER, DIMENSION(:,:), ALLOCATABLE :: DATE_start_LOG_previous
INTEGER, DIMENSION(:,:), ALLOCATABLE :: DATE_finish_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff1_value_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff2_value_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  ff3_value_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_mass_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  variation_sup_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  min_u_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  max_u_LOG_previous

REAL*8,  DIMENSION(:), ALLOCATABLE ::  tF_LOG_previous
REAL*8,  DIMENSION(:), ALLOCATABLE ::  elapsed_time_LOG_previous

REAL*8   variation_mass_0
REAL*8   variation_sup_0

! Local variables:

INTEGER  i, j, numLoops

REAL*8   h_to_pm2, q

LOGICAL  file_present

INTEGER  DATE_start(8)
INTEGER  DATE_finish(8)

REAL*8   time_start, time_finish   ! <--- used to measure elapsed time
!        of each iteration cycle
REAL*8   tF_new, tF_check

INTEGER  previous_length
INTEGER  new_length

INTEGER  no_dumps

INTEGER  new_count

REAL*8   FP_count, FLOPs

! Intrinsic functions used:

REAL*8,  INTRINSIC :: sum
REAL*8,  INTRINSIC :: abs
!  REAL*8,  INTRINSIC :: min, max
REAL*8,  INTRINSIC :: minval, maxval
INTEGER, INTRINSIC :: size, nint    ! <--- NINT rounds REALs to nearest INTEGERs

INTRINSIC DATE_AND_TIME
INTRINSIC CPU_TIME

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!*****************************************************************************
!*****************************************************************************
!**              EXECUTION of  u24h_serial.f90  STARTS HERE:                **
!*****************************************************************************
!*****************************************************************************
do numLoops = 1, 2000

	CALL DATE_AND_TIME(VALUES = DATE_start)

	!--------------------------------------------------------------------------
	!
	!  READING INPUT DATA FILE  u24h_serial_INPUT.dat:
	!   (created by program previous run of "u24h_serial.f90"
	!    or by program "initialize.f90" on first run)
	!
	!--------------------------------------------------------------------------

	INQUIRE( FILE = 'u24h_serial_INPUT.dat', EXIST = file_present )

	IF ( file_present .EQV. .FALSE. ) THEN

		OPEN(11, FILE = 'u24h_serial_ERROR.txt', ACTION = 'WRITE', &
			STATUS = 'REPLACE')
		write(11,*)
		write(11,*) '*******************************************************'
		write(11,*) ' ERROR ***** ERROR ***** ERROR ***** ERROR ***** ERROR '
		write(11,*) '*******************************************************'
		write(11,*)
		write(11,51) DATE_start(3), DATE_start(2), DATE_start(1)
		51     FORMAT(1X,'DAY :  ',I2.2,' / ',I2.2,' / ',I4)
		write(11,52) DATE_start(5:7)
		52     FORMAT(1X,'TIME:  ',I2.2,' : ',I2.2,' : ',I2.2)
		write(11,*)
		write(11,*) '*******************************************************'
		write(11,*) ' ERROR CONDITION:'
		write(11,*) ' FILE  u24h_serial_INPUT.dat  NOT FOUND'
		write(11,*)
		write(11,*) ' EXECUTION WILL TERMINATE'
		write(11,*)
		write(11,*) '*******************************************************'

		CLOSE(11)

		write(6,*) '*****************************************'
		write(6,*) ' ERROR: see file "u24h_serial_ERROR.txt" '
		write(6,*) '        for explanation                  '
		write(6,*) ' ***** EXECUTION HAS BEEN ABORTED *****  '
		write(6,*) '*****************************************'

		STOP

	ENDIF

	OPEN(10, FILE = 'u24h_serial_INPUT.dat', ACTION = 'READ')

	510 FORMAT( 10(1X,I10) )
	515 FORMAT( 10(1X,I12) )
	520 FORMAT( 10(1X,E25.16))

	read(10,510) M1, M2
	read(10,510) N1, N2

	read(10,510) n

	ALLOCATE( x(M1:M2), y(N1:N2) )                ! <--- spatial meshgrid (grid points)
	ALLOCATE( u(M1-1:M2+1,N1-1:N2+1) )            ! <--- extended solution
	ALLOCATE( v(M1-1:M2+1,N1-1:N2+1) )            ! <--- (temporary) solution at new time level
	ALLOCATE( u_previous(M1-1:M2+1,N1-1:N2+1) )   ! <--- solution at the previous statistics dumping

	ALLOCATE( F(M1:M2+1,N1:N2) )
	ALLOCATE( G(M1:M2,N1:N2+1) )

	read(10,520) p
	read(10,520) h

	read(10,520) cfl

	read(10,520) x_min, x_max
	read(10,520) y_min, y_max

	read(10,520) t0, tF, dt_dump

	read(10,520) x1, x2, x3, x4
	read(10,520) y1, y2, y3, y4
	read(10,510) i1, i2, i3, i4
	read(10,510) j1, j2, j3, j4

	read(10,520) b1, b2, b3, b4
	read(10,520) b

	read(10,520) refv1, refv2, refv3

	read(10,510) count
	read(10,510) no_runs

	read(10,510) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
	read(10,510) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
	read(10,510) length_i1, length_i2, length_i3
	read(10,510) length_j1, length_j2, length_j3

	read(10,520) ff1_value_u0, ff2_value_u0, ff3_value_u0

	read(10,520) min_u0, max_u0
	read(10,520) mass_u0

	read(10,520) t1_min, t2_min, t3_min, t4_min
	read(10,520) t1_max, t2_max, t3_max, t4_max


	read(10,520)

	!  Reallocating and READING previous LOG history variables:

	ALLOCATE( DATE_start_LOG_previous (0:no_runs,1:8) )
	ALLOCATE( DATE_finish_LOG_previous(0:no_runs,1:8) )
	ALLOCATE( tF_LOG_previous(0:no_runs) )
	ALLOCATE( ff1_value_LOG_previous(0:no_runs) )
	ALLOCATE( ff2_value_LOG_previous(0:no_runs) )
	ALLOCATE( ff3_value_LOG_previous(0:no_runs) )
	ALLOCATE( variation_mass_LOG_previous(0:no_runs) )
	ALLOCATE( variation_sup_LOG_previous (0:no_runs) )
	ALLOCATE( min_u_LOG_previous(0:no_runs) )
	ALLOCATE( max_u_LOG_previous(0:no_runs) )
	ALLOCATE( elapsed_time_LOG_previous(0:no_runs) )

	DO i = 0, no_runs
		read(10,510) (DATE_start_LOG_previous (i,j), j = 1,8)
	ENDDO
	DO i = 0, no_runs
		read(10,510) (DATE_finish_LOG_previous(i,j), j = 1,8)
	ENDDO

	DO i = 0, no_runs
		read(10,520) tF_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		read(10,520) ff1_value_LOG_previous(i)
	ENDDO
	DO i = 0, no_runs
		read(10,520) ff2_value_LOG_previous(i)
	ENDDO
	DO i = 0, no_runs
		read(10,520) ff3_value_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		read(10,520) variation_mass_LOG_previous(i)
	ENDDO
	DO i = 0, no_runs
		read(10,520) variation_sup_LOG_previous (i)
	ENDDO

	DO i = 0, no_runs
		read(10,520) min_u_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		read(10,520) max_u_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		read(10,520) elapsed_time_LOG_previous(i)
	ENDDO

	read(10,520)

	!  READING the latest computed solution and grid points:

	read(10,520) (x(i), i = M1, M2)
	read(10,520)
	read(10,520) (y(j), j = N1, N2)

	read(10,520)

	DO i = M1, M2
		read(10,520) (u(i,j), j = N1, N2)  ! <- serves as the initial state for next running
	ENDDO                                  !                (performed here)

	!  RECREATING and READING previous SOLUTION STATISTICS data:

	read(10,510)
	read(10,510) previous_length    ! <--- previous_length should
	read(10,510)                    !      always be equal to: count

	if ( previous_length .NE. count ) then

		OPEN(11, FILE = 'u24h_serial_ERROR.txt', ACTION = 'WRITE', &
			STATUS = 'REPLACE')
		write(11,*)
		write(11,*) '*******************************************************'
		write(11,*) ' ERROR ***** ERROR ***** ERROR ***** ERROR ***** ERROR '
		write(11,*) '*******************************************************'
		write(11,*)
		write(11,51) DATE_start(3), DATE_start(2), DATE_start(1)
		!51    FORMAT(1X,'DAY :  ',I2.2,' / ',I2.2,' / ',I4)
		write(11,52) DATE_start(5:7)
		!52    FORMAT(1X,'TIME:  ',I2.2,' : ',I2.2,' : ',I2.2)
		write(11,*)
		write(11,*) '*******************************************************'
		write(11,*) ' ERROR CONDITION:'
		write(11,*) ' variable "count" should be equal to '
		write(11,*) ' variable "previous_length"'
		write(11,*)
		write(11,*) ' count = ', count
		write(11,*) ' previous_length = ', previous_length
		write(11,*)
		write(11,*) ' EXECUTION WILL TERMINATE'
		write(11,*)
		write(11,*) '*******************************************************'

		CLOSE(11)

		write(6,*) '*****************************************'
		write(6,*) ' ERROR: see file "u24h_serial_ERROR.txt" '
		write(6,*) '        for explanation '
		write(6,*)
		write(6,*) 'count = ', count
		write(6,*) 'previous_length = ', previous_length
		write(6,*)
		write(6,*) ' ***** EXECUTION HAS BEEN ABORTED *****  '
		write(6,*) '*****************************************'

		STOP

	endif

	ALLOCATE( ff1_value_u_previous(0:previous_length) )
	ALLOCATE( ff2_value_u_previous(0:previous_length) )
	ALLOCATE( ff3_value_u_previous(0:previous_length) )
	ALLOCATE( min_u_previous(0:previous_length) )
	ALLOCATE( max_u_previous(0:previous_length) )
	ALLOCATE( variation_mass_previous(0:previous_length) )
	ALLOCATE( variation_sup_previous (0:previous_length) )
	ALLOCATE( time_per_cycle_previous(0:previous_length) )

	DO i = 0, previous_length
		read(10, 520) ff1_value_u_previous(i), ff2_value_u_previous(i), &
			ff3_value_u_previous(i)
	ENDDO

	read(10,520)
	DO i = 0, previous_length
		read(10, 520) variation_mass_previous(i), variation_sup_previous(i)
	ENDDO

	read(10,520)
	DO i = 0, previous_length
		read(10, 520) min_u_previous(i), max_u_previous(i)
	ENDDO

	read(10,520)
	DO i = 0, previous_length
		read(10, 520) time_per_cycle_previous(i)
	ENDDO

	read(10,520)
	read(10,515) Qt_Atualizacoes_u

	CLOSE(10)

	!-------------------------------------------------------------------------------
	!  CREATING SAFE COPY  u24h_serial_input_previous.dat  OF THE INPUT FILE ABOVE:
	!-------------------------------------------------------------------------------

	OPEN(20, FILE = 'u24h_serial_input_previous.dat', ACTION = 'WRITE', &
		STATUS = 'REPLACE')

	write(20,510) M1, M2
	write(20,510) N1, N2

	write(20,510) n

	write(20,520) p
	write(20,520) h

	write(20,520) cfl

	write(20,520) x_min, x_max
	write(20,520) y_min, y_max

	write(20,520) t0, tF, dt_dump

	write(20,520) x1, x2, x3, x4
	write(20,520) y1, y2, y3, y4
	write(20,510) i1, i2, i3, i4
	write(20,510) j1, j2, j3, j4

	write(20,520) b1, b2, b3, b4
	write(20,520) b

	write(20,520) refv1, refv2, refv3

	write(20,510) count
	write(20,510) no_runs

	write(20,510) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
	write(20,510) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
	write(20,510) length_i1, length_i2, length_i3
	write(20,510) length_j1, length_j2, length_j3

	write(20,520) ff1_value_u0, ff2_value_u0, ff3_value_u0

	write(20,520) min_u0, max_u0
	write(20,520) mass_u0

	write(20,520) t1_min, t2_min, t3_min, t4_min
	write(20,520) t1_max, t2_max, t3_max, t4_max

	write(20,520)

	DO i = 0, no_runs
		write(20,510) (DATE_start_LOG_previous (i,j), j = 1,8)
	ENDDO
	DO i = 0, no_runs
		write(20,510) (DATE_finish_LOG_previous(i,j), j = 1,8)
	ENDDO

	DO i = 0, no_runs
		write(20,520) tF_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		write(20,520) ff1_value_LOG_previous(i)
	ENDDO
	DO i = 0, no_runs
		write(20,520) ff2_value_LOG_previous(i)
	ENDDO
	DO i = 0, no_runs
		write(20,520) ff3_value_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		write(20,520) variation_mass_LOG_previous(i)
	ENDDO
	DO i = 0, no_runs
		write(20,520) variation_sup_LOG_previous (i)
	ENDDO

	DO i = 0, no_runs
		write(20,520) min_u_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		write(20,520) max_u_LOG_previous(i)
	ENDDO

	DO i = 0, no_runs
		write(20,520) elapsed_time_LOG_previous(i)
	ENDDO

	write(20,520)

	write(20,520) (x(i), i = M1, M2)
	write(20,520)
	write(20,520) (y(j), j = N1, N2)

	write(20,520)

	DO i = M1, M2
		write(20,520) (u(i,j), j = N1, N2)   ! <--- current initial state
	ENDDO

	write(20, 510)
	write(20, 510) previous_length
	write(20, 510)

	DO i = 0, previous_length
		write(20,520) ff1_value_u_previous(i), ff2_value_u_previous(i), &
			ff3_value_u_previous(i)
	ENDDO

	write(20,520)
	DO i = 0, previous_length
		write(20,520) variation_mass_previous(i), variation_sup_previous(i)
	ENDDO

	write(20,520)
	DO i = 0, previous_length
		write(20,520) min_u_previous(i), max_u_previous(i)
	ENDDO

	write(20,520)
	DO i = 0, previous_length
		write(20, 520) time_per_cycle_previous(i)
	ENDDO

	write(20,520)
	write(20,515) Qt_Atualizacoes_u

	CLOSE(20)

	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------

	!----------------------------------------------------------
	!              Setting new t0, tF:
	!----------------------------------------------------------


	tF_new = tF + 100D0;
	if ( tF < 50.1D0 ) tF_new = 100D0;
	if ( tF < 20.1D0 ) tF_new =  50D0;
	if ( tF < 10.1D0 ) tF_new =  20D0;
	if ( tF <  8.1D0 ) tF_new =  10D0;
	if ( tF <  6.1D0 ) tF_new =  8D0;
	if ( tF <  4.1D0 ) tF_new =  6D0;
	if ( tF <  2.1D0 ) tF_new =  4D0;
	if ( tF <  1.1D0 ) tF_new =  2D0;    ! <------ eliminate
	if ( tF <  0.6D0 ) tF_new =  1D0;    ! <------ eliminate
	if ( tF <  0.3D0 )  tF_new = 0.5D0;
	if ( tF <  0.15D0 ) tF_new = 0.2D0;  ! <------ eliminate
	if ( tF <  0.05D0 ) tF_new = 0.1D0;  ! <------ eliminate

	t0 = tF;
	tF = tF_new;

	!----------------------------------------------------------
	!     ALLOCATING the new solution-statistics arrays:
	!----------------------------------------------------------

	no_dumps = NINT( (tF-t0)/dt_dump ) + 10;  ! <--- 10 is added here for CAUTION!!
	! (correct value will be given to "no_dumps"
	!  AFTER all dumps have been done and counted)
	new_length = previous_length + no_dumps;  ! <--- this estimate will NOT actually be used

	ALLOCATE( ff1_value_u(1:no_dumps) )
	ALLOCATE( ff2_value_u(1:no_dumps) )
	ALLOCATE( ff3_value_u(1:no_dumps) )

	ALLOCATE( variation_mass(1:no_dumps) )
	ALLOCATE( variation_sup (1:no_dumps) )

	ALLOCATE( min_u(1:no_dumps) )
	ALLOCATE( max_u(1:no_dumps) )

	ALLOCATE( time_per_cycle(1:no_dumps) )


	!----------------------------------------------------------
	!      Extending solution values to extended grid:
	!----------------------------------------------------------

	u(M1-1,N1:N2) = u(M1,N1:N2);
	u(M2+1,N1:N2) = u(M2,N1:N2);
	u(M1-1:M2+1,N1-1) = u(M1-1:M2+1,N1);
	u(M1-1:M2+1,N2+1) = u(M1-1:M2+1,N2);

	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------

	!
	! Computing steady-state limit of solution u(.,t):
	!

	dt = cfl*h*h;  ! <--- time step length

	! --------------------------------------------------------------------
	write(6,*) ' '
	write(6,*) ' '
	write(6,*) ' '
	write(6,*) '***************************************************************'
	write(6,*) '                 Running "u24h_serial4P.f90"... '
	write(6,*) '***************************************************************'
	write(6,'(2x,a,f8.1)') "Current value of tF: ", tF
	write(6,'(2x,a,i9)') "Total amount of u updates since t=0: ", Qt_Atualizacoes_u
	write(6,*) '------------------------------------------------'
	write(6,*) '    Computing ... '

	q = (p - 2D0)/2D0;

	h_to_pm2 = h**(p-2);

	time_next_dump = t0 + dt_dump;

	t = t0;

	u_previous = u;

	new_count = 0;

	no_runs = no_runs + 1;

	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------
	!                    MAIN COMPUTATION SECTION
	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------

	tF_check = tF - dt/4D0;

	DO WHILE ( t < tF_check )

		!CALL DATE_AND_TIME(VALUES = DATE_values)
		!time_start = dble(DATE_values(7)) + dble(DATE_values(8))/1000D0;
		CALL CPU_TIME(time_start)


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
				cfl*( F(M1+1:M2+1,N1:N2)*(u(M1+1:M2+1,N1:N2)-u(M1:M2,N1:N2)) - &
				F(M1:M2,N1:N2)*(u(M1:M2,N1:N2)-u(M1-1:M2-1,N1:N2)) ) + &
				cfl*( G(M1:M2,N1+1:N2+1)*(u(M1:M2,N1+1:N2+1)-u(M1:M2,N1:N2)) - &
				G(M1:M2,N1:N2)*(u(M1:M2,N1:N2)-u(M1:M2,N1-1:N2-1)) );

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

			Qt_Atualizacoes_u = Qt_Atualizacoes_u + 1

		ENDDO ! DO WHILE ( t < time_next_dump - dt/4 )

		!--- time to save solution statistics:

		!CALL DATE_AND_TIME(VALUES = DATE_values)
		!time_finish = dble(DATE_values(7)) + dble(DATE_values(8))/1000D0;
		CALL CPU_TIME(time_finish)

		new_count = new_count + 1;

		time_per_cycle(new_count) = time_finish - time_start;

		!--- solution far-field values:

		ff1_value_u(new_count) = ( sum( u(i_ff1a:i_ff1a+10,j_ff1a:j_ff1b) )/(11*length_j1) &
			+ sum( u(i_ff1b-10:i_ff1b,j_ff1a:j_ff1b) )/(11*length_j1) &
			+ sum( u(i_ff1a:i_ff1b,j_ff1a:j_ff1a+10) )/(11*length_i1) &
			+ sum( u(i_ff1a:i_ff1b,j_ff1b-10:j_ff1b) )/(11*length_i1) )/4D0;

		ff2_value_u(new_count) = ( sum( u(i_ff2a:i_ff2a+10,j_ff2a:j_ff2b) )/(11*length_j2) &
			+ sum( u(i_ff2b-10:i_ff2b,j_ff2a:j_ff2b) )/(11*length_j2) &
			+ sum( u(i_ff2a:i_ff2b,j_ff2a:j_ff2a+10) )/(11*length_i2) &
			+ sum( u(i_ff2a:i_ff2b,j_ff2b-10:j_ff2b) )/(11*length_i2) )/4D0;

		ff3_value_u(new_count) = ( sum( u(i_ff3a:i_ff3a+10,j_ff3a:j_ff3b) )/(11*length_j3) &
			+ sum( u(i_ff3b-10:i_ff3b,j_ff3a:j_ff3b) )/(11*length_j3) &
			+ sum( u(i_ff3a:i_ff3b,j_ff3a:j_ff3a+10) )/(11*length_i3) &
			+ sum( u(i_ff3a:i_ff3b,j_ff3b-10:j_ff3b) )/(11*length_i3) )/4D0;

		!--------------------------------------------------

		! saving minimum and maximum solution values at current time:

		min_u(new_count) = minval( u );
		max_u(new_count) = maxval( u );

		! computing supnorm of solution variation from previous u:

		variation_sup(new_count) = maxval( abs( u - u_previous ) );

		! computing mass of solution variation from previous u:

		variation_mass(new_count) = sum( u - u_previous )*h**2;

		! current solution values become previous solution values
		! for the next iteration cycle:

		u_previous = u;

		! printing out some statistics before going to next cycle:
		!
		! write(6,50) (t-t0)/(tF-t0)
		! 50 FORMAT( 3E15.6 )
		! write(6,50) variation_mass(new_count), variation_sup(new_count), ff2_value_u(new_count)
		!
		! done

		! updating time of next dumping for solution statistics:

		time_next_dump = time_next_dump + dt_dump;

	ENDDO !DO WHILE ( t < tF_check )

	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------
	!                  END OF MAIN COMPUTATION SECTION
	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------

	total_elapsed_time = sum( time_per_cycle );

	count = count + new_count;   ! <--- updating "count" (= size of new total statistics arrays)
	no_dumps = new_count;        ! <--- this puts to "no_dumps" the final, CORRECT value of the
	!      actual number of dumps (which was ESTIMATED from above
	!                           by the previous value of  no_dumps)

	CALL DATE_AND_TIME(VALUES = DATE_finish)

	write(6,*) '*******************************************'
	write(6,'(2x,a)') 'Done! '
	write(6,*) '*******************************************'
	write(6,'(2x,a)') 'Saving data to output files ...'

	duracao = 3600*(DATE_finish(5)-DATE_start(5))+ 60*(DATE_finish(6)-DATE_start(6)) + &
		(DATE_finish(7)-DATE_start(7))

	write(6,*) '------------------------------------------------------------'
	write(6,'(2x,a,i0,a)') 'TOTAL elapsed time: ',  duracao, ' segundos'
	write(6,*) '------------------------------------------------------------'
	write(6,*) '             SUCCESSFUL EXECUTION                           '
	write(6,*) '------------------------------------------------------------'


	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------
	!                       SAVING OUTPUT DATA TO DISK
	!      (this data will serve as new input data for the program
	!                "u24h_serial.f90" on its next execution)
	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------

	OPEN(50, FILE = 'u24h_serial_INPUT.dat', ACTION = 'WRITE', &
		STATUS = 'REPLACE')
	530 FORMAT( 10(1X,I10))
	535 FORMAT( 10(1X,I12) )
	540 FORMAT( 10(1X,E25.16) )

	write(50,530) M1, M2  ! ind=2
	write(50,530) N1, N2

	write(50,530) n

	write(50,540) p
	write(50,540) h

	write(50,540) cfl

	write(50,540) x_min, x_max !ind x_max =10
	write(50,540) y_min, y_max

	write(50,540) t0, tF, dt_dump

	write(50,540) x1, x2, x3, x4
	write(50,540) y1, y2, y3, y4
	write(50,530) i1, i2, i3, i4
	write(50,530) j1, j2, j3, j4 ! ind j4 =31

	write(50,540) b1, b2, b3, b4
	write(50,540) b

	write(50,540) refv1, refv2, refv3

	write(50,530) count
	write(50,530) no_runs  ! ind no_runs = 41

	write(50,530) i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b
	write(50,530) j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b
	write(50,530) length_i1, length_i2, length_i3
	write(50,530) length_j1, length_j2, length_j3

	write(50,540) ff1_value_u0, ff2_value_u0, ff3_value_u0

	write(50,540) min_u0, max_u0
	write(50,540) mass_u0

	write(50,540) t1_min, t2_min, t3_min, t4_min
	write(50,540) t1_max, t2_max, t3_max, t4_max

	write(50,540)

	! Saving LOG variables: VVVVV

	DO i = 0, no_runs-1
		write(50,530) (DATE_start_LOG_previous (i,j), j = 1,8)
	ENDDO
	write(50,530) (DATE_start (j), j = 1,8)

	DO i = 0, no_runs-1
		write(50,530) (DATE_finish_LOG_previous(i,j), j = 1,8)
	ENDDO
	write(50,530) (DATE_finish(j), j = 1,8)

	DO i = 0, no_runs-1
		write(50,540) tF_LOG_previous(i)
	ENDDO
	write(50,540) tF

	DO i = 0, no_runs-1
		write(50,540) ff1_value_LOG_previous(i)
	ENDDO
	write(50,540) ff1_value_u(new_count)

	DO i = 0, no_runs-1
		write(50,540) ff2_value_LOG_previous(i)
	ENDDO
	write(50,540) ff2_value_u(new_count)

	DO i = 0, no_runs-1
		write(50,540) ff3_value_LOG_previous(i)
	ENDDO
	write(50,540) ff3_value_u(new_count)

	DO i = 0, no_runs-1
		write(50,540) variation_mass_LOG_previous(i)
	ENDDO
	write(50,540) variation_mass(new_count)

	DO i = 0, no_runs-1
		write(50,540) variation_sup_LOG_previous (i)
	ENDDO
	write(50,540) variation_sup(new_count)

	DO i = 0, no_runs-1
		write(50,540) min_u_LOG_previous(i)
	ENDDO
	write(50,540) min_u(new_count)

	DO i = 0, no_runs-1
		write(50,540) max_u_LOG_previous(i)
	ENDDO
	write(50,540) max_u(new_count)

	DO i = 0, no_runs-1
		write(50,540) elapsed_time_LOG_previous(i)
	ENDDO
	write(50,540) total_elapsed_time

	! LOG variables saved ^^^^^

	write(50,540)

	write(50,540) (x(i), i = M1, M2)
	write(50,540)
	write(50,540) (y(j), j = N1, N2)

	write(50,540)

	DO i = M1, M2
		write(50,540) (u(i,j), j = N1, N2)   ! <--- current final solution
	ENDDO

	! Saving current SOLUTION STATISTICS variables:

	new_length = previous_length + new_count;  ! <--- same as newly updated value of "count"

	write(50,530)
	write(50,530) new_length   ! <--- new_length = previous_length + new_count
	write(50,530)

	DO i = 0, previous_length
		write(50,540) ff1_value_u_previous(i), ff2_value_u_previous(i), &
			ff3_value_u_previous(i)
	ENDDO
	DO i = 1, new_count
		write(50,540) ff1_value_u(i), ff2_value_u(i), ff3_value_u(i)
	ENDDO

	write(50,540)
	DO i = 0, previous_length
		write(50,540) variation_mass_previous(i), variation_sup_previous(i)
	ENDDO
	DO i = 1, new_count
		write(50,540) variation_mass(i), variation_sup(i)
	ENDDO

	write(50,540)
	DO i = 0, previous_length
		write(50,540) min_u_previous(i), max_u_previous(i)
	ENDDO
	DO i = 1, new_count
		write(50,540) min_u(i), max_u(i)
	ENDDO

	write(50,540)
	DO i = 0, previous_length
		write(50,540) time_per_cycle_previous(i)
	ENDDO
	DO i = 1, new_count
		write(50,540) time_per_cycle(i)
	ENDDO

	write(50,540)
	write(50,535) Qt_Atualizacoes_u

	CLOSE(50)

	!-------------------------------------------------------
	!    Generating  LOG FILE  for current run:
	!-------------------------------------------------------

	580 FORMAT((2x,a,f5.2) )
	590 FORMAT((2x,a,f5.2,5x,a,f5.2))
	OPEN(70, FILE = 'u24h_serial_LOG.txt', ACTION = 'WRITE', &
		STATUS = 'REPLACE')
	write(70,*)
	write(70,*)
	write(70,*) '***************************************************************************'
	write(70,*) '*                     LOG DATA for u24h_serial.f90                        *'
	write(70,*) '***************************************************************************'
	write(70,*)
	write(70,'(2x,a,f7.1)') "tF = ", tF
	write(70,*)
	write(70,'(2x,a,i0,5x,a,i0)')  "M1 = ", M1, "M2 = ", M2
	write(70,'(2x,a,i0,5x,a,i0)')  "N1 = ", N1, "N2 = ", N2
	write(70,*)
	write(70,*)
	write(70,'(2x,a,f5.1,5x,a,f5.1)') "x_min = ", x_min, "   x_max = ", x_max
	write(70,'(2x,a,f5.1,5x,a,f5.1)') "y_min = ", y_min, "   y_max = ", y_max
	write(70,*)
	write(70,'(2x,a,f4.1)') "refv1 = ", refv1
	write(70,'(2x,a,f4.1)') "refv2 = ", refv2
	write(70,'(2x,a,f4.1)') "refv3 = ", refv3
	write(70,*)
	write(70,*) " Dimension n = ", n
	write(70,*)
	write(70,'(2x,a,f3.1)') "p = ", p
	write(70,*)
	write(70,580) "h = ", h
	write(70,'(2x,a,f7.5)') "dt = ", dt
	write(70,'(2x,a,f7.4)') "dt_dump = ", dt_dump
	write(70,'(2x,a,f7.4)') "cfl = ", cfl
	write(70,*)
	write(70,*) " Singularities:"
	write(70,590) "x1 = ", x1, "y1 = ", y1
	write(70,590) "x2 = ", x2, "y2 = ", y2
	write(70,590) "x3 = ", x3, "y3 = ", y3
	write(70,590) "x4 = ", x4, "y4 = ", y4
	write(70,580) "b1 = ", b1
	write(70,580) "b2 = ", b2
	write(70,580) "b3 = ", b3
	write(70,580) "b4 = ", b4
	write(70,*)
	write(70,580) "Value of b used: ", b
	write(70,580) "Expected value of b: ", (b1 + b2 + b3 + b4)/4D0
	write(70,*)
	write(70,'(2x,a,f7.1)') "Current value of tF: ", tF
	write(70,'(2x,a,i0)') "Atualizacoes de u desde t=0: ", Qt_Atualizacoes_u
	write(70,*)
	write(70,'(2x,a,f10.2)') "Elapsed CPU TIME (master thread)  : ", total_elapsed_time
	write(70,580) "Menor valor em u: ", min_u(new_count)
	write(70,580) "Maior valor em u: ", max_u(new_count)
	write(70,*)
	write(70,*) '***************************************************************************'
	write(70,*)
	write(70,*)
	!-----------------------------------------------------------------------
	!  TABLE 1: Run no., DATE/TIME (start & finish), tF, FP count, FLOP/s:
	!-----------------------------------------------------------------------
	701 FORMAT(T2,'|-----------','|---------------------',4('|----------------'),'|')
	702 FORMAT(T2,'|', T14,'|', T36,'|', T53,'|', T70,'|', T87, '|', T104, '|')       ! <--- blank line
	703 FORMAT(T2,'|', T5,'Run no.', T14,'|', T20,'DATE / TIME', T36,'|', &
		T44,'tF', T53,'|', T58,'FP count', T70,'|', T76,'FLOP/s', T87,'|', T89,'Duracao(segs)',T104, '|')

	704 FORMAT(T2,'|', T14,'|', T16,I2.2,'/',I2.2,'/',I4,1X,I2.2,':',I2.2,':',I2.2,&
		T36,'|', T53,'|', T70,'|', T87,'|', T104, '|')
	705 FORMAT(T2,'|', T6,I3, T14,'|', T16,I2.2,'/',I2.2,'/',I4,1X,I2.2,':',I2.2,':',I2.2,&
		T36,'|', T38,E13.6, T53,'|', T55,E13.6, T70,'|', T72,E13.6, T87,'|', T89, I8, T104, '|')

	write(70,701)
	write(70,702)
	write(70,703)
	write(70,702)
	write(70,701)
	i = 0;
	FP_count = 20*(M2-M1+1)*(N2-N1+1)* dt_dump/dt;
	FLOPs = FP_count / elapsed_time_LOG_previous(0);
	duracao = 3600*(DATE_finish_LOG_previous(i,5)-DATE_start_LOG_previous(i,5))+&
		60*(DATE_finish_LOG_previous(i,6)-DATE_start_LOG_previous(i,6)) + &
		(DATE_finish_LOG_previous(i,7)-DATE_start_LOG_previous(i,7))
	write(70,704) DATE_start_LOG_previous(i,3), DATE_start_LOG_previous(i,2), &
		DATE_start_LOG_previous(i,1), DATE_start_LOG_previous(i,5:7)
	write(70,705) i, DATE_finish_LOG_previous(i,3), DATE_finish_LOG_previous(i,2), &
		DATE_finish_LOG_previous(i,1), DATE_finish_LOG_previous(i,5:7), &
		tF_LOG_previous(i), FP_count, FLOPs,duracao
	write(70,701)
	IF ( no_runs > 1 ) THEN
		do i = 1, no_runs-1
			FP_count = 20*(M2-M1+1)*(N2-N1+1)*(tF_LOG_previous(i) - tF_LOG_previous(i-1))/dt;
			FLOPs = FP_count / elapsed_time_LOG_previous(i);
			duracao = 3600*(DATE_finish_LOG_previous(i,5)-DATE_start_LOG_previous(i,5))+&
				60*(DATE_finish_LOG_previous(i,6)-DATE_start_LOG_previous(i,6)) + &
				(DATE_finish_LOG_previous(i,7)-DATE_start_LOG_previous(i,7))

			write(70,704) DATE_start_LOG_previous(i,3), DATE_start_LOG_previous(i,2), &
				DATE_start_LOG_previous(i,1), DATE_start_LOG_previous(i,5:7)
			write(70,705) i, DATE_finish_LOG_previous(i,3), DATE_finish_LOG_previous(i,2), &
				DATE_finish_LOG_previous(i,1), DATE_finish_LOG_previous(i,5:7), &
				tF_LOG_previous(i), FP_count, FLOPs, duracao
			write(70,701)
		enddo
	ENDIF
	i = no_runs;
	FP_count = 20*(M2-M1+1)*(N2-N1+1)*(tF - tF_LOG_previous(i-1))/dt;
	FLOPs = FP_count / total_elapsed_time;
	duracao = 3600*(DATE_finish(5)-DATE_start(5))+ 60*(DATE_finish(6)-DATE_start(6)) + &
		(DATE_finish(7)-DATE_start(7))

	write(70,704) DATE_start(3), DATE_start(2), DATE_start(1), DATE_start(5:7)
	write(70,705) i, DATE_finish(3), DATE_finish(2), DATE_finish(1), DATE_finish(5:7), &
		tF, FP_count, FLOPs,duracao
	write(70,701)

	write(70,*)
	write(70,*)
	write(70,*)

	!-----------------------------------------------------------------------
	!  TABLE 2: Run no., ff1_value, ff2_value, ff3_value, elapsed time:
	!-----------------------------------------------------------------------

	711 FORMAT(T2,'|-----------',4('|----------------'),'|')
	712 FORMAT(T2,'|', T14,'|', T31,'|', T48,'|', T65,'|', T82,'|')       ! <--- blank line
	713 FORMAT(T2,'|', T5,'Run no.', T14,'|', T21,'ffv1', T31,'|', &
		T38,'ffv2', T48,'|', T55,'ffv3', T65,'|', T68,'elapsed time', T82,'|')
	714 FORMAT(T2,'|', T14,'|', T31,'|', T48,'|', T65,'|', T69, '(seconds)', T82,'|')
	715 FORMAT(T2,'|', T6,I3, T14,'|', T16,E13.6, T31,'|', T33,E13.6, T48,'|', &
		T50,E13.6, T65,'|', T67,E13.6, T82,'|')

	write(70,711)
	write(70,712)
	write(70,713)
	write(70,714)
	write(70,711)

	do i = 0, no_runs-1
		write(70,715) i, ff1_value_LOG_previous(i), ff2_value_LOG_previous(i), &
			ff3_value_LOG_previous(i), elapsed_time_LOG_previous(i)
		write(70,711)
	enddo
	i = no_runs;
	write(70,715) i, ff1_value_u(new_count), ff2_value_u(new_count), &
		ff3_value_u(new_count), total_elapsed_time
	write(70,711)

	write(70,*)
	write(70,*)
	write(70,*)

	!--------------------------------------------------------------------------
	!  TABLE 3: Run no., variation_mass, variation_sup, mean_du/dt, max_du/dt:
	!--------------------------------------------------------------------------

	721 FORMAT(T2,'|-----------',4('|----------------'),'|')
	722 FORMAT(T2,'|', T14,'|', T31,'|', T48,'|', T65,'|', T82,'|')       ! <--- blank line
	723 FORMAT(T2,'|', T5,'Run no.', T14,'|', T16,'variation_mass', T31,'|', &
		T34,'variation_sup', T48,'|', T52,'mean_du/dt', T65,'|', T70,'max_du/dt', T82,'|')
	724 FORMAT(T2,'|', T6,I3, T14,'|', T18,'----------', T31,'|', T35,'----------', &
		T48,'|', T50,E13.6, T65,'|', T67,E13.6, T82,'|')
	725 FORMAT(T2,'|', T6,I3, T14,'|', T16,E13.6, T31,'|', T33,E13.6, T48,'|', &
		T50,E13.6, T65,'|', T67,E13.6, T82,'|')

	write(70,721)
	write(70,722)
	write(70,723)
	write(70,722)
	write(70,721)

	i = 0;
	! write(70,724) i, min_u_LOG_previous(0), max_u_LOG_previous(0)
	! write(70,721)
	do i = 0, no_runs-1
		write(70,725) i, variation_mass_LOG_previous(i), variation_sup_LOG_previous(i), &
			variation_mass_LOG_previous(i)/(dt_dump *(M2-M1+1)*(N2-N1+1)), &
			variation_sup_LOG_previous(i)/dt_dump
		write(70,721)
	enddo
	i = no_runs;
	write(70,725) i, variation_mass(new_count), variation_sup(new_count), &
		variation_mass(new_count)/(dt_dump *(M2-M1+1)*(N2-N1+1)), &
		variation_sup(new_count)/dt_dump
	write(70,721)

	CLOSE(70)

	!--------------------------------------------------------------------------
	!--------------------------------------------------------------------------

	!-----------------------------------------
	!  Deallocating variables from memory:
	!-----------------------------------------

	DEALLOCATE( x, y )
	DEALLOCATE( u, v )
	DEALLOCATE( u_previous )
	DEALLOCATE( F, G )
	DEALLOCATE( ff1_value_u, ff2_value_u, ff3_value_u )
	DEALLOCATE( min_u, max_u )
	DEALLOCATE( variation_mass, variation_sup )
	DEALLOCATE( time_per_cycle )

	DEALLOCATE( ff1_value_u_previous )
	DEALLOCATE( ff2_value_u_previous )
	DEALLOCATE( ff3_value_u_previous )
	DEALLOCATE( min_u_previous )
	DEALLOCATE( max_u_previous )
	DEALLOCATE( variation_mass_previous )
	DEALLOCATE( variation_sup_previous )
	DEALLOCATE( time_per_cycle_previous )

	DEALLOCATE( DATE_start_LOG_previous )
	DEALLOCATE( DATE_finish_LOG_previous )
	DEALLOCATE( ff1_value_LOG_previous )
	DEALLOCATE( ff2_value_LOG_previous )
	DEALLOCATE( ff3_value_LOG_previous )
	DEALLOCATE( variation_mass_LOG_previous )
	DEALLOCATE( variation_sup_LOG_previous )
	DEALLOCATE( min_u_LOG_previous )
	DEALLOCATE( max_u_LOG_previous )
	DEALLOCATE( tF_LOG_previous )
	DEALLOCATE( elapsed_time_LOG_previous )

	!-----------------------------------------
	!       Deallocation completed
	!-----------------------------------------

end do !do numLoops = 1, 2000; + or - line 170

END PROGRAM u24h_serial4P
