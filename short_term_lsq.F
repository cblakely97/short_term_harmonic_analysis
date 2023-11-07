!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!     This program serves as a place to test implementing a short-    !
! term (1-3 day) least squares harmonic analysis of a time history    !
! signal and then use that analysis to predict the next couple of     !
! days. It is a relatively simple program wherein two constituents    !
! (a diurnal and semidiurnal) are used. The harmonic equation is of   !
! the form:                                                           !
!                                                                     !
!       y(t) = sum_{i=1}^n ( c_i*cos(a_i*t) + s_i*sin(a_i*t) )        !
!                                                                     !
! Where,                                                              !
!                                                                     !
!       y = the signal being solved for (\eta, u, v, etc.)            !
!       t = time of interest                                          !
!       n = number of constituents (2 in this case)                   !
!       c_i, s_i = amplitudes being solved for                        !
!       a_i = frequencies of interest                                 !
!                                                                     !
! By grouping nodal factors, equilibrium arguments, phase lags, and   !
! amplitudes into c_i and s_i, we greatly simplify the problem and    !
! are able to directly use them to solve for the next couple of days  !
! without worrying about long-term considerations.                    !
!                                                                     !
! Written by: Coleman Blakely 11/3/2023                               !
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      PROGRAM SHORT_LSQ_HA_program
      IMPLICIT NONE
      REAL(8),ALLOCATABLE :: u(:,:), v(:,:), t(:)       ! signal and times
      REAL(8),ALLOCATABLE :: xp(:), yp(:)               ! for output
      REAL(8),ALLOCATABLE :: u_out(:,:), v_out(:,:)         ! output
      REAL(8),ALLOCATABLE :: UAV(:,:),VAV(:,:)          ! used to store samples
                                                ! for 25-hour filter
      REAL(8),ALLOCATABLE :: freq(:)       ! frequencies
      ! vectors for storing c_i and s_i for resynthesis
      REAL(8),ALLOCATABLE,DIMENSION(:,:,:) :: c_i, s_i
      INTEGER :: NP     ! number of "nodes"
      INTEGER :: nconst ! number of constituents
      REAL(8) :: DTDP   ! "timestep" of times (seconds)
      REAL(8) :: LSQ_L  ! length of analysis (days)
      REAL(8) :: Fs     ! sampling frequency (seconds)
      REAL(8) :: SimL   ! length of "simulation" (days)
      INTEGER :: TimeStep ! current timestep
      INTEGER :: NT     ! number of time steps
      INTEGER :: NSCREEN ! print to screen every nscreen time steps 
      ! variables for least squares analysis (Ax=b)
      ! note that b will be of size 2 by 2*nconst by NP (u and v dirs)
      REAL(8),ALLOCATABLE :: A(:,:), b(:,:,:)
      ! parameters
      REAL(8),PARAMETER :: pi = 3.141592653589793d0
      REAL(8),PARAMETER :: cpd2radpsec = 2d0*pi/86400d0
      ! other variables
      INTEGER :: it, iconst, kk, counter ! for loops
      ! output estimate
      REAL(8),ALLOCATABLE :: UBar(:), VBar(:), UTide(:), VTide(:)
      ! filename
      CHARACTER(LEN=16) :: f62="Baroclinic.62.nc"
      CHARACTER(LEN=11) :: fn_out="file_out.nc"
      CHARACTER(LEN=40) ::screenmessage
      REAL(8) :: tmp
      !
      ! Read in velocities from a netcdf file
      !
      CALL READ_f62(f62)
      !
      ! Define input values
      !
      nconst = 2        ! I am writing it this way to make it more
                        ! portable if we want to expand later.
      ALLOCATE( freq(nconst), c_i(2,nconst,NP), s_i(2,nconst,NP) )
      freq = (/1.9323, 1.0027/)  ! diurnal (M2) and semidiurnal (K1)
                                 ! frequencies (cycles per day)
      freq = freq*cpd2radpsec    ! convert to radians per second
      LSQ_L = 2d0       ! days
      Fs = 360d0       ! seconds
      SimL = t(NT) - t(1)       ! seconds
      DTDP = t(2) - t(1)        ! seconds
      ! put everything in seconds
      LSQ_L = LSQ_L*86400d0
      NSCREEN = FLOOR(0.1d0*DBLE(NT))
      !
      counter = 0
      DO it = 1,NT
         CALL SHORT_LSQ_HA(u(:,it),v(:,it),it)
         counter = counter + 1
         IF (counter.EQ.NSCREEN) THEN
            tmp = (DBLE(IT)/DBLE(NT))*100d0
            WRITE(screenmessage,'(F6.2,A6)') tmp,"% done" 
            WRITE(*,*) trim(screenmessage)
            counter = 0
         ENDIF
      ENDDO
      CALL write_output(fn_out)
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      CONTAINS
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE write_output(filename)
!---------------------------------------------------------------------!
!     Calculates a time series y(t) based on input amplitudes c and s !
! and timestep and number of time steps.                              !
!---------------------------------------------------------------------!
      USE netcdf
      IMPLICIT NONE
      !
      ! arguments
      !
      CHARACTER(11),INTENT(IN) :: filename
      !
      ! local
      !
      INTEGER :: ncid ! file id
      INTEGER :: ierr ! err code
      INTEGER :: time_dim_id, NP_dim_id
      INTEGER :: time_id, u_id, v_id, x_id, y_id
      CALL NCERROR(NF90_CREATE(filename,NF90_CLOBBER,ncid))
      ! define dimensions
      CALL NCERROR(NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED,time_dim_id))
      CALL NCERROR(NF90_DEF_DIM(ncid,"NP",NP,NP_dim_id))
      ! define variables
      CALL NCERROR(NF90_DEF_VAR(ncid,"x",NF90_DOUBLE,NP_dim_id, x_id))
      CALL NCERROR(NF90_PUT_ATT(ncid,x_id,'long_name','longitude'))
      CALL NCERROR(NF90_PUT_ATT(ncid,x_id,'units','degrees east'))
      CALL NCERROR(NF90_DEF_VAR(ncid,"y",NF90_DOUBLE,NP_dim_id, y_id))
      CALL NCERROR(NF90_PUT_ATT(ncid,y_id,'long_name','latitude'))
      CALL NCERROR(NF90_PUT_ATT(ncid,y_id,'units','degrees north'))
      CALL NCERROR(NF90_DEF_VAR(ncid,"time",NF90_DOUBLE,time_dim_id,
     &                          time_id))
      CALL NCERROR(NF90_PUT_ATT(ncid,time_id,'long_name','model time'))
      CALL NCERROR(NF90_PUT_ATT(ncid,time_id,'units',
     &                     'seconds since 2016-10-03 00:00:00 UTC'))
      ! velocities
      CALL NCERROR(NF90_DEF_VAR(ncid,"UTide",NF90_DOUBLE,
     &                    (/NP_dim_id,time_dim_id/),u_id))
      CALL NCERROR(NF90_PUT_ATT(ncid,u_id,'long_name',
     &             'Zonal Tidal Velocity'))
      CALL NCERROR(NF90_PUT_ATT(ncid,u_id,'units','m s-1'))
      CALL NCERROR(NF90_DEF_VAR(ncid,"VTide",NF90_DOUBLE,
     &                    (/np_dim_id,time_dim_id/),v_id))
      CALL NCERROR(NF90_PUT_ATT(ncid,v_id,'long_name',
     &                     'Meridional Tidal Velocity'))
      CALL NCERROR(NF90_PUT_ATT(ncid,v_id,'units','m s-1'))
      CALL NCERROR(NF90_ENDDEF(ncid))
      ! put on x, y
      CALL NCERROR(NF90_PUT_VAR(ncid,x_id,xp))
      CALL NCERROR(NF90_PUT_VAR(ncid,y_id,yp))
      ! put on time
      CALL NCERROR(NF90_PUT_VAR(ncid,time_id,t))
      ! u and v
      CALL NCERROR(NF90_PUT_VAR(ncid,u_id,u_out))
      CALL NCERROR(NF90_PUT_VAR(ncid,v_id,v_out))
      CALL NCERROR(NF90_CLOSE(ncid))
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      END SUBROUTINE write_output
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE READ_f62(filename)
!---------------------------------------------------------------------!
!     Calculates a time series y(t) based on input amplitudes c and s !
! and timestep and number of time steps.                              !
!---------------------------------------------------------------------!
      USE netcdf
      IMPLICIT NONE
      !
      ! arguments
      !
      CHARACTER(16),INTENT(IN) :: filename
      !
      ! local
      !
      INTEGER :: ncid ! file id
      INTEGER :: ierr ! err code
      INTEGER :: tempid ! for grabbing variables, etc
      CALL NCERROR(NF90_OPEN(filename,NF90_NOWRITE,ncid))
      ! read in number of stations to NP
      CALL NCERROR(NF90_INQ_DIMID(ncid,'station',tempid))
      CALL NCERROR(NF90_INQUIRE_DIMENSION(ncid,tempid,len=NP))
      ! read in number of time snaps
      CALL NCERROR(NF90_INQ_DIMID(ncid,'time',tempid))
      CALL NCERROR(NF90_INQUIRE_DIMENSION(ncid,tempid,len=NT))
      !
      ! allocate and read in time, u-vel, v-vel
      !
      ALLOCATE( t(NT), u(NP,NT), v(NP,NT) )
      !
      ! allocate x and y for outputting
      !
      ALLOCATE( xp(NP), yp(NP) )
      !
      ! get time, u-vel, v-vel, x, and y from f62
      !
      CALL NCERROR(NF90_INQ_VARID(ncid,'time',tempid))
      CALL NCERROR(NF90_GET_VAR(ncid,tempid,t))
      CALL NCERROR(NF90_INQ_VARID(ncid,'u-vel',tempid))
      CALL NCERROR(NF90_GET_VAR(ncid,tempid,u))
      CALL NCERROR(NF90_INQ_VARID(ncid,'v-vel',tempid))
      CALL NCERROR(NF90_GET_VAR(ncid,tempid,v))
      CALL NCERROR(NF90_INQ_VARID(ncid,'x',tempid))
      CALL NCERROR(NF90_GET_VAR(ncid,tempid,xp))
      CALL NCERROR(NF90_INQ_VARID(ncid,'y',tempid))
      CALL NCERROR(NF90_GET_VAR(ncid,tempid,yp))
      ! 
      CALL NCERROR(NF90_CLOSE(ncid))
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      END SUBROUTINE READ_f62
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE NCERROR(ierr)
!---------------------------------------------------------------------!
!     Calculates a time series y(t) based on input amplitudes c and s !
! and timestep and number of time steps.                              !
!---------------------------------------------------------------------!
      USE netcdf
      IMPLICIT NONE
      !
      ! arguments
      !
      INTEGER,INTENT(IN) :: ierr
      REAL(8),ALLOCATABLE :: dmy(:)
      IF (ierr.NE.NF90_NOERR) THEN
         write(*,*) trim(NF90_STRERROR(ierr))
         dmy(1) = 999d0
         CALL ABORT('huh?')
      ENDIF
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      END SUBROUTINE NCERROR
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE SHORT_LSQ_HA(Uin,Vin, TimeStep)
!---------------------------------------------------------------------!
!     subroutine that takes the current y value and timetep and, if   !
! it is time, performs the least squares analysis.                    !
!---------------------------------------------------------------------!
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: Uin(:), Vin(:)
      INTEGER,INTENT(IN) :: TimeStep
      INTEGER :: L, Lm ! for checking if it is time to update/solve
      ! first_call is my hacky-ass way of initializing. I hopefully will
      ! fix it before actually implementing in ADCIRC
      LOGICAL,SAVE :: first_call = .true.
      ! This next flag is also a bit hacky but I use it to tell the
      ! program when the first HA has been performed. 
      LOGICAL,SAVE :: first_solve = .false.
      ! these next two flags serve to check if it is time to update the
      ! mass matrix and load vector (Ax=b) and also if it is time to
      ! solve the harmonic analysis
      LOGICAL :: update_flag = .false., solve_flag = .false.
      INTEGER,SAVE :: solve_counter = 0
      ! To avoid roundoff issues in calculating sines/cosines we will
      ! calculate a "harmonic analysis time" that is in units of seconds
      ! since the last harmonic analysis. We will calculate this using
      ! solve_counter
      REAL(8) :: HA_Time
      INTEGER :: ii
      !
      ! Calculate harmonic time before solve_counter is messed with
      !
      HA_time = TimeStep*DTDP
      !
      ! Check if it is time to update/solve
      !
      ! to check if it is time to update check if this timestep and the
      ! previous timestep are in different sampling periods
      L = FLOOR(TimeStep*DTDP/Fs)
      Lm = FLOOR((TimeStep-1)*DTDP/Fs)
      IF ( L.GT.Lm ) THEN
         update_flag = .true.
      ELSE
         update_flag = .false.
      ENDIF
      ! to check if it is time to solve we just use a counter to track
      ! how long it has been since we last solved the analysis
      L = solve_counter*DTDP
      IF ( L.EQ.LSQ_L ) THEN
         first_solve = .true.
         solve_flag = .true.
         solve_counter = 0
      ELSE
         solve_flag = .false.
      ENDIF

      IF ( first_call ) THEN
         ALLOCATE( A(2*nconst,2*nconst), b(2,2*nconst,NP) )
         A = 0d0
         b = 0d0
         ALLOCATE( UAV(25,NP), VAV(25,NP), 
     &             UBar(NP), VBar(NP),
     &             UTide(NP), VTide(NP) )
         ALLOCATE( u_out(NP,NT), v_out(NP,NT) )
         UAV = 0d0
         VAV = 0d0
         UBar = 0d0
         VBar = 0d0
         UTide = 0d0
         VTide = 0d0
         u_out = 0d0
         v_out = 0d0
         first_call = .false.
      ENDIF
      !
      ! Update LA25 filter if we haven't gotten to the LSQHA
      !
      IF (.NOT.first_solve) THEN
         CALL UNTIDE(TimeStep,Uin, Vin)
      ENDIF
      ! 
      ! Update A matrix and b vector
      !
      CALL UPDATE_LSQ_HA(Uin,Vin,HA_Time,update_flag)
      !
      ! Solve the least squares analysis
      !
      CALL SOLVE_LSQ_HA(solve_flag)
      !
      ! Estimate output
      !
      CALL ESTIMATE_TIDAL_VELO(Uin,Vin,HA_time,first_solve)
      DO ii = 1,NP
         u_out(ii,TimeStep) = Utide(ii)
         v_out(ii,TimeStep) = Vtide(ii)
      ENDDO
      solve_counter = solve_counter + 1
      return
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDSUBROUTINE SHORT_LSQ_HA
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE UNTIDE(TimeStep,Uin,Vin)
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: Uin(:),Vin(:)
      INTEGER,INTENT(IN) :: TimeStep
      INTEGER :: L, Lm ! for checking if it is time to update/solve
      INTEGER,SAVE :: ISTA = 1
      INTEGER :: NS=25
      INTEGER :: ii,kk
      ! see if it is time to update the samples
      L = FLOOR(TimeStep*DTDP/Fs)
      Lm = FLOOR((TimeStep-1)*DTDP/Fs)
      IF (L.GT.Lm) THEN
         IF (ISTA.GT.NS) THEN
            DO ii = 1,NP
               DO kk = 1,NS-1
                  UAV(ii,kk) = UAV(ii,kk+1)
                  VAV(ii,kk) = VAV(ii,kk+1)
               ENDDO
               UAV(ii,NS) = Uin(ii)
               VAV(ii,NS) = Vin(ii)
            ENDDO
            DO ii = 1,NP
               UBar(ii) = SUM(UAV(ii,1:NS))/DBLE(NS)
               VBar(ii) = SUM(VAV(ii,1:NS))/DBLE(NS)
            ENDDO
         ELSE
            DO ii = 1,NP
               UAV(ii,ISTA) = Uin(ii)
               VAV(ii,ISTA) = Vin(ii)
               UBar(ii) = SUM(UAV(ii,1:ISTA))/DBLE(ISTA)
               VBar(ii) = SUM(VAV(ii,1:ISTA))/DBLE(ISTA)
            ENDDO
         ENDIF
      ENDIF
      RETURN
      ENDSUBROUTINE UNTIDE
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE ESTIMATE_TIDAL_VELO(Uin,Vin,HA_time,first_solve)
!---------------------------------------------------------------------!
!     Returns an estimated tidal output. If we have performed a       !
!     harmonic analysis already then it estimates the output at       !
!     HA_time. If not, it estimates the tidal output using a lagged   !
!     25-hour average filter.
!---------------------------------------------------------------------!
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: Uin(:),Vin(:)
      REAL(8),INTENT(IN) :: HA_Time
      LOGICAL,INTENT(IN) :: first_solve
      ! local
      INTEGER :: kk, iconst, nconst
      REAL(8) :: cost,sint
      UTide = 0d0
      VTide = 0d0
      nconst = size(freq)
      IF (.NOT.first_solve) THEN
         ! use LA25 filter
         DO kk = 1,NP
            UTide(kk) = Uin(kk) - UBar(kk)
            VTide(kk) = Vin(kk) - VBar(kk)
         ENDDO
         RETURN
      ENDIF
      ! resynthesize
      DO kk = 1,NP
         DO iconst = 1,nconst
            cost = c_i(1,iconst,kk)*cos(freq(iconst)*HA_Time)
            sint = s_i(1,iconst,kk)*sin(freq(iconst)*HA_Time)
            UTide(kk) = UTide(kk) + cost + sint
            cost = c_i(2,iconst,kk)*cos(freq(iconst)*HA_Time)
            sint = s_i(2,iconst,kk)*sin(freq(iconst)*HA_Time)
            VTide(kk) = VTide(kk) + cost + sint
         ENDDO
      ENDDO
      RETURN
!---------------------------------------------------------------------!
      ENDSUBROUTINE ESTIMATE_TIDAL_VELO
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE UPDATE_LSQ_HA(Uin,Vin,HA_Time,flag)
!---------------------------------------------------------------------!
!     Checks if it is time to update the mass matrix and load vector  !
!     for the harmonic analysis. If it is it then updates
!---------------------------------------------------------------------!
      IMPLICIT NONE
      REAL(8),INTENT(IN) :: Uin(:), Vin(:)
      REAL(8),INTENT(IN) :: HA_Time
      LOGICAL,INTENT(IN) :: flag
      INTEGER :: ii, jj ! for loops
      INTEGER :: i1, i2, j1, j2 ! for indices
      REAL(8) :: cosii, cosjj, sinii, sinjj
      ! return if it is not time to update
      IF ( .NOT.flag ) THEN
         RETURN
      ENDIF
      ! update LHS matrix
      DO ii = 1,nconst
         cosii = COS(freq(ii)*HA_Time)
         sinii = SIN(freq(ii)*HA_Time)
         i1 = 2*ii-1
         i2 = i1+1
         DO jj = ii,nconst
            cosjj = COS(freq(jj)*HA_Time)
            sinjj = SIN(freq(jj)*HA_Time)
            j1 = 2*jj-1
            j2 = j1+1
            A(i1,j1) = A(i1,j1) + cosii*cosjj
            A(i1,j2) = A(i1,j2) + cosii*sinjj
            A(i2,j2) = A(i2,j2) + sinii*sinjj
            IF (i2.LE.j1) THEN
               A(i2,j1) = A(i2,j1) + sinii*cosjj
            ENDIF
         ENDDO
         ! load vectors
         DO jj = 1,NP
            b(1,i1,jj) = b(1,i1,jj) + Uin(jj)*cosii
            b(1,i2,jj) = b(1,i2,jj) + Uin(jj)*sinii
            b(2,i1,jj) = b(2,i1,jj) + Vin(jj)*cosii
            b(2,i2,jj) = b(2,i2,jj) + Vin(jj)*sinii
         ENDDO
      ENDDO

      RETURN
!---------------------------------------------------------------------!
      ENDSUBROUTINE UPDATE_LSQ_HA
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE SOLVE_LSQ_HA(flag)
!---------------------------------------------------------------------!
!     Checks if it is time to solve for the harmonic constituents.    !
!     If it is then decomposes mass matrix using LU decomp and solves !
!     for all points of interest. 
!---------------------------------------------------------------------!
      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: flag
      INTEGER :: ii, jj, i1,i2 ! for loops
      INTEGER,ALLOCATABLE :: pivot(:)
      ! return if it is not time to solve
      IF ( .NOT.flag ) THEN
         RETURN
      ENDIF
      ! 
      ! Populate the lower left part of the mass matrix.
      !
      DO jj = 1,2*nconst
         DO ii = jj,2*nconst
            A(ii,jj) = A(jj,ii)
         ENDDO
      ENDDO
      ! 
      ! We will solve using LU Decomposition
      !
      ! Pivot matrix allows for non-diagonal pivot points. This ensures
      ! that the LU decomposition converges.
      ALLOCATE( pivot(2*nconst) )
      ! decompose matrix
      CALL LU_DECOMP(A,pivot)
      !
      ! Solve for each point of interest
      !
      DO jj = 1,NP
         ! u-direction
         CALL LU_SOLVE(A,b(1,:,jj),pivot)
         DO ii = 1,nconst
            i1 = 2*ii-1
            i2 = i1+1
            c_i(1,ii,jj) = b(1,i1,jj)
            s_i(1,ii,jj) = b(1,i2,jj)
         ENDDO
         ! v-direction
         CALL LU_SOLVE(A,b(2,:,jj),pivot)
         DO ii = 1,nconst
            i1 = 2*ii-1
            i2 = i1+1
            c_i(2,ii,jj) = b(2,i1,jj)
            s_i(2,ii,jj) = b(2,i2,jj)
         ENDDO
      ENDDO
      A = 0d0
      b = 0d0
      DEALLOCATE( pivot )
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDSUBROUTINE SOLVE_LSQ_HA
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE LU_DECOMP(A,Pivot)
!---------------------------------------------------------------------!
!     Performs the LU decomposition of mass matrix A such that        !
!     L*U = A. Note that the matrices LU are stored in A such that:   !
!        A = | u_{1,1}, u_{1,2}, u_{1,3}, u_{1,4} |                   !
!            | l_{2,1}, u_{2,2}, u_{2,3}, u_{2,4} |                   !
!            | l_{3,1}, l_{3,2}, u_{3,3}, u_{3,4} |                   !
!            | l_{4,1}, l_{4,2}. l_{4,3}, u_{4,4} |                   !
!     and l_{i,i} = 1                                                 !
!     This decomposition can then be used to solve for every node.    !
!---------------------------------------------------------------------!
      IMPLICIT NONE
      ! 
      ! arguments
      ! 
      ! input/output matrix
      REAL(8),DIMENSION(:,:),INTENT(INOUT) :: A
      ! permutations
      INTEGER,DIMENSION(:),INTENT(OUT) :: Pivot
      !
      ! local variables
      !
      INTEGER :: i,j,k,n
      INTEGER,DIMENSION(:),ALLOCATABLE :: p
      REAL,DIMENSION(:),ALLOCATABLE :: Tmp ! temporary row
      n = size(A,1)
      ALLOCATE( Tmp(n), p(n) )
      !
      DO j = 1,n-1
         p = MAXLOC(ABS(A(J:n,j)))+j-1
         IF (p(1).NE.j) THEN
            Tmp(:) = A(j,:)
            A(j,:) = A(p(1),:)
            A(p(1),:) = Tmp(:)
            Pivot(j) = p(1)
         ELSE
            Pivot(j) = j
         ENDIF
         DO i = j+1,n
            A(i,j) = A(i,j)/A(j,j)
            A(i,j+1:n) = A(i,j+1:n) - A(i,j)*A(j,j+1:n)
         ENDDO
      ENDDO
      DEALLOCATE( Tmp, p )
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      END SUBROUTINE LU_DECOMP
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE LU_SOLVE(A,b,pivot)
!---------------------------------------------------------------------!
!     Solves the system of linear equations Ax=b given:               !
!        - LU decomposition of A (stored compactly)                   !
!        - RHS load vector b                                          !
!        - Vector of pivots if LU decomp did not pivot on diagonal    !
!     Returns:                                                        !
!        - Solution to system stored in b                             !
!     First solves Ly=b by forward substitution and then solves       !
!     Ux=y by backwards substitution.                                 !
!---------------------------------------------------------------------!
      IMPLICIT NONE
      ! 
      ! arguments
      ! 
      ! A-matrix decomposed into LU and stored compactly
      REAL(8),DIMENSION(:,:),INTENT(IN) :: A
      ! vector with stored pivots if they are not the diagonal
      INTEGER,DIMENSION(:),INTENT(IN) :: Pivot
      ! RHS Load vector
      REAL(8),DIMENSION(:),INTENT(INOut) :: b
      !
      ! local variables
      !
      INTEGER :: i,j,n
      REAL(8) :: Tmp ! temporary row
      n = size(A,1)
      !
      DO j = 1,n-1
         IF (Pivot(j).NE.j) THEN
            Tmp = b(j)
            b(j) = b(Pivot(j))
            b(Pivot(j)) = Tmp
         ENDIF
      ENDDO
      !
      ! Forward substitution
      !
      DO i = 2,n
         b(i) = b(i) - DOT_PRODUCT(A(i,1:i-1),b(1:i-1))
      ENDDO
      !
      ! Backward substitution
      !
      b(n) = b(n)/A(n,n)
      DO i = n-1,1,-1
         b(i) = (b(i) - DOT_PRODUCT(A(i,i+1:n),b(i+1:n)))/A(i,i)
      ENDDO
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDSUBROUTINE LU_SOLVE
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      SUBROUTINE CREATE_INPUT_SIGNAL(y,t,freq,c,s,timestep,numsteps,
     & numnodes)
!---------------------------------------------------------------------!
!     Calculates a time series y(t) based on input amplitudes c and s !
! and timestep and number of time steps.                              !
!---------------------------------------------------------------------!
      IMPLICIT NONE
      REAL(8),INTENT(OUT) :: y(:,:), t(:)
      REAL(8),INTENT(IN) :: freq(:), c(:,:), s(:,:), timestep
      INTEGER,INTENT(IN) :: numsteps
      INTEGER,INTENT(IN) :: numnodes
      INTEGER :: it, iconst,kk ! for iterations
      INTEGER :: nconst ! number of input constituents
      REAL(8) :: cost, sint ! makes it more readable
      ! get number of constituents
      nconst = SIZE(freq)
      ! loop through times and calculate signal
      DO it = 1,numsteps
         t(it) = (it-1)*timestep
         DO kk = 1,numnodes
            DO iconst = 1,nconst
               cost = c(iconst,kk)*COS(freq(iconst)*t(it))
               sint = s(iconst,kk)*SIN(freq(iconst)*t(it))
               y(it,kk) = y(it,kk) + cost + sint
            ENDDO
         ENDDO
      ENDDO
      RETURN
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      END SUBROUTINE CREATE_INPUT_SIGNAL
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      ENDPROGRAM SHORT_LSQ_HA_Program
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
