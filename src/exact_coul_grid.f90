
      PROGRAM TREEDRIVER
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(15,307)

! runtime parameters

      INTEGER :: numparsS, numparsT
      INTEGER :: ftype, pot_type
      REAL(KIND=r8) :: kappa, eta, eps, T

! arrays for coordinates, charge, force calculations 

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xS,yS,zS,qS !source particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy

      REAL(KIND=r8),DIMENSION(6) :: xyzminmax
      INTEGER,DIMENSION(3) :: xyzdim

! variables for potential energy computation

      REAL(KIND=r8) :: dpeng

! variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2 
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec
      CHARACTER (LEN=512) :: sampin1,sampin2,sampout
      REAL(KIND=r8)      :: timedirect

! local variables

      INTEGER :: i,j,err, stat
      CHARACTER (LEN=10) :: c1, c2, c3, c4, c5, c6, c7
      REAL(KIND=r8) :: a1, a2, a3, a4

! EXECUTABLE STATEMENTS
      WRITE(6,*) 'Enter the name of input file 1 (xyzq source data)'
      READ(5,*) sampin1
      WRITE(6,*) 'Enter the type of input file 1 (0 is pqr, 1 is flat)'
      READ(5,*) ftype
      WRITE(6,*) 'Enter the name of output file'
      READ(5,*) sampout
      WRITE(6,*) 'Enter particle number for source'
      READ(5,*) numparsS

      WRITE(6,*) 'Enter kappa value'
      READ(5,*) kappa 
      WRITE(6,*) 'Enter eta value'
      READ(5,*) eta 
      WRITE(6,*) 'Enter eps value'
      READ(5,*) eps 
      WRITE(6,*) 'Enter T value'
      READ(5,*) T 

      WRITE(6,*) 'Enter potential type (0 total corr asym; 1 direct corr asym)'
      READ(5,*) pot_type 

      WRITE(6,*) 'Enter xyzminmax grid limits'
      READ(5,*) xyzminmax(1),xyzminmax(2),xyzminmax(3)
      READ(5,*) xyzminmax(4),xyzminmax(5),xyzminmax(6)

      WRITE(6,*) 'Enter xyzdim grid dimensions'
      READ(5,*) xyzdim(1),xyzdim(2),xyzdim(3)

      numparsT=xyzdim(1)*xyzdim(2)*xyzdim(3)


      ALLOCATE(xS(numparsS),yS(numparsS),zS(numparsS),qS(numparsS),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for xS, yS, zS and qS! '
          STOP
      END IF

      ALLOCATE(denergy(numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

      OPEN(unit=82,file=sampin1,status='old',action='read')
 
      print *, "Reading in sources..."
      IF (ftype == 0) THEN
          i = 1
          DO
              READ(82, 100, iostat=stat) c1, c2, c3, c4, c5, a1, a2, a3, a4, c6, c7
              IF (stat /= 0) EXIT
              IF (c1(1:4) == 'ATOM') THEN
                  WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                        char(13), " Reading in source ", i, " of ", numparsS
                  xS(i) = a1
                  yS(i) = a2
                  zS(i) = a3
                  qS(i) = a4
                  i = i + 1
              END IF
          END DO

      ELSE
          DO i=1,numparsS
              WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                    char(13), " Reading in source ", i, " of ", numparsS
              READ(82,*) xS(i),yS(i),zS(i),qS(i)
!             READ(82,'(F15.10,3F16.10)') xS(i),yS(i),zS(i),qS(i)
          END DO
      END IF
     
      CLOSE(82)

      OPEN(unit=85,file=sampout,status='replace',action='write')

      CALL DATE_AND_TIME(datec,timec,zonec,time1)

      print *
      print *, "Calculating potentials by direct summation..."
      CALL DIRECT_ENG(xS,yS,zS,qS,numparsS,numparsT,denergy,dpeng, &
                      pot_type, kappa, eta, eps, T, xyzminmax, xyzdim)

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,timedirect)

      print *, "Writing direct results to file..."

      WRITE(85,13)timedirect
      DO j=1,numparsT
         WRITE(85,14) j, denergy(j)
      END DO

      CLOSE(unit=85)

      print *
      print *, "Total energy sum: ", dpeng
      print *,   "Total time (s): ", timedirect
      print *

 13   FORMAT(E24.16)
 14   FORMAT(I12,2X,E24.16)
100   FORMAT(A7, A5, A5, A7, A6, F8.3, F8.3, F8.3, F8.4, A8, A9) 

      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!

      SUBROUTINE DIRECT_ENG(xS,yS,zS,qS,numparsS,numparsT,denergy,dpeng, &
                            pot_type, kappa, eta, eps, T, xyzminmax, xyzdim) 

      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      REAL(KIND=r8),PARAMETER :: kb = 0.001987215873_r8
      REAL(KIND=r8),PARAMETER :: coulomb = 332.0637790571_r8

      INTEGER,INTENT(IN) :: numparsS,numparsT 
      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: xS,yS,zS,qS

      REAL(KIND=r8),INTENT(IN) :: kappa, eta, eps, T
      INTEGER,INTENT(IN) :: pot_type

      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzminmax

      REAL(KIND=r8),DIMENSION(numparsT),INTENT(INOUT) :: denergy
      REAL(KIND=r8),INTENT(INOUT) :: dpeng

! local variables 
   
      INTEGER :: i,j,k,nn,m
      REAL(KIND=r8) :: tx,ty,tz,xi,yi,zi,teng, rad
      REAL(KIND=r8),DIMENSION(3) :: xyzdd

      xyzdd(1) = (xyzminmax(2)-xyzminmax(1)) / (xyzdim(1)-1)
      xyzdd(2) = (xyzminmax(4)-xyzminmax(3)) / (xyzdim(2)-1)
      xyzdd(3) = (xyzminmax(6)-xyzminmax(5)) / (xyzdim(3)-1)

      dpeng=0.0_r8; denergy=0.0_r8

      IF (pot_type == 0) THEN
          nn=0
          DO i=0,xyzdim(1)-1
              DO j=0,xyzdim(2)-1
                  DO k=0,xyzdim(3)-1

                      nn=nn+1
                      xi=xyzminmax(1)+i*xyzdd(1)
                      yi=xyzminmax(3)+j*xyzdd(2)
                      zi=xyzminmax(5)+k*xyzdd(3)
                      teng=0.0_r8

                      DO m=1,numparsS
                          tx=xi-xS(m)
                          ty=yi-yS(m)
                          tz=zi-zS(m)
                          rad = SQRT(tx*tx + ty*ty + tz*tz)
                          teng = teng + qS(m) / rad &
                               * (exp(-kappa*rad) * erfc(kappa*eta/2 - rad/eta) &
                               - exp(kappa*rad) * erfc(kappa*eta/2 + rad/eta))
                      END DO

                      denergy(nn) = -teng * exp((kappa*eta)**2 / 4) / (2*eps) &
                                  * sqrt(coulomb/kb/T)
                  END DO
              END DO
          END DO

      ELSE IF (pot_type == 1) THEN
              
          nn=0
          DO i=0,xyzdim(1)-1
              DO j=0,xyzdim(2)-1
                  DO k=0,xyzdim(3)-1

                      nn=nn+1
                      xi=xyzminmax(1)+i*xyzdd(1)
                      yi=xyzminmax(3)+j*xyzdd(2)
                      zi=xyzminmax(5)+k*xyzdd(3)
                      teng=0.0_r8

                      DO m=1,numparsS
                          tx=xi-xS(m)
                          ty=yi-yS(m)
                          tz=zi-zS(m)
                          rad = SQRT(tx*tx + ty*ty + tz*tz)
                          teng = teng + qS(m) / rad * erf(rad/eta)
                      END DO

                      denergy(nn) = -teng * sqrt(coulomb/kb/T)
                  END DO
              END DO
          END DO

      ELSE IF (pot_type == 2) THEN
          nn=0
          DO i=0,xyzdim(1)-1
              DO j=0,xyzdim(2)-1
                  DO k=0,xyzdim(3)-1

                      nn=nn+1
                      xi=xyzminmax(1)+i*xyzdd(1)
                      yi=xyzminmax(3)+j*xyzdd(2)
                      zi=xyzminmax(5)+k*xyzdd(3)
                      teng=0.0_r8

                      DO m=1,numparsS
                          tx=xi-xS(m)
                          ty=yi-yS(m)
                          tz=zi-zS(m)
                          rad = SQRT(tx*tx + ty*ty + tz*tz)
                          teng = teng + qS(m) / rad
                      END DO

                      denergy(nn) = -teng * sqrt(coulomb/kb/T)
                  END DO
              END DO
          END DO
      END IF

      print *, "Summing up energies..."
      dpeng=SUM(denergy)
      
      END SUBROUTINE DIRECT_ENG  

!!!!!!!!!!!!!!!!!!
      SUBROUTINE TTIME(timebeg,timeend,totaltime)
      IMPLICIT NONE
!
! TTIME computes the time difference in seconds between
! the timestamps TIMEBEG and TIMEEND returned by the 
! f90 intrinsic DATE_AND_TIME
!
      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER,DIMENSION(8),INTENT(INOUT) :: timebeg,timeend
      REAL(KIND=r8),INTENT(OUT) :: totaltime

! TIMEEND is modifed by borrowing in case each of its fields
! are not .GE. to the corresponding field in TIMEBEG (up to
! and including days) 

      IF (timeend(8) .LT. timebeg(8)) THEN
          timeend(8)=timeend(8)+1000
          timeend(7)=timeend(7)-1
      END IF
      IF (timeend(7) .LT. timebeg(7)) THEN
          timeend(7)=timeend(7)+60
          timeend(6)=timeend(6)-1
      END IF
      IF (timeend(6) .LT. timebeg(6)) THEN
          timeend(6)=timeend(6)+60
          timeend(5)=timeend(5)-1
      END IF
      IF (timeend(5) .LT. timebeg(5)) THEN
          timeend(5)=timeend(5)+24
          timeend(3)=timeend(3)-1
      END IF

      totaltime=  REAL(timeend(8)-timebeg(8),KIND=r8) +          &
            1000.0_r8*( REAL(timeend(7)-timebeg(7),KIND=r8) +    &
              60.0_r8*( REAL(timeend(6)-timebeg(6),KIND=r8) +    &
              60.0_r8*( REAL(timeend(5)-timebeg(5),KIND=r8) +    &
              24.0_r8*( REAL(timeend(3)-timebeg(3),KIND=r8)))))
      totaltime=totaltime/1000.0_r8


      RETURN
      END SUBROUTINE TTIME
   
