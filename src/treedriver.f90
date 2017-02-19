!
!   Author:  Henry A. Boateng  (boateng@umich.edu)
!   Department of Mathematics
!   University of Michigan, Ann Arbor
!
!   Copyright (c) 2013. The Regents of the University of Michigan.
!   All Rights Reserved.
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM TREEDRIVER
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numparsS,numparsT,order,maxparnodeT
      INTEGER :: shrinkT,treelevelT,iflagT
      REAL(KIND=r8) :: theta 

      INTEGER :: ftype, pot_type
      REAL(KIND=r8) :: kappa, eta, eps, T

! arrays for coordinates and charges and energy of target particles

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xS,yS,zS,qS  !source particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xT,yT,zT     !target particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy,tenergy !exact energy and energy via treecode

! variables for potential energy computation

      REAL(KIND=r8) :: tpeng,dpeng

! variables needed for f90 DATE_AND_TIME intrinsic

      CHARACTER (LEN=25) :: sampin1,sampin2,sampin3,sampout
      REAL(KIND=r8)      :: timedirect,timetree

! variables for error calculations

      REAL(KIND=r8) :: inferr,relinferr
      REAL(KIND=r8) :: n2err,reln2err

! local variables

      INTEGER :: i,err,stat
      CHARACTER (LEN=10) :: c1, c2, c3, c4, c5, c6, c7
      REAL(KIND=r8) :: a1, a2, a3, a4

! EXECUTABLE STATEMENTS
      WRITE(6,*) 'Enter the name of input file 1 (xyzq source data)'
      READ(5,*) sampin1

      WRITE(6,*) 'Enter the type of input file 1 (0 is pqr, 1 is flat)'
      READ(5,*) ftype
      
      WRITE(6,*) 'Enter the name of input file 2 (xyz target positions)'
      READ(5,*) sampin2
      
      WRITE(6,*) 'Enter the name of input file 2 (xyz target positions)'
      READ(5,*) sampin3
      
      !WRITE(6,*) 'Enter the name of output file'
      !READ(5,*) sampout
      sampout = 'out.txt'
      
      WRITE(6,*) 'Enter particle number for source and target'
      READ(5,*) numparsS, numparsT

      WRITE(6,*) 'Enter THETA : ' ! The multipole acceptability criterion
      READ(5,*) theta         

      WRITE(6,*) 'Enter ORDER : ' ! The order of the approximation 
      READ(5,*) order

      WRITE(6,*) 'Enter MAXPARNODE : '  ! maximum number of particles in a leaf
      READ(5,*) maxparnodeT ! maxparnodeT is for leaves of target tree

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
      
      !WRITE(6,*) 'Enter TREELEVEL : '   ! maximum number of levels of tree
      !READ(5,*) treelevelT   ! treelevelT corresponds to the level for target tree
      treelevelT = 5
      
      !WRITE(6,*) 'Enter IFLAG (0-divide tree till number of particles &
      ! in leaf is less or equal to maxparnode)'
      !WRITE(6,*) '            (1-divide tree till maxlevel attained )'
      !READ(5,*) iflagT 
      ! Flag to decide whether tree divides till number of particles in a leaf is at 
      ! most maxparnode or till number of levels in a tree is maxlevel
      iflagT = 0            

      ALLOCATE(xS(numparsS),yS(numparsS),zS(numparsS),qS(numparsS),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for xS, yS, zS and qS! '
          STOP
      END IF

      ALLOCATE(xT(numparsT),yT(numparsT),zT(numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for xT, yT, zT ! '
          STOP
      END IF

      ALLOCATE(tenergy(numparsT),denergy(numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

! Read in coordinates and charges for the sources
      OPEN(unit=82,file=sampin1,status='old',action='read')
 
      WRITE(6,*) ' ' 
      WRITE(6,*) "Reading in sources..."

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

! Read in coordinates for the targets
      OPEN(unit=83,file=sampin2,status='old',action='read')

      WRITE(6,*) ' ' 
      WRITE(6,*) "Reading in targets..."
      DO i=1,numparsT
         READ(83,*) xT(i),yT(i),zT(i)
         IF ( MOD(i, 100000) == 0) THEN
            WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                   char(13), " Reading in target ", i, " of ", numparsT
         END IF
      END DO

      CLOSE(83)

! Read in the values for the exact energy at each target       
      OPEN(unit=84,file=sampin3,status='old',action='read')

      WRITE(6,*) ' ' 
      WRITE(6,*) "Reading in direct energy results..."

      READ(84,13) timedirect
      DO i=1,numparsT
         READ(84,17) denergy(i)
         IF ( MOD(i, 100000) == 0) THEN
            WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                   char(13), " Reading in direct ", i, " of ", numparsT
         END IF
      END DO

      CLOSE (84)


      OPEN(unit=85,file=sampout,status='replace',action='write')

! Calling main subroutine to approximate the energy

      CALL TREECODE(xS,yS,zS,qS,xT,yT,zT,numparsS,numparsT, &
                      tenergy,tpeng,order,theta, &
                      maxparnodeT,timetree, treelevelT,iflagT, &
                      pot_type, kappa, eta, eps, T)


      dpeng=SUM(denergy)

      WRITE(6,*) ' '
      WRITE(6,'("                      Direct time (s): ",ES15.8)') timedirect
      WRITE(6,'("                        Tree time (s): ",ES15.8)') timetree
      WRITE(6,*) ' '
      WRITE(6,'("              Direct Potential Energy: ",ES15.8)') dpeng 
      WRITE(6,'("                Tree Potential Energy: ",ES15.8)') tpeng
      WRITE(6,*) ' '
      WRITE(6,'("    Absolute error in total potential: ",ES15.8)') ABS(tpeng-dpeng) 
      WRITE(6,'("    Relative error in total potential: ",ES15.8)') ABS((tpeng-dpeng)/dpeng) 
      WRITE(6,*) ' '


! compute energy errors

      inferr = MAXVAL(ABS(denergy-tenergy))
      relinferr = inferr / MAXVAL(ABS(denergy))
      n2err = SQRT(DOT_PRODUCT(denergy-tenergy,denergy-tenergy))
      reln2err = n2err / SQRT(DOT_PRODUCT(denergy,denergy))

! output errors to standard out

      WRITE(6,*) ' '
      WRITE(6,'(" Absolute inf norm error in potential: ",ES15.8)') inferr
      WRITE(6,'(" Relative inf norm error in potential: ",ES15.8)') relinferr
      WRITE(6,*) ' ' 
      WRITE(6,'("   Absolute 2 norm error in potential: ",ES15.8)') n2err
      WRITE(6,'("   Relative 2 norm error in potential: ",ES15.8)') reln2err
      WRITE(6,*) ' '
         
      write(85,15)numparsS,numparsT,iflagT,maxparnodeT,treelevelT,order,theta
      write(85,16)ABS((tpeng-dpeng)/dpeng),reln2err,timedirect,timetree

      CLOSE(unit=85)

 13   FORMAT(E24.16)
 15   FORMAT(I8,2X,I8,2X,I3,2X,I4,2X,I3,2X,I3,2X,F12.8)
 16   FORMAT(E24.16,2X,E24.16,2X,F24.16,2X,F24.16) 
 17   FORMAT(14X,E24.16)
100   FORMAT(A7, A5, A5, A7, A6, F8.3, F8.3, F8.3, F8.4, A8, A9) 

      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!

