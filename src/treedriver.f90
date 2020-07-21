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
      USE mpi
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(15,307)

! runtime parameters

      INTEGER :: numparsS,numparsT,order,maxparnodeT,numparsTglob
      REAL(KIND=r8) :: theta

      INTEGER :: ftype, pot_type, direct_calc, write_to_file
      REAL(KIND=r8) :: kappa, eta, eps, T

! arrays for coordinates and charges and energy of target particles

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xS,yS,zS,qS  !source particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy,tenergy !exact energy and energy via treecode

      REAL(KIND=r8),DIMENSION(6) :: xyzminmax, xyzminmaxglob
      INTEGER,DIMENSION(3) :: xyzdim, xyzdimglob

! variables for potential energy computation

      REAL(KIND=r8) :: tpeng,dpeng,tpeng_glob,dpeng_glob,ddz

! variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2 
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec

      CHARACTER (LEN=512) :: sampin1,sampin3,sampout
      REAL(KIND=r8)      :: timedirect,timetree,timedirect_max, timetree_max
      REAL(KIND=r8)      :: t1,timetreempi,timetreempi_max

! Fortran MPI IO
      INTEGER :: direct_out, finfo, direct_in
      INTEGER(KIND=MPI_OFFSET_KIND) :: foffset

! variables for error calculations

      REAL(KIND=r8) :: inferr,relinferr, inferr_glob, relinferr_glob
      REAL(KIND=r8) :: n2err,reln2err, n2err_glob, reln2err_glob

! local variables

      INTEGER :: i,err,stat, rank, numProcs, provided
      CHARACTER (LEN=10) :: c1, c2, c3, c4, c5, c6, c7
      REAL(KIND=r8) :: a1, a2, a3, a4, junk_read

      CALL MPI_Init_thread(MPI_THREAD_FUNNELED, provided, err)

      CALL MPI_Comm_size(MPI_COMM_WORLD,numProcs,err)
      CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,err)


! EXECUTABLE STATEMENTS
      !WRITE(6,*) 'Enter the name of input file 1 (xyzq source data)'
      READ(5,*) sampin1

      !WRITE(6,*) 'Enter the type of input file 1 (0 is pqr, 1 is flat)'
      READ(5,*) ftype

      !WRITE(6,*) 'Enter the name of input file 2 (direct in)'
      READ(5,*) sampin3
      
      !WRITE(6,*) 'Enter the name of output file'
      !READ(5,*) sampout
      sampout = 'out.csv'
      
      !WRITE(6,*) 'Enter particle number for source'
      READ(5,*) numparsS

      !WRITE(6,*) 'Enter THETA : ' ! The multipole acceptability criterion
      READ(5,*) theta         

      !WRITE(6,*) 'Enter ORDER : ' ! The order of the approximation 
      READ(5,*) order

      !WRITE(6,*) 'Enter MAXPARNODE : '  ! maximum number of particles in a leaf
      READ(5,*) maxparnodeT ! maxparnodeT is for leaves of target tree

      !WRITE(6,*) 'Enter kappa value'
      READ(5,*) kappa 
      
      !WRITE(6,*) 'Enter eta value'
      READ(5,*) eta 
      
      !WRITE(6,*) 'Enter eps value'
      READ(5,*) eps 
      
      !WRITE(6,*) 'Enter T value'
      READ(5,*) T 
      
      !WRITE(6,*) 'Enter potential type (0 total corr asym; 1 direct corr asym; 2 coulomb)'
      READ(5,*) pot_type

      !WRITE(6,*) 'Enter xyzminmax grid limits'
      READ(5,*) xyzminmaxglob(1),xyzminmaxglob(2),xyzminmaxglob(3)
      READ(5,*) xyzminmaxglob(4),xyzminmaxglob(5),xyzminmaxglob(6)

      !WRITE(6,*) 'Enter xyzdim grid dimensions'
      READ(5,*) xyzdimglob(1),xyzdimglob(2),xyzdimglob(3)

      !WRITE(6,*) 'Direct calculation (0-compute, 1-readin, 2-neither)'
      READ(5,*) direct_calc

      !WRITE(6,*) 'Write to file (0-no, 1-yes)'
      READ(5,*) write_to_file


      numparsTglob=xyzdimglob(1)*xyzdimglob(2)*xyzdimglob(3)
      ddz = (xyzminmaxglob(6)-xyzminmaxglob(5))/REAL(xyzdimglob(3)-1,r8)


      xyzdim(1) = xyzdimglob(1)
      xyzdim(2) = xyzdimglob(2)
      xyzdim(3) = xyzdimglob(3) / numProcs
      
      numparsT=xyzdim(1)*xyzdim(2)*xyzdim(3)

      xyzminmax(1) = xyzminmaxglob(1)
      xyzminmax(2) = xyzminmaxglob(2)
      xyzminmax(3) = xyzminmaxglob(3)
      xyzminmax(4) = xyzminmaxglob(4)

      xyzminmax(5) = xyzminmaxglob(5) + REAL(rank*xyzdim(3),r8)*ddz
      xyzminmax(6) = xyzminmax(5) + REAL(xyzdim(3)-1,r8)*ddz


      ALLOCATE(xS(numparsS),yS(numparsS),zS(numparsS),qS(numparsS),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for xS, yS, zS and qS! '
          STOP
      END IF

      ALLOCATE(tenergy(numparsT),denergy(numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

! Read in coordinates and charges for the sources
      OPEN(unit=82,file=sampin1,status='old',action='read')
 
      !WRITE(6,*) ' ' 
      !WRITE(6,*) "Reading in sources..."

      IF (ftype == 0) THEN
          i = 1
          DO
              READ(82, 100, iostat=stat) c1, c2, c3, c4, c5, a1, a2, a3, a4, c6, c7
              IF (stat /= 0) EXIT
              IF (c1(1:4) == 'ATOM') THEN
!                  WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
!                        char(13), " Reading in source ", i, " of ", numparsS
                  xS(i) = a1
                  yS(i) = a2
                  zS(i) = a3
                  qS(i) = a4
                  i = i + 1
              END IF
          END DO

      ELSE
          DO i=1,numparsS
!              WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
!                    char(13), " Reading in source ", i, " of ", numparsS
              READ(82,*) xS(i),yS(i),zS(i),qS(i)
!             READ(82,'(F15.10,3F16.10)') xS(i),yS(i),zS(i),qS(i)
          END DO
      END IF
     
      CLOSE(82)


! Read in the values for the exact energy at each target       
      IF (direct_calc == 1) THEN 

         IF (rank == 0) THEN
            WRITE(6,*) ' ' 
            WRITE(6,*) "Reading in direct energy results..."
         END IF

         CALL MPI_File_open(MPI_COMM_WORLD, sampin3, MPI_MODE_RDONLY, finfo, direct_in, err) 
         foffset = 8*numparsT*rank
         CALL MPI_File_read_at_all(direct_in, foffset, denergy, numparsT, MPI_DOUBLE, MPI_STATUS_IGNORE, err)
         CALL MPI_File_close(direct_in, err)

         !OPEN(unit=84,file=sampin3,status='old',action='read')
         ! Throw out everything before this rank
         !DO i=1,numparsT*rank
         !   READ(84,13) junk_read
         !END DO

         !DO i=1,numparsT
         !   READ(84,13) denergy(i)
         !END DO
         !CLOSE (84)

         dpeng = SUM(denergy)

      END IF


! Calling main subroutine to approximate the energy

      t1 = MPI_Wtime()
      CALL TREECODE(xS,yS,zS,qS,xyzminmax,xyzdim,numparsS,numparsT, &
                      tenergy,tpeng,order,theta, &
                      maxparnodeT,timetree, &
                      pot_type, kappa, eta, eps, T)
      timetreempi = MPI_Wtime() - t1


      IF (direct_calc == 0) THEN
         CALL DATE_AND_TIME(datec,timec,zonec,time1)
         CALL DIRECT_ENG(xS,yS,zS,qS,numparsS,numparsT,denergy,dpeng, &
                         pot_type, kappa, eta, eps, T, xyzminmax, xyzdim) 
         CALL DATE_AND_TIME(datec,timec,zonec,time2)
         CALL TTIME(time1,time2,timedirect)
      END IF

      CALL MPI_Reduce(dpeng, dpeng_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, err)
      CALL MPI_Reduce(tpeng, tpeng_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, err)

      CALL MPI_Reduce(timetree, timetree_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, err)
      CALL MPI_Reduce(timetreempi, timetreempi_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, err)
      CALL MPI_Reduce(timedirect, timedirect_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, err)


      IF (rank == 0) THEN 
          WRITE(6,*) ' '
          IF (direct_calc == 0) THEN
             WRITE(6,'("                      Direct time (s): ",ES15.8)') timedirect_max
          END IF
          WRITE(6,'("                        Tree time (s): ",ES15.8)') timetree_max
          WRITE(6,*) ' '
          IF (direct_calc == 0 .OR. direct_calc == 1) THEN
             WRITE(6,'("              Direct Potential Energy: ",ES15.8)') dpeng_glob 
          END IF
          WRITE(6,'("                Tree Potential Energy: ",ES15.8)') tpeng_glob
          WRITE(6,*) ' '
          IF (direct_calc == 0 .OR. direct_calc == 1) THEN
             WRITE(6,'("    Absolute error in total potential: ",ES15.8)') ABS(tpeng_glob-dpeng_glob) 
             WRITE(6,'("    Relative error in total potential: ",ES15.8)') ABS((tpeng_glob-dpeng_glob)/dpeng_glob) 
             WRITE(6,*) ' '
          END IF
      END IF


! compute energy errors

      IF (direct_calc == 0 .OR. direct_calc == 1) THEN
          inferr = MAXVAL(ABS(denergy-tenergy))
          relinferr = MAXVAL(ABS(denergy))
          n2err = DOT_PRODUCT(denergy-tenergy,denergy-tenergy)
          reln2err = DOT_PRODUCT(denergy,denergy)
    
          CALL MPI_Reduce(inferr, inferr_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, err)
          CALL MPI_Reduce(relinferr, relinferr_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, err)
          CALL MPI_Reduce(n2err, n2err_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, err)
          CALL MPI_Reduce(reln2err, reln2err_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, err)
    
          relinferr_glob = inferr_glob / relinferr_glob
          n2err_glob = SQRT(ABS(n2err_glob))
          reln2err_glob = n2err_glob / SQRT(ABS(reln2err_glob))

! output errors to standard out

          IF (rank == 0) THEN 
              WRITE(6,*) ' '
              WRITE(6,'(" Absolute inf norm error in potential: ",ES15.8)') inferr_glob
              WRITE(6,'(" Relative inf norm error in potential: ",ES15.8)') relinferr_glob
              WRITE(6,*) ' ' 
              WRITE(6,'("   Absolute 2 norm error in potential: ",ES15.8)') n2err_glob
              WRITE(6,'("   Relative 2 norm error in potential: ",ES15.8)') reln2err_glob
              WRITE(6,*) ' '
                 
              OPEN(unit=85,file=sampout,status='replace',action='write')
              write(85,15) numparsS, ", ", numparsT, ", ", numProcs, ", ", theta, ", ", order, ", ",  maxparnodeT, &
                           ", ", ABS((tpeng_glob-dpeng_glob)/dpeng_glob), ", ", relinferr_glob, ", ", reln2err_glob, &
                           ", ", timetree_max, ", ", timetreempi_max
              CLOSE(unit=85)                                                                         
          END IF
      END IF 

! write to file
      IF (write_to_file == 1 .AND. direct_calc == 0) THEN

          CALL MPI_File_open(MPI_COMM_WORLD, "direct_pot.bin", MPI_MODE_CREATE + MPI_MODE_WRONLY, finfo, direct_out, err) 
          foffset = 8*numparsT*rank
          CALL MPI_File_write_at_all(direct_out, foffset, denergy, numparsT, MPI_DOUBLE, MPI_STATUS_IGNORE, err)
          CALL MPI_File_close(direct_out, err)
          !WRITE(direct_out,17) rank
          !OPEN(unit=85,file=direct_out,status='replace',action='write')
          !DO i=1,numparsT 
          !    WRITE(85,13) denergy(i)
          !END DO
          !CLOSE(85)
      END IF

                                                                                             


 13   FORMAT(E24.16)                                                                         
 15   FORMAT(I9,A,I9,A,I4,A,F12.8,A,I3,A,I4,A,E24.16,A,E24.16,A,E24.16,A,F24.16,A,F24.16)
100   FORMAT(A7, A5, A5, A7, A6, F8.3, F8.3, F8.3, F8.4, A8, A9)  

      CALL MPI_Finalize(err)

      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!

