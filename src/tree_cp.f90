!
!   Authors:  Leighton W. Wilson (lwwilson@umich.edu)
!             Henry A. Boateng  (boateng@umich.edu)
!
!   Department of Mathematics
!   University of Michigan, Ann Arbor
!
!   Copyright (c) 2013-2017. The Regents of the University of Michigan.
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
      MODULE treecode_procedures
      IMPLICIT NONE

! r8 is 8-byte (double precision) real

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(15,307)
      REAL(KIND=r8),PARAMETER :: pi = 4_r8*ATAN(1.D0)

! global variables for taylor expansions

      INTEGER :: torder, torderlim, torder3

! global array for Chebyshev points.
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cheb_pts
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cheb_wts
      
! global variables used when computing potential/force

      REAL(KIND=r8),DIMENSION(3) :: tarpos
      REAL(KIND=r8) :: thetasq,tarposq

! global variables for postition and charge storage

      REAL(KIND=r8),DIMENSION(3) :: xyz_ddglob
      INTEGER,DIMENSION(3) :: xyz_dimglob

! global variables for tree level tracking

      INTEGER :: maxlevel, minlevel

! node pointer and node type declarations

      TYPE tnode_pointer
           TYPE(tnode), POINTER :: p_to_tnode
      END TYPE tnode_pointer
      TYPE tnode
           INTEGER          :: numpar
           REAL(KIND=r8),DIMENSION(3) :: xyz_min, xyz_max, xyz_mid

           REAL(KIND=r8)    :: radius,sqradius,aspect
           INTEGER          :: level,num_children,exist_ms
           REAL(KIND=r8),DIMENSION(:,:),POINTER :: ms
           REAL(KIND=r8),DIMENSION(:,:),POINTER :: tinterp
           TYPE(tnode_pointer), DIMENSION(8) :: child

           INTEGER,DIMENSION(3) :: xyz_dim, xyz_lowindex, xyz_highindex
      END TYPE tnode

      CONTAINS


!!!!!!!!!!!!!!!


      SUBROUTINE SETUP(xyzminmax,xyzdim,xyzind,order,theta)
      IMPLICIT NONE
!
! SETUP allocates and initializes arrays needed for the Taylor expansion.
! Also, global variables are set and the Cartesian coordinates of
! the smallest box containing the particles is determined.
!
      INTEGER,INTENT(IN) :: order
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzminmax
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      INTEGER,DIMENSION(6),INTENT(INOUT) :: xyzind
      REAL(KIND=r8),INTENT(IN) :: theta

! local variables

      INTEGER :: err,i
      REAL(KIND=r8) :: t1, xx

! global integers and reals:  TORDER, TORDERLIM and THETASQ

      torder=order
      torderlim=torder+1
      torder3=torderlim*torderlim*torderlim
      thetasq=theta*theta

      xyz_dimglob = xyzdim

      xyz_ddglob(1) = (xyzminmax(2)-xyzminmax(1)) / (xyz_dimglob(1)-1)
      xyz_ddglob(2) = (xyzminmax(4)-xyzminmax(3)) / (xyz_dimglob(2)-1)
      xyz_ddglob(3) = (xyzminmax(6)-xyzminmax(5)) / (xyz_dimglob(3)-1)

! allocate global Taylor expansion variables

      ALLOCATE(cheb_pts(0:torder), cheb_wts(0:torder), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

      DO i = 0, torder
         xx = i * pi / torder
         cheb_pts(i) = COS(xx)
         cheb_wts(i) = -COS(xx) / (2.0_r8 * SIN(xx) * SIN(xx))
      END DO

      cheb_wts(0) = 0.25_r8 * (torder * torder / 3.0_r8 + 1.0_r8/6.0_r8)
      cheb_wts(torder) = -cheb_wts(0)

! find bounds of Cartesian box enclosing the particles

      xyzind(1) = 0
      xyzind(2) = xyz_dimglob(1)-1
      xyzind(3) = 0
      xyzind(4) = xyz_dimglob(2)-1
      xyzind(5) = 0
      xyzind(6) = xyz_dimglob(3)-1

      RETURN
      END SUBROUTINE SETUP


!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE CREATE_TREE_N0(p,maxparnode,xyzmm,xyzdim,xyzind,level)
      IMPLICIT NONE
!
! CREATE_TREE_N0 recursively creates the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box. The division of a cluster terminates
! when the number of particles in a cluster are is less or equal to maxparnode

      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: level,maxparnode
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      INTEGER,DIMENSION(6),INTENT(IN) :: xyzind

! local variables
      REAL(KIND=r8), DIMENSION(3) :: xyz_len
      REAL(KIND=r8) :: lmax,t2
      INTEGER :: i, err, loclev, numposchild

      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER, DIMENSION(3,8) :: xyzdims
      INTEGER, DIMENSION(6,8) :: xyzinds

      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm
      INTEGER, DIMENSION(3) :: lxyzdim
      INTEGER, DIMENSION(6) :: lxyzind

      xyzmms = 0.0_r8
      xyzdims = 0
      xyzinds = 0
      lxyzmm = 0.0_r8
      lxyzdim = 0
      lxyzind = 0
     
! allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF

! set node fields: number of particles, exist_ms
! and xyz bounds 

      p%numpar=product(xyzdim)
      p%exist_ms=0

      p%xyz_min=xyzmm(1:5:2)
      p%xyz_max=xyzmm(2:6:2)

      p%xyz_dim=xyzdim

      p%xyz_lowindex=xyzind(1:5:2)
      p%xyz_highindex=xyzind(2:6:2)

! compute aspect ratio

      xyz_len=p%xyz_max-p%xyz_min

      lmax=MAXVAL(xyz_len)
      t2=MINVAL(xyz_len)

      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=lmax/t2
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS 

      p%xyz_mid=(p%xyz_max+p%xyz_min)/2.0_r8
      p%sqradius=SUM(xyz_len**2)/4.0_r8
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%level=level
      p%num_children=0

      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF (p%numpar .GT. maxparnode) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(:,1)=xyzmm
         xyzdims(:,1)=xyzdim
         xyzinds(:,1)=xyzind

         CALL PARTITION_8(xyzmms,xyzdims,xyzinds,xyz_len,lmax,numposchild)
!
! create children if indicated and store info in parent
!
         loclev=level+1

         DO i=1,numposchild
            IF (((xyzinds(1,i) .LE. xyzinds(2,i)) .AND. &
                 (xyzinds(3,i) .LE. xyzinds(4,i))) .AND. &
                 (xyzinds(5,i) .LE. xyzinds(6,i))) THEN

               p%num_children=p%num_children+1

               lxyzmm=xyzmms(:,i)
               lxyzdim=xyzdims(:,i)
               lxyzind=xyzinds(:,i)

               CALL CREATE_TREE_N0(p%child(p%num_children)%p_to_tnode, &
                                   maxparnode,lxyzmm,lxyzdim,lxyzind,loclev)
            END IF
         END DO

      END IF   

      END SUBROUTINE CREATE_TREE_N0


!!!!!!!!!!!!!!!


      SUBROUTINE PARTITION_8(xyzmms,xyzdims,xyzinds,xyz_len,lmax,numposchild)

      IMPLICIT NONE
!
! PARTITION_8 determines the particle indices of the eight sub boxes
! containing the particles after the box defined by particles I_BEG
! to I_END is divided by its midpoints in each coordinate direction.
! The determination of the indices is accomplished by the subroutine
! PARTITION. A box is divided in a coordinate direction as long as the
! resulting aspect ratio is not too large. This avoids the creation of
! "narrow" boxes in which Talyor expansions may become inefficient.
! On exit the INTEGER array IND (dimension 8 x 2) contains
! the indice limits of each new box (node) and NUMPOSCHILD the number 
! of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
! that box J is empty.
!
      REAL(KIND=r8),DIMENSION(6,8),INTENT(INOUT) :: xyzmms
      INTEGER,DIMENSION(3,8),INTENT(INOUT) :: xyzdims
      INTEGER,DIMENSION(6,8),INTENT(INOUT) :: xyzinds

      REAL(KIND=r8),DIMENSION(3),INTENT(IN) :: xyz_len
      REAL(KIND=r8),INTENT(IN) :: lmax
      INTEGER,INTENT(INOUT) :: numposchild

! local variables

      INTEGER :: i
      REAL(KIND=r8) :: critlen
      INTEGER :: xdim,ydim,zdim,xn,yn,zn
      INTEGER :: xlowind,xhighind,ylowind,yhighind,zlowind,zhighind
      REAL(KIND=r8) :: xlowmid,xhighmid,ylowmid,yhighmid,zlowmid,zhighmid


      numposchild=1
      critlen=lmax/sqrt(2.0_r8)

      xdim=xyzdims(1,1)
      ydim=xyzdims(2,1)
      zdim=xyzdims(3,1)

      xn=xdim/2
      yn=ydim/2
      zn=zdim/2

      IF (xyz_len(1) .GE. critlen) THEN

         xlowmid=xyzmms(1,1)+(xn-1)*xyz_ddglob(1)
         xhighmid=xyzmms(2,1)-(xdim-xn-1)*xyz_ddglob(1)

         xlowind=xyzinds(1,1)+(xn-1)
         xhighind=xyzinds(2,1)-(xdim-xn-1)

         xyzmms(:,2)=xyzmms(:,1)
         xyzinds(:,2)=xyzinds(:,1)

         xyzmms(2,1)=xlowmid
         xyzmms(1,2)=xhighmid

         xyzinds(2,1)=xlowind
         xyzinds(1,2)=xhighind

         numposchild=2*numposchild

      END IF 

      IF (xyz_len(2) .GE. critlen) THEN

         ylowmid=xyzmms(3,1)+(yn-1)*xyz_ddglob(2)
         yhighmid=xyzmms(4,1)-(ydim-yn-1)*xyz_ddglob(2)

         ylowind=xyzinds(3,1)+(yn-1)
         yhighind=xyzinds(4,1)-(ydim-yn-1)

         DO i=1,numposchild
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzinds(:,numposchild+i)=xyzinds(:,i)

            xyzmms(4,i)=ylowmid
            xyzmms(3,numposchild+i)=yhighmid

            xyzinds(4,i)=ylowind
            xyzinds(3,numposchild+i)=yhighind
         END DO

         numposchild=2*numposchild

      END IF

      IF (xyz_len(3) .GE. critlen) THEN

         zlowmid=xyzmms(5,1)+(zn-1)*xyz_ddglob(3)
         zhighmid=xyzmms(6,1)-(zdim-zn-1)*xyz_ddglob(3)

         zlowind=xyzinds(5,1)+(zn-1)
         zhighind=xyzinds(6,1)-(zdim-zn-1)

         DO i=1,numposchild

            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzinds(:,numposchild+i)=xyzinds(:,i)

            xyzmms(6,i)=zlowmid
            xyzmms(5,numposchild+i)=zhighmid

            xyzinds(6,i)=zlowind
            xyzinds(5,numposchild+i)=zhighind

         END DO

         numposchild=2*numposchild

      END IF

      xyzdims(1,:)=xyzinds(2,:)-xyzinds(1,:)+1
      xyzdims(2,:)=xyzinds(4,:)-xyzinds(3,:)+1
      xyzdims(3,:)=xyzinds(6,:)-xyzinds(5,:)+1

      RETURN 
      END SUBROUTINE PARTITION_8


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE CP_TREECODE_TCF(p, xS,yS,zS,qS, EnP, &
                                 numparsS, numparsT, kappa, eta, eps)
      IMPLICIT NONE
!
! CP_TREECODE is the driver routine which calls COMPUTE_CP1 for each
! source particle, setting the global variable TARPOS before the call. After
! the calls to COMPUTE_CP1, CP_TREECODE calls COMPUTE_CP2 to compute
! the energy using power series. 
! P is the root node of the target tree, xT,yT,zT are the coordinates of
! the target particles while xS,yS,zS,qS are the coordinates and charges of the
! source particles. The energy at target 'i', is stored in EnP(i).
!
 
      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p  
      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: xS,yS,zS,qS
      REAL(KIND=r8),INTENT(IN) :: kappa, eta, eps
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(INOUT) :: EnP
 
! local variables

      INTEGER :: i,j

      EnP=0.0_r8
      DO i=1,numparsS
         tarpos(1)=xS(i)
         tarpos(2)=yS(i)
         tarpos(3)=zS(i)
         tarposq=qS(i)

         DO j=1,p%num_children
            CALL COMPUTE_CP1_TCF(p%child(j)%p_to_tnode,EnP, &
                                 numparsT, kappa, eta)
         END DO
      END DO

      CALL COMPUTE_CP2(p,EnP,numparsT)
 
      EnP = EnP * exp((kappa*eta)**2 / 4_r8) / (2_r8*eps)

      RETURN
      END SUBROUTINE CP_TREECODE_TCF


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE CP_TREECODE_DCF(p, xS,yS,zS,qS, EnP, &
                                 numparsS, numparsT, eta)
      IMPLICIT NONE
!
! CP_TREECODE is the driver routine which calls COMPUTE_CP1 for each
! source particle, setting the global variable TARPOS before the call. After
! the calls to COMPUTE_CP1, CP_TREECODE calls COMPUTE_CP2 to compute
! the energy using power series. 
! P is the root node of the target tree, xT,yT,zT are the coordinates of
! the target particles while xS,yS,zS,qS are the coordinates and charges of the
! source particles. The energy at target 'i', is stored in EnP(i).
!
 
      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p  
      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: xS,yS,zS,qS
      REAL(KIND=r8),INTENT(IN) :: eta
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(INOUT) :: EnP
 
! local variables

      INTEGER :: i,j

      EnP=0.0_r8
      DO i=1,numparsS
         tarpos(1)=xS(i)
         tarpos(2)=yS(i)
         tarpos(3)=zS(i)
         tarposq=qS(i)

         DO j=1,p%num_children
            CALL COMPUTE_CP1_DCF(p%child(j)%p_to_tnode,EnP, &
                                 numparsT, eta)
         END DO
      END DO

      CALL COMPUTE_CP2(p,EnP,numparsT)

      RETURN
      END SUBROUTINE CP_TREECODE_DCF


!!!!!!!!!!!!!!


      SUBROUTINE CP_TREECODE_COULOMB(p, xS,yS,zS,qS, EnP, &
                                     numparsS, numparsT)
      IMPLICIT NONE
!
! CP_TREECODE is the driver routine which calls COMPUTE_CP1 for each
! source particle, setting the global variable TARPOS before the call. After
! the calls to COMPUTE_CP1, CP_TREECODE calls COMPUTE_CP2 to compute
! the energy using power series. 
! P is the root node of the target tree, xT,yT,zT are the coordinates of
! the target particles while xS,yS,zS,qS are the coordinates and charges of the
! source particles. The energy at target 'i', is stored in EnP(i).
!
 
      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p  
      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: xS,yS,zS,qS
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(INOUT) :: EnP
 
! local variables

      INTEGER :: i,j

      EnP=0.0_r8
      DO i=1,numparsS
         tarpos(1)=xS(i)
         tarpos(2)=yS(i)
         tarpos(3)=zS(i)
         tarposq=qS(i)

         DO j=1,p%num_children
            CALL COMPUTE_CP1_COULOMB(p%child(j)%p_to_tnode,EnP, &
                                     numparsT)
         END DO
      END DO

      CALL COMPUTE_CP2(p,EnP,numparsT)

      RETURN
      END SUBROUTINE CP_TREECODE_COULOMB


!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP1_TCF(p, EnP, arrdim, kappa, eta)

      IMPLICIT NONE

! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p      
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP
      REAL(KIND=r8),INTENT(IN) :: kappa, eta

! local variables

      REAL(KIND=r8),DIMENSION(3) :: xyz_t
      REAL(KIND=r8) :: distsq
      INTEGER :: i, err

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8) .AND. (p%numpar > torderlim*torderlim*torderlim)) THEN

         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(8,torderlim*torderlim*torderlim),STAT=err)
             ALLOCATE(p%tinterp(3,0:torder),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF

             CALL COMP_INTERP(p)

             p%ms=0.0_r8
             p%exist_ms=1
         END IF

         CALL COMP_CMS_TCF(p, kappa, eta)
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_TCF(p, EnP, arrdim, kappa, eta)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1_TCF(p%child(i)%p_to_tnode, EnP, arrdim, kappa, eta)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1_TCF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP1_DCF(p,EnP,arrdim, eta)

      IMPLICIT NONE

! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p      
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP
      REAL(KIND=r8),INTENT(IN) :: eta

! local variables

      REAL(KIND=r8),DIMENSION(3) :: xyz_t
      REAL(KIND=r8) :: distsq
      INTEGER :: i, err

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8) .AND. (p%numpar > torderlim*torderlim*torderlim)) THEN

         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(8,torderlim*torderlim*torderlim),STAT=err)
             ALLOCATE(p%tinterp(3,0:torder),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF

             CALL COMP_INTERP(p)

             p%ms=0.0_r8
             p%exist_ms=1
         END IF

         CALL COMP_CMS_DCF(p, eta)
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_DCF(p, EnP, arrdim, eta)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1_DCF(p%child(i)%p_to_tnode, EnP, arrdim, eta)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1_DCF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP1_COULOMB(p,EnP,arrdim)

      IMPLICIT NONE

! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p      
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP

! local variables

      REAL(KIND=r8),DIMENSION(3) :: xyz_t
      REAL(KIND=r8) :: distsq
      INTEGER :: i, err

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8) .AND. (p%numpar > torderlim*torderlim*torderlim)) THEN

         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(8,torderlim*torderlim*torderlim),STAT=err)
             ALLOCATE(p%tinterp(3,0:torder),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF

             CALL COMP_INTERP(p)

             p%ms=0.0_r8
             p%exist_ms=1
         END IF

         CALL COMP_CMS_COULOMB(p)
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_COULOMB(p, EnP, arrdim)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1_COULOMB(p%child(i)%p_to_tnode, EnP, arrdim)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1_COULOMB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_CP2(ap,EnP,arrdim)

        IMPLICIT NONE

! COMPUTE_CP2 is a recursive routine that evaluates the power series  
! approximation of the potential at the targets in a cluster via 
! a 3-D Horner's rule.  
!
        INTEGER,INTENT(IN) :: arrdim
        TYPE(tnode),POINTER,INTENT(IN) :: ap
        REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP

! local variables

        REAL(KIND=r8) :: peng,xl,yl,zl,dx,dy,dz,temp11,temp12,temp21,temp22
        INTEGER :: xlind, ylind, zlind, xhind, yhind, zhind, yzhind
        INTEGER :: i,nn,j,k,k1,k2,k3,kk

        REAL(KIND=r8), DIMENSION(0:torder) :: dj, wx, wy, wz
        REAL(KIND=r8), DIMENSION(0:torder) :: a1i, a2j, a3k, b1i, b2j, b3k 
        REAL(KIND=r8), DIMENSION(torder3) :: sum1, sum2, sum3, sum4 
        REAL(KIND=r8), DIMENSION(torder3) :: sum5, sum6, sum7, sum8 
        REAL(KIND=r8) :: sumA1, sumA2, sumA3, xx, yy, zz, Dd
        INTEGER :: a1exactind, a2exactind, a3exactind

        IF (ap%exist_ms==1) THEN

          xl=ap%xyz_min(1)
          yl=ap%xyz_min(2)
          zl=ap%xyz_min(3)

          xlind=ap%xyz_lowindex(1)
          ylind=ap%xyz_lowindex(2)
          zlind=ap%xyz_lowindex(3)

          xhind=ap%xyz_highindex(1)
          yhind=ap%xyz_highindex(2)
          zhind=ap%xyz_highindex(3)

          sum1 = 0.0_r8;
          sum2 = 0.0_r8;
          sum3 = 0.0_r8;
          sum4 = 0.0_r8;
          sum5 = 0.0_r8;
          sum6 = 0.0_r8;
          sum7 = 0.0_r8;
          sum8 = 0.0_r8;

          dj = 1.0_r8
          dj(0) = 0.25_r8
          dj(torder) = 0.25_r8

          DO k1=0,torder
             wx(k1) = -4.0 * cheb_wts(k1) / (ap%xyz_max(1) - ap%xyz_min(1))
             wy(k1) = -4.0 * cheb_wts(k1) / (ap%xyz_max(2) - ap%xyz_min(2))
             wz(k1) = -4.0 * cheb_wts(k1) / (ap%xyz_max(3) - ap%xyz_min(3))
          END DO

          yzhind=xyz_dimglob(2)*xyz_dimglob(3)
          DO i=xlind,xhind
             xx=xl+(i-xlind)*xyz_ddglob(1)
             DO j=ylind,yhind
                yy=yl+(j-ylind)*xyz_ddglob(2)
                DO k=zlind,zhind

                   sumA1 = 0.0_r8
                   sumA2 = 0.0_r8
                   sumA3 = 0.0_r8

                   a1exactind = -1
                   a2exactind = -1
                   a3exactind = -1

                   peng = 0.0_r8
                   kk = 0

                   nn=(i*yzhind)+(j*xyz_dimglob(3))+k+1
                   zz=zl+(k-zlind)*xyz_ddglob(3)

                   DO k1 = 0, torder
                   
                      dx = xx - ap%tinterp(1,k1)
                      dy = yy - ap%tinterp(2,k1)
                      dz = zz - ap%tinterp(3,k1)

                      a1i(k1) = wx(k1)/dx + dj(k1)/(dx*dx)
                      a2j(k1) = wy(k1)/dy + dj(k1)/(dy*dy)
                      a3k(k1) = wz(k1)/dz + dj(k1)/(dz*dz)

                      b1i(k1) = dj(k1)/dx
                      b2j(k1) = dj(k1)/dy
                      b3k(k1) = dj(k1)/dz

                      sumA1 = sumA1 + a1i(k1)
                      sumA2 = sumA2 + a2j(k1)
                      sumA3 = sumA3 + a3k(k1)

                      IF (ABS(xx - ap%tinterp(1,k1)) < TINY(1.0_r8)) THEN
                         a1exactind = k1
                      END IF

                      IF (ABS(yy - ap%tinterp(2,k1)) < TINY(1.0_r8)) THEN
                         a2exactind = k1
                      END IF

                      IF (ABS(zz - ap%tinterp(3,k1)) < TINY(1.0_r8)) THEN
                         a3exactind = k1
                      END IF
                   END DO

                   IF (a1exactind > -1) THEN
                      sumA1 = 1.0_r8
                      a1i = 0.0_r8
                      b1i = 0.0_r8
                      a1i(a1exactind) = 1.0_r8
                   END IF

                   IF (a2exactind > -1) THEN
                      sumA2 = 1.0_r8
                      a2j = 0.0_r8
                      b2j = 0.0_r8
                      a2j(a2exactind) = 1.0_r8
                   END IF

                   IF (a3exactind > -1) THEN
                      sumA3 = 1.0_r8
                      a3k = 0.0_r8
                      b3k = 0.0_r8
                      a3k(a3exactind) = 1.0_r8
                   END IF

                   Dd = 1.0_r8 / (sumA1 * sumA2 * sumA3);

                   DO k1 = 0, torder
                      DO k2 = 0, torder
                         DO k3 = 0, torder
                            kk = kk + 1
                            temp11 = a1i(k1) * a2j(k2) * Dd
                            temp21 = b1i(k1) * a2j(k2) * Dd
                            temp12 = a1i(k1) * b2j(k2) * Dd
                            temp22 = b1i(k1) * b2j(k2) * Dd
                            peng = peng + ap%ms(1,kk) * temp11 * a3k(k3) &
                                        + ap%ms(2,kk) * temp21 * a3k(k3) &
                                        + ap%ms(3,kk) * temp12 * a3k(k3) &
                                        + ap%ms(4,kk) * temp11 * b3k(k3) &
                                        + ap%ms(5,kk) * temp22 * a3k(k3) &
                                        + ap%ms(6,kk) * temp12 * b3k(k3) &
                                        + ap%ms(7,kk) * temp21 * b3k(k3) &
                                        + ap%ms(8,kk) * temp22 * b3k(k3)
                         END DO
                      END DO
                   END DO

                   EnP(nn)=EnP(nn)+peng

                END DO
             END DO
          END DO

        END IF

        DO j=1,ap%num_children
           CALL COMPUTE_CP2(ap%child(j)%p_to_tnode,EnP,arrdim)
        END DO

        RETURN
      END SUBROUTINE COMPUTE_CP2


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_INTERP(p)
      IMPLICIT NONE

      TYPE(tnode),POINTER :: p
      INTEGER :: i

      DO i = 0, torder
          p%tinterp(:,i) = (p%xyz_min) + (cheb_pts(i) + 1.0_r8)/2.0_r8 &
                         * (p%xyz_max - p%xyz_min)
      END DO

      END SUBROUTINE COMP_INTERP


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_CMS_TCF(p, kappa, eta)
      IMPLICIT NONE
!
! COMP_CMS computes the moments for node P needed in the Taylor approximation
!
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),INTENT(IN)  :: kappa, eta

! local variables

      INTEGER :: k1,k2,k3,kk
      REAL(KIND=r8), DIMENSION(0:torder) :: dx, dy, dz, dx2, dy2, dz2
      REAL(KIND=r8) :: rad, invR, invRk, invR2, invR4, pot, d1, dpot1, dpot2, dpot3
      REAL(KIND=r8) :: kappa2

      kappa2 = kappa * kappa

      kk=0

      DO k1=0,torder
          dx(k1) = (tarpos(1) - p%tinterp(1,k1))
          dy(k1) = (tarpos(2) - p%tinterp(2,k1))
          dz(k1) = (tarpos(3) - p%tinterp(3,k1))
          dx2(k1) = dx(k1)**2
          dy2(k1) = dy(k1)**2
          dz2(k1) = dz(k1)**2
      END DO

      DO k1=0,torder
         DO k2=0,torder
            DO k3=0,torder
               kk=kk+1
               rad = SQRT(dx2(k1) + dy2(k2) + dz2(k3))
               invR = 1.0_r8 / rad
               invRk = kappa * invR
               invR2 = invR * invR
               invR4 = invR2 * invR2

               pot = 2.0_r8 * tarposq * invR * exp(-kappa * rad)

               d1 = invRk + invR2
               dpot1 = d1 * pot
               dpot2 = (invR2 * (kappa2 + 3.0_r8 * d1)) * pot
               dpot3 = (invRk * invRk * (invRk + 6.0_r8 * invR2) + 15 * invR4 * d1) * pot

               p%ms(1,kk) = p%ms(1,kk) - pot
               p%ms(2,kk) = p%ms(2,kk) - dpot1 * dx(k1)
               p%ms(3,kk) = p%ms(3,kk) - dpot1 * dy(k2)
               p%ms(4,kk) = p%ms(4,kk) - dpot1 * dz(k3)
               p%ms(5,kk) = p%ms(5,kk) - dpot2 * dx(k1) * dy(k2)
               p%ms(6,kk) = p%ms(6,kk) - dpot2 * dy(k2) * dz(k3)
               p%ms(7,kk) = p%ms(7,kk) - dpot2 * dx(k1) * dz(k3)
               p%ms(8,kk) = p%ms(8,kk) - dpot3 * dx(k1) * dy(k2) * dz(k3)

           END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE COMP_CMS_TCF


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_CMS_DCF(p, eta)
      IMPLICIT NONE
!
! COMP_CMS computes the moments for node P needed in the Taylor approximation
!
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),INTENT(IN)  :: eta

! local variables

      INTEGER :: k1,k2,k3,kk
      REAL(KIND=r8), DIMENSION(0:torder) :: dx, dy, dz, dx2, dy2, dz2
      REAL(KIND=r8) :: invR, invR2, invR4, pot, auxpot, auxpot_eta
      REAL(KIND=r8) :: rad, rad_eta, dpot1, dpot2, dpot3
      REAL(KIND=r8) :: twoinvEta2, teninvEta2, fourinvEta4

      twoinvEta2 = 2.0_r8 / (eta * eta)
      teninvEta2 = 5.0_r8 * twoinvEta2
      fourinvEta4 = twoinvEta2 * twoinvEta2

      kk=0

      DO k1=0,torder
          dx(k1) = (tarpos(1) - p%tinterp(1,k1))
          dy(k1) = (tarpos(2) - p%tinterp(2,k1))
          dz(k1) = (tarpos(3) - p%tinterp(3,k1))
          dx2(k1) = dx(k1)**2
          dy2(k1) = dy(k1)**2
          dz2(k1) = dz(k1)**2
      END DO

      DO k1=0,torder
         DO k2=0,torder
            DO k3=0,torder
               kk=kk+1

               rad = SQRT(dx2(k1) + dy2(k2) + dz2(k3))
               rad_eta = rad / eta
               invR = 1.0_r8 / rad
               invR2 = invR * invR
               invR4 = invR2 * invR2

               pot = tarposq * erf(rad_eta) * invR
               auxpot = tarposq * 2.0_r8 / sqrt(pi) * exp(-rad_eta * rad_eta)
               auxpot_eta = auxpot / eta

               dpot1 = invR2 * (pot - auxpot)
               dpot2 = invR2 * (3.0_r8 * pot * invR2 - auxpot_eta &
                          * (3.0_r8 * invR2 + twoinvEta2))
               dpot3 = invR2 * (15.0_r8 * pot * invR4 - auxpot_eta &
                          * (15.0_r8 * invR4 + teninvEta2 * invR2 + fourinvEta4))

               p%ms(1,kk) = p%ms(1,kk) - pot
               p%ms(2,kk) = p%ms(2,kk) - dpot1 * dx(k1)
               p%ms(3,kk) = p%ms(3,kk) - dpot1 * dy(k2)
               p%ms(4,kk) = p%ms(4,kk) - dpot1 * dz(k3)
               p%ms(5,kk) = p%ms(5,kk) - dpot2 * dx(k1) * dy(k2)
               p%ms(6,kk) = p%ms(6,kk) - dpot2 * dy(k2) * dz(k3)
               p%ms(7,kk) = p%ms(7,kk) - dpot2 * dx(k1) * dz(k3)
               p%ms(8,kk) = p%ms(8,kk) - dpot3 * dx(k1) * dy(k2) * dz(k3)
           END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE COMP_CMS_DCF


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_CMS_COULOMB(p)
      IMPLICIT NONE
!
! COMP_CMS computes the moments for node P needed in the Taylor approximation
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: k1,k2,k3,kk
      REAL(KIND=r8), DIMENSION(0:torder) :: dx, dy, dz, dx2, dy2, dz2
      REAL(KIND=r8) :: invR, invRq, invR2, invR3, invR5, invR7

      kk=0

      DO k1=0,torder
          dx(k1) = (tarpos(1) - p%tinterp(1,k1))
          dy(k1) = (tarpos(2) - p%tinterp(2,k1))
          dz(k1) = (tarpos(3) - p%tinterp(3,k1))
          dx2(k1) = dx(k1)**2
          dy2(k1) = dy(k1)**2
          dz2(k1) = dz(k1)**2
      END DO


      DO k1=0,torder
         DO k2=0,torder
            DO k3=0,torder
               kk=kk+1
               invR = 1.0_r8 / SQRT(dx2(k1) + dy2(k2) + dz2(k3)) 
               invRq = invR * tarposq
               invR2 = invR * invR
               invR3 = invRq * invR2
               invR5 = invR3 * invR2
               invR7 = invR5 * invR2

               invR5 = invR5 * 3.0_r8

               p%ms(1,kk) = p%ms(1,kk) - invRq
               p%ms(2,kk) = p%ms(2,kk) - invR3 * dx(k1)
               p%ms(3,kk) = p%ms(3,kk) - invR3 * dy(k2)
               p%ms(4,kk) = p%ms(4,kk) - invR3 * dz(k3)
               p%ms(5,kk) = p%ms(5,kk) - invR5 * dx(k1) * dy(k2)
               p%ms(6,kk) = p%ms(6,kk) - invR5 * dy(k2) * dz(k3)
               p%ms(7,kk) = p%ms(7,kk) - invR5 * dx(k1) * dz(k3)
               p%ms(8,kk) = p%ms(8,kk) - invR7 * dx(k1) * dy(k2) * dz(k3) * 15.0_r8
           END DO
         END DO
      END DO
         
      RETURN
      END SUBROUTINE COMP_CMS_COULOMB


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT_TCF(p, EnP, arrdim, kappa, eta)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP
      REAL(KIND=r8),INTENT(IN) :: kappa, eta

! local variables

      INTEGER :: i,j,k,nn
      REAL(KIND=r8) :: tx,ty,tz,xl,yl,zl,rad
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,yzhind
      REAL(KIND=r8) :: kap_eta_2, kap_rad, rad_eta

      kap_eta_2 = kappa * eta / 2_r8

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      yzhind=xyz_dimglob(2)*xyz_dimglob(3)
      DO i=xlind,xhind
         tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)
         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            DO k=zlind,zhind

               tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)
               nn=(i*yzhind)+(j*xyz_dimglob(3))+k+1

               rad = SQRT(tx*tx + ty*ty + tz*tz)
               kap_rad = kappa * rad
               rad_eta = rad / eta

               EnP(nn) = EnP(nn) - tarposq / rad &
                                 * (exp(-kap_rad) * erfc(kap_eta_2 - rad_eta) &
                                 -  exp( kap_rad) * erfc(kap_eta_2 + rad_eta))

            END DO
         END DO
      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT_TCF


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT_DCF(p, EnP, arrdim, eta)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP
      REAL(KIND=r8),INTENT(IN) :: eta

! local variables

      INTEGER :: i,j,k,nn
      REAL(KIND=r8) :: tx,ty,tz,xl,yl,zl,rad
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,yzhind
      REAL(KIND=r8) :: rad_eta

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      yzhind=xyz_dimglob(2)*xyz_dimglob(3)
      DO i=xlind,xhind
         tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)
         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            DO k=zlind,zhind

               nn=(i*yzhind)+(j*xyz_dimglob(3))+k+1
               tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)

               rad = SQRT(tx*tx + ty*ty + tz*tz)
               rad_eta = rad / eta

               EnP(nn) = EnP(nn) - tarposq / rad * erf(rad_eta)

            END DO
         END DO
      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT_DCF


!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT_COULOMB(p, EnP, arrdim)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP

! local variables

      INTEGER :: i,j,k,nn
      REAL(KIND=r8) :: tx,ty,tz,xl,yl,zl
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,yzhind

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      yzhind=xyz_dimglob(2)*xyz_dimglob(3)
      DO i=xlind,xhind
         tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)
         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            DO k=zlind,zhind

               nn=(i*yzhind)+(j*xyz_dimglob(3))+k+1
               tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)

               EnP(nn) = EnP(nn) - tarposq / SQRT(tx*tx + ty*ty + tz*tz)

            END DO
         END DO
      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT_COULOMB


!!!!!!!!!!!!!!!!!


      SUBROUTINE CLEANUP(p)
      IMPLICIT NONE
!
! CLEANUP deallocates allocated global variables and then
! calls recursive routine REMOVE_NODE to delete the tree.
!
      TYPE(tnode),POINTER :: p      

! local variables
  
      INTEGER :: err

      CALL REMOVE_NODE(p)
      DEALLOCATE(p, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating root node! '
         STOP
      END IF 
      NULLIFY(p)         

      RETURN
      END SUBROUTINE CLEANUP


!!!!!!!!!!!


      RECURSIVE SUBROUTINE REMOVE_NODE(p)
      IMPLICIT NONE
!
! REMOVE_NODE recursively removes each node from the
! tree and deallocates its memory for MS array if it
! exits.
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: i,err

      IF (p%exist_ms .EQ. 1) THEN
         DEALLOCATE(p%ms,STAT=err)
         IF (err .NE. 0) THEN
            WRITE(6,*) 'Error deallocating node MS! '
            STOP
         END IF               
      END IF

      IF (p%num_children .GT. 0) THEN
          DO i=1,p%num_children
            CALL REMOVE_NODE(p%child(i)%p_to_tnode)
            DEALLOCATE(p%child(i)%p_to_tnode,STAT=err)
            IF (err .NE. 0) THEN
               WRITE(6,*) 'Error deallocating node child! '
               STOP
            END IF                           
          END DO
      END IF 

      RETURN                
      END SUBROUTINE REMOVE_NODE      

      END MODULE treecode_procedures

!!!!!!!!!!!!!!!

      SUBROUTINE TREECODE(xS,yS,zS,qS,xyzminmax,xyzdim,numparsS,numparsT, &
                      tEn,tpeng,order,theta, &
                      maxparnodeT,timetree, &
                      pot_type, kappa, eta, eps, T)

      USE treecode_procedures
      IMPLICIT NONE

!====================================================================
!                                                                   
! xS,yS,zS,qS    :: x,y,z coordinates and charges of sources
! xyzdim         :: target grid dimensions
! xyzminmax      :: min, max target grid limits
! numparS        :: number of sources
! numparT        :: number of sources
! tEn            :: array of dimension numparT for storing potential
!                   at each target
! tpeng          :: total potential
! maxparnodeT    :: maximum number of particles in a leaf
!                   (only employed in the cluster-particle version)
! timetree       :: The total time for the treecode computation
!=====================================================================

      INTEGER,INTENT(IN) :: numparsS,numparsT,order,maxparnodeT
      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: xS,yS,zS,qS

      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzminmax
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim

      REAL(KIND=r8),DIMENSION(numparsT),INTENT(OUT) :: tEn
      REAL(KIND=r8),INTENT(IN) :: theta
      REAL(KIND=r8),INTENT(OUT) :: tpeng,timetree

      INTEGER,INTENT(IN) :: pot_type
      REAL(KIND=r8),INTENT(IN) :: kappa, eta, eps, T

      REAL(KIND=r8),PARAMETER :: kb = 0.001987215873_r8
      REAL(KIND=r8),PARAMETER :: coulomb = 332.0637790571_r8

! local variables

      TYPE(tnode),POINTER :: trootS,trootT
      INTEGER :: level
      INTEGER,DIMENSION(6) :: xyzind

! variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2 
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec
      REAL(KIND=r8)      :: totaltime


! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables. 

      CALL SETUP(xyzminmax,xyzdim,xyzind,order,theta)

! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(trootS,trootT)  

      CALL DATE_AND_TIME(datec,timec,zonec,time1)

! set global variables to track tree levels during construction

      level=0
      minlevel=50000

      WRITE(6,*) ' '
      WRITE(6,*) 'Creating tree...'

      maxlevel=0
      CALL CREATE_TREE_N0(trootT,maxparnodeT,xyzminmax,xyzdim,xyzind,level)

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,totaltime)
      timetree=totaltime

! print tree information to stdout 

         WRITE(6,*) ' '
         WRITE(6,*) 'Tree created. '
         WRITE(6,*) 'Tree parameters: '
         WRITE(6,*) ' '
         WRITE(6,*) '         numpar: ',trootT%numpar
         WRITE(6,*) '          x_mid: ',trootT%xyz_mid(1)
         WRITE(6,*) '          y_mid: ',trootT%xyz_mid(2)
         WRITE(6,*) '          z_mid: ',trootT%xyz_mid(3)
         WRITE(6,*) '         radius: ',trootT%radius   
         WRITE(6,*) '         torder: ',torder
         WRITE(6,*) '          theta: ',theta
         WRITE(6,*) '     maxparnode: ',maxparnodeT
         WRITE(6,*) ' '
 
      CALL DATE_AND_TIME(datec,timec,zonec,time1)

!Call driver routine for cluster-particle
      IF (pot_type == 0) THEN
          CALL CP_TREECODE_TCF(trootT, xS,yS,zS,qS, tEn, &
                               numparsS, numparsT, kappa, eta, eps)
      ELSE IF (pot_type == 1) THEN
          CALL CP_TREECODE_DCF(trootT, xS,yS,zS,qS, tEn, &
                               numparsS, numparsT, eta)
      ELSE IF (pot_type == 2) THEN
          CALL CP_TREECODE_COULOMB(trootT, xS,yS,zS,qS, tEn, &
                                   numparsS, numparsT)
      END IF

      tEn = tEn * sqrt(coulomb/kb/T)
      tpeng = SUM(tEn)

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,totaltime)
      timetree = timetree + totaltime

         WRITE(6,*) ' '
         WRITE(6,*) '   Finished calculation. '
         WRITE(6,*) '    Tree timing results: '
         WRITE(6,*) ' '
         WRITE(6,*) ' Tree creation time (s): ', timetree-totaltime
         WRITE(6,*) '  Treecode run time (s): ', totaltime
         WRITE(6,*) 'Treecode total time (s): ', timetree

! Call CLEANUP to deallocate global variables and tree structure.

         WRITE(6,*) ' '
         WRITE(6,*) 'Deallocating tree structure...'
         WRITE(6,*) ' '

      CALL CLEANUP(trootT)

      END SUBROUTINE TREECODE
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

