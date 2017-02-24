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
      MODULE treecode_procedures
      IMPLICIT NONE

! r8 is 8-byte (double precision) real

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      REAL(KIND=r8),PARAMETER :: pi = 4_r8*ATAN(1.D0)

      REAL(KIND=r8),PARAMETER :: kb = 0.001987215873_r8
      REAL(KIND=r8),PARAMETER :: coulomb = 332.0637790571_r8

! global variables for taylor expansions

      INTEGER :: torder,torderlim
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2,cf3
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:,:) :: a, b

! global variables to track tree levels 
 
      INTEGER :: minlevel,maxlevel
      
! global variables used when computing potential/force

      INTEGER :: orderoffset
      REAL(KIND=r8),DIMENSION(3) :: tarpos
      REAL(KIND=r8) :: thetasq,tarposq

! global variables for postition and charge storage

      INTEGER,ALLOCATABLE,DIMENSION(:)  :: orderarr

! node pointer and node type declarations

      TYPE tnode_pointer
           TYPE(tnode), POINTER :: p_to_tnode
      END TYPE tnode_pointer
      TYPE tnode
           INTEGER          :: numpar,ibeg,iend
           REAL(KIND=r8)    :: x_min,y_min,z_min
           REAL(KIND=r8)    :: x_max,y_max,z_max
           REAL(KIND=r8)    :: x_mid,y_mid,z_mid
           REAL(KIND=r8)    :: radius,sqradius,aspect
           INTEGER          :: level,num_children,exist_ms
           REAL(KIND=r8),DIMENSION(:,:,:),POINTER :: ms
           TYPE(tnode_pointer), DIMENSION(8) :: child
      END TYPE tnode

      CONTAINS
!!!!!!!!!!!!!!!
      SUBROUTINE SETUP(x,y,z,numpars,order,theta,xyzminmax) 
      IMPLICIT NONE
!
! SETUP allocates and initializes arrays needed for the Taylor expansion.
! Also, global variables are set and the Cartesian coordinates of
! the smallest box containing the particles is determined.
!
      INTEGER,INTENT(IN) :: numpars,order
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),INTENT(INOUT),DIMENSION(6) :: xyzminmax
      REAL(KIND=r8),INTENT(IN) :: theta

! local variables

      INTEGER :: err,i
      REAL(KIND=r8) :: t1

! global integers and reals:  TORDER, TORDERLIM and THETASQ

      torder=order
      orderoffset=1
      torderlim=torder+orderoffset
      thetasq=theta*theta

! allocate global Taylor expansion variables

      ALLOCATE(cf(0:torder), cf1(torderlim), cf2(torderlim), cf3(torderlim), &
               a(-2:torderlim, -2:torderlim, -2:torderlim), &
               b(-2:torderlim, -2:torderlim, -2:torderlim), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

! initialize arrays for Taylor sums and coeffs
      DO i = 0, torder
         cf(i) = REAL(i,KIND=r8) + 1.0_r8
      END DO

      DO i = 1, torderlim
         t1 = 1.0_r8 / REAL(i,KIND=r8)
         cf1(i) = t1;
         cf2(i) = 1.0_r8 - 0.5_r8 * t1
         cf3(i) = 1.0_r8 - t1
      END DO

! find bounds of Cartesian box enclosing the particles

      xyzminmax(1)=MINVAL(x(1:numpars))
      xyzminmax(2)=MAXVAL(x(1:numpars))
      xyzminmax(3)=MINVAL(y(1:numpars))
      xyzminmax(4)=MAXVAL(y(1:numpars))
      xyzminmax(5)=MINVAL(z(1:numpars))
      xyzminmax(6)=MAXVAL(z(1:numpars))

      ALLOCATE(orderarr(numpars),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating orderarr'
         STOP
      END IF  

      DO i=1,numpars
         orderarr(i)=i
      END DO  

      RETURN
      END SUBROUTINE SETUP
!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE CREATE_TREE_N0(p,ibeg,iend,x,y,z,shrink, &
                                      maxparnode,xyzmm,level,arrdim)
      IMPLICIT NONE
!
! CREATE_TREE_N0 recursively creates the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box. The division of a cluster terminates
! when the number of particles in a cluster are is less or equal to maxparnode

      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: ibeg,iend,shrink,level,maxparnode,arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm

! local variables

      REAL(KIND=r8) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax,t1,t2,t3
      INTEGER, DIMENSION(8,2) :: ind
      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER :: i,err,loclev,numposchild
      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm
     
! allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF

! set node fields: number of particles, exist_ms
! and xyz bounds 

      p%numpar=iend-ibeg+1
      p%exist_ms=0

      IF (shrink .EQ. 1) THEN   
         p%x_min=MINVAL(x(ibeg:iend))
         p%x_max=MAXVAL(x(ibeg:iend))
         p%y_min=MINVAL(y(ibeg:iend))
         p%y_max=MAXVAL(y(ibeg:iend))
         p%z_min=MINVAL(z(ibeg:iend))
         p%z_max=MAXVAL(z(ibeg:iend))
      ELSE
         p%x_min=xyzmm(1)
         p%x_max=xyzmm(2)
         p%y_min=xyzmm(3)
         p%y_max=xyzmm(4)
         p%z_min=xyzmm(5)
         p%z_max=xyzmm(6)        
      END IF

! compute aspect ratio

      xl=p%x_max-p%x_min
      yl=p%y_max-p%y_min
      zl=p%z_max-p%z_min

      lmax=MAX(xl,yl,zl)
      t1=lmax
      t2=MIN(xl,yl,zl)

      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=t1/t2
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS 

      p%x_mid=(p%x_max+p%x_min)/2.0_r8
      p%y_mid=(p%y_max+p%y_min)/2.0_r8
      p%z_mid=(p%z_max+p%z_min)/2.0_r8
      t1=p%x_max-p%x_mid
      t2=p%y_max-p%y_mid
      t3=p%z_max-p%z_mid
      p%sqradius=t1*t1+t2*t2+t3*t3
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%ibeg=ibeg
      p%iend=iend
      p%level=level

      IF (maxlevel .LT. level) THEN
         maxlevel=level
      END IF
      p%num_children=0
      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF (p%numpar .GT. maxparnode) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(1,1)=p%x_min
         xyzmms(2,1)=p%x_max
         xyzmms(3,1)=p%y_min
         xyzmms(4,1)=p%y_max
         xyzmms(5,1)=p%z_min
         xyzmms(6,1)=p%z_max
         ind(1,1)=ibeg
         ind(1,2)=iend
         x_mid=p%x_mid
         y_mid=p%y_mid
         z_mid=p%z_mid

         CALL PARTITION_8(x,y,z,xyzmms,xl,yl,zl,lmax,numposchild, &
                         x_mid,y_mid,z_mid,ind,arrdim)
!
! create children if indicated and store info in parent
!
         loclev=level+1

         DO i=1,numposchild
            IF (ind(i,1) .LE. ind(i,2)) THEN
               p%num_children=p%num_children+1
               lxyzmm=xyzmms(:,i)

               CALL CREATE_TREE_N0(p%child(p%num_children)%p_to_tnode, &
                               ind(i,1),ind(i,2),x,y,z,shrink, &
                               maxparnode,lxyzmm,loclev,arrdim)
            END IF
         END DO

      ELSE
         IF (level .LT. minlevel) THEN
            minlevel=level
         END IF
      END IF   

      END SUBROUTINE CREATE_TREE_N0      
!!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE CREATE_TREE_LV(p,ibeg,iend,x,y,z,shrink, &
                                          treelevel,xyzmm,level,arrdim)
      IMPLICIT NONE
!
! CREATE_TREE_LV recursively creates the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box. The subdivision terminates
! when the number of levels equals treelevel
!
      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: ibeg,iend,shrink,level,treelevel,arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm

! local variables

      REAL(KIND=r8) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax,t1,t2,t3
      INTEGER, DIMENSION(8,2) :: ind
      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER :: i,err,loclev,numposchild
      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm
     
! allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF
! set node fields: number of particles, exist_ms
! and xyz bounds 

      p%numpar=iend-ibeg+1
      p%exist_ms=0

      IF (shrink .EQ. 1) THEN   
         p%x_min=MINVAL(x(ibeg:iend))
         p%x_max=MAXVAL(x(ibeg:iend))
         p%y_min=MINVAL(y(ibeg:iend))
         p%y_max=MAXVAL(y(ibeg:iend))
         p%z_min=MINVAL(z(ibeg:iend))
         p%z_max=MAXVAL(z(ibeg:iend))
      ELSE
         p%x_min=xyzmm(1)
         p%x_max=xyzmm(2)
         p%y_min=xyzmm(3)
         p%y_max=xyzmm(4)
         p%z_min=xyzmm(5)
         p%z_max=xyzmm(6)        
      END IF

! compute aspect ratio

      xl=p%x_max-p%x_min
      yl=p%y_max-p%y_min
      zl=p%z_max-p%z_min

      lmax=MAX(xl,yl,zl)
      t1=lmax
      t2=MIN(xl,yl,zl)
      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=t1/t2
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS 

      p%x_mid=(p%x_max+p%x_min)/2.0_r8
      p%y_mid=(p%y_max+p%y_min)/2.0_r8
      p%z_mid=(p%z_max+p%z_min)/2.0_r8
      t1=p%x_max-p%x_mid
      t2=p%y_max-p%y_mid
      t3=p%z_max-p%z_mid
      p%sqradius=t1*t1+t2*t2+t3*t3
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%ibeg=ibeg
      p%iend=iend
      p%level=level
      p%num_children=0
      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF(level .LT. treelevel) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(1,1)=p%x_min
         xyzmms(2,1)=p%x_max
         xyzmms(3,1)=p%y_min
         xyzmms(4,1)=p%y_max
         xyzmms(5,1)=p%z_min
         xyzmms(6,1)=p%z_max
         ind(1,1)=ibeg
         ind(1,2)=iend
         x_mid=p%x_mid
         y_mid=p%y_mid
         z_mid=p%z_mid

         CALL PARTITION_8(x,y,z,xyzmms,xl,yl,zl,lmax,numposchild, &
                         x_mid,y_mid,z_mid,ind,arrdim)
!
! create children if indicated and store info in parent
!
         loclev=level+1
         DO i=1,numposchild
            IF (ind(i,1) .LE. ind(i,2)) THEN
               p%num_children=p%num_children+1
               lxyzmm=xyzmms(:,i)
               CALL CREATE_TREE_LV(p%child(p%num_children)%p_to_tnode, &
                               ind(i,1),ind(i,2),x,y,z,shrink, &
                               treelevel,lxyzmm,loclev,arrdim)
            END IF
            
         END DO
      ELSE
         IF (level .LT. minlevel) THEN
            minlevel=level
         END IF
      END IF   

      END SUBROUTINE CREATE_TREE_LV      


!!!!!!!!!!!!!!!
      SUBROUTINE PARTITION_8(x,y,z,xyzmms,xl,yl,zl,lmax,numposchild, &
                            x_mid,y_mid,z_mid,ind,arrdim)
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
      INTEGER, INTENT(IN) :: arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z
      INTEGER, DIMENSION(8,2),INTENT(INOUT) :: ind
      REAL(KIND=r8),DIMENSION(6,8),INTENT(INOUT) :: xyzmms
      REAL(KIND=r8), INTENT(IN) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax
      INTEGER,INTENT(INOUT) :: numposchild

! local variables

      INTEGER :: temp_ind,i
      REAL(KIND=r8) :: critlen

      numposchild=1
      critlen=lmax/sqrt(2.0_r8)

      IF (xl .GE. critlen) THEN

         CALL PARTITION(x,y,z,orderarr,ind(1,1),ind(1,2), &
                       x_mid,temp_ind,arrdim)

         ind(2,1)=temp_ind+1
         ind(2,2)=ind(1,2)
         ind(1,2)=temp_ind

         xyzmms(:,2)=xyzmms(:,1)
         xyzmms(2,1)=x_mid
         xyzmms(1,2)=x_mid
         numposchild=2*numposchild
      END IF 

      IF (yl .GE. critlen) THEN
         
         DO i=1,numposchild
            CALL PARTITION(y,x,z,orderarr,ind(i,1),ind(i,2), &
                          y_mid,temp_ind,arrdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind

            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(4,i)=y_mid
            xyzmms(3,numposchild+i)=y_mid
         END DO
         numposchild=2*numposchild

      END IF

      IF (zl .GE. critlen) THEN

         DO i=1,numposchild
            CALL PARTITION(z,x,y,orderarr,ind(i,1),ind(i,2), &
                          z_mid,temp_ind,arrdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind

            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(6,i)=z_mid
            xyzmms(5,numposchild+i)=z_mid
         END DO
         numposchild=2*numposchild

      END IF

      RETURN 
      END SUBROUTINE PARTITION_8
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CP_TREECODE_TCF(p,xS,yS,zS,qS,xT,yT,zT,tpeng,EnP,&
                                 numparsS,numparsT, kappa, eta, eps, T)
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
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(IN) :: xT,yT,zT
      REAL(KIND=r8),INTENT(IN) :: kappa, eta, eps, T

      REAL(KIND=r8),DIMENSION(numparsT),INTENT(INOUT) :: EnP
      REAL(KIND=r8),INTENT(INOUT) :: tpeng
 
! local variables

      INTEGER :: i,j
      REAL(KIND=r8) :: peng

      EnP=0.0_r8
      DO i=1,numparsS
         peng=0.0_r8
         tarpos(1)=xS(i)
         tarpos(2)=yS(i)
         tarpos(3)=zS(i)
         tarposq=qS(i)

         DO j=1,p%num_children
            CALL COMPUTE_CP1_TCF(p%child(j)%p_to_tnode,EnP, &
                                 xT,yT,zT,numparsT, kappa, eta)
         END DO
      END DO

      CALL COMPUTE_CP2(p,xT,yT,zT,EnP,numparsT)
 
      EnP = EnP * exp((kappa*eta)**2 / 4_r8) / (2_r8*eps) * sqrt(coulomb/kb/T) 
      tpeng = SUM(EnP)


      RETURN
      END SUBROUTINE CP_TREECODE_TCF
!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE COMPUTE_CP1_TCF(p,EnP,x,y,z,arrdim, kappa, eta)

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
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z
      REAL(KIND=r8),INTENT(IN) :: kappa, eta

! local variables

      REAL(KIND=r8) :: tx,ty,tz,distsq
      INTEGER :: i,err

! determine DISTSQ for MAC test

      tx=tarpos(1)-p%x_mid
      ty=tarpos(2)-p%y_mid
      tz=tarpos(3)-p%z_mid
      distsq=tx*tx+ty*ty+tz*tz

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8)) THEN

         a=0.0_r8
         b=0.0_r8
         !CALL COMP_TCOEFF_TCF(tx,ty,tz, kappa, eta)
         CALL COMP_TCOEFF_RECURSE(tx,ty,tz, kappa)
         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(0:torder,0:torder,0:torder),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF
             p%ms=0.0_r8
             p%exist_ms=1
         END IF
         CALL COMP_CMS(p)   
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT_TCF(EnP,p%ibeg,p%iend, &
                                 x,y,z,arrdim, kappa, eta)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1_TCF(p%child(i)%p_to_tnode,EnP, &
                                    x,y,z,arrdim, kappa, eta)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1_TCF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE COMPUTE_CP2(ap,x,y,z,EnP,arrdim)

        IMPLICIT NONE

! COMPUTE_CP2 is a recursive routine that evaluates the power series  
! approximation of the potential at the targets in a cluster via 
! a 3-D Horner's rule.  
!
        INTEGER,INTENT(IN) :: arrdim
        TYPE(tnode),POINTER :: ap
        REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP 
        REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z

! local variables

        REAL(KIND=r8) :: tx,ty,peng
        REAL(KIND=r8) :: xm,ym,zm,dx,dy,dz
        INTEGER :: i,nn,j,k1,k2,k3,porder,porder1

        porder=torder
        porder1=porder-1 
        IF(ap%exist_ms==1)THEN
          xm=ap%x_mid; ym=ap%y_mid; zm=ap%z_mid

         DO i=ap%ibeg,ap%iend
           nn=orderarr(i)
           dx=x(i)-xm; dy=y(i)-ym; dz=z(i)-zm
           peng=ap%ms(0,0,porder)           

           DO k3=porder1,0,-1

              ty=ap%ms(0,porder-k3,k3)

              DO k2=porder1-k3,0,-1

                 tx=ap%ms(porder-k3-k2,k2,k3)

                 DO k1=porder1-k3-k2,0,-1
                 
                   tx  = dx * tx + ap%ms(k1,k2,k3) 
                   
                 END DO

                 ty  = dy * ty + tx 
              END DO

               peng = dz * peng + ty  
             END DO

             EnP(nn)=EnP(nn)+peng
          END DO
        END IF

         DO j=1,ap%num_children
            CALL COMPUTE_CP2(ap%child(j)%p_to_tnode,x,y,z,EnP,arrdim)
         END DO

       RETURN
      END SUBROUTINE COMPUTE_CP2 

!!!!!!!!
      SUBROUTINE COMP_TCOEFF_TCF(dx,dy,dz, kappa, eta)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! directly.  The center of the expansion is the midpoint of the node P.  
! TARPOS and TORDERLIM are globally defined.
!
      REAL(KIND=r8),INTENT(IN) :: dx,dy,dz
      REAL(KIND=r8),INTENT(IN) :: kappa, eta

! local variables

      REAL(KIND=r8) :: rad, fp, g1p, g2p, h1p, h2p, h1arg, h2arg
      REAL(KIND=r8) :: kap_rad, kap_eta_2, rad_eta
      REAL(KIND=r8),DIMENSION(0:4, 0:2, 0:2, 0:2) :: fgh
      REAL(KIND=r8),DIMENSION(0:2) :: dd

! setup variables

      rad = SQRT(dx*dx + dy*dy + dz*dz)
      dd(0) = -dx
      dd(1) = -dy
      dd(2) = -dz
      kap_rad = kappa * rad

!      kap_eta_2 = kappa * eta / 2_r8
!      rad_eta = rad / eta

! 0th coeff or function val 

      !b1(0,0,0) = - 1_r8 / rad &
      !            * (exp(-kap_rad) * erfc(kap_eta_2 - rad_eta) &
      !            -  exp( kap_rad) * erfc(kap_eta_2 + rad_eta))

      a(0,0,0) = - 2_r8 / rad * exp(-kap_rad)

      IF (torder > 0) THEN
          fgh(0,0,0,0) = -1_r8 / rad
          fgh(1,0,0,0) = exp(-kap_rad)
          fgh(3,0,0,0) = exp( kap_rad)

          fp  = -1_r8 * fgh(0,0,0,0) / rad**2
          g1p = -kappa * fgh(1,0,0,0) / rad
          g2p =  kappa * fgh(3,0,0,0) / rad

          fgh(0,1,0,0) = fp * dd(0)
          fgh(0,0,1,0) = fp * dd(1)
          fgh(0,0,0,1) = fp * dd(2)

          fgh(1,1,0,0) = g1p * dd(0)
          fgh(1,0,1,0) = g1p * dd(1)
          fgh(1,0,0,1) = g1p * dd(2)

          fgh(3,1,0,0) = g2p * dd(0)
          fgh(3,0,1,0) = g2p * dd(1)
          fgh(3,0,0,1) = g2p * dd(2)

!----- Simplified code based on several approximation -----!
!----- fgh2000 = 2, fgh4000 = 0, h1p = h2p = 0 -----!

          a(1,0,0) = (fgh(0,1,0,0) * fgh(1,0,0,0)  &
                     + fgh(0,0,0,0) * fgh(1,1,0,0)) * 2_r8

          a(0,1,0) = (fgh(0,0,1,0) * fgh(1,0,0,0)  &
                     + fgh(0,0,0,0) * fgh(1,0,1,0)) * 2_r8

          a(0,0,1) = (fgh(0,0,0,1) * fgh(1,0,0,0)  &
                     + fgh(0,0,0,0) * fgh(1,0,0,1)) * 2_r8


!----- Original full precision code -----!
!
!          fgh(2,0,0,0) = erfc(kappa*eta/2_r8 - rad/eta) !Evaluates to approx 2
!          fgh(4,0,0,0) = erfc(kappa*eta/2_r8 + rad/eta) !Evaluates to approx 0
!
!         !h1p and h2p are approx 0
!          h1p = -2_r8 / eta / sqrt(pi) * exp(-(kappa*eta / 2_r8 - rad/eta)**2) / rad
!          h2p =  2_r8 / eta / sqrt(pi) * exp(-(kappa*eta / 2_r8 + rad/eta)**2) / rad
!
!          fgh(2,1,0,0) = h1p * dd(0)
!          fgh(2,0,1,0) = h1p * dd(1)
!          fgh(2,0,0,1) = h1p * dd(2)
!
!          fgh(4,1,0,0) = h2p * dd(0)
!          fgh(4,0,1,0) = h2p * dd(1)
!          fgh(4,0,0,1) = h2p * dd(2)
!
!          a(1,0,0) = fgh(0,1,0,0) * (fgh(1,0,0,0) * fgh(2,0,0,0)   &
!                                   - fgh(3,0,0,0) * fgh(4,0,0,0))  &
!
!                   + fgh(0,0,0,0) * (fgh(1,1,0,0) * fgh(2,0,0,0)   &
!                                   - fgh(3,1,0,0) * fgh(4,0,0,0)   &
!
!                                   + fgh(1,0,0,0) * fgh(2,1,0,0)   &
!                                   - fgh(3,0,0,0) * fgh(4,1,0,0))
!
!          a(0,1,0) = fgh(0,0,1,0) * (fgh(1,0,0,0) * fgh(2,0,0,0)   &
!                                   - fgh(3,0,0,0) * fgh(4,0,0,0))  &
!
!                   + fgh(0,0,0,0) * (fgh(1,0,1,0) * fgh(2,0,0,0)   &
!                                   - fgh(3,0,1,0) * fgh(4,0,0,0)   &
!
!                                   + fgh(1,0,0,0) * fgh(2,0,1,0)   &
!                                   - fgh(3,0,0,0) * fgh(4,0,1,0))
!
!          a(0,0,1) = fgh(0,0,0,1) * (fgh(1,0,0,0) * fgh(2,0,0,0)   &
!                                   - fgh(3,0,0,0) * fgh(4,0,0,0))  &
!                                                                    
!                   + fgh(0,0,0,0) * (fgh(1,0,0,1) * fgh(2,0,0,0)   &
!                                   - fgh(3,0,0,1) * fgh(4,0,0,0)   &
!                                                                    
!                                   + fgh(1,0,0,0) * fgh(2,0,0,1)   &
!                                   - fgh(3,0,0,0) * fgh(4,0,0,1))
!
!----- End of original full precision code -----!

      END IF

      RETURN
      END SUBROUTINE COMP_TCOEFF_TCF
!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!
      SUBROUTINE COMP_TCOEFF_RECURSE(dx,dy,dz,kappa)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.
!
      REAL(KIND=r8),INTENT(IN) :: dx,dy,dz
      REAL(KIND=r8),INTENT(IN)  :: kappa

! local varaibles

      REAL(KIND=r8) :: ddx,ddy,ddz,dist,fac
      REAL(KIND=r8) :: kappax,kappay,kappaz
      INTEGER :: i,j,k

! setup variables

      ddx=2.0_r8*dx
      ddy=2.0_r8*dy
      ddz=2.0_r8*dz

      kappax=kappa*dx
      kappay=kappa*dy
      kappaz=kappa*dz

      dist=dx*dx+dy*dy+dz*dz
      fac=1.0_r8/dist
      dist=SQRT(dist)

! 0th coeff or function val

      b(0,0,0)=-2_r8*EXP(-kappa*dist)
      a(0,0,0)=b(0,0,0)/dist

! 2 indices are 0

      b(1,0,0)=kappax*a(0,0,0)
      b(0,1,0)=kappay*a(0,0,0)
      b(0,0,1)=kappaz*a(0,0,0)

      a(1,0,0)=fac*dx*(a(0,0,0)+kappa*b(0,0,0))
      a(0,1,0)=fac*dy*(a(0,0,0)+kappa*b(0,0,0))
      a(0,0,1)=fac*dz*(a(0,0,0)+kappa*b(0,0,0))


      DO i=2,torderlim
         b(i,0,0)=cf1(i)*kappa*(dx*a(i-1,0,0)-a(i-2,0,0))
         b(0,i,0)=cf1(i)*kappa*(dy*a(0,i-1,0)-a(0,i-2,0))
         b(0,0,i)=cf1(i)*kappa*(dz*a(0,0,i-1)-a(0,0,i-2))

         a(i,0,0)=fac*(ddx*cf2(i)*a(i-1,0,0)-cf3(i)*a(i-2,0,0)+&
                  cf1(i)*kappa*(dx*b(i-1,0,0)-b(i-2,0,0)))
         a(0,i,0)=fac*(ddy*cf2(i)*a(0,i-1,0)-cf3(i)*a(0,i-2,0)+&
                  cf1(i)*kappa*(dy*b(0,i-1,0)-b(0,i-2,0)))
         a(0,0,i)=fac*(ddz*cf2(i)*a(0,0,i-1)-cf3(i)*a(0,0,i-2)+&
                  cf1(i)*kappa*(dz*b(0,0,i-1)-b(0,0,i-2)))
      END DO

! 1 index 0, 1 index 1, other >=1

      b(1,1,0)=kappax*a(0,1,0)
      b(1,0,1)=kappax*a(0,0,1)
      b(0,1,1)=kappay*a(0,0,1)

      a(1,1,0)=fac*(dx*a(0,1,0)+ddy*a(1,0,0)+kappax*b(0,1,0))
      a(1,0,1)=fac*(dx*a(0,0,1)+ddz*a(1,0,0)+kappax*b(0,0,1))
      a(0,1,1)=fac*(dy*a(0,0,1)+ddz*a(0,1,0)+kappay*b(0,0,1))

      DO i=2,torderlim-1
         b(1,0,i)=kappax*a(0,0,i)
         b(0,1,i)=kappay*a(0,0,i)
         b(0,i,1)=kappaz*a(0,i,0)
         b(1,i,0)=kappax*a(0,i,0)
         b(i,1,0)=kappay*a(i,0,0)
         b(i,0,1)=kappaz*a(i,0,0)

         a(1,0,i)=fac*(dx*a(0,0,i)+ddz*a(1,0,i-1)-a(1,0,i-2)+&
                  kappax*b(0,0,i))
         a(0,1,i)=fac*(dy*a(0,0,i)+ddz*a(0,1,i-1)-a(0,1,i-2)+&
                  kappay*b(0,0,i))
         a(0,i,1)=fac*(dz*a(0,i,0)+ddy*a(0,i-1,1)-a(0,i-2,1)+&
                  kappaz*b(0,i,0))
         a(1,i,0)=fac*(dx*a(0,i,0)+ddy*a(1,i-1,0)-a(1,i-2,0)+&
                  kappax*b(0,i,0))
         a(i,1,0)=fac*(dy*a(i,0,0)+ddx*a(i-1,1,0)-a(i-2,1,0)+&
                  kappay*b(i,0,0))
         a(i,0,1)=fac*(dz*a(i,0,0)+ddx*a(i-1,0,1)-a(i-2,0,1)+&
                  kappaz*b(i,0,0))
      END DO

! 1 index 0, others >= 2

      DO i=2,torderlim-2
         DO j=2,torderlim-i
            b(i,j,0)=cf1(i)*kappa*(dx*a(i-1,j,0)-a(i-2,j,0))
            b(i,0,j)=cf1(i)*kappa*(dx*a(i-1,0,j)-a(i-2,0,j))
            b(0,i,j)=cf1(i)*kappa*(dy*a(0,i-1,j)-a(0,i-2,j))

            a(i,j,0)=fac*(ddx*cf2(i)*a(i-1,j,0)+ddy*a(i,j-1,0) &
                     -cf3(i)*a(i-2,j,0)-a(i,j-2,0)+&
                     cf1(i)*kappa*(dx*b(i-1,j,0)-b(i-2,j,0)))
            a(i,0,j)=fac*(ddx*cf2(i)*a(i-1,0,j)+ddz*a(i,0,j-1)&
                     -cf3(i)*a(i-2,0,j)-a(i,0,j-2)+&
                     cf1(i)*kappa*(dx*b(i-1,0,j)-b(i-2,0,j)))
            a(0,i,j)=fac*(ddy*cf2(i)*a(0,i-1,j)+ddz*a(0,i,j-1)&
                     -cf3(i)*a(0,i-2,j)-a(0,i,j-2)+&
                     cf1(i)*kappa*(dy*b(0,i-1,j)-b(0,i-2,j)))
         END DO
      END DO

! 2 indices 1, other >= 1
! b(1,1,1) is correct, but a little tricky!
!      b(1,1,1)=5.0*dz*fac*b(1,1,0)

      b(1,1,1)=kappax*a(0,1,1)
      a(1,1,1)=fac*(dx*a(0,1,1)+ddy*a(1,0,1)+ddz*a(1,1,0)+&
               kappax*b(0,1,1))

      DO i=2,torderlim-2
         b(1,1,i)=kappax*a(0,1,i)
         b(1,i,1)=kappax*a(0,i,1)
         b(i,1,1)=kappay*a(i,0,1)

         a(1,1,i)=fac*(dx*a(0,1,i)+ddy*a(1,0,i)+ddz*a(1,1,i-1)&
                 -a(1,1,i-2)+kappax*b(0,1,i))
         a(1,i,1)=fac*(dx*a(0,i,1)+ddy*a(1,i-1,1)+ddz*a(1,i,0)&
                 -a(1,i-2,1)+kappax*b(0,i,1))
         a(i,1,1)=fac*(dy*a(i,0,1)+ddx*a(i-1,1,1)+ddz*a(i,1,0)&
                 -a(i-2,1,1)+kappay*b(i,0,1))
      END DO

! 1 index 1, others >=2

      DO i=2,torderlim-3
         DO j=2,torderlim-i
            b(1,i,j)=kappax*a(0,i,j)
            b(i,1,j)=kappay*a(i,0,j)
            b(i,j,1)=kappaz*a(i,j,0)

            a(1,i,j)=fac*(dx*a(0,i,j)+ddy*a(1,i-1,j)+ddz*a(1,i,j-1)&
                    -a(1,i-2,j)-a(1,i,j-2)+kappax*b(0,i,j))
            a(i,1,j)=fac*(dy*a(i,0,j)+ddx*a(i-1,1,j)+ddz*a(i,1,j-1)&
                    -a(i-2,1,j)-a(i,1,j-2)+kappay*b(i,0,j))
            a(i,j,1)=fac*(dz*a(i,j,0)+ddx*a(i-1,j,1)+ddy*a(i,j-1,1)&
                    -a(i-2,j,1)-a(i,j-2,1)+kappaz*b(i,j,0))

         END DO
      END DO

! all indices >=2

      DO k=2,torderlim-4
         DO j=2,torderlim-2-k
            DO i=2,torderlim-k-j
               b(i,j,k)=cf1(i)*kappa*(dx*a(i-1,j,k)-a(i-2,j,k))

               a(i,j,k)=fac*(ddx*cf2(i)*a(i-1,j,k)+ddy*a(i,j-1,k)&
                       +ddz*a(i,j,k-1)-cf3(i)*a(i-2,j,k)&
                       -a(i,j-2,k)-a(i,j,k-2)+&
                       cf1(i)*kappa*(dx*b(i-1,j,k)-b(i-2,j,k)))
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE COMP_TCOEFF_RECURSE
!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!
      SUBROUTINE COMP_CMS(p)
      IMPLICIT NONE
!
! COMP_CMS computes the moments for node P needed in the Taylor approximation
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: k1,k2,k3
        
      DO k3=0,torder
          
         DO k2=0,torder-k3
                
            DO k1=0,torder-k3-k2
               p%ms(k1,k2,k3)=p%ms(k1,k2,k3)+tarposq*a(k1,k2,k3)
                  
           END DO
             
         END DO
         
      END DO
         
      RETURN
      END SUBROUTINE COMP_CMS
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE COMP_DIRECT_TCF(EnP,ibeg,iend,x,y,z,arrdim, kappa, eta)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: ibeg,iend,arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z
      REAL(KIND=r8),INTENT(IN) :: kappa, eta

! local variables

      INTEGER :: i,nn
      REAL(KIND=r8) :: tx,ty,tz, rad
      REAL(KIND=r8) :: kap_eta_2, kap_rad, rad_eta

      DO i=ibeg,iend
         nn=orderarr(i)
         tx=x(i)-tarpos(1)
         ty=y(i)-tarpos(2)
         tz=z(i)-tarpos(3)

         rad = SQRT(tx*tx + ty*ty + tz*tz)
         kap_eta_2 = kappa * eta / 2_r8
         kap_rad = kappa * rad
         rad_eta = rad / eta

         EnP(nn) = EnP(nn) - tarposq / rad &
                           * (exp(-kap_rad) * erfc(kap_eta_2 - rad_eta) &
                           -  exp( kap_rad) * erfc(kap_eta_2 + rad_eta))

      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT_TCF
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

      DEALLOCATE(cf,cf1,cf2,cf3,a,b, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating Taylor variables! '
         STOP
      END IF      

      DEALLOCATE(orderarr,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating orderarr variables! '
         STOP
      END IF  

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
      SUBROUTINE TREECODE(xS,yS,zS,qS,xT,yT,zT,numparsS,numparsT, &
                      tEn,tpeng,order,theta, &
                      maxparnodeT,timetree,treelevelT,iflagT, &
                      pot_type, kappa, eta, eps, T)

      USE treecode_procedures
      IMPLICIT NONE

!====================================================================
!                                                                   
! xS,yS,zS,qS    :: x,y,z coordinates and charges of sources        
! xT,yT,zT       :: x,y,z coordinates of targets                    
! numparS        :: number of sources 
! numparT        :: number of targets
! tEn            :: array of dimension numparT for storing potential
!                   at each target
! tpeng          :: total potential
! shrinkT        :: flag for determining whether to shrink a source
!                   cluster to the smallest Cartesian box that
!                   encloses the particles in the cluster.
!                   No longer a parameter; all clusters are shrunk
! maxparnodeT    :: maximum number of particles in a leaf
!                   (only employed in the cluster-particle version)
! timetree       :: The total time for the treecode computation
! treelevelT     :: maximum number of levels (cluster-particle only)
! iflagT         :: if iflagT=0, the division of the target tree
!                   terminates when the number of particles in a leaf
!                   less than or equal to maxparnodeT. If iflagT is
!                   not equal to zero, then the divison terminates
!                   when the number of levels of the tree is equal
!                   to treelevelT (for cluster-particle only)
!=====================================================================

      INTEGER,INTENT(IN) :: numparsS,numparsT,order, &
                            maxparnodeT,treelevelT,iflagT
      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: xS,yS,zS,qS
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(INOUT) :: xT,yT,zT
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(OUT) :: tEn
      REAL(KIND=r8),INTENT(IN) :: theta
      REAL(KIND=r8),INTENT(OUT) :: tpeng,timetree

      INTEGER,INTENT(IN) :: pot_type
      REAL(KIND=r8),INTENT(IN) :: kappa, eta, eps, T

! local variables

      TYPE(tnode),POINTER :: trootS,trootT
      INTEGER :: level
      REAL(KIND=r8), DIMENSION(6) :: xyzminmax

! variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2 
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec
      REAL(KIND=r8)      :: totaltime


! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables. 

      CALL SETUP(xT,yT,zT,numparsT,order,theta,xyzminmax)

! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(trootS,trootT)  

      CALL DATE_AND_TIME(datec,timec,zonec,time1)

! set global variables to track tree levels during construction

      level=0
      minlevel=50000

         WRITE(6,*) ' '
         WRITE(6,*) 'Creating tree...'

      IF(iflagT .EQ. 0) THEN
         maxlevel=0
         CALL CREATE_TREE_N0(trootT,1,numparsT,xT,yT,zT,1, &
                      maxparnodeT,xyzminmax,level,numparsT)
      ELSE
         maxlevel=treelevelT
         CALL CREATE_TREE_LV(trootT,1,numparsT,xT,yT,zT,1, &
                      treelevelT,xyzminmax,level,numparsT)
      END IF
      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,totaltime)
      timetree=totaltime

! print tree information to stdout 

         WRITE(6,*) ' '
         WRITE(6,*) 'Tree created. '
         WRITE(6,*) 'Tree parameters: '
         WRITE(6,*) ' '
         WRITE(6,*) '         numpar: ',trootT%numpar
         WRITE(6,*) '          x_mid: ',trootT%x_mid
         WRITE(6,*) '          y_mid: ',trootT%y_mid
         WRITE(6,*) '          z_mid: ',trootT%z_mid
         WRITE(6,*) '         radius: ',trootT%radius   
         WRITE(6,*) '         torder: ',torder
         WRITE(6,*) '          theta: ',theta
         WRITE(6,*) '         shrink: ',1
         WRITE(6,*) '     maxparnode: ',maxparnodeT
         WRITE(6,*) '          iflag: ',iflagT
         WRITE(6,*) '  tree maxlevel: ',treelevelT
         WRITE(6,*) ' '
 
      CALL DATE_AND_TIME(datec,timec,zonec,time1)

!Call driver routine for cluster-particle
      IF (pot_type == 0) THEN
          CALL CP_TREECODE_TCF(trootT,xS,yS,zS,qS,xT,yT,zT,tpeng,tEn,&
                               numparsS,numparsT, kappa, eta, eps, T) 
      ELSE IF (pot_type == 1) THEN
          !CALL CP_TREECODE_DCF(trootT,xS,yS,zS,qS,xT,yT,zT,tpeng,tEn,&
          !                     numparsS,numparsT, kappa, eta, eps, T) 
      END IF

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
!!!!!!!!!!!!!
      SUBROUTINE PARTITION(a,b,c,indarr,ibeg,iend,val,midind,arrdim)
      IMPLICIT NONE
!
! PARTITION determines the index MIDIND, after partitioning
! in place the  arrays A,B,C and Q,  such that 
! A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL. 
! If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
! is returned as IBEG-1. 
! 
      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER, INTENT(IN) :: arrdim,ibeg,iend
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: a,b,c
      INTEGER,DIMENSION(arrdim),INTENT(INOUT) :: indarr   
      INTEGER, INTENT(INOUT) :: midind   
      REAL(KIND=r8) val

! local variables

      REAL(KIND=r8) ta,tb,tc
      INTEGER lower,upper,tind

      IF (ibeg .LT. iend) THEN

! temporarily store IBEG entries and set A(IBEG)=VAL for 
! the partitoning algorithm.  

         ta=a(ibeg)
         tb=b(ibeg)
         tc=c(ibeg)
         tind=indarr(ibeg)
         a(ibeg)=val 
         upper=ibeg
         lower=iend

         DO WHILE (upper .NE. lower)
            DO WHILE ((upper .LT. lower) .AND. (val .LT. a(lower)))
                  lower=lower-1
            END DO
            IF (upper .NE. lower) THEN
               a(upper)=a(lower)
               b(upper)=b(lower)
               c(upper)=c(lower)
               indarr(upper)=indarr(lower)
            END IF
            DO WHILE ((upper .LT. lower) .AND. (val .GE. a(upper)))
                  upper=upper+1
            END DO
            IF (upper .NE. lower) THEN
               a(lower)=a(upper)
               b(lower)=b(upper)
               c(lower)=c(upper)
               indarr(lower)=indarr(upper)
            END IF
         END DO
         midind=upper

! replace TA in position UPPER and change MIDIND if TA > VAL 

         IF (ta .GT. val) THEN
            midind=upper-1
         END IF
         a(upper)=ta
         b(upper)=tb
         c(upper)=tc
         indarr(upper)=tind

      ELSEIF (ibeg .EQ. iend) THEN
         IF (a(ibeg) .LE. val) THEN
            midind=ibeg
         ELSE
            midind=ibeg-1
         END IF
      ELSE
         midind=ibeg-1
      END IF

      RETURN
      END SUBROUTINE PARTITION
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

