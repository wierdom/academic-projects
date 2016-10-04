!-----------!
 MODULE bmsr

 USE hyzer, ONLY: nn;

 IMPLICIT NONE

 INTEGER ele(-1:1,-1:1,4)
 REAL(8),SAVE :: avu,avk,umag,sxu,ssa,sxa,rhx,rhy
 
 END MODULE bmsr
!---------------!

!-----------!
 MODULE bvxp

 USE hyzer, ONLY: nvx,maxp;

 REAL(8),SAVE :: vxp(4,4,nvx,maxp)

 END MODULE bvxp
!---------------!

!-----------!
 MODULE bwgp

 USE hyzer, ONLY: maxp;

 REAL(8),SAVE :: wgt(-1:1,-1:1,4,maxp)

 END MODULE bwgp
!---------------!

!-----------------------------------------------------!
 SUBROUTINE simulation (istep,nruns,mstep,tfr,np,pnum)

 USE hyzer; USE mpi;

 IMPLICIT NONE
      
 INTEGER i,j,l1,ie,np,tfr,pnum,istep,nruns,mstep
 INTEGER ic,oc,vx
 INTEGER stat(MPI_STATUS_SIZE),iarr(2)
 REAL*8  rarr(2)

 CALL lattice
 CALL pvect0(np)
 CALL pvect1
 CALL pvect2(pnum)
! CALL initvrtx(np,pnum)
 CALL initvrtx_dirloop(np,pnum)
 CALL pvect3(pnum)

 IF (istep.NE.0) THEN
   lopers=0.d0
   nloops=0.d0
   DO i=1,istep
       CALL mcstep
       CALL checkl
       IF (MOD(i,tfr).EQ.0) THEN
          CALL sendres(0,tfr)
          CALL MPI_RECV(iarr,2,MPI_INTEGER,0,7,MPI_COMM_WORLD,stat,ie)
          CALL MPI_RECV(rarr,2,MPI_DOUBLE_PRECISION,0,8,MPI_COMM_WORLD,stat,ie)
          pnum=iarr(1)
          nl=iarr(2)
          nloops=rarr(1)
          lopers=rarr(2)               
          CALL pvect2(pnum)               
          CALL pvect3(pnum)               
       ENDIF
       IF (MOD(i,istep/20).EQ.0) CALL adjnl
    ENDDO
    CALL sendconf
    CALL MPI_RECV(l,1,MPI_INTEGER,0,9,MPI_COMM_WORLD,stat,ie)
    CALL pvect1
 ENDIF

 DO j=1,nruns
    CALL zerodat1
    DO i=1,mstep
       CALL mcstep
       IF (MOD(i,tfr).EQ.0) THEN
          CALL sendres(1,tfr)
          CALL MPI_RECV(iarr,2,MPI_INTEGER,0,7,MPI_COMM_WORLD,stat,ie)
          pnum=iarr(1)
          nl=iarr(2)
          CALL pvect2(pnum)               
          CALL pvect3(pnum)               
       ENDIF
    ENDDO
    CALL sendconf
 ENDDO

 END SUBROUTINE simulation
!-------------------------!

!-------------------!
 SUBROUTINE SENDCONF

 USE hyzer; USE mpi;

 IMPLICIT NONE

 INTEGER arr(3),ie

 arr(1)=l
 arr(2)=nh
 arr(3)=nl
 CALL MPI_SSEND(arr,3,MPI_INTEGER,0,5,MPI_COMM_WORLD,ie)
 CALL MPI_SSEND(spn,nn,MPI_INTEGER,0,6,MPI_COMM_WORLD,ie)      
 CALL MPI_SSEND(str,2*ll,MPI_INTEGER,0,7,MPI_COMM_WORLD,ie)      

 END SUBROUTINE SENDCONF
!-----------------------!

!-----------------------------!
 SUBROUTINE sendres (ires,tfr)

 USE bmsr; USE hyzer; USE mpi;

 IMPLICIT NONE

 INTEGER ie,tfr,ires
 REAL*8  res(8),rarr(2)

 IF (ires.NE.0) THEN
    res(1)=avu/DFLOAT(tfr)
    res(2)=avk/DFLOAT(tfr)
    res(3)=umag/DFLOAT(tfr)
    res(4)=sxu/DFLOAT(tfr)
    res(5)=ssa/DFLOAT(tfr)
    res(6)=sxa/DFLOAT(tfr)
    res(7)=rhx/DFLOAT(tfr)
    res(8)=rhy/DFLOAT(tfr)
    CALL MPI_SSEND(res,8,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,ie)      
 ELSE
    CALL MPI_SSEND(nl,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ie)  
    rarr(1)=nloops
    rarr(2)=lopers
    CALL MPI_SSEND(rarr,2,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ie)  
 ENDIF
 CALL MPI_SSEND(ele,36,MPI_INTEGER,0,4,MPI_COMM_WORLD,ie)      
 CALL zerodat1

 END SUBROUTINE sendres
!----------------------!

!-------------------!
 SUBROUTINE zerodat1

 USE bmsr;

 IMPLICIT NONE

 INTEGER i

 avu=0.d0
 avk=0.d0
 umag=0.d0
 sxu=0.d0
 ssa=0.d0
 sxa=0.d0
 rhx=0.d0
 rhy=0.d0

 END SUBROUTINE zerodat1
!-----------------------!

!-----------------!
 SUBROUTINE checkl

 USE hyzer; USE random;

 IMPLICIT NONE

 INTEGER i,j,p,dl,l1,tmp(2,ll)
 LOGICAL lstr(ll)

 dl=l/10+2
 IF (nh.LT.l-dl/2) RETURN
 l1=l+dl
 DO i=1,l1
    lstr(i)=.TRUE.
 ENDDO
 DO i=1,dl
    DO
       p=MIN(INT(ranmz()*DFLOAT(l1))+1,l1)
       IF (lstr(p)) THEN
          lstr(p)=.FALSE.
          EXIT
       ENDIF
    ENDDO
 ENDDO
 j=0
 DO i=1,l1
    IF (lstr(i)) THEN
       j=j+1
       tmp(1,i)=str(1,j)
       tmp(2,i)=str(2,j)
    ELSE
       tmp(1,i)=0
       tmp(2,i)=0
    ENDIF            
 ENDDO
 l=l1
 DO i=1,l
    str(1,i)=tmp(1,i)
    str(2,i)=tmp(2,i)
 ENDDO
 CALL pvect1

 END SUBROUTINE checkl
!---------------------!

!------------------!
 SUBROUTINE lattice

 USE hyzer;

 IMPLICIT NONE

 INTEGER i,j,is,i1,i2,ix1,iy1,ix2,iy2,b
 INTEGER xy(2,nn),xyi(0:nx-1,0:ny-1)

 is=0
 DO j=0,ny-1
    DO i=0,nx-1
       is=is+1
       xy(1,is)=i
       xy(2,is)=j
       xyi(i,j)=is
    ENDDO
 ENDDO
 DO i=1,nn
    ix1=xy(1,i)
    iy1=xy(2,i)
    ix2=MOD(ix1+1,nx)
    iy2=iy1
    bst(1,i)=i
    bst(2,i)=xyi(ix2,iy2)
    btyp(i)=1
    ix2=ix1
    iy2=MOD(iy1+1,ny)
    bst(1,i+nn)=i
    bst(2,i+nn)=xyi(ix2,iy2)
    btyp(i+nn)=1
 ENDDO
 b=2*nn
 DO i=1,nn
    ix1=xy(1,i)
    iy1=xy(2,i)
    IF (MOD(iy1,2).EQ.0) THEN
       IF (MOD((ix1+iy1),2).EQ.0) THEN
          ix2=MOD(ix1+1,nx)
          iy2=MOD(iy1+1,ny)
          b=b+1
          bst(1,b)=i
          bst(2,b)=xyi(ix2,iy2)
          btyp(b)=2
       ENDIF
    ELSE
       IF (MOD((ix1+iy1),2).EQ.1) THEN
          ix2=MOD(ix1-1+nx,nx)
          iy2=MOD(iy1+1,ny)
          b=b+1
          bst(1,b)=i
          bst(2,b)=xyi(ix2,iy2)
          btyp(b)=2
       ENDIF
    ENDIF
 ENDDO
 b=2*nn+nn/2
 DO i=1,nn
    ix1=xy(1,i)
    iy1=xy(2,i)
    IF (MOD(iy1,2).EQ.0) THEN
       IF (MOD((ix1+iy1),2).EQ.0) THEN
          ix2=MOD(ix1-1+nx,nx)
          iy2=MOD(iy1+1,ny)
          b=b+1
          bst(1,b)=i
          bst(2,b)=xyi(ix2,iy2)
          btyp(b)=3
       ELSE
          ix2=MOD(ix1+1,nx)
          iy2=MOD(iy1+1,ny)
          b=b+1
          bst(1,b)=i
          bst(2,b)=xyi(ix2,iy2)
          btyp(b)=3
       ENDIF
    ELSE
       IF (MOD((ix1+iy1),2).EQ.1) THEN
          ix2=MOD(ix1+1,nx)
          iy2=MOD(iy1+1,ny)
          b=b+1
          bst(1,b)=i
          bst(2,b)=xyi(ix2,iy2)
          btyp(b)=3
       ELSE
          ix2=MOD(ix1-1+nx,nx)
          iy2=MOD(iy1+1,ny)
          b=b+1
          bst(1,b)=i
          bst(2,b)=xyi(ix2,iy2)
          btyp(b)=3
       ENDIF
    ENDIF
 ENDDO
 b=3*nn+nn/2
 DO i=1,nn
    ix1=xy(1,i)
    iy1=xy(2,i)
    ix2=MOD(ix1+2,nx)
    iy2=iy1
    b=b+1
    bst(1,b)=i
    bst(2,b)=xyi(ix2,iy2)
    btyp(b)=4
    ix2=ix1
    iy2=MOD(iy1+2,ny)
    b=b+1
    bst(1,b)=i
    bst(2,b)=xyi(ix2,iy2)
    btyp(b)=4
 ENDDO

 DO i=1,nn
    ph(i)=(-1)**(xy(1,i)+xy(2,i))
 ENDDO

 END SUBROUTINE lattice
!----------------------!

!----------------------!
 SUBROUTINE pvect0 (np)

 USE bwgp; USE hyzer;

 IMPLICIT NONE

 INTEGER i,j,s1,s2,np
 REAL*8 z

 z=11.d0
 DO i=1,np
 DO j=1,4
 DO s1=-1,1,2
 DO s2=-1,1,2
    wgt(s1,s2,j,i)=0.d0
 ENDDO
 ENDDO
 ENDDO
 ENDDO

 DO i=1,np
 DO j=1,4
    amax(j,i)=0.d0
 ENDDO
 ENDDO

 DO i=1,np
 DO j=1,4
 DO s1=-1,1,2
 DO s2=-1,1,2
    wgt(s1,s2,j,i)=jz(j)*0.25d0*s1*s2
    wgt(s1,s2,j,i)=wgt(s1,s2,j,i)-(1.d0/z)*hz(i)*0.5d0*(s1+s2)
    IF (wgt(s1,s2,j,i).GT.amax(j,i)) amax(j,i)=wgt(s1,s2,j,i)
 ENDDO
 ENDDO
 ENDDO
 ENDDO

 DO i=1,np
 DO j=1,4
    amax(j,i)=amax(j,i)+1.d0
    DO s1=-1,1,2
    DO s2=-1,1,2
       wgt(s1,s2,j,i)=amax(j,i)-wgt(s1,s2,j,i)
    ENDDO
    ENDDO
 ENDDO
 ENDDO

 END SUBROUTINE pvect0
!---------------------!

!-----------------!
 SUBROUTINE pvect1

 USE hyzer;

 IMPLICIT NONE

 INTEGER i

 DO i=0,l-1
    addp(i)=beta*DFLOAT(nb)/DFLOAT(l-i)
 ENDDO
 DO i=1,l
    delp(i)=DFLOAT(l-i+1)/(beta*DFLOAT(nb))
 ENDDO

 RETURN
 END SUBROUTINE pvect1
!---------------------!

!------------------------!
 SUBROUTINE pvect2 (pnum)

 USE bwgp; USE hyzer;

 IMPLICIT NONE

 INTEGER i,s1,s2,pnum

 DO s1=-1,1,2
    DO s2=-1,1,2
       DO i=1,4
          awgt(s1,s2,i)=wgt(s1,s2,i,pnum)
          IF (awgt(s1,s2,i).GT.1.d-6) THEN
             dwgt(s1,s2,i)=1.d0/awgt(s1,s2,i)
          ELSE
             dwgt(s1,s2,i)=1.d6
          ENDIF
       ENDDO
    ENDDO
 ENDDO

 END SUBROUTINE pvect2
!---------------------!

!------------------------!
 SUBROUTINE pvect3 (pnum)

 USE bvxp; USE hyzer;

 IMPLICIT NONE

 INTEGER i,vx,ic,is,oc,pnum

 DO vx=1,nvx
 DO ic=1,4
 DO oc=1,4
    vxprb(oc,ic,vx)=vxp(oc,ic,vx,pnum)
 END DO
 END DO
 END DO

 END SUBROUTINE pvect3
!---------------------!

!-----------------------------!
 SUBROUTINE initvrtx (np,pnum)

 USE bvxp; USE bwgp; USE hyzer;

 IMPLICIT NONE

 INTEGER np,v1,v2,v3,v4,v5,v6,v0,ii,pnum
 INTEGER i,j,op,s1,s2,t1,t2,vx,ic,oc,is,os,vxn
 INTEGER st(4)
 INTEGER vxtyp(nvx)
 REAL*8 vxwgt(nvx,maxp)

 legvx=0 ! legvx(-1:1,-1:1,-1:1,-1:1,4)=0 !

 DO j=1,4
    v0=6*(j-1)
    v1=v0+1
    v2=v0+2
    v3=v0+3
    v4=v0+4
    v5=v0+5
    v6=v0+6

    legvx(1,-1,-1,1,j)=v1
    legvx(-1,1,1,-1,j)=v2
    legvx(1,-1,1,-1,j)=v3
    legvx(-1,1,-1,1,j)=v4
    legvx(1,1,1,1,j)=v5
    legvx(-1,-1,-1,-1,j)=v6

    DO ii=1,6
       vxtyp(v0+ii)=j
    ENDDO

    vxoper(v1)=2
    vxoper(v2)=2
    vxoper(v3)=1
    vxoper(v4)=1
    vxoper(v5)=1
    vxoper(v6)=1

    DO ii=1,np
       vxwgt(v1,ii)=0.5*jx(j)
       vxwgt(v2,ii)=0.5*jx(j)
       vxwgt(v3,ii)=wgt(1,-1,j,ii)
       vxwgt(v4,ii)=wgt(-1,1,j,ii)
       vxwgt(v5,ii)=wgt(1,1,j,ii)
       vxwgt(v6,ii)=wgt(-1,-1,j,ii)
    ENDDO

    DO t2=-1,1,2
       DO t1=-1,1,2
          DO s2=-1,1,2
             DO s1=-1,1,2
                vx=legvx(s1,s2,t1,t2,j)
                IF (vx.NE.0) then
                   vxleg(1,vx)=s1
                   vxleg(2,vx)=s2
                   vxleg(3,vx)=t1
                   vxleg(4,vx)=t2
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

 ENDDO

 vxnew=0 ! vxnew(4,4,nvx)=0 !
 vxp=0.d0 ! vxp(4,4,nvx,maxp)=0.d0 !

 DO vx=1,nvx
    j=vxtyp(vx)
    DO ic=1,4
       DO oc=1,4
          st(1)=vxleg(1,vx)
          st(2)=vxleg(2,vx)
          st(3)=vxleg(3,vx)
          st(4)=vxleg(4,vx)
          st(ic)=-st(ic)
          st(oc)=-st(oc)
          vxn=legvx(st(1),st(2),st(3),st(4),j)
          IF (vxn.NE.0) THEN
             vxnew(oc,ic,vx)=vxn
             DO i=1,np
                vxp(oc,ic,vx,i)=vxwgt(vxn,i)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
 ENDDO

 DO i=1,np
    DO vx=1,nvx
       DO ic=1,4
          DO oc=2,4
             vxp(oc,ic,vx,i)=vxp(oc,ic,vx,i)+vxp(oc-1,ic,vx,i)
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 DO i=1,np
    DO vx=1,nvx
       DO ic=1,4
          DO oc=1,4
             IF (vxp(4,ic,vx,i).LT.1.d-6) THEN
                vxp(oc,ic,vx,i)=0.d0
             ELSE
                vxp(oc,ic,vx,i)=vxp(oc,ic,vx,i)/vxp(4,ic,vx,i)
                IF (vxp(oc,ic,vx,i).LT.1.d-6) vxp(oc,ic,vx,i)=-1.d0
             ENDIF
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 END SUBROUTINE initvrtx
!-----------------------!

!-------------------------------------!
 SUBROUTINE initvrtx_dirloop (np,pnum)

 USE bvxp; USE bwgp; USE hyzer;

 IMPLICIT NONE

 INTEGER np,v1,v2,v3,v4,v5,v6,v0,ii,pnum
 INTEGER i,j,k,op,s1,s2,t1,t2,vx,ic,oc,is,os,vxn
 INTEGER st(4)
 INTEGER vxtyp(nvx)
 REAL*8 vxwgt(nvx,maxp),legwgt(4,4,nvx,maxp)

 INTEGER jn(0:3),jo(0:3),jv(0:3)
 REAL*8 fac,jw(0:3),mwgt(3,3)

 INTEGER jdog(3),kdog(3),icn,nj,nk,vxk
 DATA jdog/1,2,1/
 DATA kdog/2,3,2/

 legvx=0 ! legvx(-1:1,-1:1,-1:1,-1:1,4)=0 !

 !Define legvx, vxoper, vxwgt, and vxleg
 DO j=1,4
    v0=6*(j-1)
    v1=v0+1
    v2=v0+2
    v3=v0+3
    v4=v0+4
    v5=v0+5
    v6=v0+6

    legvx(1,-1,-1,1,j)=v1
    legvx(-1,1,1,-1,j)=v2
    legvx(1,-1,1,-1,j)=v3
    legvx(-1,1,-1,1,j)=v4
    legvx(1,1,1,1,j)=v5
    legvx(-1,-1,-1,-1,j)=v6

    DO ii=1,6
       vxtyp(v0+ii)=j
    ENDDO

    vxoper(v1)=2
    vxoper(v2)=2
    vxoper(v3)=1
    vxoper(v4)=1
    vxoper(v5)=1
    vxoper(v6)=1

    DO ii=1,np
       vxwgt(v1,ii)=0.5*jx(j)
       vxwgt(v2,ii)=0.5*jx(j)
       vxwgt(v3,ii)=wgt(1,-1,j,ii)
       vxwgt(v4,ii)=wgt(-1,1,j,ii)
       vxwgt(v5,ii)=wgt(1,1,j,ii)
       vxwgt(v6,ii)=wgt(-1,-1,j,ii)
    ENDDO

    DO t2=-1,1,2
       DO t1=-1,1,2
          DO s2=-1,1,2
             DO s1=-1,1,2
                vx=legvx(s1,s2,t1,t2,j)
                IF (vx.NE.0) then
                   vxleg(1,vx)=s1
                   vxleg(2,vx)=s2
                   vxleg(3,vx)=t1
                   vxleg(4,vx)=t2
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

 ENDDO

 legwgt=0.d0 ! legwgt(4,4,nvx,maxp)=0.d0 !
 vxnew=0 ! vxnew(4,4,nvx)=0 !
 vxp=0.d0 ! vxp(4,4,nvx,maxp)=0.d0 !

 !Create vxnew links
 DO vx=1,nvx
    j=vxtyp(vx)
    DO ic=1,4
       DO oc=1,4
          st(1)=vxleg(1,vx)
          st(2)=vxleg(2,vx)
          st(3)=vxleg(3,vx)
          st(4)=vxleg(4,vx)
          st(ic)=-st(ic)
          st(oc)=-st(oc)
          vxn=legvx(st(1),st(2),st(3),st(4),j)
          IF (vxn.NE.0) THEN
             vxnew(oc,ic,vx)=vxn
          ENDIF
       ENDDO
    ENDDO
 ENDDO

 !Choose a vertex and in channel
 DO vx=1,nvx
    DO ic=1,4
       !Now we have list of vertices in Directed Loop equations
       DO i=1,np
          !Initialize vertex order
          j=0
          DO oc=1,4
             vxn=vxnew(oc,ic,vx)
             IF (vxn.NE.0) THEN
                j=j+1
                jn(j)=j
                jo(j)=oc
                jv(j)=vxn
                jw(j)=vxwgt(vxn,i)
             ENDIF
          ENDDO
          IF (j.NE.3) STOP 'DL FAIL 1'
          !Three moves to relative order
          DO oc=1,3
             j=jdog(oc)
             k=kdog(oc)
             IF (jw(j).LT.jw(k)) THEN
                jn(0)=jn(j)
                jn(j)=jn(k)
                jn(k)=jn(0)
                jw(0)=jw(j)
                jw(j)=jw(k)
                jw(k)=jw(0)
             ENDIF
          ENDDO
          !Implement Eqs. 22 and 23 of Sylju\r{a}sen
          mwgt(1,1)=MAX(0.d0,jw(1)-jw(2)-jw(3))
          mwgt(1,2)=MIN(jw(2),0.5d0*(jw(1)+jw(2)-jw(3)))
          mwgt(1,3)=MIN(jw(3),0.5d0*(jw(1)-jw(2)+jw(3)))
          mwgt(2,1)=MIN(jw(2),0.5d0*(jw(1)+jw(2)-jw(3)))
          mwgt(2,2)=0.d0
          mwgt(2,3)=MAX(0.d0,0.5d0*(-jw(1)+jw(2)+jw(3)))
          mwgt(3,1)=MIN(jw(3),0.5d0*(jw(1)-jw(2)+jw(3)))
          mwgt(3,2)=MAX(0.d0,0.5d0*(-jw(1)+jw(2)+jw(3)))
          mwgt(3,3)=0.d0
          !Finally, assign weights to Directed Loop vertex elements
          DO j=1,3
             nj=jn(j)
             icn=jo(nj)
             vxn=jv(nj)
             DO k=1,3
                nk=jn(k)
                vxk=jv(nk)
                ! Cycle through to find correct out channel
                ! (the one that changes vxj to vxk)
                DO oc=1,4
                   IF (vxnew(oc,icn,vxn).EQ.vxk) EXIT
                   IF (oc.EQ.4) STOP 'DL FAIL 2'
                ENDDO
                legwgt(oc,icn,vxn,i)=mwgt(j,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 !Convert matrix weights to probabilities
 DO i=1,np
    DO vx=1,nvx
       fac=1.d0/vxwgt(vx,i)
       DO ic=1,4
          DO oc=1,4
             vxp(oc,ic,vx,i)=fac*legwgt(oc,ic,vx,i)
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 !Create cumulative probabilities
 DO i=1,np
    DO vx=1,nvx
       DO ic=1,4
          DO oc=2,4
             vxp(oc,ic,vx,i)=vxp(oc,ic,vx,i)+vxp(oc-1,ic,vx,i)
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 !Normalize cumulative probabilities
 DO i=1,np
    DO vx=1,nvx
       DO ic=1,4
          DO oc=1,4
             IF (vxp(4,ic,vx,i).LT.1.d-6) THEN
                vxp(oc,ic,vx,i)=0.d0
             ELSE
                vxp(oc,ic,vx,i)=vxp(oc,ic,vx,i)/vxp(4,ic,vx,i)
                IF (vxp(oc,ic,vx,i).LT.1.d-6) vxp(oc,ic,vx,i)=-1.d0
             ENDIF
          ENDDO
       ENDDO
    ENDDO
 ENDDO

 END SUBROUTINE initvrtx_dirloop
!-------------------------------!
