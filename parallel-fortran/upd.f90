!------------!
 MODULE blink

 USE hyzer, ONLY: ll,nn;

 IMPLICIT NONE

 INTEGER,SAVE :: no
 INTEGER,SAVE :: frst(nn),fspn(nn),lpos(ll),lvtx(ll),plnk(4,ll),slnk(4,ll)

 END MODULE blink
!----------------!

!-----------------!
 SUBROUTINE mcstep

 IMPLICIT NONE

 LOGICAL :: passed

 DO
    CALL dupdate
    CALL linkoper
    CALL updloop(passed)
    IF (passed) EXIT
 ENDDO

 END SUBROUTINE mcstep
!---------------------!

!------------------!
 SUBROUTINE dupdate

 USE bmsr; USE hyzer; USE random;

 IMPLICIT NONE

 INTEGER i,j,o,b,s1,s2,jjx,jjy,cu,ca
 INTEGER ss1,ss2,nu,nnu(4),nk,nh1,nmc,last,typ
 REAL*8  ssa1,sxa1,p,dl

 cu=0
 ca=0
 jjx=0
 jjy=0
 nmc=0
 DO s2=-1,1
 DO s1=-1,1
    ele(s1,s2,1)=0
    ele(s1,s2,2)=0
    ele(s1,s2,3)=0
    ele(s1,s2,4)=0
 ENDDO
 ENDDO
 DO i=1,nn
    cu=cu+spn(i)
    ca=ca+phase(i)*spn(i)
 ENDDO      
 ssa1=DFLOAT(ca*ca)
 sxa1=0.d0
 nu=0
 nk=0
 do i=1,4
    nnu(i)=0
 enddo
 nh1=0
 last=0
 DO i=1,l
    o=str(1,i)
    IF (o.EQ.0) THEN
       b=MIN(INT(ranmz()*FLOAT(nb))+1,nb)
       typ=btyp(b)
       ss1=spn(bst(1,b))
       ss2=spn(bst(2,b))
       p=awgt(ss1,ss2,typ)*addp(nh)
       IF (p.GE.1.) THEN
          str(1,i)=1
          str(2,i)=b
          nh=nh+1
          nnu(typ)=nnu(typ)+1
          ele(ss1,ss2,typ)=ele(ss1,ss2,typ)+1    
       ELSEIF (ranmz().LE.p) THEN
          str(1,i)=1
          str(2,i)=b
          nh=nh+1
          nnu(typ)=nnu(typ)+1
          ele(ss1,ss2,typ)=ele(ss1,ss2,typ)+1    
       ENDIF
    ELSEIF (o.EQ.1) THEN            
       b=str(2,i)
       typ=btyp(b)
       ss1=spn(bst(1,b))
       ss2=spn(bst(2,b))
       p=dwgt(ss1,ss2,typ)*delp(nh)
       IF (p.GE.1.) THEN
          str(1,i)=0
          str(2,i)=0
          nh=nh-1
       ELSEIF (ranmz().LE.p) THEN
          str(1,i)=0
          str(2,i)=0
          nh=nh-1
       ELSE
          ele(ss1,ss2,typ)=ele(ss1,ss2,typ)+1    
       ENDIF
       nu=nu+1
       nh1=nh1+1
    ELSE
       b=str(2,i)
       typ=btyp(b)
       s1=bst(1,b)
       s2=bst(2,b)
       nh1=nh1+1
       dl=DFLOAT(nh1-last)
       ssa1=ssa1+dl*DFLOAT(ca*ca)
       sxa1=sxa1+dl*DFLOAT(ca)
       spn(s1)=-spn(s1)
       spn(s2)=-spn(s2)
       IF (b.LE.nn) THEN
          jjx=jjx+spn(s2)-spn(s1)
       ELSEIF (b.LE.2*nn) THEN
          jjy=jjy+spn(s2)-spn(s1)
       ENDIF
       ca=ca+phase(s1)*spn(s1)+phase(s2)*spn(s2)
       ssa1=ssa1+DFLOAT(ca*ca)
       sxa1=sxa1+DFLOAT(ca)
       nk=nk+1
       last=nh1
    ENDIF
 ENDDO

 dl=DFLOAT(nh1-last)
 ssa1=ssa1+0.25*DFLOAT(ca*ca)
 sxa1=sxa1+0.5*DFLOAT(ca)
 ssa1=ssa1/DFLOAT(nh1+1)
 IF (nh1.NE.0) THEN
    sxa1=(sxa1**2+DFLOAT(nh1)*ssa1)/(DFLOAT(nh1)*DFLOAT(nh1+1))
    sxa1=sxa1/DFLOAT(nn)
 ELSE
    sxa1=ssa1/DFLOAT(nn)
 ENDIF
 ssa1=ssa1/DFLOAT(nn)
 avu=avu+DFLOAT(nu)/DFLOAT(nn)
 avk=avk+DFLOAT(nk)/DFLOAT(nn)
 rhx=rhx+(DFLOAT(jjx)**2)/(DFLOAT(nn))
 rhy=rhy+(DFLOAT(jjy)**2)/(DFLOAT(nn))
 ssa=ssa+ssa1
 sxa=sxa+sxa1
 sxu=sxu+0.25*DFLOAT(cu**2)/DFLOAT(nn)
 umag=umag+0.5*DFLOAT(cu)/DFLOAT(nn)

 END SUBROUTINE dupdate
!----------------------!

!-------------------!
 SUBROUTINE linkoper

 USE blink; USE hyzer;

 IMPLICIT NONE

 INTEGER i,b,o,s1,s2,p1,p2,ss1,ss2,tt1,tt2,typ
 INTEGER last(nn),lspn(2,ll)

 DO i=1,nn
    frst(i)=0
    last(i)=0
 ENDDO
 no=0
 DO i=1,l
    o=str(1,i)
    IF (o.NE.0) THEN
       no=no+1
       b=str(2,i)
       typ=btyp(b)
       s1=bst(1,b)
       s2=bst(2,b)
       ss1=spn(s1)
       ss2=spn(s2)
       IF (o.EQ.2) THEN
          spn(s1)=-spn(s1)
          spn(s2)=-spn(s2)
       ENDIF
       tt1=spn(s1)
       tt2=spn(s2)
       lpos(no)=i
       lspn(1,no)=s1
       lspn(2,no)=s2
       lvtx(no)=legvx(ss1,ss2,tt1,tt2,typ)
       p1=last(s1)
       p2=last(s2)
       IF (p1.NE.0) THEN
          IF (lspn(1,p1).EQ.s1) THEN
             plnk(3,p1)=no
             slnk(3,p1)=1
             plnk(1,no)=p1
             slnk(1,no)=3
          ELSE
             plnk(4,p1)=no
             slnk(4,p1)=1
             plnk(1,no)=p1
             slnk(1,no)=4
          ENDIF
       ELSE
          frst(s1)=no
          fspn(s1)=1
       ENDIF
       IF (p2.NE.0) THEN
          IF (lspn(1,p2).EQ.s2) THEN
             plnk(3,p2)=no
             slnk(3,p2)=2
             plnk(2,no)=p2
             slnk(2,no)=3
          ELSE
             plnk(4,p2)=no
             slnk(4,p2)=2
             plnk(2,no)=p2
             slnk(2,no)=4
          ENDIF
       ELSE
          frst(s2)=no
          fspn(s2)=2
       ENDIF
       last(s1)=no
       last(s2)=no
    ENDIF
 END DO

 DO s1=1,nn
    p1=last(s1)
    IF (p1.NE.0) THEN
       p2=frst(s1)
       IF (s1.EQ.lspn(1,p1)) THEN
          plnk(3,p1)=p2
          IF (s1.EQ.lspn(1,p2)) THEN
             slnk(3,p1)=1
          ELSE
             slnk(3,p1)=2
          ENDIF
       ELSE
          plnk(4,p1)=p2
          IF (s1.EQ.lspn(1,p2)) THEN
             slnk(4,p1)=1
          ELSE
             slnk(4,p1)=2
          ENDIF
       ENDIF
       IF (s1.EQ.lspn(1,p2)) THEN
          plnk(1,p2)=p1
          IF (s1.EQ.lspn(1,p1)) THEN
             slnk(1,p2)=3
          ELSE
             slnk(1,p2)=4
          ENDIF
       ELSE
          plnk(2,p2)=p1
          IF (s1.EQ.lspn(1,p1)) THEN
             slnk(2,p2)=3
          ELSE
             slnk(2,p2)=4
          ENDIF
       ENDIF
    ENDIF
 ENDDO

 END SUBROUTINE linkoper
!-----------------------!

!---------------------------!
 SUBROUTINE updloop (passed)

 USE blink; USE hyzer; USE random;

 IMPLICIT NONE

 INTEGER :: i,j,k,p0,p1,p2,vx,ic,ic0,oc,nop,nop1
 REAL(8) :: r
 LOGICAL :: passed

 nop=0 
 mloop=1000*l
 DO j=1,nl
    p0=MIN(INT(ranmz()*FLOAT(no))+1,no)
    vx=lvtx(p0)
    ic0=MIN(INT(ranmz()*4.d0)+1,4)
    p1=p0
    ic=ic0
    nop1=1
    r=ranmz()
    DO k=1,4            
       IF (r.LE.vxprb(k,ic,vx)) THEN
          lvtx(p0)=vxnew(k,ic,vx)
          oc=k
          EXIT
       ENDIF
    ENDDO
    IF (lvtx(p0).NE.vx) THEN ! skip if first move is a bounce
       DO
          p2=plnk(oc,p1)
          ic=slnk(oc,p1)
          IF (p2.EQ.p0.AND.ic.EQ.ic0) EXIT
          vx=lvtx(p2)
          r=ranmz()
          DO k=1,4            
             IF (r.LE.vxprb(k,ic,vx)) THEN
                lvtx(p2)=vxnew(k,ic,vx)
                oc=k
                EXIT
             ENDIF
          ENDDO
          p1=p2
          IF (p1.EQ.p0.AND.oc.EQ.ic0) EXIT
          nop1=nop1+1
          IF (nop1.GT.mloop) THEN
             passed=.FALSE.
             RETURN
          ENDIF
       ENDDO
    ENDIF
    nop=nop+nop1
 ENDDO
 DO i=1,no
    str(1,lpos(i))=vxoper(lvtx(i))
 ENDDO
 DO i=1,nn
    IF (frst(i).NE.0) THEN
       spn(i)=vxleg(fspn(i),lvtx(frst(i)))
    ENDIF
 ENDDO
 lopers=lopers+DFLOAT(nop)
 nloops=nloops+DFLOAT(nl)
 passed=.TRUE.

 END SUBROUTINE updloop
!----------------------!

!----------------!
 SUBROUTINE adjnl

 USE hyzer;

 IMPLICIT NONE

 INTEGER :: nl1

 lopers=lopers/nloops
 nl1=1+INT(DFLOAT(2*l)/lopers)
 nl=(nl+nl1)/2
 lopers=0.d0
 nloops=0.d0

 RETURN
 END SUBROUTINE ADJNL
!--------------------!
