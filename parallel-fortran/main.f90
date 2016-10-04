!------------!
 MODULE hyzer

 IMPLICIT NONE

 INTEGER,PARAMETER :: nx=6
 INTEGER,PARAMETER :: ny=nx
 INTEGER,PARAMETER :: nn=nx*ny
 INTEGER,PARAMETER :: nb=5*nn+nn/2
 INTEGER,PARAMETER :: ll=6000000
 INTEGER,PARAMETER :: nvx=24
 INTEGER,PARAMETER :: maxp=100

 REAL*8 jz(4),jx(4)
 REAL*8 beta
 REAL*8 hz(maxp)
 INTEGER nlvv(maxp)

 INTEGER l
 INTEGER mloop
 INTEGER nh
 INTEGER str(2,ll)
 INTEGER spn(nn)
 INTEGER ph(0:nn)
 INTEGER bst(2,nb)
 INTEGER btyp(nb)
 REAL*8 awgt(-1:1,-1:1,4)
 REAL*8 addp(0:ll)
 REAL*8 dwgt(-1:1,-1:1,4)
 REAL*8 delp(0:ll)
 REAL*8 amax(4,maxp)


 INTEGER vxoper(nvx)
 INTEGER legvx(-1:1,-1:1,-1:1,-1:1,4)
 INTEGER vxleg(4,nvx)
 INTEGER vxnew(4,4,nvx)
 REAL*8 vxprb(4,4,nvx)

 INTEGER nl
 REAL*8 nloops
 REAL*8 lopers

 INTEGER phase(nn)

 END MODULE hyzer
!----------------!

!-----------!
 MODULE bacc

 USE hyzer, ONLY: maxp

 IMPLICIT NONE

 REAL(8),SAVE :: acc(maxp)

 END MODULE bacc
!---------------!

!-------------!
 MODULE bmtemp

 USE hyzer, ONLY: maxp,nn

 IMPLICIT NONE

 REAL(8),SAVE :: avu(maxp),avk(maxp)
 REAL(8),SAVE :: rhx(maxp),rhy(maxp)
 REAL(8),SAVE :: ssa(maxp),sxa(maxp),sxu(maxp),umag(maxp)

 END MODULE bmtemp
!-----------------!

!-----------!
 MODULE bwgt

 USE hyzer, ONLY: maxp

 IMPLICIT NONE

 REAL(8),SAVE :: wrat(-1:1,-1:1,4,maxp)

 END MODULE bwgt
!---------------!

!------------!
 PROGRAM main

 USE hyzer; USE mpi; USE random;

 IMPLICIT NONE

 INTEGER i,j,k,id,ie,lm,s1,s2,np,equ
 INTEGER tfr,pnum,istep,mstep,nruns
 INTEGER mord(maxp),imord(maxp)
 INTEGER stat(MPI_STATUS_SIZE),ipar(9),seed(4),arr(2)
 REAL*8  rpar(11),h0,dh,del(4),jjj(4)

 CALL MPI_INIT(ie)
 IF (ie.NE.MPI_SUCCESS) THEN
    CALL openlog
    WRITE(12,*)'MPI_INIT failure '
    CALL closelog
    STOP
 ENDIF
 CALL MPI_COMM_RANK(MPI_COMM_WORLD,id,ie)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,ie)
 np=np-1

 IF (id.EQ.0) THEN         
    CALL openlog
    WRITE(12,*)'Number of parameters ',np
    CALL closelog
    OPEN (UNIT=10,FILE='read.in',STATUS='old')
    READ(10,*)del(1)
    READ(10,*)del(2)
    READ(10,*)del(3)
    READ(10,*)del(4)
    READ(10,*)jjj(1)
    READ(10,*)jjj(2)
    READ(10,*)jjj(3)
    READ(10,*)jjj(4)
    READ(10,*)h0
    READ(10,*)dh
    READ(10,*)beta
    READ(10,*)istep
    READ(10,*)nruns
    READ(10,*)mstep
    READ(10,*)equ
    READ(10,*)tfr
    CLOSE(10)
    DO i=1,4
       jx(i)=del(i)*ABS(jjj(i))
       jz(i)=jjj(i)
    ENDDO
    ipar(1)=np
    ipar(2)=istep
    ipar(3)=nruns
    ipar(4)=mstep
    ipar(5)=tfr
    CALL initran0
    IF (equ.EQ.0) THEN
       CALL initconf
    ELSE            
       OPEN (UNIT=20,FILE='conf',STATUS='old')
    ENDIF
    DO i=1,np
       IF (equ.NE.0) CALL readconf
       hz(i)=h0+DFLOAT(i-1)*dh
       mord(i)=i
       ipar(6)=l
       ipar(7)=nh
       ipar(8)=nl
       ipar(9)=i
       rpar(1)=beta
       rpar(2)=jx(1)
       rpar(3)=jx(2)
       rpar(4)=jx(3)
       rpar(5)=jx(4)
       rpar(6)=jz(1)
       rpar(7)=jz(2)
       rpar(8)=jz(3)
       rpar(9)=jz(4)
       rpar(10)=h0
       rpar(11)=dh
       nlvv(i)=nl
       CALL getseed(seed)
       CALL MPI_SSEND(ipar,9,MPI_INTEGER,i,1,MPI_COMM_WORLD,ie)
       CALL MPI_SSEND(rpar,11,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,ie)
       CALL MPI_SSEND(spn,nn,MPI_INTEGER,i,3,MPI_COMM_WORLD,ie)
       CALL MPI_SSEND(str,2*ll,MPI_INTEGER,i,4,MPI_COMM_WORLD,ie)
       CALL MPI_SSEND(seed,4,MPI_INTEGER,i,5,MPI_COMM_WORLD,ie)
       CALL openlog
       WRITE(12,*)'Sent configuration data to ',i,nh,nl
       CALL closelog
       CALL MPI_RECV(arr,2,MPI_INTEGER,i,1,MPI_COMM_WORLD,stat,ie)
       CALL openlog
       WRITE(12,*)'Configuration received by ',i,arr(1),arr(2)
       CALL closelog
    ENDDO
    IF (equ.NE.0) CLOSE(20)
    CALL initwgt(np)
    IF (istep.NE.0) THEN
       CALL openlog
       WRITE(12,*)'Starting Equilibration'
       CALL closelog
       DO j=1,istep/tfr
          CALL tempering(0,np,mord,imord)
       ENDDO
       lm=l
       DO i=1,np
          CALL recveconf(imord(i))
          IF (l.GT.lm) lm=l
          nlvv(i)=nl
          CALL openlog
          WRITE(12,*)i,imord(i),' equilibrated;  L,nl = ',l,nl
          CALL closelog
       ENDDO
       DO i=1,np
          CALL MPI_SSEND(lm,1,MPI_INTEGER,i,9,MPI_COMM_WORLD,ie)
       ENDDO
       CALL openlog
       WRITE(12,*)'Completed equilibration, L = ',lm
       CALL closelog
       l=lm
    ENDIF
    DO k=1,nruns
       CALL zerodat0
       CALL openlog
       WRITE(12,*)'Starting run ',k
       CALL closelog
       DO j=1,mstep/tfr
          CALL tempering(1,np,mord,imord)
       ENDDO
       IF (MOD(k,2).EQ.0) THEN
          OPEN(UNIT=20,FILE='conf',STATUS='unknown')
       ELSE
          OPEN(UNIT=20,FILE='conf2',STATUS='unknown')
       ENDIF
       DO i=1,np
          CALL recveconf(imord(i))
!          CALL openlog
!          WRITE(12,*)'Received configuration ',i,imord(i),nl
!          CALL closelog
          CALL writeconf
       ENDDO
       CLOSE(20)
       CALL writeres(np,mstep/tfr)
       CALL openlog
       WRITE(12,*)'Completed run ',k
       CALL closelog
    ENDDO
 ELSE
    CALL MPI_RECV(ipar,9,MPI_INTEGER,0,1,MPI_COMM_WORLD,stat,ie)
    CALL MPI_RECV(rpar,11,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,stat,ie)
    CALL MPI_RECV(spn,nn,MPI_INTEGER,0,3,MPI_COMM_WORLD,stat,ie)
    CALL MPI_RECV(str,2*ll,MPI_INTEGER,0,4,MPI_COMM_WORLD,stat,ie)
    CALL MPI_RECV(seed,4,MPI_INTEGER,0,5,MPI_COMM_WORLD,stat,ie) 
    np=ipar(1)
    istep=ipar(2)
    nruns=ipar(3)
    mstep=ipar(4)
    tfr=ipar(5)
    l=ipar(6)
    nh=ipar(7)
    nl=ipar(8)
    pnum=ipar(9)
    beta=rpar(1)
    jx(1)=rpar(2)
    jx(2)=rpar(3)
    jx(3)=rpar(4)
    jx(4)=rpar(5)
    jz(1)=rpar(6)
    jz(2)=rpar(7)
    jz(3)=rpar(8)
    jz(4)=rpar(9)
    h0=rpar(10)
    dh=rpar(11)
    DO i=1,np
       hz(i)=h0+DFLOAT(i-1)*dh
    ENDDO
    arr(1)=nh
    arr(2)=nl
    CALL MPI_SSEND(arr,2,MPI_INTEGER,0,1,MPI_COMM_WORLD,ie)
    CALL initran1(seed)
    CALL simulation(istep,nruns,mstep,tfr,np,pnum)
 ENDIF

 CALL MPI_FINALIZE(ie)

 END PROGRAM main
!----------------!

!-----------------------------------------!
 SUBROUTINE tempering (ires,np,mord,imord)

 USE bacc; USE hyzer; USE mpi; USE random; USE bwgt;

 IMPLICIT NONE

 INTEGER i,j,k,s1,s2,i1,i2,v1,v2,np,ie,ires,mord(maxp),imord(maxp)
 INTEGER elem(-1:1,-1:1,4,maxp),ele(-1:1,-1:1,4)
 INTEGER stat(MPI_STATUS_SIZE),iarr(2)
 REAL*8 res(8),rarr(2),nloopa(maxp),lopera(maxp),rat

 DO i=1,np
    IF (ires.EQ.0) THEN
       CALL MPI_RECV(nlvv(mord(i)),1,MPI_INTEGER,i,2,MPI_COMM_WORLD,stat,ie)
       CALL MPI_RECV(rarr,2,MPI_DOUBLE_PRECISION,i,3,MPI_COMM_WORLD,stat,ie)
       nloopa(mord(i))=rarr(1)
       lopera(mord(i))=rarr(2)
    ELSE
       CALL MPI_RECV(res,8,MPI_DOUBLE_PRECISION,i,2,MPI_COMM_WORLD,stat,ie)
       CALL addres(mord(i),res)
    ENDIF
    CALL MPI_RECV(ele,36,MPI_INTEGER,i,4,MPI_COMM_WORLD,stat,ie)
    DO j=1,4
    DO s2=-1,1,2
    DO s1=-1,1,2
       elem(s1,s2,j,i)=ele(s1,s2,j)
    ENDDO
    ENDDO
    ENDDO
 ENDDO
 DO i=1,np
    imord(mord(i))=i
 ENDDO
 DO j=1,1
    DO v1=1,np-1
       v2=v1+1
       i1=imord(v1)
       i2=imord(v2)
       rat=1.d0
       DO s2=-1,1,2
       DO s1=-1,1,2
       DO k=1,4
          rat=rat*wrat(s1,s2,k,v1)**(elem(s1,s2,k,i2)-elem(s1,s2,k,i1))
       ENDDO
       ENDDO
       ENDDO
       IF (ranmz().LT.rat) THEN
          imord(v1)=i2
          imord(v2)=i1
          acc(v1)=acc(v1)+1.d0
       ENDIF
    ENDDO
 ENDDO
 DO i=1,np
    mord(imord(i))=i
 ENDDO
 DO i=1,np
    iarr(1)=mord(i)
    iarr(2)=nlvv(mord(i))
    CALL MPI_SSEND(iarr,2,MPI_INTEGER,i,7,MPI_COMM_WORLD,ie)
    IF (ires.EQ.0) THEN
       rarr(1)=nloopa(mord(i))
       rarr(2)=lopera(mord(i))
       CALL MPI_SSEND(rarr,2,MPI_DOUBLE_PRECISION,i,8,MPI_COMM_WORLD,ie)
    ENDIF
 ENDDO

 END SUBROUTINE tempering
!------------------------!

!---------------------------!
 SUBROUTINE recveconf (pnum)

 USE hyzer; USE mpi;

 IMPLICIT NONE

 INTEGER ie,pnum,arr(3),stat(MPI_STATUS_SIZE)

 CALL MPI_RECV(arr,3,MPI_INTEGER,pnum,5,MPI_COMM_WORLD,stat,ie)
 CALL MPI_RECV(spn,nn,MPI_INTEGER,pnum,6,MPI_COMM_WORLD,stat,ie)
 CALL MPI_RECV(str,2*ll,MPI_INTEGER,pnum,7,MPI_COMM_WORLD,stat,ie)               
 l=arr(1)
 nh=arr(2)
 nl=arr(3)

 END SUBROUTINE recveconf
!------------------------!

!-------------------------!
 SUBROUTINE addres (k,res)

 USE bmtemp;

 IMPLICIT NONE

 INTEGER k
 REAL*8 res(8)

 avu(k)=avu(k)+res(1)
 avk(k)=avk(k)+res(2)
 umag(k)=umag(k)+res(3)
 sxu(k)=sxu(k)+res(4)
 ssa(k)=ssa(k)+res(5)
 sxa(k)=sxa(k)+res(6)
 rhx(k)=rhx(k)+res(7)
 rhy(k)=rhy(k)+res(8)

 END SUBROUTINE addres
!---------------------!

!-----------------------------!
 SUBROUTINE writeres (np,nmsr)

 USE bacc; USE bmtemp; USE hyzer;

 IMPLICIT NONE
      
 INTEGER i,np,nmsr

 DO i=1,np
    umag(i)=umag(i)/DFLOAT(nmsr)      
    avu(i)=-avu(i)/(DBLE(beta)*DFLOAT(nmsr))
    avu(i)=avu(i)+amax(1,i)+amax(2,i)+amax(3,i)+amax(4,i)
    avk(i)=-avk(i)/(DBLE(beta)*DFLOAT(nmsr))
    sxu(i)=sxu(i)/DFLOAT(nmsr)
    ssa(i)=ssa(i)/DFLOAT(nmsr)
    sxa(i)=sxa(i)/DFLOAT(nmsr)
    rhx(i)=rhx(i)/(DBLE(beta)*DFLOAT(nmsr))
    rhy(i)=rhy(i)/(DBLE(beta)*DFLOAT(nmsr))
    acc(i)=acc(i)/(DFLOAT(nmsr))
 ENDDO

 OPEN(UNIT=10,FILE='enr.dat',STATUS='unknown',ACCESS='append')
 DO i=1,np
    WRITE(10,*)i,avu(i),avk(i)
 ENDDO
 CLOSE(10)

 OPEN(UNIT=10,FILE='uni.dat',STATUS='unknown',ACCESS='append')
 DO i=1,np
    WRITE(10,*)i,sxu(i),umag(i)
 ENDDO
 CLOSE(10)

 OPEN(UNIT=10,FILE='stg.dat',STATUS='unknown',ACCESS='append')
 DO i=1,np
    WRITE(10,*)i,ssa(i),sxa(i)
 ENDDO
 CLOSE(10)

 OPEN(UNIT=10,FILE='rho.dat',STATUS='unknown',ACCESS='append')
 DO i=1,np
    WRITE(10,*)i,rhx(i),rhy(i)
 ENDDO
 CLOSE(10)

 OPEN(UNIT=10,FILE='acc.dat',STATUS='unknown',ACCESS='append')
 DO i=1,np-1
    WRITE(10,*)i,acc(i)
 ENDDO
 CLOSE(10)

 END SUBROUTINE writeres
!-----------------------!

!-------------------!
 SUBROUTINE zerodat0

 USE bacc; USE bmtemp;

 IMPLICIT NONE

 INTEGER i

 DO i=1,maxp
    avu(i)=0.d0
    avk(i)=0.d0
    umag(i)=0.d0
    sxu(i)=0.d0
    ssa(i)=0.d0
    sxa(i)=0.d0
    rhx(i)=0.d0
    rhy(i)=0.d0
    acc(i)=0.d0
 ENDDO

 END SUBROUTINE zerodat0
!-----------------------!

!-----------------------!
 SUBROUTINE initwgt (np)

 USE hyzer; USE bwgt;

 IMPLICIT NONE

 INTEGER i,j,s1,s2,np
 REAL*8 wgt(-1:1,-1:1,4,maxp),z

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
 DO i=1,np-1
    DO j=1,4
    DO s2=-1,1,2
    DO s1=-1,1,2
       wrat(s1,s2,j,i)=wgt(s1,s2,j,i)/wgt(s1,s2,j,i+1)
    ENDDO
    ENDDO
    ENDDO
 ENDDO

 END SUBROUTINE initwgt
!----------------------!

!-------------------!
 SUBROUTINE initconf

 USE hyzer;

 IMPLICIT NONE

 INTEGER i

 l=88
 nh=0
 nl=5
 DO i=1,nn
    spn(i)=(-1)**i
 ENDDO
 DO i=1,ll
    str(1,i)=0
    str(2,i)=0
 ENDDO

 END SUBROUTINE initconf
!-----------------------!

!-------------------!
 SUBROUTINE readconf

 USE hyzer;

 IMPLICIT NONE

 INTEGER i

 READ(20,*)l,nh,nl
 DO i=1,nn
    READ(20,*)spn(i)
 ENDDO
 DO i=1,l
    READ(20,*)str(1,i),str(2,i)
 ENDDO
 DO i=l+1,ll
    str(1,i)=0
    str(2,i)=0
 ENDDO

 END SUBROUTINE readconf
!-----------------------!

!--------------------!
 SUBROUTINE readconf2

 USE hyzer;

 IMPLICIT NONE

 INTEGER i

 READ(20,*)l,nh,nl
 DO i=1,nn
    READ(20,*)spn(i)
 ENDDO
 DO i=1,l
    READ(20,*)str(1,i),str(2,i)
    str(1,i+l)=str(1,i)
    str(2,i+l)=str(2,i)
 ENDDO
 l=2*l
 nh=2*nh
 DO i=l+1,ll
    str(1,i)=0
    str(2,i)=0
 ENDDO

 END SUBROUTINE readconf2
!------------------------!

!--------------------!
 SUBROUTINE writeconf

 USE hyzer;

 IMPLICIT NONE

 INTEGER i

 WRITE(20,*)l,nh,nl
 DO i=1,nn
    WRITE(20,*)spn(i)
 ENDDO
 DO i=1,l
    WRITE(20,*)str(1,i),str(2,i)
 ENDDO

 END SUBROUTINE writeconf
!------------------------!

!------------------!
 SUBROUTINE openlog

 IMPLICIT NONE

 OPEN(UNIT=12,FILE='log.txt',STATUS='unknown',ACCESS='append')

 END SUBROUTINE openlog
!----------------------!

!-------------------!
 SUBROUTINE closelog

 IMPLICIT NONE

 CLOSE(12)

 END SUBROUTINE closelog
!-----------------------!
