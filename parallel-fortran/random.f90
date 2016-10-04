!-------------!
 MODULE random

 IMPLICIT NONE

 INTEGER, SAVE :: iir=521288629,jjr=362436069,kkr=16163801,nnr=1131199299

 CONTAINS

  !---------------------!
   SUBROUTINE initran0()

   OPEN(UNIT=10,FILE='rand.in',STATUS='old')
   READ(10,*)iir
   READ(10,*)jjr
   READ(10,*)kkr
   READ(10,*)nnr
   CLOSE(10)
   iir=1+IABS(iir)
   jjr=1+IABS(jjr)
   kkr=1+IABS(kkr)
   OPEN(UNIT=10,FILE='rand.in',STATUS='unknown')
   WRITE(10,*)mzran()
   WRITE(10,*)mzran()
   WRITE(10,*)mzran()
   WRITE(10,*)mzran()
   CLOSE(10)

   END SUBROUTINE initran0
  !-----------------------!

  !-------------------------!
   SUBROUTINE initran1(seed)

   INTEGER :: seed(4)

   iir=1+IABS(seed(1))
   jjr=1+IABS(seed(2))
   kkr=1+IABS(seed(3))
   nnr=seed(4)

   END SUBROUTINE initran1
  !-----------------------!

  !------------------------!
   SUBROUTINE getseed(seed)

   INTEGER :: seed(4)

   seed(1)=mzran()
   seed(2)=mzran()
   seed(3)=mzran()
   seed(4)=mzran()

   END SUBROUTINE getseed
  !----------------------!

  !------------------------!
   INTEGER FUNCTION mzran()

  !----------------------------!
  ! G. Marsaglia and A. Zaman, !
  ! Comp. Phys. 8, 117 (1994). !
  !----------------------------!

   mzran=iir-kkr
   IF (mzran.LT.0) mzran=mzran+2147483579
   iir=jjr
   jjr=kkr
   kkr=mzran
   nnr=69069*nnr+1013904243
   mzran=mzran+nnr

   END FUNCTION mzran
  !------------------!

  !------------------------!
   REAL FUNCTION ranmz()

   ranmz=.5+.2328306E-9*mzran()

   END FUNCTION ranmz
  !------------------!

 END MODULE random
!-----------------!
