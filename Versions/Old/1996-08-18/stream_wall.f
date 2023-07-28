* Program: Stream
* Programmer: John D. McCalpin
* Revision: 3.0, August 2, 1994
*
* This program measures memory transfer rates in MB/s for simple
* computational kernels coded in Fortran.  These numbers reveal the
* quality of code generation for simple uncacheable kernels as well
* as showing the cost of floating-point operations relative to memory
* accesses.
*
* INSTRUCTIONS:
*       1) Stream requires a cpu timing function called myclock().
*          A sample is in myclock.c.  This is unfortunately rather
*          system dependent.  It helps to know the granularity of the
*          timing.  The code below assumes that the granularity is
*          1/100 seconds.
*       2) Stream requires a good bit of memory to run.
*          Adjust the Parameter 'N' in the second line of the main
*          program to give a 'timing calibration' of at least 20 clicks.
*          This will provide rate estimates that should be good to
*          about 5% precision.
*       3) Compile the code with full optimization.  Many compilers
*          generate unreasonably bad code before the optimizer tightens
*          things up.  If the results are unreasonable good, on the
*          other hand, the optimizer might be too smart for me!
*       4) Mail the results to mccalpin@cs.virginia.edu
*          Be sure to include:
*               a) computer hardware model number and software revision
*               b) the compiler flags
*               c) all of the output from the test case.
*
* Thanks!
*
      PROGRAM stream
C     .. Parameters ..
      INTEGER n,ntimes
      PARAMETER (n=2 000 000,ntimes=10)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION t,t0
      INTEGER j,k,nbpw
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION a(n),b(n),c(n),maxtime(4),mintime(4),rmstime(4),
     $                 times(4,ntimes),sum(3)
      INTEGER bytes(4)
      CHARACTER label(4)*11
C     ..
C     .. External Functions ..
      DOUBLE PRECISION myclock
      INTEGER realsize
      EXTERNAL myclock,realsize
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,max,min,sqrt
C     ..
C     .. Data statements ..
      DATA rmstime/4*0.0D0/,mintime/4*1.0D+36/,maxtime/4*0.0D0/
      DATA label/'Assignment:','Scaling   :','Summing   :',
     $     'SAXPYing  :'/
      DATA bytes/2,2,3,3/
C     ..

*       --- SETUP --- determine precision and check timing ---

      nbpw = realsize()

      t = myclock(t0)
      t = myclock(t0)
      DO 10 j = 1,n
          a(j) = 1.0D0
          b(j) = 2.0D0
          c(j) = 0.0D0
   10 CONTINUE
      t = myclock(t0) - t
      PRINT *,'Timing calibration ; time = ',t*100,' hundredths',
     $  ' of a second'
      PRINT *,'Increase the size of the arrays if this is <30 ',
     $  ' and your clock precision is =<1/100 second'
      PRINT *,'---------------------------------------------------'

*       --- MAIN LOOP --- repeat test cases NTIMES times ---
      DO 60 k = 1,ntimes

          t = myclock(t0)
          DO 20 j = 1,n
              c(j) = a(j)
   20     CONTINUE
          t = myclock(t0) - t
          times(1,k) = t

          t = myclock(t0)
          DO 30 j = 1,n
              b(j) = 3.0D0*c(j)
   30     CONTINUE
          t = myclock(t0) - t
          times(2,k) = t

          t = myclock(t0)
          DO 40 j = 1,n
              c(j) = a(j) + b(j)
   40     CONTINUE
          t = myclock(t0) - t
          times(3,k) = t

          t = myclock(t0)
          DO 50 j = 1,n
              a(j) = b(j) + 3.0D0*c(j)
   50     CONTINUE
          t = myclock(t0) - t
          times(4,k) = t
   60 CONTINUE

*       --- SUMMARY ---
      DO 80 k = 1,ntimes
          DO 70 j = 1,4
              rmstime(j) = rmstime(j) + times(j,k)**2
              mintime(j) = min(mintime(j),times(j,k))
              maxtime(j) = max(maxtime(j),times(j,k))
   70     CONTINUE
   80 CONTINUE
      WRITE (*,FMT=9000)
      DO 90 j = 1,4
          rmstime(j) = sqrt(rmstime(j)/dble(ntimes))
          WRITE (*,FMT=9010) label(j),n*bytes(j)*nbpw/mintime(j)/1.0D6,
     $      rmstime(j),mintime(j),maxtime(j)
   90 CONTINUE
      sum(1) = 0.0
      sum(2) = 0.0
      sum(3) = 0.0
      DO 100 j=1,n
         sum(1) = sum(1) + a(j)
         sum(2) = sum(2) + b(j)
         sum(3) = sum(3) + c(j)
  100 CONTINUE
      PRINT *,'Sum of a is : ',sum(1)
      PRINT *,'Sum of b is : ',sum(2)
      PRINT *,'Sum of c is : ',sum(3)

 9000 FORMAT ('Function',5x,'Rate (MB/s)  RMS time   Min time  Max time'
     $       )
 9010 FORMAT (a,4 (f10.4,2x))
      END

*-------------------------------------
* INTEGER FUNCTION dblesize()
*
* A semi-portable way to determine the precision of DOUBLEPRECISION
* in Fortran.
* Here used to guess how many bytes of storage a DOUBLEPRECISION 
* number occupies.
*
      INTEGER FUNCTION realsize()

C     .. Local Scalars ..
      DOUBLE PRECISION result,test
      INTEGER j,ndigits
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ref(30)
C     ..
C     .. External Subroutines ..
      EXTERNAL dummy
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,acos,log10,sqrt
C     ..

C       Test #1 - compare single(1.0d0+delta) to 1.0d0

   10 DO 20 j = 1,30
          ref(j) = 1.0d0 + 10.0d0**(-j)
   20 CONTINUE

      DO 30 j = 1,30
          test = ref(j)
          ndigits = j
          CALL dummy(test,result)
          IF (test.EQ.1.0D0) THEN
              GO TO 50
          END IF
   30 CONTINUE
      GOTO 60

   50 WRITE (*,FMT='(a)') '--------------------------------------'
      WRITE (*,FMT='(1x,a,i2,a)') 'Double precision appears to have ',
     $  ndigits,' digits of accuracy'
      IF (ndigits.LE.8) THEN
          realsize = 4
      ELSE
          realsize = 8
      END IF
      WRITE (*,FMT='(1x,a,i1,a)') 'Assuming ',realsize,
     $  ' bytes per DOUBLEPRECISION word'
      WRITE (*,FMT='(a)') '--------------------------------------'
      RETURN

   60 PRINT *,'Hmmmm.  I am unable to determine the size of a REAL'
      PRINT *,'Please enter the number of Bytes per DOUBLEPRECISION',
     $  ' number : '
      READ (*,FMT=*) realsize
      IF (realsize.NE.4 .AND. realsize.NE.8) THEN
          PRINT *,'Your answer ',realsize,' does not make sense!'
          PRINT *,'Try again!'
          PRINT *,'Please enter the number of Bytes per ',
     $      'REAL number : '
          READ (*,FMT=*) realsize
      END IF
      PRINT *,'You have manually entered a size of ',realsize,
     $  ' bytes per REAL number'
      WRITE (*,FMT='(a)') '--------------------------------------'
      END

      SUBROUTINE dummy(q,r)
C     .. Scalar Arguments ..
      DOUBLE PRECISION q,r
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cos
C     ..
      r = cos(q)
      RETURN
      END
