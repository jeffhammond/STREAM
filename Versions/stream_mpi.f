* Program: STREAM
* Programmer: John D. McCalpin, P. Vezolle
*
*-----------------------------------------------------------------------
* Copyright 1991-2003: John D. McCalpin
*-----------------------------------------------------------------------
* License:
*  1. You are free to use this program and/or to redistribute
*     this program.
*  2. You are free to modify this program for your own use,
*     including commercial use, subject to the publication
*     restrictions in item 3.
*  3. You are free to publish results obtained from running this
*     program, or from works that you derive from this program,
*     with the following limitations:
*     3a. In order to be referred to as "STREAM benchmark results",
*         published results must be in conformance to the STREAM
*         Run Rules, (briefly reviewed below) published at
*         http://www.cs.virginia.edu/stream/ref.html
*         and incorporated herein by reference.
*         As the copyright holder, John McCalpin retains the
*         right to determine conformity with the Run Rules.
*     3b. Results based on modified source code or on runs not in
*         accordance with the STREAM Run Rules must be clearly
*         labelled whenever they are published.  Examples of
*         proper labelling include:
*         "tuned STREAM benchmark results"
*         "based on a variant of the STREAM benchmark code"
*         Other comparable, clear and reasonable labelling is
*         acceptable.
*     3c. Submission of results to the STREAM benchmark web site
*         is encouraged, but not required.
*  4. Use of this program or creation of derived works based on this
*     program constitutes acceptance of these licensing restrictions.
*  5. Absolutely no warranty is expressed or implied.
*-----------------------------------------------------------------------
*                            MPI VERSION
*
*  Latest Modification: February 14, 2014 (John D. McCalpin)
*     Minor update to output formatting to prevent overflows.
*  Latest Modification: May 17, 2002 (P. Vezolle)
*  Based on STREAM version 5.0, July 30, 2000
*
* This version has been shown to work under the MPI environment
* of IBM's AIX 5.1 --- comments or suggestions on how to improve
* portability are welcome!    mailto:john@mccalpin.com
* 
* This program measures memory transfer rates in MB/s for simple
* computational kernels coded in Fortran.  
* The intent is to demonstrate the extent to which ordinary user
* code can exploit the main memory bandwidth of the system under
* test.
*=========================================================================
* The STREAM web page is at:
*          http://www.streambench.org
*
* Most of the content is currently hosted at:
*          http://www.cs.virginia.edu/stream/
*
* BRIEF INSTRUCTIONS: 
*       0) See http://www.cs.virginia.edu/stream/ref.html for details
*       1) STREAM requires a timing function called mysecond().
*          Several examples are provided in this directory.
*          "CPU" timers are only allowed for uniprocessor runs.
*          "Wall-clock" timers are required for all multiprocessor runs.
*       2) The STREAM array sizes must be set to size the test.
*          The value "N" must be chosen so that each of the three
*          arrays is at least 4x larger than the sum of all the last-
*          level caches used in the run, or 1 million elements, which-
*          ever is larger.
*          ------------------------------------------------------------
*          Note that you are free to use any array length and offset
*          that makes each array 4x larger than the last-level cache.
*          The intent is to determine the *best* sustainable bandwidth
*          available with this simple coding.  Of course, lower values
*          are usually fairly easy to obtain on cached machines, but 
*          by keeping the test to the *best* results, the answers are
*          easier to interpret.
*          You may put the arrays in common or not, at your discretion.
*          There is a commented-out COMMON statement below.
*          Fortran90 "allocatable" arrays are fine, too.
*          ------------------------------------------------------------
*       3) Compile the code with full optimization.  Many compilers
*          generate unreasonably bad code before the optimizer tightens
*          things up.  If the results are unreasonably good, on the
*          other hand, the optimizer might be too smart for me
*          Please let me know if this happens.
*       4) Mail the results to mccalpin@cs.virginia.edu
*          Be sure to include:
*               a) computer hardware model number and software revision
*               b) the compiler flags
*               c) all of the output from the test case.
*          Please let me know if you do not want your name posted along
*          with the submitted results.
*       5) See the web page for more comments about the run rules and
*          about interpretation of the results.
*
* Thanks,
*   Dr. Bandwidth
*=========================================================================
*
      PROGRAM stream_mpi
*     IMPLICIT NONE

      include 'mpif.h'

C     .. Parameters ..
      INTEGER  n,offset,ndim,ntimes
      PARAMETER (n=2000000,offset=0,ndim=n+offset,ntimes=10)

*     .. MPI Variables ..
      integer numtask, rank
      integer rc
      DOUBLE PRECISION rc1

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION scalar,t
      INTEGER j,k,nbpw,quantum
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION maxtime(4),mintime(4),avgtime(4),
     $                 times(4,ntimes)
      INTEGER bytes(4)
      CHARACTER*11 label(4)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION mysecond
      INTEGER checktick,realsize
      EXTERNAL mysecond,checktick,realsize
C     ..
C     .. Intrinsic Functions ..
C
      INTRINSIC dble,max,min,nint,sqrt
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION a(ndim),b(ndim),c(ndim)
C     ..
C     .. Common blocks ..
*     COMMON a,b,c
C     ..
C     .. Data statements ..
      DATA avgtime/4*0.0D0/,mintime/4*1.0D+36/,maxtime/4*0.0D0/
      DATA label/'Copy:      ','Scale:     ','Add:       ',
     $     'Triad:     '/
      DATA bytes/2,2,3,3/
C     ..

*       --- SETUP --- determine precision and check timing ---

*     .. MPI Initialization ..

      call MPI_INIT ( rc )
      if ( rc .ne. 0 ) then
         WRITE(*,*) ' MPI Initialization problem, error code: ',rc
         stop
      endif
      call MPI_COMM_RANK ( MPI_COMM_WORLD, rank, rc )
      call MPI_COMM_SIZE ( MPI_COMM_WORLD, numtask, rc )

      if ( rank .eq. 0 ) then
        nbpw = realsize()
        WRITE (*,*)        'Number of processors = ', numtask
        WRITE (*,FMT=9010) 'Array size = ',n
        WRITE (*,FMT=9010) 'Offset     = ',offset
        WRITE (*,FMT=9060) 'The total memory requirement is ',
     $     3.0*nbpw*n*numtask/(1024*1024),' MB',
     $     ' (',3.0*nbpw*n/(1024*1024),'MB/task)'
        WRITE (*,FMT=9030) 'You are running each test ',ntimes,' times'
        WRITE (*,FMT=9030) '--'
        WRITE (*,FMT=9030) 'The *best* time for each test is used'
        WRITE (*,FMT=9030) '*EXCLUDING* the first and last iterations'
      endif

      DO 10 j = 1,n
          a(j) = 2.0d0
          b(j) = 0.5D0
          c(j) = 0.0D0
   10 CONTINUE
      t = mysecond()
      DO 20 j = 1,n
          a(j) = 0.5d0*a(j)
   20 CONTINUE
      t = mysecond() - t
      if ( rank .eq. 0 ) then
        PRINT *,'----------------------------------------------------'
        quantum = checktick()
        WRITE (*,FMT=9000)
     $    'Your clock granularity/precision appears to be ',quantum,
     $    ' microseconds'
        PRINT *,'----------------------------------------------------'
      endif

      call MPI_BARRIER( MPI_COMM_WORLD, rc)

*       --- MAIN LOOP --- repeat test cases NTIMES times ---
      scalar = 0.5d0*a(1)
      DO 70 k = 1,ntimes

!********** COPY

        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond()
          a(1) = a(1) + t
          DO 30 j = 1,n
              c(j) = a(j)
   30     CONTINUE
        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond() - t
          c(n) = c(n) + t
          times(1,k) = t

!********** SCALE

        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond()
          c(1) = c(1) + t
          DO 40 j = 1,n
              b(j) = scalar*c(j)
   40     CONTINUE
        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond() - t
          b(n) = b(n) + t
          times(2,k) = t

!********** ADD

        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond()
          a(1) = a(1) + t
          DO 50 j = 1,n
              c(j) = a(j) + b(j)
   50     CONTINUE
        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond() - t
          c(n) = c(n) + t
          times(3,k) = t

!********** TRIAD

        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond()
          b(1) = b(1) + t
          DO 60 j = 1,n
              a(j) = b(j) + scalar*c(j)
   60     CONTINUE
        call MPI_BARRIER( MPI_COMM_WORLD, rc)
          t = mysecond() - t
          a(n) = a(n) + t
          times(4,k) = t
   70 CONTINUE

*       --- SUMMARY ---
      if ( rank .eq. 0 ) then
        DO 90 k = 2,ntimes-1
          DO 80 j = 1,4
              avgtime(j) = avgtime(j) + times(j,k)
              mintime(j) = min(mintime(j),times(j,k))
              maxtime(j) = max(maxtime(j),times(j,k))
   80     CONTINUE
   90   CONTINUE
        WRITE (*,FMT=9040)
        DO 100 j = 1,4
          avgtime(j) = avgtime(j)/dble(ntimes-2)
          rc1 = 1.*n*bytes(j)*nbpw/mintime(j)/1.0D6
          WRITE (*,FMT=9050) label(j),rc1*numtask,
     $      avgtime(j),mintime(j),maxtime(j)
  100   CONTINUE
       endif
       CALL checksums (a,b,c,n,ntimes,rank,numtask)

 1001 continue

      call MPI_FINALIZE ( rc )
      stop

 9000 FORMAT (1x,a,i6,a)
 9010 FORMAT (1x,a,i10)
 9020 FORMAT (1x,a,i4,a)
 9030 FORMAT (1x,a,i3,a,a)
 9040 FORMAT ('Function',5x,'Rate (MB/s)  Avg time   Min time  Max time'
     $       )
 9050 FORMAT (a,f12.1,3 (f11.6,2x))
 9060 FORMAT (1x,a,f9.1,a,a,f9.1,a)
      END

*-------------------------------------
* INTEGER FUNCTION dblesize()
*
* A semi-portable way to determine the precision of DOUBLE PRECISION
* in Fortran.
* Here used to guess how many bytes of storage a DOUBLE PRECISION
* number occupies.
*
      INTEGER FUNCTION realsize()
*     IMPLICIT NONE

C     .. Local Scalars ..
      DOUBLE PRECISION result,test
      INTEGER j,ndigits
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ref(30)
C     ..
C     .. External Subroutines ..
      EXTERNAL confuse
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,acos,log10,sqrt
C     ..

C       Test #1 - compare single(1.0d0+delta) to 1.0d0

   10 DO 20 j = 1,30
          ref(j) = 1.0d0 + 10.0d0** (-j)
   20 CONTINUE

      DO 30 j = 1,30
          test = ref(j)
          ndigits = j
          CALL confuse(test,result)
          IF (test.EQ.1.0D0) THEN
              GO TO 40
          END IF
   30 CONTINUE
      GO TO 50

   40 WRITE (*,FMT='(a)')
     $  '----------------------------------------------'
      WRITE (*,FMT='(1x,a,i2,a)') 'Double precision appears to have ',
     $  ndigits,' digits of accuracy'
      IF (ndigits.LE.8) THEN
          realsize = 4
      ELSE
          realsize = 8
      END IF
      WRITE (*,FMT='(1x,a,i1,a)') 'Assuming ',realsize,
     $  ' bytes per DOUBLE PRECISION word'
      WRITE (*,FMT='(a)')
     $  '----------------------------------------------'
      RETURN

   50 PRINT *,'Hmmmm.  I am unable to determine the size.'
      PRINT *,'Please enter the number of Bytes per DOUBLE PRECISION',
     $  ' number : '
      READ (*,FMT=*) realsize
      IF (realsize.NE.4 .AND. realsize.NE.8) THEN
          PRINT *,'Your answer ',realsize,' does not make sense.'
          PRINT *,'Try again.'
          PRINT *,'Please enter the number of Bytes per ',
     $      'DOUBLE PRECISION number : '
          READ (*,FMT=*) realsize
      END IF
      PRINT *,'You have manually entered a size of ',realsize,
     $  ' bytes per DOUBLE PRECISION number'
      WRITE (*,FMT='(a)')
     $  '----------------------------------------------'
      END

      SUBROUTINE confuse(q,r)
*     IMPLICIT NONE
C     .. Scalar Arguments ..
      DOUBLE PRECISION q,r
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cos
C     ..
      r = cos(q)
      RETURN
      END

* A semi-portable way to determine the clock granularity
* Adapted from a code by John Henning of Digital Equipment Corporation
*
      INTEGER FUNCTION checktick()
*     IMPLICIT NONE

C     .. Parameters ..
      INTEGER n
      PARAMETER (n=20)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION t1,t2
      INTEGER i,j,jmin
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION timesfound(n)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION mysecond
      EXTERNAL mysecond
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC max,min,nint
C     ..
      i = 0
      t1 = mysecond()

   10 t2 = mysecond()
      IF (t2.EQ.t1) GO TO 10

      t1 = t2
      i = i + 1
      timesfound(i) = t1
      IF (i.LT.n) GO TO 10

      jmin = 1000000
      DO 20 i = 2,n
          j = nint((timesfound(i)-timesfound(i-1))*1d6)
          jmin = min(jmin,max(j,0))
   20 CONTINUE

      IF (jmin.GT.0) THEN
          checktick = jmin
      ELSE
          PRINT *,'Your clock granularity appears to be less ',
     $      'than one microsecond'
          checktick = 1
      END IF
      RETURN

*      PRINT 14, timesfound(1)*1d6
*      DO 20 i=2,n
*         PRINT 14, timesfound(i)*1d6,
*     &       nint((timesfound(i)-timesfound(i-1))*1d6)
*   14    FORMAT (1X, F18.4, 1X, i8)
*   20 CONTINUE

      END


      SUBROUTINE checksums(a,b,c,n,ntimes,rank,numtask)
*     IMPLICIT NONE
      include 'mpif.h'
C     ..
C     .. Arguments ..
      DOUBLE PRECISION a(*),b(*),c(*)
      INTEGER n,ntimes
      INTEGER rank,numtask,rc, Stat(MPI_STATUS_SIZE)
      DOUBLE PRECISION Status(2), ival
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION aa,bb,cc,scalar,suma,sumb,sumc,epsilon
      INTEGER k
C     ..

      ival = 0

C     Repeat the main loop, but with scalars only.
C     This is done to check the sum & make sure all
C     iterations have been executed correctly.

      aa = 2.0D0
      bb = 0.5D0
      cc = 0.0D0
      aa = 0.5D0*aa
      scalar = 0.5d0*aa
      DO k = 1,ntimes
          cc = aa
          bb = scalar*cc
          cc = aa + bb
          aa = bb + scalar*cc
      END DO
      aa = aa*DBLE(n-2)
      bb = bb*DBLE(n-2)
      cc = cc*DBLE(n-2)

C     Now sum up the arrays, excluding the first and last
C     elements, which are modified using the timing results
C     to confuse aggressive optimizers.

      suma = 0.0d0
      sumb = 0.0d0
      sumc = 0.0d0
!$OMP PARALLEL DO REDUCTION(+:suma,sumb,sumc)
      DO 110 j = 2,n-1
          suma = suma + a(j)
          sumb = sumb + b(j)
          sumc = sumc + c(j)
  110 CONTINUE

      epsilon = 1.D-6
    
*     .. Gather results by process  0
      IF (ABS(suma-aa)/suma .GT. epsilon) THEN
          Status(1) = aa
          Status(2) = suma
      ELSE
           Status(1) = 0
      ENDIF
      if ( rank .ne. 0 ) then
          call MPI_SEND(Status, 2, MPI_DOUBLE_PRECISION, 0, rank,
     $                  MPI_COMM_WORLD, rc)
      else
         ival = ival + Status(1)
         if ( Status(1) .ne. 0 ) then
            PRINT *,'Failed Validation on array a(), Process ',rank
            PRINT *,'Target   Sum of a is = ',Status(1)
            PRINT *,'Computed Sum of a is = ',Status(2)
         endif
         do i=1,numtask-1
          call MPI_RECV(Status, 2, MPI_DOUBLE_PRECISION, i, i,
     s                  MPI_COMM_WORLD, Stat, rc);
          ival = ival + Status(1)
          if ( Status(1) .ne. 0 ) then
            PRINT *,'Failed Validation on array a(), Process ',i
            PRINT *,'Target   Sum of a is = ',Status(1)
            PRINT *,'Computed Sum of a is = ',Status(2)
          endif
         enddo
      endif

      IF (ABS(sumb-bb)/sumb .GT. epsilon) THEN
          Status(1) = bb
          Status(2) = sumb
      ELSE
           Status(1) = 0
      ENDIF
      if ( rank .ne. 0 ) then
          call MPI_SEND(Status, 2, MPI_DOUBLE_PRECISION, 0, 2*rank,
     $                  MPI_COMM_WORLD, rc)
      else
         ival = ival + Status(1)
         if ( Status(1) .ne. 0 ) then
            PRINT *,'Failed Validation on array b(), Process ',rank
            PRINT *,'Target   Sum of b is = ',Status(1)
            PRINT *,'Computed Sum of b is = ',Status(2)
         endif
         do i=1,numtask-1
          call MPI_RECV(Status, 2, MPI_DOUBLE_PRECISION, i, 2*i,
     s                  MPI_COMM_WORLD, Stat, rc);
          ival = ival + Status(1)
          if ( Status(1) .ne. 0 ) then
            PRINT *,'Failed Validation on array b(), Process ',i
            PRINT *,'Target   Sum of b is = ',Status(1)
            PRINT *,'Computed Sum of b is = ',Status(2)
          endif
         enddo
      endif

      IF (ABS(sumc-cc)/sumc .GT. epsilon) THEN
          Status(1) = cc
          Status(2) = sumc
      ELSE
           Status(1) = 0
      ENDIF
      if ( rank .ne. 0 ) then
          call MPI_SEND(Status, 2, MPI_DOUBLE_PRECISION, 0, 3*rank,
     $                  MPI_COMM_WORLD, rc)
      else
         ival = ival + Status(1)
         if ( Status(1) .ne. 0 ) then
            PRINT *,'Failed Validation on array c(), Process ',rank
            PRINT *,'Target   Sum of c is = ',Status(1)
            PRINT *,'Computed Sum of c is = ',Status(2)
         endif
         do i=1,numtask-1
          call MPI_RECV(Status, 2, MPI_DOUBLE_PRECISION, i, 3*i,
     s                  MPI_COMM_WORLD, Stat, rc);
          ival = ival + Status(1)
          if ( Status(1) .ne. 0 ) then
            PRINT *,'Failed Validation on array c(), Process ',i
            PRINT *,'Target   Sum of c is = ',Status(1)
            PRINT *,'Computed Sum of c is = ',Status(2)
          endif
         enddo
        
         if ( ival .eq. 0. ) then
           PRINT *,'-----------------------------------------------'
           PRINT*,'Solution Validates!'
           PRINT *,'-----------------------------------------------'
         endif
      endif

      END
