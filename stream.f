*=======================================================================
* Program: STREAM
* Programmer: John D. McCalpin
* RCS Revision: $Id: stream.f,v 5.6 2005/10/04 00:20:48 mccalpin Exp mccalpin $
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
* This program measures sustained memory transfer rates in MB/s for
* simple computational kernels coded in FORTRAN.
*
* The intent is to demonstrate the extent to which ordinary user
* code can exploit the main memory bandwidth of the system under
* test.
*=======================================================================
* The STREAM web page is at:
*          http://www.streambench.org
*
* Most of the content is currently hosted at:
*          http://www.cs.virginia.edu/stream/
*
* BRIEF INSTRUCTIONS:
*       0) See http://www.cs.virginia.edu/stream/ref.html for details
*       1) "CPU" timers are only allowed for uniprocessor runs.
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
      PROGRAM stream

      use, intrinsic :: iso_fortran_env, only : int64
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER n,offset,ndim,ntimes
      PARAMETER (n=20000000,offset=0,ndim=n+offset,ntimes=10)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION :: scalar
      integer(int64) :: t64, tic, toc
      integer(int64) :: tick_rate
      INTEGER j,k,nbpw
      integer(int64) :: quantum
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION maxtime(4),mintime(4),avgtime(4),
     $                 times(4,ntimes)
      INTEGER bytes(4)
      CHARACTER label(4)*11
C     ..

!$    INTEGER, external :: omp_get_num_threads
C     ..
C     .. Intrinsic Functions ..
C
      INTRINSIC dble,max,min,nint,sqrt
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION, allocatable, dimension(:) :: a, b, c
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

      allocate(a(ndim), b(ndim), c(ndim))

      call system_clock(COUNT_RATE=tick_rate)
C     set timing to max precision, typically sub-microsecond

      nbpw = storage_size(a)/8

      PRINT *,'----------------------------------------------'
      PRINT *,'STREAM Version $Revision: 5.6 $'
      PRINT *,'----------------------------------------------'
      WRITE (*,FMT=9010) 'Array size = ',n
      WRITE (*,FMT=9010) 'Offset     = ',offset
      WRITE (*,FMT=9020) 'The total memory requirement is ',
     $  3*nbpw*n/ (1024*1024),' MB'
      WRITE (*,FMT=9030) 'You are running each test ',ntimes,' times'
      WRITE (*,FMT=9030) '--'
      WRITE (*,FMT=9030) 'The *best* time for each test is used'
      WRITE (*,FMT=9030) '*EXCLUDING* the first and last iterations'

!$OMP PARALLEL
!$OMP MASTER
      PRINT *,'----------------------------------------------'
!$    PRINT *,'Number of Threads = ',OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

      PRINT *,'----------------------------------------------'
!$OMP PARALLEL
      PRINT *,'Printing one line per active thread....'
!$OMP END PARALLEL

!$OMP PARALLEL DO
      DO 10 j = 1,n
          a(j) = 2.0d0
          b(j) = 0.5D0
          c(j) = 0.0D0
   10 CONTINUE
      call system_clock(count=tic)
!$OMP PARALLEL DO
      DO 20 j = 1,n
          a(j) = 0.5d0*a(j)
   20 CONTINUE
      call system_clock(count=toc)
      t64 = toc - tic
      PRINT *,'----------------------------------------------------'

      print '(a,f10.3)','Clock granularity/precision (microseconds):',
     &   1/dble(tick_rate) * 1e6
      PRINT *,'----------------------------------------------------'

*       --- MAIN LOOP --- repeat test cases NTIMES times ---
      scalar = 0.5d0*a(1)
      DO 70 k = 1,ntimes

          call system_clock(count=tic)
          a(1) = a(1) + tic
!$OMP PARALLEL DO
          DO 30 j = 1,n
              c(j) = a(j)
   30     CONTINUE
          call system_clock(count=toc)
          t64 = toc - tic
          c(n) = c(n) + t64
          times(1,k) = t64 / dble(tick_rate)

          call system_clock(count=tic)
          c(1) = c(1) + tic
!$OMP PARALLEL DO
          DO 40 j = 1,n
              b(j) = scalar*c(j)
   40     CONTINUE
          call system_clock(count=toc)
          t64 = toc - tic
          b(n) = b(n) + t64
          times(2,k) = t64 / dble(tick_rate)

          call system_clock(count=tic)
          a(1) = a(1) + tic
!$OMP PARALLEL DO
          DO 50 j = 1,n
              c(j) = a(j) + b(j)
   50     CONTINUE
          call system_clock(count=toc)
          t64 = toc - tic
          c(n) = c(n) + t64
          times(3,k) = t64 / dble(tick_rate)

          call system_clock(count=tic)
          b(1) = b(1) + tic
!$OMP PARALLEL DO
          DO 60 j = 1,n
              a(j) = b(j) + scalar*c(j)
   60     CONTINUE
          call system_clock(count=toc)
          t64 = toc - tic
          a(n) = a(n) + t64
          times(4,k) = t64 / dble(tick_rate)
   70 CONTINUE

*       --- SUMMARY ---
      DO 90 k = 2,ntimes
          DO 80 j = 1,4
              avgtime(j) = avgtime(j) + times(j,k)
              mintime(j) = min(mintime(j),times(j,k))
              maxtime(j) = max(maxtime(j),times(j,k))
   80     CONTINUE
   90 CONTINUE
      WRITE (*,FMT=9040)
      DO 100 j = 1,4
          avgtime(j) = avgtime(j)/dble(ntimes-1)
          WRITE (*,FMT=9050) label(j),n*bytes(j)*nbpw/mintime(j)/1.0D6,
     $      avgtime(j),mintime(j),maxtime(j)
  100 CONTINUE
      PRINT *,'----------------------------------------------------'
      CALL checksums (a,b,c,n,ntimes)
      PRINT *,'----------------------------------------------------'

 9010 FORMAT (1x,a,i10)
 9020 FORMAT (1x,a,i4,a)
 9030 FORMAT (1x,a,i3,a,a)
 9040 FORMAT ('Function',5x,'Rate (MB/s)  Avg time   Min time  Max time'
     $       )
 9050 FORMAT (a,4 (f10.4,2x))

      contains


      SUBROUTINE checksums(a,b,c,n,ntimes)
*     IMPLICIT NONE
C     ..
C     .. Arguments ..
      DOUBLE PRECISION a(*),b(*),c(*)
      INTEGER n,ntimes
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION aa,bb,cc,scalar,suma,sumb,sumc,epsilon
      INTEGER k
C     ..

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

      IF (ABS(suma-aa)/suma .GT. epsilon) THEN
          PRINT *,'Failed Validation on array a()'
          PRINT *,'Target   Sum of a is = ',aa
          PRINT *,'Computed Sum of a is = ',suma
      ELSEIF (ABS(sumb-bb)/sumb .GT. epsilon) THEN
          PRINT *,'Failed Validation on array b()'
          PRINT *,'Target   Sum of b is = ',bb
          PRINT *,'Computed Sum of b is = ',sumb
      ELSEIF (ABS(sumc-cc)/sumc .GT. epsilon) THEN
          PRINT *,'Failed Validation on array c()'
          PRINT *,'Target   Sum of c is = ',cc
          PRINT *,'Computed Sum of c is = ',sumc
      ELSE
          PRINT *,'Solution Validates!'
      ENDIF

      END

      END program stream
