*-------------------------------------
* Sample timing routine
*       This code works on Sun and Silicon Graphics machines.
*       DOUBLE PRECISION function mysecond()
*       real arg(2)
*       mysecond = etime(arg)
*       end
* Sample timing routine
*       This code works on IBM RS/6000 machines
      DOUBLE PRECISION FUNCTION mysecond()
C     ..
C     .. External Functions ..
      INTEGER mclock
      EXTERNAL mclock
C     ..
      mysecond = mclock()*0.01D0
      END

