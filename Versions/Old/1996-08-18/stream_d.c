/*
* Program: Stream
* Programmer: John D. McCalpin
* Revision: 2.1, August 30, 1995
*
* This program measures memory transfer rates in MB/s for simple 
* computational kernels coded in Fortran.  These numbers reveal the
* quality of code generation for simple uncacheable kernels as well
* as showing the cost of floating-point operations relative to memory
* accesses.
*
* INSTRUCTIONS:
*	1) (fortran-specific, omitted.)
*	2) Stream requires a good bit of memory to run.
*	   Adjust the Parameter 'N' in the second line of the main
*	   program to give a 'timing calibration' of at least 20 clicks.
*	   This will provide rate estimates that should be good to 
*	   about 5% precision.
*	3) Compile the code with full optimization.  Many compilers
*	   generate unreasonably bad code before the optimizer tightens
*	   things up.  If the results are unreasonable good, on the
*	   other hand, the optimizer might be too smart for me!
*	4) Mail the results to mccalpin@cs.virginia.edu
*	   Be sure to include:
*		a) computer hardware model number and software revision
*		b) the compiler flags
*		c) all of the output from the test case.
* Thanks!
*
* this version was ported from fortran to c by mark hahn, hahn+@pitt.edu.
*/

#define N 1000000
#define NTIMES 10

#ifdef __hpux
#define _HPUX_SOURCE 1
#else
#define _INCLUDE_POSIX_SOURCE 1
#endif
#include <limits.h>
#include <sys/time.h>
#include <math.h>
#include <stdio.h>

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

struct timeval tvStart;

void utimeStart() {
    struct timezone tz;
    gettimeofday(&tvStart,&tz);
}

float utime() {
    struct timeval tv;
    struct timezone tz;
    float utime;
    gettimeofday(&tv,&tz);
    utime = 1e6 * (tv.tv_sec - tvStart.tv_sec) + tv.tv_usec - tvStart.tv_usec;
    return utime;
}

typedef double real;
static real a[N],b[N],c[N];

int main() {
    int j,k;
    float times[4][NTIMES];
    static float rmstime[4] = {0};
    static float mintime[4] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};
    static float maxtime[4] = {0};
    static char *label[4] = {"Assignment:",
			     "Scaling   :",
			     "Summing   :",
			     "SAXPYing  :"};
    static float bytes[4] = { 2 * sizeof(real) * N,
			      2 * sizeof(real) * N,
			      3 * sizeof(real) * N,
			      3 * sizeof(real) * N};

    /* --- SETUP --- determine precision and check timing --- */
    utimeStart();
    for (j=0; j<N; j++) {
	a[j] = 1.0;
	b[j] = 2.0;
	c[j] = 0.0;
    }
    printf("Timing calibration ; time = %f usec.\n",utime());
    printf("Increase the size of the arrays if this is < 300000\n"
	   "and your clock precision is =< 1/100 second.\n");
    printf("---------------------------------------------------\n");
    
    /*	--- MAIN LOOP --- repeat test cases NTIMES times --- */
    for (k=0; k<NTIMES; k++) {
	utimeStart();
	for (j=0; j<N; j++)
	    c[j] = a[j];
	times[0][k] = utime();
	
	utimeStart();
	for (j=0; j<N; j++)
	    c[j] = 3.0e0*a[j];
	times[1][k] = utime();
	
	utimeStart();
	for (j=0; j<N; j++)
	    c[j] = a[j]+b[j];
	times[2][k] = utime();
	
	utimeStart();
	for (j=0; j<N; j++)
	    c[j] = a[j]+3.0e0*b[j];
	times[3][k] = utime();
    }
    
    /*	--- SUMMARY --- */
    for (k=0; k<NTIMES; k++) {
	for (j=0; j<4; j++) {
	    rmstime[j] = rmstime[j] + (times[j][k] * times[j][k]);
	    mintime[j] = MIN(mintime[j], times[j][k]);
	    maxtime[j] = MAX(maxtime[j], times[j][k]);
	}
    }
    
    printf("Function Rate   (MB/s)   RMS time     Min time     Max time\n");
    for (j=0; j<4; j++) {
	rmstime[j] = sqrt(rmstime[j]/(float)NTIMES);

	printf("%s%11.3f  %11.3f  %11.3f  %11.3f\n",
	       label[j],
	       bytes[j]/mintime[j],
	       rmstime[j],
	       mintime[j],
	       maxtime[j]);
    }
    return 0;
}
