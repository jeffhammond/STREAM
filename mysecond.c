/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.

   This version defines two entry points -- with 
   and without appended underscores, so it *should*
   automagically link with FORTRAN */

#ifndef _WIN32
# include <sys/time.h>
#else
# define WIN32_LEAN_AND_MEAN
# include <Windows.h>
#endif

double mysecond()
{
/* struct timeval { long        tv_sec;
            long        tv_usec;        };

struct timezone { int   tz_minuteswest;
             int        tz_dsttime;      };     */

#ifndef _WIN32
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
#else
        __int64 t;
        GetSystemTimeAsFileTime((FILETIME*)&t);
        return ((double)(t - 116444736000000000LL)) / 10000000.0;
#endif
}

double mysecond_() {return mysecond();}

