===============================================

STREAM is the de facto industry standard benchmark
for measuring sustained memory bandwidth.

Documentation for STREAM is on the web at:
   http://www.cs.virginia.edu/stream/ref.html

===============================================
NEWS
===============================================
UPDATE: October 28 2014:

"stream_mpi.c" released in the Versions directory.

Based on Version 5.10 of stream.c, stream_mpi.c
brings the following new features:
* MPI implementation that *distributes* the arrays
  across all MPI ranks. (The older Fortran version
  of STREAM in MPI *replicates* the arrays across
  all MPI ranks.)
* Data is allocated using "posix_memalign" 
  rather than using static arrays.  Different
  compiler flags may be needed for both portability
  and optimization.
  See the READ.ME file in the Versions directory
  for more details.
* Error checking and timing done by all ranks and
  gathered by rank 0 for processing and output.
* Timing code uses barriers to ensure correct
  operation even when multiple MPI ranks run on
  shared memory systems.

NOTE: MPI is not a preferred implementation for
  STREAM, which is intended to measure memory
  bandwidth in shared-memory systems.  In stream_mpi,
  the MPI calls are only used to properly synchronize
  the timers (using MPI_Barrier) and to gather
  timing and error data, so the performance should 
  scale linearly with the size of the cluster.
  But it may be useful, and was an interesting 
  exercise to develop and debug.

===============================================
UPDATE: January 17 2013:

Version 5.10 of stream.c is finally available!

There are no changes to what is being measured, but
a number of long-awaited improvements have been made:

* Updated validation code does not suffer from 
  accumulated roundoff error for large arrays.
* Defining the preprocessor variable "VERBOSE"
  when compiling will (1) cause the code to print the
  measured average relative absolute error (rather than
  simply printing "Solution Validates", and (2) print
  the first 10 array entries with relative error exceeding
  the error tolerance.
* Array index variables have been upgraded from
  "int" to "ssize_t" to allow arrays with more
  than 2 billion elements on 64-bit systems.
* Substantial improvements to the comments in 
  the source on how to configure/compile/run the
  benchmark.
* The proprocessor variable controlling the array
  size has been changed from "N" to "STREAM_ARRAY_SIZE".
* A new preprocessor variable "STREAM_TYPE" can be
  used to override the data type from the default
  "double" to "float".
  This mechanism could also be used to change to 
  non-floating-point types, but several "printf"
  statements would need to have their formats changed
  to accomodate the modified data type.
* Some small changes in output, including printing
  array sizes is GiB as well as MiB.
* Change to the default output format to print fewer
  decimals for the bandwidth and more decimals for
  the min/max/avg execution times.


===============================================
UPDATE: February 19 2009:

The most recent "official" versions have been renamed
"stream.f" and "stream.c" -- all other versions have
been moved to the "Versions" subdirectory and should be
considered obsolete.

The "official" timer (was "second_wall.c") has been
renamed "mysecond.c".   This is embedded in the C version
("stream.c"), but still needs to be externally linked to
the FORTRAN version ("stream.f").  The new version defines
entry points both with and without trailing underscores,
so it *should* link automagically with any Fortran compiler.

===============================================

STREAM is a project of "Dr. Bandwidth":
	John D. McCalpin, Ph.D.
	john@mccalpin.com

===============================================

The STREAM web and ftp sites are currently hosted at
the Department of Computer Science at the University of
Virginia under the generous sponsorship of Professor Bill
Wulf and Professor Alan Batson.

===============================================
