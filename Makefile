CC = gcc
FC = gfortran
CFLAGS = -O2 -fopenmp \
		-DSTREAM_ARRAY_SIZE=$(size) \
		-DNTIMES=$(ntime) \
		-DOFFSET=$(offset)
FFLAGS = -O2 -fopenmp

size ?= 10000000
ntime ?= 20
offset ?= 0
TARGET = stream_f stream_c

all: $(TARGET)

stream.o: stream.f 
	$(FC) $(FFLAGS) -c stream.f

mysecond.o: mysecond.c
	$(CC) $(CFLAGS) -c mysecond.c

stream_f: stream.o mysecond.o
	$(FC) $(FFLAGS) stream.o mysecond.o -o $@

stream_c: stream.c
	$(CC) $(CFLAGS) stream.c -o $@

clean:
	rm -f $(TARGET) *.o

# an example of a more complex build line for the Intel icc compiler
stream.icc: stream.c
	icc -O3 -xCORE-AVX2 -ffreestanding -qopenmp -DSTREAM_ARRAY_SIZE=80000000 -DNTIMES=20 stream.c -o stream.omp.AVX2.80M.20x.icc
