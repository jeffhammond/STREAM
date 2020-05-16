CC = icc
CFLAGS = -O3 -xCORE-AVX512 \
	 -qopt-zmm-usage=low \
         -qopenmp \
         -mcmodel=medium \
         -shared-intel \
         -DSTREAM_ARRAY_SIZE=134217728 \
         -DOFFSET=0 \
         -DNTIMES=100

#CC = gcc
#CFLAGS = -O3 -march=skylake-avx512 \
	 -fopenmp \
         -mcmodel=medium \
         -DSTREAM_ARRAY_SIZE=134217728 \
         -DOFFSET=0 \
         -DNTIMES=100

TEMPORAL = -qopt-streaming-stores never
NONTEMPORAL = -qopt-streaming-stores always

FC = ifort
FFLAGS = -O3 -xCORE-AVX512 -qopenmp -qopt-zmm-usage=low

all: stream_f.exe stream_c.exe stream_c.temporal stream_c.nontemporal
#all: stream_c.exe

stream_f.exe: stream.f mysecond.o
	$(CC) $(CFLAGS) -c mysecond.c
	$(FC) $(FFLAGS) -c stream.f
	$(FC) $(FFLAGS) stream.o mysecond.o -o stream_f.exe

stream_c.temporal: stream.c
	$(CC) $(CFLAGS) $(TEMPORAL) stream.c -o stream_c.temporal

stream_c.nontemporal: stream.c
	$(CC) $(CFLAGS) $(NONTEMPORAL) stream.c -o stream_c.nontemporal

stream_c.exe: stream.c
	$(CC) $(CFLAGS) stream.c -o stream_c.exe

clean:
	rm -f stream_f.exe stream_c.exe stream_c.temporal stream_c.nontemporal *.o
