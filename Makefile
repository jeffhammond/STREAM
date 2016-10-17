C99FLAG = -std=c99

CC = icc
CFLAGS = -O3 \
         -qopenmp \
         -mcmodel=medium \
         -shared-intel \
         -DSTREAM_ARRAY_SIZE=134217728 \
         -DOFFSET=0 \
         -DNTIMES=100 \
	 -xMIC-AVX512

TEMPORAL = -qopt-streaming-stores never
NONTEMPORAL = -qopt-streaming-stores always

FC = ifort
FFLAGS = -O3 -qopenmp -xMIC-AVX512

ifdef SORT_TIMES
    CFLAGS += -DSORT_TIMES
endif

all: stream_f.exe stream_c.exe stream_c.temporal stream_c.nontemporal test_qsort

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

test_qsort: test_qsort.c
	$(CC) $(CFLAGS) $(C99FLAG) $< -o $@

clean:
	rm -f *.o
	rm -f stream_f.exe stream_c.exe
	rm -f stream_c.temporal stream_c.nontemporal
	rm -f test_qsort

