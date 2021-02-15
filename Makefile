CC = clang
CFLAGS = -O3 \
         -DSTREAM_ARRAY_SIZE=134217728 \
         -DOFFSET=0 \
         -DNTIMES=100

FC = gfortran
FFLAGS = -O3

all: stream_f.exe stream_c.exe
#all: stream_c.exe

stream_f.exe: stream.f mysecond.o
	$(CC) $(CFLAGS) -c mysecond.c
	$(FC) $(FFLAGS) -c stream.f
	$(FC) $(FFLAGS) stream.o mysecond.o -o stream_f.exe

stream_c.exe: stream.c
	$(CC) $(CFLAGS) stream.c -o stream_c.exe

clean:
	rm -f stream_f.exe stream_c.exe *.o
