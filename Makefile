CC = gcc
CFLAGS = $(level) -fopenmp \
		-DSTREAM_ARRAY_SIZE=$(size) \
		-DNTIMES=$(ntime) \
		-DOFFSET=$(offset)
TARGET := stream

level ?= -02
size ?= 10000000
ntime ?= 20
offset ?= 0

all: clean $(TARGET)

$(TARGET): stream.c
	$(CC) $(CFLAGS) stream.c -o $@

clean:
	rm -f $(TARGET)
