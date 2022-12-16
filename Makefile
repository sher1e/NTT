CC = /usr/bin/gcc
CFLAGS = -Wall -Wextra -g -O3 -fomit-frame-pointer -march=native

HEADERS = params.h poly.h ntt.h reduce.h
SOURCES = poly.c reduce.c ntt.c precomp.c

all:	test_ntt speed_ntt


test_ntt: $(HEADERS) $(SOURCES) test_ntt.c
	$(CC) $(CFLAGS) $(SOURCES) test_ntt.c -o $@

speed_ntt: $(HEADERS) $(SOURCES) speed.c  cpucycles.h cpucycles.c
	$(CC) $(CFLAGS) $(SOURCES) cpucycles.c speed.c -o $@

.PHONY: clean

clean:
	-rm test_ntt
	-rm speed_ntt
