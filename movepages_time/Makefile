#CFLAGS=-g -O2
CC=gcc-4.8
#remove moveall to disable move of pages
CFLAGS=-O3 -fopenmp -march=native -std=gnu99 -lnuma -g


all: distgen

distgen:
	$(CC) distgen.c -o distgen  $(CFLAGS)
clean:
	rm -f *.o distgen
