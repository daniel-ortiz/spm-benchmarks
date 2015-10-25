#CFLAGS=-g -O2
CC=gcc
#remove moveall to disable move of pages

ifndef MOVEALL
	CFLAGS=-O3 -fopenmp -march=native -std=gnu99 -lnuma -g3
	OUTNAME=distgen-allth-nomov
else
	CFLAGS=-O3 -fopenmp -march=native -std=gnu99 -lnuma -g3 -DMOVEALL=1
	OUTNAME=distgen-allth-movall
endif

all: distgen

distgen:
	$(CC)  distgen.c -o $(OUTNAME)  $(CFLAGS)
clean:
	rm -f *.o distgen-allth*

