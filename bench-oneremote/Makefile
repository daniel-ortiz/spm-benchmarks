#CFLAGS=-g -O2
CC=gcc
#remove moveall to disable move of pages

ifndef MOVEALL
	CFLAGS=-O3 -fopenmp -march=native -std=gnu99 -lnuma 
	OUTNAME=distgen-twoth-nomov
else
	CFLAGS=-O3 -fopenmp -march=native -std=gnu99 -lnuma -DMOVEALL=1
	OUTNAME=distgen-twoth-movall
endif

all: distgen

distgen:
	$(CC) -c distgen.c -o $(OUTNAME)  $(CFLAGS)
clean:
	rm -f *.o distgen-twoth
