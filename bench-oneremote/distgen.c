#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>
#include <sys/time.h>
#include <pthread.h>
#include <numaif.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXDISTCOUNT 10
#define BLOCKLEN 64
#define MAXTHREADS 64
//how much time to wait before moving pages
#define PHAS1_TIME 45

struct entry {
  double v;
  struct entry *next;
};
typedef struct {
	double* start[2];
	long size;
} memory_region;


int distSize[MAXDISTCOUNT];
int distIter[MAXDISTCOUNT];
int distBlocks[MAXDISTCOUNT];
int distsUsed = 0;
int verbose = 0;
int trasv_count=10;
long count_through=1;
double wtime()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec+1e-6*tv.tv_usec;
}

void * call_move_pages(void * arg){
	memory_region *mr=(memory_region*)arg;
	void** pages;
	int pages_size=(8*mr->size)/4096, begin=0,j=0;
	
	pages=malloc(pages_size*sizeof(void *));
	int * status=malloc(sizeof(int) * pages_size);
	int * nodes=malloc(sizeof(int) * pages_size);
	
	sleep(PHAS1_TIME);
	for(int jj=0; jj<2; jj++){
		for(int i=0; i<pages_size; i++){
				double * offset= mr->start[jj]+ i*512;
				pages[i]=(void *) offset;
				nodes[i]=1;
			}
			
		for (j=0; j<9; j++ ){
			begin = j* (pages_size/10);
			
			int ret= 	move_pages(getpid(), pages_size/10, pages+begin, nodes+begin, status+begin,0);
		
			 printf("move pages res %d \n", ret);
		}
	j++;
	begin = 9* (pages_size/10);
	int ret= 	move_pages(getpid(), pages_size-(9*pages_size/10), pages+begin, nodes+begin, status+begin,0);
	}
		
}

void addDist(int size)
{
  int d, dd;

  // distances are sorted
  for(d=0; d<distsUsed; d++) {
    if (distSize[d] == size) return; // one dist only once
    if (distSize[d] < size) break;
  }

  assert(distsUsed < MAXDISTCOUNT);
  for(dd = distsUsed; dd >= d; dd--)
    distSize[dd] = distSize[dd-1];

  distSize[d] = size;
  distsUsed++;
}

void initBufs(int blocks)
{
  int d, i;

  for(d=0; d<distsUsed; d++) {
    // each memory block of cacheline size gets accessed
    distBlocks[d] = (distSize[d] + BLOCKLEN - 1) / BLOCKLEN;
    distIter[d] = distSize[0] / distSize[d];
  }

  if (verbose) {
    fprintf(stderr, "Number of distances: %d\n", distsUsed);
    for(d=0; d<distsUsed; d++)
      fprintf(stderr, "  D%2d: size %d (iter %d)\n",
      	d+1, distSize[d], distIter[d]);
  }
}

void usage(char* argv0)
{
  fprintf(stderr,
	  "Distance Generator\n"
	  "Usage: %s [Options] [-<iter>] [<dist1> [<dist2> ... ]]\n"
	  "\nParameters:\n"
	  "  <iter>       number of times accessing arrays (def. 1000)\n"
	  "  <dist1>, ... different reuse distances (def. 1 dist with 1MB)\n"
	  "\nOptions:\n"
	  "  -h           show this help\n"
	  "  -p           use pseudo-random access pattern\n"
	  "  -d           travers by dependency chain\n"
	  "  -c <MHz>     provide clock frequency to show cycles per access\n"
	  "  -v           be verbose\n", argv0);
  fprintf(stderr,
	  "\nNumbers can end in k/m/g for Kilo/Mega/Giga factor\n");
  exit(1);
}

// helper for adjustSize
int gcd(int a, int b)
{
  if (b == 0) return a;
  return gcd(b, a % b);
}

// make sure that gcd(size,diff) is 1 by increasing size, return size
int adjustSize(int size, int diff)
{
  while(gcd(size,diff) >1) size++;
  return size;
}

void runBench(struct entry* buffer,
	      int iter, int blocks, int blockDiff, int depChain,
	      double* sum, unsigned long* aCount)
{
  int i, d, k, j, idx, max;
  double lsum = *sum;
  int idxIncr = blockDiff * BLOCKLEN/sizeof(struct entry);
  int idxMax = blocks * BLOCKLEN/sizeof(struct entry);
	double time1,diff;
	time1=wtime();
  for(i=0; i<iter; i++) {
    lsum += buffer[0].v;
    for(d=0; d<distsUsed; d++) {
      lsum += buffer[0].v;
      for(k=0; k<distIter[d]; k++) {
	//fprintf(stderr, "D %d, I %d\n", d, k);
		*aCount += distBlocks[d];
		max = distBlocks[d];
		count_through++;
		
			if (!depChain) {
			  idx = 0;
			  for(j=0; j<max; j++) {
				lsum += buffer[idx].v;
				idx += idxIncr;
				if (idx >= idxMax) idx -= idxMax;
				//fprintf(stderr, " Off %d\n", idx);
			  }
			}
			else {
			  struct entry* p = buffer;
			  for(j=0; j<max; j++) {
				lsum += p->v;
				p = p->next;
				//fprintf(stderr, " POff %d\n", (int)(p - buffer));
			  }
			}
			if(count_through==trasv_count){
				diff=wtime()-time1;
				diff=(trasv_count*max)/diff;
				time1=wtime();
				count_through=1;
				printf("%f \n",diff);
			}
      }
    }
  }
  printf("counts %ld -", count_through);

  *sum = lsum;
}

char* prettyVal(char *s, unsigned long v)
{
  static char str[20];

  if (!s) s = str;
  if (v > 1000000000)
    sprintf(s,  "%.1f G", 1.0 / 1024.0 / 1024.0 / 1024.0 * v);
  else if (v > 1000000)
    sprintf(s,  "%.1f M", 1.0 / 1024.0 / 1024.0 * v);
  else if (v > 1000)
    sprintf(s,  "%.1f K", 1.0 / 1024.0 * v);
  else
    sprintf(s,  "%lu", v);

  return s;
}

int toInt(char* s, int isSize)
{
  char* end;
  int d;
  int f = isSize ? 1024 : 1000;

  d = strtol(s, &end, 10);
  if ((*end == 'k') || (*end == 'K')) d = d * f;
  else if ((*end == 'm') || (*end == 'M')) d = d * f * f;
  else if ((*end == 'g') || (*end == 'G')) d = d * f * f * f;
  return d;
}

int main(int argc, char* argv[])
{
  int arg, d, i, j, k, idx;
  int iter = 0;
  int pseudoRandom = 0;
  int depChain = 0;
  int clockFreq = 2400;
  unsigned long aCount = 0;
  int blocks, blockDiff;
  double sum = 0.0;
  int t, tcount;
  struct entry* buffer[MAXTHREADS];
  double tt;
  int thr=0;
  memory_region mr;
  pthread_t move_thread;
  
  verbose = 1;  
  for(arg=1; arg<argc; arg++) {
    if (argv[arg][0] == '-') {
      if (argv[arg][1] == 'h') usage(argv[0]);
      if (argv[arg][1] == 'v') { verbose++; continue; }
      if (argv[arg][1] == 'p') { pseudoRandom = 1; continue; }
      if (argv[arg][1] == 'd') { depChain = 1; continue; }
      if (argv[arg][1] == 'c') {
	if (arg+1<argc) {
	  clockFreq = atoi(argv[arg+1]);
	  arg++;
	}
	continue;
      }
      iter = toInt(argv[arg]+1, 0);
      if (iter == 0) usage(argv[0]);
      continue;
    }
    d = toInt(argv[arg], 1);
    if (d <= 0) usage(argv[0]);
    addDist(d);
  }

  if (distsUsed == 0) addDist(1024*1024);
  if (iter == 0) iter = 1000;

  blocks = (distSize[0] + BLOCKLEN - 1) / BLOCKLEN;  
  blockDiff = pseudoRandom ? (blocks * 7/17) : 1;
  blocks = adjustSize(blocks, blockDiff);
  initBufs(blocks);

#pragma omp parallel
#ifdef _OPENMP
  tcount = omp_get_num_threads();
#else
  tcount = 1;
#endif

  // calculate expected number of accesses
  aCount = 0;
  for(d=0; d<distsUsed; d++)
    aCount += distIter[d] * distBlocks[d];

  if (verbose) {
    char sBuf[20], tsBuf[20], acBuf[20], tacBuf[20];
    prettyVal(sBuf, BLOCKLEN * blocks);
    prettyVal(tsBuf, BLOCKLEN * blocks * tcount);
    prettyVal(acBuf, aCount);
    prettyVal(tacBuf, aCount * tcount * iter);

    fprintf(stderr, "Buffer size per thread %sB (total %sB), address diff %d\n",
	    sBuf, tsBuf, BLOCKLEN * blockDiff);
    fprintf(stderr, "Accesses per iteration and thread: %s, total %s\n",
	    acBuf, tacBuf);
    fprintf(stderr, "Iterations: %d, threads: %d\n",
	    iter, tcount);
  }

  assert(tcount < MAXTHREADS);
  assert(sizeof(struct entry) == 16);
  
 //Modified so that initialization is done only in the master
	#pragma omp master
	{
    struct entry *next, *buf;
    int idx, blk, nextIdx;
    int idxMax = blocks * BLOCKLEN/sizeof(struct entry);
    int idxIncr = blockDiff * BLOCKLEN/sizeof(struct entry);

    // allocate and initialize used memory
    buffer[t] = (struct entry*) memalign(64, blocks * BLOCKLEN);
    mr.size=blocks * BLOCKLEN;
    mr.start[t]= buffer[t];
    
    buf = buffer[t];
    for(idx=0; idx < idxMax; idx++) {
      buf[idx].v = (double) idx;
      buf[idx].next = 0;
    }
    // generate dependency chain
    idx = 0;
    for(blk=0; blk < blocks; blk++) {
      nextIdx = idx + idxIncr;
      if (nextIdx >= idxMax) nextIdx -= idxMax;
      //fprintf(stderr, " Blk %d, POff %d\n", blk, nextIdx);
      assert(buf[idx].next == 0);
      buf[idx].next = buf + nextIdx;
      idx = nextIdx;
    }
  } 
   #ifdef MOVEALL	
    thr = pthread_create(&move_thread, NULL,&call_move_pages, &mr);
   #endif
  if (verbose)
    fprintf(stderr, "Running ...\n");

  aCount = 0;
  tt = wtime();
//Modified so that access is done only by other threads than the master
//it is advisable to use pinning to place this threads on the other numa node
#pragma omp parallel
{ 
    if (omp_get_thread_num()!=0){
		double tsum = 0.0;
		unsigned long taCount = 0;

		runBench(buffer[t], iter, blocks, blockDiff, depChain,
			 &tsum, &taCount);

		sum += tsum;
		aCount += taCount;
	}
  }

  tt = wtime() - tt;

  if (verbose) {
    double avg = tt * tcount / aCount * 1000000000.0;
    double cTime = 1000.0 / clockFreq;

    fprintf(stderr, "Finished (ACount: %lu, sum: %g)\n", aCount, sum);
    fprintf(stderr, "Elapsed: %.3fs => %.3f GB/s, %.3f GF/s"
	    " (per core: %.3f GB/s, %.3f GF/s)\n",
	    tt,
	    aCount * 64.0 / tt / 1000000000.0,
	    aCount / tt / 1000000000.0,
	    aCount * 64.0 / (tt * tcount) / 1000000000.0,
	    aCount / (tt * tcount) / 1000000000.0 );
    fprintf(stderr, " avg. time per access: %.3f ns (%.1f cycles @ %.1f GHz)\n",
	    avg, avg/cTime, 1.0 / 1000.0 * clockFreq);
  }

  return 0;
}

