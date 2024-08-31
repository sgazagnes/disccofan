
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <sys/types.h>
#include <sys/times.h>
#include <stdarg.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <libgen.h>
#include <ctype.h>
#include <assert.h>


typedef uint8_t  ubyte;
typedef uint16_t ushort;
typedef uint32_t uint;
typedef uint64_t ulong;
typedef int64_t  idx;
typedef uint     value_t;



/* If not floating point : */

//#define FLOAT_TYPE 0
//typedef ubyte value; 			 // 8-bit (not working with flooding 1)
//typedef ushort value;			 // 16-bit
//typedef uint value; 			 // 32-bit
//typedef ulong value; 			 // 64-bit

/* If floating point:     */

#define FLOAT_TYPE 1
typedef float value;

#define MXT_HISTO_SZ_LOG2 8
#define MXT_HISTO_SZ (1U << MXT_HISTO_SZ_LOG2)
#define MXT_HISTO_MASK (MXT_HISTO_SZ - 1)
#define NUM_DIGITS (sizeof(value) * CHAR_BIT + MXT_HISTO_SZ_LOG2 - 1) / MXT_HISTO_SZ_LOG2
#define NO_INLINE __attribute__((noinline))

/* Teeninga sort item Structure */
#pragma pack(push, 1)
typedef struct _SortItem {
  value_t val;				/* Unsigned value */
  ulong rank;				/* Rank of value */
} SortItem;
#pragma pack(pop)


void read_fits(value **img,  const char *fname, ulong dims[3], int *bitpix){
  /* +++++++++++++++++++++++++++ */
  /*          Read Fits          */
  /* +++++++++++++++++++++++++++ */

  fitsfile 	*fptr;   			/* FITS file pointer */
  int		n_dims;
  int 		type;
  int 		status        = 0;   		/* CFITSIO status value MUST be initialized to zero! */
  long 		inc[4]        = {1,1,1,1};	/* Increment for fits read function */
  long 		naxes[4]      = {1,1,1,1};	/* Dimensions (4D max) */
  long 		counts[4]     = {1,1,1,1};	/* Pixels to read in fits function */
  long 		offsets[4]    = {1,1,1,1};	/* Pixel index to start in fist function */
  
  printf("Reading FITS Image %s \n", fname);

  if(FLOAT_TYPE)
    type = TFLOAT;
  else
    type = (int) (log2(sizeof(value))+1)*10 + ttype;

  if (!fits_open_file(&fptr, fname, READONLY, &status)) {	/* Open file */
    if (!fits_get_img_param(fptr, 4,  bitpix,  &n_dims, naxes, &status)) {  /* Get image parameter */
      if (n_dims > 4 || n_dims == 0 || (n_dims == 4 && naxes[3] > 1)) { /* 3D Max */	   	      
	printf("Only 2D and 3D images are supported\n");
	exit(0);
      }
      else {	
	if(n_dims == 4) (n_dims)--;
	  
	for (int i = n_dims; i-- ; ) 
	  counts[i]  = dims[i] = naxes[i];

	
	printf("Fits reading: \n Tile dimensions: %ld by %ld by %ld \n Offsets %ld, %ld, %ld \n Counts %ld, %ld, %ld \n", dims[0], dims[1], dims[2],offsets[0],offsets[1],offsets[2],counts[0],counts[1],counts[2]);
	*img = malloc(dims[0]*dims[1]*dims[2]* sizeof(value));
	fits_read_subset(fptr, type, offsets, counts, inc, NULL, *img, NULL, &status);
      }
    }
    
    fits_close_file(fptr, &status);        
  } else {
    printf("Cannot open the file %s, wrong name ? (FITS) \n", fname);
    exit(0);
  }    
}


void gen_histogram(value *gvals, ulong histos[][MXT_HISTO_SZ],  ulong lwb, ulong upb, value_t (*functor)(const value)) {
  memset(histos, 0, NUM_DIGITS * MXT_HISTO_SZ * sizeof(histos[0][0]));
  
  for (ulong i = lwb; i != upb; ++i) {
    int shift = 0;
    value_t gval = functor(gvals[i]);
    for (int digit_nr = 0; digit_nr != NUM_DIGITS; ++digit_nr) {
      int digit = (gval >> shift) & MXT_HISTO_MASK;
      ++histos[digit_nr][digit];
      shift += MXT_HISTO_SZ_LOG2;
    }
  }
}

void exclusive_sum(ulong *it, ulong *it_end){
  ulong sum = 0;
  while (it != it_end) {
    ulong next = *it;
    *it++ = sum;
    sum += next;
  }
}

void scatter_first_digit(value *gvals, SortItem *pair_it, ulong* histo, ulong lwb, ulong upb, value_t (*functor)(const value)  ){
  for (ulong i = lwb; i != upb; ++i) {
    value_t gval; 
    gval = functor(gvals[i]);
    int digit = gval & MXT_HISTO_MASK;
    pair_it[histo[digit]++] = (SortItem) {gval, i-lwb};
  }
}

void scatter_digit_bit1(value *gvals, ulong *out, ulong* histo, ulong lwb, ulong upb, value_t (*functor)(const value))
{   
  for (ulong i = lwb; i != upb; ++i)
  {
    int digit = functor(gvals[i]) & MXT_HISTO_MASK;
    out[histo[digit]++] = i-lwb;
  }
}

void scatter_digit(SortItem *in, SortItem *out, int digit_nr, ulong* histo, ulong size) {  
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i) {
    const SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair;
  }
}

void scatter_last_digit(SortItem *in, ulong *out, int digit_nr, ulong* histo, ulong size){
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i)  {
    const SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair.rank;
  }
}


void create_ranks_inv(value *gvals, ulong **ranks_inv, ulong histos[][MXT_HISTO_SZ], ulong lwb, ulong upb, value_t (*functor)(const value)) {

  const int num_digits = NUM_DIGITS;
  //  size_t alloc_pos =allocator_pos(work_alloc);
  // *ranks_inv = (ulong *)  allocator_allocate(work_alloc, upb-lwb, sizeof(ulong));
  //  size_t alloc_pos_rank_inv = allocator_pos(work_alloc);

  if (num_digits == 1) {
    *ranks_inv = calloc(upb-lwb, sizeof(ulong));
    scatter_digit_bit1(gvals, *ranks_inv, histos[0], lwb, upb, functor);
    return;
  }

  //allocator_set_pos(work_alloc, alloc_pos);
  SortItem *pairs1 =  malloc((upb-lwb)*sizeof(SortItem));//*ranks_inv;
  //allocator_allocate(work_alloc, upb-lwb, sizeof(SortItem));
  SortItem *pairs2 =  malloc((upb-lwb)*sizeof(SortItem));//*ranks_inv + (upb-lwb)*sizeof(SortItem);
  //allocator_allocate(work_alloc, upb-lwb, sizeof(SortItem));

  SortItem *pairswap = NULL;
  scatter_first_digit(gvals, pairs2, histos[0], lwb, upb, functor);

  int digit_nr = 1;
  for (; digit_nr != NUM_DIGITS - 1; ++digit_nr)
  {
    pairswap = pairs1;
    pairs1 = pairs2;
    pairs2 = pairswap;
    
    scatter_digit(pairs1, pairs2, digit_nr, histos[digit_nr], (upb-lwb));
  }
  
  free(pairs1);
  *ranks_inv = calloc(upb-lwb, sizeof(ulong));
  scatter_last_digit(pairs2, (ulong *) *ranks_inv, digit_nr, histos[digit_nr], (upb-lwb));
  free(pairs2);
  // allocator_set_pos(work_alloc, alloc_pos_rank_inv);
}

NO_INLINE void create_mappings(value *gvals, ulong **ranks_inv, ulong lwb, ulong upb) {
  const int   num_digits = NUM_DIGITS;
  value_t (*foo)(const value);

  if(FLOAT_TYPE == 1)
    foo = &transform_float;
  else
    foo = &transform_dum;

  
  ulong histos[num_digits][MXT_HISTO_SZ];
  gen_histogram(gvals, histos, lwb, upb, foo);

  for (int digit_nr = 0; digit_nr != num_digits; digit_nr++) 
    exclusive_sum(histos[digit_nr], histos[digit_nr] + MXT_HISTO_SZ);
  
  create_ranks_inv(gvals, ranks_inv, histos, lwb, upb, foo);

  
}


int main(int argc, char** argv) {


  value *gval;
  ulong *ranks 	= NULL;
  ulong 	dims[3]    = {1, 1, 1};	
  int bitpix;
  read_fits(&img,  "test.fits", dims, &bitpix);
  nthreads = 1;
  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			/* Thread number */
    ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads) 
      : dims[0]*dims[1]*((id*dims[2])/nthreads) ;              /* Lower bound for current thread */
    ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/nthreads)
      : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);       /* Upper bound for current thread */

    // Create mappings
    
    ulong *ranks_inv  	= NULL;
    create_mappings(local_tree->gval, &ranks_inv, lwb, upb);

    #pragma omp single
    {
      ranks 	       		= malloc(dims[0]*dims[1]*dims[2] * sizeof(ulong));
    }
    
    for (ulong i = 0; i != upb-lwb; ++i){
      ranks[ranks_inv[i]+lwb] = i;
    }

  }

  
  free(gval);
  free(ranks);
}
