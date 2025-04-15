#include "types.h"
#include "refine.h"
#include "../misc/kmeans.c"

/*static int bg_find_usable_tile_size(image* img,  int *tile_width,  int *tile_height, float significance_level)
{        
   int current_size = BG_TILE_SIZE_START;
    
  if (bg_available_tiles(img, current_size, current_size,
    significance_level) == TRUE)
  {
    current_size *= 2;
    
    while (current_size <= BG_TILE_SIZE_MAX &&
      bg_available_tiles(img, current_size, current_size,
        significance_level) == TRUE)
    {
      current_size *= 2;
    }
    
    current_size /= 2;
  }
  else
  {
    current_size /= 2;

    if (current_size <  BG_TILE_SIZE_MIN)
    {
      return FALSE;
    }
    
    while (bg_available_tiles(img, current_size, current_size,
      significance_level) == FALSE)
    {
      current_size /= 2;
      
      if (current_size <  BG_TILE_SIZE_MIN)
      {
        return FALSE;
      }
    }    
  }
  
  *tile_width = current_size;
  *tile_height = current_size;
  
  return TRUE;
  }*/


/*void bg_info( value* gval, float *mean, float *variance)
{
  int tile_width = BG_TILE_SIZE_START;
  int tile_height = BG_TILE_SIZE_START;
  
  if (bg_find_usable_tile_size(img, &tile_width, &tile_height,
    BG_REJECTION_RATE) == FALSE)
  {
    bg_error("could not find usable tiles.");
  }
  
  if (verbosity_level)
  {
    printf("Using a tile size of %dx%d in the background estimation.\n",
      tile_height, tile_width);
  }
    
  bg_collect_info(img, tile_width, tile_height, mean,
    variance, BG_REJECTION_RATE, verbosity_level); 
    }*/


void bg_subtract(value *gvals, ulong size, value mean)
{
  #pragma omp parallel
  for (int i = 0; i < size; ++i)
  {
    gvals[i] -= mean;
  }
}

void bg_truncate(value *gvals, ulong size)
{
  int i;
  #pragma omp parallel for private(i)
  for (i = 0; i < size; ++i)
  {
    if (gvals[i] < 0.0)
      gvals[i] = 0;
  }
}


void find_extremas(value *gvals, ulong size, double extrema[2]){
  // extrema[0] = extrema[1] = (double) gvals[0];
  /*for (ulong i = 0; i < size; i++) {
    if (gvals[i] < (value) extrema[0])
      extrema[0] = (double) gvals[i];
    if (gvals[i] > (value) extrema[1])
      extrema[1] = (double) gvals[i];
      }*/
  double min = (double) gvals[0], max = (double) gvals[0];
  #pragma omp parallel for reduction(max:max) reduction(min:min) 
  for (ulong i = 0; i < size; i++) {
    if (gvals[i] < (value) min)
      min = (double) gvals[i];
    if (gvals[i] > (value) max)
      max = (double) gvals[i];
  }
  extrema[0] = min;
  extrema[1] = max;
}
    
void apply_log(value *gvals, ulong size){
  #pragma omp parallel for
  for (ulong i = 0; i <size; i++) 
    gvals[i] = (value) log(gvals[i]);
}

void apply_log10(value *gvals, ulong size){
  #pragma omp parallel for
  for (ulong i = 0; i < size; i++) 
    gvals[i] = (value) log10(gvals[i]);
}
	 
ulong *histogram_tile(value *gvals, int samples, double maxval, double minval, ulong size){
  ulong *hist = calloc(samples, sizeof(ulong));
  double range = (double) maxval - minval;
  double size_sample = range / (double) samples;
  for(ulong i = 0; i < size; i++){
    if(gvals[i] != gvals[i]){
      debug("NaN value spotted");
      gvals[i] = 0.0;
    }
    int subsample = (((int) gvals[i] - (int) minval)/(int) size_sample);
    hist[subsample]++;
  }
  return hist;
}




double *refine(Arguments *args, value *gvals, double extrema[2], ulong size){
  double *factors = calloc(2,sizeof(double));
  factors[0] = 0; factors[1] = 1;
  find_extremas(gvals, size, extrema);
  MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  info("Min value found: %3.3lf", extrema[0]);
  info("Max value found: %3.3lf", extrema[1]);

  g_max_greyval = extrema[1]+extrema[0];
  g_max_levels  = args->bpp_arg < 0 ? size : pow(2, args->bpp_arg);
 
  if(args->refine_arg == NULL)
    info("Running on original dataset");
  else{
    int init_size = strlen(args->refine_arg);
    char delim[] = "_";

    char *ptr = strtok(args->refine_arg, delim);
    while(ptr != NULL) {
	
      if (!strcmp(ptr, "rmean")){
	info("Removing mean double from dataset");
	double tot = 0;
	ulong  totsize = 0;
	for(ulong i = 0; i < size; i++)
	  tot += (double) gvals[i];

	MPI_Allreduce(&tot, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&size, &totsize, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	double mean =  tot/(double)totsize;
	info("Average in data is %lf", mean);
	#pragma omp parallel for
	for(ulong i = 0; i < size; i++)
	  gvals[i] -= (double) mean;

	find_extremas(gvals, size, extrema);
	MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	info("New min value found: %3.3lf", extrema[0]);
	info("New max value found: %3.3lf", extrema[1]);
	g_max_greyval = extrema[1]+extrema[0];
      }
      else if (!strcmp(ptr, "dstd")){
	info("Dividing by std");
	double tot = 0, tot2 = 0;
	ulong totsize = 0;
	for(ulong i = 0; i < size; i++){
	  tot += (double) gvals[i];
	  tot2 += (double) gvals[i]*gvals[i];
	}
    	MPI_Allreduce(&tot, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&tot2, &tot2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(&size, &totsize, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	double std = sqrt((tot2 - (tot*tot)/totsize)/(totsize-1));
	info("STD in data is %lf", std);

	//	#pragma omp parallel for
	for(ulong i = 0; i < size; i++)
	  gvals[i] /= std;
    
	find_extremas(gvals, size, extrema);
	MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	info("New min value found: %3.3lf", extrema[0]);
	info("New max value found: %3.3lf", extrema[1]);
	g_max_greyval = extrema[1]+extrema[0];
      }
      else if (!strcmp(ptr, "thresh")){
	info("Thresholding with limits: [%lf, %lf]", args->lims_arg[0], args->lims_arg[1]);
	#pragma omp parallel for
	for(ulong i = 0; i < size; i++)
	  gvals[i] = gvals[i] >= args->lims_arg[0] && gvals[i] <= args->lims_arg[1] ? 255:  0;

    
	extrema[1] = 255;
	extrema[0] = 0;
	g_max_greyval = extrema[1]+extrema[0];
	g_max_levels = 256;
	args->bpp_arg = 8; 

      }

      else if( !strcmp(ptr, "kmeans")){
	double center[2] = {0, 10};
	if( args->lims_arg != NULL){
	  center[0] = args->lims_arg[0];
	  center[1] = args->lims_arg[1];
	}
	info("Thresholding with Kmeans and 2 clusters, initial center [%lf, %lf]", center[0], center[1]);
	int *labels       = malloc(size * sizeof(int));
	double *gvals_cop = malloc(size * sizeof(double));
	#pragma omp parallel for
	for(ulong i =0; i< size; i++)
	  gvals_cop[i]= (double) gvals[i];
    
	kmeans(1, gvals_cop, size, 2, center, labels);
    
	#pragma omp parallel for
	for(ulong i = 0; i < size; i++)
	  gvals[i] = 255 - 255*labels[i];
	free(gvals_cop);

	extrema[1] = 255;
	extrema[0] = 0;
	g_max_greyval = extrema[1]+extrema[0];
	g_max_levels = 256;
	args->bpp_arg = 8; 

      }

 

      else if( !strcmp(ptr, "16bits")){    
	info("Refining the volume into 16-bits per pixel");
	double range       =  extrema[1] - extrema[0];
	double size_refine = (double) pow(2,16) / range;
	factors[1]        = size_refine;
	factors[0]        = extrema[0];
	//args->bpp_arg     = 16;
	g_max_greyval     = trunc((extrema[1]-factors[0])*factors[1]);
	g_max_levels      = pow(2,16);
	#pragma omp parallel for
	for (ulong i = 0; i < size; i++) 
	  gvals[i] = floor((gvals[i]- (value) factors[0])* (value)factors[1]);
	find_extremas(gvals, size, extrema);
	MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	info("New min value found: %3.3lf", extrema[0]);
	info("New max value found: %3.3lf", extrema[1]);
	g_max_greyval = extrema[1]+extrema[0];
      }
      else if( !strcmp(ptr, "12bits")){    
	info("Refining the volume into 16-bits per pixel");
	double range       =  extrema[1] - extrema[0];
	double size_refine = (double) pow(2,12) / range;
	factors[1]        = size_refine;
	factors[0]        = extrema[0];
	//args->bpp_arg     = 12;
	args->flood_arg = 0;
	g_max_greyval     = trunc((extrema[1]-factors[0])*factors[1]);
	g_max_levels      = pow(2,12);
	#pragma omp parallel for
	for (ulong i = 0; i < size; i++) 
	  gvals[i] = floor((gvals[i]- (value) factors[0])* (value)factors[1]);
	find_extremas(gvals, size, extrema);
	MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	info("New min value found: %3.3lf", extrema[0]);
	info("New max value found: %3.3lf", extrema[1]);
	g_max_greyval = extrema[1]+extrema[0];
      }
      else if( !strcmp(ptr, "8bits")){    
	info("Refining the volume into 16-bits per pixel");
	double range       =  extrema[1] - extrema[0];
	double size_refine = (double) pow(2,8) / range;
	factors[1]        = size_refine;
	factors[0]        = extrema[0];
	//args->bpp_arg     = 8;
	args->flood_arg = 0;
	g_max_greyval     = trunc((extrema[1]-factors[0])*factors[1]);
	g_max_levels      = pow(2,8);
	#pragma omp parallel for
	for (ulong i = 0; i < size; i++) 
	  gvals[i] = floor((gvals[i]- (value) factors[0])* (value)factors[1]);
	find_extremas(gvals, size, extrema);
	MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	info("New min value found: %3.3lf", extrema[0]);
	info("New max value found: %3.3lf", extrema[1]);
	g_max_greyval = extrema[1]+extrema[0];
      }
      else if( !strcmp(ptr, "nlev")){    
	info("Refining the volume into given levels per pixel");
	double range       =  extrema[1] - extrema[0];
	double size_refine = (double) (args->nlev_arg-1) / range;
	factors[1]        = size_refine;
	factors[0]        = extrema[0];
	//args->bpp_arg     = 8;
	args->flood_arg = 0;
	g_max_greyval     = trunc((extrema[1]-factors[0])*factors[1]);
	g_max_levels      = args->nlev_arg;
	#pragma omp parallel for
	for (ulong i = 0; i < size; i++) 
	  gvals[i] = (value) trunc((gvals[i]- (value) factors[0])* (value)factors[1]);
	find_extremas(gvals, size, extrema);
	MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	info("New min value found: %3.3lf", extrema[0]);
	info("New max value found: %3.3lf", extrema[1]);
	g_max_greyval = extrema[1]+extrema[0];
      }
      else if( !strncmp(ptr, "prec",4)){    
	info("Cutting precision of float or double");
	int prec           = atoi(ptr+4);
	info("NEW prec %d", prec);
	for (ulong i = 0; i < size; i++) 
	  gvals[i] = floor(pow(10,prec)*gvals[i])/pow(10,prec); 

	find_extremas(gvals, size, extrema);
	MPI_Allreduce(extrema, extrema, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(extrema+1, extrema+1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	info("New min value found: %3.3lf", extrema[0]);
	info("New max value found: %3.3lf", extrema[1]);
	g_max_greyval = extrema[1]+extrema[0];
      }
      
      ptr = strtok(NULL, delim);
    }
  }

    
  return factors;
}



