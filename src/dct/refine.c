
#include "types.h"
#include "refine.h"


float* find_extremas(value *gvals, ulong size){
  float *extrema = calloc(2, sizeof(float));
  extrema[0] = extrema[1] = gvals[0];
  for (ulong i = 0; i < size; i++) {
    if (gvals[i] < extrema[0])
      extrema[0] = gvals[i];
    if (gvals[i] > extrema[1])
      extrema[1] = gvals[i];
  }
  return extrema;
}
    
void apply_log(value *gvals, ulong size){
  for (ulong i = 0; i <size; i++) 
    gvals[i] = log(gvals[i]);
}

void apply_log10(value *gvals, ulong size){
  for (ulong i = 0; i < size; i++) 
    gvals[i] = log10(gvals[i]);
}
	 
ulong *histogram_tile(value *gvals, int samples, float maxval, float minval, ulong size){
  ulong *hist = calloc(samples, sizeof(ulong));
  float range = (float) maxval - minval;
  float size_sample = range / samples;
  for(ulong i = 0; i < size; i++){
    if(gvals[i] != gvals[i]){
      debug("NaN value spotted");
      gvals[i] = 0.0;
    }
    int subsample = (int)((gvals[i] - minval)/size_sample);
    hist[subsample]++;
  }
  return hist;
}

void apply_refine(value *gvals, ulong size, float *factors){
  #pragma omp parallel for
  for (ulong i = 0; i < size; i++) 
    gvals[i] = (gvals[i]-factors[0])*factors[1];

   float *extrema_allS = find_extremas(gvals, size);
   float *extrema_allR = find_extremas(gvals, size);

   MPI_Allreduce(extrema_allS, extrema_allR, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
   debug("New minimum value: %f", extrema_allR[0]);
   MPI_Allreduce(extrema_allS+1, extrema_allR+1, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   debug("New maximum value: %f", extrema_allR[1]);
}



void check_refine(Arguments *args, value *gvals, ulong size){
  float factors[2];
  float *extrema_allS = find_extremas(gvals, size);
  float *extrema_allR = find_extremas(gvals, size);

  MPI_Allreduce(extrema_allS, extrema_allR, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(extrema_allS+1, extrema_allR+1, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  info("Min and Max values: [%f, %f]", extrema_allR[0],extrema_allR[1]);
  
  if ((args->bpp_arg > 0 && args->bpp_arg <= 16) || args->refine_arg == 0){
    if(extrema_allR[0] < 0){
      info("Negative values spotted, transforming the set in positive values");
      factors[1]        = 1;
      factors[0]        = extrema_allR[0];
      //  apply_refine(gvals, size, factors);
      g_max_greyval = extrema_allR[1];//-extrema_all[0];
      info("New min and Max values in the full volume: [%f, %f]", 0, g_max_greyval);
    } else {
      info("No refinement");
      g_max_greyval = extrema_allR[1];
    }
    g_max_levels  = args->bpp_arg < 0? 2*FLT_MAX : pow(2, args->bpp_arg);
  } else {    
    info("Refining the volume into 16-bits per pixel");
    float range       = (float) extrema_allR[1] - extrema_allR[0];
    float size_refine = (float) pow(2,16) / range;
    factors[1]        = size_refine;
    factors[0]        = extrema_allR[0];
    args->bpp_arg     = 16;
    g_max_greyval     = trunc((extrema_allR[1]-factors[0])*factors[1]);
    g_max_levels      = pow(2,16);
    apply_refine(gvals, size, factors);
  }
}



