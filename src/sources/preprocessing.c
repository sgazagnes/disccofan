
#include <math.h>
#include "types.h"
#include "preprocessing.h"


void find_extremas(value *gvals, ulong size, float *extremas, float *extremas_local){
  extremas[0] = extremas[1] = gvals[0];
  for (ulong i = 0; i < size; i++) {
    if (gvals[i] < extremas[0])
      extremas[0] = gvals[i];
    if (gvals[i] > extremas[1])
      extremas[1] = gvals[i];
  }
  extremas_local[0] = extremas[0];
  extremas_local[1] = extremas[1];
}

void MPI_extremas(float *extremas, float *extremas_local){
  if(np()>1){
    MPI_Allreduce(&extremas_local[0], &extremas[0], 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&extremas_local[1], &extremas[1], 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    debug("Min and Max values in entire data set: [%f, %f]", extremas[0],extremas[1]);
    debug("Min and Max values in this process: [%f, %f]", extremas_local[0],extremas_local[1]);
  }else{
    debug("Min and Max values in data set: [%f, %f]", extremas[0],extremas[1]);
  }
  data_properties.g_max_gval = extremas[1];
  data_properties.g_min_gval = extremas[0];
}


void apply_log(value *gvals, ulong size){
  for (ulong i = 0; i <size; i++) 
    gvals[i] = log(gvals[i]);
}

void apply_log10(value *gvals, ulong size){
  for (ulong i = 0; i < size; i++) 
    gvals[i] = log10(gvals[i]);
}
	 
void apply_sqrt(value *gvals, ulong size){
  for (ulong i = 0; i < size; i++) 
    gvals[i] = sqrt(gvals[i]);
}

void apply_exp(value *gvals, ulong size){
  for (ulong i = 0; i < size; i++) 
    gvals[i] = exp(gvals[i]);
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
}


void preprocessing(Arguments *args, value *gvals, ulong size){

  float extremas[2], extremas_local[2];
  float factors[2];
  find_extremas(gvals, size, extremas, extremas_local);
  MPI_extremas(extremas, extremas_local);

  if (!strcmp(args->tree_type, "min"))
  {
    info("Min tree: inverting the image");
    value max_val = data_properties.g_max_gval;
#pragma omp parallel for
    for (ulong i = 0; i < size; i++)
      gvals[i] = max_val - gvals[i];
  }
  if (!strcmp(args->preprocessing, "ubyte"))
  {
    info("Transforming into 8 bit-per-pixel");
    float range       = (float) extremas[1] - extremas[0];
    float size_refine = (float) (pow(2,8)-1) / range;
    factors[1]        = size_refine;
    factors[0]        = extremas[0];
    args->bpp         = 8;
    data_properties.g_max_gval     = trunc((extremas[1]-factors[0])*factors[1]);
    g_max_levels      = pow(2,8);
    apply_refine(gvals, size, factors);
    find_extremas(gvals, size, extremas, extremas_local);
    MPI_extremas(extremas, extremas_local);
  }
  else if (!strcmp(args->preprocessing, "ushort"))
  {
    info("Transforming into 16 bit-per-pixel");
    float range       = (float) extremas[1] - extremas[0];
    float size_refine = (float) (pow(2,16)-1) / range;
    factors[1]        = size_refine;
    factors[0]        = extremas[0];
    args->bpp         = 16;
    data_properties.g_max_gval     = trunc((extremas[1]-factors[0])*factors[1]);
    g_max_levels      = pow(2,16);
    apply_refine(gvals, size, factors);
    find_extremas(gvals, size, extremas, extremas_local);
    MPI_extremas(extremas, extremas_local);
  }
  else if (!strcmp(args->preprocessing, "uint"))
  {
    info("Transforming into 32 bit-per-pixel");
    float range       = (float) extremas[1] - extremas[0];
    float size_refine = (float) (pow(2,32)-1) / range;
    factors[1]        = size_refine;
    factors[0]        = extremas[0];
    args->bpp         = 32;
    data_properties.g_max_gval     = trunc((extremas[1]-factors[0])*factors[1]);
    g_max_levels      = pow(2,32);
    apply_refine(gvals, size, factors);
    find_extremas(gvals, size, extremas, extremas_local);
    MPI_extremas(extremas, extremas_local);
  }
  else if (!strcmp(args->preprocessing, "log"))
  {
    if(FLOAT_TYPE == 0){
      error("Log preprocessing but floating type not activated, exitinng");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    info("Applying log transformation");
    apply_log(gvals, size);
    args->bpp         = -32;
    find_extremas(gvals, size, extremas, extremas_local);
    MPI_extremas(extremas, extremas_local);
  }
  else if (!strcmp(args->preprocessing, "log10"))
  {
    if(FLOAT_TYPE == 0){
      error("Log10 preprocessing but floating type not activated, exitinng");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    info("Applying log10 transformation");
    apply_log10(gvals, size);   
     args->bpp         = -32;
    find_extremas(gvals, size, extremas, extremas_local);
    MPI_extremas(extremas, extremas_local);
  }
  else if (!strcmp(args->preprocessing, "sqrt"))
  {
    if(FLOAT_TYPE == 0){
      error("sqrt preprocessing but floating type not activated, exitinng");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    info("Applying sqrt transformation");
    apply_sqrt(gvals, size);
    args->bpp         = -32;
    find_extremas(gvals, size, extremas, extremas_local);
    MPI_extremas(extremas, extremas_local);
  }
  else if (!strcmp(args->preprocessing, "exp"))
  {
    if(FLOAT_TYPE == 0){
      error("Exp preprocessing but floating type not activated, exitinng");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    info("Applying exp transformation");
    apply_exp(gvals, size);
    args->bpp         = -32;
    find_extremas(gvals, size, extremas, extremas_local);
    MPI_extremas(extremas, extremas_local);
  }
  else
    info("No pre-processing of the data");
}


