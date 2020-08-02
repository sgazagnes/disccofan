#include "misc.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>
#include <gsl/gsl_randist.h>

#if RAND_MAX/256 >= 0xFFFFFFFFFFFFFF
  #define LOOP_COUNT 1
#elif RAND_MAX/256 >= 0xFFFFFF
  #define LOOP_COUNT 2
#elif RAND_MAX/256 >= 0x3FFFF
  #define LOOP_COUNT 3
#elif RAND_MAX/256 >= 0x1FF
  #define LOOP_COUNT 4
#else
  #define LOOP_COUNT 5
#endif



size_t mod_fwrite (const void *array, unsigned long long size, unsigned long long count, FILE *stream){
  unsigned long long pos, tot_size, pos_ct;
  const unsigned long long block_size = 4*512*512*512;

  tot_size = size*count; // total size of buffer to be written

  //check if the buffer is smaller than our pre-defined block size
  if (tot_size <= block_size)
    return fwrite(array, size, count, stream);
  else{
    pos = 0;
    pos_ct=0;
    // buffer is large, so write out in increments of block_size
    while ( (pos+block_size) < tot_size ){
      fwrite((char *)array+pos, block_size, 1, stream);
      pos_ct++;
      pos = block_size*pos_ct;
    }
    return fwrite((char *)array+pos, tot_size-pos, 1, stream);
  }
}


size_t mod_fread(void * array, unsigned long long size, unsigned long long count, FILE * stream){
  unsigned long long pos, tot_size, pos_ct;
  const unsigned long long block_size = 512*512*512*4;

  tot_size = size*count; // total size of buffer to be written
  //check if the buffer is smaller than our pre-defined block size
  if (tot_size <= block_size)
    return fread(array, size, count, stream);
  else{
    pos = 0;
    pos_ct=0;
    // buffer is large, so write out in increments of block_size
    while ( (pos+block_size) < tot_size ){
      fread((char *)array + pos, block_size, 1, stream);
      pos_ct++;
      pos = block_size*pos_ct;
    }
    return fread((char *)array+pos, tot_size-pos, 1, stream);
  }
}

int compare_floats (const void *a, const void *b)
{
  const float *da = (const float *) a;
  const float *db = (const float *) b;

  return (*da > *db) - (*da < *db);
}

double keithRandom() {
    // Random number function based on the GNU Scientific Library
    // Returns a random float between 0 and 1, exclusive; e.g., (0,1)
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_default; // Generator setup
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, mySeed);
    double u = gsl_rng_uniform(r); // Generate it!
    gsl_rng_free (r);
    return u;
}

int find_scale_double(double *thresh, double val, int nthresh){
  int upper = nthresh-1, lower = 0, mid;

  if (val >=   thresh[upper])
    return upper;


  mid = (upper + lower) / 2;
  while (mid!=lower) {
    if(val >=  thresh[mid])
      lower = mid;
    else
      upper = mid;

    mid = (upper + lower) / 2;
  }
  return lower;
} /* find_scale */


double find_max(double *gvals, ulong size){
  double max = gvals[0];
  for(size_t i = 0; i<size;i++){
    if(gvals[i] > max)
      max = gvals[i];
  }
  return max;
}

double find_min(double *gvals, ulong size){
  double max = gvals[0];
  for(size_t i = 0; i<size;i++){
    if(gvals[i] < max)
      max = gvals[i];
  }
  return max;
}
   
double* find_extremas(double *gvals, ulong size){
  double *extrema = calloc(2, sizeof(double));
  extrema[0] = extrema[1] = gvals[0];
  for (ulong i = 0; i < size; i++) {
    if (gvals[i] < extrema[0])
      extrema[0] = gvals[i];
    if (gvals[i] > extrema[1])
      extrema[1] = gvals[i];
  }
  return extrema;
}
    
void apply_log(double *gvals, ulong size){
  for (ulong i = 0; i <size; i++) 
    gvals[i] = log(gvals[i]);
}

void apply_log10(double *gvals, ulong size){
  for (ulong i = 0; i < size; i++) 
    gvals[i] = log10(gvals[i]);
}

float as_float(uint32_t i){
  union {
    uint32_t i;
    float f;
  } pun = { i };
  return pun.f;
}

double as_double(uint64_t i){
  union  {
    uint64_t i;
    double f;
  } pun = { i };
  return pun.f;
}

double randn (double mu, double sigma){
  
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

double *create_random(ulong size, int bitpix, int neg){
  double *gvals_out = calloc(size, sizeof(double));
  double max_gval_out, min_gval_out;
  srand(time(NULL));
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  struct timeval tv; // Seed generation based on time
  gettimeofday(&tv,0);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;
  T = gsl_rng_default; // Generator setup
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, mySeed);

  
  max_gval_out = bitpix < 0 ? 1 : pow(2, bitpix);

  if(neg){
    max_gval_out /= 2;
    min_gval_out = - max_gval_out;
  } else
    min_gval_out = 0;
  
  printf("Creating a random image/volume of size %lu and %d bit-per-pixel depth. \n Min and max doubles:[%lf, %lf) \n", size, bitpix, min_gval_out,  max_gval_out);

  for(ulong i = 0; i < size; i++){
    double tt = gsl_ran_flat( r, min_gval_out, max_gval_out);
    gvals_out[i] =  bitpix > 0? trunc(tt):tt;//gsl_ran_flat( r, min_gval_out, max_gval_out);
  }
  return gvals_out;
  }

/*double *create_random(ulong size, int bitpix, int neg){
  double *gvals_out = calloc(size, sizeof(double));
  double max_gval_out, min_gval_out;
  srand(time(NULL));

  max_gval_out = bitpix < 0 ? 1 : pow(2, bitpix)-1;

  if(neg){
    max_gval_out /= 2;
    min_gval_out = - max_gval_out-1;
  } else
    min_gval_out = 0;
  
  printf("Creating a random image/volume of size %lu and %d bit-per-pixel depth. \n Min and max doubles:[%lf, %lf) \n", size, bitpix, min_gval_out,  max_gval_out);

  for(ulong i = 0; i < size; i++){
    ulong maxrand = 0;//=   rand()*(RAND_MAX)*(RAND_MAX) + rand()*(RAND_MAX) +rand();
    for (int j=LOOP_COUNT; j > 0; j--) {
      maxrand = maxrand*(RAND_MAX + (uint64_t)1) + rand();
    }
    float maxrandd =  (float)((double)rand()/(double)(RAND_MAX));
      gvals_out[i] = bitpix > 0 ? min_gval_out + maxrand % (ulong) max_gval_out : min_gval_out + max_gval_out*maxrandd;
  }
  return gvals_out;
  }*/


double *create_sequence(ulong size, int bitpix, int neg){
  double *gvals_out = calloc(size, sizeof(double));
  double max_gval_out, min_gval_out;

  max_gval_out = bitpix < 0 ? 1 : MIN(size, pow(2, bitpix)-1);

  if(neg){
    max_gval_out /= 2;
    min_gval_out = - max_gval_out-1;
  } else
    min_gval_out = 0;
  

  printf(" Creating an image/volume with a linear sequence of pixel intensities of %lu pixels and %d bit-per-pixel depth.\n Min and max doubles: [%f, %f) \n", size, bitpix, (float) min_gval_out, (float) max_gval_out);

  for(ulong i = 0; i < size; i++){
    if(bitpix > 0)
      gvals_out[i] = i < max_gval_out ?  (double) min_gval_out + i :  (double) max_gval_out;
    else
      gvals_out[i] = (double) (max_gval_out * i)/size ;
  }
  return gvals_out;
}

double *create_ioniz(ulong size, float fract, ulong dims[3]){
  double *gvals_out = calloc(size, sizeof(double));
  srand(time(NULL));
  
  ulong npixels   = (ulong) (fract * size);
  ulong npixels_c = 0;

  while(npixels_c != npixels){
    ulong maxrand = 0;
    for (int j=LOOP_COUNT; j > 0; j--) {
      maxrand = maxrand*(RAND_MAX + (uint64_t)1) + rand();
    }
    double myRand = rand()/(1.0 + RAND_MAX); 
    double range = size + 1;
    ulong myRand_scaled = myRand * range;
    double z = myRand_scaled / (dims[0]*dims[1]);
    double y = myRand_scaled % (dims[0]*dims[1]) / dims[0];
    double x = myRand_scaled % dims[0];

    if	(gvals_out[(ulong) (x+y*dims[0]+z*dims[0]*dims[1])] == 1) continue;
    else{
      gvals_out[(ulong) (x+y*dims[0]+z*dims[0]*dims[1])]= 1;
      npixels_c++;
    }
      
  }
  return gvals_out;
}


double *create_bubbles_ov(ulong dim, ulong size, float fract, int periodic){
  double *gvals_out = calloc(size, sizeof(double));
  srand(time(NULL));
  ulong npixels = (ulong) (fract * size);
  ulong npixels_c = 0;
  while(npixels_c < npixels){
    double myRand = rand()/(1.0 + RAND_MAX); 
    int range = size + 1;
    int myRand_scaled = myRand * range;
    double r = round(randn(5, 3));
    // double r = 10;
    ulong x0 = myRand_scaled%dim;
    ulong z0 = myRand_scaled/(dim*dim);
    ulong y0 = myRand_scaled%(dim*dim)/dim;
    ulong x,y,z,rr=r*r,dx,dy,dz,xp,yp,zp;
    if ((x0 - r < 0 || x0 + r >= dim) && !periodic) continue;
    if ((y0 - r < 0 || y0 + r >= dim) && !periodic) continue;
    if ((z0 - r < 0 || z0 + r >= dim) && !periodic) continue;

    for (z=z0-r;z<=z0+r;z++){
      dz = z-z0;
      zp = z;
      if (z < 0 && periodic) zp = dim + z;
      else if (z >= dim && periodic) zp = z - dim;
      for (y=y0-r;y<=y0+r;y++){
	dy = y-y0;
	yp = y;
	if (y < 0 && periodic) yp = dim + y;
	else if (y >= dim && periodic) yp = y - dim;
	for(x=x0-r; x<=x0+r;x++){
	  dx = x-x0;
	  xp = x;
	  if (x < 0 && periodic) xp = dim + x;
	  else if (x >= dim && periodic) xp = x - dim;
	  if(dx*dx+dy*dy+dz*dz <= r*r){
	    npixels_c++;
	    gvals_out[zp*dim*dim+yp*dim+xp]=255;
	  }
	}
      }
    }
  }
  return gvals_out;
}



double *create_bubbles_nov(ulong dim, ulong size,  float fract, int periodic){
  double *gvals_out = calloc(size, sizeof(double));
  //memset(gvals_out, 1, size * sizeof(double));
  srand(time(NULL));
  ulong npixels = (ulong) (fract * size);
  ulong npixels_c = 0;
  while(npixels_c < npixels){
    double myRand = rand()/(1.0 + RAND_MAX); 
    int range = size + 1;
    int myRand_scaled = myRand * range;
    double r = round(randn(5, 3));
    // double r = 10;
    long x0 = myRand_scaled%dim;
    long z0 = myRand_scaled/(dim*dim);
    long y0 = myRand_scaled%(dim*dim)/dim;
    long x,y,z,rr=r*r,dx,dy,dz,xp,yp,zp;
    bool ok = 1;
    if ((x0 - r < 0 || x0 + r >= dim) && !periodic) continue;
    if ((y0 - r < 0 || y0 + r >= dim) && !periodic) continue;
    if ((z0 - r < 0 || z0 + r >= dim) && !periodic) continue;
    if ((x0 - r -2< 0 || x0 + r +2>= dim) && !periodic) continue;
    if ((y0 - r -2< 0 || y0 + r +2>= dim) && !periodic) continue;
    if ((z0 - r -2< 0 || z0 + r +2>= dim) && !periodic) continue;
    for (z=z0-r-2;z<=z0+r+2;z++){
      dz = z-z0;
      zp = z;
      if (z < 0 && periodic) zp = dim + z;
      else if (z >= dim && periodic) zp = z - dim;
      for (y=y0-r-2;y<=y0+r+2;y++){
	dy = y-y0;
	yp = y;
	if (y < 0 && periodic) yp = dim + y;
	else if (y >= dim && periodic) yp = y - dim;
	for(x=x0-r-2; x<=x0+r+2;x++){
	  dx = x-x0;
	  xp = x;
	  if (x < 0 && periodic) xp = dim + x;
	  else if (x >= dim && periodic) xp = x - dim;
	  if(gvals_out[zp*dim*dim+yp*dim+xp]) ok = 0;
	}
      }
    }
    if(ok){
      for (z=z0-r;z<=z0+r;z++){
	dz = z-z0;
	zp = z;
	if (z < 0 && periodic) zp = dim + z;
	else if (z >= dim && periodic) zp = z - dim;
	for (y=y0-r;y<=y0+r;y++){
	  dy = y-y0;
	  yp = y;
	  if (y < 0 && periodic) yp = dim + y;
	  else if (y >= dim && periodic) yp = y - dim;
	  for(x=x0-r; x<=x0+r;x++){
	    dx = x-x0;
	    xp = x;
	    if (x < 0 && periodic) xp = dim + x;
	    else if (x >= dim && periodic) xp = x - dim;
	    if(dx*dx+dy*dy+dz*dz <= r*r){
	      npixels_c++;
	      gvals_out[zp*dim*dim+yp*dim+xp]=255;
	    }
	  }
	}
      }
    }
  }
  return gvals_out;
}



double *modify_bpp(double *gvals_in, ulong size, int bitpix, int nextbitpix){
  srand(time(NULL));
  double *gvals_out = calloc(size, sizeof(double));
  double nlevels_out = nextbitpix < 0 ? FLT_MAX : pow(2, nextbitpix);
  printf("Number levels in new image %lf\n",(float) nlevels_out);
  
  double max_gval_in, min_gval_in, min_gval_out = 0, max_gval_out = 0;
  max_gval_in = find_max(gvals_in, size);
  min_gval_in = find_min(gvals_in, size);
  double nlevels_in = bitpix < 0 ? max_gval_in - min_gval_in: pow(2, bitpix);
  max_gval_out = nextbitpix > 0 ? pow(2, nextbitpix)-1 : FLT_MAX;
  printf("Previous image: size %lu, bit-per-pixel %d, min and max doubles = [%lf,  %lf] \n", size, bitpix, min_gval_in, max_gval_in);
  printf("New image: size %lu, bit-per-pixel %d, min and max doubles = [%lf,  %lf] \n", size, nextbitpix,min_gval_out, max_gval_out);
  
  for(size_t i = 0; i<size; i++){
    gvals_out[i] = (double) (gvals_in[i]-min_gval_in)*(nlevels_out-1)/(nlevels_in);
    if(bitpix > 0 && nlevels_out > nlevels_in)
      gvals_out[i] += (rand()*RAND_MAX*RAND_MAX*RAND_MAX + rand()*RAND_MAX*RAND_MAX + rand()*RAND_MAX + rand())% (uint) (nlevels_out/nlevels_in);
  }

  return gvals_out;
}

    
int compare_txt_files(FILE *fp1, FILE *fp2)
{
  char ch1 = getc(fp1);
  char ch2 = getc(fp2);
  int error = 0, pos = 0, line = 1;
  while (ch1 != EOF && ch2 != EOF)
    {
      pos++;

      if (ch1 == '\n' && ch2 == '\n')
        {
	  line++;
	  pos = 0;
        }

      if (ch1 != ch2)
        {
	  error++;
	  printf("Line Number : %d \t Error "
		 "Position : %d \n", line, pos);
        }
      ch1 = getc(fp1);
      ch2 = getc(fp2);
    }

  if(error)
    printf("\x1B[91m Errors spotted \033[0m\n" );
  else
    printf("\x1B[92m No differences spotted \033[0m\n");

  return error;
}

double *cross_corr(double *imgA, double *imgB, ulong size){
  double mx,my,sx,sy,sxy,denom,r;
  double *corr_map = calloc(size, sizeof(double));
  /* Calculate the mean of the two series x[], y[] */
  mx = 0;
  my = 0;   
  for (ulong i=0;i<size;i++) {
    mx += imgA[i];
    my += imgB[i];
  }
  mx /= size;
  my /= size;

  /* Calculate the denominator */
  sx = 0;
  sy = 0;
  
  for (ulong i=0;i<size;i++) {
    sx += (imgA[i] - mx) * (imgA[i] - mx);
    sy += (imgB[i] - my) * (imgB[i] - my);
  }
  
  denom = sqrt(sx*sy);

  /* Calculate the correlation series */
  //  for (delay=-maxdelay;delay<maxdelay;delay++) {
  sxy = 0;
  for (ulong i=0;i<size;i++) {
    // j = i;
    //  if (j < 0 || j >= n)
    //	continue;
    //  else
    sxy += (imgA[i] - mx) * (imgB[i] - my);
    corr_map[i] = (imgA[i] - mx) * (imgB[i] - my);
    /* Or should it be (?)
       if (j < 0 || j >= n)
       sxy += (x[i] - mx) * (-my);
       else
       sxy += (x[i] - mx) * (y[j] - my);
    */
  }
  r = sxy / denom;
  printf("Correlaton coefficient if %lf \n", r);
  return corr_map;
  /* r is the correlation coefficient at "delay" */
    
}

double *gaussian_filter(value *gvals, ulong *dims, double fhwm, ulong size){
  // set standard deviation to 1.0
  double sigma = fhwm/2.355;
  double r, s = 2.0 * sigma * sigma;
  // sum is for normalization
  double sum = 0.0;

  double kernel[size][size][size];
  long width, height, depth;
  width = height  = size/2;
  depth = dims[2] > 1 ? size/2 : 0;
  // generate 5x5 kernel

  for (long x = -depth; x <= depth; x++)
    {
      for(long y  = -height; y <= height; y++)
        {
	  for(long z = -width; z<= width; z++)
	    {
	      r = sqrt(x*x + y*y + z*z);
	      //   printf("IN%lf, %lf, %lf \n", x + width, y + height, z + depth);

	      kernel[x + depth][y + height][z + width] = (exp(-(r*r)/s))/(M_PI * s);
	      //  printf("%lf \n", kernel[x + width][y + height][z + depth]);

	      sum += kernel[x + depth][y + height][z + width];
	    }
        }
    }

  // normalize the Kernel
   for (long x = 0; x < (long) size; x++)
    for(long y = 0; y <(long) size; y++)
      for(long z = 0; z < (long)size; z++)
            kernel[x][y][z] /= sum;

   for(long i = 0; i <(long) size; ++i){
     for (long j = 0; j <(long) size; ++j){
      for (long z = 0; z < (long)size; ++z)
	printf(" %lf \t",  kernel[i][j][z]);
      printf("\n");
     }
     printf("\n\n");
   }
   printf("\n\n");

   double *newdata =  calloc(dims[0]*dims[1]*dims[2], sizeof(double));
   //printf("ERROR ?\n");
   #pragma omp parallel for
   for(long i = 0; i < dims[2]; i++){
     // printf("%ld\n", i);
     for(long j = 0; j < dims[1]; j++){
       for(long k = 0; k < dims[0]; k++){
	 long startz = (long) i - depth < 0 ? 0 : i - depth;
	 // printf("ERROR ? %ld\n", startz);

	 long starty = (long) j - height < 0 ? 0 : j - height;
	 // printf("ERROR ? %ld\n", starty);

   	 long startx = (long) k - width < 0 ? 0 : k - width;
	 // printf("ERROR ? %ld \n", startx);

	 long endz = (long) i + depth >= dims[2] ? (long) dims[2] : i + (long)depth;
	 // printf("ERROR ?%ld\n", endz);

	 long endy = (long) j + height >= dims[1] ? (long) dims[1] : j + (long)height;
	 //printf("ERROR ? %ld\n", endy);

	 long endx = (long) k + width >= dims[0] ? (long) dims[0] : k + (long)width;
	 //printf("ERROR ?%ld \n", endx);

	 //	/ printf("%ld %ld %ld\n",  startx, starty, startz);
	 long idxfix = i *dims[0]*dims[1] + j *dims[0] + k;

	 for(long z = startz; z < endz; z++){
	   for(long y = starty; y < endy; y++){
	     for(long x = startx; x < endx; x++){
	       long idx = z *dims[0]*dims[1] + y *dims[0] + x;
	       //  printf("%ld, %ld %ld \n", x,y,z);
	       newdata[idxfix] += (double) gvals[idx] * kernel[z-startz][y-starty][x-startx];
	     }
	   }
	 }

       }
     }
   }

   return newdata;

}

double *hessian(value *gvals, ulong *dims, ulong p){
  printf("%ld %lf \n", p, gvals[p]);

  double hxx, hyy, hzz, hxy, hxz, hyz;
  long x = p %dims[0];
  long y = (p %(dims[0]*dims[1]))/dims[0];
  long z = (p /(dims[0]*dims[1]));
  // info("x %ld, y %ld z %ld", x, y, z);

  long xbef = (x - 1) >= 0? x-1: (long) dims[0]-1;
  long xaf  = (x + 1) < (long) dims[0]? x+1: 0;
  long ybef = (y - 1) >= 0? y-1: (long) dims[1]-1;
  long yaf  = (y + 1) < (long) dims[1]?  y+1: 0;
  long zbef = (z - 1) >= 0? z-1: (long) dims[2]-1;
  long zaf  = (z + 1) < (long) dims[2]?  z+1: 0;
  // info("x %ld, y %ld z %ld", xbef, ybef, zbef);
  // info("x %ld, y %ld z %ld", xaf, yaf, zaf);

  hxx = (double) (gvals[z*(dims[0]*dims[1])+y*dims[0]+xbef] + gvals[z*(dims[0]*dims[1])+y*dims[0]+xaf] - 2*gvals[p]);
  // printf("%lf %lf %lf\n", gvals[z*(dims[0]*dims[1])+y*dims[0]+xbef], gvals[z*(dims[0]*dims[1])+y*dims[0]+xaf], 2*gvals[p]);

  hyy = (double) (gvals[z*(dims[0]*dims[1])+ybef*dims[0]+x] + gvals[z*(dims[0]*dims[1])+yaf*dims[0]+x] - 2*gvals[p]);
  hzz = (double) (gvals[zbef*(dims[0]*dims[1])+y*dims[0]+x] + gvals[zaf*(dims[0]*dims[1])+y*dims[0]+x] - 2*gvals[p]);

  hxy = (double) (gvals[z*(dims[0]*dims[1])+ybef*dims[0]+xbef] + gvals[z*(dims[0]*dims[1])+yaf*dims[0]+xaf] - gvals[z*(dims[0]*dims[1])+yaf*dims[0]+xbef] -gvals[z*(dims[0]*dims[1])+ybef*dims[0]+xaf]);
  hxz = (double) (gvals[zbef*(dims[0]*dims[1])+y*dims[0]+xbef] + gvals[zaf*(dims[0]*dims[1])+y*dims[0]+xaf] - gvals[zaf*(dims[0]*dims[1])+y*dims[0]+xbef] -gvals[zbef*(dims[0]*dims[1])+y*dims[0]+xaf]);
  hyz = (double) (gvals[zbef*(dims[0]*dims[1])+ybef*dims[0]+x] + gvals[zaf*(dims[0]*dims[1])+yaf*dims[0]+x] - gvals[zaf*(dims[0]*dims[1])+ybef*dims[0]+x] -gvals[zbef*(dims[0]*dims[1])+yaf*dims[0]+x]);
  double *arr = calloc(6, sizeof(double));
  arr[0]=hxx;
  arr[1]=hyy;
  arr[2]=hzz;
  arr[3]=hxy;
  arr[4]=hyz;
  arr[5]=hxz;
  return arr;
}
