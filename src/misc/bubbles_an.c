#include "common.h"
#include "eispack.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <sys/time.h>


char *remove_ext(char* mystr);
double *inertia_analysis(char *prefix, double *matrix, ulong *volume, double *intens, double *eigvec, double *eigval, ulong nbubbles);
void orientation_analysis(char *prefix, double *vec, double *center, ulong *volume, ulong *volume_th, ulong *thresh, ulong *dims, double *intens, ulong nbubbles, ulong nvol, double maxint[2]);
void density_bubbles(char *prefix, double *matrix, double *matrix_aux, double *vec_main, double *vec_aux,  ulong nbubbles);
int find_scale_double(double *thresh, double val, int nthresh);
void twoptcorr(char *prefix, ulong *volume, ulong volume_th, double *intens, int nbins, int boots, double maxint[2]);
double frob(double matrix[3][3], int size1, int size2);

int main(int argc, char** argv) {
  char cmnd[10000];
  srand(time(NULL)); 

  if(argc <= 1){
    printf("HII regions analysis function:  %s <analysis> : {run, all, inertia, orientation, twoptcorr, density} \n", argv[0]);
    exit(0);
  }else{

    /**************************************************************************************/

    /******************************  Run the program  *************************************/

    /**************************************************************************************/
    if(!strcmp(argv[1], "run")){
      if(argc < 9){
	printf("Running disccoman: <input name> <input type> <output name> <dim size> <lambda> <tree> <periodic> [bpp] [density] \n");
	exit(0);
      }
      ulong dims[3] = {atoi(argv[5]),atoi(argv[5]),atoi(argv[5])};
      double lambda = atof(argv[6]);
      sprintf(cmnd, "mpirun -np 1 ./disccoman -g 1,1,1 -f eor -l %d --refine 0 --inprefix %s -v info --intype %s --attribute 4 -c 26 --outprefix out/%s --dims %d,%d,%d --tree %s --periodic 1 ",  lambda, argv[2], argv[3], argv[4],dims[0],dims[1],dims[2], argv[7], atoi(argv[8]));
      if(argc == 10){
	int bpp = atoi(argv[9]);
	sprintf(cmnd, "%s --bpp %d ", cmnd, bpp);
      }
      if(argc == 11)
	sprintf(cmnd, "%s --density %d ", cmnd, argv[10]);
      
      printf("Running disccoman with the following parameters: \n %s \n", cmnd);
      system(cmnd);
    }
    /**************************************************************************************/

    /******************************  Inertia   *************************************/

    /**************************************************************************************/
    else if(!strcmp(argv[1], "inertia")){
    
      if(argc < 4){
	printf("Running inertia computation: <input name> <dims> <orientation> <twoptcorr> <density> [maxint]\n");
	exit(0);
      }
      printf("/*****************************************/ \n");
      printf("/*** Computing inertia attributes      ***/ \n");
      printf("/*****************************************/ \n");
      char *fname, *prefix;
      ulong nbubbles;
      ulong dims[3] = {atoi(argv[3]),atoi(argv[3]),atoi(argv[3])};
      double maxint[2] = {0};
      if(argc == 9){
	maxint[0] = atof(argv[7]);
	maxint[1] = atof(argv[8]);
      }
      asprintf(&fname, "%s_ine.bin", argv[2]);
      printf("Opening inertia file %s \n", fname);
      FILE *infile = fopen(fname, "rb");
      if(infile == NULL)
	printf("File not opened, issue in name ? \n");
		      
      fseek(infile, 0, SEEK_END);
      nbubbles=ftell(infile)/(6*8);
      fseek(infile, 0, SEEK_SET);
  
      printf("Number of objects detected : %lu\n", nbubbles);

      double *matrix  = malloc(6*nbubbles  * sizeof(double));
      fread(matrix, nbubbles*6*sizeof(double), 1, infile);
      fclose(infile);
  
      printf("File %s closed \n", fname);
    
      ulong  *volume  = malloc(nbubbles * sizeof(ulong));
      asprintf(&fname, "%s_vol.bin", argv[2]);
      printf("Opening volume file %s \n", fname);
      infile = fopen(fname, "rb");
      if(infile == NULL)
	printf("File not opened, issue in name ? \n");
      fread(volume, nbubbles*sizeof(ulong), 1, infile);
      fclose(infile);
      printf("File %s closed \n", fname);


      double *intens  = malloc(2*nbubbles  * sizeof(double));
      printf("Opening intensity file %s \n", fname);
      asprintf(&fname, "%s_int.bin", argv[2]);
      infile = fopen(fname, "rb");
      if(infile == NULL)
	printf("File not opened, issue in name ? \n");
      fread(intens, 2*nbubbles*sizeof(double), 1, infile);
      fclose(infile);
      printf("File %s closed \n", fname);

      double *eigval  = malloc(3*nbubbles  * sizeof(double));
      double *eigvec  = malloc(9*nbubbles  * sizeof(double));
  
      #pragma omp parallel for
      for (ulong i = 0; i < nbubbles; i++) {
	double tens_mat[9] = {matrix[i*6], matrix[i*6+3], matrix[i*6+5], matrix[i*6+3], matrix[i*6+1],matrix[i*6+4], matrix[i*6+5], matrix[i*6+4],matrix[i*6+2]};
	int ierr = rs (3, tens_mat, eigval+i*3, 1, eigvec+i*9);
      }
      printf("Eigenvalue and eigenvectors computed \n");

      asprintf(&prefix, "%s", argv[2]);
      double *vec_main = inertia_analysis(prefix, matrix, volume, intens, eigvec, eigval, nbubbles);
      printf("Inertia attributes of the objects computed \n");

      if(atoi(argv[4])){
	double *center  = malloc(3*nbubbles  * sizeof(double));
	printf("Opening position file %s \n", fname);
	asprintf(&fname, "%s_cen.bin", argv[2]);
	infile = fopen(fname, "rb");
	if(infile == NULL)
	  printf("File not opened, issue in name ? \n");
	fread(center, 3*nbubbles*sizeof(double), 1, infile);
	fclose(infile);
	printf("File %s closed \n", fname);

	ulong volume_th[3]  = {10,100,1000};
	ulong thresh[3]  = {20,30,50};
	orientation_analysis(prefix, vec_main, center,volume,volume_th, thresh, dims, intens, nbubbles, 3, maxint);
	printf("Cross correlations between orientation vectors done \n");
      }

      if(atoi(argv[5])){
	ulong volume_thc[3]  = {10,100,1000};
	for(int i = 0; i <3; i++)
	  twoptcorr(prefix, volume, volume_thc[i], intens, 100, 20, maxint);
	printf("Two pt correlation function done \n");
      }
      
      if(atoi(argv[6])){
	printf("\n Running and comparing the density file \n");

	asprintf(&fname, "%s_dens_ine.bin", argv[2]);
	FILE *infile = fopen(fname, "rb");
	if(infile == NULL)
	  printf("File not opened, issue in name ? \n");
	double *matrix_aux  = malloc(6*nbubbles  * sizeof(double));
	fread(matrix_aux, nbubbles*6*sizeof(double), 1, infile);
	fclose(infile);
	printf("File %s closed \n", fname);

	double *intens_dens  = malloc(2*nbubbles  * sizeof(double));
	asprintf(&fname, "%s_dens.bin", argv[2]);
	infile = fopen(fname, "rb");
	if(infile == NULL)
	  printf("File not opened, issue in name ? \n");
	fread(intens, 2*nbubbles*sizeof(double), 1, infile);
	fclose(infile);
	printf("File %s closed \n", fname);

  
	#pragma omp parallel for
	for (ulong i = 0; i < nbubbles; i++) {
	  double tens_mat[9] = {matrix_aux[i*6], matrix_aux[i*6+3], matrix_aux[i*6+5], matrix_aux[i*6+3], matrix_aux[i*6+1],matrix_aux[i*6+4], matrix_aux[i*6+5], matrix_aux[i*6+4],matrix_aux[i*6+2]};
	  //  printf("%lf, %lf, %lf, %lf, %lf, %lf\n", matrix_aux[i*6], matrix_aux[i*6+3], matrix_aux[i*6+5], matrix_aux[i*6+1],matrix_aux[i*6+4], matrix_aux[i*6+2]);
	  int ierr = rs (3, tens_mat, eigval+i*3, 1, eigvec+i*9);
	}
    
	printf("[Density] Eigenvalue and eigenvectors computed  \n");

	asprintf(&prefix, "%s_dens", argv[2]);
    
	double *vec_aux = inertia_analysis(prefix, matrix_aux, volume, intens_dens, eigvec, eigval, nbubbles);
	printf("[Density] Inertia attributes of the objects computed\n");

	density_bubbles(prefix, matrix, matrix_aux, vec_main, vec_aux,  nbubbles);
	printf("[Density] Cross correlations computed \n");
	free(vec_aux);
    
	asprintf(&fname, "%s_dens_rand_ine.bin", argv[2]);
	/*	infile = fopen(fname, "rb");
	fread(matrix_aux, nbubbles*6*sizeof(double), 1, infile);
	fclose(infile);

	double *intens_dens_rand  = malloc(2*nbubbles  * sizeof(double));
	asprintf(&fname_int, "out/%s_dens_rand.bin", argv[2]);
	infile = fopen(fname_int, "rb");
	fread(intens, 2*nbubbles*sizeof(double), 1, infile);
	fclose(infile);
    
	#pragma omp parallel for
	for (ulong i = 0; i < nbubbles; i++) {
	  double tens_mat[9] = {matrix_aux[i*6], matrix_aux[i*6+3], matrix_aux[i*6+5], matrix_aux[i*6+3], matrix_aux[i*6+1],matrix_aux[i*6+4], matrix_aux[i*6+5], matrix_aux[i*6+4],matrix_aux[i*6+2]};
	  int ierr = rs (3, tens_mat, eigval+i*3, 1, eigvec+i*9);
	}
    
	printf("Eigenvalue and eigenvectors computed 2\n");

	asprintf(&prefix, "out/%s_dens_rand", argv[2]);
    
	vec_aux = inertia_analysis(prefix, matrix_aux, volume, intens_dens_rand, eigvec, eigval, nbubbles);
	printf("Inertia attributes for the aux cube written 2\n");

	printf("Comparing orientation in main and aux files  2\n");
	density_bubbles(prefix, matrix, matrix_aux, vec_main, vec_aux, nbubbles);
	printf("Density correlations done 2\n");

	free(vec_aux);
	free(matrix_aux);*/

      }
  
      free(vec_main);
      free(matrix);
  
      printf("/***********************/\n");
      printf("/**** Analysis done ****/\n");
      printf("/***********************/\n");

    }
  }
  exit(0);
}


double *inertia_analysis(char *prefix, double *matrix, ulong *volume, double *intens, double *eigvec, double *eigval, ulong nbubbles){
  char *fname;
  FILE *outelong = malloc(1 * sizeof(FILE*));
  asprintf(&fname, "%s_elong.bin", prefix);
  outelong = fopen(fname, "wb");
  
  FILE *outflat   = malloc(1 * sizeof(FILE*));
  asprintf(&fname, "%s_flat.bin", prefix);
  outflat = fopen(fname, "wb");

  FILE *outspars  = malloc(1 * sizeof(FILE*));
  asprintf(&fname, "%s_spars.bin", prefix);
  outspars = fopen(fname, "wb");

  FILE *outncomp  = malloc(1 * sizeof(FILE*));
  asprintf(&fname, "%s_ncomp.bin", prefix);
  outncomp = fopen(fname, "wb");

  FILE *outvec   = malloc(1 * sizeof(FILE*));
  asprintf(&fname, "%s_vec.bin", prefix);
  outvec = fopen(fname, "wb");


  if (outelong==NULL || outflat == NULL || outspars == NULL || outncomp == NULL)  exit(-1);  

  double *elong  = malloc(nbubbles  * sizeof(double));
  double *flat   = malloc(nbubbles  * sizeof(double));
  double *spars  = malloc(nbubbles  * sizeof(double));
  double *ncomp  = malloc(nbubbles  * sizeof(double));
  double *vec    = malloc(nbubbles*3* sizeof(double));

  #pragma omp parallel for
 for (ulong i = 0; i < nbubbles; i++) {
   ncomp[i]    = (matrix[i*6]+matrix[i*6+1]+matrix[i*6+2])/pow(intens[i*2], 5.0/3.0);
   elong[i]    = eigval[i*3 + 2]/eigval[i*3 + 1];
   flat[i]     = eigval[i*3 + 1]/eigval[i*3 + 0];
   //  if ( eigval[i*3 + 0] < 0 || eigval[i*3 + 1] <0 || eigval[i*3 + 2] < 0) 
   // printf("%lf\n", eigval[i*3 + 0], eigval[i*3 + 1], eigval[i*3 + 2]);
   double ax_len[3] = {sqrt(20*eigval[i*3 + 2]/intens[i*2]), sqrt(20*eigval[i*3 + 1]/intens[i*2]), sqrt(20*eigval[i*3 + 0]/intens[i*2])};
   spars[i]     = intens[i*2]/(double) volume[i] * 3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*intens[i*2]);
   int pic = 2;
   if(eigval[i*3+2] == eigval[i*3+1])
    pic = eigval[i*3+1] == eigval[i*3] ? rand() % 3: rand() % 2 + 1;
   vec[3*i] = eigvec[i*9+pic*3];
   vec[3*i+1] = eigvec[i*9+pic*3+1];
   vec[3*i+2] = eigvec[i*9+pic*3+2]; 
 }

 #pragma omp parallel
 #pragma omp single
 {
   #pragma omp task
   fwrite(elong,nbubbles*sizeof(double),1,outelong);
   #pragma omp task
   fwrite(flat,nbubbles*sizeof(double),1,outflat);
   #pragma omp task
   fwrite(spars,nbubbles*sizeof(double),1,outspars);
   #pragma omp task
   fwrite(ncomp,nbubbles*sizeof(double),1,outncomp);
   #pragma omp task
   fwrite(vec,3*nbubbles*sizeof(double),1,outvec);
   #pragma omp taskwait
 }

 free(elong);
 free(flat);
 free(spars);
 free(ncomp);
 fclose(outelong);
 fclose(outflat);
 fclose(outspars);
 fclose(outncomp);
 fclose(outvec);
 return vec;
}

int find_scale(float *thresh, float val, int nthresh){
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

void orientation_analysis(char *prefix, double *vec, double *center, ulong *volume, ulong *volume_th, ulong *thresh, ulong *dims, double *intens, ulong nbubbles, ulong nvol, double maxint[2]){
  
  FILE **outfile = malloc(nvol*sizeof(FILE*)); 
  char *temp ;
  uint nbins    = 20;

  printf("Maximum number of cross correlations to perform: %lu\n", nbubbles*(nbubbles-1)/2);

  float ***cos_arr = calloc(nbins, sizeof(float*));
  uint  ***sca_arr = calloc(nbins, sizeof(uint*));
  for(int i = 0; i < nbins; i++){
    
  ulong *alloc    = calloc(nvol, sizeof(ulong));
  float **dist_bins = calloc(nvol, sizeof(float*));
  ulong *nbub     = calloc(nvol, sizeof(ulong));

  
  for(int i = 0; i < nvol; i++){
    asprintf(&temp, "%s_orien_V%lu.csv", prefix, volume_th[i]);
    outfile[i] = fopen(temp, "wb");
    if (outfile[i]==NULL)  exit(-1);
    alloc[i]    = 1000000;
    dist_bins[i] = calloc(nbins, sizeof(float));
    cos_arr[i]   = calloc(alloc[i], sizeof(float));
    sca_arr[i]   = calloc(alloc[i], sizeof(uint));
    dist_bins[i][1] =  thresh[i];
    dist_bins[i][2] =  thresh[i] + thresh[i]/2;
    for(uint j = 3; j < nbins; j++){
      double dr = dist_bins[i][j-1] - dist_bins[i][j-2];
      dist_bins[i][j] = dist_bins[i][j-1] + pow(dist_bins[i][j-1],2)*dr/(pow(dist_bins[i][j-1]+dr,2));
    }
  }

  for (ulong i = 0; i < nbubbles; i++) {
    for(ulong j = i+1; j < nbubbles; j++) {
      float deltax = fabs(center[3*i]-center[3*j]) > (float) dims[0]/2.0 ?
	(float)	dims[0] - fabs(center[3*i]-center[3*j]):
	center[3*i]-center[3*j];
      float deltay = fabs(center[3*i+1]-center[3*j+1]) > (float) dims[1]/2.0 ?
	(float)	dims[1] - fabs(center[3*i+1]-center[3*j+1]):
	center[3*i]-center[3*j];
      float deltaz = fabs(center[3*i+2]-center[3*j+2]) > (float) dims[2]/2.0 ?
	(float)	dims[2] - fabs(center[3*i+2]-center[3*j+2]):
	center[3*i+2]-center[3*j+2];
      float dist = sqrt(pow(deltax,2)+pow(deltay,2)+pow(deltaz,2));
      
      for(int k = 0; k < nvol; k++){
	if(volume[i] >= volume_th[k] && volume[j] >= volume_th[k]){
	  int scale = find_scale(dist_bins[k], dist, nbins);
	  //  printf("Total max cross correlations to perform: %lu\n", nbub[k]);
	  if(scale == nbins - 1) continue;
	  if(nbub[k] == alloc[k] -1){
	    alloc[k] *= 1.5;
	    cos_arr[k] = realloc(cos_arr[k], alloc[k]*sizeof(float));
	    sca_arr[k] = realloc(sca_arr[k], alloc[k]*sizeof(uint));
	  }
	  cos_arr[k][nbub[k]] = 2*pow(vec[i*3]*vec[j*3] + vec[i*3+1]*vec[j*3+1] + vec[i*3+2]*vec[j*3+2],2) - 1;
	  sca_arr[k][nbub[k]++] = scale;
	}
      }
    }
  }
  printf("Cross corr computed, starting bootstrapping\n");
  ulong *pairs    ;// = calloc(nbins, sizeof(ulong));
  float *cos_sum  ;//= calloc(nbins, sizeof(float));

      
  for(int k = 0; k < nvol; k++){
    pairs     = calloc(nbins, sizeof(ulong));
    cos_sum   = calloc(nbins, sizeof(float));
    for(ulong j = 0; j < nbub[k]; j++){
	cos_sum[sca_arr[k][j]] += cos_arr[k][j];
	pairs[sca_arr[k][j]]++;
    }
    
    for (ulong i = 0; i < nbins; i++) {
      if (pairs[i] < 100 && i != nbins - 1){
	pairs[i+1]   += pairs[i];
	cos_sum[i+1] += cos_sum[i];
      } else if (pairs[i] > 0){
	fprintf(outfile[k], "%10.2f,%10lu,%10.3f,", dist_bins[k][i],pairs[i],cos_sum[i]/pairs[i]);
	#pragma omp parallel
	{
	  const gsl_rng_type * T;
	  gsl_rng * r;
	  gsl_rng_env_setup();
	  struct timeval tv; // Seed generation based on time
	  gettimeofday(&tv,0);
	  unsigned long mySeed = tv.tv_sec + tv.tv_usec + omp_get_thread_num();
	  T = gsl_rng_default; // Generator setup
	  r = gsl_rng_alloc (T);
	  gsl_rng_set(r, mySeed);

	  uint seed = 25234 + 17*omp_get_thread_num();
	  #pragma omp for 
	  for(ulong j = 0; j < 100; j++){
	    float sum = 0;
	    for(ulong t = 0; t  < pairs[i]; t++){
      
	      ulong id = gsl_rng_uniform_int(r, nbub[k]);//rand_64(&seed) % nbub[k];
	      sum += cos_arr[k][id];
	    }
	    fprintf(outfile[k], "%10.3f,", sum / pairs[i]);
	  }
	}
	fprintf(outfile[k], "\n");
      }
    }
    fclose(outfile[k]);
    free(pairs);
    free(cos_sum);
  }
  free(alloc);
  for(int i = 0; i < nvol; i++){
    free(cos_arr[i]);
    free(sca_arr[i]);
    free(dist_bins[i]);
  }
  free(sca_arr);
  free(dist_bins);
  free(cos_arr);
}


 
void density_bubbles(char *prefix, double *matrix, double *matrix_aux, double *vec_main, double *vec_aux, ulong nbubbles){
  char *fname;
  asprintf(&fname, "%s_cos.bin", prefix);
  FILE *outfile = fopen(fname, "wb"); 
  asprintf(&fname, "%s_dcorr.bin", prefix);
  FILE *outfileb = fopen(fname, "wb"); 

  printf("Number of vectors to compare: %lu\n", nbubbles);

  float *cos_arr     = calloc(nbubbles+1, sizeof(float));
  float *dcorr_arr   = calloc(nbubbles+1, sizeof(float));
  
  for (ulong i = 0; i < nbubbles; i++){
    cos_arr[i] = 2*pow(vec_main[i*3]*vec_aux[i*3] + vec_main[i*3+1]*vec_aux[i*3+1] + vec_main[i*3+2]*vec_aux[i*3+2],2) - 1;
    double ma[3][3] = {{matrix[i*6],  matrix[i*6+3],matrix[i*6+5]},
		      {matrix[i*6+3], matrix[i*6+1],matrix[i*6+4]},
		      {matrix[i*6+5], matrix[i*6+4],matrix[i*6+2]}};
    double mb[3][3] = {{matrix_aux[i*6],   matrix_aux[i*6+3], matrix_aux[i*6+5]},
		       {matrix_aux[i*6+3], matrix_aux[i*6+1], matrix_aux[i*6+4]},
		       {matrix_aux[i*6+5], matrix_aux[i*6+4], matrix_aux[i*6+2]}};
    double mc[3][3] = {0};
    // printf("%lf %lf %lf\n", ma[0],  ma[1], ma[2]);
    // printf("%lf %lf %lf\n", mb[0], mb[1], mb[2]);

    double sum = 0;
    for (ulong c = 0; c < 3; c++) {
      for (ulong d = 0; d < 3; d++) {
	for (ulong k = 0; k < 3; k++) {
	  sum = sum + ma[c][k]*mb[k][d];
	}
	mc[c][d] = sum;
	sum = 0;
      }
    }

    double trace = 0;
    for (ulong c = 0; c < 3; c++) 
      trace += mc[c][c];

    dcorr_arr[i] = (float) 1-trace/(frob(ma, 3,3)*frob(mb,3,3));
    // fprintf(outfile, "%10.3f,", cos_arr[i]);
    // printf("%lf\n", dcorr_arr[i]);
  }
  // fprintf(outfile, "\n");
  fwrite(cos_arr, (nbubbles+1) * sizeof(float), 1, outfile);
  fwrite(dcorr_arr, (nbubbles+1) * sizeof(float), 1, outfileb);

  fclose(outfile);
  fclose(outfileb);

  free(cos_arr);
  free(dcorr_arr);
  printf("Comparison done \n");
}

void twoptcorr(char *prefix, ulong *volume, ulong volume_th, double *intens, int nbins, int boots, double maxint[2]){
  printf("Two points correlation for bubbles larger than %ld \n", volume_th);
  char *fname;
  asprintf(&fname, "%s_Position.txt", prefix);
  FILE *fp1 = fopen(fname, "r");
  if(fp1 == NULL)
    printf("Positions file not opened \n");
 
 
  ulong dims[3];
  ulong nval, nvalr, count=0;
  fscanf(fp1, "%ld\n", &nval);
  printf("%ld positions to read \n", nval);
  fscanf(fp1, "0 %ld 0 %ld 0 %ld \n", dims, dims+1, dims+2);
  printf("Full dimensions: %ld  %ld  %ld \n", dims[0],dims[1],dims[2]);

  double tempx, tempy, tempz, temp;
  double *x = calloc(nval, sizeof(double));
  double *y = calloc(nval, sizeof(double));
  double *z = calloc(nval, sizeof(double));

  for(ulong i = 0; i < nval; i++){
    fscanf(fp1,"%lf %lf %lf %lf \n",&tempx, &tempy, &tempz, &temp);
    if(volume[i] >= volume_th){
      x[count] = tempx;
      y[count] = tempy;
      z[count] = tempz;
      count++;
    }
  }
  fclose(fp1);
  if(count < 500){
    printf("Not enough bubbles \n");
    return;
  }
  nval = count;
  nvalr = 2*nval;
  asprintf(&fname, "%s_twoptcorr_V%ld.csv", prefix, volume_th);
  FILE *fp2 = fopen(fname, "w");
  if(fp2 == NULL)
    printf("Output file not opened \n");
  ulong size = dims[0]*dims[1]*dims[2];
  double *count_dd = calloc(nbins, sizeof(double));
  double **count_dr = calloc(boots, sizeof(double *));
  double **count_rr = calloc(boots, sizeof(double *));
  double **corr = calloc(boots, sizeof(double *));
  for(ulong i = 0; i<boots; i++){
    corr[i] = calloc(nbins, sizeof(double));
    count_dr[i] = calloc(nbins, sizeof(double));
    count_rr[i] = calloc(nbins, sizeof(double));
  }

  ulong max_dist = sqrt(dims[0]*dims[0] + dims[1]*dims[1] + dims[2]*dims[2]);
  double bin_val = (double) max_dist / (double) nbins;

  double *dist_bins = calloc(nbins, sizeof(double));
  for(ulong i = 0; i < nbins; i++)
    dist_bins[i] = i*bin_val;
  for (ulong i = 0; i < nval; i++) {
    for(ulong j = i+1; j < nval; j++) {
      double deltax, deltay, deltaz;
      deltax = fabs(x[i]-x[j]) > (double) dims[0]/2.0 ? (double) dims[0] -fabs(x[i] - x[j]): x[i]-x[j];
      deltay = fabs(y[i]-y[j]) > (double) dims[1]/2.0 ?(double) dims[1] - fabs(y[i] - y[j]): y[i]-y[j];
      deltaz = fabs(z[i]-z[j]) > (double) dims[2]/2.0 ? (double) dims[2] - fabs(z[i] - z[j]):	z[i]-z[j];
      double dist = sqrt(pow(deltax,2)+pow(deltay,2)+pow(deltaz,2));
      int scale = find_scale_double(dist_bins, dist, nbins);
      count_dd[scale]++;
    }
  }
  #pragma omp parallel for
  for(int k =0; k< boots; k++){
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec + omp_get_thread_num();
    T = gsl_rng_default; // Generator setup
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, mySeed);
    double *x_rand = calloc(nvalr, sizeof(double));
    double *y_rand = calloc(nvalr, sizeof(double));
    double *z_rand = calloc(nvalr, sizeof(double));
    for(ulong l = 0; l < nvalr; l++){
      x_rand[l] = gsl_rng_uniform(r)*dims[0];
      y_rand[l] = gsl_rng_uniform(r)*dims[1];
      z_rand[l] = gsl_rng_uniform(r)*dims[2];
    }
    for (ulong i = 0; i < nvalr; i++) {
      for(ulong j = 0; j < nvalr; j++) {
	double dist, deltax, deltay, deltaz;
	int scale;
	if(j > i){
	  deltax = fabs(x_rand[i]-x_rand[j]) > (double) dims[0]/2.0 ?  (double) dims[0] - fabs(x_rand[i] - x_rand[j]):   x_rand[i]-x_rand[j];
	  deltay = fabs(y_rand[i]-y_rand[j]) > (double) dims[1]/2.0 ?  (double) dims[1] - fabs(y_rand[i] - y_rand[j]):   y_rand[i]-y_rand[j];
	  deltaz = fabs(z_rand[i]-z_rand[j]) > (double) dims[2]/2.0 ?  (double) dims[2] - fabs(z_rand[i] - z_rand[j]):   z_rand[i]-z_rand[j];
	  dist = sqrt(pow(deltax,2)+pow(deltay,2)+pow(deltaz,2));
      
	  scale = find_scale_double(dist_bins, dist, nbins);
	  count_rr[k][scale]++;
	}
	if( i<nval){
	  deltax = fabs(x[i]-x_rand[j]) > (double) dims[0]/2.0 ? (double) dims[0] - fabs(x[i] - x_rand[j]):    x[i]-x_rand[j];
	  deltay = fabs(y[i]-y_rand[j]) > (double) dims[1]/2.0 ? (double) dims[1] - fabs(y[i] - y_rand[j]):   y[i]-y_rand[j];
	  deltaz = fabs(z[i]-z_rand[j]) > (double) dims[2]/2.0 ? (double) dims[2] - fabs(z[i] - z_rand[j]):  z[i]-z_rand[j];
	  dist = sqrt(pow(deltax,2)+pow(deltay,2)+pow(deltaz,2));
      	  scale = find_scale_double(dist_bins, dist, nbins);
	  count_dr[k][scale]++;
	}
      }
    }
    for(int s = 0; s < nbins; s++)
      corr[k][s] = (count_dd[s]/(nval*(nval-1)/2) - 2*count_dr[k][s]/(nvalr*nval) + count_rr[k][s]/(nvalr*(nvalr-1)/2)) / (count_rr[k][s]/(nvalr*(nvalr-1)/2));
    gsl_rng_free (r);
  }

  double *corr_mean = calloc(nbins, sizeof(double));
  for(int s = 0; s < nbins; s++){
    for(int k = 0; k < boots; k++){
      corr_mean[s] += corr[k][s];
    }
    corr_mean[s] /= boots;
  }

  double *corr_std = calloc(nbins, sizeof(double));
  for(int s = 0; s < nbins; s++){
    for(int k = 0; k < boots; k++){
      corr_std[s] += pow(corr[k][s]-corr_mean[s],2);
    }
    corr_std[s] = sqrt(1.0/(boots)*corr_std[s]);
  }
     
   for(int j = 0; j < nbins; j++){
     // printf("Corr in bin %d is %lf +- %lf\n", j, corr_mean[j], corr_std[j]);
    fprintf(fp2, "%lf, %lf, %lf\n", dist_bins[j], corr_mean[j], corr_std[j]);
  }
  
  fclose(fp2);
  free(corr_std);
  free(corr_mean);
  free(dist_bins);
  free(x);free(y);free(z);
    free(count_dd);
  for(int i = 0; i < boots; i++){
    free(count_dr[i]);
    free(count_rr[i]);
    free(corr[i]);
    }
  free(count_dr); free(count_rr), free(corr);
  
}





/*void orientation_analysis(char *prefix, double *vec, double *center, ulong *volume, ulong *volume_th, ulong *thresh, ulong nbubbles, ulong nvol){
  
  FILE **outfile = malloc(nvol*sizeof(FILE*)); 
  char *temp ;
  uint nbins    = 50;

  printf("Total max cross correlations to perform: %lu\n", nbubbles*(nbubbles-1)/2);

  float *cos_arr  = malloc(nbubbles*(nbubbles-1)/2*sizeof(float)); 
  float *dist_arr = malloc(nbubbles*(nbubbles-1)/2*sizeof(float));
  char  *vol_arr  = malloc(nbubbles*(nbubbles-1)/2*sizeof(char));
  memset(vol_arr, -1, nbubbles*(nbubbles-1)/2*sizeof(char));
  ulong s = 0;
  float **dist_bins = calloc(nvol, sizeof(float*));
  
  for(int i = 0; i < nvol; i++){
   asprintf(&temp, "%s_orien_V%lu.csv", prefix, volume_th[i]);
   outfile[i] = fopen(temp, "wb");
   if (outfile[i]==NULL)  exit(-1);  
   dist_bins[i] = calloc(nbins, sizeof(float));
   dist_bins[i][1] =  thresh[i];
   dist_bins[i][2] =  thresh[i] + thresh[i]/2;
   for(uint j = 3; j < nbins; j++){
     double dr = dist_bins[i][j-1] - dist_bins[i][j-2];
     dist_bins[i][j] = dist_bins[i][j-1] + pow(dist_bins[i][j-1],2)*dr/(pow(dist_bins[i][j-1]+dr,2));
   }
 }
  // maxdist = dist_bins[nvol-1][nbins-1];

  for (ulong i = 0; i < nbubbles; i++) {
    for(ulong j = i+1; j < nbubbles; j++, s++) {
      dist_arr[s]= sqrt(pow(center[3*i]-center[3*j],2)+pow(center[3*i+1]-center[3*j+1],2)+pow(center[3*i+2]-center[3*j+2],2));
      // if(dist_curr > maxdist) continue;
      cos_arr[s] = 2*pow(vec[i*3]*vec[j*3] + vec[i*3+1]*vec[j*3+1] + vec[i*3+2]*vec[j*3+2],2) - 1;
      for(int k = 0; k < nvol; k++)
	if(volume[i] >= volume_th[k] && volume[j] >= volume_th[k])
	  vol_arr[s] = k;
    }
  }
	  /*for(int k = 0; k < nvol; k++){
	if(volume[i] >= volume_th[k] && volume[j] >= volume_th[k]){
	  int scale = find_scale(dist_bins[k], dist_curr, nbins);
	  if(scale == nbins -1) continue;


	 nbub_bins[k]++;
	 vol_arr[s] = k;
	 // if(id_cur[k] == alloc_size-1)
	 //  id_arr[k] = realloc(id_arr[k]
	 // printf("volume %d ,k %d, scale  %d\n", vol_arr[s], k, scale);
	 // cos_sum[k][scale] += cos_arr[s];
	 // pairs[k][scale]++;
       }
       }*/
     // s++;
     // dist_arr[s] = 
     // if(25 > dist_arr[s]) continue;
  // int scale = find_scale(dist_bins, dist, nbins);


     /* double ma[3][3] = {{matrix[i*9],matrix[i*9+1],matrix[i*9+2]},
	{matrix[i*9+3], matrix[i*9+4],matrix[i*9+5]},
	{matrix[i*9+6], matrix[i*9+7],matrix[i*9+8]}};
	double mb[3][3] = {{matrix[j*9],matrix[j*9+1],matrix[j*9+2]},
	{matrix[j*9+3], matrix[j*9+4],matrix[j*9+5]},
	{matrix[j*9+6], matrix[j*9+7],matrix[j*9+8]}};
	double mc[3][3] = {0};
	double sum = 0;
	for (ulong c = 0; c < 3; c++) {
	for (ulong d = 0; d < 3; d++) {
	for (ulong k = 0; k < 3; k++) {
	sum = sum + ma[c][k]*mb[k][d];
	}
	mc[c][d] = sum;
	sum = 0;
	}
	}

	double trace = 0;
	for (ulong c = 0; c < 3; c++) 
	trace += mc[c][c];

	double dcorr = 1-trace/(frob(ma, 3,3)*frob(mb,3,3));*/
     // cos_arr[scale][pairs[scale]++] = cos;
     // sumcos[scale] += cos;
     // sumdcorr[scale] += dcorr;
     // pairs[scale]++;
     //double dcorr =0;

/* printf("%lu cross corr computed, starting bootstrapping\n", s);
  ulong size_alloc = 10000000;
  float *cos_temp ;// = calloc(size_alloc, sizeof(float));
  float *dist_temp;// = calloc(size_alloc, sizeof(float));
  ulong *pairs    ;// = calloc(nbins, sizeof(ulong));
  float *cos_sum  ;//= calloc(nbins, sizeof(float));
  float *cos_boot = calloc(100, sizeof(float));
  for(int k = 0; k < nvol; k++){
    cos_temp  = calloc(size_alloc, sizeof(float));
    dist_temp = calloc(size_alloc, sizeof(float));
    pairs     = calloc(nbins, sizeof(ulong));
    cos_sum   = calloc(nbins, sizeof(float));
    ulong nbub_k = 0;
    //cos_temp[k]  = malloc(nbub_bins[k] * sizeof(float));
    //dist_temp[k] = malloc(nbub_bins[k] * sizeof(float));
    for(ulong j = 0; j < s; j++){
      // printf("%lu \n", cur);
      int scale = find_scale(dist_bins[k], dist_arr[j], nbins);
      if(vol_arr[j] > -1 && vol_arr[j] <= k && scale != nbins -1){
	cos_temp[nbub_k] = cos_arr[j];
	dist_temp[nbub_k++] = dist_arr[j];
	cos_sum[scale] += cos_arr[j];
	pairs[scale]++;
      }
    }
    // for(ulong j = 0; j< nbub_k; j++)
    //  printf("DOne %f\n", cos_arr[j] );

    for (ulong i = 0; i < nbins; i++) {
      if (pairs[i] < 100 && i != nbins - 1){
	pairs[i+1]   += pairs[i];
	cos_sum[i+1] += cos_sum[i];
      } else if (pairs[i] > 0){
	fprintf(outfile[k], "%10.2f,%10lu,%10.3f,", dist_bins[k][i],pairs[i],cos_sum[i]/pairs[i]);
	#pragma omp parallel
	{
	  uint seed = 25234 + 17*omp_get_thread_num();
	  #pragma omp for 
	  for(ulong j = 0; j < 100; j++){
	    float sum = 0;
	    for(ulong t = 0; t  < pairs[i]; t++){
	      ulong id = rand_64(&seed) % nbub_k;
	      // printf("DOne %f\n", cos_temp[id] );
	      sum += cos_temp[id];
	    }
	    fprintf(outfile[k], "%10.3f,", sum / pairs[i]);
	    //   cos_boot[j] = sum / pairs[i];
	  }
	}
	//	for (ulong j = 0; j < 100; j++)
	//  fprintf(outfile[k], "%10.3f,", cos_boot[j]);
	//	fprintf(outfile[k], "\n");
      }
    }
    fclose(outfile[k]);
    free(cos_temp);
    free(dist_temp);
    free(pairs);
    free(cos_sum);
  }

  
    
  free(cos_arr);
  free(vol_arr);
}
*/

void orientbis_analysis(char *fname, char *fname_out, ulong volumeth, float elongth){
  ulong count;
  FILE *infile = fopen(fname, "r");
  FILE *outfile = fopen(fname_out, "w");
  if (infile==NULL)  exit(-1);
  if (outfile==NULL)  exit(-1);  

  fscanf(infile, "count %lu\n", &count);
  printf("Number of bubbles to analyse: %lu\n", count);
  //printf("Thresh %ld %lf, %lf\n", volumeth, elongth, ncompth);

 float *center = malloc(3*count*sizeof(float));
 ulong *volume = malloc(count*sizeof(ulong));
 float *matrix = malloc(9*count*sizeof(float));
 ulong *id = malloc(count*sizeof(ulong));
 memset(id, -1, sizeof(ulong)*count);
 double *vec_main = malloc(count*3*sizeof(double));

 float *elong = malloc(count*sizeof(float));
 double *dist = malloc(count*sizeof(double));
 for (ulong i = 0; i < count; i++)
   dist[i]=-1;

 double *cos = malloc(count*sizeof(double));

 for (ulong i = 0; i < count; i++) {
   ulong lol;
   fscanf(infile, "id %lu\n",&lol);
   fscanf(infile, "center %f %f %f\n", center+i*3,  center+i*3+1,center+i*3+2 );
   fscanf(infile, "volume %lu\n", volume+i);
   fscanf(infile, "matrix %f %f %f %f %f %f %f %f %f\n",matrix+i*9,  matrix+i*9+1,matrix+i*9+2, matrix+i*9+3, matrix+i*9+4,matrix+i*9+5, matrix+i*9+6, matrix+i*9+7,matrix+i*9+8);
 }

 fprintf(outfile, "idi,idj,volume,dist,elong,cos\n");
 #pragma omp parallel for
 for (ulong i = 0; i < count; i++) {
   double tens_mat[9] = {matrix[i*9],matrix[i*9+1],matrix[i*9+2], matrix[i*9+3], matrix[i*9+4],matrix[i*9+5], matrix[i*9+6], matrix[i*9+7],matrix[i*9+8]};
   double eigval[3] = {0};
   double eigvec[9] = {0};
  rs (3, tens_mat,eigval, 1, eigvec );
   // ncomp[i] = (matrix[i*9]+matrix[i*9+4]+matrix[i*9+8])/pow(volume[i], 5/3);
   elong[i] = eigval[2]/eigval[1];
   for(int j = 0; j <3; j++)
     vec_main[i*3+j] = eigvec[6+j];
   // fprintf(outfile, "%ld,%ld,%10.1lf,%10.1lf,%10.1lf,%10.1lf,%10.1lf,%10.1lf,%10.1lf\n", id[i], volume[i],elong,flat,spars,ncomp,ax_len[0],ax_len[1],ax_len[2]);
 }

 for (ulong i = 0; i < count; i++) {

   if(volume[i] < volumeth || elong[i] < elongth) continue;

   for (ulong j = 0; j < count; j++) {
     if(i == j || volume[j] < volumeth || elong[j] < elongth) continue;
     double ddist = sqrt(pow(center[3*i]-center[3*j],2)+pow(center[3*i+1]-center[3*j+1],2)+pow(center[3*i+2]-center[3*j+2],2));
     if((dist[i] != -1 && ddist > dist[i]) || id[j] == i) continue;
     dist[i] = ddist;
     id[i] = j;
     cos[i] = 2*pow(vec_main[i*3]*vec_main[j*3] + vec_main[i*3+1]*vec_main[j*3+1] + vec_main[i*3+2]*vec_main[j*3+2],2) - 1;
   }
 }

 for (ulong i = 0; i < count; i++) {
   if(dist[i] != -1)
     fprintf(outfile, "%ld,%ld,%ld,%10.0lf,%10.1lf,%10.3lf\n", i,id[i],volume[i],dist[i],elong[i],cos[i]);
 }
   fclose(infile);
   fclose(outfile);
   printf("FILES CLOSED\n");
}

void direction_analysis(char *fname, char *fname_out){
  ulong count;
  FILE *infile = fopen(fname, "r");
  FILE *outfile = fopen(fname_out, "w");
  if (infile==NULL)  exit(-1);
  if (outfile==NULL)  exit(-1);  

   fscanf(infile, "count %lu\n", &count);
  printf("Number of bubbles to analyse: %lu\n", count);

 double *center = malloc(3*count*sizeof(double));
 ulong *volume = malloc(count*sizeof(ulong));
 double *matrix = malloc(9*count*sizeof(double));
 ulong *id = malloc(count*sizeof(ulong));
 double *vec_main = malloc(count*3*sizeof(double));
 
 for (ulong i = 0; i < count; i++) {
   fscanf(infile, "id %lu\n", id+i);
   fscanf(infile, "center %lf %lf %lf\n", center+i*3, center+i*3+1,center+i*3+2);
   fscanf(infile, "volume %lu\n", &volume[i]);
   fscanf(infile, "matrix %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",matrix+i*9,  matrix+i*9+1,matrix+i*9+2, matrix+i*9+3, matrix+i*9+4,matrix+i*9+5, matrix+i*9+6, matrix+i*9+7,matrix+i*9+8);
   double tens_mat[9] = {matrix[i*9],matrix[i*9+1],matrix[i*9+2], matrix[i*9+3], matrix[i*9+4],matrix[i*9+5], matrix[i*9+6], matrix[i*9+7],matrix[i*9+8]};
   double eigval[3] = {0};
   double eigvec[9] = {0};
   rs (3, tens_mat,eigval, 1, eigvec );
   vec_main[i*3] = eigvec[6];
   vec_main[i*3+1] = eigvec[7];
   vec_main[i*3+2] = eigvec[8];
 }
 
 fprintf(outfile, "id,ida,dista,cosa,idb,distb,cosb,idc,distc,cosc\n");
 #pragma omp parallel for
 for (ulong i = 0; i < count; i++) {
   double dist_cu[3] = {FLT_MAX,FLT_MAX,FLT_MAX};
   double cos_cu[3] = {-2,-2,-2};
   ulong id_cu[3] = {-1,-1,-1};

   for(ulong j = 0; j < count; j++) {
     if(j == i) continue;
     double dir[3] = {center[3*i]-center[3*j], center[3*i+1]-center[3*j+1], center[3*i+2]-center[3*j+2]};
     double dist = sqrt(pow(dir[0],2)+pow(dir[1],2)+pow(dir[2],2));
     dir[0]/= dist;
     dir[1]/= dist;
     dir[2]/= dist;
     if(dist < dist_cu[0]){
       cos_cu[0] = 2*pow(vec_main[i*3]*dir[0] + vec_main[i*3+1]*dir[1] + vec_main[i*3+2]*dir[2],2) - 1;
       dist_cu[0] = dist;
       id_cu[0] = j;
     }
     if(dist < dist_cu[1] && volume[j] >= 100){
       cos_cu[1] = 2*pow(vec_main[i*3]*dir[0] + vec_main[i*3+1]*dir[1] + vec_main[i*3+2]*dir[2],2) - 1;
       dist_cu[1] = dist;
       id_cu[1] = j;
     }
     if(dist < dist_cu[2] && volume[j] >= 1000){
       cos_cu[2] = 2*pow(vec_main[i*3]*dir[0] + vec_main[i*3+1]*dir[1] + vec_main[i*3+2]*dir[2],2) - 1;
       dist_cu[2] = dist;
       id_cu[2] = j;
     }
   }
   fprintf(outfile, "%ld,%ld,%10.3lf,%10.3lf,%ld,%10.3lf,%10.3lf,%ld,%10.3lf,%10.3lf\n", id[i],id_cu[0],dist_cu[0],cos_cu[0],id_cu[1],dist_cu[1],cos_cu[1],id_cu[2],dist_cu[2],cos_cu[2]);
 }
 fclose(infile);
 fclose(outfile);
 printf("FILES CLOSED\n");
}


void direction_analysisbis(char *fname, char *fname_out, ulong threshold){
  ulong count;
  FILE *infile = fopen(fname, "r");
  FILE *outfile = fopen(fname_out, "w");
  if (infile==NULL)  exit(-1);
  if (outfile==NULL)  exit(-1);  

  fscanf(infile, "count %lu\n", &count);
  printf("Number of bubbles to analyse: %lu\n", count);

 double *center = malloc(3*count*sizeof(double));
 ulong *volume = malloc(count*sizeof(ulong));
 double *matrix = malloc(9*count*sizeof(double));
 ulong *id = malloc(count*sizeof(ulong));
 double *vec_main = malloc(count*3*sizeof(double));
 
 for (ulong i = 0; i < count; i++) {
   fscanf(infile, "id %lu\n", id+i);
   fscanf(infile, "center %lf %lf %lf\n", center+i*3, center+i*3+1,center+i*3+2);
   fscanf(infile, "volume %lu\n", &volume[i]);
   fscanf(infile, "matrix %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",matrix+i*9,  matrix+i*9+1,matrix+i*9+2, matrix+i*9+3, matrix+i*9+4,matrix+i*9+5, matrix+i*9+6, matrix+i*9+7,matrix+i*9+8);
   double tens_mat[9] = {matrix[i*9],matrix[i*9+1],matrix[i*9+2], matrix[i*9+3], matrix[i*9+4],matrix[i*9+5], matrix[i*9+6], matrix[i*9+7],matrix[i*9+8]};
   double eigval[3] = {0};
   double eigvec[9] = {0};
   rs (3, tens_mat,eigval, 1, eigvec );
   vec_main[i*3] = eigvec[6];
   vec_main[i*3+1] = eigvec[7];
   vec_main[i*3+2] = eigvec[8];
 }
 
 fprintf(outfile, "id,ida,idb,distmean,cos\n");
 #pragma omp parallel for
 for (ulong i = 0; i < count; i++) {
   double dist_cu[2] = {FLT_MAX,FLT_MAX};
   ulong id_cu[2] = {-1,-1};
   if(volume[i] < threshold) continue;
   for(ulong j = 0; j < count; j++) {
     if(j == i || volume[j] < threshold) continue;
     double dir[3] = {center[3*i]-center[3*j], center[3*i+1]-center[3*j+1], center[3*i+2]-center[3*j+2]};
     double dist = sqrt(pow(dir[0],2)+pow(dir[1],2)+pow(dir[2],2));
     dir[0]/= dist;
     dir[1]/= dist;
     dir[2]/= dist;
     if(dist < dist_cu[0]){
       dist_cu[1] = dist_cu[0];
       dist_cu[0] = dist;
       id_cu[1] = id_cu[0];
       id_cu[0] = j;
     } else if(dist < dist_cu[1]){
       dist_cu[1] = dist;
       id_cu[1] = j;
     }
   }
   double dira[3] = {center[3*i]-center[3*id_cu[0]], center[3*i+1]-center[3*id_cu[0]+1], center[3*i+2]-center[3*id_cu[0]+2]};
   double dista = sqrt(pow(dira[0],2)+pow(dira[1],2)+pow(dira[2],2));
   dira[0]/= dista;
   dira[1]/= dista;
   dira[2]/= dista;
   double dirb[3] = {center[3*i]-center[3*id_cu[1]], center[3*i+1]-center[3*id_cu[1]+1], center[3*i+2]-center[3*id_cu[1]+2]};
   double distb = sqrt(pow(dirb[0],2)+pow(dirb[1],2)+pow(dirb[2],2));
   dirb[0]/= distb;
   dirb[1]/= distb;
   dirb[2]/= distb;
   double dist_mean = (dista+distb)/2;
   double cos = 2*pow(dira[0]*dirb[0] + dira[1]*dirb[1] + dira[2]*dirb[2],2) - 1;
 
   fprintf(outfile, "%ld,%ld,%ld,%10.3lf,%10.3lf\n", id[i],id_cu[0],id_cu[1], dist_mean, cos);
 }
 fclose(infile);
 fclose(outfile);
 printf("FILES CLOSED\n");
}

const char *get_filename_ext(const char *filename) {
  const char *dot = strrchr(filename, '.');
  if(!dot || dot == filename) return "";
  return dot + 1;
}




void test_or(char *fname_inA,  ulong volth){
  ulong count;
  FILE *fp1 = fopen(fname_inA, "r");
  fscanf(fp1, "count %lu\n", &count);
  printf("Number of bubbles to analyse: %ld\n", count);

 ulong *volume = malloc(count*sizeof(ulong));
 double *vec_main = malloc(count*3*sizeof(double));

 // double *cos = malloc(count*sizeof(double));

 for (ulong i = 0; i < count; i++) {
   fscanf(fp1, "volume %lu\n", volume+i);
   fscanf(fp1, "vector %lf %lf %lf\n",vec_main+i*3,  vec_main+i*3+1, vec_main+i*3+2 );
 }
 // #pragma omp parallel for

 for (int j = 0; j< 3; j++){
   double cos = 0;
   ulong num = 0;
   for (ulong i = 0; i < count; i++) {
     if(volume[i] > volth && (volume[i] < volth*10 && volume[i] != 8000000)){
       double temp =2*pow(vec_main[3*i+j],2) - 1;
       cos+=temp;
       num++;
     }
   }
   printf("Number corr: %ld \n COS %lf\n\n",num, cos/num);
 }

 fclose(fp1);
}




char *remove_ext(char* mystr) {
    char *retstr;
    char *lastdot;
    if (mystr == NULL)
         return NULL;
    if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
        return NULL;
    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}

double frob(double matrix[3][3], int size1, int size2)
{
    double result = 0.0;
    for(int i = 0; i < size1; ++i)
    {
        for(int j = 0; j < size2; ++j)
        {
	  //double val = *(matrix + (i*size2) + j);
	  double val = matrix[i][j];
            result += val * val;
        }
    }
    // printf("%lf \n", sqrt(result));
    if(result == 0) return 10000000;
    else
      return sqrt(result);
}
