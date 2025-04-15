#include "types.h"
#include "attributes.h"
#include "lambdavec.h"
#include "tree_filt.h"
#include "tree_flood.h"
#include "eor_analysis.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "eispack.h"
struct tms 	tstruct;				
clock_t 	start;

ulong long_rand(ulong size){
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  struct timeval tv; // Seed generation based on time
  gettimeofday(&tv,0);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;
  T = gsl_rng_default; // Generator setup
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, mySeed);
  ulong u =  gsl_rng_uniform_int(r, size); // Generate it!
  gsl_rng_free (r);
  return u;
}

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*   Analysis  Moment inertia      */
/*				   */
/* +++++++++++++++++++++++++++++++ */

int find_scale_misc(double *bins, double attribute, int nbins){
  int upper = nbins-1, lower = 0, mid;

  if (attribute >=  bins[upper]){
    //  warn("One value if higher than the highest bin");
    return upper;
  }

  mid = (upper + lower) / 2;
  while (mid!=lower) {
    if(attribute >= bins[mid])
      lower = mid;
    else
      upper = mid;

    mid = (upper + lower) / 2;
  }
  return lower;
} /* find_scale */



ulong nobj = 1;
void write_inertia_json(Node *tree, char *filename, LambdaVec *lvec, double *(*attribute)(void *), double lambda, ulong size){
  
  FILE *fp = malloc(1 * sizeof(FILE*));
  char *fname;
  
  asprintf(&fname, "%s_Inertia.json",filename);
  fp = fopen(fname, "wb");
  fprintf(fp, "{\"bubbles\":\n[");
  
  double *inertia;
  ulong  count = 0;
  bool   comma = 0;
  
  for (ulong v = 0; v < size; v++) { // To change with border
    if(tree->parent[v] != BOTTOM && is_levelroot(tree, v)) {
      if(!tree->attribute[v]) error("No attributes ? ");
      inertia = (*attribute)(tree->attribute[v]);
      if(inertia[0] < lambda) continue;

      double c_x = inertia[1] / inertia[0];
      double c_y = inertia[2] / inertia[0];
      double c_z = inertia[3] / inertia[0];
      double ixx = inertia[4] - inertia[1]*inertia[1]/inertia[0] + inertia[0]/12;
      double iyy = inertia[5] - inertia[2]*inertia[2]/inertia[0] + inertia[0]/12;
      double izz = inertia[6] - inertia[3]*inertia[3]/inertia[0] + inertia[0]/12;
      double ixy = inertia[7] - inertia[1]*inertia[2]/inertia[0];
      double iyz = inertia[8] - inertia[2]*inertia[3]/inertia[0];
      double ixz = inertia[9] - inertia[1]*inertia[3]/inertia[0];

      if(comma)
	fprintf(fp, ",\n");
      else
	comma= true;

      fprintf(fp, "{\"id\": %ld, \n", count);
      fprintf(fp, "\"center\": [%3.2lf, %3.2lf, %3.2lf], \n", c_x, c_y, c_z);
      fprintf(fp, "\"volume\": %.0lf,\n", inertia[0]);
      fprintf(fp, "\"t_matrix\":  [[%10.5lf, %10.5lf, %10.5lf], [%10.5lf, %10.5lf, %10.5lf], [%10.5lf, %10.5lf, %10.5lf]]}", ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz);
      count++;

    }
  }
  
  fprintf(fp, "],\n \"count\": %lu}", count);
  fclose(fp);
}



void write_inertia_txt(Node *tree, char *filename, double *(*attribute)(void *), double lambda, ulong size, ulong *dims, int attrib){
  FILE *fpine = malloc(1 * sizeof(FILE*));
  char *fnameine;
  asprintf(&fnameine, "%s_Inertia.txt", filename);
  fpine = fopen(fnameine, "wb");
  if (fpine==NULL)  error("File not opened");  

  FILE *fppos = malloc(1 * sizeof(FILE*));
  char *fnamepos;
  asprintf(&fnamepos, "%s_Position.txt", filename);
  fppos = fopen(fnamepos, "wb");
  if (fppos==NULL)  error("File not opened");  
  //fuse_periodic_3D(tree, tree->store, dims, 26);

  double *inertia;
  ulong  count = 0;
  // bool   comma = 0;

  // #pragma omp parallel for
  for (ulong v = 0; v < size; v++) { // To change with border
    if(tree->parent[v] != BOTTOM && is_levelroot(tree, v)){
      inertia = (*attribute)(tree->attribute[v]);
      if(inertia[0] < lambda) continue;
      count++;
    }
  }

  info("Number of bubbles detected: %lu", count);
  fprintf(fpine, "count %lu\n", count);
  fprintf(fppos, "%lu\n", count);
  fprintf(fppos, "0 %lu 0 %lu 0 %lu\n", dims[0], dims[1], dims[2]);

  srand(time(NULL));
  count = 0;

  for (ulong v = 0; v < size; v++) { // To change with border
    if(tree->parent[v] != BOTTOM && is_levelroot(tree, v)){

      inertia = (*attribute)(tree->attribute[v]);
      if(inertia[0] < lambda) continue;

      double c_x = inertia[1]/inertia[0];
      double c_y = inertia[2]/inertia[0];
      double c_z = inertia[3]/inertia[0];
      double ixx = inertia[4] - inertia[1]*inertia[1]/inertia[0] + inertia[0]/12;
      double iyy = inertia[5] - inertia[2]*inertia[2]/inertia[0] + inertia[0]/12;
      double izz = inertia[6] - inertia[3]*inertia[3]/inertia[0] + inertia[0]/12;
      double ixy = inertia[7] - inertia[1]*inertia[2]/inertia[0];
      double iyz = inertia[8] - inertia[2]*inertia[3]/inertia[0];
      double ixz = inertia[9] - inertia[1]*inertia[3]/inertia[0];

      //  fprintf(fpine, "id %lu \n", count);
      //   fprintf(fpine, "center %3.2lf %3.2lf %3.2lf \n", c_x < 0 ? c_x + dims[0]: c_x, c_y <0 ? c_y + dims[1]:c_y, c_z<0 ? c_z +dims[2]:c_z);
      //    fprintf(fppos, "%3.0lf %3.0lf %3.0lf %3.2lf\n", c_x < 0 ? c_x + dims[0]: c_x, c_y <0 ? c_y + dims[1]:c_y, c_z<0 ? c_z +dims[2]:c_z , pow(inertia[0], 1.0/3.0));
      //  fprintf(fpine, "volume %.0lf \n", inertia[0]);
      fprintf(fpine, "matrix %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf\n", ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz);
     
      count++;

    }
  }
  fclose(fpine);
  fclose(fppos);
}


void write_inertia_bin(Node *tree, char *filename, double *(*attribute)(void *), double lambda, ulong size, ulong *dims, ulong *attr_off, int attrib){

  
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  struct timeval tv; // Seed generation based on time
  gettimeofday(&tv,0);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;
  T = gsl_rng_default; // Generator setup
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, mySeed);
  // ulong u =  gsl_rng_uniform_int(r, size); // Generate it!

  /* Opening files, declaring variables */
  
  FILE *fpine = malloc(1 * sizeof(FILE*));
  char *fnameine;
  asprintf(&fnameine, "%s_ine.bin", filename);
  fpine = fopen(fnameine, "wb");
  if (fpine==NULL)  error("Inertia file not opened");  

  FILE *fpvol = malloc(1 * sizeof(FILE*));
  char *fnamevol;
  asprintf(&fnamevol, "%s_vol.bin", filename);
  fpvol = fopen(fnamevol, "wb");
  if (fpvol==NULL)  error("Volume file not opened");

  FILE *fpint = malloc(1 * sizeof(FILE*));
  char *fnameint;
  asprintf(&fnameint, "%s_int.bin", filename);
  fpint = fopen(fnameint, "wb");
  if (fpint==NULL)  error("Intensity file not opened");

  FILE *fpcurv = malloc(1 * sizeof(FILE*));
  char *fnamecurv;
  asprintf(&fnamecurv, "%s_curv.bin", filename);
  fpcurv = fopen(fnamecurv, "wb");
  if (fpcurv==NULL)  error("Intensity file not opened");

  FILE *fpcen = malloc(1 * sizeof(FILE*));
  char *fnamecen;
  asprintf(&fnamecen, "%s_cen.bin", filename);
  fpcen = fopen(fnamecen, "wb");
  if (fpcen==NULL)  error("Centroids file not opened");
  
  FILE *fppos = malloc(1 * sizeof(FILE*));
  char *fnamepos;
  asprintf(&fnamepos, "%s_Position.txt", filename);
  fppos = fopen(fnamepos, "wb");
  if (fppos==NULL)  error("Position file not opened");  

  FILE *fpdens;
  char *fnamedens;

  FILE *fpdensine;
  char *fnamedensine;

  FILE *fpdenscen;
  char *fnamedenscen;

  FILE *fpdensr;
  char *fnamedensr;

  FILE *fpdensiner;
  char *fnamedensiner;

  FILE *fpdenscenr;
  char *fnamedenscenr;
  
  double *inertia;
  void   **attr_dens, **attr_dens_rand;
  AuxDataStore *store_dens, *store_dens_rand;  

  if(tree->dens){
    //  sum_dens   = calloc(size,sizeof(double));
    //   sum_dens_2 = calloc(size,sizeof(double));
    attr_dens  = calloc(size,sizeof(void*));
    store_dens = malloc(sizeof(AuxDataStore)); check_alloc(store_dens, 0);
    init_aux_data_store(store_dens, AttribsArray[attrib].size, size);
    
    // sum_dens_rand   = calloc(size,sizeof(double));
    //   sum_dens_rand_2 = calloc(size,sizeof(double));
    attr_dens_rand  = calloc(size,sizeof(void*));
    store_dens_rand = malloc(sizeof(AuxDataStore)); check_alloc(store_dens_rand, 0);
    init_aux_data_store(store_dens_rand, AttribsArray[attrib].size, size);
    
    fpdens = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedens, "%s_dens.bin", filename);

    fpdens = fopen(fnamedens, "wb");
    if (fpdens==NULL)  error("Density average file not opened");

    fpdensine = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedensine, "%s_dens_ine.bin", filename);
    fpdensine = fopen(fnamedensine, "wb");
    if (fpdensine==NULL)  error("Density inertia file not opened");

    fpdenscen = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedenscen, "%s_dens_cen.bin", filename);
    fpdenscen = fopen(fnamedenscen, "wb");
    if (fpdenscen==NULL)  error("Density centroids file not opened");

    fpdensr = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedensr, "%s_dens_rand.bin", filename);
    fpdensr = fopen(fnamedensr, "wb");
    if (fpdensr==NULL)  error("Density r average file not opened");

    fpdensiner = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedensiner, "%s_dens_rand_ine.bin", filename);
    fpdensiner = fopen(fnamedensiner, "wb");
    if (fpdensiner==NULL)  error("Density r inertia file not opened");

    fpdenscenr = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedenscenr, "%s_dens_rand_cen.bin", filename);
    fpdenscenr = fopen(fnamedenscenr, "wb");
    if (fpdenscenr==NULL)  error("Density r centroids file not opened");
    
    info("Finding the number of bubbles and the average density in the cube");
  } else
    info("Finding the number of bubbles in the cube");

  double *init_attr    = calloc(4, sizeof(double));
  double mean_dens_tot = 0;
  double mean_dens_out = 0;
  ulong  count = 0;
  
  for (ulong v = 0; v < size; v++) { // To change with MPI
    ulong v_lr = get_levelroot(tree, v);
    
    if(tree->dens){
      mean_dens_tot  += (double) tree->gval_dens[v];
      init_attr[0]    = (double) ((v % (dims[0] * dims[1])) % dims[0]  + attr_off[0]);
      init_attr[1]    = (double) ((v % (dims[0] * dims[1])) / dims[0]  + attr_off[1]);
      init_attr[2]    = (double) (v / (dims[0] * dims[1]) + attr_off[2]);
      init_attr[3]    = (double) tree->gval_dens[v];
      if( attr_dens[v_lr])
	add_to_aux_data(attr_dens[v_lr], init_attr);
      else
	attr_dens[v_lr] = new_aux_data(store_dens, init_attr);

      ulong v_rand =  gsl_rng_uniform_int(r, size);
      init_attr[3]         = (double) tree->gval_dens[v_rand] ;
      if(attr_dens_rand[v_lr])
	add_to_aux_data(attr_dens_rand[v_lr], init_attr);
      else
	attr_dens_rand[v_lr] =	new_aux_data(store_dens_rand, init_attr);

    }
    
    if(tree->parent[v] != BOTTOM && is_levelroot(tree, v)){
      inertia = (*attribute)(tree->attribute[v]);
      if(inertia[0] < lambda) continue;
      count++;   
    }
    
  }

  info("Number of individual bubbles: %lu", count);
  fprintf(fppos, "%lu\n", count);
  fprintf(fppos, "0 %lu 0 %lu 0 %lu\n", dims[0], dims[1], dims[2]);
  
  if(tree->dens){
    mean_dens_tot /= (double) size;
    info("Average density in the cube is %lf", mean_dens_tot);
  }
  
  info("Writing the results in the file");

  for (ulong v = 0; v < size; v++) { // To change with border

    if(is_levelroot(tree, v)){
      if(tree->parent[v] != BOTTOM ){
      
	if(!tree->attribute[v]) error("One node is not defined ?");
	
	inertia = (*attribute)(tree->attribute[v]);
	if(inertia[0] < lambda) continue;
	
	ulong volume     = (ulong) inertia[0];
	double intensities[2] =   {inertia[10], inertia[11]};
	double c_xyz[3]  = { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
	double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[10] + inertia[10]/12,
			     inertia[5] - inertia[2]*inertia[2]/inertia[10] + inertia[10]/12,
			     inertia[6] - inertia[3]*inertia[3]/inertia[10] + inertia[10]/12,
			     inertia[7] - inertia[1]*inertia[2]/inertia[10],
			     inertia[8] - inertia[2]*inertia[3]/inertia[10],
			     inertia[9] - inertia[1]*inertia[3]/inertia[10] };

	fwrite(&volume, 1 * sizeof(ulong),  1, fpvol);
	fwrite(intensities, 2 * sizeof(double),  1, fpint);
	fwrite(inertia + 12, 1 * sizeof(double),  1, fpcurv);

	fwrite(c_xyz,   3 * sizeof(double), 1, fpcen);
	fwrite(matrix,  6 * sizeof(double), 1, fpine);	
	fprintf(fppos, "%3.0lf %3.0lf %3.0lf %3.2lf\n", c_xyz[0], c_xyz[1], c_xyz[2], pow(inertia[0], 1.0/3.0));

	if(tree->dens){

	  inertia = (*attribute)(attr_dens[v]);
	  if((ulong) inertia[0] != volume) warn("IMPOSSIBLE ?");
	  
	  // double temp_mean  = inertia[10]/inertia[0];
	  //  double temp_std   = (inertia[11]-inertia[10]*inertia[10]/inertia[0])/inertia[0];
	  double intensitiesr[2] =   {inertia[10], inertia[11]};

	  double c_xyzd[3]  =  { inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
	  double matrixd[6] =  { inertia[4] - inertia[1]*inertia[1]/inertia[10] + inertia[10]/12,
				 inertia[5] - inertia[2]*inertia[2]/inertia[10] + inertia[10]/12,
				 inertia[6] - inertia[3]*inertia[3]/inertia[10] + inertia[10]/12,
				 inertia[7] - inertia[1]*inertia[2]/inertia[10],
				 inertia[8] - inertia[2]*inertia[3]/inertia[10],
				 inertia[9] - inertia[1]*inertia[3]/inertia[10] };
	  fwrite(intensitiesr,   2 * sizeof(double), 1, fpdens);
	  // fwrite(&temp_std, 1 * sizeof(double), 1, fpdens);
	  fwrite(c_xyzd,  3 * sizeof(double), 1, fpdenscen);      
	  fwrite(matrixd, 6 * sizeof(double), 1, fpdensine);

	  inertia = (*attribute)(attr_dens_rand[v]);
	  if((ulong) inertia[0] != volume) warn("IMPOSSIBLE ?");
	  
	  // temp_mean   	     =  inertia[10]/inertia[0];
	  // temp_std   	     =  (inertia[11]-inertia[10]*inertia[10]/inertia[0])/inertia[0];
	  double intensitiesdr[2] =   {inertia[10], inertia[11]};
	  double c_xyzdr[3]  =  {inertia[1]/inertia[10],inertia[2]/inertia[10],inertia[3]/inertia[10]} ;
	  double matrixdr[6] =  { inertia[4] - inertia[1]*inertia[1]/inertia[10] + inertia[10]/12,
				 inertia[5] - inertia[2]*inertia[2]/inertia[10] + inertia[10]/12,
				 inertia[6] - inertia[3]*inertia[3]/inertia[10] + inertia[10]/12,
				 inertia[7] - inertia[1]*inertia[2]/inertia[10],
				 inertia[8] - inertia[2]*inertia[3]/inertia[10],
				 inertia[9] - inertia[1]*inertia[3]/inertia[10] };
	  fwrite(intensitiesdr,    2 * sizeof(double), 1, fpdensr);
	  //fwrite(&temp_std,    1 * sizeof(double), 1, fpdensr);	

	  fwrite(c_xyzdr,  3 * sizeof(double), 1, fpdenscenr);      
	  fwrite(matrixdr, 6 * sizeof(double), 1, fpdensiner);
	}
	
      }

    }

  }

  gsl_rng_free (r);

  /*...*/
  fclose(fpine);
  fclose(fpvol);
  fclose(fpcen);
  fclose(fppos);
  fclose(fpint);
  fclose(fpcurv);

  if(tree->dens) {
    //   free(sum_dens);
    //free(sum_dens_rand);
    //free(sum_dens_2);
    //free(sum_dens_rand_2);;
    fclose(fpdens);
    fclose(fpdensine);
    fclose(fpdenscen);
    fclose(fpdensr);
    fclose(fpdensiner);
    fclose(fpdenscenr);
  }
}


void write_inertiaall_bin(Arguments *args, Node *tree, char *filename, double *(*attribute)(void *), double lambda, ulong size, ulong *dims, ulong *attr_off, int attrib, double *extrema){

  
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  struct timeval tv; // Seed generation based on time
  gettimeofday(&tv,0);
  unsigned long mySeed = tv.tv_sec + tv.tv_usec;
  T = gsl_rng_default; // Generator setup
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, mySeed);
  
  FILE *fpbins = malloc(1 * sizeof(FILE*));
  char *fname;
  asprintf(&fname, "%s_bins.bin", filename);
  fpbins = fopen(fname, "wb");
  if (fpbins==NULL)  error("Bins file not opened");
  int nbins = 101;
  double bin = (extrema[1]-extrema[0])/(double) nbins;
  
  double *int_bins = calloc(nbins, sizeof(double));
  for (int i = 0; i < nbins; i++)
    int_bins[i] = extrema[0] + i*bin;
  fwrite(int_bins, nbins * sizeof(double),  1, fpbins);
  fclose(fpbins);
  info("First bin is [%lf, %lf)", int_bins[0], int_bins[1]);

  /* Opening files, declaring variables */
  
  FILE *fpine = malloc(1 * sizeof(FILE*));
  char *fnameine;
  asprintf(&fnameine, "%s_ine.bin", filename);
  fpine = fopen(fnameine, "wb");
  if (fpine==NULL)  error("Inertia file not opened");  

  FILE *fpvol = malloc(1 * sizeof(FILE*));
  char *fnamevol;
  asprintf(&fnamevol, "%s_vol.bin", filename);
  fpvol = fopen(fnamevol, "wb");
  if (fpvol==NULL)  error("Volume file not opened");

  FILE *fpint = malloc(1 * sizeof(FILE*));
  char *fnameint;
  asprintf(&fnameint, "%s_int.bin", filename);
  fpint = fopen(fnameint, "wb");
  if (fpint==NULL)  error("Intensity file not opened");

  FILE *fpcurv = malloc(1 * sizeof(FILE*));
  char *fnamecurv;
  asprintf(&fnamecurv, "%s_curv.bin", filename);
  fpcurv = fopen(fnamecurv, "wb");
  if (fpcurv==NULL)  error("Intensity file not opened");

  FILE *fpcen = malloc(1 * sizeof(FILE*));
  char *fnamecen;
  asprintf(&fnamecen, "%s_cen.bin", filename);
  fpcen = fopen(fnamecen, "wb");
  if (fpcen==NULL)  error("Centroids file not opened");
  
  FILE *fppos = malloc(1 * sizeof(FILE*));
  char *fnamepos;
  asprintf(&fnamepos, "%s_Position.txt", filename);
  fppos = fopen(fnamepos, "wb");
  if (fppos==NULL)  error("Position file not opened");  

  FILE *fpdens;
  char *fnamedens;

  FILE *fpdensine;
  char *fnamedensine;

  FILE *fpdenscen;
  char *fnamedenscen;

  FILE *fpdensr;
  char *fnamedensr;

  FILE *fpdensiner;
  char *fnamedensiner;

  FILE *fpdenscenr;
  char *fnamedenscenr;
  
  double *inertia; 

  if(tree->dens){

    fpdens = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedens, "%s_dens.bin", filename);

    fpdens = fopen(fnamedens, "wb");
    if (fpdens==NULL)  error("Density average file not opened");

    fpdensine = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedensine, "%s_dens_ine.bin", filename);
    fpdensine = fopen(fnamedensine, "wb");
    if (fpdensine==NULL)  error("Density inertia file not opened");

    fpdenscen = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedenscen, "%s_dens_cen.bin", filename);
    fpdenscen = fopen(fnamedenscen, "wb");
    if (fpdenscen==NULL)  error("Density centroids file not opened");

    fpdensr = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedensr, "%s_dens_rand.bin", filename);
    fpdensr = fopen(fnamedensr, "wb");
    if (fpdensr==NULL)  error("Density r average file not opened");

    fpdensiner = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedensiner, "%s_dens_rand_ine.bin", filename);
    fpdensiner = fopen(fnamedensiner, "wb");
    if (fpdensiner==NULL)  error("Density r inertia file not opened");

    fpdenscenr = malloc(1 * sizeof(FILE*));
    asprintf(&fnamedenscenr, "%s_dens_rand_cen.bin", filename);
    fpdenscenr = fopen(fnamedenscenr, "wb");
    if (fpdenscenr==NULL)  error("Density r centroids file not opened");
    
    info("Finding the number of bubbles and the average density in the cube");
  } else
    info("Finding the number of bubbles in the cube");

  ulong  count = 0;

  bool *visited = calloc(size, sizeof(bool));
  idx  *list = calloc(size, sizeof(idx));
  
  for (ulong v = 0; v < size; v++) { 
    ulong v_lr = get_levelroot(tree, v);
   
    if(is_levelroot(tree, v) && !visited[v]){

      idx p = v;

      int scale = find_scale_misc(int_bins, (double) tree->gval[p], nbins);

      idx parent = get_parent(tree, p);
      if(parent != BOTTOM ){
	int newscale = find_scale_misc(int_bins, (double) tree->gval[parent], nbins);
	while(newscale == scale && get_parent(tree, p) != BOTTOM && !visited[p]){
	  visited[p] = true;
	  p = parent;
	  parent = get_parent(tree, p);
	  newscale = find_scale_misc(int_bins, (double) tree->gval[parent], nbins);
	}
      }
      inertia = (*attribute)(tree->attribute[p]);

      if(inertia[0] < lambda || visited[p]) continue;

      list[count++] = p;
      visited[p] = true;
    }
  }

  info("Number of individual bubbles: %lu", count);
  fprintf(fppos, "%lu\n", count);
  fprintf(fppos, "0 %lu 0 %lu 0 %lu\n", dims[0], dims[1], dims[2]);
  
  /*if(tree->dens){
    mean_dens_tot /= (double) size;
    info("Average density in the cube is %lf", mean_dens_tot);
    }*/
  
  info("Writing the results");
  for (ulong p = 0; p < count; p++) { // To change with border
    ulong v = list[p];

    if(!tree->attribute[v]) error("One node is not defined ?");
    

    inertia = (*attribute)(tree->attribute[v]);
    //	if(inertia[0] < lambda) continue;
    ulong volume     = (ulong) inertia[0];
    double intensities[2] =   {inertia[10], inertia[11]};
    double c_xyz[3]  = { inertia[1]/inertia[0],inertia[2]/inertia[0],inertia[3]/inertia[0]} ;
    double matrix[6] = { inertia[4] - inertia[1]*inertia[1]/inertia[0] + inertia[0]/12,
			 inertia[5] - inertia[2]*inertia[2]/inertia[0] + inertia[0]/12,
			 inertia[6] - inertia[3]*inertia[3]/inertia[0] + inertia[0]/12,
			 inertia[7] - inertia[1]*inertia[2]/inertia[0],
			 inertia[8] - inertia[2]*inertia[3]/inertia[0],
			 inertia[9] - inertia[1]*inertia[3]/inertia[0] };

    double intensity = (double) tree->gval[v];
    fwrite(&volume, 1 * sizeof(ulong),  1, fpvol);
    fwrite(intensities, 2 * sizeof(double),  1, fpint);
    fwrite(&intensity, 1 * sizeof(double),  1, fpcurv);
    fwrite(c_xyz,   3 * sizeof(double), 1, fpcen);
    fwrite(matrix,  6 * sizeof(double), 1, fpine);	
    fprintf(fppos, "%3.0lf %3.0lf %3.0lf %3.2lf\n", c_xyz[0], c_xyz[1], c_xyz[2], pow(inertia[0], 1.0/3.0));

    if(tree->dens){

      //  inertia = (*attribute)(attr_dens[v]);
      // if((ulong) inertia[0] != volume) warn("IMPOSSIBLE ?");
	  
      // double temp_mean  = inertia[10]/inertia[0];
      //  double temp_std   = (inertia[11]-inertia[10]*inertia[10]/inertia[0])/inertia[0];
      double intensitiesr[2] =   {inertia[21], inertia[22]};
      double c_xyzd[3]  =  { inertia[12]/inertia[21],inertia[13]/inertia[21],inertia[14]/inertia[21]} ;
      double matrixd[6] =  { inertia[15] - inertia[12]*inertia[12]/inertia[21] + inertia[21]/12,
			     inertia[16] - inertia[13]*inertia[13]/inertia[21] + inertia[21]/12,
			     inertia[17] - inertia[14]*inertia[14]/inertia[21] + inertia[21]/12,
			     inertia[18] - inertia[12]*inertia[13]/inertia[21],
			     inertia[19] - inertia[13]*inertia[14]/inertia[21],
			     inertia[20] - inertia[12]*inertia[14]/inertia[21] };
      fwrite(intensitiesr,   2 * sizeof(double), 1, fpdens);
      fwrite(c_xyzd,  3 * sizeof(double), 1, fpdenscen);      
      fwrite(matrixd, 6 * sizeof(double), 1, fpdensine);
    }
	
  }
   

  gsl_rng_free (r);

  /*...*/
  fclose(fpine);
  fclose(fpvol);
  fclose(fpcen);
  fclose(fppos);
  fclose(fpint);
  fclose(fpcurv);

  if(tree->dens) {
    fclose(fpdens);
    fclose(fpdensine);
    fclose(fpdenscen);
    fclose(fpdensr);
    fclose(fpdensiner);
    fclose(fpdenscenr);
  }
}
 

void inertia_attributes_all(Arguments *args, Node *tree,  double *(*attribute)(void *), ulong *dims, ulong *attr_off, double *extrema){

  double lambda =  args->lambda_arg;

  MPI_File  fpine, fpvol, fpelo, fpfla, fpspa, fpnco, fpcurv;
  char *filename = args->outprefix_arg;
  char *fname;
  double lims[6] = {attr_off[0]+tree->border[0],
		    attr_off[0]+dims[0]-tree->border[1],
		    attr_off[1]+tree->border[2],
		    attr_off[1]+dims[1]-tree->border[3],
		    attr_off[2]+tree->border[4],
		    attr_off[2]+dims[2]-tree->border[5]};

    FILE *fpbins = malloc(1 * sizeof(FILE*));
    if(rank() ==0){
      asprintf(&fname, "%s_bins.bin", filename);
      fpbins = fopen(fname, "wb");
      if (fpbins==NULL)  error("Bins file not opened");
    } 

  
  asprintf(&fname, "%s_ine.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpine);


  asprintf(&fname, "%s_vol.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpvol);
 
  asprintf(&fname, "%s_elong.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpelo);

  asprintf(&fname, "%s_flat.bin", filename);
    if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpfla);

  asprintf(&fname, "%s_spars.bin", filename);
    if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpspa);

  
  asprintf(&fname, "%s_ncomp.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpnco);

  
  asprintf(&fname, "%s_curv.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpcurv);

  
  
  
  int nbins = args->nlev_arg+2;
  double bin = (extrema[1]-extrema[0])/(double) args->nlev_arg;
  double *int_bins = calloc(nbins, sizeof(double));
  for (int i = 0; i < nbins; i++)
    int_bins[i] = extrema[0] - bin/2.0 + i*bin;
  
  if(rank() ==0){
    fwrite(int_bins, nbins * sizeof(double),  1, fpbins);
    fclose(fpbins);
  }

  info("Number of bins is %d, first bin is [%lf, %lf), last bin is [%lf, %lf)",nbins, int_bins[0], int_bins[1],int_bins[nbins-2], int_bins[nbins-1] );


  bool *visited = calloc(tree->size, sizeof(bool));
  ulong   count_tot = 0;
  ulong   count_loc = 0;
  ulong   *count    = calloc(args->threads_arg, sizeof(ulong));
  ulong   *volume;
  double  *matrix;
  double  *intensity;
  double  *elong;
  double  *ncomp;
  double  *spars;
  double  *flat;

  info("Finding the number of bubbles in the cube");


  
  #pragma omp parallel 
  {
    int nthreads = omp_get_num_threads();
    int id	   = omp_get_thread_num();
    ulong lwb 	   = id*tree->size/nthreads;
    ulong upb	   = (id+1)*tree->size/nthreads;
    idx  *list     = calloc(upb-lwb, sizeof(idx));
    double *inertia;

    for (ulong v = lwb; v < upb; v++) { // To change with MPI
      if(is_levelroot(tree, v) && !visited[v]){
	ulong p = v;
	int scale = find_scale_misc(int_bins, (double) tree->gval[p], nbins);
	idx parent = get_parent(tree, p);
	if(parent != BOTTOM ){
	  int newscale = find_scale_misc(int_bins, (double) tree->gval[parent], nbins);
	  while(newscale == scale && get_parent(tree, p) != BOTTOM && !visited[p] ){
	    if (p >=lwb && p < upb) visited[p] = true;
	    p = parent;
	    parent = get_parent(tree, p);
	    newscale = find_scale_misc(int_bins, (double) tree->gval[parent], nbins);
	  }
	}
	if(tree->attribute[p] == NULL || *(*attribute)(tree->attribute[p])< lambda || visited[p] || (p <lwb && p >= upb) ) continue;

	inertia = (*attribute)(tree->attribute[p]);
	double cen[3]  =  {inertia[1]/inertia[0],inertia[2]/inertia[0],inertia[3]/inertia[0]} ;
	if( cen[0] >= lims[0] && cen[0] < lims[1] && cen[1] >= lims[2] && cen[1] < lims[3] && cen[2] >= lims[4] && cen[2] < lims[5]) {
	  list[count[id]++] = p;
	  visited[p] = true;
	}
      }
    }

    ulong copycount = count[id];  
    #pragma omp barrier

    if(id == 0){
      for(int i = 1; i < nthreads; i++)
	count[i] += count[i-1];
      count_loc = count[nthreads-1];
      free(visited);
      MPI_Reduce(&count_loc, &count_tot, 1, MPI_UINT64_T, MPI_SUM, 0,  MPI_COMM_WORLD);
      info("Number of individual bubbles among all nodes: %lu (node 0: %lu)", count_tot,count_loc );

      volume = calloc( count_loc, sizeof(ulong));
      matrix = calloc( count_loc*6, sizeof(double));
      intensity = calloc( count_loc, sizeof(double));
      elong = calloc( count_loc, sizeof(double));
      ncomp = calloc( count_loc, sizeof(double));
      spars = calloc( count_loc, sizeof(double));
      flat  = calloc( count_loc, sizeof(double));

    }

    #pragma omp barrier

    ulong offset = id == 0 ? 0 :count[id-1];
    for (ulong p = 0; p <  copycount; p++) { // To change with border
      ulong v = list[p];
      if(!tree->attribute[v]) info("One node (v %ld) is not defined ?", v);
    
      inertia = (*attribute)(tree->attribute[v]);
      volume[p+offset]     = (ulong) inertia[0];
      matrix[p*6+offset]   = inertia[4] - inertia[1]*inertia[1]/inertia[0] + inertia[0]/12;
      matrix[p*6+1+offset] = inertia[5] - inertia[2]*inertia[2]/inertia[0] + inertia[0]/12;
      matrix[p*6+2+offset] = inertia[6] - inertia[3]*inertia[3]/inertia[0] + inertia[0]/12;
      matrix[p*6+3+offset] = inertia[7] - inertia[1]*inertia[2]/inertia[0];
      matrix[p*6+4+offset] = inertia[8] - inertia[2]*inertia[3]/inertia[0];
      matrix[p*6+5+offset] = inertia[9] - inertia[1]*inertia[3]/inertia[0];

      
      double eigval[3];
      double tens_mat[9] = {matrix[p*6+offset], matrix[p*6+3+offset], matrix[p*6+5+offset], matrix[p*6+3+offset], matrix[p*6+1+offset],matrix[p*6+4+offset], matrix[p*6+5+offset], matrix[p*6+4+offset],matrix[p*6+2+offset]};
      rs (3, tens_mat, eigval, 0, NULL);
      double ax_len[3]        = {sqrt(20*eigval[2]/volume[p+offset]), sqrt(20*eigval[1]/volume[p+offset]), sqrt(20*eigval[0]/volume[p+offset])};
	
      intensity[p+offset]  = (double) tree->gval[v];
      ncomp[p+offset]      = (matrix[p*6+0+offset]+matrix[p*6+1+offset]+matrix[p*6+2+offset])/pow(volume[p+offset], 5.0/3.0);
      elong[p+offset]      = eigval[2]/eigval[1];
      flat[p+offset]       = eigval[1]/eigval[0];
      spars[p+offset]      =  3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*volume[p+offset]);

    }
  }


  MPI_File_write_ordered(fpvol, volume,  count_loc,  MPI_UINT64_T, NULL);
  MPI_File_write_ordered(fpelo, elong,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpfla, flat,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpcurv, intensity,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpspa, spars,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpnco, ncomp,  count_loc,  MPI_DOUBLE, NULL);  
  MPI_File_close(&fpelo);
  MPI_File_close(&fpvol);
  MPI_File_close(&fpfla);
  MPI_File_close(&fpcurv);
  MPI_File_close(&fpspa);
  MPI_File_close(&fpnco);


}




void inertia_attributes_bins(Arguments *args, Node *tree,  double *(*attribute)(void *), ulong *dims, ulong *attr_off, double *extrema){

  double lambda =  args->lambda_arg;

  MPI_File  fpine, fpvol, fpelo, fpfla, fpspa, fpnco, fpcurv, fpsca;
  char *filename = args->outprefix_arg;
  char *fname;
  double lims[6] = {attr_off[0]+tree->border[0],
		    attr_off[0]+dims[0]-tree->border[1],
		    attr_off[1]+tree->border[2],
		    attr_off[1]+dims[1]-tree->border[3],
		    attr_off[2]+tree->border[4],
		    attr_off[2]+dims[2]-tree->border[5]};

  FILE *fpbins = malloc(1 * sizeof(FILE*));
  if(rank() ==0){
    asprintf(&fname, "%s_bins.bin", filename);
    fpbins = fopen(fname, "wb");
    if (fpbins==NULL)  error("Bins file not opened");
  } 

  
  asprintf(&fname, "%s_ine.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpine);


  asprintf(&fname, "%s_vol.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpvol);
 
  asprintf(&fname, "%s_elong.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpelo);

  asprintf(&fname, "%s_flat.bin", filename);
    if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpfla);

  asprintf(&fname, "%s_spars.bin", filename);
    if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpspa);

  
  asprintf(&fname, "%s_ncomp.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpnco);

  
  asprintf(&fname, "%s_curv.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpcurv);

  asprintf(&fname, "%s_sca.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpsca);

  
  int nbins = args->nlev_arg+1;
  double bin = (extrema[1]-extrema[0])/(double) args->nlev_arg;
  double *int_bins = calloc(nbins, sizeof(double));
  for (int i = 0; i < nbins; i++)
    int_bins[i] = extrema[0] + i*bin;
  
  if(rank() ==0){
    fwrite(int_bins, nbins * sizeof(double),  1, fpbins);
    fclose(fpbins);
  }

  info("Number of bins is %d, first bin is [%lf, %lf), last bin is [%lf, %lf)",nbins, int_bins[0], int_bins[1],int_bins[nbins-2], int_bins[nbins-1] );


  ulong   count_tot = 0;
  ulong   count_loc = 0;
  bool    *visited = calloc(tree->size, sizeof(bool));
  ulong   *count = calloc(args->threads_arg, sizeof(ulong));
  ulong   *volume;
  double  *matrix;
  double  *intensity;
  double  *elong;
  double  *ncomp;
  double  *spars;
  double  *flat;
  int  *scale_int;
  info("Finding the number of bubbles in the cube");


  
  #pragma omp parallel 
  {
    int nthreads   = omp_get_num_threads();
    int id	   = omp_get_thread_num();
    ulong lwb 	   = id*tree->size/nthreads;
    ulong upb	   = (id+1)*tree->size/nthreads;
    idx  *list    = calloc(upb-lwb, sizeof(idx));
    double *inertia;

    for (ulong v = lwb; v < upb; v++) { // To change with MPI
      if(is_levelroot(tree, v) && !visited[v]){
	ulong p = v;
	int scale = find_scale_misc(int_bins, (double) tree->gval[p], nbins);
	idx parent = get_parent(tree, p);
	if(parent != BOTTOM ){
	  int newscale = find_scale_misc(int_bins, (double) tree->gval[parent], nbins);
	  while(newscale == scale && get_parent(tree, p) != BOTTOM && !visited[p] ){
	    if (p >=lwb && p < upb) visited[p] = true;
	    p = parent;
	    parent = get_parent(tree, p);
	    newscale = find_scale_misc(int_bins, (double) tree->gval[parent], nbins);
	  }
	}
	if(tree->attribute[p] == NULL || *(*attribute)(tree->attribute[p])< lambda || visited[p] || (p <lwb && p >= upb) ) continue;

	inertia = (*attribute)(tree->attribute[p]);
	double cen[3]  =  {inertia[1]/inertia[0],inertia[2]/inertia[0],inertia[3]/inertia[0]} ;
	if( cen[0] >= lims[0] && cen[0] < lims[1] && cen[1] >= lims[2] && cen[1] < lims[3] && cen[2] >= lims[4] && cen[2] < lims[5]) {
	  list[count[id]++] = p;
	  visited[p] = true;
	}
      }
    }

    ulong copycount = count[id];  
    #pragma omp barrier

    if(id == 0){
      for(int i = 1; i < nthreads; i++)
	count[i] += count[i-1];
      count_loc = count[nthreads-1];
      free(visited);
      MPI_Reduce(&count_loc, &count_tot, 1, MPI_UINT64_T, MPI_SUM, 0,  MPI_COMM_WORLD);
      info("Number of individual bubbles among all nodes: %lu (node 0: %lu)", count_tot,count_loc );

      volume = calloc( count_loc, sizeof(ulong));
      matrix = calloc( count_loc*6, sizeof(double));
      intensity = calloc( count_loc, sizeof(double));
      elong = calloc( count_loc, sizeof(double));
      ncomp = calloc( count_loc, sizeof(double));
      spars = calloc( count_loc, sizeof(double));
      flat  = calloc( count_loc, sizeof(double));
      scale_int  = calloc( count_loc*2, sizeof(int));

    }

    #pragma omp barrier

    ulong offset = id == 0 ? 0 :count[id-1];
    for (ulong p = 0; p <  copycount; p++) { // To change with border
      ulong v = list[p];
      if(!tree->attribute[v]) info("One node (v %ld) is not defined ?", v);

      int scale = find_scale_misc(int_bins, (double) tree->gval[v], nbins);
      idx parent = get_parent(tree, v);
      int newscale = parent != -1 ? find_scale_misc(int_bins, (double) tree->gval[parent], nbins): nbins;

      inertia = (*attribute)(tree->attribute[v]);
      volume[p+offset]     = (ulong) inertia[0];
      matrix[p*6+offset]   = inertia[4] - inertia[1]*inertia[1]/inertia[0] + inertia[0]/12;
      matrix[p*6+1+offset] = inertia[5] - inertia[2]*inertia[2]/inertia[0] + inertia[0]/12;
      matrix[p*6+2+offset] = inertia[6] - inertia[3]*inertia[3]/inertia[0] + inertia[0]/12;
      matrix[p*6+3+offset] = inertia[7] - inertia[1]*inertia[2]/inertia[0];
      matrix[p*6+4+offset] = inertia[8] - inertia[2]*inertia[3]/inertia[0];
      matrix[p*6+5+offset] = inertia[9] - inertia[1]*inertia[3]/inertia[0];

      
      double eigval[3];
      double tens_mat[9] = {matrix[p*6+offset], matrix[p*6+3+offset], matrix[p*6+5+offset], matrix[p*6+3+offset], matrix[p*6+1+offset],matrix[p*6+4+offset], matrix[p*6+5+offset], matrix[p*6+4+offset],matrix[p*6+2+offset]};
      rs (3, tens_mat, eigval, 0, NULL);
      double ax_len[3]        = {sqrt(20*eigval[2]/volume[p+offset]), sqrt(20*eigval[1]/volume[p+offset]), sqrt(20*eigval[0]/volume[p+offset])};
	
      intensity[p+offset]  = (double) tree->gval[v];
      ncomp[p+offset]      = (matrix[p*6+0+offset]+matrix[p*6+1+offset]+matrix[p*6+2+offset])/pow(volume[p+offset], 5.0/3.0);
      elong[p+offset]      = eigval[2]/eigval[1];
      flat[p+offset]       = eigval[1]/eigval[0];
      spars[p+offset]      =  3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*volume[p+offset]);
      scale_int[(p+offset)*2] = scale;
      scale_int[(p+offset)*2+1] = newscale;
    }
  }


  MPI_File_write_ordered(fpvol, volume,  count_loc,  MPI_UINT64_T, NULL);
  MPI_File_write_ordered(fpelo, elong,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpfla, flat,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpcurv, intensity,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpspa, spars,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpnco, ncomp,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpsca, scale_int,  count_loc*2,  MPI_INT, NULL);  

  MPI_File_close(&fpelo);
  MPI_File_close(&fpvol);
  MPI_File_close(&fpfla);
  MPI_File_close(&fpcurv);
  MPI_File_close(&fpspa);
  MPI_File_close(&fpnco);
  MPI_File_close(&fpsca);


}



void inertia_attributes_hii(Arguments *args, Node *tree,  double *(*attribute)(void *), ulong *dims){

  double lambda =  args->lambda_arg;

  MPI_File  fpine, fpvol, fpelo, fpfla, fpspa, fpnco, fpcurv, fpsca;
  char *filename = args->outprefix_arg;
  char *fname;


  
  asprintf(&fname, "%s_ine.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpine);


  asprintf(&fname, "%s_vol.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpvol);
 
  asprintf(&fname, "%s_elong.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpelo);

  asprintf(&fname, "%s_flat.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpfla);

  asprintf(&fname, "%s_spars.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpspa);

  
  asprintf(&fname, "%s_ncomp.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpnco);

  
  asprintf(&fname, "%s_curv.bin", filename);
  if(rank() == 0)
    MPI_File_delete(fname,MPI_INFO_NULL);
  err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE|MPI_MODE_EXCL|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpcurv);
  


  ulong   count_tot = 0;
  ulong   count_loc = 0;
  ulong   *count = calloc(args->threads_arg, sizeof(ulong));
  ulong   *volume;
  double  *matrix;
  double  *intensity;
  double  *elong;
  double  *ncomp;
  double  *spars;
  double  *flat;
  int  *scale_int;
  
  info("Finding the number of bubbles in the cube");


  
  #pragma omp parallel 
  {
    int nthreads   = omp_get_num_threads();
    int id	   = omp_get_thread_num();
    ulong lwb 	   = id*tree->size/nthreads;
    ulong upb	   = (id+1)*tree->size/nthreads;
    idx  *list    = calloc(upb-lwb, sizeof(idx));
    double *inertia;

    for (ulong v = lwb; v < upb; v++) { // To change with MPI
      //info("%f", tree->gval[v]);
      if(is_levelroot(tree, v) && tree->gval[v] == 1.){
	ulong p = v;
	if(tree->attribute[v] == NULL || *(*attribute)(tree->attribute[v])< lambda) continue;
	inertia = (*attribute)(tree->attribute[v]);
	double cen[3]  =  {inertia[1]/inertia[0],inertia[2]/inertia[0],inertia[3]/inertia[0]} ;
	list[count[id]++] = v;
      }
    }
    

    ulong copycount = count[id];  
    #pragma omp barrier

    if(id == 0){
      for(int i = 1; i < nthreads; i++)
	count[i] += count[i-1];
      count_loc = count[nthreads-1];
      MPI_Reduce(&count_loc, &count_tot, 1, MPI_UINT64_T, MPI_SUM, 0,  MPI_COMM_WORLD);
      info("Number of individual bubbles among all nodes: %lu (node 0: %lu)", count_tot,count_loc );

      volume = calloc( count_loc, sizeof(ulong));
      matrix = calloc( count_loc*6, sizeof(double));
      elong = calloc( count_loc, sizeof(double));
      ncomp = calloc( count_loc, sizeof(double));
      spars = calloc( count_loc, sizeof(double));
      flat  = calloc( count_loc, sizeof(double));
    }

    #pragma omp barrier

    ulong offset = id == 0 ? 0 :count[id-1];
    for (ulong p = 0; p <  copycount; p++) { // To change with border
      ulong v = list[p];
      if(!tree->attribute[v]) info("One node (v %ld) is not defined ?", v);

      inertia = (*attribute)(tree->attribute[v]);
      volume[p+offset]     = (ulong) inertia[0];
      matrix[p*6+offset]   = inertia[4] - inertia[1]*inertia[1]/inertia[0] + inertia[0]/12;
      matrix[p*6+1+offset] = inertia[5] - inertia[2]*inertia[2]/inertia[0] + inertia[0]/12;
      matrix[p*6+2+offset] = inertia[6] - inertia[3]*inertia[3]/inertia[0] + inertia[0]/12;
      matrix[p*6+3+offset] = inertia[7] - inertia[1]*inertia[2]/inertia[0];
      matrix[p*6+4+offset] = inertia[8] - inertia[2]*inertia[3]/inertia[0];
      matrix[p*6+5+offset] = inertia[9] - inertia[1]*inertia[3]/inertia[0];

      
      double eigval[3];
      double tens_mat[9] = {matrix[p*6+offset], matrix[p*6+3+offset], matrix[p*6+5+offset], matrix[p*6+3+offset], matrix[p*6+1+offset],matrix[p*6+4+offset], matrix[p*6+5+offset], matrix[p*6+4+offset],matrix[p*6+2+offset]};
      rs (3, tens_mat, eigval, 0, NULL);
      double ax_len[3]        = {sqrt(20*eigval[2]/volume[p+offset]), sqrt(20*eigval[1]/volume[p+offset]), sqrt(20*eigval[0]/volume[p+offset])};

      ncomp[p+offset]      = (matrix[p*6+0+offset]+matrix[p*6+1+offset]+matrix[p*6+2+offset])/pow(volume[p+offset], 5.0/3.0);
      elong[p+offset]      = eigval[2]/eigval[1];
      flat[p+offset]       = eigval[1]/eigval[0];
      spars[p+offset]      =  3.14*ax_len[0]*ax_len[1]*ax_len[2]/(6*volume[p+offset]);

    }
  }


  MPI_File_write_ordered(fpvol, volume,  count_loc,  MPI_UINT64_T, NULL);
  MPI_File_write_ordered(fpelo, elong,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpfla, flat,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpspa, spars,  count_loc,  MPI_DOUBLE, NULL);
  MPI_File_write_ordered(fpnco, ncomp,  count_loc,  MPI_DOUBLE, NULL);

  MPI_File_close(&fpelo);
  MPI_File_close(&fpvol);
  MPI_File_close(&fpfla);
  MPI_File_close(&fpspa);
  MPI_File_close(&fpnco);

}
 
/* +++++++++++++++++++++++++++++++ */
/*				   */
/*         Annexes Functions       */
/*				   */
/* +++++++++++++++++++++++++++++++ */


double *hessian(value *gvals, ulong *dims, ulong p){
  // info("%ld ", p);

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

void tree_seg_dir(Node *tree, value *out, bool *reached, ulong *rank, ulong lwb, ulong upb, double max_var, double (*attribute)(void *), double lambda) {
  value  val;
  idx	 parent;
  ulong   u, w;
  srand(time(NULL));
  for (long v = upb-1; v >= 0; v-- ) {
    idx p = rank[v];
    //  info("v %ld, gval %f", p, tree->gval[p]);

    if (!reached[p]) { /* not filtered yet */
      // info("IN REACHED v %ld, gval %f", p, tree->gval[p]);

      w = get_levelroot(tree, p);
      parent = get_parent(tree, w);

      
      //  info("Start %ld, %f, %f", w, tree->gval[w], tree->gval[parent]);
      /* repeat while we're not at the bottom, th */
      while ((parent != BOTTOM) && (!reached[w]) && ((tree->gval[w]-tree->gval[parent])/tree->gval[w] < max_var) ) {
        w = parent;
        parent = tree->parent[w];	      
      }
      //   if(parent != BOTTOM){
      //	out[parent] = 0;
      //	reached[parent] = true;
      while (parent != BOTTOM && !reached[parent]) {
	out[parent] = 0;
	reached[parent] = true;
	parent = tree->parent[parent];
	//  info("Start %ld, %f, %f", w, tree->gval[w], tree->gval[parent]);
      }
      // if((tree->gval[w]-tree->gval[parent])/tree->gval[w]> 0.0)
      if (reached[w]) {
        /* criterion satisfied at level tree[w].filter */
	val = out[w];
      }
      else if (tree->attribute[w] && (*attribute)(tree->attribute[w]) >= lambda)
        /* criterion cannot be satisfied */
        val = rand()%255;//nobj++;
      else
	val = 0;

      //if((*attribute)(tree->attribute[w]) < lambda && val == 1) info("%lf, %lf", val);
		
      // info("AFTER REACHED");
    /* set filt along par-path from v to w */
    u = p;
    while (u != w) {
      if ((lwb <= u) && (upb > u)){
	out[u] 	= val;
	reached[u] 	= true;
      }
      u = tree->parent[u];
    }
    if ( (lwb <= w) && (upb > w)){
      out[w] 		= val;
      reached[w]      = true;
    }
    
    }
  }
} /* tree_filter_direct */


