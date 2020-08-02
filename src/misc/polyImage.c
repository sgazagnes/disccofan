#include "common.h"
#include "misc.h"
#include "kmeans.h"
#include "io.h"
#include "eispack.h"
#include <fftw3.h>
#include <complex.h>


double *im_out, *im_in, *im_inb;
#define PI (double) (3.14159265358979323846264338327)

int main(int argc, char** argv) {
  if(argc <= 1){
    printf("Miscelaneous functions:  %s <transformation> : {create, mod_type, mod_depth, print_info, print_all, compare, bin_to_fits, slice_plot, binary, lofar} \n", argv[0]);
    exit(0);
  } else {
    
    /**************************************************************************************/

    /******************************  Generate images  *************************************/

    /**************************************************************************************/    
      
    if(!strcmp(argv[1], "create")){
      if (argc <= 9 || (atoi(argv[2]) < 0 || atoi(argv[2]) > 4) ) {
	printf("Usage <create>: generate images with a bit-per-pixel between -32 and 32 bits.\n Arguments: <type> : 0 for random, 1 for linear sequence, 2 for random binary ionization field \n <output name without type suffix> \n <output type suffix> \n <width> \n <height> \n <depth> \n <bit depth> \n <negative values>: 1 to generate negative values, 0 otherwise \n [ionization fraction] \n [periodic]");
	exit(0);
      }
            
      ulong dims[3]      = {atoi(argv[5]),atoi(argv[6]),atoi(argv[7])};
      int   bitpix       = atoi(argv[8]);
      int   negative     = atoi(argv[9]);
      int   n_dims       = dims[2] > 1 ? 3 : 2;
      ulong size         = dims[0] * dims[1] * dims[2];
      
      if(atoi(argv[2]) == 0)       im_out = create_random(size, bitpix, negative);
      else if(atoi(argv[2]) == 1)  im_out = create_sequence(size, bitpix, negative);
      else if(atoi(argv[2]) == 2)  im_out = create_ioniz(size, atof(argv[10]), dims);
      else if(atoi(argv[2]) == 3)  im_out = create_bubbles_ov(dims[0], size, atof(argv[10]),atof(argv[11]));
      else if(atoi(argv[2]) == 4)  im_out = create_bubbles_nov(dims[0], size, atof(argv[10]), atof(argv[11]));

      write_file(argv[3], argv[4], "data", dims, im_out, bitpix, n_dims);
    }




    /**************************************************************************************/

    /*****************************  Change image type *************************************/

    /**************************************************************************************/ 

    else if (!strcmp(argv[1], "mod_type")){
       
      if (argc<=4 || (argc<=5 && !strcmp(argv[3], "h5"))){
	printf("Usage <mod_type>: modify data format.\n %s <input without type suffix> <input type suffixe> <output type suffix> [h5 input dataset if hdf5 file]\n", argv[1]);
	exit(0);
      }
      
      int     bitpix;
      ulong   dims[3] = {1,1,1};

      im_in = read_file(argv[2], argv[3], argv[5], dims, &bitpix);

      int n_dims = dims[2] > 1 ? 3 : 2;

      write_file(argv[2], argv[4], "data", dims, im_in, bitpix, n_dims);

    } 




    
    /**************************************************************************************/

    /*************************  Change image dynamic range  *******************************/

    /**************************************************************************************/ 

    else if (!strcmp(argv[1], "mod_depth")){
      if (argc<=6 || (argc<=7 && !strcmp(argv[3], "h5"))) {
	printf("Usage <mod_depth>: Modify the dynamic range of an image.\n Arguments: <input without type prefix> \n <input type suffixe> \n <output name without type suffix> \n <output type suffix> \n <new bit depth> \n  [dataset if hdf5 file] \n");
	exit(0);
      }

      ulong dims[3] = {1,1,1};
      int   bitpix;      
      im_in = read_file(argv[2], argv[3], argv[7], dims, &bitpix);
      int   n_dims = dims[2] > 1 ? 3 : 2;
      int   newbitpix = atoi(argv[6]);
      
      im_out = modify_bpp(im_in, dims[0]*dims[1]*dims[2], bitpix, newbitpix);      
      bitpix = newbitpix;
      write_file(argv[4], argv[5], "data", dims, im_out, bitpix, n_dims);

    }



    
    /**************************************************************************************/

    /***************************  Print image info/values  ********************************/

    /**************************************************************************************/ 

    else if( !strcmp(argv[1], "print_info") || !strcmp(argv[1], "print_all") ){
      if (argc<=3 || (argc<=4 && !strcmp(argv[3], "h5"))) {
	printf("Function <%s>: print image parameters. Arguments: <input name without type suffix> \n <input prefix> \n [h5 dataset name]", argv[1]);
	exit(0);
      }

      ulong dims[3] = {1,1,1};
      int   bitpix;
      im_in = read_file(argv[2], argv[3], argv[4], dims, &bitpix);
      int   n_dims = dims[2] > 1 ? 3 : 2;

      double max_gval_in = find_max(im_in, dims[0]*dims[1]*dims[2]);
      double min_gval_in = find_min(im_in, dims[0]*dims[1]*dims[2]);
      
      printf("Image paramaters: Width %lu \n Height %lu \n Depth %lu \n bit_per_pixel %d\n [min, max] values [%lf, %lf] \n", dims[0], dims[1], dims[2], bitpix,min_gval_in, max_gval_in);

      if(!strcmp(argv[1], "print_all")){
	printf("Pixel intensities: \n");
	for(size_t k = 0; k < dims[2]; k++){
	  for(size_t i = 0; i < dims[1]; i++){
	    for(size_t j = 0; j < dims[0]; j++){
	      printf("%5.2lf   ", im_in[k*(dims[0]*dims[1]) + i*dims[0] + j]);
	    }
	    printf("\n");
	  }
	  printf("\n\n");
	}
      }
    }


    
    /**************************************************************************************/

    /*********************************  Compare image   ***********************************/

    /**************************************************************************************/ 

  
    else if(!strcmp(argv[1], "compare") ){
      if (argc<=6) {
	printf("Function <compare>: Compare images or txt files.\n Arguments: <input name of image A without suffix> \n  <input type suffix image A> \n <input name of image B without suffix> \n <input type suffix image A> \n <number of images> \n [dataset if image A or B is hdf5] \n [dataset if image B and A is hdf5] \n");
	exit(0);
      }

      if(!strcmp(argv[3], "txt")){
	char *filename_in, *filename_out;
	asprintf(&filename_in, "%s.%s", argv[2], argv[3]);
	FILE *fp1 = fopen(filename_in, "r");
	asprintf(&filename_in, "%s.%s", argv[4], argv[5]);
	FILE *fp2 = fopen(filename_in, "r");
	if (fp1 == NULL || fp2 == NULL)
	  {
	    printf("Error : Files not open\n");
	    exit(0);
	  }
 
	int error = compare_txt_files(fp1, fp2);
	fclose(fp1);
	fclose(fp2);
      } else{
    
	int   n_images    = atoi(argv[6]);
	ulong dims_A[3] = {1,1,1};
	ulong dims_B[3] = {1,1,1};
	int   bitpix_A, bitpix_B;
	float diff = 0, maxdiff = 0;
	char *prefix;
	for (int i = 0; i < n_images; i++) {
	  n_images > 1 ? asprintf(&prefix, "%s-%d", argv[2],i) :
	    asprintf(&prefix, "%s", argv[2]); 
	  im_in = read_file(prefix, argv[3], argv[7], dims_A, &bitpix_A);	  		
	  int n_dimsA = dims_A[2] > 1 ? 3 : 2;
	  n_images > 1 ? asprintf(&prefix, "%s-%d", argv[4],i) :
	    asprintf(&prefix, "%s", argv[4]); 
	  im_inb = read_file(prefix, argv[5], argv[8], dims_B, &bitpix_B);	  	      
	  int n_dimsB = dims_B[2] > 1 ? 3 : 2;	
	  if(n_dimsA != n_dimsB){
	    printf("Error: N dimensions not matching in A (%d) and B (%d)",n_dimsA, n_dimsB);
	    exit(0);
	  }
	  if(bitpix_A != bitpix_B)
	    printf("Warn: different bit depth in the A (%d) and B (%d)",bitpix_A, bitpix_B);
	  if(dims_A[0] != dims_B[0] || dims_A[1] != dims_B[1] || dims_A[2] != dims_B[2] ){
	    printf("Error: dimensions not matching in A (W:%lu, H:%lu, D:%lu) and B (W:%lu, H:%lu, D:%lu)", dims_A[0], dims_A[1], dims_A[2], dims_B[0], dims_B[1], dims_B[2]);
	    exit(0);
	  }

	  ulong size = dims_A[0]*dims_A[1]*dims_A[2];
	  im_out = malloc(size * sizeof(double));
	
	  for(size_t k = 0; k < size; k++){
	    im_out[k] = bitpix_A < 0 ? (double) fabs( im_in[k] -  im_inb[k]) : (double) abs( im_in[k] -  im_inb[k]);
	    if(im_out[k] > 0){
	      diff = 1;
	      maxdiff = im_out[k] > maxdiff ? im_out[k] : maxdiff;
	    }
	  }

	  n_images > 1 ? asprintf(&prefix, "Compare-%d", i) :
	    asprintf(&prefix, "Compare"); 
	  write_file(prefix, "fits", "data", dims_A, im_out, -32, n_dimsA);
	}
	if(diff)
	  printf("\x1B[91m Error: differences spotted (Max Diff is %f) \033[0m\n", maxdiff);
	else
	  printf("\x1B[92m No differences spotted \033[0m\n");
      }


    /**************************************************************************************/

    /***********************************  bintofits   *************************************/

    /**************************************************************************************/ 
      
    }else if( !strcmp(argv[1], "bintofits") ){
      if (argc<=6) {
	printf("Function <bin_to_fits>: convert binary files to fits. Arguments: <input binary file name> <output name without type suffix> <output type suffix>  <input bpp> <dimension of the sides>");
	exit(0);
      }

      ulong dims[3];
      if(argc == 6){
	dims[0] = dims[1] = dims[2] = atoi(argv[6]);
      }else{
	dims[0] = atoi(argv[6]);
	dims[1] = atoi(argv[7]);
	dims[2] = atoi(argv[8]);  
      }

      int bitpix = atoi(argv[5]);
      FILE *F = fopen(argv[2], "rb");
      ulong size = dims[0]*dims[1]*dims[2];
      
      im_out = calloc(size, sizeof(double));

      //mod_fread(gvals_out, size * sizeof(value), 1, F);
      for(ulong i = 0; i < dims[0]; i++){
	for(ulong j = 0; j < dims[1]; j++){
	  for(ulong k = 0; k < dims[2]; k++){
	    ulong index = k*dims[2]*dims[1]+(dims[0]-i-1)*dims[2]+j;
	    // ulong index = i*hdims[0]*hdims[1]+hdims[0]*j+k;
	    if(bitpix == -64) {
	      double val;
	      fread(&val, sizeof(double), 1, F);
	      im_out[index] =(double) val;

	    } else if (bitpix == -32){
	      float val;
	      fread(&val, sizeof(float), 1, F);
	      im_out[index] =(double) val;

	    }else if (bitpix == 8){
	      ubyte val;
	      fread(&val, sizeof(ubyte), 1, F);
	      im_out[index] =(double) val;

	    }else if (bitpix == 16){
	      ushort val;
	      fread(&val, sizeof(ushort), 1, F);
	      im_out[index] =(double) val;

	    }else if (bitpix == 32){
	      uint val;
	      fread(&val, sizeof(uint), 1, F);
	      im_out[index] =(double) val;

	    }else{
	      ulong val;
	      fread(&val, sizeof(ulong), 1, F);
	      im_out[index] =(double) val;

	    }  
	  }
	}
      }
      fclose(F);
      
      int n_dims = dims[2] > 1 ? 3 : 2;

      write_file(argv[3], "fits", "data", dims, im_out, bitpix, n_dims);
	
    }


    

    /**************************************************************************************/

    /***********************************  slice_plot   ************************************/

    /**************************************************************************************/ 

    else if( !strcmp(argv[1], "slice_plot") ){
      if (argc<=9 || (argc<=10 && !strcmp(argv[3], "h5"))) {
	printf("Function <%s>: cut a smaller tile in the image. Arguments: <input name without type suffix> \n <input prefix> \n <start x> \n <end x> \n <start y> \n <end y> \n <start z> \n <end z> \n [h5 dataset name]", argv[1]);
	exit(0);
      }
      
      char *prefix;
      ulong dims[3] = {1,1,1};
      ulong ndims[3] = {1,1,1};
      int   bitpix;
      im_in = read_file(argv[2], argv[3], argv[10], dims, &bitpix);

      ulong x_s = atoi(argv[4]);
      ulong x_e = atoi(argv[5]);
      ulong y_s = atoi(argv[6]);
      ulong y_e = atoi(argv[7]);
      ulong z_s = atoi(argv[8]);
      ulong z_e = atoi(argv[9]);
      ndims[0] = (x_e-x_s);
      ndims[1] = (y_e-y_s);
      ndims[2] = (z_e-z_s);      
      im_out = malloc((x_e-x_s)*(y_e-y_s)*(z_e-z_s)*sizeof(double));
      ulong p = 0;
      for(size_t k = z_s; k < z_e; k++)
	for(size_t i = y_s; i < y_e; i++)
	  for(size_t j = x_s; j < x_e; j++,p++)
	    im_out[p] = im_in[k*dims[0]*dims[1]+i*dims[0]+j];

      asprintf(&prefix, "%s_slice", argv[2]);
      int n_dims = ndims[2] > 1 ? 3 : 2;
      write_file(prefix, argv[3], "data", ndims, im_out, bitpix, n_dims);
    }

    
    /**************************************************************************************/

    /*****************************************  bw  ***************************************/

    /**************************************************************************************/

    
    else if( !strcmp(argv[1], "bw") ){
      if (argc<=6 || (argc<=7 && !strcmp(argv[3], "h5"))) {
	printf("Function <%s>: create a black and white data with thresholding method. Arguments:<input name without type suffix> \n <input suffix> \n  <Thresholding method: 0 use a threshold value , 2 kmeans \n <threshold value min> \n <threshold value max> \n   [h5 dataset name]", argv[1]);
	exit(0);
      }

      char *prefix;
      ulong dims[3] = {1,1,1};
      int   bitpix;
      im_in = read_file(argv[2], argv[3], argv[9], dims, &bitpix);
      ulong size = dims[0]*dims[1]*dims[2];
      int   n_dims = dims[2] > 1 ? 3 : 2;
      
      im_out = malloc(size*sizeof(double));
      if(atoi(argv[4]) == 0){
	for(ulong i = 0; i < size; i++)
	  im_out[i] = im_in[i] >= atof(argv[5]) &&  im_in[i] <= atof(argv[6]) ? 255 : 0;
      } else {
	double center[2] = {0, 10};
	int *labels = malloc(size * sizeof(int));
        double *gvals_cop = malloc(size*sizeof(double));
	for(ulong i =0; i< size; i++)
	  gvals_cop[i]= (double) im_in[i];
	kmeans(1, gvals_cop, size, 2, center, labels);
	for(ulong i = 0; i < size; i++)
	  im_out[i] = 255 - 255*labels[i];
	free(gvals_cop);
      }

      bitpix = 8;
      asprintf(&prefix, "%s_thr", argv[2]);
      write_file(prefix, "fits", "data", dims, im_out, bitpix, n_dims);
  }



    
    /**************************************************************************************/

    /**************************************  Lofar  ***************************************/

    /**************************************************************************************/
    
    else if( !strcmp(argv[1], "lofar") ){
      if (argc<=3) {
	printf("Function <%s>: combine lofar channels. < obs id min > < obs id max > < sky dim >", argv[1]);
	exit(0);
      }

      char *prefix;
      ulong dims[3] = {1,1,1};
      int   bitpix;
      int   interv = atoi(argv[3]) - atoi(argv[2]);
      ulong dims_out[3] = {atoi(argv[4]), atoi(argv[4]), interv*3};
      ulong size;
      im_out = malloc(dims_out[0]*dims_out[1]*dims_out[2]*sizeof(double));
      ulong slice = 0;
      for(int i =0; i < interv ; i++){
	for(int j = 0; j < 3; j++, slice++){
	  asprintf(&prefix, "L254871-SB%03d-UV50_250_natural-000%d-I-image_reweighted", atoi(argv[2]),j);
	  im_in = read_file(prefix, "fits", NULL, dims, &bitpix);
	  for(ulong ii = 0; ii < MAX(dims_out[1], dims[1]); ii++)
	    for(ulong jj = 0; jj < MAX(dims_out[0], dims[0]); jj++)
	      im_out[slice*dims_out[0]*dims_out[1] + ii*dims_out[0] +jj] = im_in[dims[0]*(ii+dims_out[0])+jj+dims_out[0]];
	  free(im_in);
	}
      }
      write_file("combined_lofar", "fits", "data", dims, im_out, bitpix, 3);
    }

        
    /**************************************************************************************/

    /******************************  Read density from DFTE  ******************************/

    /**************************************************************************************/

    
    else if( !strcmp(argv[1], "density") ){
      if (argc<=4) {
	printf("Function <%s>: Read density files, <input binary>, <output name (automatic fits)>, <dimensions> \n", argv[1]);
	exit(0);
      }
      ulong dims[3];

      if(argc == 5){
	dims[0] = dims[1] = dims[2] = atoi(argv[4]);
      }else{
	dims[0] = atoi(argv[4]);
	dims[1] = atoi(argv[5]);
	dims[2] = atoi(argv[6]);  
      }
      FILE *fp1 = fopen(argv[2], "r");
      int   bitpix = -32;    
      ulong size = dims[0]*dims[1]*dims[2];
      im_out     = malloc(dims[0]*dims[1]*dims[2]*sizeof(double));
      int n_dims = dims[2] > 1 ? 3 : 2;
      double x, y, z, val;
      for(ulong i = 0; i < size; i++){
	fscanf(fp1,"%lf %lf %lf %lf",&x, &y, &z, &val);
	im_out[((ulong)z)*dims[0]*dims[1]+((ulong)y)*dims[0]+((ulong)x)] = (double) val;
      }
      write_file(argv[3], "fits", "data", dims, im_out, bitpix, n_dims);
    }

    /**************************************************************************************/

    /**********************************  Read position   **********************************/

    /**************************************************************************************/

    else if( !strcmp(argv[1], "position") ){
      if (argc<=3) {
	printf("Function <%s>: Read position files , <input binary/txt>, <output name (automatic fits)> \n", argv[1]);
	exit(0);
      }

      FILE *fp1 = fopen(argv[2], "r");
      int   bitpix = -32;    
      ulong dims[3];
      double x, y, z;
      ulong nval;
      fscanf(fp1, "%ld\n", &nval);
      printf("%ld positions to read \n", nval);
      fscanf(fp1, "0 %ld 0 %ld 0 %ld \n", dims, dims+1, dims+2);
      printf("Full dimensions: %ld  %ld  %ld \n", dims[0],dims[1],dims[2]);
      ulong size = dims[0]*dims[1]*dims[2];
      im_out     = malloc(dims[0]*dims[1]*dims[2]*sizeof(double));
      int n_dims = dims[2] > 1 ? 3 : 2;
      double temp;
      for(ulong i = 0; i < nval; i++){
	fscanf(fp1,"%lf %lf %lf %lf \n",&x, &y, &z, &temp);
	im_out[((ulong)z)*dims[0]*dims[1]+((ulong)y)*dims[0]+((ulong)x)] = 255;
      }

      write_file(argv[3], "fits", "data", dims, im_out, bitpix, n_dims);
    }


    /**************************************************************************************/
    
    /*******************************  Cross-correlating   *********************************/

    /**************************************************************************************/
    
    else if( !strcmp(argv[1], "cross_corr") ){
      if (argc<=5) {
	printf("Function <%s>: Compute cross corr between two images: <input name of image A without suffix> \n  <input type suffix image A> \n <input name of image B without suffix> \n <input type suffix image A> \n", argv[1]);
	exit(0);
      }

      ulong dims_A[3] = {1,1,1};
      ulong dims_B[3] = {1,1,1};
      int bitpix_A, bitpix_B;
      im_in = read_file(argv[2], argv[3], argv[7], dims_A, &bitpix_A);	  		
      int n_dimsA = dims_A[2] > 1 ? 3 : 2;
      im_inb = read_file(argv[4], argv[5], argv[8], dims_B, &bitpix_B);	  		
      im_out   = cross_corr(im_in, im_inb, dims_A[0]*dims_A[1]*dims_A[2]);
      write_file("cross_corr", "fits", "data", dims_A, im_out, -32, n_dimsA);
    }
    
     /**************************************************************************************/
    
    /*******************************  Hessian   *********************************/

    /**************************************************************************************/
    
    else if( !strcmp(argv[1], "hessian") ){
      if (argc<=5) {
	printf("Function <%s>: Compute hessian and eigenvec at point x,y,z: <input name of image A without suffix> \n  <input type suffix image A> \n <x>, <y>, <z> \n", argv[1]);
	exit(0);
      }


      ulong dims[3] = {1,1,1};
      int bitpix;
      im_in = read_file(argv[2], argv[3], NULL, dims, &bitpix);	  		
      int n_dims = dims[2] > 1 ? 3 : 2;

      im_out     = calloc(dims[0]*dims[1]*dims[2],sizeof(double));
      for(ulong i =0; i< dims[0]*dims[1]*dims[2]; i++){
	ulong x = i % dims[0];
	ulong y = (i %(dims[0]*dims[1]))/dims[0];
	ulong z = i/(dims[0]*dims[1]);
	double *arr = hessian(im_in, dims, i);
	//	printf("Hessian at %ld, %ld, %ld: hxx =  %lf, hyy = %lf, hzz = %lf, hxy = %lf, hyz = %lf, hxz = %lf\n",x,y,z,arr[0],arr[1], arr[2], arr[3], arr[4], arr[5]);
	double tens_mat[9] = {arr[0], arr[3], arr[5], arr[3], arr[1],arr[4], arr[5], arr[4],arr[2]};
	//double tens_mat[4] = {arr[0], arr[3], arr[3], arr[1]};

	double eigval[3] = {0};
	double eigvec[9] = {0};
	rs (3, tens_mat,eigval, 1, eigvec );
	im_out[i] = eigval[2];
      }
      write_file("maineig", "fits", "data", dims, im_out, -32, n_dims);
    }

    /**************************************************************************************/
    
    /*******************************  Percentile   *********************************/
    
    /**************************************************************************************/

    else if( !strcmp(argv[1], "percentile") ){
      if (argc<=4) {
	printf("Function <%s>: Gets x percentile value: <input name of image A without suffix> \n  <input type suffix image A> \n <x>, <y>, <z> \n", argv[1]);
	exit(0);
      }

      ulong dims[3] = {1,1,1};
      int bitpix;
      im_in = read_file(argv[2], argv[3], NULL, dims, &bitpix);	  		
      int n_dims = dims[2] > 1 ? 3 : 2;
      ulong percentile = (ulong)  atoi(argv[4]) * dims[0]*dims[1]*dims[2] / 100;
      double inmax = find_max(im_in, dims[0]*dims[1]*dims[2]);
      double inmin = find_min(im_in,  dims[0]*dims[1]*dims[2]);
      int nbins = 100;
      ulong histos[200] = {0};
      double inc = (inmax-inmin) / nbins;
      for (ulong i = 0; i < dims[0]*dims[1]*dims[2]; i++) {
	int bin = (int) ((im_in[i]-inmin) / inc);
	histos[bin]++;
      }

      ulong ncum = 0;
      int j =0;
      for(j = 0; j < dims[0]*dims[1]*dims[2]; j++){
	ncum += histos[j];
	if(ncum >= percentile) break;
      }
      printf("Bin for percentile %ld is [%lf, %lf) (Cumul is %ld)\n", percentile, (j-1)*inc+inmin, (j)*inc+inmin, ncum);

    }

     /**************************************************************************************/
    
    /*******************************  Two point corr   *********************************/
    
    /**************************************************************************************/

    else if( !strcmp(argv[1], "twoptcorr") ){
      if (argc<=3) {
	printf("Function <%s>: Get two pts correlation: <input name of file with positions> \n  <outputname> \n <nbins> \n <nbootstraps>", argv[1]);
	exit(0);
      }
      srand(time(NULL));
      FILE *fp1 = fopen(argv[2], "r");
      FILE *fp2 = fopen(argv[3], "w");

      int   bitpix = -32;    
      ulong dims[3];
      ulong nval, nvalr;
      fscanf(fp1, "%ld\n", &nval);
      printf("%ld positions to read \n", nval);
      fscanf(fp1, "0 %ld 0 %ld 0 %ld \n", dims, dims+1, dims+2);
      printf("Full dimensions: %ld  %ld  %ld \n", dims[0],dims[1],dims[2]);

      nvalr = 3*nval;
	
      double tempx, tempy, tempz, temp;
      double *x = calloc(nval, sizeof(double));
      double *y = calloc(nval, sizeof(double));
      double *z = calloc(nval, sizeof(double));

      for(ulong i = 0; i < nval; i++){
       	fscanf(fp1,"%lf %lf %lf %lf \n",&tempx, &tempy, &tempz, &temp);
	x[i] = tempx;
	y[i] = tempy;
	z[i] = tempz;
      }
      ulong size = dims[0]*dims[1]*dims[2];
      int nbins = atoi(argv[4]);
      int boots = atoi(argv[5]);
      double *count_dd = calloc(nbins, sizeof(double));
      double **count_dr = calloc(boots, sizeof(double *));
      double **count_rr = calloc(boots, sizeof(double *));
      double **corr = calloc(boots, sizeof(double *));
      for(ulong i = 0; i<nbins; i++){
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
	  double dist = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2));
	  int scale = find_scale_double(dist_bins, dist, nbins);
	  count_dd[scale]++;
	}
      }

      #pragma omp parallel for
      for(int k =0; k< boots; k++){
	double *x_rand = calloc(nvalr, sizeof(double));
	double *y_rand = calloc(nvalr, sizeof(double));
	double *z_rand = calloc(nvalr, sizeof(double));
	for(ulong l = 0; l < nvalr; l++){
	  x_rand[l] = keithRandom()*dims[0];
	  y_rand[l] = keithRandom()*dims[1];
	  z_rand[l] = keithRandom()*dims[2];
	}
	for (ulong i = 0; i < nvalr; i++) {
	  for(ulong j = 0; j < nvalr; j++) {
	    double dist;
	    int scale;
	    if(j > i){
	      dist = sqrt(pow(x_rand[i]-x_rand[j],2)+pow(y_rand[i]-y_rand[j],2)+pow(z_rand[i]-z_rand[j],2));
	      scale = find_scale_double(dist_bins, dist, nbins);
	      count_rr[k][scale]++;
	    }
	    if( j!= i && i<nval){
	      dist = sqrt(pow(x[i]-x_rand[j],2)+pow(y[i]-y_rand[j],2)+pow(z[i]-z_rand[j],2));
	      scale = find_scale_double(dist_bins, dist, nbins);
	      count_dr[k][scale]++;
	    }
	  }
	}
	for(int s = 0; s < nbins; s++)
	  corr[k][s] = (count_dd[s]/(nval*(nval-1)/2) - 2*count_dr[k][s]/(nvalr*nval) + count_rr[k][s]/(nvalr*(nvalr-1)/2)) / (count_rr[k][s]/(nvalr*(nvalr-1)/2));
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
	printf("Corr in bin %d is %lf +- %lf\n", j, corr_mean[j], corr_std[j]);
	fprintf(fp2, "%lf, %lf, %lf\n", dist_bins[j], corr_mean[j], corr_std[j]);
      }
      fclose(fp2);
    }
    
    /**************************************************************************************/
    
    /*******************************  Gaussian Filter  *********************************/
    
    /**************************************************************************************/

    
    else if( !strcmp(argv[1], "gaussian") ){
      if (argc<=4) {
	printf("Function <%s> Gaussian filter: <input name of image A without suffix> \n  <input type suffix image A> \n <x>, <y>, <z> \n", argv[1]);
	exit(0);
      }

      double fwhm = atof(argv[4]);
      ulong size = atoi(argv[5]);
      if(size % 2  == 0)
	size++;

      ulong dims[3] = {1,1,1};
      int bitpix;

      im_in = read_file(argv[2], argv[3], NULL, dims, &bitpix);	  		
      int n_dims = dims[2] > 1 ? 3 : 2;
      
      double *im_out = gaussian_filter(im_in, dims, fwhm, size);
      char *fname;
      asprintf(&fname, "%s_filtered", argv[2]);
      write_file(fname, "fits", "data", dims, im_out, -32, n_dims);

    }

    /**************************************************************************************/
    
    /*******************************  Power spectrum  *********************************/
    
    /**************************************************************************************/


    else if( !strcmp(argv[1], "ps") ){
      if (argc<=4) {
	printf("Function <%s>: Compute power spectrum: <input name of image A without suffix> \n  <input type suffix image A> \n <nfiles> \n  \n", argv[1]);
	exit(0);
      }

      ulong dims[3] = {256,256,256};
      int bitpix;

      im_in = read_file(argv[2], argv[3], NULL, dims, &bitpix);

      ulong size = dims[0]*dims[1]*dims[2];
      float *temp = calloc(size, sizeof(float));
      for(ulong i = 0; i< size; i++)
	temp[i] = (float) im_in[i];
      double tot = 0;
      for(ulong i = 0; i< dims[0]*dims[1]*dims[2]; i++)
	tot += (double) temp[i];
      tot /= (double) (dims[0]*dims[1]*dims[2]);
      float ave = (float) tot;
      printf("AVE in bof %lf ( box size %ld^3)\n", ave, dims[0]);
      int n_dims = dims[2] > 1 ? 3 : 2;
      float dims_box = atof(argv[4]);
      float volume = dims_box*dims_box*dims_box;
      float k_factor = 1.5;
      float delta_k = ((float)(2.0*PI))/((float)dims_box);
      float k_first_bin_ceil = delta_k;
      float k_max = delta_k*dims[0];
      int HII_DIM = (int) (dims[0]);
      int HII_MIDDLE = (int) (dims[0]/2);

      unsigned long long HII_D =  (unsigned long long) (dims[0]);
      unsigned long long HII_MID = ((unsigned long long) HII_D/2);

      unsigned long long kpixels=  (unsigned long long)(HII_D*HII_D*(HII_MID+1llu));


      int num_bins = 0;
      float k_floor = 0;
      float k_ceil = k_first_bin_ceil;
      while (k_ceil < k_max){
	num_bins++;
	k_floor=k_ceil;
	k_ceil*=k_factor;
      }

      double *p_box =  (double *)malloc(sizeof(double)*num_bins);
      double *k_ave =  (double *)malloc(sizeof(double)*num_bins);
      unsigned long long *in_bin_ct = (unsigned long long *)malloc(sizeof(unsigned long long)*num_bins);

      for (int ct=0; ct<num_bins; ct++){
	p_box[ct] = k_ave[ct] = 0;
	in_bin_ct[ct] = 0;
      }
      
      fftwf_plan plan;
      fftwf_complex *deldel_T;
      deldel_T = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*kpixels);
      // fill-up the real-space of the deldel box
      for (int i=0; i<HII_D; i++){
	for (int j=0; j<HII_D; j++){
	  for (int k=0; k<HII_D; k++){
	    unsigned long long HII_R_FFT_INDEX = ((unsigned long long)((k)+2llu*(HII_MID+1llu)*((j)+HII_D*(i)))) ;
	      unsigned long long HII_R_INDEX = ((unsigned long long)((k)+HII_D*((j)+HII_D*(i)))); // for 3D real array with no padding
	      *((float *)deldel_T + HII_R_FFT_INDEX) = ((float)(temp[ HII_R_INDEX])/ave - 1)*volume/(((unsigned long long)(HII_D*HII_D*HII_D))+0.0);
	      *((float *)deldel_T + HII_R_FFT_INDEX) *= ave;
	      // Note: we include the V/N factor for the scaling after the fft
	    }
	  }
	}

	// transform to k-space
	plan = fftwf_plan_dft_r2c_3d(HII_DIM, HII_DIM, HII_DIM, (float *)deldel_T, (fftwf_complex *)deldel_T, FFTW_ESTIMATE);
	fftwf_execute(plan);
	fftwf_destroy_plan(plan);
	fftwf_cleanup();
	float k_x, k_y, k_z;
	// now construct the power spectrum file
	for (int n_x=0; n_x<HII_DIM; n_x++){
	  if (n_x>HII_MIDDLE)
	    k_x =(n_x-HII_DIM) * delta_k;  // wrap around for FFT convention
	  else
	    k_x = n_x * delta_k;

	  for (int n_y=0; n_y<HII_DIM; n_y++){
	    if (n_y>HII_MIDDLE)
	      k_y =(n_y-HII_DIM) * delta_k;
	    else
	      k_y = n_y * delta_k;

	    for (int n_z=0; n_z<=HII_MIDDLE; n_z++){ 
	      k_z = n_z * delta_k;
	
	      float k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);

	      // now go through the k bins and update
	      int ct = 0;
	      k_floor = 0;
	      k_ceil = k_first_bin_ceil;
	      while (k_ceil < k_max){
		// check if we fal in this bin
		if ((k_mag>=k_floor) && (k_mag < k_ceil)){

		  in_bin_ct[ct]++;

		  unsigned long long HII_C_INDEX = ((unsigned long long)((n_z)+(HII_MID+1llu)*((n_y)+HII_D*(n_x))));// for 3D complex array;
		  p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T[HII_C_INDEX]), 2)/(2.0*PI*PI*volume);

		  // note the 1/VOLUME factor, which turns this into a power density in k-space

		  k_ave[ct] += k_mag;
		  break;
		}

		ct++;
		k_floor=k_ceil;
		k_ceil*=k_factor;
	      }
	    }
	  }
	} // end looping through k box
	char *fname;
	asprintf(&fname, "%s_ps", argv[2]);
	FILE *fp2 = fopen(fname, "w");
       

	for (int ct=1; ct<num_bins; ct++){
	  if (in_bin_ct[ct]>0){
	    printf("%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
	    fprintf(fp2, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
	  }
	}
	fclose(fp2);
	free(im_in);
	free(p_box);
	free(temp);
	//    }
    }
    else{
      printf("Incorrect function chosen, choices are:  create, mod_type, mod_depth, print_info, print_all, compare, bintofits, slice_plot, bw, lofar, density, position, cross_corr \n", argv[0]);
      exit(0);
    }
  }
  return 0;
}


