#include <hdf5.h>
#include <fitsio.h>
#include <FreeImage.h>
#include "types.h"
#include "nifti1.h"
#include "workspace.h"
#include "image.h"
#include "flood.h"

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352


FIBITMAP* freeimage_generic_loader(const char* lpszPathName, int flag);
const char *get_filename_ext(const char *filename);
void get_distributed_borders(value **img, int grid[3], ulong dims[3], ulong dims_T[3], bool border[6], int overlap);
void read_fits(Arguments *args, value **img, const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix);
void read_hdf5(Arguments *args, value **img, const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix);
void read_basic(Arguments *args,value **img, const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix);
void read_nifti_file(Arguments *args,value **img,  const char *fname, ulong dims[3], ulong dims_T[3],  int *bitpix);
void write_hdf5(Arguments *args, const char* fname, const char *dataset_out, value *out, ulong dims_T[3], ulong dims[3], bool border[6]);
void write_fits(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3], bool border[6]);
void write_basic(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3],  bool border[6]);


/*+++++++++++++++++++++++++++++++++++++++*/
/*     					 */
/*	     Functions Read     	 */
/*                                       */
/*+++++++++++++++++++++++++++++++++++++++*/

void read_input(Arguments *args, value **img,  ulong dims[3], ulong dims_T[3], ulong offsets[3], bool border[6]){

  /* +++++++++++++++++++++++++++ */
  /*       Main Read Function    */
  /* +++++++++++++++++++++++++++ */
  
  char 		*fname;
  int 		bitpix;
  int		myrank        = rank();
  int	 	grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  int 		tile   	      = args->infile_arg == 1 ?  myrank_arr[0] + (grid[1]-myrank_arr[1]-1)*grid[0] +myrank_arr[2]*grid[1]*grid[0]: myrank;

  // INIT BORDER ARRAY //

  border[0] = myrank_arr[0] % grid[0] != 0 ? true : false;
  border[1] = myrank_arr[0] % grid[0] < grid[0]-1 ? true : false;
  border[2] = myrank_arr[1] % grid[1] != 0 ? true : false;
  border[3] = myrank_arr[1] % grid[1] < grid[1]-1 ? true : false;
  border[4] = myrank_arr[2] % grid[2] != 0 ? true : false;
  border[5] = myrank_arr[2] % grid[2] < grid[2]-1 ? true: false;

  // Classic image reading  //
  
  if(!strcmp(args->image_arg, "classic")){
    if (args->infile_arg == 0)
      asprintf(&fname, "%s.%s", args->inprefix_arg, args->intype_arg);
    else if (args->infile_arg == 1)
      asprintf(&fname, "%s-%d.%s", args->inprefix_arg, tile, args->intype_arg);
    
      
    if (!strcmp(args->intype_arg, "h5")) 
      read_hdf5(args, img, fname, dims, dims_T, &bitpix);     
    else if (!strcmp(args->intype_arg, "fits") || !strcmp(args->intype_arg, "fit"))     
      read_fits(args, img, fname, dims, dims_T, &bitpix);
    else if (!strcmp(args->intype_arg, "nii"))
      read_nifti_file(args, img, fname, dims, dims_T, &bitpix);
    else  
      read_basic(args, img, fname, dims, dims_T, &bitpix);
  

    if (args->bpp_orig == NULL)
      args->bpp_arg = bitpix;
    else if (bitpix != args->bpp_arg && myrank == 0)
	warn("The user bit depth (%d) does not match the image header (%d)", args->bpp_arg, bitpix);
    
    free(fname);
    
  }

  // Reading RGB images //

  
  else if (!strcmp(args->image_arg, "rgb")){
    char 	  rgb[3]    = {'R','G','B'};
    float         coef[3]   = {0.2126, 0.7152, 0.0722};
    value         *img_curr = NULL;

    info("Reading RGB channels");
  
    if(!FLOAT_TYPE) {
      error("Need to activate the floating point values ! ");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    
    for(int i =3; i--;)
      {
	if (args->infile_arg == 0)
	  asprintf(&fname, "%s%d.%s", args->inprefix_arg, rgb[0], args->intype_arg);
	else if (args->infile_arg == 1)
	  asprintf(&fname, "%s%s-%d.%s", args->inprefix_arg, rgb[0], tile, args->intype_arg);
      
	if (!strcmp(args->intype_arg, "h5"))
	  read_hdf5(args, &img_curr, fname, dims, dims_T, &bitpix);     
	else if (!strcmp(args->intype_arg, "fits") || !strcmp(args->intype_arg, "fit"))     
	  read_fits(args, &img_curr, fname, dims, dims_T, &bitpix);   
	else   
	  read_basic(args, &img_curr, fname, dims, dims_T, &bitpix);  

	if(i == 2) *img = calloc(dims[0]*dims[1]*dims[2], sizeof(value));

	#pragma omp parallel for
	for(ulong j = 0; j < dims[0]*dims[1]*dims[2]; j++)
	  (*img)[j] += (value) coef[i]*img_curr[j];
      
	free(img_curr); 
      }
    
    args->bpp_arg = -32;
  }

  // Reading Lofar data //

  else if (!strcmp(args->image_arg, "lofar")){
    if(grid[0] > 1 || grid[1] > 1){
      error("Take 3D sequence, can't divide it in space");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    
    args->infile_arg 	= 1;

    value 	*img_curr 	= NULL;
    int 	nchannels 	= 2;
    int 	slices 		= 2;
    int 	nfiles 		= nchannels*slices;
    int 	offset[2] 	= {0};
    
    if(border[4]) offset[0]++;
    if(border[5]) offset[1]++;
    
    if(offset[0]) {
      // asprintf(&fname, "L254871-SB%3d-UV50_250_natural000%d-l-image_reweighted.fits", slices*rank()-1, nchannels-1);
      asprintf(&fname, "test%03d-000%d.fits",  slices*rank()-1, nchannels-1);
      read_fits(args, &img_curr, fname, dims, dims_T, &bitpix);   
      *img = malloc((nfiles+offset[0]+offset[1])*dims[0]*dims[1]*dims[2]*sizeof(value));	  
      memcpy(*img, img_curr, dims[0]*dims[1]*dims[2]*sizeof(value));
      free(img_curr);
    }
    
    for(int i = 0; i < slices; i++){
      for(int j =0; j< nchannels; j++){
	asprintf(&fname, "test%03d-000%d.fits",  slices*rank()+i,j);
	//asprintf(&fname, "L254871-SB%3d-UV50_250_natural000%d-l-image_reweighted.fits", slices*rank()+i,j);

	read_fits(args, &img_curr, fname, dims, dims_T, &bitpix);   
	if(!i && !j&& !offset[0])
	  *img = malloc((nfiles+offset[0]+offset[1])*dims[0]*dims[1]*dims[2]*sizeof(value));	  
	memcpy(*img+dims[0]*dims[1]*dims[2]*(i*nchannels+j+offset[0]), img_curr, dims[0]*dims[1]*dims[2]*sizeof(value));
	free(img_curr);
      }
    }
    
    if(offset[1]){
      //asprintf(&fname, "L254871-SB%3d-UV50_250_natural0000-l-image_reweighted.fits", slices*(rank()+1));
      asprintf(&fname, "test%03d-0000.fits",  slices*(rank()+1));
      read_fits(args, &img_curr, fname, dims, dims_T, &bitpix);   
      memcpy(*img+dims[0]*dims[1]*dims[2]*(nfiles+offset[0]), img_curr, dims[0]*dims[1]*dims[2]*sizeof(value));
      free(img_curr);
    }
    
    dims[2] = nfiles+offset[0]+offset[1];
     args->bpp_arg = bitpix;

     
  }


  if(args->infile_arg)
    get_distributed_borders(img, grid, dims, dims_T, border, args->overlap_arg);

  // Init offsets for attributes //
  
  ulong		counts[3]     = {0};
  
    for (int i = 3; i-- ; ) {
    counts[i]  = dims_T[i] / grid[i];
    offsets[i] = myrank_arr[i]*counts[i];
    if (myrank_arr[i] < dims_T[i]%grid[i]) {
      counts[i]++;
      offsets[i] += myrank_arr[i];
    } else {
      offsets[i] += dims_T[i]%grid[i];
    }
    counts[i] += offsets[i] - 1;
    if((offsets[i] > 1)){
      offsets[i]--;
    }
  }

  
  // Init connectivity //
  
  if (args->connectivity_orig == NULL)
    args->connectivity_arg = dims[2] > 1 ? 6 : 4;
  else if(args->connectivity_arg != 4 && args->connectivity_arg != 6 && args->connectivity_arg != 8 && args->connectivity_arg != 26){
    if(rank() == 0){
      error("Wrong Connectivity");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
  } else if( args->connectivity_arg == 8 && dims[2] > 1){
    if(rank() == 0) warn("Wrong connectivity (8) for data dimension (3D), changing to 26");
    args->connectivity_arg = 26;
  } else if( args->connectivity_arg == 26 && dims[2] == 1){
    if(rank() == 0) warn("Wrong connectivity (26) for data dimension (2D), changing to 8");
    args->connectivity_arg = 8;
  }

  info("Data initialized");

} /* read_input */


void read_fits(Arguments *args, value **img,  const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix){
  /* +++++++++++++++++++++++++++ */
  /*          Read Fits          */
  /* +++++++++++++++++++++++++++ */

  fitsfile 	*fptr;   			/* FITS file pointer */
  int		n_dims;
  int 		type;
  int 		status        = 0;   		/* CFITSIO status value MUST be initialized to zero! */
  int	 	grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank        = rank();		/* Rank of process */
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  long 		inc[4]        = {1,1,1,1};	/* Increment for fits read function */
  long 		naxes[4]      = {1,1,1,1};	/* Dimensions (4D max) */
  long 		counts[4]     = {1,1,1,1};	/* Pixels to read in fits function */
  long 		offsets[4]    = {1,1,1,1};	/* Pixel index to start in fist function */
  int  		overlap       = args->overlap_arg;
  int 		negative      = args->negative_arg;
  int 		infile        = args->infile_arg;
  int           ttype         = sizeof(value) == 1 ? 1-negative: negative;
  
  info("Reading FITS Image %s", fname);

  if(FLOAT_TYPE)
    type = TFLOAT;
  else
    type = (int) (log2(sizeof(value))+1)*10 + ttype;

  if (!fits_open_file(&fptr, fname, READONLY, &status)) {	/* Open file */
    if (!fits_get_img_param(fptr, 4,  bitpix,  &n_dims, naxes, &status)) {  /* Get image parameter */
      if (n_dims > 4 || n_dims == 0 || (n_dims == 4 && naxes[3] > 1)) { /* 3D Max */	   	      
	error("Only 2D and 3D images are supported");
	MPI_Abort(MPI_COMM_WORLD, 701);
      }
      else {	
	if(n_dims == 4) (n_dims)--;
	
	if(!infile){
	  if(n_dims == 2 && grid[2] > 1){
	    error("Ask distributionin depth but data is 2D");
	    MPI_Abort(MPI_COMM_WORLD, 701);
	  }
	  for (int i = n_dims; i-- ; ) {
	    counts[i]  = naxes[i] / grid[i];
	    offsets[i] = myrank_arr[i]*counts[i] + 1;
	    if (myrank_arr[i] < naxes[i]%grid[i]) {
	      counts[i]++;
	      offsets[i] += myrank_arr[i];
	    } else {
	      offsets[i] += naxes[i]%grid[i];
	    }
	    dims_T[i]  = naxes[i];
	    dims[i]    = counts[i];
	    counts[i] += offsets[i] - 1;
	    if((offsets[i] > 1) && overlap){
	      offsets[i]--;
	      dims[i]++;
	    }
	    if(((ulong) counts[i] != dims_T[i]) && overlap) {
	      counts[i]++;				       			       
	      dims[i]++;
	    }
	  }
	}
	else {	  
	  for (int i = n_dims; i-- ; ) 
	    counts[i]  = dims[i] = naxes[i];

	}
	debug("Fits reading: \n Tile dimensions: %ld by %ld by %ld \n Offsets %ld, %ld, %ld \n Counts %ld, %ld, %ld", dims[0], dims[1], dims[2],offsets[0],offsets[1],offsets[2],counts[0],counts[1],counts[2]);
	*img = malloc(dims[0]*dims[1]*dims[2]* sizeof(value));
	fits_read_subset(fptr, type, offsets, counts, inc, NULL, *img, NULL, &status);
      }
    }
    
    fits_close_file(fptr, &status);        
  } else {
    error("Cannot open the file %s, wrong name ? (FITS)", fname);
    MPI_Abort(MPI_COMM_WORLD, 703);
  }    
}


void read_hdf5(Arguments *args, value **img, const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix){
  
  /* +++++++++++++++++++++++++++ */
  /*      Read file HDF5        */
  /* +++++++++++++++++++++++++++ */
  
  H5T_class_t   t_class;
  hid_t       	file_id, dataset_id, dataspace;  	    
  hsize_t 	counts[3]     = {1,1,1};
  hsize_t 	offsets[3]    = {0};
  hsize_t 	*hdims;
  hsize_t 	hn_dims;
  herr_t	err;
  int	 	grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank        = rank();		/* Rank of process */
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  int  		overlap       = args->overlap_arg;
  int 		negative      = args->negative_arg;
  int 		infile        = args->infile_arg;
  
  info("Reading HDF5 file %s", fname);

  /*      Open HDF5 file and dataset     */
  
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);	
  if (file_id < 0) {
    error("Could not open file %s, wrong name ? (HDF5) ", fname);
    MPI_Abort(MPI_COMM_WORLD, 707);
  }
  
  dataset_id = H5Dopen2(file_id, args->dataset_arg, H5P_DEFAULT);	
  if (dataset_id < 0) {
    error("Could not open dataset %s, wrong name ? (HDF5) ", args->dataset_arg);
    MPI_Abort(MPI_COMM_WORLD, 708);
  }

  
  dataspace = H5Dget_space(dataset_id);
  hid_t tid = H5Dget_type(dataset_id);
  t_class   = H5Tget_class(tid);
  int ord   = H5Tget_precision(tid);
  
  if(t_class == H5T_FLOAT)
    *bitpix = -32;
  else {  
    if(ord == 8)
      *bitpix = 8;
    else if(ord == 16)
      *bitpix = 16;
    else 
      *bitpix = 32;
  }
  
  hn_dims = H5Sget_simple_extent_ndims(dataspace);
  
  /*       Get number of dimensions  and dimensions    */
  
  if (hn_dims > 3) {
    error("Only handle 2D or 3D data");
    MPI_Abort(MPI_COMM_WORLD, 709);
  }

  hdims = calloc(hn_dims, sizeof(hsize_t));			
  H5Sget_simple_extent_dims(dataspace, hdims, NULL);

  if(infile){
    for (hsize_t i = 0; i < hn_dims; i++)
      counts[i] = dims[hn_dims - i - 1] = hdims[i];

  } else {
    if(hn_dims == 2 && grid[2] > 1){
      error("Ask distribution in depth but data is 2D");
      MPI_Abort(MPI_COMM_WORLD, 701);
    }
    for (hsize_t i = 0; i < hn_dims; i++) {
      counts[i]  = hdims[i]/grid[hn_dims - i - 1];
      offsets[i] = myrank_arr[hn_dims - i - 1] * counts[i];
      if (myrank_arr[hn_dims - i - 1] < (int) (hdims[i]%grid[hn_dims - i - 1])) {
	counts[i]++;
	offsets[i] += myrank_arr[hn_dims - i - 1];
      } else {
	offsets[i] += (hdims[i]%grid[hn_dims - i - 1]);
      }
      if((offsets[i] > 0) && overlap){
	offsets[i]--;
	counts[i]++;
      }
      if((counts[i] + offsets[i] != hdims[i]) && overlap) {
	counts[i]++;				       			       
      }
      dims_T[hn_dims - i - 1] = hdims[i];
      dims[hn_dims - i - 1]   = counts[i];
    }

  }
  
  debug("HDF5: tile dimensions: depth %ld, height %ld, width %ld.\n First pixel offset: depth %ld, height %ld, width %ld", counts[0], counts[1], counts[2], offsets[0], offsets[1], offsets[2]);
  
  *img = malloc(counts[0]*counts[1]*counts[2] * sizeof(value));   check_alloc(img,710);
  
  hid_t memory_window = H5Screate_simple(hn_dims, counts, NULL); /* Allocate memory to read the data */
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL); /* Select the data subset */
  if(FLOAT_TYPE)
    err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 8)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_LONG, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_ULONG, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 4)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_INT, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_UINT, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 2)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_SHORT, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_USHORT, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 1)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_CHAR, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_UCHAR, memory_window, dataspace, H5P_DEFAULT, *img);

  if (err < 0){ 
    error("Could not read data from the dataset (HDF5)");
    MPI_Abort(MPI_COMM_WORLD, 711);
  }

  free(hdims);
  
  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);

  //  return img;
}


void read_basic(Arguments *args,value **img, const char *fname, ulong dims[3], ulong dims_T[3],  int *bitpix){
  
  /* +++++++++++++++++++++++++++ */
  /*       Read FreeImage        */
  /* +++++++++++++++++++++++++++ */
  
  info("Reading image %s with FreeImage", fname);
  
  ulong 	i   	      = 0;
  FIBITMAP 	*dib 	      = freeimage_generic_loader(fname, 0);  /* Generic Loader from FreeImage */
  int	 	grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank        = rank();		/* Rank of process */
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  ulong 	offsets[2]    = {0};
  int  		overlap       = args->overlap_arg;
  int 		infile        = args->infile_arg;
  
  if (dib == NULL) {
    error("Can't read the image %s (FreeImage)", fname);
    MPI_Abort(MPI_COMM_WORLD, 717);
  }

  FREE_IMAGE_COLOR_TYPE nc = FreeImage_GetColorType(dib);

  /*  if(nc != 1){
      FIBITMAP* hImage= dib;
      dib = FreeImage_ConvertTo32Bits( hImage );
      hImage= dib;
      dib = FreeImage_ConvertTo8Bits( hImage );
      FreeImage_Unload( hImage );
      }*/
  
  if(nc != 1)
    dib =  FreeImage_ConvertToGreyscale(dib);
  *bitpix = FreeImage_GetBPP(dib);				/* Get dynamic range */

  
  if(!infile){
    dims_T[0] = FreeImage_GetWidth(dib);			    
    dims_T[1] = FreeImage_GetHeight(dib);
    if(grid[2] > 1){
      error("Ask distributionin depth but data is 2D");
      MPI_Abort(MPI_COMM_WORLD, 701);
    }
    for (i = 0; i < 2; i++){
      dims[i]  = dims_T[i]/grid[i];
      offsets[i] = myrank_arr[i] * dims[i];
      if (myrank_arr[i] < (int) (dims_T[i]%grid[i])) {
	dims[i]++;
	offsets[i] += myrank_arr[i];
      } else {
	offsets[i] += (dims_T[i]%grid[i]);
      }
      if((offsets[i] > 0) && overlap){
	offsets[i]--;
	dims[i]++;
      }
      if((dims[i] + offsets[i] != dims_T[i]) && overlap) {
	dims[i]++;
      }
    }
  } else {
    dims[0] = FreeImage_GetWidth(dib);				/* Get width */
    dims[1] = FreeImage_GetHeight(dib);
  }

  *img = malloc(dims[0]*dims[1] * sizeof(value));	check_alloc(img,718);

  i = 0;
  if (*bitpix <= 8) {
    for(ulong y = offsets[1]; y < dims[1] +offsets[1]; y++) {
      BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, y);
      for(ulong x = offsets[0]; x < dims[0] + offsets[0]; x++,i++){
	(*img)[i] = (value) bits[x];
      }
    }
    FreeImage_Unload(dib);
    //   return img;
    
  } else if (*bitpix <= 16) {
    for(ulong y = offsets[1]; y < dims[1] + offsets[1]; y++) {
      ushort *bits = (ushort *)FreeImage_GetScanLine(dib, y);
      for(ulong x = offsets[0]; x < dims[0] + offsets[0]; x++,i++)
	(*img)[i] = (value) bits[x];
    }
    FreeImage_Unload(dib);
    return img;  
  } else if (*bitpix <= 32) {
    for(ulong y = offsets[1]; y < dims[1] + offsets[1]; y++) {
      uint *bits = (uint *)FreeImage_GetScanLine(dib, y);
      for(ulong x = offsets[0]; x < dims[0] + offsets[0]; x++,i++)
	(*img)[i] = (value) bits[x];
    }
    FreeImage_Unload(dib);
    //   return img;  
  } else{
    error("Floating point or 64 bits depth are not supported with classical images (Freeimage)");
    MPI_Abort(MPI_COMM_WORLD, 719);
    FreeImage_Unload(dib);
    //  return NULL;
  }
}

void read_nifti_file(Arguments *args, value **img, const char *fname, ulong dims[3], ulong dims_T[3],  int *bitpix){
  nifti_1_header hdr;
  FILE *fp;
  int ret,i;
  double total;
  value *data=NULL;


  /********** open and read header */
  fp = fopen(fname,"r");
  if (fp == NULL) {
    fprintf(stderr, "\nError opening header file %s\n",fname);
    exit(1);
  }
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  if (ret != 1) {
    fprintf(stderr, "\nError reading header file %s\n",fname);
    exit(1);
  }
  fclose(fp);


  /********** print a little header information */
  fprintf(stderr, "\n%s header information:",fname);
  fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
  fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
  fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
  fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
  fprintf(stderr, "\n");

  *bitpix = hdr.bitpix;
  
  /********** open the datafile, jump to data offset */
  fp = fopen(fname,"r");
  if (fp == NULL) {
    error("Error opening data file %s",fname);
    exit(1);
  }

  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
  if (ret != 0) {
    error("Error doing fseek() to %ld in data file %s",(long)(hdr.vox_offset), fname);
    exit(1);
  }


/********** allocate buffer and read first 3D volume from data file */
  *img = (value *) malloc(sizeof(value) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
  if (*img == NULL) {
    error("Error allocating data buffer for %s",fname);
    exit(1);
  }
  if(*bitpix < 0)
    ret = fread(*img, sizeof(float), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  else if(*bitpix <=8)
    ret = fread(*img, sizeof(ubyte), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  else if(*bitpix <=16)
    ret = fread(*img, sizeof(ushort), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  else if(*bitpix <=32)
    ret = fread(*img, sizeof(uint), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
    error("Error reading volume 1 from %s (%d)");
    exit(1);
  }
  fclose(fp);


  /********** scale the data buffer  */
 if (hdr.scl_slope != 0) {
   for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
     (*img)[i] = ((*img)[i] * hdr.scl_slope) + hdr.scl_inter;
 }


/********** print mean of data */
 total = 0;
 for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
   total += (*img)[i];
 total /= (hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
 fprintf(stderr, "\nMean of volume 1 in %s is %.3f\n",fname,total);

 dims[0] = dims_T[0] =  hdr.dim[1];
 dims[1] = dims_T[1] =  hdr.dim[2];
 dims[2] = dims_T[2] =  hdr.dim[3];

 //return(data);
}

/*+++++++++++++++++++++++++++++++++++++++*/
/*     					 */
/*	     Functions Write      	 */
/*                                       */
/*+++++++++++++++++++++++++++++++++++++++*/

void write_output(Arguments *args, value *img, const char *attr_name, ulong dims_T[3], ulong dims[3],  bool border[6]){
  
  /* +++++++++++++++++++++++++++ */
  /*     Main Write Function     */
  /* +++++++++++++++++++++++++++ */
  
  char 		*fname_out;
  char		*fname_in;
  char 		*dataset_out;
  int	 	grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank     = rank();
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};   int   	tile = args->infile_arg || args->outfile_arg ?  myrank_arr[0] + (grid[1]-myrank_arr[1]-1)*grid[0] + myrank_arr[2]*grid[2]*grid[1] : myrank;
  
  if (args->outfile_arg == 0)
    asprintf(&fname_out, "%s.%s", args->outprefix_arg, args->outtype_arg);
  else{
    asprintf(&fname_out, "%s-%d.%s", args->outprefix_arg, tile, args->outtype_arg);
    if(rank() == 0)
      warn("Output without overlapping pixels !");
  }

  if (!strcmp(args->outtype_arg, "h5")){
    asprintf(&dataset_out, "%s_%s_%1.0lf", args->output_arg, attr_name, args->lambda_arg);
    write_hdf5(args, fname_out, dataset_out, img, dims_T, dims, border);
  }  else if (!strcmp(args->outtype_arg, "fits") || !strcmp(args->outtype_arg, "fit"))     
    write_fits(args, fname_out, img, dims_T, dims, border);      
  else 
    write_basic(args, fname_out, img, dims_T, dims, border);   
  
} /* write_output */


void write_differential(Arguments *args, value *outOrig, value *outDH, value *outScale, const char *attr_name, ulong dims_T[3], ulong dims[3],  bool border[6]) {
  /* Specific  differential profile function     */
  
  char		*fname_out;
  char 		*dataset_out;
  int	 	grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank     = rank();
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  int   	tile = args->infile_arg || args->outfile_arg ?  myrank_arr[0] + (grid[1]-myrank_arr[1]-1)*grid[0] + myrank_arr[2]*grid[2]*grid[1] : myrank;
  char 		step[3]    = {'C','S','L'};
  value 	*out;
  
  for(int i = 0; i<3; i++ ){
    if (args->outfile_arg == 0)
      asprintf(&fname_out, "%s-%c.%s", args->outprefix_arg, step[i], args->outtype_arg);
    else{
      asprintf(&fname_out, "%s-%s-%d.%s", args->outprefix_arg, step[i], tile, args->outtype_arg);
      if(rank() == 0 && i == 3)
	warn("Output without overlapping pixels !");
    }
    
    if(i)
      out = i == 1 ? outScale : outOrig ;
    else
      out = outDH;

    if (!strcmp(args->outtype_arg, "h5")) {
      asprintf(&dataset_out, "%s_%s", args->output_arg, attr_name);
      write_hdf5(args, fname_out, dataset_out, out,  dims_T, dims, border);     
    } else if (!strcmp(args->outtype_arg, "fits") || !strcmp(args->outtype_arg, "fit"))
      write_fits(args, fname_out, out,  dims_T, dims, border);
    else 
      write_basic(args, fname_out, out, dims_T, dims, border);   
    
    free(fname_out);

    MPI_Barrier(MPI_COMM_WORLD);
  }
} /* write_differential */




void write_hdf5(Arguments *args, const char* fname, const char *dataset_out, value *out, ulong dims_T[3], ulong dims[3],  bool border[6]) {
  
  /* +++++++++++++++++++++++++++ */
  /*     HDF5 Write Function     */
  /* +++++++++++++++++++++++++++ */  
  hid_t       	file_id, dataset_id, dataspace, type;  /* identifiers */
  hsize_t	hdims[3]      = {1,1,1};
  hsize_t 	counts[3]     = {1,1,1};
  hsize_t 	offsets[3]    = {0};

  int         	message;
  int 		grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank        = rank();		
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  int 		n_dims 	      = dims[2] > 1 ? 3 : 2;
  int 		negative      = args->negative_arg;
  int 		outfile       = args->outfile_arg;
  int 		bitpix        = args->bpp_arg;
  
  debug("Writing HDF5 file %s", fname);

  if(bitpix < 0)
    type = H5T_NATIVE_FLOAT;
  else if(bitpix <= 8)
    type =  negative ? H5T_NATIVE_CHAR: H5T_NATIVE_UCHAR;
  else if(bitpix <= 16)
    type =  negative ? H5T_NATIVE_SHORT: H5T_NATIVE_USHORT;
  else if(bitpix <= 32)
    type =  negative ? H5T_NATIVE_INT: H5T_NATIVE_UINT;
  else{
    error("Non supported data type");
    MPI_Abort(MPI_COMM_WORLD, 726);
  }

  if(outfile){
    for (int i = 0; i < n_dims; i++){
      counts[i] = dims[n_dims - i - 1]-border[2*(n_dims - i -1)]-border[2*(n_dims - i -1)+1];
      hdims[i]   = dims[n_dims - 1 - i]-border[2*(n_dims - i -1)]-border[2*(n_dims - i -1)+1];
    }
  } else {   
    for(int i = 0; i < n_dims; i++){
      hdims[i]   = dims_T[n_dims - 1 - i];
      counts[i]  = hdims[i] / grid[n_dims - i - 1];
      offsets[i] = myrank_arr[n_dims - i - 1] * counts[i];
      if(myrank_arr[n_dims-i-1] < (hdims[i] % grid[n_dims-i-1])){
	counts[i]++;
	offsets[i] += myrank_arr[n_dims - i - 1];
      } else
	offsets[i] += (hdims[i] % grid[n_dims - i - 1]);
    }
  }
  
  if(myrank == 0 || outfile) {
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);    
    hid_t global_memory_space = H5Screate_simple(n_dims, hdims, NULL);
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    dataset_id = H5Dopen(file_id, dataset_out, H5P_DEFAULT);		/* Open dataset */
    if (dataset_id < 0) 
      dataset_id = H5Dcreate(file_id, dataset_out, type, global_memory_space,
			     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	/* Create dataset */
    
    dataspace = H5Dget_space(dataset_id);			/* Copy dataset */
  } else {    
    MPI_Recv(&message, 1, MPI_INT, myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id,dataset_out, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset_id);   
  }
   
  if(FLOAT_TYPE == 1)
    type =  H5T_NATIVE_FLOAT;
  else if(sizeof(value) == 1)
    type =  negative ? H5T_NATIVE_CHAR: H5T_NATIVE_UCHAR;
  else if(sizeof(value) == 2)
    type =  negative ? H5T_NATIVE_SHORT: H5T_NATIVE_USHORT;
  else if(sizeof(value) == 4)
    type =  negative ? H5T_NATIVE_INT: H5T_NATIVE_UINT;
  
  ulong off_out  = border[0] + border[2]*dims[0] + border[4]*dims[0]*dims[1];
  ulong off_save = offsets[1];
  ulong copy[3]  = {counts[0], counts[1], counts[2]};
  ulong s 	   = 0;

  for(ulong i = 0; i < copy[0]; i++){
    counts[0] = 1;
    if(n_dims == 2) {
      hid_t memory_window = H5Screate_simple((hsize_t)n_dims, counts, NULL);
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
      H5Dwrite(dataset_id, type, memory_window, dataspace, H5P_DEFAULT, out+i*dims[0]+off_out);
    } else {
      offsets[1]= off_save;
      for(ulong j = 0; j < copy[1]; j++, s++){
	counts[1] = 1;
	hid_t memory_window = H5Screate_simple((hsize_t)n_dims, counts, NULL);
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
	H5Dwrite(dataset_id, type, memory_window, dataspace, H5P_DEFAULT, out+off_out+(s+i*(border[2]+border[3]))*dims[0]);
	offsets[1]++;
      }
    }
    offsets[0]++;
  }
  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);
  if((myrank < np()-1) && !outfile)
    MPI_Send(&myrank, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);
} /* write_hdf5 */


 
void write_fits(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3], bool border[6]){
  /* +++++++++++++++++++++++++++ */
  /*     FITS Write Function     */
  /* +++++++++++++++++++++++++++ */
  
  fitsfile 	*outfptr; 	/* FITS file pointers */
  int 		status = 0; 		/* CFITSIO status value MUST be initialized to zero! */
  int 		message;
  int 		type;
  int 		grid[3]       = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		myrank        = rank();	
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  long 		naxes[3]      = {1,1,1};      
  long 		counts[3]     = {1,1,1};
  long 		offsets[3]    = {1,1,1};
  int 		n_dims 	      = dims[2] > 1 ? 3 : 2;
  int 		negative      = args->negative_arg;
  int 		outfile       = args->outfile_arg;
  int 		bitpix        = args->bpp_arg;
  int           ttype         = sizeof(value) == 1 ? 1-negative: negative;

  char 		str1[100];
  strcpy(str1, "!");   		/* '!' symbol makes the output image, if existing, to be overwritten */
  if(outfile) {
    for (int i = 0; i < n_dims; i++) 
      naxes[i] = counts[i] = dims[i]-border[2*i]-border[2*i+1];     //-border[2*i]-border[2*i+1]
  } else {  
    for (int i = 0; i < n_dims; i++) {
      naxes[i]   = dims_T[i];
      counts[i]  = naxes[i] / grid[i];
      offsets[i] = myrank_arr[i]*counts[i] + 1;
      if (myrank_arr[i] < naxes[i]%grid[i]) {
	counts[i]++;
	offsets[i] += myrank_arr[i];
      } else {
	offsets[i] += naxes[i]%grid[i];
      }
      counts[i] += offsets[i] - 1;
    }
  }
  if(FLOAT_TYPE || args->bpp_arg == -32){
    type = TFLOAT;
    out = (float *) out;
  }else
    type = (int) (log2(sizeof(value))+1)*10 + ttype;
	
  if ((myrank == 0) || outfile) {
    if (!fits_create_file(&outfptr, strcat(str1, fname), &status)) {
      if (bitpix < 0)
	fits_create_img(outfptr, FLOAT_IMG, n_dims, naxes, &status);
      else if (bitpix <= 8)
	negative ? fits_create_img(outfptr, SBYTE_IMG, n_dims, naxes, &status):
	  fits_create_img(outfptr, BYTE_IMG, n_dims, naxes, &status);
      else if (bitpix <= 16)
	negative ? fits_create_img(outfptr, SHORT_IMG, n_dims, naxes, &status):
	  fits_create_img(outfptr, USHORT_IMG, n_dims, naxes, &status);
      else if (bitpix <= 32)
	negative ? fits_create_img(outfptr, LONG_IMG, n_dims, naxes, &status):
	  fits_create_img(outfptr, ULONG_IMG, n_dims, naxes, &status);
      else 
	fits_create_img(outfptr, LONGLONG_IMG, n_dims, naxes, &status);
      
      fits_update_key(outfptr, TSTRING , "FILTER", args->output_arg, "Type of filter used", &status);
      fits_update_key(outfptr, TSTRING , "ATTRIBUTE", AttribsArray[args->attribute_arg].name, "Attribute name", &status);
      fits_update_key(outfptr, TSTRING , "DECISION", Decisions[args->decision_arg].name, "Pruning rule", &status);
      fits_update_key(outfptr, TINT , "CONNECTIVITY",  &args->connectivity_arg, "Neighbour connectivity", &status);
      if(!strcmp(args->output_arg, "filter"))
	fits_update_key(outfptr, TULONG , "LAMBDA", &args->lambda_arg, "Threshold value", &status);
      else
	fits_update_key(outfptr, TSTRING , "LVEC", args->lvec_arg, "Name of threshold file", &status);
      fits_update_key(outfptr, TSTRING , "COMMENT", "Filtering using connected filter method DISCCOMAN", NULL, &status);
     
    } else {
      error("Cannot create the file");
      MPI_Abort(MPI_COMM_WORLD, 730);
    }
  } else {
    MPI_Recv(&message, 1, MPI_INT, myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    fits_open_file(&outfptr, fname, READWRITE, &status); // open input images
  }

  //  if(!outfile){
  ulong off_out   = border[0] + border[2]*dims[0] + border[4]*dims[0]*dims[1];
  ulong off_save  = offsets[1];
  ulong s	  = 0;
  counts[1] = off_save;
  counts[2] = offsets[2];

  for(ulong i = 0; i < dims[2]-border[5]-border[4]; i++){
    offsets[1] = off_save;
    counts[1] = off_save;
    for(ulong j = 0; j < dims[1]-border[3]-border[2]; j++,s++){
      fits_write_subset(outfptr, type, offsets, counts, out+off_out+(s+i*(border[2]+border[3]))*dims[0], &status);
      offsets[1]++;
      counts[1]++;
    }
    counts[2]++;
    offsets[2]++;
  }
  // }else{
  //   fits_write_subset(outfptr, type, offsets, counts, out, &status);
    // }
  fits_close_file(outfptr, &status);
  if((myrank < np()-1) && !outfile)
    MPI_Send(&status, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);
} /* write_fits */




void write_basic(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3],  bool border[6]){
  
  /* +++++++++++++++++++++++++++ */
  /*        Write FreeImage      */
  /* +++++++++++++++++++++++++++ */

  FIBITMAP 	*outmap;
  ulong 	width  	      = dims[0];
  ulong 	height 	      = dims[1];
  int 		negative      = args->negative_arg;
  int 		outfile       = args->outfile_arg;
  int 		bitpix        = args->bpp_arg;

  if(dims[2] > 1 || (dims_T[2] > 1 && !outfile)){
    if(rank() == 0) error("FreeImage does not handle writing 3D images, using FITS");
    args->outfile_arg == 0 ? asprintf(&fname, "%s.fits", args->outprefix_arg):
      asprintf(&fname, "%s-%d.fits", args->outprefix_arg, rank());
    write_fits(args, fname, out, dims_T, dims, border);
    return;
  }
  if(!outfile && np() > 1){
    if(rank() == 0)error("FreeImage does not handle writing 3D images, using FITS");
    args->outfile_arg == 0 ? asprintf(&fname, "%s.fits", args->outprefix_arg):
      asprintf(&fname, "%s-%d.fits", args->outprefix_arg, rank());
    write_fits(args, fname, out, dims_T, dims, border);
    return;
  }
  if(bitpix < 0){
    if(rank() == 0) error("We do not handle floating point writing with FreeImage, using FITS instead");
    args->outfile_arg == 0 ? asprintf(&fname, "%s.fits", args->outprefix_arg):
      asprintf(&fname, "%s-%d.fits", args->outprefix_arg, rank());
    write_fits(args, fname, out, dims_T, dims, border);
    return;
  }
  
  debug("Writing image %s using FreeImage", fname);	  


  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname);
  if (fif == FIF_UNKNOWN) 
    fif = FreeImage_GetFIFFromFilename(fname);
  if ((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsWriting(fif)) {  
    if (bitpix <= 8) {   
      ubyte *imagebuf = malloc(width*height* sizeof(ubyte));;
      outmap =  FreeImage_AllocateT(FIT_BITMAP, width-border[0]-border[1], height-border[2]-border[3], 8, 0xFF, 0xFF, 0xFF);
      for (ulong y = border[2]; y < height - border[3]; y++) {
	imagebuf = (BYTE *) FreeImage_GetScanLine(outmap, y - border[2]);      
	for (ulong x = border[0]; x < width-border[1]; x++)
	  imagebuf[x-border[0]] = out[width*y + x];	
      }
    } else if (bitpix <= 16) {
      if(!negative){
	ushort *imagebuf;
	outmap = FreeImage_AllocateT(FIT_UINT16, width-border[0]-border[1], height-border[2]-border[3], 16,0xFFFF,0xFFFF,0xFFFF);
	for (ulong y = border[2]; y < height - border[3]; y++) {
	  imagebuf = (ushort *) FreeImage_GetScanLine(outmap, y - border[2]);  
	  for (ulong x = border[0]; x < width-border[1]; x++)
	    imagebuf[x-border[0]] = out[width*y + x];
	}
      } else {
	short *imagebuf;
	outmap = FreeImage_AllocateT(FIT_INT16,width-border[0]-border[1], height-border[2]-border[3], 16,0xFFFF,0xFFFF,0xFFFF);
	for (ulong y = border[2]; y < height - border[3]; y++) {
	  imagebuf = (short *) FreeImage_GetScanLine(outmap, y - border[2]);  
	  for (ulong x = border[0]; x < width-border[1]; x++)
	    imagebuf[x-border[0]] = out[width*y + x];
	}
      }
    } else if (bitpix <= 32) {
      if(!negative){
	uint *imagebuf;
	outmap = FreeImage_AllocateT(FIT_UINT32, width-border[0]-border[1], height-border[2]-border[3], 32,0xFFFF,0xFFFF,0xFFFF);
	for (ulong y = border[2]; y < height - border[3]; y++) {
	  imagebuf = (uint *) FreeImage_GetScanLine(outmap, y - border[2]);  
	  for (ulong x = border[0]; x < width-border[1]; x++)
	    imagebuf[x-border[0]] = out[width*y + x];
	}
      } else {
	int *imagebuf;
	outmap = FreeImage_AllocateT(FIT_INT32, width-border[0]-border[1], height-border[2]-border[3], 32,0xFFFF,0xFFFF,0xFFFF);
	for (ulong y = border[2]; y < height - border[3]; y++) {
	  imagebuf = (int *) FreeImage_GetScanLine(outmap, y - border[2]);  
	  for (ulong x = border[0]; x < width-border[1]; x++)
	    imagebuf[x-border[0]] = out[width*y + x];
	}
      }
    }      
    FreeImage_Save(fif,outmap,fname,0); 
    FreeImage_Unload(outmap);
  }
  else{
    error("FreeImage couldn't write the output");
    MPI_Abort(MPI_COMM_WORLD, 725);
  }
}

 

void write_pattern_spectra(Arguments *args, double* spectrum,  LambdaVec *lvec){
  char *filename;
  int numscales = lvec->num_lambdas;
  if(args->outprefix_orig != NULL)
    asprintf(&filename, "%s.txt", args->outprefix_arg);
  else
    asprintf(&filename, "pattern.txt");
	       
  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }
  fprintf(f, "#lambdas \t spectra \n");
  //fprintf(f, "#lambdas \t spectra \n");

  for (ulong i = 0; i < (ulong) numscales; i++) 
    //fprintf(f, "Spec[%ld]= %lf \n", i, spectrum[i]);
    fprintf(f, "%f \t %lf \n", lvec->lambdas[i], spectrum[i]);
  
  fclose(f);
  info("2D Pattern spectrum file closed (%s)", filename);
  free(filename);
}

void write_pattern_spectra2d(Arguments *args, double* spectrum, LambdaVec *lvec_attr1, LambdaVec *lvec_attr2){
  char *filename;
  int numscales_attr1 = lvec_attr1->num_lambdas;
  int numscales_attr2 = lvec_attr2->num_lambdas;

  if(args->outprefix_orig != NULL)
    asprintf(&filename, "%s.txt", args->outprefix_arg);
  else
    asprintf(&filename, "pattern2d.txt");
	       
  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }
  fprintf(f, "#lambdas_attr1 \t #lambda_attr2 \t spectra \n");

  for (ulong i = 0; i < (ulong) numscales_attr1*numscales_attr2; i++) 
    //fprintf(f, "Spec[%ld]= %lf \n", i, spectrum[i]);
    fprintf(f, "%f \t %f \t %lf \n", lvec_attr1->lambdas[i%numscales_attr1], lvec_attr2->lambdas[i/numscales_attr1],  spectrum[i]);
  
  fclose(f);
  info("Pattern spectrum file closed (%s)", filename);
  free(filename);
}

void write_tree_file_txt(Arguments *args, Node *tree,ulong *dims ) {
  char *filename;
  ulong width  = dims[0];
  ulong height = dims[1];
  ulong depth  = dims[2];
  ulong size2D = width * height;

  
  if(args->outprefix_orig != NULL)
    asprintf(&filename, "%s.txt", args->outprefix_arg);
  else
    asprintf(&filename, "ComponentTree.txt");

  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }
  
  fprintf(f, "#Pixel list, -1 in column Area means pixel does NOT represent a connected component\n");
  fprintf(f, "#index \t x \t y \t z \t gval \t parent \t flux \t Area|Volume \t");
  
  if(args->attribute_arg == 1)
    fprintf(f, "Area of Min Enclosing rectangle");
  else if(args->attribute_arg == 2)
    fprintf(f, "Diag square of Min Enclosing rectangle");
  else if(args->attribute_arg >= 3 && args->attribute_arg <= 7 )
    fprintf(f, "Inertia tensor trace \t Inertia tensor trace / area ^2 \t Mean X \t Mean Y \t Mean Z");
  else if( args->attribute_arg > 7 )
    fprintf(f, "Elongation \t Flatness \t Sparseness \t Non-compactness \t cm_x \t cm_y \t cm_z \t cm_wx \t cm_wy \t cm_wz");

   fprintf(f, "\n");

  for (ulong i = 0; i < tree->size_curr; i++) {
    long size = -1;
    float dh = 0;
    ulong x, y, z;
    x = (i % size2D) % width + tree->offsets[0];
    y = (i % size2D) / width +  tree->offsets[1];
    z = i / size2D +  tree->offsets[2];
    
    if(is_levelroot(tree,i) ){
      size = (long) (*AttribsArray[args->attribute_arg].area)(tree->attribute + i*tree->size_attr);
      dh = tree->parent[i] != BOTTOM ? (tree->gval[i] - tree->gval[tree->parent[i]]): tree->gval[i];
    } 

    fprintf(f, "%ld \t %ld \t %ld \t %ld \t %f \t %ld \t %f \t %ld \t", i, x, y, z, tree->gval[i], get_levelroot(tree, tree->parent[i]), dh * (float) size, size);
    if(args->attribute_arg == 1){
      if(is_levelroot(tree,i) ){
	float encRec = (float) (*AttribsArray[args->attribute_arg].attribute)(tree->attribute + i*tree->size_attr);
	fprintf(f, "%f", encRec);
      } else {
	fprintf(f, "%f", -1);
      }
    } else if(args->attribute_arg == 2){
      if(is_levelroot(tree,i) ){
	float diagRec = (float) (*AttribsArray[args->attribute_arg].attribute)(tree->attribute + i*tree->size_attr);
	fprintf(f, "%f", diagRec);
      } else {
	fprintf(f, "%f", -1);
      }
    }  else if(args->attribute_arg >= 3 && args->attribute_arg <= 7 ){
      if(is_levelroot(tree,i) ){
	float in1 = (float) (*AttribsArray[3].attribute)(tree->attribute + i*tree->size_attr);
	float in2 = (float) (*AttribsArray[4].attribute)(tree->attribute + i*tree->size_attr);
	float meanx = (float) (*AttribsArray[5].attribute)(tree->attribute + i*tree->size_attr);
	float meany = (float) (*AttribsArray[6].attribute)(tree->attribute + i*tree->size_attr);
	float meanz = (float) (*AttribsArray[7].attribute)(tree->attribute + i*tree->size_attr);

	fprintf(f, "%f \t %f \t %f \t %f \t %f", in1, in2, meanx, meany, meanz);
      } else {
	fprintf(f, "%f \t %f \t %f \t %f \t %f", -1, -1, -1, -1, -1);
      }
    } else if(args->attribute_arg > 7){
      if(is_levelroot(tree,i) ){
	double *attrArr = (double*) inertiafull_attribute_arr(tree->attribute + i*tree->size_attr);
	float elong = (float) attrArr[8];//(*AttribsArray[9].attribute)(tree->attribute + i*tree->size_attr);
	float flat  = (float) attrArr[9]; //(*AttribsArray[10].attribute)(tree->attribute + i*tree->size_attr);
	float spars = (float) attrArr[10];//(*AttribsArray[11].attribute)(tree->attribute + i*tree->size_attr);
	float ncomp = (float) attrArr[11];//(*AttribsArray[12].attribute)(tree->attribute + i*tree->size_attr);
	float c_x   = (float) attrArr[2];
	float c_y   = (float) attrArr[3];
	float c_z   = (float) attrArr[4];
	float c_xw   = (float) attrArr[5];
	float c_yw   = (float) attrArr[6];
	float c_zw   = (float) attrArr[7];	
	fprintf(f, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f ", elong, flat, spars, ncomp, c_x, c_y, c_z, c_xw, c_yw, c_zw);
	
      } else {
	fprintf(f, "%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f ", -1, -1, -1, -1,-1, -1, -1,-1, -1, -1);
      }
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");

  fclose(f);
  info("Component tree file closed (%s)", filename);
  free(filename);
} /* write_maxtree_file_binary */


FIBITMAP* freeimage_generic_loader(const char* lpszPathName, int flag) {
  FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
  //  check the file signature and deduce its format
  fif = FreeImage_GetFileType(lpszPathName, 0);
  if(fif == FIF_UNKNOWN) {
    // no signature ? try to guess the file format from the file extension
    fif = FreeImage_GetFIFFromFilename(lpszPathName);
  }
  // check that the plugin has reading capabilities ...
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
    FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
    // unless a bad file format, we are done !
    return dib;
  }
  else{
    error("FreeImage couldn't read the input \n");
    MPI_Abort(MPI_COMM_WORLD, 735);
  }
  return NULL;
} /* freeimage_generic_loader */

const char *get_filename_ext(const char *filename) {
  const char *dot = strrchr(filename, '.');
  if(!dot || dot == filename) return "";
  return dot + 1;
}

void get_distributed_borders(value **img, int grid[3],  ulong dims[3], ulong dims_T[3], bool border[6], int overlap){

  value 	*img_curr ;   
  int 		myrank	      = rank();
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};


  MPI_Comm row_comm, col_comm, slice_comm;
  MPI_Comm_split(MPI_COMM_WORLD, myrank_arr[1]+myrank_arr[2]*grid[1], myrank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_arr[0]+myrank_arr[2]*grid[0], myrank, &col_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_2D, myrank, &slice_comm);
  if (!overlap){
    img_curr     = malloc(dims[0]*dims[1]*dims[2]*sizeof(value));
    memcpy(img_curr, *img, dims[0]*dims[1]*dims[2]*sizeof(value));

    ulong new_dims[3] = {dims[0]+border[0]+border[1], dims[1]+border[2]+border[3], dims[2]+border[4]+border[5]};
    free(*img);
    *img = malloc(new_dims[0]*new_dims[1]*new_dims[2] * sizeof(value));
    MPI_Datatype gdep_s, gdep_r, gcol_s, gcol_r, gpla_s, gpla_r, grow_s, grow_r, gang_s, gang_r;
    MPI_Type_vector(dims[1], dims[0], dims[0], mpi_value_type, &gdep_s);
    MPI_Type_commit(&gdep_s);
    MPI_Type_vector(dims[1], dims[0], new_dims[0], mpi_value_type, &gdep_r);
    MPI_Type_commit(&gdep_r);
    MPI_Type_vector(dims[1], 1, dims[0], mpi_value_type, &gcol_s);
    MPI_Type_commit(&gcol_s);
    MPI_Type_vector(dims[1], 1, new_dims[0], mpi_value_type, &gcol_r);
    MPI_Type_commit(&gcol_r);
    MPI_Type_hvector(dims[2], 1, dims[0]*dims[1]*sizeof(value), gcol_s, &gpla_s);
    MPI_Type_commit(&gpla_s);
    MPI_Type_hvector(dims[2], 1, new_dims[0]*new_dims[1]*sizeof(value), gcol_r, &gpla_r);
    MPI_Type_commit(&gpla_r);
    MPI_Type_vector(dims[2], dims[0], dims[0]*dims[1], mpi_value_type, &grow_s);
    MPI_Type_commit(&grow_s);
    MPI_Type_vector(dims[2], dims[0], new_dims[0]*new_dims[1], mpi_value_type, &grow_r);
    MPI_Type_commit(&grow_r);
    MPI_Type_vector(dims[2], 1, dims[0]*dims[1], mpi_value_type, &gang_s);
    MPI_Type_commit(&gang_s);
    MPI_Type_vector(dims[2], 1, new_dims[0]*new_dims[1], mpi_value_type, &gang_r);
    MPI_Type_commit(&gang_r);
    int row_rank, col_rank, dep_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_rank(slice_comm, &dep_rank);

    if(rank() == 2)
      info("ROW %d, Col %d, dep %d", row_rank, col_rank, dep_rank);
    /* Send planes*/
    if(row_rank % 2 == 0){
      if(border[1]){
	MPI_Send(img_curr+dims[0]-1, 1, gpla_s, row_rank+1, 1, row_comm);
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+(1+border[2])*new_dims[0]-1, 1, gpla_r, row_rank+1, 1, row_comm, NULL);
      }
      if(border[0]){
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+border[2]*new_dims[0], 1, gpla_r, row_rank-1, 1, row_comm, NULL);
	MPI_Send(img_curr,1, gpla_s, row_rank-1,1,  row_comm);
      }
    } else {
      if(border[0]){
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+border[2]*new_dims[0], 1, gpla_r, row_rank-1, 1, row_comm, NULL);
	MPI_Send(img_curr,1, gpla_s, row_rank-1,1,  row_comm);
      }
      if(border[1]){
	MPI_Send(img_curr+dims[0]-1, 1, gpla_s, row_rank+1, 1, row_comm);
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+(1+border[2])*new_dims[0]-1, 1, gpla_r, row_rank+1, 1, row_comm, NULL);
      }
    }

    if(col_rank % 2 == 0){
      if(border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, grow_s, col_rank+1, 1, col_comm);
	MPI_Recv(*img+(border[4]+1)*new_dims[0]*new_dims[1]-new_dims[0]+border[0], 1, grow_r, col_rank+1, 1, col_comm, NULL);
      }
      if(border[2]){
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+border[0],1, grow_r, col_rank-1,1,  col_comm, NULL);
	MPI_Send(img_curr,1, grow_s, col_rank-1,1,  col_comm);
      }
    } else {
      if(border[2]){
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+border[0],1, grow_r, col_rank-1,1,  col_comm, NULL);
	MPI_Send(img_curr,1, grow_s, col_rank-1,1,  col_comm);

      }
      if(border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, grow_s, col_rank+1, 1, col_comm);
	MPI_Recv(*img+(border[4]+1)*new_dims[0]*new_dims[1]-new_dims[0]+border[0], 1, grow_r, col_rank+1, 1, col_comm, NULL);
      }
    }

    if(dep_rank %2 == 0){
      if(border[5]){
	MPI_Send(img_curr+(dims[2]-1)*dims[1]*dims[0], 1, gdep_s, dep_rank+1, 1, slice_comm);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]*border[2]+border[0], 1, gdep_r, dep_rank+1, 1, slice_comm, NULL);
      }
      if(border[4]){
	MPI_Recv(*img+new_dims[0]*border[2]+border[0],1, gdep_r, dep_rank-1,1,  slice_comm, NULL);
	MPI_Send(img_curr,1,gdep_s, dep_rank-1,1,  slice_comm);
      }
    } else {
      if(border[4]){
	MPI_Recv(*img+new_dims[0]*border[2]+border[0],1, gdep_r, dep_rank-1,1,  slice_comm, NULL);
	MPI_Send(img_curr,1,gdep_s, dep_rank-1,1,  slice_comm);
      }
      if(border[5]){
	MPI_Send(img_curr+(dims[2]-1)*dims[1]*dims[0], 1, gdep_s, dep_rank+1, 1, slice_comm);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]*border[2]+border[0], 1, gdep_r, dep_rank+1, 1, slice_comm, NULL);
      }
    }

    if(row_rank % 2 == 0){
      if(border[1] && border[2]){
	MPI_Send(img_curr+dims[0]-1,1, gang_s, rank()+1-grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+new_dims[0]-1,1, gang_r, rank()+1-grid[0],1,  MPI_COMM_WORLD, NULL);
      }
      if(border[1] && border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]-1,1, gang_s, rank()+1+grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(1+border[4])*new_dims[0]*new_dims[1]-1,1, gang_r, rank()+1+grid[0],1,  MPI_COMM_WORLD, NULL);
      }
      if(border[0] && border[2]){
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1],1, gang_r, rank()-1-grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr,1, gang_s, rank()-1-grid[0],1,  MPI_COMM_WORLD);
      }
      if(border[0] && border[3]){
	MPI_Recv(*img+(1+border[4])*new_dims[0]*new_dims[1]-new_dims[0],1, gang_r, rank()-1+grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+(dims[1]-1)*dims[0],1, gang_s, rank()-1+grid[0],1,  MPI_COMM_WORLD);
      }
    } else {
      if(border[0] && border[3]){
	MPI_Recv(*img+(1+border[4])*new_dims[0]*new_dims[1]-new_dims[0],1, gang_r, rank()-1+grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+(dims[1]-1)*dims[0],1, gang_s, rank()-1+grid[0],1,  MPI_COMM_WORLD);
      }
      if(border[0] && border[2]){
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1],1, gang_r, rank()-1-grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr,1, gang_s, rank()-1-grid[0],1,  MPI_COMM_WORLD);
      }
      if(border[1] && border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]-1,1, gang_s, rank()+1+grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(1+border[4])*new_dims[0]*new_dims[1]-1,1, gang_r, rank()+1+grid[0],1,  MPI_COMM_WORLD, NULL);
      }
      if(border[1] && border[2]){
	MPI_Send(img_curr+dims[0]-1,1, gang_s, rank()+1-grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+border[4]*new_dims[0]*new_dims[1]+new_dims[0]-1,1, gang_r, rank()+1-grid[0],1,  MPI_COMM_WORLD, NULL);
      }
    }
      
    if(dep_rank %2 == 0){
      if(border[4] && border[2]){
	MPI_Send(img_curr,dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+border[0],dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
      }
      if(border[4] && border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1),dims[0], mpi_value_type, rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(new_dims[1]-1)*new_dims[0]+border[0],dims[0], mpi_value_type,  rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
      }		
      if(border[5] && border[2]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+border[0],dims[0], mpi_value_type,rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1),dims[0], mpi_value_type,  rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
      }
      if(border[5] && border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2]-new_dims[0]+border[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*dims[2]-dims[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
      }

    } else {
      if(border[5] && border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2]-new_dims[0]+border[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*dims[2]-dims[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
      }	
      if(border[5] && border[2]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+border[0],dims[0], mpi_value_type,rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1),dims[0], mpi_value_type,  rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
      }
      if(border[4] && border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1),dims[0], mpi_value_type, rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(new_dims[1]-1)*new_dims[0]+border[0],dims[0], mpi_value_type,  rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
      }
      if(border[4] && border[2]){
	MPI_Send(img_curr,dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+border[0],dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
      }
    }

    if(dep_rank %2 == 0){
      if(border[4] && border[0]){
	MPI_Send(img_curr,1, gcol_s, rank()-grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+border[2]*new_dims[0],1, gcol_r, rank()-grid[0]*grid[1] -1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[4] && border[1]){
	MPI_Send(img_curr+dims[0]-1,1, gcol_s, rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(1+border[2])*new_dims[0]-1,1, gcol_r,  rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[5] && border[0]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, gcol_s,  rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+border[2]*new_dims[0],1, gcol_r,rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[5] && border[1]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+(border[2]+1)*new_dims[0]-1,1, gcol_r, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
      }
    } else {
      if(border[5] && border[1]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+(border[2]+1)*new_dims[0]-1,1, gcol_r, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
      }
      if(border[5] && border[0]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+border[2]*new_dims[0],1, gcol_r,rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, gcol_s,  rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
      }
      if(border[4] && border[1]){
	MPI_Recv(*img+(1+border[2])*new_dims[0]-1,1, gcol_r,  rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]-1,1, gcol_s, rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
      }
      if(border[4] && border[0]){
	MPI_Recv(*img+border[2]*new_dims[0],1, gcol_r, rank()-grid[0]*grid[1] -1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr,1, gcol_s, rank()-grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
      }
    }

    if(dep_rank %2 == 0){
      if(border[4] && border[0] && border[2]){
	MPI_Send(img_curr, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[4] && border[1] && border[2]){
	MPI_Send(img_curr+dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[4] && border[0] && border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*(new_dims[1]-1),1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[4] && border[1] && border[3]){
	MPI_Send(img_curr+dims[1]*dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }
	
      if(border[5] && border[0] && border[2]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1),1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[5] && border[1] && border[2]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }

      if(border[5] && border[0] && border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2] - new_dims[0],1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(border[5] && border[1] && border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[1]*dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }
    } else{
      if(border[5] && border[1] && border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[1]*dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(border[5] && border[0] && border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2] - new_dims[0],1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
      }
      if(border[5] && border[1] && border[2]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(border[5] && border[0] && border[2]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1),1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
      }
      if(border[4] && border[1] && border[3]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[1]*dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(border[4] && border[0] && border[3]){
	MPI_Recv(*img+new_dims[0]*(new_dims[1]-1),1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
      }
      if(border[4] && border[1] && border[2]){
	MPI_Recv(*img+new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(border[4] && border[0] && border[2]){
	MPI_Recv(*img,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
      }
    }

    ulong offset = border[4]*new_dims[0]*new_dims[1] + border[2]*new_dims[0]+ border[0];
    
    for(ulong i = 0; i < dims[0]*dims[1]*dims[2]; i++){
      ulong x = i % dims[0];
      ulong y = i %(dims[0]*dims[1]) / dims[0];
      ulong z = i / (dims[0]*dims[1]);
      (*img)[offset+z*(new_dims[0]*new_dims[1])+y*new_dims[0]+x] = img_curr[i];
    }

    dims[0] = new_dims[0];
    dims[1] = new_dims[1];
    dims[2] = new_dims[2];
    free(img_curr);
    MPI_Type_free( &grow_r );
    MPI_Type_free( &grow_s );
    MPI_Type_free( &gcol_r );
    MPI_Type_free( &gcol_s );
    MPI_Type_free( &gpla_r );
    MPI_Type_free( &gpla_s );
    MPI_Type_free( &gdep_r );
    MPI_Type_free( &gdep_s );

  }

   
  ulong send = dims[0]-border[0]-border[1];
  MPI_Allreduce(&send, dims_T, 1, MPI_UNSIGNED_LONG, MPI_SUM, row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_arr[0]+myrank_arr[2]*grid[0], myrank, &col_comm);
  send = dims[1]-border[2]-border[3];
  MPI_Allreduce(&send, dims_T+1, 1, MPI_UNSIGNED_LONG, MPI_SUM, col_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_2D, myrank, &slice_comm);
  send = dims[2]-border[4]-border[5];
  MPI_Allreduce(&send, dims_T+2, 1, MPI_UNSIGNED_LONG, MPI_SUM, slice_comm);
  MPI_Comm_free(&row_comm);
  MPI_Comm_free(&col_comm);
  MPI_Comm_free(&slice_comm);

}

