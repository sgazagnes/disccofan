#include "io.h"

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
    printf("FreeImage error");
    exit(0);
  }
  return NULL;
} /* freeimage_generic_loader */


double *read_fits_file(const char *fname,  ulong *dims,  int *bitpix){
  /* +++++++++++++++++++++++++++ */
  /*          Read Fits          */
  /* +++++++++++++++++++++++++++ */
  
  fitsfile 	*fptr;   		/* FITS file pointer */
  int 		status        = 0;   	/* CFITSIO status value MUST be initialized to zero! */
  int		n_dims 	      = 0;
  long 		inc[4]        = {1,1,1,1};	/* Increment for fits read function */
  long 		naxes[4]      = {1,1,1,1};	/* Dimensions (4D max) */
  long 		counts[4]     = {1,1,1,1};	/* Pixels to read in fits function */
  long 		offsets[4]    = {1,1,1,1};	/* Pixel index to start in fist function */

  double 	*img;
  
  if (!fits_open_file(&fptr, fname, READONLY, &status)) {	/* Open file */
    if (!fits_get_img_param(fptr, 4,  bitpix,  &n_dims, naxes, &status)) {  /* Get image parameter */
      if (n_dims > 4 || n_dims == 0 || (n_dims == 4 && naxes[3] > 1)) { /* 3D Max */
	printf("Error reading fits \n");
	exit(0);
      } else {
	if(n_dims == 4) (n_dims)--;
	
	for (int i = 0; i < n_dims; i++) 
	  counts[i]  = dims[i] = naxes[i];
	
	img  =  malloc(dims[0]*dims[1]*dims[2] * sizeof(double));  
	fits_read_subset(fptr, TDOUBLE , offsets, counts, inc, NULL, img, NULL, &status);	
      }
    }
    fits_close_file(fptr, &status);
  } else {
    printf("Error reading fits \n");
    exit(0);
  }    
  return img;
}


double *read_hdf5_file(const char *fname, const char *dataset,  ulong *dims, int *bitpix){
  
  /* +++++++++++++++++++++++++++ */
  /*      Read file HDF5        */
  /* +++++++++++++++++++++++++++ */
  
  H5T_class_t   t_class;
  hid_t       	file_id, dataset_id, dataspace;  		/* identifiers */
  hsize_t 	counts[3]     = {1,1,1};
  hsize_t 	offsets[3]    = {0};
  hsize_t 	*hdims;
  hsize_t 	hn_dims;
  herr_t	err;
    
  /*      Open HDF5 file and dataset  */
  
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);	
  if (file_id < 0) {
    printf("Could not open file %s, wrong name ? (HDF5) \n", fname);
    exit(0);
  }
  
  dataset_id = H5Dopen2(file_id,dataset, H5P_DEFAULT);	
  if (dataset_id < 0) {
    printf("Could not open dataset %s, wrong name ? (HDF5) ", dataset);
    exit(0);
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
    else if(ord == 32)
      *bitpix = 32;
    else
      *bitpix = 64;
  }
  
  hn_dims = H5Sget_simple_extent_ndims(dataspace);
  
  /*       Get number of dimensions  and dimensions    */
  
  if (hn_dims > 3) {
    printf("Only handle 2D or 3D data \n");
    exit(0);
  }
  
  hdims = calloc(hn_dims, sizeof(hsize_t));			
  H5Sget_simple_extent_dims(dataspace, hdims, NULL);

  for (hsize_t i = 0; i < hn_dims; i++)
    counts[i] = dims[hn_dims - i - 1] = hdims[i];
  
  double *img = malloc(counts[0]*counts[1]*counts[2] * sizeof(double));  
  
  hid_t memory_window = H5Screate_simple(hn_dims, counts, NULL); /* Allocate memory to read the data */
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL); 
  err = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memory_window, dataspace, H5P_DEFAULT, img);
 
  if (err < 0){ 
    printf("Could not read data from the dataset (HDF5) \n");
    exit(0);
  }

  free(hdims);
  
  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);
  return img;
}


double *read_basic_file(const char *fname, ulong *dims, int *bitpix) {
  /* +++++++++++++++++++++++++++ */
  /*       Read FreeImage        */
  /* +++++++++++++++++++++++++++ */
  
  // printf("Reading image %s with FreeImage \n", fname);	

  ulong 	i    = 0; 
  FIBITMAP 	*dib = freeimage_generic_loader(fname, 0);  /* Generic Loader from FreeImage */
  
  if (dib == NULL) {
    printf("Can't read the image %s (FreeImage) \n", fname);
    exit(0);
  }
  
  *bitpix = FreeImage_GetBPP(dib);				/* Get dynamic range */
  dims[0] = FreeImage_GetWidth(dib);				/* Get width */
  dims[1] = FreeImage_GetHeight(dib);				/* Get height */
  FREE_IMAGE_COLOR_TYPE nc = FreeImage_GetColorType(dib);

  /* if(nc != 1){
     FIBITMAP* hImage= dib;
     dib = FreeImage_ConvertTo32Bits( hImage );
     hImage= dib;
     dib = FreeImage_ConvertTo8Bits( hImage );
     FreeImage_Unload( hImage );
     }*/
  if(nc != 1)
    dib =  FreeImage_ConvertToGreyscale(dib);

  *bitpix = FreeImage_GetBPP(dib);				/* Get dynamic range */

  double *img = malloc(dims[0]*dims[1] * sizeof(double));	
   
  if (*bitpix <= 8) {
    for(ulong y = 0; y < dims[1]; y++) {
      BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, y);
      for(ulong x = 0; x < dims[0]; x++,i++){
	img[i] = (double) bits[x];
      }
    }
    FreeImage_Unload(dib);
    return img;
    
  } else if (*bitpix <= 16) {
    for(ulong y = 0; y < dims[1]; y++) {
      ushort *bits = (ushort *)FreeImage_GetScanLine(dib, y);
      for(ulong x = 0; x < dims[0]; x++,i++) 
	img[i] = (double) bits[x];
    }
    FreeImage_Unload(dib);
    return img;  
  } else if (*bitpix <= 32) {
    for(ulong y = 0; y < dims[1]; y++) {
      uint *bits = (uint *)FreeImage_GetScanLine(dib, y);
      for(ulong x = 0; x < dims[0]; x++,i++) 
	img[i] = (double) bits[x];
    }
    FreeImage_Unload(dib);
    return img;  
  } else{
    printf("Floating point or 64 bits depth are not supported with classical images (Freeimage) \n");
    exit(0);
    FreeImage_Unload(dib);
    return NULL;
  }
}



void write_basic_file(const char *fname_out, double *im, ulong *dims, int bitpix){
  
  /* +++++++++++++++++++++++++++ */
  /*        Write FreeImage      */
  /* +++++++++++++++++++++++++++ */
  
  printf("Writing image %s using FreeImage \n", fname_out);	  

  if(dims[2] > 1){
    printf("FreeImage does not handle writing 3D images \n");
    exit(0);
  }

  FIBITMAP 	*outmap;
  long		j,y;
  ulong 	width  = dims[0];
  ulong 	height = dims[1];
 
  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname_out);
  if (fif == FIF_UNKNOWN) 
    fif = FreeImage_GetFIFFromFilename(fname_out);
  if ((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsWriting(fif)) {
    
    if (bitpix <= 8) {
      ubyte *imagebuf = malloc(width*height* sizeof(ubyte));;
      outmap = FreeImage_AllocateT(FIT_BITMAP, width, height, 8, 0xFF, 0xFF, 0xFF);
      for (j=height-1,  y=0; j>=0; j--, y++) {
	imagebuf = FreeImage_GetScanLine(outmap,j);      
	for (ulong x=0;x<width;x++)
	  imagebuf[x]=im[width*j + x];
      }	 
    } else if (bitpix <= 16) {
      ushort *imagebuf;
      outmap = FreeImage_AllocateT(FIT_UINT16,width,height,16,0xFFFF,0xFFFF,0xFFFF);
      for (j=height-1, y=0; j>=0; j--, y++){      
	imagebuf = (ushort *)FreeImage_GetScanLine(outmap,j);
	for (ulong x=0; x<width ;x++)
	  imagebuf[x]=im[(width)*j + x];	
      }
    } else if (bitpix <= 32) {
      unsigned int *imagebuf;
      outmap = FreeImage_AllocateT(FIT_UINT32,width,height,32,0xFFFF,0xFFFF,0xFFFF);
      for (j=height-1, y=0; j>=0; j--, y++){      
	imagebuf = (uint *)FreeImage_GetScanLine(outmap,j);
	for (ulong x=0; x<width ;x++)
	  imagebuf[x]=im[(width)*j + x];	
      }
    }else{
      printf("not handling more than 32, 64 bits or floating point\n");
      exit(0);
    }       
    FreeImage_Save(fif,outmap,fname_out,0); 
    FreeImage_Unload(outmap);
  }
  else{
    printf("FreeImage couldn't write the output\n");
    exit(0);
  }
}


void write_hdf5_file(const char* fname_out, const char *dataset_out, double *out, ulong *dims, int bitpix, int n_dims){
  
  /* +++++++++++++++++++++++++++ */
  /*     HDF5 Write Function     */
  /* +++++++++++++++++++++++++++ */
  
  hid_t       	file_id, dataset_id, dataspace, type;  /* identifiers */
  hsize_t	hdims[3]      = {0};
  hsize_t 	offsets[3]    = {0,0,0};

  if(n_dims > 2){
    hdims[0] = dims[2];hdims[1] = dims[1]; hdims[2] = dims[0];
  } else {
    hdims[0] = dims[1];hdims[1] = dims[0];
  }

  if(bitpix < 0)
    type = H5T_NATIVE_FLOAT;
  else if(bitpix <= 8)
    type =  H5T_NATIVE_UCHAR;
  else if(bitpix <= 16)
    type =  H5T_NATIVE_USHORT;
  else if(bitpix <= 32)
    type =  H5T_NATIVE_UINT;
  else{
    type =  H5T_NATIVE_ULONG;
    }

  file_id = H5Fcreate(fname_out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
  hid_t global_memory_space = H5Screate_simple(n_dims, hdims, NULL);
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  dataset_id = H5Dopen(file_id,dataset_out, H5P_DEFAULT);		/* Open dataset */
  if (dataset_id < 0) {
    dataset_id = H5Dcreate(file_id,dataset_out, type, global_memory_space,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	/* Create dataset */
  }
  dataspace = H5Dget_space(dataset_id);			/* Copy dataset */
  type =  H5T_NATIVE_DOUBLE; 
  
  hid_t memory_window = H5Screate_simple((hsize_t)n_dims, hdims, NULL);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, hdims, NULL);
  H5Dwrite(dataset_id, type, memory_window, dataspace, H5P_DEFAULT, out);

  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);
}



void write_fits_file(char *fnameout,  double *out, ulong *dims, int bitpix, int n_dims){
  /* +++++++++++++++++++++++++++ */
  /*     FITS Write Function     */
  /* +++++++++++++++++++++++++++ */
  
  fitsfile 	*outfptr, *infptr; 	/* FITS file pointers */
  int 		status = 0; 		/* CFITSIO status value MUST be initialized to zero! */
  int 		type;
  long 		naxes[3]      = {dims[0], dims[1], dims[2]};      
  long 		counts[3]     = {dims[0], dims[1], dims[2]};
  long 		offsets[3]    = {1,1,1};
  char 		str1[100];
  int           negative = 0;
  strcpy(str1, "!");   		/* '!' symbol makes the output image, if existing, to be overwritten */


  
  //  if(FLOAT_TYPE == 1)
    type = TDOUBLE;
    /*  else if (sizeof(value) == 1)					
    type = TBYTE;
  else if (sizeof(value) == 2)
    type = TUSHORT;
  else if (sizeof(value) == 4)
    type = TUINT;
  else {
    printf("Do not handle > 32 bits \n");
    exit(0);
  }*/

  if (!fits_create_file(&outfptr, strcat(str1, fnameout), &status)) {
    if(bitpix < -32)
      fits_create_img(outfptr, DOUBLE_IMG, n_dims, dims, &status);
    else if(bitpix < 0)
      fits_create_img(outfptr, FLOAT_IMG, n_dims, dims, &status);
    else if(bitpix <=8 && negative == 0)
      fits_create_img(outfptr, BYTE_IMG, n_dims, dims, &status);
    else if(bitpix <=8 && negative == 1)
      fits_create_img(outfptr, SBYTE_IMG, n_dims, dims, &status);
    else if(bitpix <= 16 && negative == 1)
      fits_create_img(outfptr, SHORT_IMG, n_dims, dims, &status);
    else if(bitpix <= 32 && negative == 1)
      fits_create_img(outfptr, LONG_IMG, n_dims, dims, &status);
    else if(bitpix <= 16 && negative == 0)
      fits_create_img(outfptr, USHORT_IMG, n_dims, dims, &status);
    else if(bitpix <= 32 && negative == 0)
      fits_create_img(outfptr, ULONG_IMG, n_dims, dims, &status);
    else
      fits_create_img(outfptr, LONGLONG_IMG, n_dims, dims, &status);

    fits_write_subset(outfptr, type, offsets, counts, out, &status);
  } else {
    printf("Cannot create the file \n");
    exit(0);
  }
 
  fits_close_file(outfptr, &status);

} /* write_fits_file */

double *read_file(char *prefix, char *suffix,  char *dataset, ulong *dims,  int *bitpix){
  char *fname;
  double *img;
  
  asprintf(&fname, "%s.%s",prefix, suffix);
  if(!strcmp(suffix, "fits"))
    img = read_fits_file(fname,  dims, bitpix);
  else if(!strcmp(suffix, "h5"))
    img = read_hdf5_file(fname,  dataset,  dims, bitpix);
  else
    img = read_basic_file(fname, dims, bitpix);

  return img;
}


void write_file(char *prefix, char *suffix,  char *dataset, ulong *dims,  double *img, int bitpix, ulong n_dims){
  char *fname;
  void *img_out;

  /*if(bitpix == 8){
    img_out = (char *) calloc(dims[0]*dims[1]*dims[2], sizeof(char));
    #pragma omp parallel for
    for(long i = 0; i < dims[0]*dims[1]*dims[2]; i++)
      img_out[i] = (char) img[i];
  } else  if(bitpix == 16){
    img_out = (unsigned short *) calloc(dims[0]*dims[1]*dims[2], sizeof(unsigned short));
    #pragma omp parallel for
    for(long i = 0; i < dims[0]*dims[1]*dims[2]; i++)
      img_out[i] = (unsigned short) img[i];
  } else  if(bitpix == 32){
    img_out = (unsigned int *) calloc(dims[0]*dims[1]*dims[2], sizeof(unsigned int));
    #pragma omp parallel for
    for(long i = 0; i < dims[0]*dims[1]*dims[2]; i++)
      img_out[i] = (unsigned int) img[i];
  } else  if(bitpix == 64){
    img_out = (unsigned long *) calloc(dims[0]*dims[1]*dims[2], sizeof(unsigned long));
    #pragma omp parallel for
    for(long i = 0; i < dims[0]*dims[1]*dims[2]; i++)
      img_out[i] = (unsigned long) img[i];
  } else if(bitpix == -32){
    img_out = (float *) calloc(dims[0]*dims[1]*dims[2], sizeof(float));
    #pragma omp parallel for
    for(long i = 0; i < dims[0]*dims[1]*dims[2]; i++)
      img_out[i] = (float) img[i];
  } else
  img_out = img;*/

  asprintf(&fname, "%s.%s",prefix, suffix);
  if(!strcmp(suffix, "fits"))
    write_fits_file(fname, img, dims, bitpix, n_dims);
  else if(!strcmp(suffix, "h5"))
    write_hdf5_file(fname, dataset,  img, dims, bitpix, n_dims);
  else{
    if((bitpix > 16 || bitpix < 0))
      printf("Warn, 16 bits and floating point might not be handled correctly with FreeImage. \n ");
    write_basic_file(fname, img, dims, bitpix);
  }
}
