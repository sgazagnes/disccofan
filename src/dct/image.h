#ifndef IMAGE_H
#define IMAGE_H

/* Read */
value *read_input(Arguments *args, char *prefix, ulong dims[3], ulong dims_T[3], bool *border);
value *read_fits(Arguments *args, const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix);
value *read_hdf5(Arguments *args, const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix);
value *read_basic(Arguments *args, const char *fname, ulong dims[3], ulong dims_T[3], int *bitpix);
value *read_nifti_file(Arguments *args, const char *fname, ulong dims[3], ulong dims_T[3],  int *bitpix);
value *read_raw_file(Arguments *args, const char *fname, ulong dims[3], ulong dims_T[3],  int *bitpix);
/* Write */
void write_output(Arguments *args, value *img, const char *attr_name, ulong dims_T[3], ulong dims[3],  bool border[6], int bitpix);
void write_differential(Arguments *args, value *outOrig, value *outDH, value *outScale, const char *attr_name, ulong dims_T[3], ulong dims[3],  bool border[6], int bitpix);
void write_hdf5(Arguments *args, const char* fname, const char *dataset_out, value *out, ulong dims_T[3], ulong dims[3], bool border[6], int bitpix);
void write_fits(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3], bool border[6], int bitpix);
void write_basic(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3],  bool border[6], int bitpix);
void write_pattern_spectra(Arguments *args, double* spectrum, int numscales);

#endif
