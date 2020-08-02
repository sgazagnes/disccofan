#ifndef IMAGE_H
#define IMAGE_H

/* Read */
void read_input(Arguments *args, value **img, ulong dims[3], ulong dims_T[3], ulong offsets[3], bool border[6]);
/* Write */
void write_output(Arguments *args, value *img, const char *attr_name, ulong dims_T[3], ulong dims[3],  bool border[6]);
void write_differential(Arguments *args, value *outOrig, value *outDH, value *outScale, const char *attr_name, ulong dims_T[3], ulong dims[3],  bool border[6]);
void write_pattern_spectra(Arguments *args, double* spectrum, int numscales);

#endif
