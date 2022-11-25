#ifndef REFINE_H
#define REFINE_H

float* find_extremas(value *gvals, ulong size);
void apply_log(value *gvals, ulong size );
void apply_log10(value *gvals,  ulong size);
ulong *histogram_tile(value *gvals, int samples, float maxval, float minval, ulong size);
void apply_refine(value *gvals, ulong size, float *factors);
void check_refine(Arguments *args, value *gvals, ulong size);

#endif
