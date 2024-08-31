#include "common.h"

void write_file(char *prefix, char *suffix,  char *dataset, ulong *dims,  double *img, int bitpix, ulong n_dims);
double *read_file(char *prefix, char *suffix,  char *dataset, ulong *dims,  int *bitpix);
