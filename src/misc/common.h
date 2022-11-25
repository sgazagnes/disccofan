#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <FreeImage.h>
#include <time.h>
#include <fitsio.h>
#include <hdf5.h>
#include <float.h> 
#include <ctype.h>
#include <time.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

typedef uint8_t ubyte;
typedef uint16_t ushort;
typedef uint32_t uint;
typedef uint64_t ulong;
typedef int64_t idx;
typedef uint value_t;

#define FLOAT_TYPE 0
#define DOUBLE_TYPE 1

typedef double value;

