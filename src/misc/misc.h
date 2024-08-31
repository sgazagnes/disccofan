#include "common.h"

/*** Some usefull math macros ***/
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double mnarg1,mnarg2;
#define FMAX(a,b) (mnarg1=(a),mnarg2=(b),(mnarg1) > (mnarg2) ?\
(mnarg1) : (mnarg2))

static double mnarg1,mnarg2;
#define FMIN(a,b) (mnarg1=(a),mnarg2=(b),(mnarg1) < (mnarg2) ?\
(mnarg1) : (mnarg2))

#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))


/*********   BEGIN PROTOTYPE DEFINITIONS  ***********/

/*** Wrapper functions for the std library functions fwrite and fread
     which should improve stability on certain 64-bit operating systems
     when dealing with large (>4GB) files ***/
size_t mod_fwrite (const void *, unsigned long long, unsigned long long, FILE *);
size_t mod_fread(void *, unsigned long long, unsigned long long, FILE *);

/* generic function to compare floats */
int compare_floats(const void *, const void *);
double find_max(double *gvals, ulong size);
double find_min(double *gvals, ulong size);
double* find_extremas(double *gvals, ulong size);
void apply_log(double *gvals, ulong size);
void apply_log10(double *gvals, ulong size);
float as_float(uint32_t i);
double as_double(uint64_t i);
double randn (double mu, double sigma);
double *create_random(ulong size, int bitpix, int neg);
double *create_sequence(ulong size, int bitpix, int neg);
double *cross_corr(double *imgA, double *imgB, ulong size);
int compare_txt_files(FILE *fp1, FILE *fp2);
double *modify_bpp(double *gvals_in, ulong size, int bitpix, int nextbitpix);
double *create_bubbles_nov(ulong dim, ulong size,  float fract, int periodic);
double *create_bubbles_ov(ulong dim, ulong size, float fract, int periodic);
double *create_ioniz(ulong size, float fract, ulong dims[3]);
double keithRandom();
double *gaussian_filter(value *gvals, ulong *dims, double fhwm, ulong size );
int find_scale_double(double *thresh, double val, int nthresh);
double *hessian(value *gvals, ulong *dims, ulong p);
