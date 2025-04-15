#ifndef REFINE_H
#define REFINE_H

#define BG_REJECT_TILE 0
#define BG_ACCEPT_TILE 1

#define BG_REJECTION_RATE 0.05

#define BG_TILE_SIZE_START 64
#define BG_TILE_SIZE_MIN 16
#define BG_TILE_SIZE_MAX 128


double *refine(Arguments *args, value *gvals, double *extrema, ulong size);
void segment(Arguments *args, value *gvals, ulong size);

#endif
