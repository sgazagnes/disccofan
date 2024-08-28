#ifndef TYPES_H
#define TYPES_H

#define _GNU_SOURCE 1

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <sys/types.h>
#include <sys/times.h>
#include <stdarg.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <libgen.h>
#include <ctype.h>
#include <assert.h>


#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))

/* +++++++++++++++++++++++++++++++ */
/*     	    Type Definition        */
/* +++++++++++++++++++++++++++++++ */

typedef uint8_t  ubyte;
typedef uint16_t ushort;
typedef uint32_t uint;
typedef uint64_t ulong;
typedef int64_t  idx;
typedef uint     value_t;

/* If not floating point : */

//#define FLOAT_TYPE 0
//typedef ubyte value; 			 // 8-bit (not working with flooding 1)
//typedef ushort value;			 // 16-bit
//typedef uint value; 			 // 32-bit
//typedef ulong value; 			 // 64-bit

/* If floating point:     */

#define FLOAT_TYPE 1
typedef float value;

#define MAX_OPERATION_SIZE 10 // Maximum size of an array parameter
#define MAX_PARAMS_SIZE 10 // Maximum size of an array parameter

/* +++++++++++++++++++++++++++++++ */
/*  	 Structure definition      */
/* +++++++++++++++++++++++++++++++ */

/* Lambda vector Structure */
typedef struct _LambdaVec{
  uint num_lambdas;
  float *lambdas;    
} LambdaVec;

/* Image Operation Structure */
typedef struct _Operation{
  char  name[256];	
  char lvec_name[256];	
  char lvec_name_2D[256];			
  int   attribute_idx;
  int   attribute_idx_2D;
  int   decision_idx;
  double lambda;
  double imscale;
  double imscale_2D;
  LambdaVec *lvec1;
  LambdaVec *lvec2;
} Operation;

typedef struct _DataParams{
    ulong dims_process[3];
    ulong dims_full[3];
    float pixdim[3];
    bool border[6];
    long offsets[3];
    value g_max_gval;
    value g_min_gval;
} DataParams;


/* Parsing structure */
typedef struct _Arguments {
    int interactive;
    char folder_run[256];
    char config_file[256];
    char input_name[256];
    char output_name[256];
    char input_prefix[256];
    char output_prefix[256];
    char input_type[256];
    char output_type[256];
    char hdf5_dataset[256];
    char image_options[256];
    char tree_type[256];
    char input_file_style[256];
    char output_file_style[256];
    char verbosity[256];
    char preprocessing[256];
    int tile_overlap;
    int mpi_grid[3];
    float pixel_dim[3];
    int threads;
    int bpp;
    int pixel_connectivity;
    int include_background;
    int save_output;
    int save_tree;
   // int morphology_choice;
    int num_operations;
    int image_stats;
    Operation ope_list[MAX_OPERATION_SIZE];
} Arguments;

/* Attribute storage Structure */
typedef struct _Allocator{
  size_t d_size;
  size_t d_pos;
  char* data;   /* Current number of attributes */
} Allocator;



/* Attribute storage Structure */
typedef struct _AuxDataStore{
  void  *data;		       		/* Attributes */ 
  size_t size_item;           		/* Size of attribute used */
  ulong size_alloc;			/* Allocation size */	
  ulong item_curr;			/* Current number of attributes */
} AuxDataStore;


/* Tree node Structure */
typedef struct _Node {
  idx 		*parent;
  idx 	 	*border_idx; 		    /* index in boundary tree */
  void 		*attribute;		      /* Pointer to attribute */
  double  *area_copy;		     /* Pointer to area copy (PS)*/
  value		*gval;			        /* Pixel intensity */
  value		*gval_par;			    /* Pixel intensity parent (PS) */
  ulong 	size_init;
  ulong		size_curr;
  ulong   size_attr;
  ulong   offsets[3];
  bool 		border[6];
} Node;

/* Output Structure */
typedef struct _Output {
  value   *filtered;
  value   *contrast;
  value   *scale;
  value   *luminance;
  double  *pattern_spec;
} Output;


/* Boundary tree node Structure */
typedef struct _BoundaryNode {
  value  	gval;			/* Pixel intensity */
  idx 		index;			/* Pixel index from local tree */
  idx 	 	border_idx; 		/* Index in boundary tree */
  //  BorderIndex 	border_par;
} BoundaryNode;

typedef struct _BorderIndex BorderIndex;

/* Boundary tree Structure */
typedef struct _Boundary {
  //  AuxDataStore  *store;			/* Attribute storage structure */
  BoundaryNode 	*array;
  void 		*attribute;		/* Pointer to attribute */
			/* Node array */
  
  BorderIndex 	*border_par;   		/* Parent in the boundary tree */
  BorderIndex 	*border_ori;   		/* Position in previous boundary tree */
  BorderIndex 	*border_lr;    		/* Levelroot at same intensity in other merged tree */
  ulong		*merge_idx;		/* Index of node to be merged */
  //  idx 		*attribute_idx;		/* Index of attribute in storage */
  bool	 	*reached;		/* Boolean array */
  ulong 	offset[7]; 		/* Number of nodes in different sides  */
  ulong		dims[3];		/* Dimensions of boundary tree */
  ulong 	size_init;		/* Initial size of tree */
  ulong 	size_curr;     		/* Current size of tree */
  ulong 	size_alloc;		/* Allocated size of tree */
  ulong 	size_attr;
} Boundary;

/* Boundary index Structure */
struct _BorderIndex {			
  Boundary 	*b;			/* Pointer to boundary tree */
  idx 		i;			/* Index in boundary */
};

/* Salembier queue Structure */
typedef struct _Queue {
  ulong *pixels;			/* Array of pixels */
  ulong head;		       		/* First pixel in list */
  ulong tail; 				/* Last pixel in list */
} Queue;

/* Wilkinson stack Structure */
typedef struct _pStack{
  ulong  pos_cur;			/* Current position in stack */
  ulong  size_max;			/* Size max of stack */
  ulong *array;  			/* Array of pixels */
} pStack;

/* Wilkinson queue Structure */
typedef struct _pQueue {
  ulong size_cur;			/* Current size */
  ulong size_max;			/* Max size */
  ulong *array;				/* Array of pixels */
} pQueue;

/* Teeninga sort item Structure */
#pragma pack(push, 1)
typedef struct _SortItem {
  value_t val;				/* Unsigned value */
  ulong rank;				/* Rank of value */
} SortItem;
#pragma pack(pop)

/* Teeninga queue Structure */
typedef struct _PrioQueue {
  ulong **m_levels;			/* Array of pixels by levels in the queue */
  ulong m_top;				/* Top pixel */
  int m_num_levels;			/* Number of levels in queue */
} PrioQueue;

/* Teeninga bitarray Structure */
typedef struct _BitArray {
  ulong *data;				/* Array */
  ulong num_words;			/* Bit word */
} BitArray;

/* Area Structure */
typedef struct _AreaData{
   ulong area;
} AreaData;

/* Extent Structure */
typedef struct _ExtentData {
   ulong area;
   ulong minX;
   ulong minY;
   ulong minZ;
   ulong maxX;
   ulong maxY;
   ulong maxZ;
} ExtentData;


/* Mean Structure */
typedef struct _MeanData {
  ulong area;
  ulong sumX, sumY, sumZ;
} MeanData;


/* Mean Structure */
typedef struct _WMeanData {
  ulong area;
  double sumGval;
  double sumXd, sumYd, sumZd;
} WMeanData;

typedef struct _InertiaData {
  ulong area;
  double sumval, sumval2;//, sumval2;//, sumvald, sumval2d;
  double sumX, sumY, sumZ, sumX2, sumY2, sumZ2, sumXY, sumYZ, sumXZ;
  double sumXd, sumYd, sumZd;//, sumX2d, sumY2d, sumZ2d, sumXYd, sumYZd, sumXZd;
} InertiaData;



/* Pruning decision Structure */
typedef struct DecisionStruct {
  const char *name;
  void (*filter)(Node *, value *, bool *, ulong, ulong, double (*attribute)(void *), double);
} DecisionStruct;

/* Attribute Structure */
typedef struct _AttribStruct{
  const char 	 *name;
  const char 	 *short_name;
  ulong  size;
  uint   group;
  void 	 *(*new_data)(AuxDataStore *);
  void   (*init_data)(void *, bool , ulong , ulong,  ulong, value);
  void   (*delete_data)(void *);
  void   (*add_to_data)(void *, void*);
  void   (*merge_data)(void *, void *);
  void   (*merge_to_data)( AuxDataStore *, void **, void *, void *);
  void   (*clone_data)( AuxDataStore *, void **, void *);
  void   (*create_mpi_data)(void);
  double (*area)(void *);
  double (*attribute)(void *);
} AttribStruct;


/* +++++++++++++++++++++++++++++++ */
/*      Constant Definition        */
/* +++++++++++++++++++++++++++++++ */


#define BOTTOM (-1)
#define MAXDIM 4
#define MAXTHREADS 512
#define NUMATTR 18
#define NUMDECISIONS 4
#define NUMBUCKETS 65536

#define PARENTS 0
#define ATTRIBUTES 1
#define PARENTS_AND_ATTRIBUTES 2

#define MXT_HISTO_SZ_LOG2 8
#define MXT_HISTO_SZ (1U << MXT_HISTO_SZ_LOG2)
#define MXT_HISTO_MASK (MXT_HISTO_SZ - 1)
#define NUM_DIGITS (sizeof(value) * CHAR_BIT + MXT_HISTO_SZ_LOG2 - 1) / MXT_HISTO_SZ_LOG2
#define NO_INLINE __attribute__((noinline))

// extern float 		g_max_greyval;		      	/* Maximum pixel intensity in data */
extern ulong 		g_max_levels;			/* Theoretical max level in data */
extern struct tms 	tstruct;			/* Timing structure */
extern clock_t 		start;		       		/* Clock structure */
extern AttribStruct 	AttribsArray[NUMATTR];		/* Attributes array */
extern DecisionStruct 	Decisions[NUMDECISIONS];	/* Pruning decision array */
// extern double           pixdim[3];
extern DataParams data_properties;

typedef enum {
  HORIZONTAL = 0,
  VERTICAL = 1,
  DEPTH = 2
} Direction;

int bit_scan_reverse(ulong val);
int bits_per_word_log2(void);
int bits_per_word(void);

typedef  void *(*f_new_aux_data)(AuxDataStore *);
typedef  void (*f_init_aux_data)(void *, bool , ulong , ulong,  ulong, value);
typedef  void (*f_delete_aux_data)(void *);
typedef  void (*f_add_to_aux_data)(void *, void *);
typedef  void (*f_merge_aux_data)(void *, void *);
typedef  void (*f_merge_to_aux_data)(AuxDataStore *, void **, void *, void *);
typedef  void (*f_clone_aux_data)(AuxDataStore *, void **, void *);
typedef  void (*f_create_mpi_aux_data)(void);

extern f_new_aux_data  new_aux_data;
extern f_init_aux_data init_aux_data;
extern f_delete_aux_data delete_aux_data;
extern f_add_to_aux_data  add_to_aux_data;  
extern f_merge_aux_data  merge_aux_data;
extern f_merge_to_aux_data merge_to_aux_data;
extern f_clone_aux_data clone_aux_data;
extern f_create_mpi_aux_data create_mpi_aux_data;


/*
extern  void *(*new_aux_data)(AuxDataStore *, void *);
extern  void (*init_aux_data)(void *,  bool , ulong , ulong,  ulong, value);
extern  void (*delete_aux_data)(void *);
extern  void (*add_to_aux_data)(void *, void *);
extern  void (*merge_aux_data)(void *, void *);
extern  void (*merge_to_aux_data)(AuxDataStore *, void **, void *, void *);
extern  void (*clone_aux_data)(AuxDataStore *, void **, void *);
extern  void (*create_mpi_aux_data)(void);

*/
/* INCLUDE */

//#include "cmdline.h"
#include "logc.h"
#include "mpihelper.h"
#include "checks.h"


#endif
