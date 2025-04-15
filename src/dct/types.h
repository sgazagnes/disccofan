#ifndef TYPES_H
#define TYPES_H

#define _GNU_SOURCE

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

/* +++++++++++++++++++++++++++++++ */
/*  	 Structure definition      */
/* +++++++++++++++++++++++++++++++ */

/* Parsing structure */
typedef struct gengetopt_args_info Arguments;

/* Attribute storage Structure */
typedef struct _AuxDataStore{
  size_t size_item;           		/* Size of attribute used */
  ulong size_alloc;			/* Allocation size */	
  ulong item_curr;			/* Current number of attributes */
  void  *data;		       		/* Attributes */ 
} AuxDataStore;

/* Tree node Structure */
typedef struct _Node {
  bool 		dens;
  ulong		size;
  bool 		border[6];
  value		*gval;			/* Pixel intensity */
  value		*gval_dens;
  idx 		*parent;
  idx 	 	*border_idx; 		/* index in boundary tree */
  void 		**attribute;		/* Pointer to attribute */
  AuxDataStore  **store;
} Node;

typedef struct mt_object_data
{  
  Node* mt;
  float bg_variance;
  float gain;
  float move_factor;
  ubyte* flags;
  idx* relevant_indices;
  idx relevant_indices_len;
  idx* closest_significant_ancestors;
  idx* main_branches;
  idx* main_power_branches;
  idx* object_ids;
  idx num_significant_nodes;
  idx num_objects;
  int (*node_significance_test)(struct mt_object_data *, idx);
  // Pointer to function that takes an mt_o and returns nothing
  void (*significant_nodes)(struct mt_object_data *);
  void *node_significance_test_data;
  void (*node_significance_test_data_free)(struct mt_object_data *);

} mt_object_data;


/* Boundary tree node Structure */
typedef struct _BoundaryNode {
  value  	gval;			/* Pixel intensity */
  idx 		index;			/* Pixel index from local tree */
  idx 	 	border_idx; 		/* Index in boundary tree */
} BoundaryNode;

typedef struct _BorderIndex BorderIndex;

/* Boundary tree Structure */
typedef struct _Boundary {
  AuxDataStore  *store;			/* Attribute storage structure */
  BoundaryNode 	*array;			/* Node array */
  BorderIndex 	*border_par;   		/* Parent in the boundary tree */
  BorderIndex 	*border_ori;   		/* Position in previous boundary tree */
  BorderIndex 	*border_lr;    		/* Levelroot at same intensity in other merged tree */
  ulong		*merge_idx;		/* Index of node to be merged */
  idx 		*attribute_idx;		/* Index of attribute in storage */
  bool	 	*reached;		/* Boolean array */
  ulong 	offset[7]; 		/* Number of nodes in different sides  */
  ulong		dims[3];		/* Dimensions of boundary tree */
  ulong 	size_curr;     		/* Current size of tree */
  ulong 	size_init;		/* Initial size of tree */
  ulong 	size_alloc;		/* Allocated size of tree */
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
typedef struct _SortItem {
  value_t val;				/* Unsigned value */
  ulong rank;				/* Rank of value */
} SortItem;

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

/* Attribute Area Structure */
typedef struct _AreaData{
   ulong area;
} AreaData;

/* Attribute Inertia Structure */
typedef struct _InertiaData {
   ulong area;
   double sumX, sumY, sumZ, sumR2;
} InertiaData;


typedef struct _InertiaDataEoR {
  double area, sumval, sumval2;//, sumvald, sumval2d;
  double sumX, sumY, sumZ, sumX2, sumY2, sumZ2, sumXY, sumYZ, sumXZ;
  //double sumXd, sumYd, sumZd, sumX2d, sumY2d, sumZ2d, sumXYd, sumYZd, sumXZd;
} InertiaDataEoR;
/*
typedef struct _DataEoR {
  double area, sumval, oriX, oriY, oriZ;
  double sumX, sumY, sumZ, sumX2, sumY2, sumZ2, sumXY, sumYZ, sumXZ;
  } DataEoR;*/
typedef struct _MtoData{
   ulong area;
  value gval;
  double power;
  double volume;
} MtoData;
/* Attribute enclosing rectangle Structure */
typedef struct _EnclRectData {
   ulong minX;
   ulong minY;
   ulong minZ;
   ulong maxX;
   ulong maxY;
   ulong maxZ;
} EnclRectData;

/* Lambda vector Structure */
typedef struct _LambdaVec{
  uint num_lambdas;
  float *lambdas;    
} LambdaVec;

/* Pruning decision Structure */
typedef struct DecisionStruct {
  const char *name;
  void (*filter)(Node *, value *, bool *, ulong, ulong, double (*attribute)(void *), double);
} DecisionStruct;

/* Attribute Structure */
typedef struct _AttribStruct{
  const char 	 *name;
  ulong  size;
  void 	 *(*new_data)( AuxDataStore *, double *);
  void   (*delete_data)(void *);
  void   (*add_to_data)(void *, double *);
  void   (*merge_data)(void *, void *);
  void   (*merge_to_data)( AuxDataStore *, void **, void *, void *);
  void   (*clone_data)( AuxDataStore *, void **, void *);
  void   (*create_mpi_data)(void);
  double (*attribute)(void *);
  double *(*attribute_arr)(void *);
  // void 	 *(*read_data_file)( AuxDataStore *, FILE *);
  // void 	 (*write_data_file)( FILE *, void *);
} AttribStruct;


/* +++++++++++++++++++++++++++++++ */
/*      Constant Definition        */
/* +++++++++++++++++++++++++++++++ */


#define BOTTOM (-1)
#define MAXDIM 4
#define MAXTHREADS 512
#define NUMATTR 10
#define NUMDECISIONS 4
#define NUMBUCKETS 65536

#define MXT_HISTO_SZ_LOG2 8
#define MXT_HISTO_SZ (1U << MXT_HISTO_SZ_LOG2)
#define MXT_HISTO_MASK (MXT_HISTO_SZ - 1)
#define NUM_DIGITS (sizeof(value) * CHAR_BIT + MXT_HISTO_SZ_LOG2 - 1) / MXT_HISTO_SZ_LOG2

extern float 		g_max_greyval;		      	/* Maximum pixel intensity in data */
extern ulong 		g_max_levels;			/* Theoretical max level in data */
extern ulong 		DIMS[3];			/* Theoretical max level in data */
extern struct tms 	tstruct;			/* Timing structure */
extern clock_t 		start;		       		/* Clock structure */
extern AttribStruct 	AttribsArray[NUMATTR];		/* Attributes array */
extern DecisionStruct 	Decisions[NUMDECISIONS];	/* Pruning decision array */

typedef enum {
  HORIZONTAL = 0,
  VERTICAL = 1,
  DEPTH = 2
} Direction;

/* INCLUDE */

#include "cmdline.h"
#include "logc.h"
#include "mpihelper.h"
#include "checks.h"

#endif
