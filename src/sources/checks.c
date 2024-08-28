#include "types.h"


void check_alloc(void *array, int code) {
  if (NULL == array) {
    if (code < 100)
      error("Allocation failed, out of memory? code: %d (seq_main.c || step_main.c)", code);	
    else if (code < 200)
      error("Allocation failed, out of memory? code: %d (attributes.c)", code);
    else if (code < 300)
      error("Allocation failed, out of memory? code: %d (boundary.c || boundary_step.c)", code);
    else if (code < 400)
      error("Allocation failed, out of memory? code: %d (communication.c)", code);
    else if (code < 500)
      error("Allocation failed, out of memory? code: %d (queue.c)", code);
    else if (code < 600)
      error("Allocation failed, out of memory? code: %d (tree_flood.c || tree_step.c)", code);
    else if (code < 700)
      error("Allocation failed, out of memory? code: %d (tree_filt.c)", code);
    else if (code < 800)
      error("Allocation failed, out of memory? code: %d (image.c)", code);
    MPI_Abort(MPI_COMM_WORLD, code);
  }
} /* check_alloc */


void check_mpi_error(int errorval, int code) {
  if (MPI_SUCCESS != errorval) {
    error("MPI communication error %d! %d", errorval, code);
    MPI_Abort(MPI_COMM_WORLD, code);
  }
} /* check_mpi_error */

void check_file_close(int errorval, const char *filename) {
  if (0 != errorval) {
    error("Error closing file %s : errorval", filename);
    MPI_Abort(MPI_COMM_WORLD, errorval);
  }
} /* check_file_close */



void check_boundary(Boundary *b, int attrib_choice) {
  info("Proc %d, Boundary exists of array %p, dims %zu %zu %zu, offset %zu %zu %zu %zu %zu %zu %zu, size %zu \n", rank(), (void *)  b->array, b->dims[0], b->dims[1], b->dims[2],
	   b->offset[0], b->offset[1], b->offset[2], b->offset[3], b->offset[4],
	   b->offset[5], b->offset[6], b->size_curr);
  
   for (ulong a = 0; a < b->size_curr; a++) {
     double attrib = (*AttribsArray[attrib_choice].attribute)((char *) b->attribute + a * b->size_attr);
     printf("%zu: index %ld gval %f border_idx %ld border_par %ld attr %lf \n",
	    a, (b->array[a]).index, (float) (b->array[a]).gval, (b->array[a]).border_idx, b->border_par[a].i, attrib);
   }
} /* check_boundary */

void check_maxtree(Node *tree, double (*attribute)(void *), ulong size) {
  // ulong attribb[6]={0};
  for (ulong a = 0; a < size; a++) {
    double attrib = (*attribute)(tree->attribute + a * tree->size_attr);
    printf("index %ld,  parent %ld, gval %f, attribute %lf\n",
	   a, tree->parent[a],  (float) tree->gval[a],  attrib);
  }
} /* check_maxtree */



