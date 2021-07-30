#define _GNU_SOURCE
#include <stddef.h>
//#include <stdio.h>
#include <errno.h>
//#include <stdlib.h>
#include <string.h>

//#include "flood.h"
//include "constants.h"
#include "writefile.h"
#include "types.h"
#include "checks.h"
#include "logc.h"
#include "attributes.h"
#include "mpihelper.h"


void write_boundary_file_ascii(Boundary *b, const char *filename) {
  FILE *fp = fopen(filename, "w");
  // check_not_null(fp, errno);

  fprintf(fp, "array: %p\n", (void*) b->array);
  fprintf(fp, "size: %zu\n", b->size_curr);
  fprintf(fp, "offsets: %zu %zu %zu %zu %zu %zu %zu\n", b->offset[0], b->offset[1], b->offset[2], b->offset[3], b->offset[4], b->offset[5], b->offset[6]);
  fprintf(fp, "+ b-index   +   index   +   parent  +  bparent + gval  \n");
  ulong size = b->size_curr;
  for (ulong i = 0; i < size; i++) {
    fprintf(fp, "|%11ld", b->array[i].border_idx);
    fprintf(fp, "|%11ld", b->array[i].index);
    fprintf(fp, "|%11ld", b->border_par[i].i);
    fprintf(fp, "|%f", (float) b->array[i].gval);
    fprintf(fp, "|\n");
  }
  fprintf(fp, "\n");

  int err = fclose(fp);
  check_file_close(err, filename);
} /* write_boundary_file_ascii */

void write_boundary_file_txt(Boundary *b, const char *fname) {
  char *fname_attr;
  FILE *fp = fopen(fname, "wb");
  //check_not_null(fp, errno);
  asprintf(&fname_attr, "%s_attr", fname);
  FILE *fp_attr = fopen(fname_attr, "wb");
  check_not_null(fp_attr, errno);

  fprintf(fp, "%zu\n", b->size_curr);
  fprintf(fp, "%zu %zu %zu %zu %zu %zu %zu\n", b->offset[0], b->offset[1], b->offset[2], b->offset[3], b->offset[4], b->offset[5], b->offset[6]);
  ulong size = b->size_curr;
  for (ulong i = 0; i < size; i++) {
    fprintf(fp, "%ld ", b->array[i].index);
    fprintf(fp, "%ld ", b->array[i].border_idx);
    fprintf(fp, "%11ld ", b->border_par[i].i);
    fprintf(fp, "%f ", (float)b->array[i].gval);
    if(b->attribute_idx[i] == BOTTOM){
      fprintf(fp, "0 ");
    }else{
      void *b_attr = (char *) b->store->data + (b->attribute_idx[i] * b->store->size_item);
      write_aux_file_binary(fp_attr, b_attr);
      fprintf(fp, "1 ");
    }
  }
  fprintf(fp, "\n");

  int err = fclose(fp);
  check_file_close(err, fname);
  err = fclose(fp_attr);
  check_file_close(err, fname_attr);
  free(fname_attr);
} /* write_boundary_file_binary */

void write_tree_file_acsii(Node *m, ulong size, const char *filename) {
  FILE *fp = fopen(filename, "w");
  //  check_not_null(fp, errno);

  fprintf(fp, "size: %zu\n", size);
  fprintf(fp, "+ b-index   +   index   +   parent   + gval +\n");

  for (ulong i = 0; i < size; i++) {
    //   fprintf(fp, "|%11ld", m[i].border_idx);
    fprintf(fp, "|%11ld", i);
    fprintf(fp, "|%11ld", m->parent[i]);
    fprintf(fp, "|%f", (float) m->gval[i]);
    fprintf(fp, "|\n");
  }
  fprintf(fp, "\n");

  int err = fclose(fp);
  check_file_close(err, filename);
} /* write_boundary_file_ascii */

void write_tree_file_txt(Arguments args, Node *tree, ulong size, ulong *dims_ti, ulong *dims_i, int bitpix) {
  char *fname_file;
  char *fname_attr;
  char *image_name;
  image_name = strrchr( args.inprefix_arg,'/');
  if(image_name == NULL) image_name = args.inprefix_arg;
  else image_name++;
  if(np() == 1) asprintf(&fname_file, "saves/Component_tree_%s.txt", image_name);
  else asprintf(&fname_file, "saves/Component_tree_%s_part_%d.txt", image_name, rank());
  FILE *fp = fopen(fname_file, "wb");
  // check_not_null(fp, errno);

  if(np() == 1) asprintf(&fname_attr, "saves/Attributes_%s.txt", image_name);
  else asprintf(&fname_attr, "saves/Attributes_%s_part_%d.txt", image_name, rank());
  FILE *fp_attr = fopen(fname_attr, "wb");
  // check_not_null(fp_attr, errno);

 
  fprintf(fp, "Image prefix: %s \n", args.inprefix_arg);
  fprintf(fp, "Image type: %s \n", args.intype_arg);
  fprintf(fp, "Image dimensions %lu x %lu x %lu \n", dims_i[0], dims_i[1], dims_i[2]);
  fprintf(fp, "Tile dimensions %lu x %lu x %lu \n", dims_ti[0], dims_ti[1], dims_ti[2]);
  fprintf(fp, "Tree size: %lu\n", size);
  fprintf(fp, "Bit-per-pixel: %d\n", bitpix);
  fprintf(fp, "Attribute choice: %d\n", args.attribute_arg);
  fprintf(fp, "Number of processes used: %d \n", np());
  fprintf(fp, "Number of threads used: %d \n", args.threads_arg);
  fprintf(fp, "Tree: %s\n", args.tree_arg);
  fprintf(fp, "Operation: %s\n", args.morphology_arg);
  fprintf(fp, "Pruning choice: %d\n", args.decision_arg);
  fprintf(fp, "Filter: %s\n", args.filter_arg);
  fprintf(fp, "Flooding function: %d\n", args.flood_arg);
  fprintf(fp, "Connectivity: %d\n", args.connectivity_arg);
  for (ulong i = 0; i < size; i++) {
    fprintf(fp, "%ld ", i);
    fprintf(fp, "%ld ", tree->parent[i]);
    if(bitpix>0)  fprintf(fp, "%u ", tree->gval[i]);
    else fprintf(fp, "%f ", (float) tree->gval[i]);
    if(is_levelroot(tree,i) && tree->attribute[i]){
      write_aux_file_binary(fp_attr, tree->attribute[i]);
      fprintf(fp, "1 \n");
    }else
      fprintf(fp, "0 \n");
  }
  fprintf(fp, "\n");

  int err = fclose(fp_attr);
  check_file_close(err, fname_attr);
  err = fclose(fp);
  check_file_close(err, fname_file);
} /* write_maxtree_file_binary */




void write_tree_binary(Arguments args, Node *tree, ulong size) {
  char *fname_file;
  char *fname_attr;
  char *image_name;
  image_name = strrchr( args.inprefix_arg,'/');
  if(image_name == NULL) image_name = args.inprefix_arg;
  else image_name++;
  if(np() == 1) asprintf(&fname_file, "Component_tree_%s", image_name);
  else asprintf(&fname_file, "Component_tree_%s_part_%d", image_name, rank());
  FILE *fp = fopen(fname_file, "wb");
  //  check_not_null(fp, errno);

  if(np() == 1) asprintf(&fname_attr, "Attributes_%s", image_name);
  else asprintf(&fname_file, "Attributes_%s_part_%d", image_name, rank());
  FILE *fp_attr = fopen(fname_attr, "wb");
  // check_not_null(fp_attr, errno);

  fwrite(tree, sizeof(Node), size, fp);

  int err = fclose(fp_attr);
  check_file_close(err, fname_attr);
  err = fclose(fp);
  check_file_close(err, fname_file);
} /* write_maxtree_file_binary */


void write_params_txt(Arguments args,  ulong *dims_i, int bitpix) {
  char *fname_file;
  char *image_name;
  image_name = strrchr( args.inprefix_arg,'/');
  if(image_name == NULL) image_name = args.inprefix_arg;
  else image_name++;
  asprintf(&fname_file, "Parameters_%s.txt", image_name);
  FILE *fp = fopen(fname_file, "wb");
  // check_not_null(fp, errno);

 
  fprintf(fp, "Image prefix: %s \n", args.inprefix_arg);
  fprintf(fp, "Image type: %s \n", args.intype_arg);
  fprintf(fp, "Image dimensions %lu x %lu x %lu \n", dims_i[0], dims_i[1], dims_i[2]);
  //fprintf(fp, "Tile dimensions %lu x %lu x %lu \n", dims_ti[0], dims_ti[1], dims_ti[2]);
  //fprintf(fp, "Tree size: %lu\n", size);
  fprintf(fp, "Bit-per-pixel: %d\n", bitpix);
  fprintf(fp, "Attribute choice: %d\n", args.attribute_arg);
  fprintf(fp, "Number of processes used: %d \n", np());
  fprintf(fp, "Number of threads used: %d \n",  args.threads_arg);
  fprintf(fp, "Tree: %s\n", args.tree_arg);
  fprintf(fp, "Operation: %s\n", args.morphology_arg);
  fprintf(fp, "Pruning choice: %d\n", args.decision_arg);
  fprintf(fp, "Filter: %s\n", args.filter_arg);
  fprintf(fp, "Flooding function: %d\n", args.flood_arg);
  fprintf(fp, "Connectivity: %d\n", args.connectivity_arg);
  
  int err = fclose(fp);
  check_file_close(err, fname_file);
} /* write_params_file_txt */

Node *read_tree_file_txt(Arguments *args, ulong* size_tree, ulong *size_tile, ulong *dims_ti, ulong *dims_i, int *bitpix) {
  char *fname_file;
  
  if(np() == 1) asprintf(&fname_file, "saves/Component_tree_%s.txt",args->inprefix_arg);
  else asprintf(&fname_file, "saves/Component_tree_%s_part_%d.txt", args->inprefix_arg, rank());
  FILE *fp = fopen(fname_file, "rb");
  // check_not_null(fp, errno);
  int numproc;
  
  fscanf(fp, "Image prefix: %s \n", args->inprefix_arg);
  fscanf(fp, "Image type: %s \n", args->intype_arg);
  fscanf(fp, "Image dimensions %lu x %lu x %lu \n", dims_i, dims_i + 1, dims_i + 2);
  fscanf(fp, "Tile dimensions %lu x %lu x %lu \n", dims_ti, dims_ti + 1, dims_ti + 2);
  fscanf(fp, "Tree size: %lu\n", size_tree);
  fscanf(fp, "Bit-per-pixel: %d\n", &args->bitsperpixel_arg);
  fscanf(fp, "Attribute choice: %d\n", &args->attribute_arg);
  fscanf(fp, "Number of processes used: %d \n", &numproc);
  if(numproc != np())     MPI_Abort(MPI_COMM_WORLD, 001);
  fscanf(fp, "Number of threads used: %d \n", &args->threads_arg);
  fscanf(fp, "Tree: %s\n", args->tree_arg);
  fscanf(fp, "Operation: %s\n", args->morphology_arg);
  fscanf(fp, "Pruning choice: %d\n", &args->decision_arg);
  fscanf(fp, "Filter: %s\n", args->filter_arg);
  fscanf(fp, "Flooding function: %d\n", &args->flood_arg);
  fscanf(fp, "Connectivity: %d\n", &args->connectivity_arg);
  *bitpix = args->bitsperpixel_arg;
  *size_tile = dims_ti[0]*dims_ti[1]*dims_ti[2];
  
  Node *tree = malloc(*size_tree * sizeof(Node));

  for (ulong i = 0; i < *size_tree; i++) {
    fscanf(fp, "%ld ", &tree->parent[i]);
    if(!FLOAT_TYPE) fscanf(fp, "%u ", (uint *) &tree->gval[i]);
    else fscanf(fp, "%f ", (float *) &tree->gval[i]);
    //  fscanf(fp, "%d \n", &tree[i].attribute);
  }
 
  int err = fclose(fp);
  check_file_close(err, fname_file);

  return tree;
} /* read_tree_file_binary */

void write_area_file_binary(FILE *fp,  void *areaattr) {

  AreaData *areadata = areaattr;
  fprintf(fp, "%lu ", areadata->area);
} /* write_area_file_binary */

void write_encl_rect_file_binary(FILE *fp, void *rectattr) {

  EnclRectData *rectdata = rectattr;
  fprintf(fp, "%lu ", rectdata->minX);
  fprintf(fp, "%lu ", rectdata->maxX);
  fprintf(fp, "%lu ", rectdata->minY);
  fprintf(fp, "%lu ", rectdata->maxY);
  fprintf(fp, "%lu ", rectdata->minZ);
  fprintf(fp, "%lu ", rectdata->maxZ);
} /* write_encl_rect_file_binary */

void write_inertia_file_binary(FILE *fp,  void *inertiaattr) {

  InertiaData *inertiadata = inertiaattr;
  fprintf(fp, "%lu ",  inertiadata->area);
  fprintf(fp, "%f ",  inertiadata->sumX);
  fprintf(fp, "%f ", inertiadata->sumY);
  fprintf(fp, "%f ", inertiadata->sumZ);
  fprintf(fp, "%f ", inertiadata->sumR2);
} /* write_inertia_file_binary */

void *read_area_file_binary(AuxDataStore *store, FILE *fp_attr){
  ulong area;
  void *attr = NULL;	
  fscanf(fp_attr, "%lu ", &area);
  attr = load_area_data(store, &area);
  return attr;
}

void *read_encl_rect_file_binary(AuxDataStore *store, FILE *fp_attr){
  ulong *initval = malloc(6 * sizeof(ulong));
  void *attr = NULL;	
  fscanf(fp_attr, "%lu ", initval);
  fscanf(fp_attr, "%lu ", initval + 1);
  fscanf(fp_attr, "%lu ", initval + 2);
  fscanf(fp_attr, "%lu ", initval + 3);
  fscanf(fp_attr, "%lu ", initval + 4);
  fscanf(fp_attr, "%lu ", initval + 5);

  attr = load_encl_rect_data(store, initval);
  free(initval);
  return attr;
}


void *read_inertia_file_binary(AuxDataStore *store, FILE *fp_attr){
  ulong *initval = malloc(5 * sizeof(ulong));
  void *attr = NULL;	
  fscanf(fp_attr, "%lu ", initval);
  fscanf(fp_attr, "%lu ", initval + 1);
  fscanf(fp_attr, "%lu ", initval + 2);
  fscanf(fp_attr, "%lu ", initval + 3);
  fscanf(fp_attr, "%lu ", initval + 4);

  attr = load_encl_rect_data(store, initval);
  free(initval);
  return attr;
}



/*Boundary *read_boundary_file_binary(Boundary *b, const char *fname) {
  char *fname_attr;
  FILE *fp = fopen(fname, "wb");
  check_not_null(fp, errno);
  asprintf(&fname_attr, "%s_attr", fname);
  FILE *fp_attr = fopen(fname_attr, "wb");
  check_not_null(fp_attr, errno);

  fprintf(fp, "%zu\n", b->size_curr);
  fprintf(fp, "%zu %zu %zu %zu %zu %zu %zu\n", b->offset[0], b->offset[1], b->offset[2], b->offset[3], b->offset[4], b->offset[5], b->offset[6]);
  size_t size = b->size_curr;
  for (size_t i = 0; i < size; i++) {
    fprintf(fp, "%ld ", b->array[i].index);
    fprintf(fp, "%ld ", b->array[i].border_idx);
    fprintf(fp, "%11ld ", b->border_par[i].i);
    fprintf(fp, "%d ", b->array[i].gval);
    if(b->array[i].attribute_idx == BOTTOM){
      fprintf(fp, "0 ");
    }else{
      void *b_attr = b->store->data + b->array[i].attribute_idx * b->store->size_item;
      write_aux_file_binary(fp_attr, b_attr);
      fprintf(fp, "1 ");
    }
  }
  fprintf(fp, "\n");

  int err = fclose(fp);
  check_file_close(err, filename);
} */

/*void queues_to_file(Queue *hq, value maxvalue, const char *filename) {
  FILE *fp = fopen(filename, "w");
  check_not_null(fp, errno);

  for (value i = 0; i < maxvalue; i++) {
    fprintf(fp, "hq[%d]: %lu pixels head: %lu tail: %lu \n", i, hq[i].pixels[hq[i].head], hq[i].head, hq[i].tail);
  }

  int err = fclose(fp);
  check_file_close(err, filename);
  }*/
