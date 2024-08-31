#ifndef WRITEFILE_H
#define WRITEFILE_H

#include <stdlib.h>
#include "types.h"
#include "cmdline.h"
void write_boundary_file_ascii(Boundary *b, const char *filename);
void write_boundary_file_txt(Boundary *b, const char *filename);
void write_tree_file_ascii(Node *b, ulong size, const char *filename);
void write_tree_file_txt(Arguments args, Node *tree, ulong size, ulong *dims_ti, ulong *dims_i, int bitpix);
void write_tree_binary(Arguments args, Node *tree, ulong size);
void write_params_txt(Arguments args,  ulong *dims_i, int bitpix);
void write_area_file_binary(FILE *fp_attr, void *areaattr);
void write_encl_rect_file_binary(FILE *fp,  void *rectattr);
void write_inertia_file_binary(FILE *fp,  void *inertiaattr);
void *read_area_file_binary(AuxDataStore *store, FILE *fp_attr);
void *read_encl_rect_file_binary(AuxDataStore *store, FILE *fp_attr);
void *read_inertia_file_binary(AuxDataStore *store, FILE *fp_attr);
Node *read_tree_file_txt(Arguments *args,  ulong* size_tree, ulong *size_tile, ulong *dims_ti, ulong *dims_i, int *bitpix);
void queues_to_file(Queue *hq, value maxvalue, const char *filename);

#endif
