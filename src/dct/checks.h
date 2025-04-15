#ifndef CHECKS_H
#define CHECKS_H

void set_border(Arguments *args, Node *tree);
ulong *attribute_offsets(Arguments *args, ulong dims[3]);
void set_flooding(Arguments *args, int bit_depth);
void set_connectivity(Arguments *args, ulong dims_z);
void check_alloc(void *array, int code);
void check_not_null(void *ptr, int code);
void check_mpi_error(int errorval, int code);
void check_file_close(int errorval, const char *filename);
void check_area_size(ulong size, ulong area);
void check_decision(Arguments *args);
void check_file_arg(Arguments *args);
void check_interactive(Arguments *args);
void check_threads(Arguments *args);
void check_operation(Arguments *args, Node *tree, value offset);
void check_boundary(Boundary *b);
void check_maxtree(Node *maxtree, double (*attribute)(void *), ulong size);
void check_bytes(void);
void check_leros(Node *tree);

#endif
