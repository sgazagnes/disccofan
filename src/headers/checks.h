#ifndef CHECKS_H
#define CHECKS_H

void check_alloc(void *array, int code);
void check_mpi_error(int errorval, int code);
void check_file_close(int errorval, const char *filename);
//void check_file_arg(Arguments *args);
//void check_operation(Arguments *args, Node *tree, ulong offset);
void check_boundary(Boundary *b, int attrib_choice);
void check_maxtree(Node *maxtree, double (*attribute)(void *), ulong size);

#endif
