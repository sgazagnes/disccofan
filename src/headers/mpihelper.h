#ifndef MPIHELPER_H
#define MPIHELPER_H

extern MPI_Datatype mpi_bound_node_type;
extern MPI_Datatype mpi_borderindex_type;
extern MPI_Datatype mpi_attribute_type;
extern MPI_Datatype mpi_value_type;

int rank(void);
int np(void);
void init_mpi(void);
void finalize_mpi(void);
void create_mpi_value_type(void);
void create_mpi_boundnode_type(void);
void create_mpi_borderindex_type(void);
void create_mpi_area_type(void);
void create_mpi_extent_type(void);
void create_mpi_mean_type(void);
void create_mpi_wmean_type(void);
void create_mpi_inertia_type(void);
void create_mpi_inertiafull_type(void);
#endif
