#include "types.h"

MPI_Datatype mpi_bound_node_type;
MPI_Datatype mpi_borderindex_type;
MPI_Datatype mpi_attribute_type;
MPI_Datatype mpi_value_type;

int rank(void) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
} /* rank */

int np(void) {
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  return numprocs;
} /* np */


void init_mpi(void) {
  int initialized, provided;

  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
  }
  if(MPI_THREAD_FUNNELED != provided)
    warn("MPI asked %d, provided: %d", MPI_THREAD_FUNNELED, provided);
} /* init_mpi */


void finalize_mpi(void) {
  int finalized;
  MPI_Finalized(&finalized);
  if (!finalized) {
    MPI_Finalize();
  }
} /* finalize_mpi */

void create_mpi_value_type(void) {
  /* Create MPI_Type for MaxTree Nodes */
  
  if(FLOAT_TYPE)
    mpi_value_type = MPI_FLOAT;
  else if(sizeof(value) == 8)
    mpi_value_type = MPI_UINT64_T;
  else if(sizeof(value) == 4)
    mpi_value_type = MPI_UINT32_T;
  else if(sizeof(value) == 2)
    mpi_value_type = MPI_UINT16_T;
  else if(sizeof(value) == 1)
    mpi_value_type = MPI_UINT8_T;

  MPI_Type_commit(&mpi_value_type);
} /* create_mpi_boundnode_type */

void create_mpi_boundnode_type(void) {
  /* Create MPI_Type for MaxTree Nodes */
  const int nitems = 3;
  int blocklengths[3] = {1, 1, 1};
  MPI_Aint offsets[3];
  offsets[0] = offsetof(BoundaryNode, index);
  offsets[1] = offsetof(BoundaryNode, gval);
  offsets[2] = offsetof(BoundaryNode, border_idx);
  /* skipped global index here */

  MPI_Datatype types[3] = {
    MPI_INT64_T,
    mpi_value_type,
    MPI_INT64_T,
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_bound_node_type);
  MPI_Type_commit(&mpi_bound_node_type);
} /* create_mpi_boundnode_type */

void create_mpi_borderindex_type(void) {
  

  // MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_borderindex_type);
  MPI_Type_vector(1, 2, 2, MPI_INT64_T, &mpi_borderindex_type);
  MPI_Type_commit(&mpi_borderindex_type);
  //
} /* create_mpi_borderindex_type */


void create_mpi_area_type(void) {
  /* Create MPI_Type for MaxTree Nodes */

  const int nitems = 1;
  int blocklengths[1] = {1};
  MPI_Aint offsets[1];
  offsets[0] = offsetof(AreaData, area);

  /* skipped global index here */

  MPI_Datatype types[1] = {
    MPI_UNSIGNED_LONG
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_attribute_type);
  MPI_Type_commit(&mpi_attribute_type);
} /* create_mpi_area_type */


void create_mpi_rect_type(void) {
  /* Create MPI_Type for MaxTree Nodes */

  const int nitems = 6;
  int blocklengths[6] = {1,1,1,1,1,1};
  MPI_Aint offsets[6];
  offsets[0] = offsetof(EnclRectData, minX);
  offsets[1] = offsetof(EnclRectData, minY);
  offsets[2] = offsetof(EnclRectData, minZ);
  offsets[3] = offsetof(EnclRectData, maxX);
  offsets[4] = offsetof(EnclRectData, maxY);
  offsets[5] = offsetof(EnclRectData, maxZ);

  /* skipped global index here */

  MPI_Datatype types[6] = {
    MPI_UINT64_T,
    MPI_UINT64_T,
    MPI_UINT64_T,
    MPI_UINT64_T,
    MPI_UINT64_T,
    MPI_UINT64_T
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_attribute_type);
  MPI_Type_commit(&mpi_attribute_type);
} /* create_mpi_rect_type */


void create_mpi_inertia_type(void) {
  /* Create MPI_Type for MaxTree Nodes */

  const int nitems = 5;
  int blocklengths[5] = {1,1,1,1,1};
  MPI_Aint offsets[5];
  offsets[0] = offsetof(InertiaData, area);
  offsets[1] = offsetof(InertiaData, sumX);
  offsets[2] = offsetof(InertiaData, sumY);
  offsets[3] = offsetof(InertiaData, sumZ);
  offsets[4] = offsetof(InertiaData, sumR2);


  /* skipped global index here */

  MPI_Datatype types[5] = {
    MPI_UINT64_T,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_attribute_type);
  MPI_Type_commit(&mpi_attribute_type);
} /* create_mpi_inertia_type */


