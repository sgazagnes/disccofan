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
  /* Create MPI_Type for area attributes */

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


void create_mpi_extent_type(void) {
  /* Create MPI_Type for extent attributes */

  const int nitems = 6;
  int blocklengths[6] = {1,1,1,1,1,1};
  MPI_Aint offsets[6];
  offsets[0] = offsetof(ExtentData, minX);
  offsets[1] = offsetof(ExtentData, minY);
  offsets[2] = offsetof(ExtentData, minZ);
  offsets[3] = offsetof(ExtentData, maxX);
  offsets[4] = offsetof(ExtentData, maxY);
  offsets[5] = offsetof(ExtentData, maxZ);

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
} /* create_mpi_extent_type */


void create_mpi_mean_type(void) {
  /* Create MPI_Type for mean attributes */

  const int nitems = 4;
  int blocklengths[4] = {1,1,1,1};
  MPI_Aint offsets[4];
  offsets[0] = offsetof(MeanData, area);
  offsets[1] = offsetof(MeanData, sumX);
  offsets[2] = offsetof(MeanData, sumY);
  offsets[3] = offsetof(MeanData, sumZ);


  /* skipped global index here */

  MPI_Datatype types[4] = {
    MPI_UINT64_T,
    MPI_UINT64_T,
    MPI_UINT64_T,
    MPI_UINT64_T
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_attribute_type);
  MPI_Type_commit(&mpi_attribute_type);
} /* create_mpi_mean_type */


void create_mpi_wmean_type(void) {
  /* Create MPI_Type for wmean attributes */

  const int nitems = 4;
  int blocklengths[4] = {1,1,1,1};
  MPI_Aint offsets[4];
  offsets[0] = offsetof(WMeanData, sumGval);
  offsets[1] = offsetof(WMeanData, sumXd);
  offsets[2] = offsetof(WMeanData, sumYd);
  offsets[3] = offsetof(WMeanData, sumZd);


  /* skipped global index here */

  MPI_Datatype types[4] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_attribute_type);
  MPI_Type_commit(&mpi_attribute_type);
} /* create_mpi_wmean_type */

void create_mpi_inertia_type(void) {
  /* Create MPI_Type for inertia attributes */

  const int nitems = 15;
  int blocklengths[15] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  MPI_Aint offsets[15];
  offsets[0] = offsetof(InertiaData, area);
  offsets[1] = offsetof(InertiaData, sumval);
  offsets[2] = offsetof(InertiaData, sumval2);
  offsets[3] = offsetof(InertiaData, sumX);
  offsets[4] = offsetof(InertiaData, sumY);
  offsets[5] = offsetof(InertiaData, sumZ);
  offsets[6] = offsetof(InertiaData, sumX2);
  offsets[7] = offsetof(InertiaData, sumY2);
  offsets[8] = offsetof(InertiaData, sumZ2);
  offsets[9] = offsetof(InertiaData, sumXY);
  offsets[10] = offsetof(InertiaData, sumYZ);
  offsets[11] = offsetof(InertiaData, sumXZ);
  offsets[12] = offsetof(InertiaData, sumXd);
  offsets[13] = offsetof(InertiaData, sumYd);
  offsets[14] = offsetof(InertiaData, sumZd);

  /* skipped global index here */

  MPI_Datatype types[15] = {
    MPI_UINT64_T,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_attribute_type);
  MPI_Type_commit(&mpi_attribute_type);
} /* create_mpi_inertia_type */


void create_mpi_morphology_type(void) {
  /* Create MPI_Type for area attributes */

  const int nitems = 3;
  int blocklengths[3] = {1,1,1};
  MPI_Aint offsets[3];
  offsets[0] = offsetof(MorphologyData, area);
  offsets[1] = offsetof(MorphologyData, sumval);
  offsets[2] = offsetof(MorphologyData, maxval);

  /* skipped global index here */

  MPI_Datatype types[3] = {
    MPI_UNSIGNED_LONG,
    MPI_DOUBLE,
    MPI_DOUBLE
  };

  MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_attribute_type);
  MPI_Type_commit(&mpi_attribute_type);
} /* create_mpi_morphology_type */
