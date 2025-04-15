#include "types.h"
#include "attributes.h"
#include "communication.h"


void send_boundary(Boundary *b, int dest, int merge_case) {

  int tag = 1; // MPI Tag is ignored 
  int err;

  if(b->size_curr >  INT_MAX)
    warn("MAXIMUM AMOUNT OF DATA REACHED FOR OPENMPI, MIGHT CAUSE ERRORS");

  err = MPI_Send(b->offset, 7, MPI_UINT64_T, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 300);
  err = MPI_Send(b->dims,   3, MPI_UINT64_T,  dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 301);
  err = MPI_Send(b->merge_idx, b->offset[6], MPI_UINT64_T, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 302);
  err = MPI_Send(b->array, b->size_curr, mpi_bound_node_type,  dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 303);
  if(merge_case){
    err = MPI_Send(b->attribute, b->size_curr, mpi_attribute_type,  dest, tag, MPI_COMM_WORLD);
    check_mpi_error(err, 304);
  }
  err = MPI_Send(b->border_par, b->size_curr,  mpi_borderindex_type, dest, tag, MPI_COMM_WORLD);
  check_mpi_error(err, 305);
} /* send_boundary */

Boundary *receive_boundary(int src, ulong size_item, int merge_case)
{
  /* Initialize MPI */
  int tag = 1, err; // MPI Tag is ignored
  MPI_Status status;
  int recv_size;

  Boundary *b = calloc(1, sizeof(Boundary));
  check_alloc(b, 300);

  err = MPI_Recv(b->offset, 7, MPI_UINT64_T, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 307);
  err = MPI_Recv(b->dims, 3, MPI_UINT64_T, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 308);

  b->merge_idx = malloc(b->offset[6] * sizeof(ulong));
  check_alloc(b->merge_idx, 302);

  err = MPI_Recv(b->merge_idx, b->offset[6], MPI_UINT64_T, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 309);

  MPI_Probe(src, tag, MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, mpi_bound_node_type, &recv_size);

  b->array = malloc(recv_size * sizeof(BoundaryNode));
  check_alloc(b->array, 303);
  if (merge_case)
  {
    b->size_attr = size_item;
    b->attribute = malloc(recv_size * size_item);
    check_alloc(b->attribute, 304);
  }
  b->border_par = malloc(recv_size * sizeof(BorderIndex));
  check_alloc(b->border_par, 305);
  if (merge_case != 1)
  {
    b->border_lr = malloc(recv_size * sizeof(BorderIndex));
    check_alloc(b->border_lr, 306);
  }
  /* Receive array, parents and idx of attributes */

  err = MPI_Recv(b->array, recv_size, mpi_bound_node_type, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 310);
  if (merge_case)
  {
    err = MPI_Recv(b->attribute, recv_size, mpi_attribute_type, src, tag, MPI_COMM_WORLD, &status);
    check_mpi_error(err, 311);
  }
  err = MPI_Recv(b->border_par, recv_size, mpi_borderindex_type, src, tag, MPI_COMM_WORLD, &status);
  check_mpi_error(err, 312);

  b->size_curr = b->size_init = b->size_alloc = (ulong)recv_size;

  /* Init border_pars and levelroots*/
  if (merge_case != 1)
  {
#pragma omp parallel for
    for (ulong i = 0; i < b->size_curr; i++)
    {
      b->border_lr[i] = (BorderIndex){.b = b, .i = BOTTOM};
      b->border_par[i].b = b;
    }
  }
  else
  {
#pragma omp parallel for
    for (ulong i = 0; i < b->size_curr; i++)
    {
      b->border_par[i].b = b;
    }
  }
  return b;
} /* receive_boundary */



void send_updated_boundary(Boundary *b, int dest, int merge_case) {
  //Send information through MPI
  int tag = 1; // MPI Tag is ignored 
  int err;

  if (b->size_curr >  INT_MAX)
    warn("MAXIMUM AMOUNT OF DATA REACHED FOR OPENMPI, MIGHT CAUSE ERRORS");
  
  if(merge_case == ATTRIBUTES){
    err = MPI_Send(b->attribute, b->size_curr, mpi_attribute_type, dest, tag, MPI_COMM_WORLD);
    check_mpi_error(err, 315);
    return;
  } else {
    err = MPI_Send(b->border_par, b->size_curr, mpi_borderindex_type, dest, tag, MPI_COMM_WORLD);
    check_mpi_error(err, 314);
    if(merge_case == PARENTS_AND_ATTRIBUTES){
      err = MPI_Send(b->attribute, b->size_curr, mpi_attribute_type,  dest, tag, MPI_COMM_WORLD);
      check_mpi_error(err, 315);
    }
    if (b->size_curr > b->size_init) {
      err = MPI_Send(b->array+b->size_init, b->size_curr - b->size_init, mpi_bound_node_type, dest, tag, MPI_COMM_WORLD);
      check_mpi_error(err, 316);
    }
  }

} /* send_updated_bobundary */

Boundary *receive_updated_boundary(Boundary *b, int src, int merge_case)
{
  /* Initialize MPI */
  int tag = 1; // MPI Tag is ignored
  MPI_Status status;
  int err;
  int recv_size;

  if (merge_case == 1)
  {
    MPI_Probe(src, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_INT64_T, &recv_size);
    err = MPI_Recv(b->attribute, recv_size, mpi_attribute_type, src, tag, MPI_COMM_WORLD, &status);
    check_mpi_error(err, 322);
  }
  else
  {
    MPI_Probe(src, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, mpi_borderindex_type, &recv_size);
    if ((ulong)recv_size > b->size_curr)
    {

      // Reallocation of b to receive the updated btree

      b->array = realloc(b->array, recv_size * sizeof(BoundaryNode));
      check_alloc(b->array, 308);
      if (merge_case == PARENTS_AND_ATTRIBUTES)
      {
        b->attribute = realloc(b->attribute, recv_size * b->size_attr);
        check_alloc(b->attribute, 309);
      }
      b->border_par = realloc(b->border_par, recv_size * sizeof(BorderIndex));
      check_alloc(b->border_par, 311);
      b->border_ori = realloc(b->border_ori, recv_size * sizeof(BorderIndex));
      check_alloc(b->border_ori, 310);
      b->border_lr = realloc(b->border_lr, recv_size * sizeof(BorderIndex));
      check_alloc(b->border_lr, 312);

      err = MPI_Recv(b->border_par, recv_size, mpi_borderindex_type, src, tag, MPI_COMM_WORLD, &status);
      check_mpi_error(err, 318);
      if (merge_case == PARENTS_AND_ATTRIBUTES)
      {
        err = MPI_Recv(b->attribute, recv_size, mpi_attribute_type, src, tag, MPI_COMM_WORLD, &status);
        check_mpi_error(err, 319);
      }
      err = MPI_Recv(b->array + b->size_curr, recv_size - b->size_curr, mpi_bound_node_type,
                     src, tag, MPI_COMM_WORLD, &status);
      check_mpi_error(err, 320);

      b->size_alloc = (ulong)recv_size;
    }
    else
    {
      err = MPI_Recv(b->border_par, recv_size, mpi_borderindex_type, src, tag, MPI_COMM_WORLD, &status);
      check_mpi_error(err, 318);
      if (merge_case == PARENTS_AND_ATTRIBUTES)
      {
        err = MPI_Recv(b->attribute, recv_size, mpi_attribute_type, src, tag, MPI_COMM_WORLD, &status);
        check_mpi_error(err, 319);
      }
    }

/* Assigning area and parents */
#pragma omp parallel for
    for (ulong i = 0; i < b->size_alloc; i++)
    {
      b->border_par[i].b = b;
      if (i >= b->size_curr)
        b->border_lr[i] = b->border_ori[i] = (BorderIndex){.b = b, .i = BOTTOM};
    }
    b->size_curr = b->size_alloc;
  }

  return b;
} /* receive_updated_boundary */
