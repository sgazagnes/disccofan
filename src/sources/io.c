#define _GNU_SOURCE 1

#include <hdf5.h>
#include <fitsio.h>
#include <FreeImage.h>
#include "types.h"
#include "nifti1.h"
#include "attributes.h"
#include "io.h"
#include "flood.h"

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

void read_grayscale_image(Arguments* args, value** img, 
                          uint grid[3], int* bitpix, uint myrank, uint myrank_2D, uint tile,
                          uint myrank_arr[3], uint file_style);
void read_rgb_channels(Arguments* args, value** img, 
                       uint grid[3], int* bitpix, uint myrank, uint myrank_2D, uint tile,
                       uint myrank_arr[3], uint file_style);
void read_lofar_data(Arguments* args, value** img, 
                     uint grid[3], int* bitpix, uint myrank, uint myrank_2D,
                     uint myrank_arr[3], uint file_style);

void read_fits(Arguments *args, value **img,  const char *fname,  uint grid[3], int *bitpix, uint myrank, uint myrank2D, uint myrank_arr[3], uint file_style);
void read_hdf5(Arguments *args, value **img, const char *fname,  uint grid[3], int *bitpix, uint myrank, uint myrank2D, uint myrank_arr[3], uint file_style);
void  read_basic(Arguments *args, value **img, const char *fname,   uint grid[3], int *bitpix, uint myrank, uint myrank2D, uint myrank_arr[3], uint file_style);
void read_nifti_file(Arguments *args,value **img,  const char *fname,  int *bitpix);
void write_hdf5(Arguments *args, const char* fname, const char *dataset_out, value *out,   ulong dims_T[3], ulong dims[3]);
void write_fits(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3]);
void write_basic(Arguments *args, const char* fname, value *out,  ulong dims_T[3], ulong dims[3]);

/*+++++++++++++++++++++++++++++*/
/*	     Side functions     	 */
/*+++++++++++++++++++++++++++++*/

void set_rank_borders(uint grid[3], uint myrank_arr[3]) {
    data_properties.border[0] = myrank_arr[0] % grid[0] != 0; //left
    data_properties.border[1] = myrank_arr[0] % grid[0] < grid[0] - 1; //right
    data_properties.border[2] = myrank_arr[1] % grid[1] != 0; //top
    data_properties.border[3] = myrank_arr[1] % grid[1] < grid[1] - 1; //bottom
    data_properties.border[4] = myrank_arr[2] % grid[2] != 0; 
    data_properties.border[5] = myrank_arr[2] % grid[2] < grid[2] - 1;
}

void set_attribute_offsets(uint grid[3], uint myrank_arr[3]) {
  ulong		counts[3]     = {0};
  for (int i = 3; i-- ; ) {
    counts[i]  = data_properties.dims_full[i] / grid[i];
    data_properties.offsets[i] = myrank_arr[i]*counts[i];
    if (myrank_arr[i] < data_properties.dims_full[i]%grid[i]) {
      counts[i]++;
      data_properties.offsets[i] += myrank_arr[i];
    } else {
      data_properties.offsets[i] += data_properties.dims_full[i]%grid[i];
    }
    counts[i] += data_properties.offsets[i] - 1;
    if((data_properties.offsets[i] > 1)){
      data_properties.offsets[i]--;
    }
  }
}

void check_connectivity(Arguments* args){
  if (data_properties.dims_process[2] > 1 && args->pixel_connectivity == 4)
    args->pixel_connectivity = 6;
  else if(args->pixel_connectivity != 4 && args->pixel_connectivity != 6 && args->pixel_connectivity != 8 && args->pixel_connectivity != 26){
    if(rank() == 0){
      error("Wrong Connectivity (must be 4 or 8 for 2D, 6 or 26 for 3D)");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
  } else if( args->pixel_connectivity == 8 && data_properties.dims_process[2] > 1){
    if(rank() == 0) warn("Wrong connectivity (8) for data dimension (3D), changing to 26");
    args->pixel_connectivity = 26;
  } else if( args->pixel_connectivity == 26 && data_properties.dims_process[2] == 1){
    if(rank() == 0) warn("Wrong connectivity (26) for data dimension (2D), changing to 8");
    args->pixel_connectivity = 8;
  }
}

FIBITMAP* freeimage_generic_loader(const char* lpszPathName, int flag) {
  FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
  //  check the file signature and deduce its format
  fif = FreeImage_GetFileType(lpszPathName, 0);
  if(fif == FIF_UNKNOWN) {
    // no signature ? try to guess the file format from the file extension
    fif = FreeImage_GetFIFFromFilename(lpszPathName);
  }
  // check that the plugin has reading capabilities ...
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
    FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
    // unless a bad file format, we are done !
    return dib;
  }
  else{
    error("FreeImage couldn't read the input \n");
    MPI_Abort(MPI_COMM_WORLD, 735);
  }
  return NULL;
} /* freeimage_generic_loader */

void get_distributed_borders(value **img, uint grid[3],  int overlap){

  value 	*img_curr ;   
  int 		myrank	      = rank();
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  ulong *dims = data_properties.dims_process;

  MPI_Comm row_comm, col_comm, slice_comm;
  MPI_Comm_split(MPI_COMM_WORLD, myrank_arr[1]+myrank_arr[2]*grid[1], myrank, &row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_arr[0]+myrank_arr[2]*grid[0], myrank, &col_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_2D, myrank, &slice_comm);
  if (!overlap){
    img_curr     = malloc(dims[0]*dims[1]*dims[2]*sizeof(value));
    memcpy(img_curr, *img, dims[0]*dims[1]*dims[2]*sizeof(value));

    ulong new_dims[3] = {dims[0]+data_properties.border[0]+data_properties.border[1], dims[1]+data_properties.border[2]+data_properties.border[3], dims[2]+data_properties.border[4]+data_properties.border[5]};
    free(*img);
    *img = malloc(new_dims[0]*new_dims[1]*new_dims[2] * sizeof(value));
    MPI_Datatype gdep_s, gdep_r, gcol_s, gcol_r, gpla_s, gpla_r, grow_s, grow_r, gang_s, gang_r;
    MPI_Type_vector(dims[1], dims[0], dims[0], mpi_value_type, &gdep_s);
    MPI_Type_commit(&gdep_s);
    MPI_Type_vector(dims[1], dims[0], new_dims[0], mpi_value_type, &gdep_r);
    MPI_Type_commit(&gdep_r);
    MPI_Type_vector(dims[1], 1, dims[0], mpi_value_type, &gcol_s);
    MPI_Type_commit(&gcol_s);
    MPI_Type_vector(dims[1], 1, new_dims[0], mpi_value_type, &gcol_r);
    MPI_Type_commit(&gcol_r);
    MPI_Type_hvector(dims[2], 1, dims[0]*dims[1]*sizeof(value), gcol_s, &gpla_s);
    MPI_Type_commit(&gpla_s);
    MPI_Type_hvector(dims[2], 1, new_dims[0]*new_dims[1]*sizeof(value), gcol_r, &gpla_r);
    MPI_Type_commit(&gpla_r);
    MPI_Type_vector(dims[2], dims[0], dims[0]*dims[1], mpi_value_type, &grow_s);
    MPI_Type_commit(&grow_s);
    MPI_Type_vector(dims[2], dims[0], new_dims[0]*new_dims[1], mpi_value_type, &grow_r);
    MPI_Type_commit(&grow_r);
    MPI_Type_vector(dims[2], 1, dims[0]*dims[1], mpi_value_type, &gang_s);
    MPI_Type_commit(&gang_s);
    MPI_Type_vector(dims[2], 1, new_dims[0]*new_dims[1], mpi_value_type, &gang_r);
    MPI_Type_commit(&gang_r);
    int row_rank, col_rank, dep_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_rank(slice_comm, &dep_rank);

    if(rank() == 2)
      info("ROW %d, Col %d, dep %d", row_rank, col_rank, dep_rank);
    /* Send planes*/
    if(row_rank % 2 == 0){
      if(data_properties.border[1]){
	MPI_Send(img_curr+dims[0]-1, 1, gpla_s, row_rank+1, 1, row_comm);
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+(1+data_properties.border[2])*new_dims[0]-1, 1, gpla_r, row_rank+1, 1, row_comm, NULL);
      }
      if(data_properties.border[0]){
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+data_properties.border[2]*new_dims[0], 1, gpla_r, row_rank-1, 1, row_comm, NULL);
	MPI_Send(img_curr,1, gpla_s, row_rank-1,1,  row_comm);
      }
    } else {
      if(data_properties.border[0]){
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+data_properties.border[2]*new_dims[0], 1, gpla_r, row_rank-1, 1, row_comm, NULL);
	MPI_Send(img_curr,1, gpla_s, row_rank-1,1,  row_comm);
      }
      if(data_properties.border[1]){
	MPI_Send(img_curr+dims[0]-1, 1, gpla_s, row_rank+1, 1, row_comm);
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+(1+data_properties.border[2])*new_dims[0]-1, 1, gpla_r, row_rank+1, 1, row_comm, NULL);
      }
    }

    if(col_rank % 2 == 0){
      if(data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, grow_s, col_rank+1, 1, col_comm);
	MPI_Recv(*img+(data_properties.border[4]+1)*new_dims[0]*new_dims[1]-new_dims[0]+data_properties.border[0], 1, grow_r, col_rank+1, 1, col_comm, NULL);
      }
      if(data_properties.border[2]){
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+data_properties.border[0],1, grow_r, col_rank-1,1,  col_comm, NULL);
	MPI_Send(img_curr,1, grow_s, col_rank-1,1,  col_comm);
      }
    } else {
      if(data_properties.border[2]){
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+data_properties.border[0],1, grow_r, col_rank-1,1,  col_comm, NULL);
	MPI_Send(img_curr,1, grow_s, col_rank-1,1,  col_comm);

      }
      if(data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, grow_s, col_rank+1, 1, col_comm);
	MPI_Recv(*img+(data_properties.border[4]+1)*new_dims[0]*new_dims[1]-new_dims[0]+data_properties.border[0], 1, grow_r, col_rank+1, 1, col_comm, NULL);
      }
    }

    if(dep_rank %2 == 0){
      if(data_properties.border[5]){
	MPI_Send(img_curr+(dims[2]-1)*dims[1]*dims[0], 1, gdep_s, dep_rank+1, 1, slice_comm);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]*data_properties.border[2]+data_properties.border[0], 1, gdep_r, dep_rank+1, 1, slice_comm, NULL);
      }
      if(data_properties.border[4]){
	MPI_Recv(*img+new_dims[0]*data_properties.border[2]+data_properties.border[0],1, gdep_r, dep_rank-1,1,  slice_comm, NULL);
	MPI_Send(img_curr,1,gdep_s, dep_rank-1,1,  slice_comm);
      }
    } else {
      if(data_properties.border[4]){
	MPI_Recv(*img+new_dims[0]*data_properties.border[2]+data_properties.border[0],1, gdep_r, dep_rank-1,1,  slice_comm, NULL);
	MPI_Send(img_curr,1,gdep_s, dep_rank-1,1,  slice_comm);
      }
      if(data_properties.border[5]){
	MPI_Send(img_curr+(dims[2]-1)*dims[1]*dims[0], 1, gdep_s, dep_rank+1, 1, slice_comm);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]*data_properties.border[2]+data_properties.border[0], 1, gdep_r, dep_rank+1, 1, slice_comm, NULL);
      }
    }

    if(row_rank % 2 == 0){
      if(data_properties.border[1] && data_properties.border[2]){
	MPI_Send(img_curr+dims[0]-1,1, gang_s, rank()+1-grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+new_dims[0]-1,1, gang_r, rank()+1-grid[0],1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[1] && data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]-1,1, gang_s, rank()+1+grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(1+data_properties.border[4])*new_dims[0]*new_dims[1]-1,1, gang_r, rank()+1+grid[0],1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[0] && data_properties.border[2]){
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1],1, gang_r, rank()-1-grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr,1, gang_s, rank()-1-grid[0],1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[0] && data_properties.border[3]){
	MPI_Recv(*img+(1+data_properties.border[4])*new_dims[0]*new_dims[1]-new_dims[0],1, gang_r, rank()-1+grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+(dims[1]-1)*dims[0],1, gang_s, rank()-1+grid[0],1,  MPI_COMM_WORLD);
      }
    } else {
      if(data_properties.border[0] && data_properties.border[3]){
	MPI_Recv(*img+(1+data_properties.border[4])*new_dims[0]*new_dims[1]-new_dims[0],1, gang_r, rank()-1+grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+(dims[1]-1)*dims[0],1, gang_s, rank()-1+grid[0],1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[0] && data_properties.border[2]){
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1],1, gang_r, rank()-1-grid[0],1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr,1, gang_s, rank()-1-grid[0],1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[1] && data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]-1,1, gang_s, rank()+1+grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(1+data_properties.border[4])*new_dims[0]*new_dims[1]-1,1, gang_r, rank()+1+grid[0],1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[1] && data_properties.border[2]){
	MPI_Send(img_curr+dims[0]-1,1, gang_s, rank()+1-grid[0],1,  MPI_COMM_WORLD);
	MPI_Recv(*img+data_properties.border[4]*new_dims[0]*new_dims[1]+new_dims[0]-1,1, gang_r, rank()+1-grid[0],1,  MPI_COMM_WORLD, NULL);
      }
    }
      
    if(dep_rank %2 == 0){
      if(data_properties.border[4] && data_properties.border[2]){
	MPI_Send(img_curr,dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+data_properties.border[0],dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[4] && data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1),dims[0], mpi_value_type, rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(new_dims[1]-1)*new_dims[0]+data_properties.border[0],dims[0], mpi_value_type,  rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
      }		
      if(data_properties.border[5] && data_properties.border[2]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+data_properties.border[0],dims[0], mpi_value_type,rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1),dims[0], mpi_value_type,  rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[5] && data_properties.border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2]-new_dims[0]+data_properties.border[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*dims[2]-dims[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
      }

    } else {
      if(data_properties.border[5] && data_properties.border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2]-new_dims[0]+data_properties.border[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*dims[2]-dims[0],dims[0], mpi_value_type, rank()+grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
      }	
      if(data_properties.border[5] && data_properties.border[2]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+data_properties.border[0],dims[0], mpi_value_type,rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1),dims[0], mpi_value_type,  rank()+grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[4] && data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1),dims[0], mpi_value_type, rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(new_dims[1]-1)*new_dims[0]+data_properties.border[0],dims[0], mpi_value_type,  rank()-grid[0]*(grid[1]-1),1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[4] && data_properties.border[2]){
	MPI_Send(img_curr,dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD);
	MPI_Recv(*img+data_properties.border[0],dims[0], mpi_value_type, rank()-grid[0]*(grid[1]+1),1,  MPI_COMM_WORLD, NULL);
      }
    }

    if(dep_rank %2 == 0){
      if(data_properties.border[4] && data_properties.border[0]){
	MPI_Send(img_curr,1, gcol_s, rank()-grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+data_properties.border[2]*new_dims[0],1, gcol_r, rank()-grid[0]*grid[1] -1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[4] && data_properties.border[1]){
	MPI_Send(img_curr+dims[0]-1,1, gcol_s, rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+(1+data_properties.border[2])*new_dims[0]-1,1, gcol_r,  rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[5] && data_properties.border[0]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, gcol_s,  rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+data_properties.border[2]*new_dims[0],1, gcol_r,rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[5] && data_properties.border[1]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+(data_properties.border[2]+1)*new_dims[0]-1,1, gcol_r, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
      }
    } else {
      if(data_properties.border[5] && data_properties.border[1]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+(data_properties.border[2]+1)*new_dims[0]-1,1, gcol_r, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[5] && data_properties.border[0]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]*(new_dims[2]-1)+data_properties.border[2]*new_dims[0],1, gcol_r,rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, gcol_s,  rank()+grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[4] && data_properties.border[1]){
	MPI_Recv(*img+(1+data_properties.border[2])*new_dims[0]-1,1, gcol_r,  rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]-1,1, gcol_s, rank()-grid[0]*grid[1]+1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[4] && data_properties.border[0]){
	MPI_Recv(*img+data_properties.border[2]*new_dims[0],1, gcol_r, rank()-grid[0]*grid[1] -1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr,1, gcol_s, rank()-grid[0]*grid[1]-1,1,  MPI_COMM_WORLD);
      }
    }

    if(dep_rank %2 == 0){
      if(data_properties.border[4] && data_properties.border[0] && data_properties.border[2]){
	MPI_Send(img_curr, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[4] && data_properties.border[1] && data_properties.border[2]){
	MPI_Send(img_curr+dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[4] && data_properties.border[0] && data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*(new_dims[1]-1),1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[4] && data_properties.border[1] && data_properties.border[3]){
	MPI_Send(img_curr+dims[1]*dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }
	
      if(data_properties.border[5] && data_properties.border[0] && data_properties.border[2]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1),1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[5] && data_properties.border[1] && data_properties.border[2]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }

      if(data_properties.border[5] && data_properties.border[0] && data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2] - new_dims[0],1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
      }
      if(data_properties.border[5] && data_properties.border[1] && data_properties.border[3]){
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[1]*dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
      }
    } else{
      if(data_properties.border[5] && data_properties.border[1] && data_properties.border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[1]*dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[5] && data_properties.border[0] && data_properties.border[3]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*new_dims[2] - new_dims[0],1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[5] && data_properties.border[1] && data_properties.border[2]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1)+new_dims[0]-1,1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1)+dims[0]-1, 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[5] && data_properties.border[0] && data_properties.border[2]){
	MPI_Recv(*img+new_dims[0]*new_dims[1]*(new_dims[2]-1),1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*dims[1]*(dims[2]-1), 1, mpi_value_type, rank()+grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[4] && data_properties.border[1] && data_properties.border[3]){
	MPI_Recv(*img+new_dims[1]*new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[1]*dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[4] && data_properties.border[0] && data_properties.border[3]){
	MPI_Recv(*img+new_dims[0]*(new_dims[1]-1),1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]*(dims[1]-1), 1, mpi_value_type, rank()-grid[0]*grid[1]+grid[0]-1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[4] && data_properties.border[1] && data_properties.border[2]){
	MPI_Recv(*img+new_dims[0]-1,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr+dims[0]-1, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]+1,1,  MPI_COMM_WORLD);
      }
      if(data_properties.border[4] && data_properties.border[0] && data_properties.border[2]){
	MPI_Recv(*img,1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD, NULL);
	MPI_Send(img_curr, 1, mpi_value_type, rank()-grid[0]*grid[1]-grid[0]-1,1,  MPI_COMM_WORLD);
      }
    }

    ulong offset = data_properties.border[4]*new_dims[0]*new_dims[1] + data_properties.border[2]*new_dims[0]+ data_properties.border[0];
    
    for(ulong i = 0; i < dims[0]*dims[1]*dims[2]; i++){
      ulong x = i % dims[0];
      ulong y = i %(dims[0]*dims[1]) / dims[0];
      ulong z = i / (dims[0]*dims[1]);
      (*img)[offset+z*(new_dims[0]*new_dims[1])+y*new_dims[0]+x] = img_curr[i];
    }

    data_properties.dims_process[0] = new_dims[0];
    data_properties.dims_process[1] = new_dims[1];
    data_properties.dims_process[2] = new_dims[2];
    free(img_curr);
    MPI_Type_free( &grow_r );
    MPI_Type_free( &grow_s );
    MPI_Type_free( &gcol_r );
    MPI_Type_free( &gcol_s );
    MPI_Type_free( &gpla_r );
    MPI_Type_free( &gpla_s );
    MPI_Type_free( &gdep_r );
    MPI_Type_free( &gdep_s );

  }

   
  ulong send = data_properties.dims_process[0]-data_properties.border[0]-data_properties.border[1];
  MPI_Allreduce(&send, data_properties.dims_full, 1, MPI_UNSIGNED_LONG, MPI_SUM, row_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_arr[0]+myrank_arr[2]*grid[0], myrank, &col_comm);
  send = data_properties.dims_process[1]-data_properties.border[2]-data_properties.border[3];
  MPI_Allreduce(&send, data_properties.dims_full+1, 1, MPI_UNSIGNED_LONG, MPI_SUM, col_comm);
  MPI_Comm_split(MPI_COMM_WORLD, myrank_2D, myrank, &slice_comm);
  send = data_properties.dims_process[2]-data_properties.border[4]-data_properties.border[5];
  MPI_Allreduce(&send, data_properties.dims_full+2, 1, MPI_UNSIGNED_LONG, MPI_SUM, slice_comm);
  MPI_Comm_free(&row_comm);
  MPI_Comm_free(&col_comm);
  MPI_Comm_free(&slice_comm);

}

void print_image_parameters(Arguments *args){
    info("\033[1m*** INPUT IMAGE PARAMETERS ***\033[0m");
    info("\033[1mImage type:\033[0m %s", args->image_options);
    info("\033[1mImage total dimensions:\033[0m  %lu, %lu, %lu", data_properties.dims_full[0],data_properties.dims_full[1],data_properties.dims_full[2]);
    info("\033[1mImage dimensions on current MPI process:\033[0m  %lu, %lu, %lu", data_properties.dims_process[0],data_properties.dims_process[1],data_properties.dims_process[2]);
    info("\033[1mImage bit per pixel:\033[0m  %d", args->bpp);
    if(FLOAT_TYPE == 0 && args->bpp < 0){
      error("The image is single or double floating point but the FLOAT_TYPE flag is 0. That won't work. Make sure to change the flag in types.h and re-compile");
      MPI_Abort(MPI_COMM_WORLD, 001);
    }
    if(FLOAT_TYPE == 1 && args->bpp > 0){
      warn("The image is %d bit-per-pixel, and the FLOAT_TYPE flag is 1 (assume the data is floating point). That will still work, but you could gain some decent memory if you set the proper flag in the file types.h",args->bpp);
    }
    info("\033[1mImage pixel connectivity:\033[0m  %d", args->pixel_connectivity);
    info("\033[1m*** ***\033[0m");
    info("");
}

/*+++++++++++++++++++++++++++++*/
/*	     Reading functions   	 */
/*+++++++++++++++++++++++++++++*/


void read_input(Arguments *args,  value **img){

  int 		bitpix;
  uint	 	grid[3]       = {args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]};
  uint	  myrank        = rank();
  uint 		myrank_2D     = myrank % (grid[1]*grid[0]);
  uint 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  // uint 		tile   	      = !strcmp(args->input_file_style,"tile") ?  myrank_arr[0] + (grid[1]-myrank_arr[1]-1)*grid[0] +myrank_arr[2]*grid[1]*grid[0]: myrank;
  uint 		tile   	      = !strcmp(args->input_file_style,"tile") ?  myrank_arr[0] + (myrank_arr[1])*grid[0] +myrank_arr[2]*grid[1]*grid[0]: myrank;

  uint    file_style    = !strcmp(args->input_file_style,"single")? 0: 1;
 

  set_rank_borders(grid, myrank_arr); // set border array

  if(!strcmp(args->image_options, "grayscale")){

    read_grayscale_image(args, img, grid, &bitpix, myrank, myrank_2D, tile,
                         myrank_arr, file_style); 
  }

  else if (!strcmp(args->image_options, "rgb_channels")){
    read_rgb_channels(args, img,  grid, &bitpix, myrank, myrank_2D, tile,
                         myrank_arr, file_style); 
  }

  else if (!strcmp(args->image_options, "lofar")){
    read_lofar_data(args, img, grid, &bitpix,  myrank, myrank_2D,
                         myrank_arr, file_style);
  }

  if(!strcmp(args->input_file_style,"tile"))
    get_distributed_borders(img, grid,  args->tile_overlap); //add an overlap of 1 pixel in each tile

  // Offsets for attributes //
  set_attribute_offsets(grid, myrank_arr);

  check_connectivity(args);

  print_image_parameters(args);
} /* read_input */

void read_grayscale_image(Arguments* args, value** img, 
                          uint grid[3], int* bitpix, uint myrank, uint myrank_2D, uint tile,
                          uint myrank_arr[3], uint file_style) {
    char* fname;
    if (!strcmp(args->input_file_style, "single")) {
        asprintf(&fname, "%s", args->input_name);
    } else if (!strcmp(args->input_file_style, "tile")) {
        asprintf(&fname, "%s-T%dT.%s", args->input_prefix, tile, args->input_type);
    }
    
    if (!strcmp(args->input_type, "h5")) {
        read_hdf5(args, img, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
    } else if (!strcmp(args->input_type, "fits") || !strcmp(args->input_type, "fit")) {
        read_fits(args, img, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
    } else if (!strcmp(args->input_type, "nii")) {
        read_nifti_file(args, img, fname,  bitpix);
    } else {
        read_basic(args, img, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
    }
    
    args->bpp = *bitpix;
    free(fname);
}

void read_rgb_channels(Arguments* args, value** img,
                       uint grid[3], int* bitpix, uint myrank, uint myrank_2D, uint tile,
                       uint myrank_arr[3], uint file_style) {

    char  rgb[3]    = {'R', 'G', 'B'};
    float coef[3]   = {0.2126, 0.7152, 0.0722};

    value* img_curr = NULL;

    for (int i = 0; i < 3; i++) {
        char* fname;
        if (!strcmp(args->input_file_style, "single")) {
            asprintf(&fname, "%s-%c.%s", args->input_prefix, rgb[i], args->input_type);
        } else if (!strcmp(args->input_file_style, "tile")) {
            asprintf(&fname, "%s-%c-T%dT.%s", args->input_prefix, rgb[i], tile, args->input_type);
        }

        if (!strcmp(args->input_type, "h5")) {
            read_hdf5(args, &img_curr, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
        } else if (!strcmp(args->input_type, "fits") || !strcmp(args->input_type, "fit")) {
            read_fits(args, &img_curr, fname,   grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
        } else if (!strcmp(args->input_type, "nii")) {
            read_nifti_file(args, &img_curr, fname,  bitpix);
        } else {
            read_basic(args, &img_curr, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
        }

        if (i == 0) {
            *img = calloc(data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2], sizeof(value));
        }

        #pragma omp parallel for
        for (ulong j = 0; j < data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2]; j++) {
            (*img)[j] += (value) coef[i] * img_curr[j];
        }

        free(img_curr);
        free(fname);
    }

    args->bpp = -32;
}

void read_lofar_data(Arguments* args, value** img, 
                     uint grid[3], int* bitpix, uint myrank, uint myrank_2D,
                     uint myrank_arr[3], uint file_style) {
    if (grid[0] > 1 || grid[1] > 1) {
        error("Take 3D sequence, can't divide it in space");
        MPI_Abort(MPI_COMM_WORLD, 001);
    }

    char* fname = NULL;
    
    strcpy(args->input_type, "tile");

    value* img_curr = NULL;
    int nchannels = 2;
    int slices = 2;
    int nfiles = nchannels * slices;
    int offset[2] = {0};
    
    if (data_properties.border[4]) {
        offset[0]++;
    }
    
    if (data_properties.border[5]) {
        offset[1]++;
    }
    
    if (offset[0]) {
        asprintf(&fname, "test%03d-000%d.fits", slices * myrank - 1, nchannels - 1);
        read_fits(args, &img_curr, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
        
        *img = malloc((nfiles + offset[0] + offset[1]) * data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * sizeof(value));	  
        memcpy(*img, img_curr, data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * sizeof(value));
        free(img_curr);
    }
    
    for (int i = 0; i < slices; i++) {
        for (int j = 0; j < nchannels; j++) {
            asprintf(&fname, "test%03d-000%d.fits", slices * myrank + i, j);
            read_fits(args, &img_curr, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
            
            if (!i && !j && !offset[0]) {
                *img = malloc((nfiles + offset[0] + offset[1]) * data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * sizeof(value));	  
            }
            
            memcpy(*img + data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * (i * nchannels + j + offset[0]), img_curr, data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * sizeof(value));
            free(img_curr);
        }
    }
    
    if (offset[1]) {
        asprintf(&fname, "test%03d-0000.fits", slices * (myrank + 1));
        read_fits(args, &img_curr, fname,  grid, bitpix, myrank, myrank_2D, myrank_arr, file_style);
        memcpy(*img + data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * (nfiles + offset[0]), img_curr, data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * sizeof(value));
        free(img_curr);
    }
    
    data_properties.dims_process[2] = nfiles + offset[0] + offset[1];
    args->bpp = *bitpix; 
    free(fname);
}



void read_fits(Arguments* args, value** img, const char* fname, 
               uint grid[3], int* bitpix, uint myrank, uint myrank_2D, uint myrank_arr[3], uint file_style) {
    fitsfile* fptr;     // FITS file pointer
    int n_dims;         // Number of dimensions
    int type;           // Data type of the FITS image
    int status = 0;     // CFITSIO status value; must be initialized to zero
    long inc[4] = {1, 1, 1, 1}; // Increment for FITS read function
    long naxes[4] = {1, 1, 1, 1}; // Dimensions (4D max) of the FITS image
    long counts[4] = {1, 1, 1, 1}; // Pixels to read in FITS read function
    long offsets[4] = {1, 1, 1, 1}; // Pixel index to start in FITS read function
    int overlap = args->tile_overlap;
    int negative = 0;
    int ttype = sizeof(value) == 1 ? 1 - negative : negative; // Determine type based on value size
    
    debug("Reading FITS Image %s", fname);

    // Determine the FITS data type based on the value type
    if (FLOAT_TYPE) {
        type = TFLOAT;
    } else {
        type = (int)(log2(sizeof(value)) + 1) * 10 + ttype;
    }

    if (!fits_open_file(&fptr, fname, READONLY, &status)) { // Open the FITS file
        if (!fits_get_img_param(fptr, 4, bitpix, &n_dims, naxes, &status)) { // Get image parameters
            if (n_dims > 4 || n_dims == 0 || (n_dims == 4 && naxes[3] > 1)) { // Check dimensionality
                error("Only 2D and 3D images are supported");
                MPI_Abort(MPI_COMM_WORLD, 701);
            } else {
                if (n_dims == 4) {
                    n_dims--;
                }

                if (!file_style) {
                    // Determine dimensions and offsets for distribution
                    if (n_dims == 2 && grid[2] > 1) {
                        error("Requested distribution in depth, but data is 2D");
                        MPI_Abort(MPI_COMM_WORLD, 701);
                    }

                    for (int i = n_dims; i-- ;) {
                        counts[i] = naxes[i] / grid[i];
                        offsets[i] = myrank_arr[i] * counts[i] + 1;
                        if (myrank_arr[i] < naxes[i] % grid[i]) {
                            counts[i]++;
                            offsets[i] += myrank_arr[i];
                        } else {
                            offsets[i] += naxes[i] % grid[i];
                        }

                        data_properties.dims_full[i] = naxes[i];
                        data_properties.dims_process[i] = counts[i];
                        counts[i] += offsets[i] - 1;

                        if ((offsets[i] > 1) && overlap) {
                            offsets[i]--;
                            data_properties.dims_process[i]++;
                        }
                        
                        if (((ulong)counts[i] != data_properties.dims_full[i]) && overlap) {
                            counts[i]++;
                            data_properties.dims_process[i]++;
                        }
                    }
                } else {
                    // Use specified dimensions
                    for (int i = n_dims; i--;) {
                        counts[i] = data_properties.dims_process[i] = naxes[i];
                    }
                }
                
                debug("Fits reading:\n Tile dimensions: %ld by %ld by %ld\n Offsets %ld, %ld, %ld\n Counts %ld, %ld, %ld", data_properties.dims_process[0], data_properties.dims_process[1], data_properties.dims_process[2], offsets[0], offsets[1], offsets[2], counts[0], counts[1], counts[2]);
                *img = malloc(data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2] * sizeof(value));
                fits_read_subset(fptr, type, offsets, counts, inc, NULL, *img, NULL, &status);
            }
        }
        
        fits_close_file(fptr, &status); // Close the FITS file
    } else {
        error("Cannot open the file %s, wrong name? (FITS)", fname);
        MPI_Abort(MPI_COMM_WORLD, 703);
    }
}

void read_hdf5(Arguments *args, value **img, const char *fname, uint grid[3], int *bitpix, uint myrank, uint myrank2D, uint myrank_arr[3], uint file_style){
  
  /* +++++++++++++++++++++++++++ */
  /*      Read file HDF5        */
  /* +++++++++++++++++++++++++++ */
  
  H5T_class_t   t_class;
  hid_t       	file_id, dataset_id, dataspace;  	    
  hsize_t 	    counts[3]     = {1,1,1};
  hsize_t 	    offsets[3]    = {0};
  hsize_t 	    *hdims;
  hsize_t 	    hn_dims;
  herr_t	      err;
  int  		overlap       = args->tile_overlap;
  int 		negative      = 0;
  
  info("Reading HDF5 file %s", fname);

  /*      Open HDF5 file and dataset     */
  
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);	
  if (file_id < 0) {
    error("Could not open file %s, wrong name ? (HDF5) ", fname);
    MPI_Abort(MPI_COMM_WORLD, 707);
  }
  
  dataset_id = H5Dopen2(file_id, args->hdf5_dataset, H5P_DEFAULT);	
  if (dataset_id < 0) {
    error("Could not open dataset %s, wrong name ? (HDF5) ", args->hdf5_dataset);
    MPI_Abort(MPI_COMM_WORLD, 708);
  }

  
  dataspace = H5Dget_space(dataset_id);
  hid_t tid = H5Dget_type(dataset_id);
  t_class   = H5Tget_class(tid);
  int ord   = H5Tget_precision(tid);
  
  if(t_class == H5T_FLOAT)
    *bitpix = -32;
  else {  
    if(ord == 8)
      *bitpix = 8;
    else if(ord == 16)
      *bitpix = 16;
    else 
      *bitpix = 32;
  }
  
  hn_dims = H5Sget_simple_extent_ndims(dataspace);
  
  /*       Get number of dimensions  and dimensions    */
  
  if (hn_dims > 3) {
    error("Only handle 2D or 3D data");
    MPI_Abort(MPI_COMM_WORLD, 709);
  }

  hdims = calloc(hn_dims, sizeof(hsize_t));			
  H5Sget_simple_extent_dims(dataspace, hdims, NULL);

  if(file_style){
    for (hsize_t i = 0; i < hn_dims; i++)
      counts[i] = data_properties.dims_process[hn_dims - i - 1] = hdims[i];

  } else {
    if(hn_dims == 2 && grid[2] > 1){
      error("Ask distribution in depth but data is 2D");
      MPI_Abort(MPI_COMM_WORLD, 701);
    }
    for (hsize_t i = 0; i < hn_dims; i++) {
      counts[i]  = hdims[i]/grid[hn_dims - i - 1];
      offsets[i] = myrank_arr[hn_dims - i - 1] * counts[i];
      if (myrank_arr[hn_dims - i - 1] < (int) (hdims[i]%grid[hn_dims - i - 1])) {
        counts[i]++;
        offsets[i] += myrank_arr[hn_dims - i - 1];
      } else {
	      offsets[i] += (hdims[i]%grid[hn_dims - i - 1]);
      }
      if((offsets[i] > 0) && overlap){
        offsets[i]--;
        counts[i]++;
      }
      if((counts[i] + offsets[i] != hdims[i]) && overlap) {
	      counts[i]++;				       			       
      }
        data_properties.dims_full[hn_dims - i - 1] = hdims[i];
        data_properties.dims_process[hn_dims - i - 1]   = counts[i];
    }

  }
  
  debug("HDF5: tile dimensions: depth %ld, height %ld, width %ld.\n First pixel offset: depth %ld, height %ld, width %ld", counts[0], counts[1], counts[2], offsets[0], offsets[1], offsets[2]);
  
  *img = malloc(counts[0]*counts[1]*counts[2] * sizeof(value));   check_alloc(img,710);
  
  hid_t memory_window = H5Screate_simple(hn_dims, counts, NULL); /* Allocate memory to read the data */
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL); /* Select the data subset */
  if(FLOAT_TYPE)
    err = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 8)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_LONG, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_ULONG, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 4)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_INT, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_UINT, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 2)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_SHORT, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_USHORT, memory_window, dataspace, H5P_DEFAULT, *img);
  else if(sizeof(value) == 1)
    err = negative ? H5Dread(dataset_id, H5T_NATIVE_CHAR, memory_window, dataspace, H5P_DEFAULT, *img):
      H5Dread(dataset_id, H5T_NATIVE_UCHAR, memory_window, dataspace, H5P_DEFAULT, *img);

  if (err < 0){ 
    error("Could not read data from the dataset (HDF5)");
    MPI_Abort(MPI_COMM_WORLD, 711);
  }

  free(hdims);
  
  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);
}

void read_basic(Arguments *args, value **img, const char *fname, uint grid[3], int *bitpix, uint myrank, uint myrank2D, uint myrank_arr[3], uint file_style){
  

  /* +++++++++++++++++++++++++++ */
  /*       Read FreeImage        */
  /* +++++++++++++++++++++++++++ */
  
  debug("Reading image %s with FreeImage library", fname);
  
  ulong 	i   	        = 0;
  FIBITMAP 	*dib 	      = freeimage_generic_loader(fname, 0);  /* Generic Loader from FreeImage */
  ulong 	offsets[2]    = {0};
  int  		overlap       = args->tile_overlap;
  
  if (dib == NULL) {
    error("Can't read the image %s (FreeImage)", fname);
    MPI_Abort(MPI_COMM_WORLD, 717);
  }

  FREE_IMAGE_COLOR_TYPE nc = FreeImage_GetColorType(dib);

  
  if(nc != 1)
    dib =  FreeImage_ConvertToGreyscale(dib);
  *bitpix = FreeImage_GetBPP(dib);				/* Get dynamic range */

  
  if(!file_style){
    data_properties.dims_full[0] = FreeImage_GetWidth(dib);			    
    data_properties.dims_full[1] = FreeImage_GetHeight(dib);
    if(grid[2] > 1){
      error("Ask distributionin depth but data is 2D");
      MPI_Abort(MPI_COMM_WORLD, 701);
    }
    for (i = 0; i < 2; i++){
      data_properties.dims_process[i]  = data_properties.dims_full[i]/grid[i];
      offsets[i] = myrank_arr[i] * data_properties.dims_process[i];
      if (myrank_arr[i] < (int) (data_properties.dims_full[i]%grid[i])) {
        data_properties.dims_process[i]++;
        offsets[i] += myrank_arr[i];
      } else {
      	offsets[i] += (data_properties.dims_full[i]%grid[i]);
      }
      if((offsets[i] > 0) && overlap){
	      offsets[i]--;
	      data_properties.dims_process[i]++;
      }
      if((data_properties.dims_process[i] + offsets[i] != data_properties.dims_full[i]) && overlap) {
	      data_properties.dims_process[i]++;
      }
    }
  } else {
    data_properties.dims_process[0] = FreeImage_GetWidth(dib);				/* Get width */
    data_properties.dims_process[1] = FreeImage_GetHeight(dib);
  }

  *img = malloc(data_properties.dims_process[0]*data_properties.dims_process[1] * sizeof(value));	check_alloc(img,718);

  i = 0;
  if (*bitpix <= 8) {
    for(ulong y = offsets[1]; y < data_properties.dims_process[1] +offsets[1]; y++) {
      BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, y);
      for(ulong x = offsets[0]; x < data_properties.dims_process[0] + offsets[0]; x++,i++){
	      (*img)[i] = (value) bits[x];
      }
    }
    FreeImage_Unload(dib);
    
  } else if (*bitpix <= 16) {
    for(ulong y = offsets[1]; y < data_properties.dims_process[1] + offsets[1]; y++) {
      ushort *bits = (ushort *)FreeImage_GetScanLine(dib, y);
      for(ulong x = offsets[0]; x < data_properties.dims_process[0] + offsets[0]; x++,i++)
	      (*img)[i] = (value) bits[x];
    }
    FreeImage_Unload(dib);
  } else if (*bitpix <= 32) {
    for(ulong y = offsets[1]; y < data_properties.dims_process[1] + offsets[1]; y++) {
      uint *bits = (uint *)FreeImage_GetScanLine(dib, y);
      for(ulong x = offsets[0]; x < data_properties.dims_process[0] + offsets[0]; x++,i++)
	      (*img)[i] = (value) bits[x];
    }
    FreeImage_Unload(dib);
  } else{
    error("Floating point or 64 bits depth are not supported with classical images (Freeimage)");
    MPI_Abort(MPI_COMM_WORLD, 719);
    FreeImage_Unload(dib);
  }
}

void read_nifti_file(Arguments *args, value **img, const char *fname, int *bitpix){
  nifti_1_header hdr;
  FILE *fp;
  int ret,i;
  double total;
  value *data=NULL;


  /********** open and read header */
  fp = fopen(fname,"r");
  if (fp == NULL) {
    fprintf(stderr, "\nError opening header file %s\n",fname);
    exit(1);
  }
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  if (ret != 1) {
    fprintf(stderr, "\nError reading header file %s\n",fname);
    exit(1);
  }
  fclose(fp);


  /********** print a little header information */
  fprintf(stderr, "\n%s header information:",fname);
  fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
  fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
  fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
  fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
  fprintf(stderr, "\n");

  *bitpix = hdr.bitpix;
  
  /********** open the datafile, jump to data offset */
  fp = fopen(fname,"r");
  if (fp == NULL) {
    error("Error opening data file %s",fname);
    exit(1);
  }

  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
  if (ret != 0) {
    error("Error doing fseek() to %ld in data file %s",(long)(hdr.vox_offset), fname);
    exit(1);
  }


/********** allocate buffer and read first 3D volume from data file */
  *img = (value *) malloc(sizeof(value) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
  if (*img == NULL) {
    error("Error allocating data buffer for %s",fname);
    exit(1);
  }
  if(*bitpix < 0)
    ret = fread(*img, sizeof(float), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  else if(*bitpix <=8)
    ret = fread(*img, sizeof(ubyte), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  else if(*bitpix <=16)
    ret = fread(*img, sizeof(ushort), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  else if(*bitpix <=32)
    ret = fread(*img, sizeof(uint), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
  if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
    error("Error reading volume 1 from %s (%d)");
    exit(1);
  }
  fclose(fp);


  /********** scale the data buffer  */
 if (hdr.scl_slope != 0) {
   for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
     (*img)[i] = ((*img)[i] * hdr.scl_slope) + hdr.scl_inter;
 }


/********** print mean of data */
 total = 0;
 for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
   total += (*img)[i];
 total /= (hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
 fprintf(stderr, "\nMean of volume 1 in %s is %.3f\n",fname,total);

 data_properties.dims_process[0] = data_properties.dims_full[0] =  hdr.dim[1];
 data_properties.dims_process[1] = data_properties.dims_full[1] =  hdr.dim[2];
 data_properties.dims_process[2] = data_properties.dims_full[2] =  hdr.dim[3];

 //return(data);
}


/*++++++++++++++++++++++++++++++*/
/*	     Functions Write      	*/
/*++++++++++++++++++++++++++++++*/

// void write_output(Arguments* args, value* img, value* outOrig, value* outDH, value* outScale, double *spectrum,
//                   Operation* ope_cur, ulong dims_T[3], ulong dims[3], bool border[6], LambdaVec *lvec_1, LambdaVec *lvec_2) {
// void write_output(Arguments* args, Node* tree, Output *output, Operation* ope_cur, ulong dims_T[3], ulong dims[3], bool border[6], LambdaVec *lvec_1, LambdaVec *lvec_2) {
//     int grid[3] = {args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]};
//     int myrank = rank();
//     int myrank_2D = myrank % (grid[1] * grid[0]);
//     int myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1] * grid[0])}; 
//     int tile = (!strcmp(args->input_file_style, "single") || !strcmp(args->output_file_style, "single")) ?
//                myrank_arr[0] + myrank_arr[1] * grid[0] + myrank_arr[2] * grid[2] * grid[1] : myrank;

//     //debug("%d, %d", tile, myrank);
//     if (!strcmp(ope_cur->name, "filter") || !strcmp(ope_cur->name, "extract")) {
//         write_filtered(args, output->filtered, ope_cur, dims_T, dims, border, tile);
//     } else if (!strcmp(ope_cur->name, "csl")) {
//         //write_dap_csv(args, outOrig, outDH, outScale, ope_cur, dims_T, dims, border);
//         write_dap_csv(args, output->contrast, output->scale, output->luminance, ope_cur, dims_T, dims, border);
//     } else if (!strcmp(ope_cur->name, "pattern")) {
//         write_pattern_spectra(args, ope_cur, NULL, lvec_1);//spectrum
//     } else if (!strcmp(ope_cur->name, "pattern2D")) {
//         write_pattern_spectra2d(args,  ope_cur, NULL, lvec_1, lvec_2);//spectrum
//     } else if(!strcmp(ope_cur->name, "check")){
//         write_check_files(args, output->filtered, ope_cur, dims_T, dims, border, tile);
//     } else {
//         write_no_operation(args, tree->gval, NULL, dims_T, dims, border, tile);
//     }
// } // write_output

void write_no_operation(Arguments* args, value* img, Operation* ope_cur) {
  char* fname_out;
  int 	outfile  = !strcmp(args->output_file_style, "single")? 0 : 1;
  ulong *dims_T = data_properties.dims_full;
  ulong *dims   = data_properties.dims_process;
  info("No operation requested, re-writing the original pixel values");

  if (!strcmp(args->output_file_style, "single")) {
      asprintf(&fname_out, "%s.%s", args->output_prefix, args->output_type);
  } else {
      asprintf(&fname_out, "%s-T%dT.%s", args->output_prefix, rank(), args->output_type);
  } 

  if (!strcmp(args->output_type, "h5")) {
      char* dataset_out = "disccofan";
      write_hdf5(args, fname_out, dataset_out, img, dims_T, dims);
  } else if (!strcmp(args->output_type, "fits") || !strcmp(args->output_type, "fit")) {
      write_fits(args, fname_out, img, dims_T, dims);      
  } else {
    if(dims[2] > 1 || (dims_T[2] > 1 && !outfile)){
      if(rank() == 0) warn("FreeImage does not handle writing 3D images, using FITS instead");
      outfile == 0 ? asprintf(&fname_out, "%s.fits", args->output_prefix):
        asprintf(&fname_out, "%s-T%dT.fits", args->output_prefix, rank());
      write_fits(args, fname_out, img, dims_T, dims);      
    } else if(!outfile && np() > 1){
      if(rank() == 0) warn("FreeImage does not handle combining tiles into a single image, using FITS instead");
      asprintf(&fname_out, "%s.fits", args->output_prefix);
      write_fits(args, fname_out, img, dims_T, dims);      
    } else if(args->bpp < 0){
      if(rank() == 0) warn("We do not handle floating point writing with FreeImage, using FITS instead");
      outfile == 0 ? asprintf(&fname_out, "%s.fits", args->output_prefix):
        asprintf(&fname_out, "%s-T%dT.fits", args->output_prefix, rank());
      write_fits(args, fname_out, img, dims_T, dims);      
    }  else
        write_basic(args, fname_out, img, dims_T, dims);   
  }

  free(fname_out);
}

void write_check_files(Arguments* args, value* img, Operation* ope_cur) {
  char* fname_out;
  ulong *dims_T = data_properties.dims_full;
  ulong *dims   = data_properties.dims_process;
  info("Check file automatically written as fits");

  if (!strcmp(args->output_file_style, "single")) {
      asprintf(&fname_out, "%s_%s_%s.fits", args->output_prefix, ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name);
  } else {
      asprintf(&fname_out, "%s-T%dT_%s_%s.fits", args->output_prefix, rank(), ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name);
  }
  int old_bpp = args->bpp;
  args->bpp = -32;
  write_fits(args, fname_out, img, dims_T, dims);      
  args->bpp = old_bpp;

  free(fname_out);
}

void write_filtered(Arguments* args, value* img, Operation* ope_cur) {
  char* fname_out;
  int 	outfile  = !strcmp(args->output_file_style, "single")? 0 : 1;
  ulong *dims_T = data_properties.dims_full;
  ulong *dims   = data_properties.dims_process;
  if (!strcmp(args->output_file_style, "single")) {
      asprintf(&fname_out, "%s_%s_%s_%.3f.%s", args->output_prefix, ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->lambda, args->output_type);
  } else {
      asprintf(&fname_out, "%s-T%dT_%s_%s_%.3f.%s", args->output_prefix, rank(), ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->lambda, args->output_type);
  }

  if (!strcmp(args->output_type, "h5")) {
      char* dataset_out = "disccofan";
      write_hdf5(args, fname_out, dataset_out, img, dims_T, dims);
  } else if (!strcmp(args->output_type, "fits") || !strcmp(args->output_type, "fit")) {
      write_fits(args, fname_out, img, dims_T, dims);      
  } else {
         if(dims[2] > 1 || (dims_T[2] > 1 && !outfile)){
      if(rank() == 0) warn("FreeImage does not handle writing 3D images, using FITS instead");
      outfile == 0 ? asprintf(&fname_out, "%s_%s_%s_%.3f.fits", args->output_prefix, ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->lambda):
        asprintf(&fname_out, "%s-T%dT_%s_%s_%.3f.fits", args->output_prefix, rank(), ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->lambda);
      write_fits(args, fname_out, img, dims_T, dims);      
    } else if(!outfile && np() > 1){
      if(rank() == 0) warn("FreeImage does not handle combining tiles into a single image, using FITS instead");
      asprintf(&fname_out, "%s_%s_%s_%.3f.fits", args->output_prefix, ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->lambda);
      write_fits(args, fname_out, img, dims_T, dims);      
    } else if(args->bpp < 0){
      if(rank() == 0) warn("We do not handle floating point writing with FreeImage, using FITS instead");
      outfile == 0 ? asprintf(&fname_out, "%s_%s_%s_%.3f.fits", args->output_prefix, ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->lambda):
        asprintf(&fname_out, "%s-T%dT_%s_%s_%.3f.fits",args->output_prefix, rank(), ope_cur->name,
                AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->lambda);
      write_fits(args, fname_out, img, dims_T, dims);      
    }  else
        write_basic(args, fname_out, img, dims_T, dims);   
  }

  free(fname_out);
}

void write_dmp(Arguments *args, value *out, int num_lambdas, Operation *ope_cur) {
  // DMP profile

  info("Writing DMP profile as fits");

  int n_dims;
  ulong dims_T[4];
  dims_T[0] =  data_properties.dims_full[0];
  dims_T[1] =  data_properties.dims_full[1];
  dims_T[2] =  data_properties.dims_full[2];
  ulong dims[4];
  dims[0] =  data_properties.dims_process[0];
  dims[1] =  data_properties.dims_process[1];
  dims[2] =  data_properties.dims_process[2];
  if(dims_T[2]>1){
    dims[3] =  num_lambdas;
    dims_T[3] =  num_lambdas;
    n_dims = 4;
  } else {
    n_dims = 3;
    dims[2] =  num_lambdas;
    dims_T[2] = num_lambdas;
  }

  char		*fname_out;

  int 	outfile  = !strcmp(args->output_file_style, "single")? 0 : 1;
  if (!strcmp(args->output_file_style, "single")) {
      asprintf(&fname_out, "%s_%s_%s_%s.fits", args->output_prefix, ope_cur->name, args->tree_type,
                AttribsArray[ope_cur->attribute_idx].short_name);
  } else {
      asprintf(&fname_out, "%s-T%dT_%s_%s_%s.fits", args->output_prefix, rank(), ope_cur->name, args->tree_type,
                AttribsArray[ope_cur->attribute_idx].short_name);
  }

  int	 	  grid[3]       = {args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]};
  int 		myrank        = rank();
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};

  fitsfile 	*outfptr; 	/* FITS file pointers */
  int 		status = 0; 		/* CFITSIO status value MUST be initialized to zero! */
  int 		message;
  int 		type;
  long 		naxes[4]      = {1,1,1,1};      
  long 		counts[4]     = {1,1,1,1};
  long 		offsets[4]    = {1,1,1,1};
  int 		negative      =  0;
  // int 		bitpix        = -32;
  // int     ttype         = sizeof(value) == 1 ? 1-negative: negative;
  char 		str1[100];
  strcpy(str1, "!");   		/* '!' symbol makes the output image, if existing, to be overwritten */

  if(outfile) {
    for (int i = 0; i < n_dims-1; i++) 
      naxes[i] = counts[i] = dims[i]-data_properties.border[2*i]-data_properties.border[2*i+1];     //-data_properties.border[2*i]-data_properties.border[2*i+1]
    naxes[n_dims-1] = counts[n_dims-1] = dims[n_dims-1];
  } else {  
    for (int i = 0; i < n_dims-1; i++) {
      naxes[i]   = dims_T[i];
      counts[i]  = naxes[i] / grid[i];
      offsets[i] = myrank_arr[i]*counts[i] + 1;
      if (myrank_arr[i] < naxes[i]%grid[i]) {
        counts[i]++;
        offsets[i] += myrank_arr[i];
      } else {
	      offsets[i] += naxes[i]%grid[i];
      }
      counts[i] += offsets[i] - 1;
    }
    naxes[n_dims-1] = counts[n_dims-1] = dims_T[n_dims-1];
  }

  type = TFLOAT;
  out = (float *) out;

  if ((myrank == 0) || outfile) {
    if (!fits_create_file(&outfptr, strcat(str1, fname_out), &status)) {
	      fits_create_img(outfptr, FLOAT_IMG, n_dims, naxes, &status);

    } else {
      error("Cannot create the file");
      MPI_Abort(MPI_COMM_WORLD, 730);
    }
  } else {
    MPI_Recv(&message, 1, MPI_INT, myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    fits_open_file(&outfptr, fname_out, READWRITE, &status); // open input images
  }

  ulong off_out   = data_properties.border[0] + data_properties.border[2]*dims[0] + data_properties.border[4]*dims[0]*dims[1];
  ulong off_save  = offsets[1];
  ulong s	  = 0;
  counts[1] = off_save;
  counts[2] = offsets[2];

    if (n_dims == 4) {
        ulong off_save_2  = offsets[2];
        counts[3] = offsets[3];
        for (ulong k = 0; k < dims[3]; k++) {
            // offsets[3] = counts[3] = k + 1;  // Adjust for 4th dimension
            offsets[2] = counts[2] = off_save_2;
            for (ulong i = 0; i < dims[2] - data_properties.border[5] - data_properties.border[4]; i++) {
                offsets[1] = off_save;
                counts[1] = off_save;
                for (ulong j = 0; j < dims[1] - data_properties.border[3] - data_properties.border[2]; j++, s++) {
                    fits_write_subset(outfptr, type, offsets, counts, out + off_out + (s + i * (data_properties.border[2] + data_properties.border[3])) * dims[0], &status);
                    offsets[1]++;
                    counts[1]++;
                }
              counts[2]++;
              offsets[2]++;
            }
            offsets[3]++;
            counts[3]++;
        }
    } else {
        for (ulong i = 0; i < dims[2] - data_properties.border[5] - data_properties.border[4]; i++) {
            offsets[1] = off_save;
            counts[1] = off_save;
            for (ulong j = 0; j < dims[1] - data_properties.border[3] - data_properties.border[2]; j++, s++) {
                fits_write_subset(outfptr, type, offsets, counts, out + off_out + (s + i * (data_properties.border[2] + data_properties.border[3])) * dims[0], &status);
                offsets[1]++;
                counts[1]++;
            }
            counts[2]++;
            offsets[2]++;
        }
    }

  fits_close_file(outfptr, &status);
  if((myrank < np()-1) && !outfile)
    MPI_Send(&status, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);

  if(np()==1 || !outfile)
    info("File %s successfully written", fname_out);
  else
    info("Files (%s on process 0) successfully written", fname_out);
  free(fname_out);

  MPI_Barrier(MPI_COMM_WORLD);
  
} /* write_dmp */

void write_csl(Arguments *args, value *contrast, value *scale, value *luminance, Operation *ope_cur) {
  // CSL segmentation / DAP profile  
  ulong *dims_T = data_properties.dims_full;
  ulong *dims   = data_properties.dims_process;
  char		*fname_out;
  char 		*dataset_out;
  int	 	  grid[3]       = {args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]};
  int 		myrank     = rank();
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};

  char 		step[3]    = {'C','S','L'};

  value 	*out;
  for(int i = 0; i<3; i++ ){
    if (!strcmp(args->output_file_style, "single"))
      asprintf(&fname_out, "%s-%c_%s_%.2f.%s", args->output_prefix, step[i], AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->imscale, args->output_type);
    else{
      asprintf(&fname_out, "%s-%c-T%dT_%s_%.2f.%s", args->output_prefix, step[i], rank(), AttribsArray[ope_cur->attribute_idx].short_name, ope_cur->imscale,  args->output_type);
    }
    
    if(i)
      out = i == 1 ? scale : luminance ;
    else
      out = contrast;

    if (!strcmp(args->output_type, "h5")) {
      asprintf(&dataset_out, "%s", "disccofan");
      write_hdf5(args, fname_out, dataset_out, out,  dims_T, dims);     
    } else if (!strcmp(args->output_type, "fits") || !strcmp(args->output_type, "fit"))
      write_fits(args, fname_out, out,  dims_T, dims);
    else 
      write_basic(args, fname_out, out, dims_T, dims);   
    
    free(fname_out);

    MPI_Barrier(MPI_COMM_WORLD);
  }
} /* write_differential */

void write_pattern_spectra(Arguments *args, Operation *ope_cur, double* spectrum, LambdaVec *lvec){
  char *filename;
  int numscales = lvec->num_lambdas;
  asprintf(&filename, "%s_PS_1D_%s.txt", args->output_prefix,AttribsArray[ope_cur->attribute_idx].short_name);
	       
  FILE *f = fopen(filename, "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }

  fprintf(f, "lambdas \t bin_label \t spectra \n");

  // Allocate buffer for each thread
  int num_threads = omp_get_max_threads();
  char **buffers = malloc(num_threads * sizeof(char *));
  size_t *buffer_sizes = calloc(num_threads, sizeof(size_t));
  for (int i = 0; i < num_threads; i++) {
      buffers[i] = NULL;
  }

  #pragma omp parallel
  {
      int tid = omp_get_thread_num();
      FILE *thread_f = open_memstream(&buffers[tid], &buffer_sizes[tid]);

      #pragma omp for schedule(static)
      for (ulong i = 0; i < (ulong) numscales; i++){ 
        fprintf(thread_f, "%f \t %lf \n", lvec->lambdas[i], spectrum[i]);
      }

      fclose(thread_f);
  }

  // Write buffers to the file sequentially
  for (int i = 0; i < num_threads; i++) {
      fwrite(buffers[i], 1, buffer_sizes[i], f);
      free(buffers[i]);
  }
  free(buffers);
  free(buffer_sizes);

  fclose(f);
    
  info("1D Pattern spectrum file closed (%s)", filename);
  free(filename);
}

void write_pattern_spectra2d(Arguments *args, Operation *ope_cur, double* spectrum, LambdaVec *lvec_attr1, LambdaVec *lvec_attr2) {
    char *filename;
    uint numscales_attr1 = lvec_attr1->num_lambdas;
    uint numscales_attr2 = lvec_attr2->num_lambdas;

    asprintf(&filename, "%s_PS_2D_%s_%s.txt", args->output_prefix, AttribsArray[ope_cur->attribute_idx].short_name, AttribsArray[ope_cur->attribute_idx_2D].short_name);

    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(f, "lambdas_attr1 \t lambda_attr2 \t spectra \n");

    // Allocate buffer for each thread
    int num_threads = omp_get_max_threads();
    char **buffers = malloc(num_threads * sizeof(char *));
    size_t *buffer_sizes = calloc(num_threads, sizeof(size_t));
    for (int i = 0; i < num_threads; i++) {
        buffers[i] = NULL;
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        FILE *thread_f = open_memstream(&buffers[tid], &buffer_sizes[tid]);

        #pragma omp for schedule(static)
        for (ulong i = 0; i < (ulong)(numscales_attr1 * numscales_attr2); i++) {
            uint lambda1_idx = i % numscales_attr1;
            uint lambda2_idx = i / numscales_attr1;
            fprintf(thread_f, "%f \t %f \t %lf \n", lvec_attr1->lambdas[lambda1_idx], lvec_attr2->lambdas[lambda2_idx], spectrum[i]);
        }

        fclose(thread_f);
    }

    // Write buffers to the file sequentially
    for (int i = 0; i < num_threads; i++) {
        fwrite(buffers[i], 1, buffer_sizes[i], f);
        free(buffers[i]);
    }
    free(buffers);
    free(buffer_sizes);

    fclose(f);
    info("2D Pattern spectrum file closed (%s)", filename);
    free(filename);
}


void write_tree(Arguments *args, Node *tree, Operation *ope_cur) {
    char *filename;
    int myrank = rank();
    ulong width  = data_properties.dims_process[0];
    ulong height = data_properties.dims_process[1];
    ulong depth  = data_properties.dims_process[2];
    ulong size2D = width * height;
    bool *border = data_properties.border;
    int     message;
    FILE *f;
    asprintf(&filename, "%s_tree.txt", args->output_prefix);
    int group = AttribsArray[ope_cur->attribute_idx].group;

    if (rank() == 0)
    {
      f = fopen(filename, "w");
      if (f == NULL)
      {
        printf("Error opening file!\n");
        exit(1);
      }
      // Write the header
      fprintf(f, "index,x,y,z,gval,parent,flux,area,");

      if (group == 1)
        fprintf(f, "area_enc_rect,diag_enc_rect,x_extent,y_extent,z_extent");
      else if (group == 2)
        fprintf(f, "x_mean,y_mean,z_mean");
      else if (group == 3)
        fprintf(f, "x_wmean,y_wmean,z_wmean");
      else if (group == 4)
        fprintf(f, "x_mean,y_mean,z_mean,x_wmean,y_wmean,z_wmean,elong,flat,spars,ncomp");

      fprintf(f, "\n");
    }
    else
    {
      MPI_Recv(&message, 1, MPI_INT, rank() - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      f = fopen(filename, "a");
      if (f == NULL)
      {
        printf("Error opening file!\n");
        exit(1);
      }
    }

    // Allocate buffer for each thread to store its output
    int num_threads = omp_get_max_threads();
    char **buffers = malloc(num_threads * sizeof(char *));
    size_t *buffer_sizes = calloc(num_threads, sizeof(size_t));
    for (int i = 0; i < num_threads; i++) {
        buffers[i] = NULL;
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        FILE *thread_f = open_memstream(&buffers[tid], &buffer_sizes[tid]);

        #pragma omp for schedule(static)
        for (ulong i = 0; i < tree->size_curr; i++) {
            if (!is_levelroot(tree, i))
              continue;

            ulong x, y, z, img_index;
            x = (i % size2D) % width;
            y = (i % size2D) / width;
            z = i / size2D;
            bool border_case = (((x == 0 && border[0]) || (x == width-1  && border[1])  ||
            (y == 0 && border[2]) || (y == height-1 && border[3])  ||
            (z == 0 && border[4]) || (z == depth-1  && border[5])) || (i >= size2D*depth));
            if(border_case)
              continue;

            long size = -1;
            float dh = 0;
            x += data_properties.offsets[0];
            y += data_properties.offsets[1];
            z += data_properties.offsets[2];
            img_index = x + y * data_properties.dims_full[0] + z * data_properties.dims_full[0] * data_properties.dims_full[1];


            size = (long)(*AttribsArray[ope_cur->attribute_idx].area)(tree->attribute + i * tree->size_attr);
            dh = tree->parent[i] != BOTTOM ? (tree->gval[i] - tree->gval[tree->parent[i]]) : tree->gval[i];

            fprintf(thread_f, "%ld,%ld,%ld,%ld,%f,%ld,%f,%ld,", img_index, x, y, z, tree->gval[i], get_levelroot(tree, tree->parent[i]), dh * (float)size, size);

            if (group == 1) {
                float encRec = (float)(*AttribsArray[1].attribute)(tree->attribute + i * tree->size_attr);
                float diagRec = (float)(*AttribsArray[2].attribute)(tree->attribute + i * tree->size_attr);
                float x_ext = (float)(*AttribsArray[3].attribute)(tree->attribute + i * tree->size_attr);
                float y_ext = (float)(*AttribsArray[4].attribute)(tree->attribute + i * tree->size_attr);
                float z_ext = (float)(*AttribsArray[5].attribute)(tree->attribute + i * tree->size_attr);
                fprintf(thread_f, "%.2f,%.2f,%.2f,%.2f,%.2f", encRec, diagRec, x_ext, y_ext, z_ext);
            }
            else if (group == 2) {
                float x_mean = (float)(*AttribsArray[6].attribute)(tree->attribute + i * tree->size_attr);
                float y_mean = (float)(*AttribsArray[7].attribute)(tree->attribute + i * tree->size_attr);
                float z_mean = (float)(*AttribsArray[8].attribute)(tree->attribute + i * tree->size_attr);
                fprintf(thread_f, "%.2f,%.2f,%.2f", x_mean, y_mean, z_mean);
            }
            else if (group == 3) {
                float x_wmean = (float)(*AttribsArray[9].attribute)(tree->attribute + i * tree->size_attr);
                float y_wmean = (float)(*AttribsArray[10].attribute)(tree->attribute + i * tree->size_attr);
                float z_wmean = (float)(*AttribsArray[11].attribute)(tree->attribute + i * tree->size_attr);
                fprintf(thread_f, "%.2f,%.2f,%.2f", x_wmean, y_wmean, z_wmean);
            }
            else if (group == 4) {
                double *attrArr = (double *)inertia_attribute_arr(tree->attribute + i * tree->size_attr);
                float elong = (float)attrArr[9];
                float flat = (float)attrArr[10];
                float spars = (float)attrArr[11];
                float ncomp = (float)attrArr[12];
                float c_x = (float)attrArr[3];
                float c_y = (float)attrArr[4];
                float c_z = (float)attrArr[5];
                float c_xw = (float)attrArr[6];
                float c_yw = (float)attrArr[7];
                float c_zw = (float)attrArr[8];
                fprintf(thread_f, "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.4f ", c_x, c_y, c_z, c_xw, c_yw, c_zw, elong, flat, spars, ncomp);
            }
            fprintf(thread_f, "\n");
        }
        fclose(thread_f);
    }

    // Write buffers to the file sequentially
    for (int i = 0; i < num_threads; i++) {
        fwrite(buffers[i], 1, buffer_sizes[i], f);
        free(buffers[i]);
    }
    free(buffers);
    free(buffer_sizes);
    fclose(f);

    if((rank() < np()-1))
      MPI_Send(&myrank, 1, MPI_INT, rank()+1, 1, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    
    info("Component tree file closed (%s)", filename);
    free(filename);
}

void write_hdf5(Arguments *args, const char* fname, const char *dataset_out, value *out, ulong dims_T[3], ulong dims[3]) {
  
  /* +++++++++++++++++++++++++++ */
  /*     HDF5 Write Function     */
  /* +++++++++++++++++++++++++++ */  
  hid_t       	file_id, dataset_id, dataspace, type;  /* identifiers */
  hsize_t	  hdims[3]      = {1,1,1};
  hsize_t 	counts[3]     = {1,1,1};
  hsize_t 	offsets[3]    = {0};

  int     message;
  int 		grid[3]       = {args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]};
  int 		myrank        = rank();		
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  int 		n_dims 	      = dims[2] > 1 ? 3 : 2;
  int 		negative      = 0;
  int 		outfile       = !strcmp(args->output_file_style, "single")? 0 : 1;
  int 		bitpix        = args->bpp;
  
  debug("Writing HDF5 file %s", fname);

  if(bitpix < 0)
    type = H5T_NATIVE_FLOAT;
  else if(bitpix <= 8)
    type =  negative ? H5T_NATIVE_CHAR: H5T_NATIVE_UCHAR;
  else if(bitpix <= 16)
    type =  negative ? H5T_NATIVE_SHORT: H5T_NATIVE_USHORT;
  else if(bitpix <= 32)
    type =  negative ? H5T_NATIVE_INT: H5T_NATIVE_UINT;
  else{
    error("Non supported data type");
    MPI_Abort(MPI_COMM_WORLD, 726);
  }

  if(outfile){
    for (int i = 0; i < n_dims; i++){
      counts[i] = dims[n_dims - i - 1]-data_properties.border[2*(n_dims - i -1)]-data_properties.border[2*(n_dims - i -1)+1];
      hdims[i]   = dims[n_dims - 1 - i]-data_properties.border[2*(n_dims - i -1)]-data_properties.border[2*(n_dims - i -1)+1];
    }
  } else {   
    for(int i = 0; i < n_dims; i++){
      hdims[i]   = dims_T[n_dims - 1 - i];
      counts[i]  = hdims[i] / grid[n_dims - i - 1];
      offsets[i] = myrank_arr[n_dims - i - 1] * counts[i];
      if(myrank_arr[n_dims-i-1] < (hdims[i] % grid[n_dims-i-1])){
        counts[i]++;
        offsets[i] += myrank_arr[n_dims - i - 1];
      } else
	      offsets[i] += (hdims[i] % grid[n_dims - i - 1]);
    }
  }
  
  if(myrank == 0 || outfile) {
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);    
    hid_t global_memory_space = H5Screate_simple(n_dims, hdims, NULL);
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    dataset_id = H5Dopen(file_id, dataset_out, H5P_DEFAULT);		// Open dataset 
    if (dataset_id < 0) 
      dataset_id = H5Dcreate(file_id, dataset_out, type, global_memory_space,
			     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	// Create dataset 
    
    dataspace = H5Dget_space(dataset_id);			// Copy dataset 
  } else {    
    MPI_Recv(&message, 1, MPI_INT, myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
    dataset_id = H5Dopen(file_id,dataset_out, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset_id);   
  }
   
  if(FLOAT_TYPE == 1)
    type =  H5T_NATIVE_FLOAT;
  else if(sizeof(value) == 1)
    type =  negative ? H5T_NATIVE_CHAR: H5T_NATIVE_UCHAR;
  else if(sizeof(value) == 2)
    type =  negative ? H5T_NATIVE_SHORT: H5T_NATIVE_USHORT;
  else if(sizeof(value) == 4)
    type =  negative ? H5T_NATIVE_INT: H5T_NATIVE_UINT;
  
  ulong off_out  = data_properties.border[0] + data_properties.border[2]*dims[0] + data_properties.border[4]*dims[0]*dims[1];
  ulong off_save = offsets[1];
  ulong copy[3]  = {counts[0], counts[1], counts[2]};
  ulong s 	   = 0;

  for(ulong i = 0; i < copy[0]; i++){
    counts[0] = 1;
    if(n_dims == 2) {
      hid_t memory_window = H5Screate_simple((hsize_t)n_dims, counts, NULL);
      H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
      H5Dwrite(dataset_id, type, memory_window, dataspace, H5P_DEFAULT, out+i*dims[0]+off_out);
    } else {
      offsets[1]= off_save;
      for(ulong j = 0; j < copy[1]; j++, s++){
          counts[1] = 1;
          hid_t memory_window = H5Screate_simple((hsize_t)n_dims, counts, NULL);
          H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
          H5Dwrite(dataset_id, type, memory_window, dataspace, H5P_DEFAULT, out+off_out+(s+i*(data_properties.border[2]+data_properties.border[3]))*dims[0]);
          offsets[1]++;
      }
    }
    offsets[0]++;
  }
  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Fclose(file_id);
  if((myrank < np()-1) && !outfile)
    MPI_Send(&myrank, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);

  if(np()==1 || !outfile)
    info("File %s successfully written", fname);
  else
    info("Files (%s on process 0) successfully written", fname);

} /* write_hdf5 */

void write_fits(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3]){
  /* +++++++++++++++++++++++++++ */
  /*     FITS Write Function     */
  /* +++++++++++++++++++++++++++ */
  
  fitsfile 	*outfptr; 	/* FITS file pointers */
  int 		status = 0; 		/* CFITSIO status value MUST be initialized to zero! */
  int 		message;
  int 		type;
  int 		grid[3]       = {args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]};
  int 		myrank        = rank();	
  int 		myrank_2D     = myrank % (grid[1]*grid[0]);
  int 		myrank_arr[3] = {myrank_2D % grid[0], myrank_2D / grid[0], myrank / (grid[1]*grid[0])};
  long 		naxes[3]      = {1,1,1};      
  long 		counts[3]     = {1,1,1};
  long 		offsets[3]    = {1,1,1};
  int 		n_dims 	      = dims[2] > 1 ? 3 : 2;
  int 		negative      =  0;
  int 		outfile       = !strcmp(args->output_file_style, "single")? 0 : 1;
  int 		bitpix        = args->bpp;
  int           ttype         = sizeof(value) == 1 ? 1-negative: negative;
  char 		str1[100];
  strcpy(str1, "!");   		/* '!' symbol makes the output image, if existing, to be overwritten */
  if(outfile) {
    for (int i = 0; i < n_dims; i++) 
      naxes[i] = counts[i] = dims[i]-data_properties.border[2*i]-data_properties.border[2*i+1];     //-data_properties.border[2*i]-data_properties.border[2*i+1]
  } else {  
    for (int i = 0; i < n_dims; i++) {
      naxes[i]   = dims_T[i];
      counts[i]  = naxes[i] / grid[i];
      offsets[i] = myrank_arr[i]*counts[i] + 1;
      if (myrank_arr[i] < naxes[i]%grid[i]) {
        counts[i]++;
        offsets[i] += myrank_arr[i];
      } else {
	      offsets[i] += naxes[i]%grid[i];
      }
      counts[i] += offsets[i] - 1;
    }
  }
  if(FLOAT_TYPE || args->bpp == -32){
    type = TFLOAT;
    out = (float *) out;
  }else
    type = (int) (log2(sizeof(value))+1)*10 + ttype;

  if ((myrank == 0) || outfile) {
    if (!fits_create_file(&outfptr, strcat(str1, fname), &status)) {
      if (bitpix < 0)
	      fits_create_img(outfptr, FLOAT_IMG, n_dims, naxes, &status);
      else if (bitpix <= 8)
	      negative ? fits_create_img(outfptr, SBYTE_IMG, n_dims, naxes, &status):
	                 fits_create_img(outfptr, BYTE_IMG, n_dims, naxes, &status);
      else if (bitpix <= 16)
	      negative ? fits_create_img(outfptr, SHORT_IMG, n_dims, naxes, &status):
	                 fits_create_img(outfptr, USHORT_IMG, n_dims, naxes, &status);
      else if (bitpix <= 32)
      	negative ? fits_create_img(outfptr, LONG_IMG, n_dims, naxes, &status):
	                 fits_create_img(outfptr, ULONG_IMG, n_dims, naxes, &status);
      else 
	      fits_create_img(outfptr, LONGLONG_IMG, n_dims, naxes, &status);
      
      // fits_update_key(outfptr, TSTRING , "FILTER", args->operation_arg, "Type of filter used", &status);
      // fits_update_key(outfptr, TSTRING , "ATTRIBUTE", AttribsArray[args->attribute_arg].name, "Attribute name", &status);
      // fits_update_key(outfptr, TSTRING , "DECISION", Decisions[args->decision_arg].name, "Pruning rule", &status);
      // fits_update_key(outfptr, TINT ,    "CONNECTIVITY",  &args->pixel_connectivity, "Neighbour connectivity", &status);
  //     if(!strcmp(args->operation_arg, "filter"))
	// fits_update_key(outfptr, TULONG , "LAMBDA", &args->lambda_arg, "Threshold value", &status);
  //     else
	// fits_update_key(outfptr, TSTRING , "LVEC", args->lvec_arg, "Name of threshold file", &status);
  //     fits_update_key(outfptr, TSTRING , "COMMENT", "Filtering using connected filter method DISCCOMAN", NULL, &status);
     
    } else {
      error("Cannot create the file");
      MPI_Abort(MPI_COMM_WORLD, 730);
    }
  } else {
    MPI_Recv(&message, 1, MPI_INT, myrank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    fits_open_file(&outfptr, fname, READWRITE, &status); // open input images
  }

  //  if(!outfile){
  ulong off_out   = data_properties.border[0] + data_properties.border[2]*dims[0] + data_properties.border[4]*dims[0]*dims[1];
  ulong off_save  = offsets[1];
  ulong s	  = 0;
  counts[1] = off_save;
  counts[2] = offsets[2];

  for(ulong i = 0; i < dims[2]-data_properties.border[5]-data_properties.border[4]; i++){
    offsets[1] = off_save;
    counts[1] = off_save;
    for(ulong j = 0; j < dims[1]-data_properties.border[3]-data_properties.border[2]; j++,s++){
      fits_write_subset(outfptr, type, offsets, counts, out+off_out+(s+i*(data_properties.border[2]+data_properties.border[3]))*dims[0], &status);
      offsets[1]++;
      counts[1]++;
    }
    counts[2]++;
    offsets[2]++;
  }
  // }else{
  //   fits_write_subset(outfptr, type, offsets, counts, out, &status);
    // }
  fits_close_file(outfptr, &status);
  if((myrank < np()-1) && !outfile)
    MPI_Send(&status, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);

  if(np()==1 || !outfile)
    info("File %s successfully written", fname);
  else
    info("Files (%s on process 0) successfully written", fname);
} /* write_fits */

void write_basic(Arguments *args, const char* fname, value *out, ulong dims_T[3], ulong dims[3]){
  
  /* +++++++++++++++++++++++++++ */
  /*        Write FreeImage      */
  /* +++++++++++++++++++++++++++ */

  FIBITMAP 	*outmap;
  ulong 	width  	      = dims[0];
  ulong 	height 	      = dims[1];
  int 		negative      = 0;
  int 		outfile       = !strcmp(args->output_file_style, "single")? 0 : 1;
  int 		bitpix        = args->bpp;


  
  debug("Writing image %s using FreeImage", fname);

  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname);
  if (fif == FIF_UNKNOWN)
    fif = FreeImage_GetFIFFromFilename(fname);
  if ((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsWriting(fif))
  {
    if (bitpix <= 8)
    {
      ubyte *imagebuf = malloc(width * height * sizeof(ubyte));
      ;
      outmap = FreeImage_AllocateT(FIT_BITMAP, width - data_properties.border[0] - data_properties.border[1], height - data_properties.border[2] - data_properties.border[3], 8, 0xFF, 0xFF, 0xFF);
      for (ulong y = data_properties.border[2]; y < height - data_properties.border[3]; y++)
      {
        imagebuf = (BYTE *)FreeImage_GetScanLine(outmap, y - data_properties.border[2]);
        for (ulong x = data_properties.border[0]; x < width - data_properties.border[1]; x++)
          imagebuf[x - data_properties.border[0]] = out[width * y + x];
      }
    }
    else if (bitpix <= 16)
    {
      if (!negative)
      {
        ushort *imagebuf;
        outmap = FreeImage_AllocateT(FIT_UINT16, width - data_properties.border[0] - data_properties.border[1], height - data_properties.border[2] - data_properties.border[3], 16, 0xFFFF, 0xFFFF, 0xFFFF);
        for (ulong y = data_properties.border[2]; y < height - data_properties.border[3]; y++)
        {
          imagebuf = (ushort *)FreeImage_GetScanLine(outmap, y - data_properties.border[2]);
          for (ulong x = data_properties.border[0]; x < width - data_properties.border[1]; x++)
            imagebuf[x - data_properties.border[0]] = out[width * y + x];
        }
      }
      else
      {
        short *imagebuf;
        outmap = FreeImage_AllocateT(FIT_INT16, width - data_properties.border[0] - data_properties.border[1], height - data_properties.border[2] - data_properties.border[3], 16, 0xFFFF, 0xFFFF, 0xFFFF);
        for (ulong y = data_properties.border[2]; y < height - data_properties.border[3]; y++)
        {
          imagebuf = (short *)FreeImage_GetScanLine(outmap, y - data_properties.border[2]);
          for (ulong x = data_properties.border[0]; x < width - data_properties.border[1]; x++)
            imagebuf[x - data_properties.border[0]] = out[width * y + x];
        }
      }
    }
    else if (bitpix <= 32)
    {
      if (!negative)
      {
        uint *imagebuf;
        outmap = FreeImage_AllocateT(FIT_UINT32, width - data_properties.border[0] - data_properties.border[1], height - data_properties.border[2] - data_properties.border[3], 32, 0xFFFF, 0xFFFF, 0xFFFF);
        for (ulong y = data_properties.border[2]; y < height - data_properties.border[3]; y++)
        {
          imagebuf = (uint *)FreeImage_GetScanLine(outmap, y - data_properties.border[2]);
          for (ulong x = data_properties.border[0]; x < width - data_properties.border[1]; x++)
            imagebuf[x - data_properties.border[0]] = out[width * y + x];
        }
      }
      else
      {
        int *imagebuf;
        outmap = FreeImage_AllocateT(FIT_INT32, width - data_properties.border[0] - data_properties.border[1], height - data_properties.border[2] - data_properties.border[3], 32, 0xFFFF, 0xFFFF, 0xFFFF);
        for (ulong y = data_properties.border[2]; y < height - data_properties.border[3]; y++)
        {
          imagebuf = (int *)FreeImage_GetScanLine(outmap, y - data_properties.border[2]);
          for (ulong x = data_properties.border[0]; x < width - data_properties.border[1]; x++)
            imagebuf[x - data_properties.border[0]] = out[width * y + x];
        }
      }
    }
    FreeImage_Save(fif, outmap, fname, 0);
    FreeImage_Unload(outmap);
  }
  else
  {
    error("FreeImage couldn't write the output");
    MPI_Abort(MPI_COMM_WORLD, 725);
  }

  if(np()==1 || !outfile)
    info("File %s successfully written", fname);
  else
    info("Files (%s on process 0) successfully written", fname);
}
