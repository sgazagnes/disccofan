/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*		       				 	 */
/*               Distributed Component Tree              */
/* Including component analtysis and attribute filtering */
/*                 Author: Simon Gazagnes                */
/*                University of Groningen                */
/*                    Compiler : C99                     */
/*                   Libraries needed:                   */
/*      OpenMPI, OpenMP, FreeImage, CFITSIO, HDF5.       */
/*		       				 	 */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "types.h"
#include "arguments.h"
#include "boundary.h"
#include "attributes.h"
#include "image.h"
#include "flood.h"
#include "filter.h"
#include "lambdavec.h"
#include "refine.h"
#include "workspace.h"

/* +++++++++++++++++++++++++++++++ */
/*     	   Global Variables        */
/* +++++++++++++++++++++++++++++++ */

float 		g_max_greyval;
ulong 		g_max_levels;
struct tms 	tstruct;				
clock_t 	start;					
f_new_aux_data  new_aux_data = NULL;
f_init_aux_data init_aux_data= NULL;
f_delete_aux_data delete_aux_data= NULL;
f_add_to_aux_data  add_to_aux_data= NULL;  
f_merge_aux_data  merge_aux_data= NULL;
f_merge_to_aux_data merge_to_aux_data= NULL;
f_clone_aux_data clone_aux_data= NULL;
f_create_mpi_aux_data create_mpi_aux_data= NULL;

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Main Function         */
/*				   */
/* +++++++++++++++++++++++++++++++ */


int main(int argc, char** argv) {
  

  /* +++++++++++++++++++++++++++ */
  /*         Initializing        */
  /* +++++++++++++++++++++++++++ */

  init_mpi();
  
  Arguments 	args;
  ulong 	dims_tile[3]    = {1, 1, 1};	
  ulong 	dims_img[3]     = {1, 1, 1}; 
  parse_args(argc, argv, &args);			// Read input arguments //

  info("******* DISCCOMAN, starting... ******* \n");

  check_bytes();
  create_mpi_value_type();				// Create MPI struct for value type //
  create_mpi_boundnode_type();                          // Create MPI struct for boundary nodes // 
  create_mpi_borderindex_type();                        // Create MPI struct for borderindex struct //

  Node          *local_tree 	=  malloc(1* sizeof(Node));
	
  read_input(&args, &(local_tree->gval), dims_tile, dims_img, local_tree->offsets, local_tree->border);
  local_tree->size_curr = dims_tile[0]*dims_tile[1]*dims_tile[2];
  check_refine(&args,  local_tree->gval, local_tree->size_curr);   
  check_operation(&args, local_tree,  g_max_greyval);    

  // Initialization of attributes functions //
  int attrib;

  if (rank() == 0){
    print_args(&args, dims_tile);
    if(args.bpp_arg < 0 && FLOAT_TYPE == 0) warn("Floating point not activated");
    if(args.bpp_arg > 0 && FLOAT_TYPE == 1) warn("Floating point activated, but data is %d bits.", args.bpp_arg);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  info("Starting clock");
  start = times(&tstruct);
 

  /* ++++++++++++++++++++++++++++++++++ */
  /*     Local tree parent building     */
  /* ++++++++++++++++++++++++++++++++++ */
 
  build_local_tree_par(&args, local_tree, dims_tile);
  MPI_Barrier(MPI_COMM_WORLD);

  timing("Parents computed: wallclock time = %0.2f",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

  
  /* ++++++++++++++++++++++++++++++++++ */
  /*        Update parents border       */
  /* ++++++++++++++++++++++++++++++++++ */
  local_tree->size_init 	= local_tree->size_curr;

  if (np() > 1){
    local_tree = correct_borders_parents(&args, local_tree, dims_tile);
    timing("Parent tree correct: wallclock time = %0.2f",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  }
  
  /* +++++++++++++++++++++++++++ */
  /*     Flood attributes in     */
  /* +++++++++++++++++++++++++++ */

  LambdaVec* lvec       = NULL;
  ulong *ranks_inv  	= NULL;
  create_mappings(local_tree->gval, &ranks_inv, 0, local_tree->size_curr);

  do{
     if(rank() == 0){
       do{
	 attrib = -2;
	 info("Choose one attribute in the following list, or write -1 to exit:");
	 for (int i=0; i<NUMATTR; i++)
	   info("\t%d - %s", i, AttribsArray[i].name);
	 char term;
	 int err = scanf(" %d%c", &attrib, &term);
	 if(err != 2 || term != '\n')
	   {
	     while (fgetc(stdin) != '\n');
	     continue;
	   }
       }while((attrib < -1 ||attrib >= NUMATTR));
       if(attrib > -1)
	 info("Attribute chosen: %s", AttribsArray[attrib].name);
     }
     MPI_Bcast(&attrib, 1, MPI_INT, 0, MPI_COMM_WORLD);
     if (attrib == -1) break;
     MPI_Barrier(MPI_COMM_WORLD);
    args.attribute_arg = attrib;

    local_tree->size_attr = AttribsArray[attrib].size;
    new_aux_data          = AttribsArray[attrib].new_data;
    init_aux_data         = AttribsArray[attrib].init_data;
    delete_aux_data       = AttribsArray[attrib].delete_data;
    add_to_aux_data       = AttribsArray[attrib].add_to_data;  
    merge_aux_data        = AttribsArray[attrib].merge_data;
    merge_to_aux_data     = AttribsArray[attrib].merge_to_data;
    clone_aux_data        = AttribsArray[attrib].clone_data;
    create_mpi_aux_data   = AttribsArray[attrib].create_mpi_data;
    create_mpi_aux_data();


    build_local_tree_att(local_tree, ranks_inv, dims_tile);
    timing("Local tree flood: total wallclock time = %0.2f",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    
  //
  double	*copy_attr	= NULL; 	       	// Attributes copy (pattern spectrum)  // 
  value 	*gvals_par	= NULL;			// Parents node intensities (pattern spec) //
  
  if (np() > 1){
    if (!strcmp(args.output_arg, "pattern")) {
     
      //Need to copy attributes and parent values for distributed pattern spectra //
      gvals_par = calloc(local_tree->size_curr, sizeof(value));   check_alloc(gvals_par, 1);
      copy_attr = calloc(local_tree->size_curr, sizeof(double));  check_alloc(copy_attr, 2);
      #pragma omp parallel for
      for (ulong i = 0; i < local_tree->size_curr; i++) {
	copy_attr[i] = is_levelroot(local_tree, i) ?
	  (*AttribsArray[args.attribute_arg].area)(local_tree->attribute + i*local_tree->size_attr) : -DBL_MAX;	

	if(local_tree->parent[i] != BOTTOM)
	  gvals_par[i] = local_tree->gval[local_tree->parent[i]];
      }    
    }
    
    correct_borders_att(&args, local_tree, dims_tile);
    timing("Tiles attributes corrected: wallclock time = %0.2f",
	     (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  }

  /* +++++++++++++++++++++++++++++++ */
  /*         Process the tree        */
  /* +++++++++++++++++++++++++++++++ */
  

  //        CASE 1: Filtering        //
  char *cmdd;
  do{       
    info("Enter the new lambda value or modify the vector file. Write -1 to change the attribute.\n New image will replace the former one.");
    double new_val;
    if(rank() == 0){
      do{
	char term2;
	info("New value: ");
	int err2 = scanf(" %lf%c", &new_val, &term2);
	if(err2 != 2 || term2 != '\n') {
	  while (fgetc(stdin) != '\n');
	  continue;
	}
      }while(new_val < -1);
    }
    MPI_Bcast(&new_val, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(new_val == -1) break;
    args.lambda_arg = new_val;
    MPI_Barrier(MPI_COMM_WORLD);

    start = times(&tstruct);
    lvec = lambda_vector_read(argv[0], args.lvec_arg, args.imscale_arg);

    if (!strcmp(args.output_arg, "filter")) {
      value *out_filter = calloc(local_tree->size_init, sizeof(value));   check_alloc(out_filter, 11);
      tree_filtering(local_tree, out_filter, local_tree->size_init, args.decision_arg, attrib, args.lambda_arg);

      if(args.saveout_arg){
	if ((!strcmp(args.tree_arg, "min") && !strcmp(args.morphology_arg, "opening")) ||
	    (!strcmp(args.tree_arg, "max") && !strcmp(args.morphology_arg, "closing"))) {
	  #pragma omp parallel for
	  for (ulong i = 0; i < local_tree->size_init; i++) 
	      out_filter[i] = g_max_greyval - out_filter[i];
	  }

	  write_output(&args, out_filter, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border);
	  if(rank() == 0){
	    asprintf(&cmdd, "display %s.%s", args.outprefix_arg, args.outtype_arg);
	    system(cmdd);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	  timing("File written, wallclock time = %0.2f",
			       (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	}
	free(out_filter);
      }

      //   CASE 2: Differential profile  //

      else if (!strcmp(args.output_arg, "csl")) {
      
	value *out_orig   = calloc(local_tree->size_curr, sizeof(value));   check_alloc(out_orig,   3);
	value *out_dh     = calloc(local_tree->size_curr, sizeof(value));   check_alloc(out_dh,     4);
	value *out_scale  = calloc(local_tree->size_curr, sizeof(value));   check_alloc(out_scale,  5);
	value *temp_scale = calloc(local_tree->size_curr, sizeof(value));   check_alloc(temp_scale, 6);
	value *temp_dh    = calloc(local_tree->size_curr, sizeof(value));   check_alloc(temp_dh,    7);
	bool *temp_valid  = calloc(local_tree->size_curr, sizeof(bool));    check_alloc(temp_valid, 8);

	tree_differential(local_tree, local_tree->size_init, lvec, out_dh, temp_dh, out_orig, out_scale, temp_scale, temp_valid,  AttribsArray[attrib].attribute);
      
	MPI_Barrier(MPI_COMM_WORLD);
	timing("CSL segmentation: wallclock time = %0.2f",
				(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	if(args.saveout_arg)
	  write_differential(&args, out_orig, out_dh, out_scale, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border);

	MPI_Barrier(MPI_COMM_WORLD);
	timing("Files written, wallclock time = %0.2f",
				(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	free(out_orig);
	free(out_dh);
	free(out_scale);
	free(temp_valid);
	free(temp_scale);
	free(temp_dh); 
      
      }
    
      //   CASE 3: Pattern spectra  //

      else if (!strcmp(args.output_arg, "pattern")) {
	double *all_spectrum  = calloc(lvec->num_lambdas, sizeof(double)); check_alloc(all_spectrum, 18);
	double *loc_spectrum  = calloc(lvec->num_lambdas, sizeof(double)); check_alloc(loc_spectrum, 19);
      
	tree_pattern_spectrum(local_tree, local_tree->size_init, lvec, copy_attr, gvals_par, loc_spectrum, args.background_arg, AttribsArray[attrib].area, AttribsArray[attrib].attribute);

	MPI_Reduce(loc_spectrum, all_spectrum, lvec->num_lambdas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank() == 0) timing("Pattern spectra built: wallclock time = %0.2f",
				(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
      
	if(args.saveout_arg && rank() == 0)
	  write_pattern_spectra(&args, all_spectrum, lvec);

	timing("Pattern spectra written, wallclock time = %0.2f",
				(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

	free(gvals_par);
	free(copy_attr);
	free(loc_spectrum);
	free(all_spectrum);
      } else 
	info("No filter chosen");

      
      lambda_vector_delete(lvec);       
    } while (args.lambda_arg != -1);
    free(local_tree->attribute);
    MPI_Type_free(&mpi_attribute_type);
    // clear_aux_data_store(local_tree->store[0]);
  }while(attrib != -1);
  
  MPI_Barrier(MPI_COMM_WORLD);
  info("\n******* End of DISCCOMAN interactive, cleaning ... ***** \n\n");

  
  /* +++++++++++++++++++++++++++++++ */
  /*              Clean Up           */
  /* +++++++++++++++++++++++++++++++ */
  free(local_tree->parent);
  free(local_tree->gval);
  free(local_tree);
  cmdline_parser_free(&args);
  MPI_Type_free( &mpi_bound_node_type );
  MPI_Type_free( &mpi_borderindex_type );
  finalize_mpi();
    
  return 0;
} /* main */
