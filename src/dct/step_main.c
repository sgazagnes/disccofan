/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*		       				 	 */
/*               Distributed Component Tree              */
/*              Including attribute filtering            */
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
#include "boundary_step.h"
#include "attributes.h"
#include "image.h"
#include "tree_flood.h"
#include "tree_filt.h"
#include "tree_step.h"
#include "lambdavec.h"
#include "refine.h"
#include "writefile.h"
//#include "moschini.h"

/* +++++++++++++++++++++++++++++++ */
/*     	   Global Variables        */
/* +++++++++++++++++++++++++++++++ */

float 		g_max_greyval;
ulong 		g_max_levels;
struct tms 	tstruct;				
clock_t 	start;					
clock_t 	inter;					

/* ... */


/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Main Function         */
/*				   */
/* +++++++++++++++++++++++++++++++ */


int main(int argc, char** argv) {
  
 
  /* +++++++++++++++++++++++++++ */
  /*        Parsing Values       */
  /* +++++++++++++++++++++++++++ */
  
  init_mpi();
  if(rank() == 0) printf("\n******* DISCCOMAN , starting... ***** \n\n");
  
  Arguments 	args;
  parse_args(argc, argv, &args);			/* Read input arguments */
  check_bytes();

  create_mpi_value_type();				/* Create MPI struct for value type */
  create_mpi_boundnode_type();                          /* Create MPI struct for boundary nodes */ 
  create_mpi_borderindex_type();                        /* Create MPI struct for borderindex struct */

  double	*copy_attr	= NULL; 	       	/* Attributes copy (pattern spectrum)  */ 
  value 	*gvals_par	= NULL;			/* Parents node intensities (pattern spec) */
  LambdaVec 	*lvec		; 
  ulong 	dims_tile[3]    = {1, 1, 1};		/* Tile dimensions */
  ulong 	dims_img[3]     = {1, 1, 1};		/* Image dimensions */
  ulong 	size_tile;				/* Tile size */
  Node          *local_tree 	= malloc(1* sizeof(Node));
  int attrib = -2;
  double *extrema;
  
  set_border(&args, local_tree);			/* Check overlapping borders */
  local_tree->gval   = read_input(&args,args.inprefix_arg, dims_tile, dims_img, local_tree->border);
  size_tile   	     = dims_tile[0]*dims_tile[1]*dims_tile[2];
  local_tree->size   = size_tile;

  refine(&args,  local_tree->gval, extrema, size_tile); /* Check is values need refinement */
  set_flooding(&args, args.bpp_arg);			/* Check flooding choice */
  set_connectivity(&args, dims_tile[2]);		/* Check connectivity choice */ 
  check_operation(&args, local_tree,  g_max_greyval); /* Check morph operation */
  ulong *attr_off = attribute_offsets(&args, dims_img);

  /* Initialization of attributes functions */

 
  if (rank() == 0){
    print_args(&args, dims_tile);
    if(args.bpp_arg < 0 && FLOAT_TYPE == 0) warn("Floating point not activated");
    if(args.bpp_arg > 0 && FLOAT_TYPE == 1) warn("Floating point activated, but data is not. Might lead to numerical errors");
  }
  MPI_Barrier(MPI_COMM_WORLD);
  start = times(&tstruct);

  /* +++++++++++++++++++++++++++ */
  /*     Local tree building     */
  /* +++++++++++++++++++++++++++ */
  local_tree->parent = malloc(local_tree->size * sizeof(idx)); 
  build_tree_parents(&args, local_tree, dims_tile);
  MPI_Barrier(MPI_COMM_WORLD);
  timing("Parent tree built: wallclock time = %0.2f",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

  
  if (np() > 1)
    local_tree = correct_borders_parents(&args, local_tree, dims_tile);

  timing("Parent tree correct: wallclock time = %0.2f",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  // ulong *ranks        = calloc(size_tile, sizeof(ulong));  check_alloc(ranks, 604);
  // ulong *ranks_inv    = calloc(size_tile, sizeof(ulong));  check_alloc(ranks_inv, 605);
  // create_mappings(local_tree->gval, ranks, ranks_inv, size_tile, 0, size_tile);

      
  ulong *ranks = calloc(local_tree->size, sizeof(ulong));
  create_mappings(local_tree->gval, NULL, ranks, local_tree->size, 0, local_tree->size);
  local_tree->store =  malloc(1 * sizeof(AuxDataStore*)); check_alloc(local_tree->store, 001);

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
    //  start = times(&tstruct);
    args.attribute_arg = attrib;

    new_aux_data          = AttribsArray[attrib].new_data;
    delete_aux_data       = AttribsArray[attrib].delete_data;
    add_to_aux_data       = AttribsArray[attrib].add_to_data;  
    merge_aux_data        = AttribsArray[attrib].merge_data;
    merge_to_aux_data     = AttribsArray[attrib].merge_to_data;
    clone_aux_data        = AttribsArray[attrib].clone_data;
    create_mpi_aux_data   = AttribsArray[attrib].create_mpi_data;
    //read_aux_file_binary  = AttribsArray[attrib].read_data_file;
    // write_aux_file_binary = AttribsArray[attrib].write_data_file;
    create_mpi_aux_data();
   
    local_tree->attribute = calloc(local_tree->size, sizeof(void*));
    local_tree->store[0] =  malloc(sizeof(AuxDataStore)); check_alloc(local_tree->store[0], 0);
    init_aux_data_store(local_tree->store[0], AttribsArray[attrib].size, local_tree->size);
    //  PrioQueue *q        = create_prio_queue(size_tile);	    
    bool *visited   	= calloc( local_tree->size, sizeof(bool));  check_alloc(visited, 603);
    //tree_flood_tee_att(local_tree, local_tree->store[0], visited, q, ranks, ranks_inv, dims_tile, attr_off, 0, size_tile, args.connectivity_arg);

    double *init_attr = calloc(4, sizeof(double));
    for( long i = local_tree->size-1; i >= 0; i--) {
      idx p = ranks[i];
      if(!local_tree->attribute[p] && !is_border(local_tree->border, dims_tile, p) && p < size_tile){
	init_attr[0] = (double) ((p % (dims_tile[0] * dims_tile[1])) % dims_tile[0]  + attr_off[0]);
	init_attr[1] = (double) ((p % (dims_tile[0] * dims_tile[1])) / dims_tile[0]  + attr_off[1]);
	init_attr[2] = (double) (p / (dims_tile[0] * dims_tile[1]) + attr_off[2]);
	init_attr[3] = (double) 1;
	local_tree->attribute[p] = new_aux_data(local_tree->store[0], init_attr);
      }
      visited[p] = true;
      idx parent = get_parent(local_tree, p);
      if(parent == BOTTOM) continue;
      if(!local_tree->attribute[parent] && !is_border(local_tree->border, dims_tile, parent) && parent < size_tile){
	init_attr[0] = (double) ((parent % (dims_tile[0] * dims_tile[1])) % dims_tile[0]  + attr_off[0]);
	init_attr[1] = (double) ((parent % (dims_tile[0] * dims_tile[1])) / dims_tile[0]  + attr_off[1]);
	init_attr[2] = (double) (parent / (dims_tile[0] * dims_tile[1]) + attr_off[2]);
	init_attr[3] = (double) 1;
	local_tree->attribute[parent] = new_aux_data(local_tree->store[0], init_attr);

      }
      if(local_tree->attribute[p]){
	if((!is_border(local_tree->border, dims_tile, parent) && parent < size_tile) || local_tree->attribute[parent])
	  merge_aux_data(local_tree->attribute[parent], local_tree->attribute[p]);
	else{
	  clone_aux_data(local_tree->store[0], &local_tree->attribute[parent], local_tree->attribute[p]);
	}
	if(local_tree->gval[p] == local_tree->gval[parent] && visited[parent]){
	  parent = get_parent(local_tree, parent);
	  if(parent != BOTTOM)
	    local_tree->attribute[parent] ? merge_aux_data(local_tree->attribute[parent], local_tree->attribute[p]):  clone_aux_data(local_tree->store[0], &local_tree->attribute[parent], local_tree->attribute[p]);
	}
      }

    }

    free(visited);
   // free(q);
    timing("Local tree flood: total wallclock time = %0.2f",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  
  
  /* +++++++++++++++++++++++++++ */
  /*     Update tree borders     */
  /* +++++++++++++++++++++++++++ */

    if (np() > 1){
      if (!strcmp(args.filter_arg, "pattern")) {
	//Need to copy attributes and parent values for distributed pattern spectra //
	gvals_par = calloc(size_tile, sizeof(value));   check_alloc(gvals_par, 1);
	copy_attr = calloc(size_tile, sizeof(double));  check_alloc(copy_attr, 2);
           
	#pragma omp parallel for
	for (ulong i = 0; i < size_tile; i++) {
	  copy_attr[i] = is_levelroot(local_tree, i) && local_tree->attribute[i] ?
	    (*AttribsArray[attrib].attribute)(local_tree->attribute[i]) : -DBL_MAX;	    
	  if(local_tree->parent[i] != BOTTOM)
	    gvals_par[i] = local_tree->gval[local_tree->parent[i]];
	  else if(FLOAT_TYPE == 1)
	    gvals_par[i] = -FLT_MAX;
	  else
	    gvals_par[i] = INT_MIN;
	}     
      }
    
      correct_borders_att(&args, local_tree, dims_tile, local_tree->store[0]);
      MPI_Barrier(MPI_COMM_WORLD);
      timing("Tiles attributes corrected: wallclock time = %0.2f",
	     (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    }


    /* +++++++++++++++++++++++++++++++ */
    /*         Process the tree        */
    /* +++++++++++++++++++++++++++++++ */
  

    //        CASE 1: Filtering        //
    char cmdd[1000];
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

      if (!strcmp(args.filter_arg, "filter")) {
	value *out_filter = calloc(size_tile, sizeof(value));   check_alloc(out_filter, 11);
	tree_filtering(local_tree, out_filter, size_tile, args.decision_arg, attrib, args.lambda_arg);
	if(rank() == 0) timing("Tree filtering: wallclock time = %0.2f",
			       (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	if(args.saveout_arg){
	  if ((!strcmp(args.tree_arg, "min") && !strcmp(args.morphology_arg, "opening")) ||
	      (!strcmp(args.tree_arg, "max") && !strcmp(args.morphology_arg, "closing"))) {
	    #pragma omp parallel for
	    for (ulong i = 0; i < size_tile; i++) 
	      out_filter[i] = g_max_greyval - out_filter[i];
	  }

	  write_output(&args, out_filter, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border, args.bpp_arg);
	  if(rank() == 0){
	    sprintf(&cmdd, "display %s.%s", args.outprefix_arg, args.outtype_arg);
	    system(cmdd);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	  timing("File written, wallclock time = %0.2f",
			       (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	}
	free(out_filter);
      }

      //   CASE 2: Differential profile  //

      else if (!strcmp(args.filter_arg, "csl")) {
      
	value *out_orig   = calloc(local_tree->size, sizeof(value));   check_alloc(out_orig,   3);
	value *out_dh     = calloc(local_tree->size, sizeof(value));   check_alloc(out_dh,     4);
	value *out_scale  = calloc(local_tree->size, sizeof(value));   check_alloc(out_scale,  5); 
	value *temp_scale = calloc(local_tree->size, sizeof(value));   check_alloc(temp_scale, 6);
	value *temp_dh    = calloc(local_tree->size, sizeof(value));   check_alloc(temp_dh,    7);
	bool *temp_valid  = calloc(local_tree->size, sizeof(bool));    check_alloc(temp_valid, 8);

	tree_differential(local_tree, size_tile, lvec, out_dh, temp_dh, out_orig, out_scale, temp_scale, temp_valid,  AttribsArray[attrib].attribute);
	//    combine_results(local_tree, size_tile, lvec, out_dh, out_dh2, out_orig, out_orig2, out_scale, out_scale2);
      
	MPI_Barrier(MPI_COMM_WORLD);
	timing("CSL segmentation: wallclock time = %0.2f",
				(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	if(args.saveout_arg)
	  write_differential(&args, out_orig, out_dh, out_scale, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border, args.bpp_arg);

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

      else if (!strcmp(args.filter_arg, "pattern")) {
	double *all_spectrum  = calloc(lvec->num_lambdas, sizeof(double)); check_alloc(all_spectrum, 18);
	double *loc_spectrum  = calloc(lvec->num_lambdas, sizeof(double)); check_alloc(loc_spectrum, 19);
      
	tree_pattern_spectrum(local_tree, size_tile, lvec, copy_attr, gvals_par, loc_spectrum, args.background_arg, AttribsArray[attrib].attribute);

	MPI_Reduce(loc_spectrum, all_spectrum, lvec->num_lambdas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank() == 0) timing("Pattern spectra built: wallclock time = %0.2f",
				(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
      
	if(args.saveout_arg && rank() == 0)
	  write_pattern_spectra(&args, all_spectrum, lvec->num_lambdas);

	timing("Pattern spectra written, wallclock time = %0.2f",
				(float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

	free(gvals_par);
	free(copy_attr);
	free(loc_spectrum);
	free(all_spectrum);
      } else {
	info("No filter chosen");
	MPI_Barrier(MPI_COMM_WORLD);

	if(args.saveout_arg)
	  write_output(&args, local_tree->gval, "None", dims_img, dims_tile, local_tree->border, args.bpp_arg);
	
	MPI_Barrier(MPI_COMM_WORLD);
	timing("File written, wallclock time = %0.2f",
	       (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
      }
      lambda_vector_delete(lvec);       
    } while (args.lambda_arg != -1);
    free(local_tree->attribute);
    MPI_Type_free(&mpi_attribute_type);
    clear_aux_data_store(local_tree->store[0]);
  }while(attrib != -1);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank() == 0) printf("\n******* End of DISCCOMAN interactive, cleaning ... ***** \n\n");

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
