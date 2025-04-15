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
#include "tree_flood.h"
#include "tree_filt.h"
#include "lambdavec.h"
#include "refine.h"
#include "writefile.h"
#include "eor_analysis.h"

//#include "moschini.h"

/* +++++++++++++++++++++++++++++++ */
/*     	   Global Variables        */
/* +++++++++++++++++++++++++++++++ */

float 		g_max_greyval;
ulong 		g_max_levels;
struct tms 	tstruct;				
clock_t 	start;					
ulong DIMS[3];

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Main Function         */
/*				   */
/* +++++++++++++++++++++++++++++++ */


int main(int argc, char** argv) {
  
 
  /* +++++++++++++++++++++++++++ */
  /*        Parsing Values       */
  /* +++++++++++++++++++++++++++ */
 printf("\n******* DISCCOMAN, starting... ***** \n\n");

 init_mpi();

  Arguments 	args;
  parse_args(argc, argv, &args);			/* Read input arguments */
  check_bytes();

  create_mpi_value_type();				/* Create MPI struct for value type */
  create_mpi_boundnode_type();                          /* Create MPI struct for boundary nodes */ 
  create_mpi_borderindex_type();                        /* Create MPI struct for borderindex struct */

  double	*copy_attr	= NULL; 	       	/* Attributes copy (pattern spectrum)  */ 
  value 	*gvals_par	= NULL;			/* Parents node intensities (pattern spec) */
  LambdaVec 	*lvec		= lambda_vector_read(argv[0], args.lvec_arg, args.imscale_arg); 
  ulong 	dims_tile[3]    = {1, 1, 1};		/* Tile dimensions */
  ulong 	dims_img[3]     = {1, 1, 1};		/* Image dimensions */
  ulong 	size_tile;				/* Tile size */
  int 		attrib 		= args.attribute_arg;
  Node          *local_tree 	= calloc(1, sizeof(Node));   check_alloc(local_tree, 000);
  double 	extrema[2] = {0};
  double 	*factors;

  set_border(&args, local_tree);	/* Check overlapping borders */
  local_tree->gval   = read_input(&args, args.inprefix_arg, dims_tile, dims_img, local_tree->border);
  size_tile   	     = dims_tile[0]*dims_tile[1]*dims_tile[2];
  DIMS[0]            = dims_tile[0];
  DIMS[1]            = dims_tile[1];
  DIMS[2]            = dims_tile[2];
  local_tree->size   = size_tile;
  ulong *attr_off    = attribute_offsets(&args, dims_img);
  set_flooding(&args, args.bpp_arg);			/* Check flooding choice */

  /*   +++   Pre-processing   +++  */
  if(!strcmp(args.eor_arg, "hi")){
    #pragma omp parallel for
    for (int i = 0; i < size_tile; ++i)
      {
	local_tree->gval[i] = 1-local_tree->gval[i];
      }
  }
  factors = refine(&args, local_tree->gval, extrema, size_tile); /* Check is values need refinement */


  /*
  for (int i = 0; i < size_tile; ++i)
  {
    if (local_tree->gval[i] < 0.0)
      local_tree->gval[i] = 0;
      }*/
  check_operation(&args, local_tree,  g_max_greyval);   /* Check morph operation */



  //args.flood_arg = 1;
  set_connectivity(&args, dims_tile[2]);		/* Check connectivity choice */ 
  
  if(args.density_arg != NULL){
    info("Density file to read: %s", args.density_arg);
    local_tree->gval_dens = read_input(&args, args.density_arg, dims_tile, dims_img, local_tree->border);

    local_tree->dens = 1;
    double tot = 0;
    for(ulong i = 0; i < size_tile; i++)
      tot += (double) local_tree->gval_dens[i];

    double mean =  tot/(double)size_tile;
    info("Density field mean is %lf", mean);
    #pragma omp parallel for
    for(ulong i = 0; i < size_tile; i++)
      local_tree->gval_dens[i] -= (value) mean;

    tot = 0;
    for(ulong i = 0; i < size_tile; i++)
      tot += (double) local_tree->gval_dens[i];

    mean =  tot/(double)size_tile;
    info("New density field mean is %lf", mean);
 
  }



  /* Initialization of attributes functions */
  new_aux_data          = AttribsArray[attrib].new_data;
  delete_aux_data       = AttribsArray[attrib].delete_data;
  add_to_aux_data       = AttribsArray[attrib].add_to_data;  
  merge_aux_data        = AttribsArray[attrib].merge_data;
  merge_to_aux_data     = AttribsArray[attrib].merge_to_data;
  clone_aux_data        = AttribsArray[attrib].clone_data;
  create_mpi_aux_data   = AttribsArray[attrib].create_mpi_data;
  //read_aux_file_binary  = AttribsArray[attrib].read_data_file;
  //write_aux_file_binary = AttribsArray[attrib].write_data_file;
  
  create_mpi_aux_data();

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

  local_tree->store =  malloc(args.threads_arg * sizeof(AuxDataStore*)); check_alloc(local_tree->store, 001);
  local_tree->attribute = calloc(local_tree->size, sizeof(void*));check_alloc(local_tree->attribute, 001);

  local_tree->parent = build_local_tree(&args, local_tree, dims_tile, attr_off);

  MPI_Barrier(MPI_COMM_WORLD);
  timing("Local tree built: wallclock time = %0.2f",
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
    
    local_tree = correct_borders(&args, local_tree, dims_tile);
 
    timing("Tiles border corrected: wallclock time = %0.2f",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));


  }

  /* +++++++++++++++++++++++++++++++ */
  /*         Process the tree        */
  /* +++++++++++++++++++++++++++++++ */
  

  //        CASE 1: Filtering        //
  
  if (!strcmp(args.filter_arg, "filter")) {
    /*value *out_filter = calloc(size_tile, sizeof(value));   check_alloc(out_filter, 3);
    tree_filtering(local_tree, out_filter, size_tile, args.decision_arg, attrib, args.lambda_arg);
    timing("Tree filtering: wallclock time = %0.2f",
			   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    if(args.saveout_arg){
       if ((!strcmp(args.tree_arg, "min") && !strcmp(args.morphology_arg, "opening")) ||
	  (!strcmp(args.tree_arg, "max") && !strcmp(args.morphology_arg, "closing"))) {
	#pragma omp parallel for
	for (ulong i = 0; i < size_tile; i++) 
	  out_filter[i] = g_max_greyval - out_filter[i];
       }
      write_output(&args, out_filter, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border, args.bpp_arg);
      
      timing("File written, wallclock time = %0.2f",
	     (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    }
    free(out_filter);*/
    float min_distance = 0.0;
    float alpha = 0.000001;
    mt_object_data mt_o;
    mt_objects_init(local_tree, &mt_o);
    mt_use_node_test_4(&mt_o, alpha, min_distance);
    mt_objects(&mt_o);
    value *out_filter = calloc(size_tile, sizeof(value));   check_alloc(out_filter, 3);
    for (ulong i = 0; i < local_tree->size; i++){
      if(mt_o.flags[i] & 8)
	out_filter[i] = 1;
    }
     write_output(&args, out_filter, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border, args.bpp_arg);
      
  }

  //   CASE 2: Differential profile  //

  else if (!strcmp(args.filter_arg, "csl")) {
      
    value *out_orig   = calloc(local_tree->size, sizeof(value));   check_alloc(out_orig,   4);
    value *out_dh     = calloc(local_tree->size, sizeof(value));   check_alloc(out_dh,     5);
    value *out_scale  = calloc(local_tree->size, sizeof(value));   check_alloc(out_scale,  6); 
    value *temp_scale = calloc(local_tree->size, sizeof(value));   check_alloc(temp_scale, 7);
    value *temp_dh    = calloc(local_tree->size, sizeof(value));   check_alloc(temp_dh,    8);
    bool *temp_valid  = calloc(local_tree->size, sizeof(bool));    check_alloc(temp_valid, 9);

    tree_differential(local_tree, size_tile, lvec, out_dh, temp_dh, out_orig, out_scale, temp_scale, temp_valid,  AttribsArray[attrib].attribute);
    //    combine_results(local_tree, size_tile, lvec, out_dh, out_dh2, out_orig, out_orig2, out_scale, out_scale2);
      
    MPI_Barrier(MPI_COMM_WORLD);
    timing("CSL segmentation: wallclock time = %0.2f",
			    (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    if(args.saveout_arg)
      write_differential(&args, out_orig, out_dh, out_scale, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border, args.bpp_arg);

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
    double *all_spectrum  = calloc(lvec->num_lambdas, sizeof(double)); check_alloc(all_spectrum, 10);
    double *loc_spectrum  = calloc(lvec->num_lambdas, sizeof(double)); check_alloc(loc_spectrum, 11);
      
    tree_pattern_spectrum(local_tree, size_tile, lvec, copy_attr, gvals_par, loc_spectrum, args.background_arg, AttribsArray[attrib].attribute);

    MPI_Reduce(loc_spectrum, all_spectrum, lvec->num_lambdas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    timing("Pattern spectra built: wallclock time = %0.2f",
			   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
      
    if(args.saveout_arg && rank() == 0)
      write_pattern_spectra(&args, all_spectrum, lvec->num_lambdas);

    timing("Pattern spectra written, wallclock time = %0.2f",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

    free(gvals_par);
    free(copy_attr);
    free(loc_spectrum);
    free(all_spectrum);
  }

  //   CASE 4: EoR/density  //

  else if (!strcmp(args.filter_arg, "eor")) {
    value *out_filter = calloc(size_tile, sizeof(value));   check_alloc(out_filter, 3);
    if (attrib == 4 && args.periodic_arg){
      info("Merging boundaries to account for periodic edges");
      fuse_periodic_3D(local_tree, local_tree->store, dims_tile, 26);
    } else if (args.periodic_arg) {
      info("Periodic boundaries is only accounted for EoR application");
    }
    if ((!strcmp(args.tree_arg, "min") && !strcmp(args.morphology_arg, "opening")) ||
	(!strcmp(args.tree_arg, "max") && !strcmp(args.morphology_arg, "closing"))) {
      #pragma omp parallel for
      for (ulong i = 0; i < local_tree->size; i++) 
	local_tree->gval[i] = g_max_greyval - local_tree->gval[i];
    }
    #pragma omp parallel for
    for (ulong i = 0; i <local_tree->size; i++){
      local_tree->gval[i] =  local_tree->gval[i]/(value)factors[1] + (value)factors[0];
    }
    extrema[0] = extrema[0]/factors[1] + factors[0];
    extrema[1] = extrema[1]/factors[1] + factors[0];

    info("New min value found: %3.3lf", extrema[0]);
    info("New max value found: %3.3lf", extrema[1]);
    
    //tree_filtering(local_tree, out_filter, size_tile, args.decision_arg, attrib, args.lambda_arg);
    //timing("Tree filtering: wallclock time = %0.2f",
    //	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    //  write_output(&args, out_filter, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border, args.bpp_arg);
    timing("Image written: wallclock time = %0.2f",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    free(out_filter);

    if(args.eor_arg == NULL || !strcmp(args.eor_arg, "hii") || !strcmp(args.eor_arg, "hi")){
      info("Assuming segmented data");
      inertia_attributes_hii(&args,local_tree, AttribsArray[attrib].attribute_arr,  dims_tile);
    }else if(!strcmp(args.eor_arg, "21all")){
      info("Getting every attributes on the 21cm field");
      write_inertiaall_bin(&args,local_tree, args.outprefix_arg, AttribsArray[attrib].attribute_arr, args.lambda_arg, size_tile, dims_img, attr_off, attrib, extrema);
    } else if(!strcmp(args.eor_arg, "21ine")){
      info("Getting the attributes from Inertia matrice");
      inertia_attributes_bins(&args,local_tree, AttribsArray[attrib].attribute_arr,  dims_tile ,attr_off, extrema);
    } else 
      info("Wrong choice ?");
    timing("Inertia file written, wallclock time = %0.2f",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  
  }else if (!strcmp(args.filter_arg, "segment")) {   
     value *out_filter = calloc(size_tile, sizeof(value));   check_alloc(out_filter, 3);
     bool *reached = calloc(local_tree->size, sizeof(value));
    ulong *ranks        = calloc(local_tree->size, sizeof(ulong));  check_alloc(ranks, 504);
    create_mappings(local_tree->gval, NULL, ranks, local_tree->size, 0, local_tree->size);
    tree_seg_dir(local_tree, out_filter, reached, ranks, 0, local_tree->size, 0.0, AttribsArray[attrib].attribute, args.lambda_arg);
    timing("Tree filtering: wallclock time = %0.2f",
			   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    if(args.saveout_arg){
      write_output(&args, out_filter, AttribsArray[attrib].name, dims_img, dims_tile, local_tree->border, args.bpp_arg);
      
      timing("File written, wallclock time = %0.2f",
	     (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    }
    free(out_filter);        
  }


  //   CASE 5: No Analysis  //

  
  else {
    info("No filter chosen");

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank() == 0) printf("Tree building: wallclock time = %0.2f \n",
			   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    if(args.saveout_arg)
      write_output(&args, local_tree->gval, "None", dims_img, dims_tile, local_tree->border, args.bpp_arg);
    timing("File written, wallclock time = %0.2f",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
  }

  if(rank() == 0) printf("\n******* End of DISCCOMAN, cleaning ... ***** \n\n");
  /* +++++++++++++++++++++++++++++++ */
  /*              Clean Up           */
  /* +++++++++++++++++++++++++++++++ */
  int nthreads = args.threads_arg;
 if(args.bpp_arg < 0 || args.bpp_arg >=16)
     nthreads = 1;
  for(int i =0; i < nthreads; i++)
    clear_aux_data_store(local_tree->store[i]);
  
  free_tree(local_tree, size_tile);
  lambda_vector_delete(lvec);
  cmdline_parser_free(&args);
  MPI_Type_free( &mpi_bound_node_type );
  MPI_Type_free( &mpi_borderindex_type );
  MPI_Type_free( &mpi_attribute_type );
  finalize_mpi();
    
  return 0;
} /* main */
