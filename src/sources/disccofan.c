/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*		       				 	                                   */
/*               Distributed Component Forest            */
/* Including component analtysis and attribute filtering */
/*                 Author: Simon Gazagnes                */
/*                    Compiler : C99                     */
/*                   Libraries needed:                   */
/*      OpenMPI, OpenMP, FreeImage, CFITSIO, HDF5.       */
/*		       				 	                                   */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "types.h"
#include "boundary.h"
#include "attributes.h"
#include "io.h"
#include "flood.h"
#include "lambdavec.h"
#include "preprocessing.h"
#include "workspace.h"
#include "config.h"
#include "operation.h"

/* +++++++++++++++++++++++++++++++ */
/*     	   Global Variables        */
/* +++++++++++++++++++++++++++++++ */

ulong 		      g_max_levels;
struct tms 	      tstruct;				
clock_t 	      start;

DataParams data_properties = {
  .dims_process = {1, 1, 1},
  .dims_full = {1, 1, 1},
  .pixdim = {1.0, 1.0, 1.0},
  .border = {false, false, false, false, false, false},
  .offsets = {0, 0, 0},
  .g_max_gval =0,
  .g_min_gval = 0
};

f_new_aux_data        new_aux_data            = NULL;
f_init_aux_data       init_aux_data           = NULL;
f_delete_aux_data     delete_aux_data         = NULL;
f_add_to_aux_data     add_to_aux_data         = NULL;  
f_merge_aux_data      merge_aux_data          = NULL;
f_merge_to_aux_data   merge_to_aux_data       = NULL;
f_clone_aux_data      clone_aux_data          = NULL;
f_create_mpi_aux_data create_mpi_aux_data     = NULL;

/* +++++++++++++++++++++++++++++++ */
/*     	   Helper functions        */
/* +++++++++++++++++++++++++++++++ */

Node *init_local_tree(){
  Node *local_tree 	=  malloc(1 * sizeof(Node));        // Creat local tree structure
  local_tree->gval = NULL;
  local_tree->parent = NULL;
  local_tree->attribute = NULL;
  local_tree->border_idx = NULL;
  local_tree->area_copy = NULL;
  local_tree->gval_par = NULL;
  return local_tree;
}
void free_tree(Node *tree) {
    if (tree == NULL) return;

    free(tree->gval);
    free(tree->attribute);
    free(tree->parent);
    free(tree->border_idx);
    free(tree->gval_par);
    free(tree->area_copy);
    // Free the node itself
    free(tree);
}
Output *init_output_struct(){
  Output *output 	=  malloc(1 * sizeof(Output));        // Creat local tree structure
  output->filtered = NULL;
  output->contrast = NULL;
  output->scale = NULL;
  output->luminance = NULL;
  return output;
}

void free_output(Output *output) {
    if (output == NULL) return;

    free(output->filtered);
    free(output->contrast);
    free(output->scale);
    free(output->luminance);
    output->filtered = NULL;
    output->contrast = NULL;
    output->scale = NULL;
    output->luminance = NULL;
    // Free the node itself
}

void initialize_attribute_functions(Node *tree, int attr_idx) {
    tree->size_attr = AttribsArray[attr_idx].size;
    new_aux_data = AttribsArray[attr_idx].new_data;
    init_aux_data = AttribsArray[attr_idx].init_data;
    delete_aux_data = AttribsArray[attr_idx].delete_data;
    add_to_aux_data = AttribsArray[attr_idx].add_to_data;
    merge_aux_data = AttribsArray[attr_idx].merge_data;
    merge_to_aux_data = AttribsArray[attr_idx].merge_to_data;
    clone_aux_data = AttribsArray[attr_idx].clone_data;
    create_mpi_aux_data = AttribsArray[attr_idx].create_mpi_data;
    create_mpi_aux_data();
}

Node *correct_relations_and_attributes(Node *tree, Arguments *args, int build_type)
{
  if (build_type == PARENTS)
  {
    tree = correct_borders(args, tree, PARENTS);
    timing("Nested relations corrected: wallclock time = %0.2f",
           (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
  }
  else if (build_type == ATTRIBUTES)
  {
    tree = correct_borders(args, tree, ATTRIBUTES);
    timing("Leaves attributes corrected: wallclock time = %0.2f",
           (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
  }
  else
  {
    tree = correct_borders(args, tree, PARENTS_AND_ATTRIBUTES);
    MPI_Barrier(MPI_COMM_WORLD);
    timing("Nested relations and leaves attributes corrected: wallclock time = %0.2f",
           (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
  }

  return tree;
}

void loop_operations(Node *tree, Arguments *args, Output *output)
  {
    for (int i = 0; i < args->num_operations; i++)
    {
      info("Starting operation #%d:\033[0m", i + 1);
      info("\tType: %s", args->ope_list[i].name);
      info("\tAttribute choice: %s", AttribsArray[args->ope_list[i].attribute_idx].name);

      if (i == 0)
      {
        initialize_attribute_functions(tree, args->ope_list[i].attribute_idx);
        build_tree_with_attributes(args, tree);
        MPI_Barrier(MPI_COMM_WORLD);
        timing("Local tree built: wallclock time = %0.2f",
               (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));

        if (np() > 1)
        {
          if (tree->gval_par == NULL)
          {
            for (int j = 0; j < args->num_operations; j++)
            {
              if (strncmp(args->ope_list[j].name, "pattern", strlen("pattern")) == 0)
              {
                tree->gval_par = calloc(tree->size_curr, sizeof(value));
                check_alloc(tree->gval_par, 1);
                tree->area_copy = malloc(tree->size_curr * sizeof(double));
                check_alloc(tree->area_copy, 2);
#pragma omp parallel for
                for (ulong k = 0; k < tree->size_curr; k++)
                {
                  tree->area_copy[k] = is_levelroot(tree, k) ? (*AttribsArray[args->ope_list[0].attribute_idx].area)(tree->attribute + k * tree->size_attr) : -10;
                  if (tree->parent[k] != BOTTOM)
                    tree->gval_par[k] = tree->gval[tree->parent[k]];
                }
                break;
              }
            }
          }
          tree = correct_relations_and_attributes(tree, args, 2);
        }
      }
      else if (AttribsArray[args->ope_list[i].attribute_idx].group != AttribsArray[args->ope_list[i - 1].attribute_idx].group)
      {
        debug("%s is from a new group, reflooding the existing tree with the new attributes", AttribsArray[args->ope_list[i].attribute_idx].name);

        free(tree->attribute);
        tree->attribute = NULL;
        MPI_Type_free(&mpi_attribute_type);
        initialize_attribute_functions(tree, args->ope_list[i].attribute_idx);

        flood_attributes(tree);
        timing("Re-building attributes: wallclock time = %0.2f",
               (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
        if (np() > 1)
        {
          tree = correct_relations_and_attributes(tree, args, 1);
        }
      }

      else
      {
        info("%s is in the same group, no need to rebuild the tree", AttribsArray[args->ope_list[i].attribute_idx].name);
      }

      perform_operation(args, tree, &(args->ope_list[i]), output);

    }
  }
  
  
  
  
  /* +++++++++++++++++++++++++++++++ */
  /*     	     Main Function         */
  /* +++++++++++++++++++++++++++++++ */

int main(int argc, char **argv)
  {

    init_mpi();

    Arguments args;


    parse_args(argc, argv, &args); // Read input arguments //
    create_mpi_value_type();       // Create MPI struct for value type //
    create_mpi_boundnode_type();   // Create MPI struct for boundary nodes //
    create_mpi_borderindex_type(); // Create MPI struct for borderindex struct //

    Node *local_tree = init_local_tree(); // Creat local tree structure
    Output *output = init_output_struct();

    read_input(&args, &(local_tree->gval));
    local_tree->size_curr = data_properties.dims_process[0] * data_properties.dims_process[1] * data_properties.dims_process[2];
    preprocessing(&args, local_tree->gval, local_tree->size_curr);
    MPI_Barrier(MPI_COMM_WORLD);

    info("Starting clock");
    start = times(&tstruct);

    if (args.interactive)
    {
      info("******* Interactive mode *******");
      info("");
      info("Well unlicky you, I have removed this feature for now...");
      // do more stuff
    }
    else if (args.num_operations > 0)
    {
      loop_operations(local_tree, &args, output);
    }
    else
    {
      if (args.save_output)
        write_no_operation(&args, local_tree->gval, NULL);
    }


  free(output);
  free_tree(local_tree);
  MPI_Type_free( &mpi_bound_node_type );
  MPI_Type_free( &mpi_borderindex_type );
  finalize_mpi();
    
  return 0;
} /* main */


/*FORMER INTERACTIVE VERSION REMOVED*/
    /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    /*         Process the tree using the interactive approach       */
    /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    /*
    else {
      /* +++++++++++++++++++++++++++ */
    /*     Flood attributes in     */
    /* +++++++++++++++++++++++++++ */
    /*
    int attrib;
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
      timing("Local tree attributes computed: total wallclock time = %0.2f",
	     (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    
      //
      copy_attr	= NULL; 	       	// Attributes copy (pattern spectrum)  // 
      gvals_par	= NULL;			// Parents node intensities (pattern spec) //
  
      if (np() > 1){
	if (!strcmp(args.operation_arg, "pattern")) {
     
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
  
/*
      //        CASE 1: Filtering   or extracting     //
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

	if (!strcmp(args.operation_arg, "filter") || !strcmp(args.operation_arg, "extract")) {
	  value *out_filter = calloc(local_tree->size_init, sizeof(value));   check_alloc(out_filter, 11);
	  tree_filtering(local_tree, out_filter, local_tree->size_init, args.decision_arg, attrib, args.lambda_arg);

	  if(args.saveout_arg){
	    if ((!strcmp(args.tree_type, "min") && !strcmp(args.morphology_choice, "opening")) ||
		(!strcmp(args.tree_type, "max") && !strcmp(args.morphology_choice, "closing"))) {
#pragma omp parallel for
	      for (ulong i = 0; i < local_tree->size_init; i++) 
		out_filter[i] = csl - out_filter[i];
	    }

	    if ( !strcmp(args.operation_arg, "extract")) {
#pragma omp parallel for
	      for (ulong i = 0; i < local_tree->size_init; i++) 
		out_filter[i] = local_tree->gval[i] - out_filter[i];
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

	else if (!strcmp(args.operation_arg, "csl")) {
      
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

	else if (!strcmp(args.operation_arg, "pattern")) {
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
	  info("No or wrong operation chosen for the interactive mode");

      
	lambda_vector_delete(lvec);       
      } while (args.lambda_arg != -1);
      free(local_tree->attribute);
      MPI_Type_free(&mpi_attribute_type);
      // clear_aux_data_store(local_tree->store[0]);
    }while(attrib != -1);
  }
  
  info("******* End of DISCCOFAN, cleaning ... *****\n\n");
*/