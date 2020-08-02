#include "types.h"
#include "arguments.h"

Arguments *parse_args(int argc, char** argv, Arguments *args) {
  if (cmdline_parser (argc, argv, args) != 0) {
    error("Error while parsing command-line arguments: %s", *argv);
    MPI_Abort(MPI_COMM_WORLD, 001);
  }

  if(rank()==0 && args->grid_arg[0]*args->grid_arg[1]*args->grid_arg[2] != np()){
    error("This program assumes an equal amount of processes to the amount of grid tiles");
    MPI_Abort(MPI_COMM_WORLD, 001);
  }

  if(rank()==0 && (args->attribute_arg < -1 || args->attribute_arg >= NUMATTR)){
    error("Wrong attribute choice. Please write -1 for interactive attribute choice, or choose one in the following list:");
    for (int i=0; i<NUMATTR; i++)
      error("\t%d - %s", i, AttribsArray[i].name);
    MPI_Abort(MPI_COMM_WORLD, 001);
  }

  if(rank()==0 && (args->decision_arg < 0 || args->decision_arg >= NUMDECISIONS)){
    error("Wrong pruning choice. Please choose one in the following list:");
    for (int i=0; i<NUMDECISIONS; i++)
      error("\t%d - %s", i, Decisions[i].name);
    MPI_Abort(MPI_COMM_WORLD, 001);
  }

  set_verbosity(args->verbosity_arg);

  if (args->outtype_arg == NULL)
    asprintf(&args->outtype_arg, "%s", args->intype_arg);
  if (args->outfile_orig == NULL){
    if(args->infile_arg == 3)
      args->outfile_arg = 0;
    else if(args->infile_arg == 4)
      args->outfile_arg = 1;
    else
      args->outfile_arg = args->infile_arg;
  }
  if(np() == 1)
    args->outfile_arg = 0;

  args->threads_arg = args->threads_orig != NULL ?
    MAX(MIN(MAXTHREADS, args->threads_arg),1): omp_get_max_threads();
  omp_set_num_threads(args->threads_arg);

    
  return args;
} /* parse_args */


void check_bytes(void){
  debug("*** CHECK VARIABLE SIZES ***");
  debug("Byte size of types:\nSize of char: %ld\n size of short: %ld\n size of int: %ld\n size of long %ld\n size of long long %d\n size of float: %ld\n size of double: %ld", sizeof(char), sizeof(short), sizeof(int), sizeof(long), sizeof(long long), sizeof(float), sizeof(double));
  debug("Max value of types:\n char: %hhu\n short: %hu\n int: %d\n long %lu\n long long %llu\n float: %e\n double: %le", CHAR_MAX, SHRT_MAX, INT_MAX, LONG_MAX, LLONG_MAX, FLT_MAX, DBL_MAX);
  debug("Min value of types:\n char: %hhi\n short: %hi\n int: %i\n long %li\n long long %lli\n float: %e\n double: %le", CHAR_MIN, SHRT_MIN, INT_MIN, LONG_MIN, LLONG_MIN, -FLT_MAX, -DBL_MAX);
  debug("*** END OF CHECK ***\n");
}

void print_args(Arguments *args, ulong* dims) {
  debug("*** PRINTING ARGUMENTS ***");
  args->infile_arg == 0 ?
    debug("Input image:  %s.%s", args->inprefix_arg, args->intype_arg):
    debug("Input images:  %s-[X].%s", args->inprefix_arg, args->intype_arg);

  if(args->saveout_arg)
    args->outfile_arg == 0 ?
      debug("Output image: %s.%s", args->outprefix_arg, args->outtype_arg):
      debug("Output images: %s-[X].%s", args->outprefix_arg, args->outtype_arg);
  else
    debug("No output files");
  
  debug("Tile (process 0): size %ld pixels (width %ld, height %ld, depth %ld)",  dims[0]*dims[1]*dims[2], dims[0], dims[1], dims[2]);
  debug("Number of processes: %d (MPI: %d, threads: %d)", np()*args->threads_arg, np(), args->threads_arg);
  debug("Grid division: %d x %d x %d", 	args->grid_arg[0], args->grid_arg[1],  args->grid_arg[2]);
  debug("Bit depth: %d bit per pixel (%s)", args->bpp_arg, args->refine_arg ? "Data was transformed internally" : "Original dataset");
  debug("Type of tree: %s-tree",  	args->tree_arg);
  debug("Morphological operation: %s",   args->morphology_arg);
  debug("Attribute choice: %s",     	AttribsArray[args->attribute_arg].name);
  debug("Pruning decision choice: %s",   Decisions[args->decision_arg].name);
  debug("Type of operation: %s",      	args->filter_arg);
  debug("Flooding algorithm: Teeninga algorithm");
  debug("Connectivity: %d neighbours",   args->connectivity_arg);
  debug("Lambda: %0.2lf",       		args->lambda_arg);
  debug("lambda vector file: %s", 	args->lvec_arg);
  debug("Scaling of the lambda vector values: %0.2lf", 	args->imscale_arg);
  debug("Including background intensity in pattern spectra: %d", args->background_arg);
  debug("Verbosity: %s",     		args->verbosity_arg);
  // debug("Interactive mode: %d",     	args->interactive_arg);
  // debug("Save local tree: %d",     	args->savetree_arg);
  debug("*** END OF ARGUMENTS PRINTING, STARTING ***\n");
} /* print_args */
