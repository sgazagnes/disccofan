#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libconfig.h>
#include "config.h"

char* get_executable_path(const char* argv0) {
    char* executable_path = (char*)malloc(strlen(argv0) + 1);
    if (executable_path == NULL) {
        perror("Memory allocation error");
        return NULL;
    }

    strcpy(executable_path, argv0);

    // Find the last occurrence of the directory separator '/'
    char* last_separator = strrchr(executable_path, '/');
    if (last_separator != NULL) {
        // Terminate the string at the last separator to remove the executable name
        *(last_separator + 1) = '\0';
    } else {
        // If no separator is found, the executable is in the current directory
        strcpy(executable_path, "./");
    }

    return executable_path;
}

// Function to set default values for parameters
void set_default_values(Arguments* args) {
    args->interactive = 0;
    strcpy(args->input_name, "input_image.png");
    strcpy(args->output_name, "output_image.png");
    strcpy(args->input_prefix, "auto");
    strcpy(args->output_prefix, "auto");
    strcpy(args->input_type, "auto");
    strcpy(args->output_type, "auto");
    strcpy(args->input_file_style, "auto");
    strcpy(args->output_file_style, "auto");
    strcpy(args->hdf5_dataset, "my_dataset");
    strcpy(args->image_options, "classic");
    args->save_output = 1;
    args->tile_overlap = 1;
    args->mpi_grid[0] = 1;
    args->mpi_grid[1] = 1;
    args->mpi_grid[2] = 1;
    args->pixel_dim[0] = 1;
    args->pixel_dim[1] = 1;
    args->pixel_dim[2] = 1;
    args->threads = 1;
    args->bpp = -32;
    strcpy(args->preprocessing, "None");
    strcpy(args->tree_type, "max");
    //args->morphology_choice = 0;
    args->pixel_connectivity = 4;
    args->include_background = 0;
    strcpy(args->verbosity, "INFO");
    args->save_tree = 0;
}

void extract_file_info(const char* filename, char* file_prefix, char* file_type, char* file_style) {
    const char* delimiter = ".";
    const char* suffix_delimiter = "-T";
    size_t filename_length = strlen(filename);

    // Find the last occurrence of the delimiter
    const char* last_dot = strrchr(filename, delimiter[0]);
    const char* last_suffix_delimiter = strstr(filename, suffix_delimiter);

    // Check if a delimiter was found and if the suffix contains "-["
    if (last_dot != NULL) {
        if(last_suffix_delimiter != NULL && last_suffix_delimiter < last_dot){
            size_t prefix_length = last_suffix_delimiter - filename;
            memset(file_prefix,0,256);
            strncpy(file_prefix, filename, prefix_length);
            strcpy(file_style, "tile");
        
            // Check if the suffix starts with "-[X]" (regardless of digits after "-")
            const char* suffix_start = last_suffix_delimiter + 1;
            const char* suffix_end = strstr(suffix_start, "T");
            if (suffix_end != NULL) {
                // Extract the file type (suffix)
                 strcpy(file_type, last_dot + 1);
            }
        } else{
            size_t prefix_length = last_dot - filename;
            memset(file_prefix,0,256);
            strncpy(file_prefix, filename, prefix_length);
            strcpy(file_style, "single");
            strcpy(file_type, last_dot + 1);
        }
    } else {
        if(last_suffix_delimiter != NULL && last_suffix_delimiter < last_dot){
            size_t prefix_length = last_suffix_delimiter - filename;
            memset(file_prefix,0,256);
            strncpy(file_prefix, filename, prefix_length);
            strcpy(file_style, "tile");
        } else{
            size_t prefix_length = last_dot - filename;
            memset(file_prefix,0,256);
            strncpy(file_prefix, filename, prefix_length);
            strcpy(file_style, "single");
        }
    }
    
}

void parse_int_array(config_setting_t* setting, int* array, int* array_size) {
    *array_size = config_setting_length(setting);
    if (*array_size > MAX_PARAMS_SIZE) {
        error("Error: Array size exceeds maximum limit.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *array_size; i++) {
        config_setting_t* element = config_setting_get_elem(setting, i);
        array[i] = config_setting_get_int(element);
    }
}

void parse_double_array(config_setting_t* setting, double* array, int* array_size) {
    *array_size = config_setting_length(setting);
    if (*array_size > MAX_PARAMS_SIZE) {
        error("Error: Array size exceeds maximum limit.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *array_size; i++) {
        config_setting_t* element = config_setting_get_elem(setting, i);
        array[i] = config_setting_get_float(element);
    }
}

void parse_double_int_array(config_setting_t* setting, double* array, int* array_size) {
    *array_size = config_setting_length(setting);
    if (*array_size > MAX_PARAMS_SIZE) {
        error("Error: Array size exceeds maximum limit.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *array_size; i++) {
        config_setting_t* element = config_setting_get_elem(setting, i);
        if (config_setting_type(element) == CONFIG_TYPE_FLOAT) {
            array[i] = config_setting_get_float(element);
        } else if (config_setting_type(element) == CONFIG_TYPE_INT) {
            array[i] = (double)config_setting_get_int(element);
        } 
    }
}

void parse_float_array(config_setting_t* setting, float* array, int* array_size) {
    *array_size = config_setting_length(setting);
    if (*array_size > MAX_PARAMS_SIZE) {
        error("Error: Array size exceeds maximum limit.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *array_size; i++) {
        config_setting_t* element = config_setting_get_elem(setting, i);
        array[i] = (float)config_setting_get_float(element);
    }
}

void parse_string_array(config_setting_t* setting, char array[MAX_PARAMS_SIZE][256], int* array_size) {
    *array_size = config_setting_length(setting);
    if (*array_size > MAX_PARAMS_SIZE) {
        error("Error: Array size exceeds maximum limit.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < *array_size; i++) {
        config_setting_t* element = config_setting_get_elem(setting, i);
        const char* str_value = config_setting_get_string(element);
        strncpy(array[i], str_value, 255);
        array[i][255] = '\0'; // Ensure null-terminated string
    }
}

// Function to parse the config.ini file
void parse_config_file(const char* filename, Arguments* args) {
    config_t cfg;
    config_set_options (&cfg,  CONFIG_OPTION_AUTOCONVERT);

    config_init(&cfg);
    if (!config_read_file(&cfg, filename)) {
        fprintf(stderr, "Error: %s:%d - %s\n", config_error_file(&cfg),
                config_error_line(&cfg), config_error_text(&cfg));
        config_destroy(&cfg);
        exit(EXIT_FAILURE);
    }

    // Get interactive
    config_lookup_bool(&cfg, "interactive", &args->interactive);

    // Get inname
    const char* inname;
    if (config_lookup_string(&cfg, "input_name", &inname)) {
        strcpy(args->input_name, inname);
        extract_file_info(args->input_name,args->input_prefix, args->input_type, args->input_file_style);
    }

    const char* outname;
    if (config_lookup_string(&cfg, "output_name", &outname)) {
        strcpy(args->output_name, outname);
        extract_file_info(args->output_name,args->output_prefix, args->output_type, args->output_file_style);
    }


    // Get dataset
    const char* dataset;
    if (config_lookup_string(&cfg, "hdf5_dataset", &dataset)) {
        strcpy(args->hdf5_dataset, dataset);
    }

    // Style image
    const char* combine;
    if (config_lookup_string(&cfg, "image_options", &combine)) {
        if(!strcmp(combine, "grayscale") || !strcmp(combine, "rgb_channels") || !strcmp(combine, "lofar") || !strcmp(combine, "sequence"))
          strcpy(args->image_options, combine);
        else{
            warn("Combine image is weirdly set up. This parameter is a mess anyway, so I'll put the default value");
            strcpy(args->image_options, "classic");
        }
    }


    // Get saveout
    config_lookup_bool(&cfg, "save_output", &args->save_output);

    // Get tile overlap
    config_lookup_int(&cfg, "tile_overlap", &args->tile_overlap);

    // Grid operation
    config_setting_t* setting = config_lookup(&cfg, "mpi_grid");
    if (setting) {
        int check_size;
        if (config_setting_is_array(setting)) {
            parse_int_array(setting, args->mpi_grid, &(check_size));
        } else {
            check_size = 1;
        }

        if(check_size != 2 && check_size != 3){
            error("Your grid array should either be with two values (2D image) or with three values (3D)");
            MPI_Abort(MPI_COMM_WORLD, 001);
        }
    }
    
    // Pixel dimension
    setting = config_lookup(&cfg, "pixel_dim");
    if (setting) {
        int check_size;
        if (config_setting_is_array(setting)) {
            parse_float_array(setting, args->pixel_dim, &(check_size));
        } else {
            check_size = 1;
            args->pixel_dim[0] = args->pixel_dim[1] = args->pixel_dim[2] = config_setting_get_float(setting);
        }
        if(check_size > 3){
            error("Your pixel dimensions should be 3 values at maximum");
            MPI_Abort(MPI_COMM_WORLD, 001);
        }
    }

    // Get threads
    config_lookup_int(&cfg, "threads", &args->threads);

    const char* preprocc;
    if (config_lookup_string(&cfg, "preprocessing", &preprocc)) 
        strcpy(args->preprocessing, preprocc);

    // Get pixel connectivity
    config_lookup_int(&cfg, "pixel_connectivity", &args->pixel_connectivity);

    // Get include background
    config_lookup_int(&cfg, "include_background", &args->include_background);

    // Get save tree option
    config_lookup_int(&cfg, "save_tree", &args->save_tree);

    // Get image stats option
   // config_lookup_int(&cfg, "image_stats", &args->image_stats);


    // Get tree type
    const char* tree_type;
    if (config_lookup_string(&cfg, "tree_type", &tree_type)) {
        strcpy(args->tree_type, tree_type);
    }

    const char* verbose;
    if (config_lookup_string(&cfg, "verbosity", &verbose)) {
        strcpy(args->verbosity, verbose);
    }

        // Image operation
    config_setting_t* operation_setting = config_lookup(&cfg, "operation");

    int num_operations = config_setting_length(operation_setting);
    if (num_operations > MAX_OPERATION_SIZE) {
        error("Error: Number of operations exceeds maximum limit.\n");
        exit(EXIT_FAILURE);
    }
    if (num_operations ==0) {
        warn("No operation selected. The program will only do the I/O\n");
        //exit(EXIT_FAILURE);
    }

    args->num_operations = num_operations;

    for (int i = 0; i < num_operations; i++) {

        setting = config_setting_get_elem(operation_setting, i);
        int num_params = config_setting_length(setting);
        if (num_params > MAX_PARAMS_SIZE) {
            error("Error: Number of parameters for an operation exceeds maximum limit.\n");
            exit(EXIT_FAILURE);
        }
        if (num_params <2 ) {
                error("Error: Number of parameters must at least 2 for all operations. Check config ini file for parameters requirements\n");
                exit(EXIT_FAILURE);
        }

        config_setting_t* param_setting = config_setting_get_elem(setting, 0);
        const char* param_value = config_setting_get_string(param_setting);     
        strcpy(args->ope_list[i].name,param_value);

        param_setting = config_setting_get_elem(setting, 1);
        int att_value = config_setting_get_int(param_setting);   
        args->ope_list[i].attribute_idx = att_value;
        if(!strcmp(args->ope_list[i].name,"filter") || !strcmp(args->ope_list[i].name,"extract")){

            if (num_params < 3 || num_params > 4) {
                error("Error: Number of parameters must be either 3 or 4 for filtering or extracting. Check config ini file for parameters requirements\n");
                exit(EXIT_FAILURE);
            }


            param_setting = config_setting_get_elem(setting, 2);
            if (config_setting_type(param_setting) == CONFIG_TYPE_FLOAT) {
                args->ope_list[i].lambda = config_setting_get_float(param_setting);
            } else if (config_setting_type(param_setting) == CONFIG_TYPE_INT) {
                args->ope_list[i].lambda  = (double)config_setting_get_int(param_setting);
            }

            if(num_params > 3){
                // Decision strategy set up
                param_setting = config_setting_get_elem(setting, 3);
                int dec_value = config_setting_get_int(param_setting);   
                args->ope_list[i].decision_idx = dec_value;
            }else{
                 args->ope_list[i].decision_idx = 0;
            }

        } else if(!strcmp(args->ope_list[i].name,"pattern") || !strcmp(args->ope_list[i].name,"csl") || !strcmp(args->ope_list[i].name,"dmp")){

            if (num_params < 3 || num_params > 4) {
                error("Error: Number of parameters must be either 3 or 4 for the dmp, csl, or 1-D pattern spectrum. Check config ini file for parameters requirements\n");
                exit(EXIT_FAILURE);
            }

            param_setting = config_setting_get_elem(setting, 2);
            const char* lvec_value = config_setting_get_string(param_setting);     

            strcpy(args->ope_list[i].lvec_name,lvec_value);

            if(num_params > 3){
                // Scaling factor lvec
                param_setting = config_setting_get_elem(setting, 3);
                if (config_setting_type(param_setting) == CONFIG_TYPE_FLOAT) {
                    args->ope_list[i].imscale = config_setting_get_float(param_setting);
                } else if (config_setting_type(param_setting) == CONFIG_TYPE_INT) {
                    args->ope_list[i].imscale  = (double)config_setting_get_int(param_setting);
                }
            }else{
                args->ope_list[i].imscale  = 1.0;
            }                       

            args->ope_list[i].decision_idx = 0;
        } else if(!strcmp(args->ope_list[i].name,"pattern2D")){


            if (num_params != 7) {
                error("Error: Number of parameters must be 7 for 2D pattern spectrum. Check config ini file for parameters requirements\n");
                exit(EXIT_FAILURE);
            }

            param_setting = config_setting_get_elem(setting, 2);
            const char* lvec_value = config_setting_get_string(param_setting);     
            strcpy(args->ope_list[i].lvec_name,lvec_value);
            args->ope_list[i].imscale_2D  = 1.0;
                // Scaling factor lvec
            param_setting = config_setting_get_elem(setting, 3);
            if (config_setting_type(param_setting) == CONFIG_TYPE_FLOAT) {
                args->ope_list[i].imscale = config_setting_get_float(param_setting);
            } else if (config_setting_type(param_setting) == CONFIG_TYPE_INT) {
                args->ope_list[i].imscale  = (double)config_setting_get_int(param_setting);
            }

            param_setting = config_setting_get_elem(setting, 4);
            args->ope_list[i].attribute_idx_2D = config_setting_get_int(param_setting);   

            param_setting = config_setting_get_elem(setting, 5);
            const char* lvec_value2 = config_setting_get_string(param_setting);     
            strcpy(args->ope_list[i].lvec_name_2D,lvec_value2);

                // Scaling factor lvec
            param_setting = config_setting_get_elem(setting, 6);
            if (config_setting_type(param_setting) == CONFIG_TYPE_FLOAT) {
                args->ope_list[i].imscale_2D = config_setting_get_float(param_setting);
            } else if (config_setting_type(param_setting) == CONFIG_TYPE_INT) {
                args->ope_list[i].imscale_2D  = (double)config_setting_get_int(param_setting);
            }
   
        }  else if(!strcmp(args->ope_list[i].name,"tree")){
            // Nothing to do

            if (num_params > 2) {
                error("Error: Number of parameters must be 2 for saving the tree. Check config ini file for parameters requirements\n");
                exit(EXIT_FAILURE);
            }

        } else if(!strcmp(args->ope_list[i].name,"check")){

            if (num_params > 2) {
                error("Error: Number of parameters must be 2 for saving the check operation. Check config ini file for parameters requirements\n");
                exit(EXIT_FAILURE);
            }
        } else{
            error("Operation not recognized. Currently available: \"filter\", \"extract\", \"dap\", \"csl\", \"pattern\", \"pattern2D\", \"tree\"");
            exit(EXIT_FAILURE);
        }
    }
    if (&cfg != NULL) {
        config_destroy(&cfg);
    }
}


void parse_operations(const char* optarg, Arguments* args) {
    char* p = (char*)optarg;
    int op_idx = 0;
    int state = 0;
    char operation[256] = {0};
    int op_char_idx = 0;

    while (*p) {
        if (*p == '(') {
            // Start parsing a new tuple, reset the state and character index
            //printf("Beg:%c %d\n", *p, state );

            state = 0;
            op_char_idx = 0;
            p++;
        } else if (*p == ',' || *p == ')') {
            // Null-terminate the current element string
            //operation[op_char_idx] = '\0';

            // Determine which field to populate based on the current state
            Operation* op = &(args->ope_list[op_idx]);
            switch (state) {
                case 0:
                    strncpy(op->name, operation, sizeof(op->name));
                    op->decision_idx = 0; //default
                    strcpy(op->lvec_name,"lvec.txt"); //default
                    op->imscale = 1.;
                    memset(operation, 0, sizeof(operation));
                    state++;
                    break;
                case 1:
                    op->attribute_idx = atoi(operation);
                    memset(operation, 0, sizeof(operation));
                    state++;
                    break;
                case 2:
                    if (!strcmp(op->name, "filter") || !strcmp(op->name, "extract")) {
                        op->lambda = atof(operation);
                    } else if (!strcmp(op->name, "pattern") || !strcmp(op->name, "csl") || !strcmp(op->name, "dmp") || !strcmp(op->name, "pattern2D")) {
                        strncpy(op->lvec_name, operation, sizeof(op->lvec_name));
                    }
                    state++;
                    break;
                case 3:
                    if (!strcmp(op->name, "filter") || !strcmp(op->name, "extract")) {
                        op->decision_idx = atoi(operation);
                    } else if (!strcmp(op->name, "pattern") || !strcmp(op->name, "csl") || !strcmp(op->name, "dmp") || !strcmp(op->name, "pattern2D")) {
                        op->imscale = atof(operation);
                    }
                    state++;
                    break;
                case 4:
                    if (!strcmp(op->name, "pattern2D")) {
                        op->attribute_idx_2D = atof(operation);
                    }
                case 5:
                    if (!strcmp(op->name, "pattern2D")) {
                        strncpy(op->lvec_name_2D, operation, sizeof(op->lvec_name));
                    }
                case 6:
                    if (!strcmp(op->name, "pattern2D")) {
                        op->imscale_2D = atof(operation);
                    }
            }
            // Move to the next state
            

            // Move past the comma or closing parenthesis

            if (*p == ')') {
                // End of the tuple, move to the next tuple
                state = 0;
                op_idx++;
            } 
            p++;
            op_char_idx = 0;
            memset(operation, 0, sizeof(operation));
        } else if (*p == ' ') {
            // Skip spaces
           // printf("Skip:%c %d\n", *p, state )
            p++;
        } else {
            // Append the current character to the current element string
            operation[op_char_idx++] = *p++;
        }
    }

    // Update the number of operations
    args->num_operations = op_idx-1;
}


void check_bytes(void){
  debug("*** CHECK VARIABLE SIZES ***");
  debug("Byte size of types:\nSize of char: %ld\n size of short: %ld\n size of int: %ld\n size of long %ld\n size of long long %d\n size of float: %ld\n size of double: %ld", sizeof(char), sizeof(short), sizeof(int), sizeof(long), sizeof(long long), sizeof(float), sizeof(double));
  debug("Max value of types:\n char: %hhu\n short: %hu\n int: %d\n long %lu\n long long %llu\n float: %e\n double: %le", CHAR_MAX, SHRT_MAX, INT_MAX, LONG_MAX, LLONG_MAX, FLT_MAX, DBL_MAX);
  debug("Min value of types:\n char: %hhi\n short: %hi\n int: %i\n long %li\n long long %lli\n float: %e\n double: %le", CHAR_MIN, SHRT_MIN, INT_MIN, LONG_MIN, LLONG_MIN, -FLT_MAX, -DBL_MAX);
  debug("*** END OF CHECK ***\n");
}

void check_args(Arguments *args){

    debug("Checking arguments");
    const char* allowed_operations[] = {"filter", "extract", "pattern", "csl", "dmp", "pattern2D", "tree", "check"};
    const int num_allowed_operations = sizeof(allowed_operations) / sizeof(allowed_operations[0]);

    if(rank()==0 && args->mpi_grid[0]*args->mpi_grid[1]*args->mpi_grid[2] != np()){
        error("This program assumes an equal amount of MPI processes to the amount of grid tiles. Currently you have %d MPI for %d tiles.",np(),args->mpi_grid[0]*args->mpi_grid[1]*args->mpi_grid[2]);
        MPI_Abort(MPI_COMM_WORLD, 001);
    }

    for(int j = 0; j<args->num_operations; j++){
        int is_allowed = 0;
        for (int i = 0; i < num_allowed_operations; i++) {
            if (strcmp(args->ope_list[j].name, allowed_operations[i]) == 0) {
                is_allowed = 1;
                break;
            }
        }
        if(!is_allowed){
            error("One operation is not recognized. Make sure all values belong to the following list:");
            for (int i = 0; i < num_allowed_operations; i++) {
                error("\t %s", allowed_operations[i]);
            }
            MPI_Abort(MPI_COMM_WORLD, 001);

        }
        int cur_attr = args->ope_list[j].attribute_idx;
        int cur_attr2D = 0;
        if(!strcmp(args->ope_list[j].name, "pattern2D")){
            cur_attr2D = args->ope_list[j].attribute_idx_2D;
        }
        if(rank()==0 && (cur_attr < -1 || cur_attr >= NUMATTR || cur_attr2D < -1 || cur_attr2D >= NUMATTR)){
            error("The attribute function choice from one operation is incorrect. Make sure all values belong to the following list:");
            for (int i=0; i<NUMATTR; i++)
                error("\t%d - %s", i, AttribsArray[i].name);
            MPI_Abort(MPI_COMM_WORLD, 001);
        }

        if(!strcmp(args->ope_list[j].name, "filter") || !strcmp(args->ope_list[j].name, "extract") || !strcmp(args->ope_list[j].name, "dmp")){
            int decision_idx = args->ope_list[j].decision_idx;
            if(rank()==0 && (decision_idx < 0 || decision_idx >= NUMDECISIONS)){
                error("The pruning decision strategy from one operation is incorrect. Make sure all values belong to the following list:");
                for (int i=0; i<NUMDECISIONS; i++)
                    error("\t%d - %s", i, Decisions[i].name);
                MPI_Abort(MPI_COMM_WORLD, 001);
            }
        }

    }
    
    if(args->threads> omp_get_max_threads()){
        warn("You choose too many threads compared to your system available processes. Make sure to set optimal values");   
    }
    omp_set_num_threads(args->threads);

    if(!(!strcmp(args->preprocessing, "None") || !strcmp(args->preprocessing, "ubyte") || !strcmp(args->preprocessing, "ushort") || !strcmp(args->preprocessing, "uint") || !strcmp(args->preprocessing, "log") || !strcmp(args->preprocessing, "log10") || !strcmp(args->preprocessing, "exp") || !strcmp(args->preprocessing, "sqrt") )){
        error("The pre-processing option you set (%s) is not available. Available options are: 'ubyte', 'ushort', 'uint', 'log', 'log10', 'exp', sqrt'.", args->preprocessing);
        MPI_Abort(MPI_COMM_WORLD, 001);
    }
    
    debug("Arguments checked, no error detected so far");



}

void *parse_args(int argc, char** argv, Arguments *args) {
    char* executable_path = get_executable_path(argv[0]);
    strcpy(args->folder_run,executable_path);
    free(executable_path);

    set_default_values(args);
    bool found_config = false;

    // Parse config file
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-config") == 0 && i + 1 < argc) {
            // Read config file and parse parameters from it
            found_config = true;
            strcpy(args->config_file, argv[i + 1]);
            parse_config_file(args->config_file, args);
            break;
        } 
    } 

    if(!found_config){
        strcpy(args->config_file, strcat(args->folder_run, "config.ini"));
        parse_config_file(args->config_file, args);
    }       

    // Overide with command line
    for (int i = 1; i < argc; i++) {
        if (((strcmp(argv[i], "-input_image") == 0) || (strcmp(argv[i], "-i") == 0)) && i + 1 < argc) {
            strncpy(args->input_name, argv[i + 1], 256); 
            extract_file_info(args->input_name, args->input_prefix, args->input_type, args->input_file_style);
            i++; 
        } else  if (((strcmp(argv[i], "-output_image") == 0) || (strcmp(argv[i], "-o") == 0)) && i + 1 < argc) {
            strncpy(args->output_name, argv[i + 1], 256); 
            extract_file_info(args->output_name, args->output_prefix, args->output_type, args->output_file_style);
            i++;     
        } else  if ((strcmp(argv[i], "-interactive") == 0) && i + 1 < argc) {
            args->interactive = atoi(argv[i + 1]); 
            i++;    
        } else  if ((strcmp(argv[i], "-input_type") == 0) && i + 1 < argc) {
            strncpy(args->input_type, argv[i + 1], 256); 
            i++; 
        }  else  if ((strcmp(argv[i], "-output_type") == 0) && i + 1 < argc) {
            strncpy(args->output_type, argv[i + 1], 256); 
            i++; 
        }  else  if ((strcmp(argv[i], "-hdf5_dataset") == 0) && i + 1 < argc) {
            strncpy(args->hdf5_dataset, argv[i + 1], 256); 
            i++; 
        }  else  if ((strcmp(argv[i], "-image_options") == 0) && i + 1 < argc) {
            strncpy(args->image_options, argv[i + 1], 256); 
            i++; 
        } else  if (((strcmp(argv[i], "-verbosity") == 0)|| (strcmp(argv[i], "-v") == 0)) && i + 1 < argc) {
            strncpy(args->verbosity, argv[i + 1], 256); 
            i++; 
        }  else  if ((strcmp(argv[i], "-save_output") == 0) && i + 1 < argc) {
            args->save_output = atoi(argv[i + 1]); 
            i++; 
        // }    else  if ((strcmp(argv[i], "-image_stats") == 0) && i + 1 < argc) {
        //     args->image_stats = atoi(argv[i + 1]); 
        //     i++; 
         }  else  if ((strcmp(argv[i], "-save_tree") == 0) && i + 1 < argc) {
            args->save_tree = atoi(argv[i + 1]); 
            i++; 
        }   else  if ((strcmp(argv[i], "-tile_overlap") == 0) && i + 1 < argc) {
            args->tile_overlap = atoi(argv[i + 1]); 
            i++; 
        }  else  if (((strcmp(argv[i], "-threads") == 0)|| (strcmp(argv[i], "-t") == 0)) && i + 1 < argc) {
            args->threads = atoi(argv[i + 1]); 
            i++; 
        } else  if (((strcmp(argv[i], "-connectivity") == 0)|| (strcmp(argv[i], "-c") == 0)) && i + 1 < argc) {
            args->pixel_connectivity = atoi(argv[i + 1]); 
            i++; 
        } else  if ((strcmp(argv[i], "-preprocessing") == 0) && i + 1 < argc) {
            strncpy(args->preprocessing, argv[i + 1], 256); 
            i++; 
        } else  if ((strcmp(argv[i], "-tree_type") == 0) && i + 1 < argc) {
            strncpy(args->tree_type, argv[i + 1], 256); 
            i++; 
        } else  if ((strcmp(argv[i], "-include_background") == 0) && i + 1 < argc) {
            args->include_background = atoi(argv[i + 1]); 
            i++; 
        }  else  if (((strcmp(argv[i], "-grid") == 0)|| (strcmp(argv[i], "-g") == 0)) && i + 1 < argc) {
            sscanf(argv[i + 1], "[%d,%d,%d]", &args->mpi_grid[0], &args->mpi_grid[1], &args->mpi_grid[2]);
            i++; 
        } else  if ((strcmp(argv[i], "-operation") == 0) && i + 1 < argc) {
            parse_operations(argv[i+1], args);
            i++; 
        } else {
            warn("Command line option not recognized!");
        }
    }

    // Remaining
    set_verbosity(args->verbosity);
    check_args(args);
    print_args(args);

} // parse_args

void print_args(Arguments *args) {
    info("\033[1m*** DISCCOFAN PARAMETERS LIST ***\033[0m");

    info("\033[1mInitial config file:\033[0m %s", args->config_file);
    info("\033[1mInteractive mode:\033[0m %s", args->interactive == 1 ? "Yes!" : "No");
    info("\033[1mInput image options:\033[0m %s", args->image_options);
    info("\033[1mInput image name:\033[0m %s", args->input_name);
    debug("\033[1mInput image prefix:\033[0m %s", args->input_prefix);
    debug("\033[1mInput image type:\033[0m %s", args->input_type);
    info("\033[1mInput image style:\033[0m %s", args->input_file_style);

    if (!(strcmp(args->input_file_style, "tile")))
    {
        info("\033[1mTiles include overlap?\033[0m %s", args->tile_overlap == 1 ? "Yes!" : "No");
    }

    if (!strcmp(args->input_type, "hdf5") || !strcmp(args->input_type, "h5"))
    {
        info("\033[1mHDF5 dataset:\033[0m %s", args->hdf5_dataset);
    }

    if (args->save_output == 1)
    {
        info("\033[1mOutput image name:\033[0m %s", args->output_name);
        debug("\033[1mOutput image prefix:\033[0m %s", args->output_prefix);
        debug("\033[1mOutput image type:\033[0m %s", args->output_type);
        info("\033[1mOutput image style:\033[0m %s", args->output_file_style);
    }
    else
    {
        info("\033[1mNot saving any outputs\033[0m");
    }

    if (args->save_tree == 1)
    {
        info("\033[1mWe are also saving the tree(s)\033[0m");
    }

    info("\033[1mTree to build:\033[0m %s", args->tree_type);
    info("\033[1mMPI grid (W, H, D):\033[0m %d, %d, %d", args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]);
    info("\033[1mNumber of threads:\033[0m %d", args->threads);

    info("\033[1mNumber of operations to perform:\033[0m %d", args->num_operations);
    for (int i = 0; i < args->num_operations; i++)
    {
        info("\033[1mOperation #%d:\033[0m", i + 1);
        info("\033[1m\tType:\033[0m %s", args->ope_list[i].name);
        info("\033[1m\tAttribute choice:\033[0m %s", AttribsArray[args->ope_list[i].attribute_idx].name);

        if (!strcmp(args->ope_list[i].name, "filter") || !strcmp(args->ope_list[i].name, "extract"))
        {
            info("\033[1m\tAttribute threshold value:\033[0m %lf", args->ope_list[i].lambda);
            info("\033[1m\tPruning the tree using decision strategy:\033[0m %s", Decisions[args->ope_list[i].decision_idx].name);
        }
        else if (!strcmp(args->ope_list[i].name, "pattern") || !strcmp(args->ope_list[i].name, "csl") || !strcmp(args->ope_list[i].name, "dmp"))
        {
            info("\033[1m\tThreshold array file:\033[0m %s", args->ope_list[i].lvec_name);
            info("\033[1m\tScaling factor:\033[0m %lf", args->ope_list[i].imscale);
        }
        else if (!strcmp(args->ope_list[i].name, "pattern2D"))
        {
            info("\033[1m\tAttribute choice 2:\033[0m %s", AttribsArray[args->ope_list[i].attribute_idx_2D].name);
            info("\033[1m\tThreshold array file 1:\033[0m %s", args->ope_list[i].lvec_name);
            info("\033[1m\tThreshold array file 2:\033[0m %s", args->ope_list[i].lvec_name_2D);
            info("\033[1m\tScaling factor 1:\033[0m %lf", args->ope_list[i].imscale);
            info("\033[1m\tScaling factor 2:\033[0m %lf", args->ope_list[i].imscale_2D);
        }
        info("");
    }

    info("\033[1mVerbosity Level:\033[0m %s", args->verbosity);

    info("\033[1m*** ***\033[0m");
    info("");

} //  print_args
