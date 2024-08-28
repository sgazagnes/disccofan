#include "types.h"
#include "operation.h"
#include "filtering.h"
#include "io.h"
#include "pattern.h"
#include "csl.h"
#include "lambdavec.h"


void perform_filter_or_extract(Node *tree, Operation *ope_cur, value *filtered, bool invert){
    info("Applying filtering");
    tree_filtering(tree, filtered, tree->size_init, ope_cur->decision_idx, ope_cur->attribute_idx, ope_cur->lambda);
    MPI_Barrier(MPI_COMM_WORLD);
    timing("Filtering: wallclock time = %0.2f",
            (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));

    if (!strcmp(ope_cur->name, "extract"))
    {
#pragma omp parallel for
        for (ulong i = 0; i < tree->size_init; i++)
            filtered[i] = tree->gval[i] - filtered[i];
    }
    if (invert)
    {
#pragma omp parallel for
        for (ulong i = 0; i < tree->size_init; i++)
            filtered[i] = data_properties.g_max_gval - filtered[i];
    }
}

void perform_csl(Node *tree, Operation *ope_cur, value *contrast, value *scale, value *luminance){
    
    info("Applying CSL segmentation");
    value *temp_scale = calloc(tree->size_curr, sizeof(value));   check_alloc(temp_scale, 7);
    value *temp_contrast    = calloc(tree->size_curr, sizeof(value));   check_alloc(temp_contrast,    8);
    bool  *temp_valid = calloc(tree->size_curr, sizeof(bool));    check_alloc(temp_valid, 9);

    tree_differential(tree, tree->size_init, ope_cur->lvec1, contrast, scale, luminance, temp_contrast, temp_scale, temp_valid,  AttribsArray[ope_cur->attribute_idx].attribute);
    MPI_Barrier(MPI_COMM_WORLD);
    timing("Part 1 CSL done: wallclock time = %0.2f",
    (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    free(temp_valid);
    free(temp_scale);
    free(temp_contrast);
}

void perform_check_operation(Node *tree, Operation *ope_cur, value *attributes){
    info("Filling pixels with attributes. This is for testing purposes. Make sure FLOAT_TYPE is 1");
    #pragma omp parallel 
    {
        int np_threads = omp_get_num_threads();
        int id	   = omp_get_thread_num();
        ulong lwb 	   = id*tree->size_init/np_threads;
        ulong upb	   = (id+1)*tree->size_init/np_threads;
        tree_attribute_check(tree, attributes, lwb, upb, AttribsArray[ope_cur->attribute_idx].attribute);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    timing("Filling image: wallclock time = %0.2f",
            (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
}

void perform_dmp_operation(Node *tree, Operation *ope_cur, value *filtered){
    info("Creating DMP");
    for (int i = 0; i < ope_cur->lvec1->num_lambdas-1; i++)
    {
        float lambda_cur = ope_cur->lvec1->lambdas[i+1];
        tree_filtering(tree, filtered + i * tree->size_init, tree->size_init, ope_cur->decision_idx, ope_cur->attribute_idx, lambda_cur);
#pragma omp parallel for
        for (int j = i*tree->size_init; j < (i+1)*tree->size_init; j++)
            filtered[j] = tree->gval[j%tree->size_init] - filtered[j];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    timing("DMP done: wallclock time = %0.2f",
            (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));

}
void perform_pattern_spectrum(Node *tree, Operation *ope_cur, double *spectrum, int dim){
    if(dim == 1){
        double *temp_pattern_spectrum  = calloc(ope_cur->lvec1->num_lambdas, sizeof(double)); check_alloc(temp_pattern_spectrum, 11);
        info("Computing the 1-D pattern spectrum");
        tree_pattern_spectrum(tree, tree->size_init, ope_cur->lvec1, tree->area_copy, tree->gval_par, temp_pattern_spectrum, 1, AttribsArray[ope_cur->attribute_idx].area, AttribsArray[ope_cur->attribute_idx].attribute);
        MPI_Reduce(temp_pattern_spectrum, spectrum,  ope_cur->lvec1->num_lambdas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        timing("Pattern spectrum built: wallclock time = %0.2f",
            (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
        free(temp_pattern_spectrum);
    } else{
        double *temp_2D_spectrum  = calloc(ope_cur->lvec1->num_lambdas*ope_cur->lvec2->num_lambdas, sizeof(double));
        check_alloc(temp_2D_spectrum, 11);

        info("Computing the 2D pattern spectrum");
        if (ope_cur->attribute_idx == 0 && ope_cur->attribute_idx_2D != 0)
            tree_pattern_spectrum2d(tree, tree->size_init, ope_cur->lvec1, ope_cur->lvec2, tree->area_copy, tree->gval_par, temp_2D_spectrum, 1, AttribsArray[ope_cur->attribute_idx].area, AttribsArray[ope_cur->attribute_idx_2D].area, AttribsArray[ope_cur->attribute_idx_2D].attribute);
        else if (ope_cur->attribute_idx_2D == 0 && ope_cur->attribute_idx != 0)
            tree_pattern_spectrum2d(tree, tree->size_init, ope_cur->lvec1, ope_cur->lvec2, tree->area_copy, tree->gval_par, temp_2D_spectrum, 1, AttribsArray[ope_cur->attribute_idx].area, AttribsArray[ope_cur->attribute_idx].attribute, AttribsArray[ope_cur->attribute_idx].area);
        else
            tree_pattern_spectrum2d(tree, tree->size_init, ope_cur->lvec1, ope_cur->lvec2, tree->area_copy, tree->gval_par, temp_2D_spectrum, 1, AttribsArray[ope_cur->attribute_idx].area, AttribsArray[ope_cur->attribute_idx].attribute, AttribsArray[ope_cur->attribute_idx_2D].attribute);
        MPI_Reduce(temp_2D_spectrum, spectrum, ope_cur->lvec1->num_lambdas*ope_cur->lvec2->num_lambdas, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        timing("2D Pattern spectrum built: wallclock time = %0.2f",
         (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
        free(temp_2D_spectrum);
    }
}


void perform_operation(Arguments *args, Node *tree, Operation *ope_cur, Output *output)
{
    bool invert = !strcmp(args->tree_type, "min") ? true : false;
    if (!strcmp(ope_cur->name, "filter") || !strcmp(ope_cur->name, "extract"))
    {
        output->filtered = malloc(tree->size_init * sizeof(value));
        check_alloc(output->filtered, 3);

        perform_filter_or_extract(tree, ope_cur, output->filtered, invert);

        if(args->save_output)
            write_filtered(args, output->filtered, ope_cur);

        free(output->filtered);
        output->filtered = NULL;
    }
    else if (!strcmp(ope_cur->name, "dmp"))
    {
        ope_cur->lvec1 = lambda_vector_read(ope_cur->lvec_name, ope_cur->imscale);
        output->filtered = malloc((ope_cur->lvec1->num_lambdas-1) * tree->size_init * sizeof(value));
        check_alloc(output->filtered, 3);

        perform_dmp_operation(tree, ope_cur, output->filtered);

        if (args->save_output)
        {
            write_dmp(args, output->filtered, ope_cur->lvec1->num_lambdas-1, ope_cur);
        }

        free(output->filtered);
        output->filtered = NULL;
        lambda_vector_delete(ope_cur->lvec1);
    }
    else if (!strcmp(ope_cur->name, "csl"))
    {
        ope_cur->lvec1     = lambda_vector_read(ope_cur->lvec_name, ope_cur->imscale);
        // print_lambda_vec(ope_cur->lvec1);
        output->scale     = calloc(tree->size_curr, sizeof(value));   check_alloc(output->scale,   4); //scale out_orig
        output->contrast  = calloc(tree->size_curr, sizeof(value));   check_alloc(output->contrast,     5); //contrast out_dh
        output->luminance  = calloc(tree->size_curr, sizeof(value));   check_alloc(output->luminance,  6); //luminance out_scale

        perform_csl(tree, ope_cur, output->contrast, output->scale, output->luminance);

        if(args->save_output){
            write_csl(args, output->contrast, output->scale, output->luminance, ope_cur);
        }

        info("WELL we miss part 2");
 
        lambda_vector_delete(ope_cur->lvec1);
    }
    else if (!strcmp(ope_cur->name, "pattern"))
    {
        ope_cur->lvec1     = lambda_vector_read(ope_cur->lvec_name, ope_cur->imscale);
        output->pattern_spec = calloc(ope_cur->lvec1->num_lambdas, sizeof(double)); check_alloc(output->pattern_spec , 10);

        perform_pattern_spectrum(tree, ope_cur, output->pattern_spec, 1);  
        if(args->save_output && rank() == 0){
            write_pattern_spectra(args, ope_cur, output->pattern_spec,  ope_cur->lvec1);
            timing("Pattern spectrum written, wallclock time = %0.2f",
                    (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
        }
        free(output->pattern_spec);
        lambda_vector_delete(ope_cur->lvec1);
    }
    else if (!strcmp(ope_cur->name, "pattern2D"))
    {
        if( AttribsArray[ope_cur->attribute_idx].group != AttribsArray[ope_cur->attribute_idx_2D].group && ope_cur->attribute_idx != 0 && ope_cur->attribute_idx_2D != 0){
            error("For a 2D pattern spectrum, both attributes must have belong to the same group. Make sure this is the case. In the meantime, we won't do anything.");
            return;
        }
        ope_cur->lvec1     = lambda_vector_read(ope_cur->lvec_name, ope_cur->imscale);
        // print_lambda_vec(ope_cur->lvec1);
        ope_cur->lvec2     = lambda_vector_read(ope_cur->lvec_name_2D, ope_cur->imscale_2D);
// print_lambda_vec(ope_cur->lvec2);
        output->pattern_spec = calloc(ope_cur->lvec1->num_lambdas*ope_cur->lvec2->num_lambdas, sizeof(double)); check_alloc(output->pattern_spec , 10);

        perform_pattern_spectrum(tree, ope_cur, output->pattern_spec, 2);  

        if(args->save_output && rank() == 0){
            write_pattern_spectra2d(args, ope_cur, output->pattern_spec, ope_cur->lvec1, ope_cur->lvec2);
            timing("2D Pattern spectrum written, wallclock time = %0.2f",
                (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
        }

        free(output->pattern_spec);
        lambda_vector_delete(ope_cur->lvec1);
        lambda_vector_delete(ope_cur->lvec2);
    }
    else if(!strcmp(ope_cur->name, "check")){
        output->filtered = malloc(tree->size_init * sizeof(value));
        check_alloc(output->filtered, 3);
        perform_check_operation(tree, ope_cur, output->filtered);

        if(args->save_output)
            write_check_files(args, output->filtered, ope_cur);

        free(output->filtered);
        output->filtered = NULL;
    } else if(!strcmp(ope_cur->name, "tree")){
        write_tree(args, tree, ope_cur);
    }
}