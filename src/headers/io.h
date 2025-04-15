#ifndef IO_H
#define IO_H

/* Read */
void read_input(Arguments *args, value **img);
/* Write */
// void write_output(Arguments* args, value* img, value* outOrig, value* outDH, value* outScale, double *spectrum,
//                   Operation* ope_cur, ulong dims_T[3], ulong dims[3], bool border[6], LambdaVec *lvec_1, LambdaVec *lvec_2);
// void write_output(Arguments* args, Node *tree, Output *output,
//                   Operation* ope_cur, ulong dims_T[3], ulong dims[3],  LambdaVec *lvec_1, LambdaVec *lvec_2);

void write_no_operation(Arguments* args, value* img, Operation* ope_cur);
void write_filtered(Arguments* args, value* img, Operation* ope_cur);
void write_pattern_spectra(Arguments *args, Operation *ope_cur, double* spectrum, LambdaVec *lvec);
void write_pattern_spectra2d(Arguments *args,  Operation *ope_cur, double* spectrum, LambdaVec *lvec_attr1, LambdaVec *lvec_attr2);
void write_check_files(Arguments* args, value* img, Operation* ope_cur);
void write_csl(Arguments *args, value *contrast, value *scale, value *luminance, Operation *ope_cur);
void write_dmp(Arguments *args, value *out, int num_lambdas, Operation *ope_cur) ;
void write_tree(Arguments *args, Node *tree, Operation *ope_cur);
#endif
