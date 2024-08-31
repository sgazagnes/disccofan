#ifndef FILTER_H
#define FILTER_H


/* Filtering  */
void tree_filtering(Node* tree, value *out, ulong size_tile,  int decision_choice, int attrib_choice, double lambda);
void tree_filter_test(Node *tree, float *out, ulong lwb, ulong upb, double (*attribute)(void *));
/* Diff profile */
void tree_differential(Node *tree, ulong size_tile, LambdaVec *lvec, value *out_dh, value *temp_dh, value *out_orig, value *out_scale, value *temp_scale, bool *temp_valid, double (*attribute)(void *));

void combine_results(Node *tree,ulong size, LambdaVec *lvec, value *out_dh, value *out_dh2, value *out_orig, value *out_orig2, value *out_scale, value *out_scale2);
/* Pattern spec */
void tree_pattern_spectrum(Node *tree,  ulong size, LambdaVec *lvec, double *copy_attr, value *gvals_par, double* spectrum, int background, double (*area)(void *), double (*attribute)(void *));

void tree_pattern_spectrum2d(Node *tree,  ulong size, LambdaVec *lvec_attr1,LambdaVec *lvec_attr2, double *copy_attr, value *gvals_par, double* spectrum, int background,  double (*area)(void *), double (*attribute1)(void *), double (*attribute2)(void *));
#endif
