#ifndef CSL_H
#define CSL_H


void tree_differential(Node *tree, ulong size_tile, LambdaVec *lvec, value *out_dh, value *temp_dh, value *out_orig, value *out_scale, value *temp_scale, bool *temp_valid, double (*attribute)(void *));

void combine_results(Node *tree,ulong size, LambdaVec *lvec, value *out_dh, value *out_dh2, value *out_orig, value *out_orig2, value *out_scale, value *out_scale2);

#endif

