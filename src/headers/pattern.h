#ifndef PATTERN_H
#define PATTERN_H

void tree_pattern_spectrum(Node *tree,  ulong size, LambdaVec *lvec, double *copy_attr, value *gvals_par, double* spectrum, bool background, double (*area)(void *), double (*attribute)(void *));

void tree_pattern_spectrum2d(Node *tree,  ulong size, LambdaVec *lvec_attr1,LambdaVec *lvec_attr2, double *copy_attr, value *gvals_par, double* spectrum, int background,  double (*area)(void *), double (*attribute1)(void *), double (*attribute2)(void *));

#endif
