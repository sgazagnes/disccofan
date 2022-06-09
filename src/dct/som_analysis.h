#ifndef SOM_ANALYSIS_H
#define SOM_ANALYSIS_H

/******************************************************************************/
/*                             Lambda Vector handling                         */
/******************************************************************************/

int *read_som_attributes(Arguments *args, ulong size,  char *prefix);
void som_filter(Node *tree, value *out, int *som_attr, int lambda);
#endif
