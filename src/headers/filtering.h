#ifndef FILTERING_H
#define FILTERING_H


/* Filtering  */
void tree_filtering(Node* tree, value *out, ulong size_tile,  int decision_choice, int attrib_choice, double lambda);
void tree_attribute_check(Node *tree, value *out, ulong lwb, ulong upb, double (*attribute)(void *));

#endif
