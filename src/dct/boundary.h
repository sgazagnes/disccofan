#ifndef BOUNDARY_H
#define BOUNDARY_H

Node *correct_borders(Arguments *args, Node *local_tree,  ulong *dims);
Node *correct_borders_parents(Arguments *args, Node *local_tree, ulong *dims);
void correct_borders_att(Arguments *args, Node *local_tree, ulong *dims);

#endif
