#ifndef BOUNDARY_H
#define BOUNDARY_H

Node *correct_borders(Arguments *args, Node *local_tree,  ulong *dims);
Node *correct_local_tree(Node *local_tree, Boundary *b);
void update_par_step(Boundary *c) ;
Boundary *update( Boundary *b, Boundary *c);
Boundary *update_branch(BorderIndex x, BorderIndex z, BorderIndex s);
void update_node_attribute(BorderIndex x, BorderIndex s);
Boundary *adding_node(BorderIndex x, BorderIndex s);
Boundary *combine(Boundary *a, Boundary *b, Direction d);
void b_add(Boundary *c,  Boundary *b, ulong s, ulong i);
void reset_border_idx(Boundary *b);
void merge(Boundary *a, Boundary *b, Direction d);
void merge_b_nodes(BorderIndex x, BorderIndex y);
BorderIndex idx_i(Boundary *b, ulong c, Direction d);
BorderIndex idx_j(Boundary *b, ulong c, Direction d);
Boundary *create_boundary(Node *local_tree, ulong *dims, int attrib_choice, int connectivity);
void add_side(Node *local_tree, Boundary *b, int side, ulong *length, ulong offset, ulong *offset_bd, long *increment, int connectivity) ;
idx *add_parent_qu(Node *tree, Boundary *b);

Boundary *add_ancestors(Node *local_tree, Boundary *b);
void node_to_bound(Node *local_tree, Boundary *b, ulong b_idx, ulong tree_idx);
void bound_to_tree(Node *local_tree, Boundary *b, ulong b_idx, ulong tree_idx);
Boundary *realloc_b(Boundary *b, ulong newsize);
void free_boundary(Boundary *b);

BorderIndex b_levelroot(BorderIndex bi);
BorderIndex b_id_levelroot(BorderIndex bi);
bool is_bottom(BorderIndex bi);
bool bi_equal(BorderIndex ai, BorderIndex bi);
BorderIndex b_parent_lr(BorderIndex bi);
BorderIndex b_parent(BorderIndex bi);
value b_gval(BorderIndex bi);
BoundaryNode b_node(BorderIndex bi);
int bi_gval_equal(value a, value b);


#endif
