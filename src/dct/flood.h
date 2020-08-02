#ifndef FLOOD_H
#define FLOOD_H

#define MXT_HISTO_SZ_LOG2 8
#define MXT_HISTO_SZ (1U << MXT_HISTO_SZ_LOG2)
#define MXT_HISTO_MASK (MXT_HISTO_SZ - 1)
#define NUM_DIGITS (sizeof(value) * CHAR_BIT + MXT_HISTO_SZ_LOG2 - 1) / MXT_HISTO_SZ_LOG2

/* Flooding */
void build_local_tree(Arguments *args, Node *tree, ulong *dims);
void build_local_tree_par(Arguments *args, Node *tree, ulong *dims);
void build_local_tree_att( Node *local_tree, ulong *ranks,  ulong *dims);
/* Misc */
NO_INLINE void create_mappings(value *gvals, ulong **ranks_inv, ulong lwb, ulong upb);

bool is_border(bool border[6], ulong *dims, ulong p);
bool is_levelroot(Node *tree, idx x);
idx get_levelroot(Node *tree, idx x);
idx levelroot(Node *tree, idx x);
idx get_parent(Node *tree, idx x);
void free_tree(Node *tree);
ulong *sort_image_pixels(Arguments *args, Node *tree);

#endif
