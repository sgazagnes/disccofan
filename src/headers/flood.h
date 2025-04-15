#ifndef FLOOD_H
#define FLOOD_H

/* Flooding */
void build_tree_with_attributes(Arguments *args, Node *local_tree);
void build_tree_without_attributes(Arguments *args, Node *local_tree);
void flood_attributes( Node *local_tree);

/* Misc */
NO_INLINE void create_mappings(value *gvals, ulong **ranks_inv, ulong lwb, ulong upb);

inline bool is_border(const bool border[6], const ulong *dims, ulong p)
{
    // Precompute dimensions products
    ulong xy_size = dims[0] * dims[1];  // Number of elements in one z-slice
    ulong x = p % dims[0];              // x-coordinate
    ulong y = (p % xy_size) / dims[0];  // y-coordinate
    ulong z = p / xy_size;              // z-coordinate

    // Check if the point is on any border
    return ((x == 0 && border[0]) || // Left border
            (x == dims[0] - 1 && border[1]) || // Right border
            (y == 0 && border[2]) || // Bottom border
            (y == dims[1] - 1 && border[3]) || // Top border
            (z == 0 && border[4]) || // Front border
            (z == dims[2] - 1 && border[5])); // Back border
}

inline bool is_levelroot(Node *tree, idx x)
{
  return ((tree->parent[x] == BOTTOM) || (tree->gval[x] != tree->gval[tree->parent[x]]));
} /* is_levelroot */
idx get_levelroot(Node *tree, idx x);
idx levelroot(Node *tree, idx x);
idx get_parent(Node *tree, idx x);

#endif
