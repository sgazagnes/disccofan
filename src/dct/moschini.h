#ifndef MOSCHINI_H
#define MOSCHINI_H

ulong *sort_image_pixels(Arguments *args, Node *tree, ulong *dims);
value_t *calculate_quantized_image(Node *tree, ulong *sorted, ulong *px_start, ulong *px_end,int nthreads);
idx *build_quant_tree(Arguments *args, Node *tree,  ulong *dims);
idx *refine_tree(Arguments *args, Node *tree, ulong *dims, ulong *px_start, ulong *px_end, ulong *sorted);
idx get_parent_qu(Node *tree, idx x);
idx get_levelroot_qu(Node *tree, idx x);

#endif
