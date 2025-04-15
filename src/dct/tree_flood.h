#ifndef TREE_FLOOD_H
#define TREE_FLOOD_H

#define MXT_HISTO_SZ_LOG2 8
#define MXT_HISTO_SZ (1U << MXT_HISTO_SZ_LOG2)
#define MXT_HISTO_MASK (MXT_HISTO_SZ - 1)
#define NUM_DIGITS (sizeof(value) * CHAR_BIT + MXT_HISTO_SZ_LOG2 - 1) / MXT_HISTO_SZ_LOG2

/* Flooding */
idx *build_local_tree(Arguments *args, Node *tree, ulong *dims, ulong *attr_off);
long tree_flood_sal(Node *tree, AuxDataStore *store, Queue *q,  idx *levelroot, bool *reached, ulong *dims, ulong *grid,  ulong lwb, ulong upb, long level, int connectivity, void **thisattr);

void tree_flood_wil(Node *tree, AuxDataStore *store, pQueue *queue, pStack *stack, ulong *dims, ulong *attr_off, ulong lwb, ulong upb,  ulong min_idx, int connectivity);
void tree_flood_tee(Node *tree, AuxDataStore *store,  bool *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims, ulong *attr_off, ulong lwb, ulong upb,  int connectivity, int periodic);

/* Threads */
void fuse_sections(Node *tree, AuxDataStore *store,  ulong *dims, uint id, uint i, uint n_t, int connectivity);
void fuse_parents(Node *tree, ulong *dims, uint id, uint i, uint n_t, int connectivity);
void fuse_attributes(Node *tree, AuxDataStore *store, ulong *dims, uint id, uint i, uint n_t, int connectivity);

void merge_nodes(Node *tree,  AuxDataStore *store, idx x, idx y);

/* Annexes */
int bit_scan_reverse(ulong value); 
int bits_per_word_log2(void);
int bits_per_word(void); 
/* Mappings */
void create_mappings(value *gvals, ulong *ranks, ulong *ranks_inv, ulong size, ulong lwb, ulong upb);
void gen_histogram(value *gvals, ulong histos[NUM_DIGITS][MXT_HISTO_SZ], ulong lwb, ulong upb);
void exclusive_sum(ulong *it, ulong *it_end);
void create_ranks_inv(value *gvals, ulong *ranks_inv, ulong histos[NUM_DIGITS][MXT_HISTO_SZ], ulong size, ulong lwb, ulong upb);
void scatter_first_digit(value *gvals, SortItem *pair_it, ulong* histo,  ulong lwb, ulong upb);
void scatter_digit(SortItem *in, SortItem *out, int digit_nr, ulong* histo, ulong size);
void scatter_last_digit(SortItem *in, ulong *out, int digit_nr, ulong* histo, ulong size);
/* Bit Array */
BitArray *create_bit_array(ulong size);
void bit_array_set(ulong *data, ulong index);
bool bit_array_get(ulong *data, ulong index);
void bit_array_free(BitArray *bit_array);
value_t transform(const float val);
/* Misc */
int get_neighbors(ulong *dims, ulong lwb, ulong upb, ulong *neighbors, ulong p, ulong x, ulong y, ulong z,  int connectivity);
bool check_neighbor(bool *visited, PrioQueue *q, ulong *ranks, ulong* index, ulong* x, ulong* y, ulong *z, ulong* rank, ulong *dims, ulong n_x, ulong n_y, ulong n_z, ulong lwb);
void remaining(bool *visited, PrioQueue *q, ulong *ranks, ulong* index, ulong* rank, ulong* dims, ulong lwb, ulong upb, int connectivity, int periodic);
bool is_border(bool border[6], ulong *dims, ulong p);
bool is_levelroot(Node *tree, idx x);
idx get_levelroot(Node *tree, idx x);
idx levelroot(Node *tree, idx x);
idx get_parent(Node *tree, idx x);
void free_tree(Node *tree, ulong size_old);
ulong *sort_image_pixels(Arguments *args, Node *tree, ulong *dims);
#endif
