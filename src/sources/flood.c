
#include "types.h"
#include "attributes.h"
#include "flood.h"
#include "queue.h"
#include "mpihelper.h"
#include "workspace.h"

#define MXT_HISTO_SZ_LOG2 8
#define MXT_HISTO_SZ (1U << MXT_HISTO_SZ_LOG2)
#define MXT_HISTO_MASK (MXT_HISTO_SZ - 1)
#define NUM_DIGITS (sizeof(value) * CHAR_BIT + MXT_HISTO_SZ_LOG2 - 1) / MXT_HISTO_SZ_LOG2

void tree_flood(idx *parent, void *attribute, BitArray *visited, PrioQueue *q, ulong *ranks, ulong *ranks_inv, ulong *dims, ulong lwb, ulong upb, int connectivity, int size_att);

void create_ranks_inv(value *gvals, ulong **ranks_inv, ulong histos[][MXT_HISTO_SZ], ulong lwb, ulong upb, value_t (*functor)(const value));

BitArray *create_bit_array(ulong size, Allocator *work_alloc);
void bit_array_set(ulong *data, ulong index);
bool bit_array_get(ulong *data, ulong index);
void bit_array_free(BitArray *bit_array);

void fuse_sections(Node *tree, ulong *dims, uint id, uint i, uint n_t, int connectivity, int with_attributes);
void merge_nodes(Node *tree, idx x, idx y, int with_attributes);

int get_neighbors(ulong *dims, ulong lwb, ulong upb, ulong *neighbors, ulong p, ulong x, ulong y, ulong z, int connectivity);
bool check_neighbor(BitArray *visited, PrioQueue *q, ulong *ranks, ulong *index, ulong *x, ulong *y, ulong *z, ulong *rank, long offset, ulong n_x, ulong n_y, ulong n_z, ulong lwb);
void remaining(BitArray *visited, PrioQueue *q, ulong *ranks, ulong *index, ulong *rank, ulong *dims, ulong lwb, ulong upb, int connectivity);
value_t transform_float(const float val);
value_t transform_dum(const value val);

/* +++++++++++++++++++++++++++++++ */
/*				                         */
/*      Flooding functions         */
/*				                         */
/* +++++++++++++++++++++++++++++++ */

void init_locks(omp_lock_t *lock, int nthreads) {
    for (int i = 0; i < nthreads; i++) {
        omp_init_lock(&(lock[i]));
    }
}

void destroy_locks(omp_lock_t *lock, int nthreads) {
    for (int i = 0; i < nthreads; i++) {
        omp_destroy_lock(&(lock[i]));
    }
}

void init_saval(int *saval, int nthreads) {
    for (int i = 0; i < nthreads; i++) {
        saval[i] = 0;
    }
}


/* Generalized function to build local trees */
void build_local_tree_generic(Arguments *args, Node *local_tree, int with_attributes) {
    info("Building the tree");
    ulong *dims = data_properties.dims_process;
    int connectivity = args->pixel_connectivity;
    int nthreads = args->threads;
    int *saval = calloc(nthreads, sizeof(int));  // Thread completion flags

    omp_lock_t lock[nthreads];
    init_locks(lock, nthreads);
    init_saval(saval, nthreads);

    if (local_tree->parent != NULL || local_tree->attribute != NULL) {
        warn("Tree not properly cleaned");
    }

    local_tree->parent = NULL;
    ulong *ranks = NULL;

    if (nthreads > 1 && args->bpp > 16)
        warn("Using multiple threads with high dynamic range dataset might not work well");

    #pragma omp parallel num_threads(nthreads)
    {
        int id = omp_get_thread_num(); /* Thread number */
        ulong lwb = dims[2] == 1 ? dims[0] * ((id * dims[1]) / nthreads)
                                 : dims[0] * dims[1] * ((id * dims[2]) / nthreads); /* Lower bound for current thread */
        ulong upb = dims[2] == 1 ? dims[0] * (((id + 1) * dims[1]) / nthreads)
                                 : dims[0] * dims[1] * (((id + 1) * dims[2]) / nthreads); /* Upper bound for current thread */

        // Create mappings
        ulong *ranks_inv = NULL;
        create_mappings(local_tree->gval, &ranks_inv, lwb, upb);

#pragma omp single
        {
            ranks = malloc(local_tree->size_curr * sizeof(ulong));
            local_tree->parent = (idx *)ranks;
            if (with_attributes) {
                local_tree->attribute = malloc(local_tree->size_curr * local_tree->size_attr);
            }
        }

        for (ulong i = 0; i != upb - lwb; ++i) {
            ranks[ranks_inv[i] + lwb] = i;
        }

        debug("Mappings: wallclock time = %0.2f ",
              (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));

        // Start construction
        Allocator work_alloc;
        allocator_init(&work_alloc, workspace_construct(upb - lwb));
        PrioQueue *q = create_prio_queue(upb - lwb - 1, &work_alloc);
        BitArray *visited = create_bit_array(upb - lwb, &work_alloc);

        if (with_attributes) {
            // Initialize attributes only if specified
            init_attrib_array(local_tree->attribute, local_tree->gval,  dims, 
                              lwb, upb, local_tree->size_attr);
        }

        // Generalized flood function (attribute-based or parent-based)
        tree_flood(local_tree->parent, local_tree->attribute, visited, q, ranks + lwb, ranks_inv,
                   dims, lwb, upb, connectivity, with_attributes ? local_tree->size_attr : 0);

        free(ranks_inv);
        allocator_free(&work_alloc);

        // Tree compression
        for (ulong x = lwb; x < upb; x++) {
            idx r = x;
            value gv = local_tree->gval[x];
            while ((local_tree->parent[r] != BOTTOM) && (gv == local_tree->gval[local_tree->parent[r]])) {
                r = local_tree->parent[r];
            }
            idx i = x;
            while (i != r) {
                idx y = local_tree->parent[i];
                local_tree->parent[i] = r;
                i = y;
            }

            local_tree->parent[r] = get_levelroot(local_tree, local_tree->parent[r]);
        }

        // Thread merging (if needed)
        if (nthreads > 1) {
            int i = 1;
            int qq = id;
            while (id + i < nthreads && qq % 2 == 0) {
                while (saval[id + i] <= 0) {
#pragma omp flush
                }
                omp_set_lock(&lock[id + i]);
                saval[id + i]--;
                omp_unset_lock(&lock[id + i]);

                fuse_sections(local_tree, dims, id, i, nthreads, connectivity, with_attributes);

                i *= 2;
                qq /= 2;
            }
            if (id > 0) {
                omp_set_lock(&lock[id]);
                saval[id]++;
                omp_unset_lock(&lock[id]);
            }
        }
    }

    free(saval);
    destroy_locks(lock, nthreads);

    local_tree->size_init = local_tree->size_curr;
} /* build_local_tree_generic */

/* Wrapper functions for compatibility */
void build_tree_with_attributes(Arguments *args, Node *local_tree) {
    build_local_tree_generic(args, local_tree, 1); // with attributes
}

void build_tree_without_attributes(Arguments *args, Node *local_tree) {
    build_local_tree_generic(args, local_tree,  0); // without attributes
}


void flood_attributes(Node *local_tree)
{
  info("Refilling attributes on existing tree structures");
  ulong *ranks = NULL;
  create_mappings(local_tree->gval, &ranks, 0, local_tree->size_curr);
  ulong *dims = data_properties.dims_process;
  bool *visited = calloc(local_tree->size_curr, sizeof(bool));
  check_alloc(visited, 603);
  if(local_tree->attribute != NULL){
    warn("Badly freed memory for the attribute array");
  }
  local_tree->attribute = malloc(local_tree->size_curr * local_tree->size_attr);

  init_attrib_array(local_tree->attribute, local_tree->gval, dims, 0, local_tree->size_curr, local_tree->size_attr);
  for (long i = local_tree->size_curr - 1; i >= 0; i--)
  {
    idx p = ranks[i];
    visited[p] = true;
    idx parent = get_parent(local_tree, p);
    if (parent == BOTTOM)
      continue;

    merge_aux_data(local_tree->attribute + parent * local_tree->size_attr, local_tree->attribute + p * local_tree->size_attr);

    if (local_tree->gval[p] == local_tree->gval[parent] && visited[parent])
    {
      parent = get_parent(local_tree, parent);
      if (parent != BOTTOM)
      {
        merge_aux_data(local_tree->attribute + parent * local_tree->size_attr, local_tree->attribute + p * local_tree->size_attr);
      }
    }
  }
  free(visited);
  free(ranks);
}


void tree_flood(idx *parents, void *attribute, BitArray *visited, PrioQueue *q, ulong *ranks, ulong *ranks_inv, ulong *dims, ulong lwb, ulong upb, int connectivity, int size_att)
{

  ulong index = lwb;
  ulong rank = ranks[index - lwb];
  bit_array_set(visited->data, index - lwb);
  char *attribute_data = (char *)attribute;
  while (true)
  {
    remaining(visited, q, ranks, &index, &rank, dims, lwb, upb, connectivity);

    if (q->m_levels[0][0] == 0)
      break;
    rank = q->m_top;
    ulong parent = ranks_inv[rank] + lwb;
    prio_queue_remove(q);
    // Merge auxiliary data if attributes are present
    if(size_att > 0)
      merge_aux_data(attribute + parent * size_att, attribute + index * size_att);
    parents[index] = parent;
    index = parent;
  }
 // Finalize the last parent

  parents[index] = BOTTOM;
} /* tree_flood */

/* ++++++++++++++++++++++++++ */
/*      Threads functions     */
/* ++++++++++++++++++++++++++ */

void fuse_sections(Node *tree, ulong *dims, uint id, uint i, uint n_t, int connectivity, int with_attributes)
{
  ulong inc = dims[2] == 1 ? dims[0] : dims[0] * dims[1];
  ulong mdb = dims[2] == 1 ? dims[0] * (((id + i) * dims[1]) / n_t) : dims[0] * dims[1] * (((id + i) * dims[2]) / n_t);

  ulong *ranks_inv = NULL;
  create_mappings(tree->gval, &ranks_inv, mdb - inc, mdb + inc);
  ulong *ranks = malloc(2 * inc * sizeof(ulong));

  for (ulong ii = 0; ii != 2 * inc; ++ii)
  {
    ranks[ranks_inv[ii]] = ii;
  }

  Allocator work_alloc;
  allocator_init(&work_alloc, workspace_construct(2 * inc));
  PrioQueue *q = create_prio_queue(2 * inc - 1, &work_alloc);
  BitArray *visited = create_bit_array(2 * inc, &work_alloc);

  ulong index = mdb - inc;
  ulong rank = ranks[index - mdb + inc];
  bit_array_set(visited->data, index - mdb + inc);
  while (true)
  {
    remaining(visited, q, ranks, &index, &rank, dims, mdb - inc, mdb + inc, connectivity);
    if (q->m_levels[0][0] == 0)
      break;
    rank = q->m_top;
    ulong parent = ranks_inv[rank] + mdb - inc;
    prio_queue_remove(q);
    merge_nodes(tree, index, parent, with_attributes);
    index = parent;
  }

  allocator_free(&work_alloc);
  free(ranks);
  free(ranks_inv);
}

void merge_nodes(Node *tree, idx x, idx y, int with_attributes)
{
  void *cor = NULL;
  void *copa = NULL;
  void *attr = NULL;
  idx h, z;

  if (tree->gval[y] > tree->gval[x])
  {
    h = x;
    x = y;
    y = h;
  }
  x = get_levelroot(tree, x);
  y = get_levelroot(tree, y);

  while ((x != y) && (y != BOTTOM))
  {
    z = get_parent(tree, x);
    if ((z != BOTTOM) && (tree->gval[z] >= tree->gval[y]))
    {
      if (cor && with_attributes)
      {
        merge_aux_data(tree->attribute + x * tree->size_attr, cor);
      }
      x = z;
    }
    else
    {
      if(with_attributes){
      if (cor)
        merge_to_aux_data(NULL, &copa, tree->attribute + x * tree->size_attr, cor);
      else
        clone_aux_data(NULL, &copa, tree->attribute + x * tree->size_attr);
      clone_aux_data(NULL, &cor, tree->attribute + x * tree->size_attr);

      if (copa)
      {
        attr = tree->attribute + x * tree->size_attr;
        clone_aux_data(NULL, &attr, copa);
      }
    }
      tree->parent[x] = y;
      x = y;
      y = z;
    }
  }

  if (y == BOTTOM && cor && with_attributes)
  {
    while (x != BOTTOM)
    {
      merge_aux_data(tree->attribute + x * tree->size_attr, cor);

      x = get_parent(tree, x);
    }
  }
  if (cor)
    delete_aux_data(cor);
  if (copa)
    delete_aux_data(copa);
} /* merge_nodes */


/* +++++++++++++++++++++++++++++++ */
/*				                         */
/*         Annexes Functions       */
/*				                         */
/* +++++++++++++++++++++++++++++++ */

inline value_t transform_float(const float val)
{
    value_t valu;
    memcpy(&valu, &val, sizeof(val)); // Safely convert float to value_t

    if (valu & 0x80000000) // Check if the sign bit is set
        return 0xFFFFFFFF - valu; // Invert the value for negative floats

    return valu | 0x80000000; // Set the sign bit for positive floats
}

inline value_t transform_dum(const value val)
{
  return val;
}

void gen_histogram(value *gvals, ulong histos[][MXT_HISTO_SZ], ulong lwb, ulong upb, value_t (*functor)(const value))
{
  memset(histos, 0, NUM_DIGITS * MXT_HISTO_SZ * sizeof(histos[0][0]));

  for (ulong i = lwb; i != upb; ++i)
  {
    int shift = 0;
    value_t gval = functor(gvals[i]);
    for (int digit_nr = 0; digit_nr != NUM_DIGITS; ++digit_nr)
    {
      int digit = (gval >> shift) & MXT_HISTO_MASK;
      ++histos[digit_nr][digit];
      shift += MXT_HISTO_SZ_LOG2;
    }
  }
}

void exclusive_sum(ulong *it, ulong *it_end)
{
  ulong sum = 0;
  while (it != it_end)
  {
    ulong next = *it;
    *it++ = sum;
    sum += next;
  }
}

NO_INLINE void create_mappings(value *gvals, ulong **ranks_inv, ulong lwb, ulong upb)
{
  const int num_digits = NUM_DIGITS;
  value_t (*foo)(const value);

  if (FLOAT_TYPE == 1)
    foo = &transform_float;
  else
    foo = &transform_dum;

  ulong histos[num_digits][MXT_HISTO_SZ];
  gen_histogram(gvals, histos, lwb, upb, foo);

  for (int digit_nr = 0; digit_nr != num_digits; digit_nr++)
    exclusive_sum(histos[digit_nr], histos[digit_nr] + MXT_HISTO_SZ);

  create_ranks_inv(gvals, ranks_inv, histos, lwb, upb, foo);
}

void scatter_first_digit(value *gvals, SortItem *pair_it, ulong *histo, ulong lwb, ulong upb, value_t (*functor)(const value))
{
  for (ulong i = lwb; i != upb; ++i)
  {
    value_t gval;
    gval = functor(gvals[i]);
    int digit = gval & MXT_HISTO_MASK;
    pair_it[histo[digit]++] = (SortItem){gval, i - lwb};
  }
}

void scatter_digit_bit1(value *gvals, ulong *out, ulong *histo, ulong lwb, ulong upb, value_t (*functor)(const value))
{
  for (ulong i = lwb; i != upb; ++i)
  {
    int digit = functor(gvals[i]) & MXT_HISTO_MASK;
    out[histo[digit]++] = i - lwb;
  }
}

void scatter_digit(SortItem *in, SortItem *out, int digit_nr, ulong *histo, ulong size)
{
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i)
  {
    const SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair;
  }
}

void scatter_last_digit(SortItem *in, ulong *out, int digit_nr, ulong *histo, ulong size)
{
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i)
  {
    const SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair.rank;
  }
}

void create_ranks_inv(value *gvals, ulong **ranks_inv, ulong histos[][MXT_HISTO_SZ], ulong lwb, ulong upb, value_t (*functor)(const value))
{
  const int num_digits = NUM_DIGITS;

  if (num_digits == 1)
  {
    *ranks_inv = calloc(upb - lwb, sizeof(ulong));
    scatter_digit_bit1(gvals, *ranks_inv, histos[0], lwb, upb, functor);
    return;
  }

  // allocator_set_pos(work_alloc, alloc_pos);
  SortItem *pairs1 = malloc((upb - lwb) * sizeof(SortItem)); //*ranks_inv;
  // allocator_allocate(work_alloc, upb-lwb, sizeof(SortItem));
  SortItem *pairs2 = malloc((upb - lwb) * sizeof(SortItem)); //*ranks_inv + (upb-lwb)*sizeof(SortItem);
  // allocator_allocate(work_alloc, upb-lwb, sizeof(SortItem));

  SortItem *pairswap = NULL;
  scatter_first_digit(gvals, pairs2, histos[0], lwb, upb, functor);

  int digit_nr = 1;
  for (; digit_nr != NUM_DIGITS - 1; ++digit_nr)
  {
    pairswap = pairs1;
    pairs1 = pairs2;
    pairs2 = pairswap;

    scatter_digit(pairs1, pairs2, digit_nr, histos[digit_nr], (upb - lwb));
  }

  free(pairs1);
  *ranks_inv = calloc(upb - lwb, sizeof(ulong));
  scatter_last_digit(pairs2, (ulong *)*ranks_inv, digit_nr, histos[digit_nr], (upb - lwb));
  free(pairs2);
  // allocator_set_pos(work_alloc, alloc_pos_rank_inv);
}

/* +++++++++++++++++++++++++++++++ */
/*				                         */
/*       Bit Array Functions       */
/*				                         */
/* +++++++++++++++++++++++++++++++ */

BitArray *create_bit_array(ulong size, Allocator *work_alloc)
{
  BitArray *bit_array = malloc(sizeof(BitArray));
  check_alloc(bit_array, 507);
  bit_array->num_words = (size + bits_per_word_log2() - 1) / bits_per_word_log2();
  bit_array->data = (ulong *) allocator_allocate(work_alloc, bit_array->num_words, sizeof(ulong));
  memset(bit_array->data, 0, bit_array->num_words * sizeof(ulong));
  return bit_array;
}

void bit_array_set(ulong *data, ulong index)
{
  ulong word_idx = index >> bits_per_word_log2();
  int bit_idx = index & (bits_per_word() - 1);
  data[word_idx] |= (ulong)(1) << bit_idx;
}

bool bit_array_get(ulong *data, ulong index)
{
  ulong word_idx = index >> bits_per_word_log2();
  int bit_idx = index & (bits_per_word() - 1);
  return !!(data[word_idx] & ((ulong)(1) << bit_idx));
}

void bit_array_free(BitArray *bit_array)
{
  free(bit_array->data);
  free(bit_array);
}

int bit_scan_reverse(ulong val)
{
  return sizeof(val) * CHAR_BIT - 1 - __builtin_clzl(val);
}
int bits_per_word_log2(void)
{
  return sizeof(unsigned) * CHAR_BIT - 1 - __builtin_clz((sizeof(ulong) * CHAR_BIT));
}
int bits_per_word(void)
{
  return sizeof(ulong) * CHAR_BIT;
}


/* Misc */

bool check_neighbor(BitArray *visited, PrioQueue *q, ulong *ranks, ulong *index, ulong *x, ulong *y, ulong *z, ulong *rank, long offset, ulong n_x, ulong n_y, ulong n_z, ulong lwb)
{
// info("In check neig");
  ulong n = *index + offset;

  if (bit_array_get(visited->data, n - lwb))
    return false;

  bit_array_set(visited->data, n - lwb);

  ulong rank_n = ranks[n - lwb];
  if (*rank > rank_n)
  {
    insert_prio_queue(q, rank_n);
    return false;
  }

  insert_prio_queue(q, *rank);

  *index = n;
  *rank = rank_n;
  *x = n_x;
  *y = n_y;
  *z = n_z;

  return true;
}
/* Clarity improved but slower....
void remaining(BitArray *visited, PrioQueue *q, ulong *ranks, ulong *index, ulong *rank, ulong *dims, ulong lwb, ulong upb, int connectivity)
{
    ulong z = *index / (dims[0] * dims[1]);
    ulong xy = *index % (dims[0] * dims[1]);
    ulong x = xy % dims[0];
    ulong y = xy / dims[0];

  // Define delta values for 6, 8, and 26 connectivity
    long dx_val[26] = { 1,  0, -1, 0, 0,  0, 1,  1,  -1, -1,  0,  0,  0,  0, 1,  1,  1,  1, -1, -1, -1, -1, -1, 1, -1,  1}; // x deltas for 4 or 8 connectivity
    long dy_val[26] = { 0, -1,  0, 1, 0,  0, 1, -1,   1, -1,  1,  1, -1, -1, 1,  1, -1, -1,  1,  1, -1, -1,  0, 0,  0,  0}; // y deltas for 4 or 8 connectivity
    long dz_val[26] = { 0,  0,  0, 0, 1, -1, 0,  0,   0,  0,  1, -1,  1, -1, 1, -1,  1, -1,  1, -1,  1, -1,  1, 1, -1, -1}; // z deltas for 4 or 8 connectivity


    bool cond = true;

    while (cond)
    {
        cond = false;

       // Check direct neighbors
        for (int i = 0; i < 6; ++i)
        {
          long dx = dx_val[i];
          long dy = dy_val[i];
          long dz = dz_val[i];
          ulong new_index = *index + dx + dy * dims[0] + dz * dims[0] * dims[1];

          if ((x + dx >= 0 && x + dx < dims[0])&&
              (y + dy >= 0 && y + dy < dims[1]) &&
              (z + dz >= 0 && z + dz < dims[2]))
          {
            cond = cond || ((new_index >= lwb) && (new_index < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dx + dy * dims[0] + dz * dims[0] * dims[1], x + dx, y + dy, z + dz, lwb));
          }
        }



   // Check diagonal neighbors for connectivity >= 8
        if (connectivity >= 8)
        {
            for (int i = 6; i < 10; ++i)
            {
                long dx = dx_val[i];
                long dy = dy_val[i];
                long dz = dz_val[i];
                ulong new_index = *index + dx + dy * dims[0] + dz * dims[0] * dims[1];
                if ((x + dx >= 0 && x + dx < dims[0]) &&
                    (y + dy >= 0 && y + dy < dims[1]) &&
                    (z + dz >= 0 && z + dz < dims[2]))
                {
                  cond = cond || ((new_index >= lwb) && (new_index < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dx + dy * dims[0] + dz * dims[0] * dims[1], x + dx, y + dy, z + dz, lwb));
                }
            }
        }

        // Check 3D diagonals for connectivity == 26
        if (connectivity == 26)
        {
            for (int i = 10; i < 26; ++i)
            {
                long dx = dx_val[i];
                long dy = dy_val[i];
                long dz = dz_val[i];
                ulong new_index = *index + dx + dy * dims[0] + dz * dims[0] * dims[1];
                if ((x + dx >= 0 && x + dx < dims[0]) &&
                    (y + dy >= 0 && y + dy < dims[1]) &&
                    (z + dz >= 0 && z + dz < dims[2]))
                {
                    cond = cond || ((new_index >= lwb) && (new_index < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long) dx + dy * dims[0] + dz * dims[0] * dims[1], x + dx, y + dy, z + dz, lwb));
                }
            }
        }
    }
}
*/
void remaining(BitArray *visited, PrioQueue *q, ulong *ranks, ulong *index, ulong *rank, ulong *dims, ulong lwb, ulong upb, int connectivity)
{
  ulong z = *index / (dims[0] * dims[1]);
  ulong x = (*index % (dims[0] * dims[1])) % dims[0];
  ulong y = (*index % (dims[0] * dims[1])) / dims[0];

  bool cond = 1;

  while (cond)
  {
    cond =
        ((x > 0) && (*index - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1, x - 1, y, z, lwb)) ||
        ((x < dims[0] - 1) && (*index + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)1, x + 1, y, z, lwb)) ||
        ((y > 0) && (*index - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, -dims[0], x, y - 1, z, lwb)) ||
        ((y < dims[1] - 1) && (*index + dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims[0], x, y + 1, z, lwb)) ||
        ((z > 0) && (*index >= lwb + dims[0] * dims[1]) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, -dims[1] * dims[0], x, y, z - 1, lwb)) ||
        ((z < dims[2] - 1) && (*index + dims[0] * dims[1] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, +dims[1] * dims[0], x, y, z + 1, lwb));

    if (connectivity >= 8)
    {
      cond = cond ||
             ((x > 0) && (y > 0) && (*index - 1 - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0], x - 1, y - 1, z, lwb)) ||
             ((x < dims[0] - 1) && (y > 0) && (*index + 1 - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)1 - dims[0], x + 1, y - 1, z, lwb)) ||
             ((x > 0) && (y < dims[1] - 1) && (*index - 1 + dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0], x - 1, y + 1, z, lwb)) ||
             ((y < dims[1] - 1) && (x < dims[0] - 1) && (*index + dims[0] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dims[0] + 1, x + 1, y + 1, z, lwb));
    }
    if (connectivity == 26)
    {
      cond = cond || (z > 0 && (((y > 0) && (x > 0) && (*index - dims[0] * dims[1] - dims[0] - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0] - dims[0] * dims[1], x - 1, y - 1, z - 1, lwb)) || ((y > 0) && (*index - dims[0] * dims[1] - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-dims[0] - dims[0] * dims[1], x, y - 1, z - 1, lwb)) || ((y > 0) && (x < dims[0] - 1) && (*index - dims[0] * dims[1] - dims[0] + 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 - dims[0] - dims[0] * dims[1], x + 1, y - 1, z - 1, lwb)) || ((x > 0) && (*index - dims[0] * dims[1] - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0] * dims[1], x - 1, y, z - 1, lwb)) || ((x < dims[0] - 1) && (*index - dims[0] * dims[1] + 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 - dims[0] * dims[1], x + 1, y, z - 1, lwb)) || ((y < dims[1] - 1) && (x > 0) && (*index - dims[0] * dims[1] + dims[0] - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0] - dims[0] * dims[1], x - 1, y + 1, z - 1, lwb)) || ((y < dims[1] - 1) && (*index - dims[0] * dims[1] + dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dims[0] - dims[0] * dims[1], x, y + 1, z - 1, lwb)) || ((y < dims[1] - 1) && (x < dims[0] - 1) && (*index - dims[0] * dims[1] + dims[0] + 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 + dims[0] - dims[0] * dims[1], x + 1, y + 1, z - 1, lwb)))) || (z < dims[2] - 1 && (((y > 0) && (x > 0) && (*index + dims[0] * dims[1] - dims[0] - 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0] + dims[0] * dims[1], x - 1, y - 1, z + 1, lwb)) || ((y > 0) && (*index + dims[0] * dims[1] - dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-dims[0] + dims[0] * dims[1], x, y - 1, z + 1, lwb)) || ((y > 0) && (x < dims[0] - 1) && (*index + dims[0] * dims[1] - dims[0] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 - dims[0] + dims[0] * dims[1], x + 1, y - 1, z + 1, lwb)) || ((x > 0) && (*index + dims[0] * dims[1] - 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0] * dims[1], x - 1, y, z + 1, lwb)) || ((x < dims[0] - 1) && (*index + dims[0] * dims[1] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 + dims[0] * dims[1], x + 1, y, z + 1, lwb)) || ((y < dims[1] - 1) && (x > 0) && (*index + dims[0] * dims[1] + dims[0] - 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0] + dims[0] * dims[1], x - 1, y + 1, z + 1, lwb)) || ((y < dims[1] - 1) && (*index + dims[0] * dims[1] + dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dims[0] + dims[0] * dims[1], x, y + 1, z + 1, lwb)) || ((y < dims[1] - 1) && (x < dims[0] - 1) && (*index + dims[0] * dims[1] + dims[0] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 + dims[0] + dims[0] * dims[1], x + 1, y + 1, z + 1, lwb))));
    }
  }
}


idx get_levelroot(Node *tree, idx x)
{
  /* index based, check whether this node is bottom */
  idx r = x;
  if (r == BOTTOM)
    return BOTTOM;
  value gv = tree->gval[x];
  while ((tree->parent[r] != BOTTOM) && (gv == tree->gval[tree->parent[r]]))
  {
    r = tree->parent[r];
  }
  /* tree compression */
  while (x != r)
  {
    idx y = tree->parent[x];
    tree->parent[x] = r;
    x = y;
  }
  return r;
} /* get_levelroot */

idx levelroot(Node *tree, idx index)
{
  return get_levelroot(tree, index);
} /* levelroot */

idx get_parent(Node *tree, idx x)
{
  return get_levelroot(tree, tree->parent[x]);
} /* get_parent */


/* ***************OLD FUNCTIONS NOT USED ANYMORE**************

void local_hist_rs(value *gvals, ulong lwb, ulong upb, ulong *hist, ulong *sorted, ulong step)
{
  ushort radix;
  memset(hist, 0, sizeof(ulong) * NUMBUCKETS);

  if (step == 0)
  {
    for (ulong i = lwb; i < upb; ++i)
    {
      radix = *((ushort *)&(gvals[i]));
      hist[radix]++;
    }
  }
  else
  {
    for (ulong i = lwb; i < upb; ++i)
    {
      radix = *((ushort *)(&(gvals[sorted[i]])) + step);
      hist[radix]++;
    }
  }
}

void create_sorted_ars(value *gvals, ulong lwb, ulong upb, ulong *hist, ulong *sorted_new, ulong *sorted_old, ulong step)
{
  ushort radix;

  if (step == 0)
  {
    for (ulong i = lwb; i < upb; ++i)
    {
      radix = *((ushort *)&(gvals[i]));
      sorted_new[hist[radix]++] = i; // sort in ascending order
    }
  }
  else
  {
    for (ulong i = lwb; i < upb; ++i)
    {
      radix = *((ushort *)(&(gvals[sorted_old[i]])) + step);
      sorted_new[hist[radix]++] = sorted_old[i]; // sort in ascending order
    }
  }
}

ulong *sort_image_pixels(Arguments *args, Node *tree)
{

  int nthreads = args->threads;
  int numstep = (int)abs(ceil((double)args->bpp / (16)));
  ulong *sortedrs[2];
  ulong *histogramrs[2];
  ulong *histogram;
  histogramrs[0] = (ulong *)malloc(NUMBUCKETS * sizeof(ulong));
  histogramrs[1] = (ulong *)malloc(NUMBUCKETS * sizeof(ulong));
  sortedrs[0] = (ulong *)malloc(tree->size_curr * sizeof(ulong));
  sortedrs[1] = (ulong *)malloc(tree->size_curr * sizeof(ulong));

  ulong *loc_hist = calloc(nthreads * NUMBUCKETS, sizeof(ulong));
  check_alloc(loc_hist, 610);

#pragma omp parallel num_threads(nthreads)
  {
    int id = omp_get_thread_num();                       // Thread number 
    ulong lwb = (id * tree->size_curr) / nthreads;       // Lower bound for current thread 
    ulong upb = ((id + 1) * tree->size_curr) / nthreads; // Upper bound for current thread 
    for (int step = 0; step < numstep; step++)
    {

      local_hist_rs(tree->gval, lwb, upb, &loc_hist[NUMBUCKETS * id], sortedrs[step % 2], step);
#pragma omp barrier

      if (id == 0)
      {
        histogram = histogramrs[step % 2];
        memset(histogram, 0, sizeof(ulong) * NUMBUCKETS);
        ulong prevhist, total = 0;
        for (int j = 0; j < NUMBUCKETS; j++)
        {
          for (int i = 0; i < nthreads; i++)
          {
            histogram[j] += loc_hist[NUMBUCKETS * i + j];
          }
          prevhist = histogram[j];
          histogram[j] = total;
          total += prevhist;
        }

        for (int j = 0; j < NUMBUCKETS; j++)
        {
          ulong sum = loc_hist[j], curr;
          for (int i = 1; i < nthreads; i++)
          {
            curr = loc_hist[NUMBUCKETS * i + j];
            loc_hist[NUMBUCKETS * i + j] = sum + histogram[j];
            sum += curr;
          }
          loc_hist[j] = histogram[j];
        }
      }
#pragma omp barrier
      create_sorted_ars(tree->gval, lwb, upb, &loc_hist[NUMBUCKETS * id], sortedrs[((step + 1) % 2)], sortedrs[step % 2], step);
    }
  }

  free(sortedrs[(numstep + 1) % 2]);
  free(histogramrs[(numstep) % 2]);
  free(loc_hist);
  return sortedrs[numstep % 2];
  // return ranks_inv;
}
*/


/* int get_neighbors(ulong *dims, ulong lwb, ulong upb, ulong *neighbors, ulong p, ulong x, ulong y, ulong z, int connectivity)
{
  int n = 0;

  if ((x < dims[0] - 1) && (p + 1 < upb))
    neighbors[n++] = p + 1;
  if ((y > 0) && (p - dims[0] >= lwb))
    neighbors[n++] = p - dims[0];
  if ((x > 0) && (p - 1 >= lwb))
    neighbors[n++] = p - 1;
  if ((y < dims[1] - 1) && (p + dims[0] < upb))
    neighbors[n++] = p + dims[0];

  if (connectivity == 8 || connectivity == 26)
  {
    if ((x < dims[0] - 1) && (y > 0) && (p + 1 - dims[0] >= lwb))
      neighbors[n++] = p + 1 - dims[0];
    if ((y > 0) && (x > 0) && (p - dims[0] - 1 >= lwb))
      neighbors[n++] = p - dims[0] - 1;
    if ((x > 0) && (y < dims[1] - 1) && (p - 1 + dims[0] < upb))
      neighbors[n++] = p - 1 + dims[0];
    if ((y < dims[1] - 1) && (x < dims[0] - 1) && (p + dims[0] + 1 < upb))
      neighbors[n++] = p + dims[0] + 1;
  }

  if (dims[2] > 1 && (connectivity == 6 || connectivity == 26))
  {
    if ((z > 0) && (p >= lwb + dims[0] * dims[1]))
      neighbors[n++] = p - dims[0] * dims[1];
    if ((z < dims[2] - 1) && (p + dims[0] * dims[1] < upb))
      neighbors[n++] = p + dims[0] * dims[1];
    if (connectivity == 26)
    {

      if ((z > 0) && (y > 0) && (x > 0) && (p - dims[0] * dims[1] - dims[0] - 1 >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] - dims[0] - 1;
      if ((z > 0) && (y > 0) && (p - dims[0] * dims[1] - dims[0] >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] - dims[0];
      if ((z > 0) && (y > 0) && (x < dims[0] - 1) && (p - dims[0] * dims[1] - dims[0] + 1 >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] - dims[0] + 1;
      if ((z > 0) && (x > 0) && (p - dims[0] * dims[1] - 1 >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] - 1;
      if ((z > 0) && (x < dims[0] - 1) && (p - dims[0] * dims[1] + 1 >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] + 1;
      if ((z > 0) && (y < dims[1] - 1) && (x > 0) && (p - dims[0] * dims[1] + dims[0] - 1 >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] + dims[0] - 1;
      if ((z > 0) && (y < dims[1] - 1) && (p - dims[0] * dims[1] + dims[0] >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] + dims[0];
      if ((z > 0) && (y < dims[1] - 1) && (x < dims[0] - 1) && (p - dims[0] * dims[1] + dims[0] + 1 >= lwb))
        neighbors[n++] = p - dims[0] * dims[1] + dims[0] + 1;

      if ((z < dims[2] - 1) && (y > 0) && (x > 0) && (p + dims[0] * dims[1] - dims[0] - 1 < upb))
        neighbors[n++] = p + dims[0] * dims[1] - dims[0] - 1;
      if ((z < dims[2] - 1) && (y > 0) && (p + dims[0] * dims[1] - dims[0] < upb))
        neighbors[n++] = p + dims[0] * dims[1] - dims[0];
      if ((z < dims[2] - 1) && (y > 0) && (x < dims[0] - 1) && (p + dims[0] * dims[1] - dims[0] + 1 < upb))
        neighbors[n++] = p + dims[0] * dims[1] - dims[0] + 1;
      if ((z < dims[2] - 1) && (x > 0) && (p + dims[0] * dims[1] - 1 < upb))
        neighbors[n++] = p + dims[0] * dims[1] - 1;
      if ((z < dims[2] - 1) && (x < dims[0] - 1) && (p + dims[0] * dims[1] + 1 < upb))
        neighbors[n++] = p + dims[0] * dims[1] + 1;
      if ((z < dims[2] - 1) && (y < dims[1] - 1) && (x > 0) && (p + dims[0] * dims[1] + dims[0] - 1 < upb))
        neighbors[n++] = p + dims[0] * dims[1] + dims[0] - 1;
      if ((z < dims[2] - 1) && (y < dims[1] - 1) && (p + dims[0] * dims[1] + dims[0] < upb))
        neighbors[n++] = p + dims[0] * dims[1] + dims[0];
      if ((z < dims[2] - 1) && (y < dims[1] - 1) && (x < dims[0] - 1) && (p + dims[0] * dims[1] + dims[0] + 1 < upb))
        neighbors[n++] = p + dims[0] * dims[1] + dims[0] + 1;
    }
  }

  return (n);
} /* get_neighbors */

/*
int get_neighbors(ulong *dims, ulong lwb, ulong upb, ulong *neighbors, ulong p, ulong x, ulong y, ulong z, int connectivity)
{
  error("connectivity", connectivity);
    int n = 0;
    ulong dx[26] = { 1,  0, -1, 0, 0,  0, 1,  1,  -1, -1,  0,  0,  0,  0, 1,  1,  1,  1, -1, -1, -1, -1, -1, 1, -1,  1}; // x deltas for 4 or 8 connectivity
    ulong dy[26] = { 0, -1,  0, 1, 0,  0, 1, -1,   1, -1,  1,  1, -1, -1, 1,  1, -1, -1,  1,  1, -1, -1,  0, 0,  0,  0}; // y deltas for 4 or 8 connectivity
    ulong dz[26] = { 0,  0,  0, 0, 1, -1, 0,  0,   0,  0,  1, -1,  1, -1, 1, -1,  1, -1,  1, -1,  1, -1,  1, 1, -1, -1}; // z deltas for 4 or 8 connectivity

    
    int base_delta = dims[0] * dims[1];
    int num_directions = connectivity == 8? 10: connectivity;
    
    // Check connectivity and add neighbors
    for (int i = 0; i < connectivity; i++) {
        int new_x = x + dx[i];
        int new_y = y + dy[i];
        int new_z = z + dz[i];
        ulong new_p = new_z * base_delta + new_y * dims[0] + new_x;
        
        // Check boundaries
        if (new_x < dims[0] && new_x >= 0 &&
            new_y < dims[1] && new_y >= 0 &&
            (dims[2] == 1 || (new_z < dims[2] && new_z >= 0)) &&
            new_p >= lwb && new_p < upb)
        {
            neighbors[n++] = new_p;
        }
    }

    return n;
}*/

/*
void remaining(BitArray *visited, PrioQueue *q, ulong *ranks, ulong *index, ulong *rank, ulong *dims, ulong lwb, ulong upb, int connectivity)
{
  ulong z = *index / (dims[0] * dims[1]);
  ulong x = (*index % (dims[0] * dims[1])) % dims[0];
  ulong y = (*index % (dims[0] * dims[1])) / dims[0];

  bool cond = 1;

  while (cond)
  {
    cond =
        ((x > 0) && (*index - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1, x - 1, y, z, lwb)) ||
        ((x < dims[0] - 1) && (*index + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)1, x + 1, y, z, lwb)) ||
        ((y > 0) && (*index - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, -dims[0], x, y - 1, z, lwb)) ||
        ((y < dims[1] - 1) && (*index + dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims[0], x, y + 1, z, lwb)) ||
        ((z > 0) && (*index >= lwb + dims[0] * dims[1]) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, -dims[1] * dims[0], x, y, z - 1, lwb)) ||
        ((z < dims[2] - 1) && (*index + dims[0] * dims[1] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, +dims[1] * dims[0], x, y, z + 1, lwb));

    if (connectivity >= 8)
    {
      cond = cond ||
             ((x > 0) && (y > 0) && (*index - 1 - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0], x - 1, y - 1, z, lwb)) ||
             ((x < dims[0] - 1) && (y > 0) && (*index + 1 - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)1 - dims[0], x + 1, y - 1, z, lwb)) ||
             ((x > 0) && (y < dims[1] - 1) && (*index - 1 + dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0], x - 1, y + 1, z, lwb)) ||
             ((y < dims[1] - 1) && (x < dims[0] - 1) && (*index + dims[0] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dims[0] + 1, x + 1, y + 1, z, lwb));
    }
    if (connectivity == 26)
    {
      cond = cond || (z > 0 && (((y > 0) && (x > 0) && (*index - dims[0] * dims[1] - dims[0] - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0] - dims[0] * dims[1], x - 1, y - 1, z - 1, lwb)) || ((y > 0) && (*index - dims[0] * dims[1] - dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-dims[0] - dims[0] * dims[1], x, y - 1, z - 1, lwb)) || ((y > 0) && (x < dims[0] - 1) && (*index - dims[0] * dims[1] - dims[0] + 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 - dims[0] - dims[0] * dims[1], x + 1, y - 1, z - 1, lwb)) || ((x > 0) && (*index - dims[0] * dims[1] - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0] * dims[1], x - 1, y, z - 1, lwb)) || ((x < dims[0] - 1) && (*index - dims[0] * dims[1] + 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 - dims[0] * dims[1], x + 1, y, z - 1, lwb)) || ((y < dims[1] - 1) && (x > 0) && (*index - dims[0] * dims[1] + dims[0] - 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0] - dims[0] * dims[1], x - 1, y + 1, z - 1, lwb)) || ((y < dims[1] - 1) && (*index - dims[0] * dims[1] + dims[0] >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dims[0] - dims[0] * dims[1], x, y + 1, z - 1, lwb)) || ((y < dims[1] - 1) && (x < dims[0] - 1) && (*index - dims[0] * dims[1] + dims[0] + 1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 + dims[0] - dims[0] * dims[1], x + 1, y + 1, z - 1, lwb)))) || (z < dims[2] - 1 && (((y > 0) && (x > 0) && (*index + dims[0] * dims[1] - dims[0] - 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 - dims[0] + dims[0] * dims[1], x - 1, y - 1, z + 1, lwb)) || ((y > 0) && (*index + dims[0] * dims[1] - dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-dims[0] + dims[0] * dims[1], x, y - 1, z + 1, lwb)) || ((y > 0) && (x < dims[0] - 1) && (*index + dims[0] * dims[1] - dims[0] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 - dims[0] + dims[0] * dims[1], x + 1, y - 1, z + 1, lwb)) || ((x > 0) && (*index + dims[0] * dims[1] - 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0] * dims[1], x - 1, y, z + 1, lwb)) || ((x < dims[0] - 1) && (*index + dims[0] * dims[1] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 + dims[0] * dims[1], x + 1, y, z + 1, lwb)) || ((y < dims[1] - 1) && (x > 0) && (*index + dims[0] * dims[1] + dims[0] - 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)-1 + dims[0] + dims[0] * dims[1], x - 1, y + 1, z + 1, lwb)) || ((y < dims[1] - 1) && (*index + dims[0] * dims[1] + dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)dims[0] + dims[0] * dims[1], x, y + 1, z + 1, lwb)) || ((y < dims[1] - 1) && (x < dims[0] - 1) && (*index + dims[0] * dims[1] + dims[0] + 1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, (long)+1 + dims[0] + dims[0] * dims[1], x + 1, y + 1, z + 1, lwb))));
    }
  }
}*/