#include "types.h"
#include "communication.h"
#include "attributes.h"
#include "flood.h"
#include "boundary.h"


BorderIndex idx_i(Boundary *b, ulong c, Direction d);
BorderIndex idx_j(Boundary *b, ulong c, Direction d);
void free_boundary(Boundary *b);
BorderIndex b_levelroot(BorderIndex bi);
BorderIndex b_id_levelroot(BorderIndex bi);
bool is_bottom(BorderIndex bi);
bool bi_equal(BorderIndex ai, BorderIndex bi);
BorderIndex b_parent_lr(BorderIndex bi);
BorderIndex b_parent(BorderIndex bi);
value b_gval(BorderIndex bi);
BoundaryNode b_node(BorderIndex bi);
void reset_border_idx(Boundary *b);


/* Merge functions */

void b_add(Boundary *c, Boundary *b, ulong s, ulong i, int merge_case)
{
  // Adding nodes from b in combined tree c 
  ulong idx_b = b->merge_idx[i];
  BorderIndex toadd = (BorderIndex){.b = b, .i = idx_b};
  toadd = b_levelroot(toadd);
  idx_b = toadd.i;
  b = toadd.b;

  if (b->array[idx_b].border_idx != BOTTOM)
    c->merge_idx[s] = b->array[idx_b].border_idx;
  else
  {
    ulong idx_c = c->size_curr++;
    b->array[idx_b].border_idx = idx_c;
    c->merge_idx[s] = idx_c;
    c->array[idx_c] = b->array[idx_b];
    c->border_ori[idx_c] = (BorderIndex){.b = b, .i = idx_b};
    if (merge_case)
    {
      void *b_attr = b->attribute + idx_b * b->size_attr;
      void *c_attr = c->attribute + idx_c * c->size_attr;
      clone_aux_data(NULL, &c_attr, b_attr);
    }
  }
}

Boundary *combine(Boundary *a, Boundary *b, Direction d, int merge_case)
{
  // Given two boundary trees a and b, return the combined tree c 

  ulong size_upb = a->size_curr + b->size_curr;
  Boundary *c = calloc(1, sizeof(Boundary));
  check_alloc(c, 203);
  c->array = calloc(size_upb, sizeof(BoundaryNode));
  check_alloc(c->array, 204);
  c->border_par = calloc(size_upb, sizeof(BorderIndex));
  check_alloc(c->border_par, 206);
  c->border_ori = calloc(size_upb, sizeof(BorderIndex));
  check_alloc(c->border_ori, 207);

  if(merge_case){
    c->size_attr = b->size_attr;
    c->attribute = malloc(size_upb * b->size_attr);
    check_alloc(c->attribute, 207);
  }

  reset_border_idx(a);
  reset_border_idx(b);

  ulong s = 0;   // points to entry point in c 
  ulong i, j, k; // index local to either a or b

  if (d == HORIZONTAL)
  {

    ulong size_merge_idx = a->offset[6] + b->offset[6] - 2 * (a->offset[2] - a->offset[1]);
    c->merge_idx = malloc(size_merge_idx * sizeof(ulong));
    check_alloc(c->merge_idx, 209);
    c->offset[0] = 0;

    // Top border 
    if (a->offset[1] != a->offset[0])
    {
      for (j = 0; j < a->dims[2]; j++)
      {
        for (k = 0; k < a->dims[0] + b->dims[0]; k++, s++)
        {
          if (k < a->dims[0])
          {
            i = j * a->dims[0] + k;
            b_add(c, a, s, a->offset[0] + i, merge_case);
          }
          else
          {
            i = j * b->dims[0] + k - a->dims[0];
            b_add(c, b, s, b->offset[0] + i, merge_case);
          }
        }
      }
    }
    c->offset[1] = s;

    // Right border 
    for (i = b->offset[1]; i < b->offset[2]; i++, s++)
      b_add(c, b, s, i, merge_case);
    c->offset[2] = s;

    // Bottom border 
    if (a->offset[3] != a->offset[2])
    {
      for (j = 0; j < a->dims[2]; j++)
      {
        for (k = 0; k < a->dims[0] + b->dims[0]; k++, s++)
        {
          if (k < a->dims[0])
          {
            i = j * a->dims[0] + k;
            b_add(c, a, s, a->offset[2] + i, merge_case);
          }
          else
          {
            i = j * b->dims[0] + k - a->dims[0];
            b_add(c, b, s, b->offset[2] + i, merge_case);
          }
        }
      }
    }
    c->offset[3] = s;

    // Left border 
    for (i = a->offset[3]; i < a->offset[4]; i++, s++)
      b_add(c, a, s, i, merge_case);
    c->offset[4] = s;

    // Front border 
    if (a->offset[5] != a->offset[4])
    {
      for (j = 0; j < a->dims[1]; j++)
      {
        for (k = 0; k < a->dims[0] + b->dims[0]; k++, s++)
        {
          if (k < a->dims[0])
          {
            i = j * a->dims[0] + k;
            b_add(c, a, s, a->offset[4] + i, merge_case);
          }
          else
          {
            i = j * b->dims[0] + k - a->dims[0];
            b_add(c, b, s, b->offset[4] + i, merge_case);
          }
        }
      }
    }
    c->offset[5] = s;

    // Back border 
    if (a->offset[6] != a->offset[5])
    {
      for (j = 0; j < a->dims[1]; j++)
      {
        for (k = 0; k < a->dims[0] + b->dims[0]; k++, s++)
        {
          if (k < a->dims[0])
          {
            i = j * a->dims[0] + k;
            b_add(c, a, s, a->offset[5] + i, merge_case);
          }
          else
          {
            i = j * b->dims[0] + k - a->dims[0];
            b_add(c, b, s, b->offset[5] + i, merge_case);
          }
        }
      }
    }

    c->offset[6] = s;
    c->dims[0] = a->dims[0] + b->dims[0];
    c->dims[1] = a->dims[1];
    c->dims[2] = a->dims[2];
  }
  else if (d == VERTICAL)
  {

    ulong size_merge_idx = a->offset[6] + b->offset[6] - 2 * (a->offset[1] - a->offset[0]);
    c->merge_idx = malloc(size_merge_idx * sizeof(ulong));
    check_alloc(c->merge_idx, 210);
    c->offset[0] = 0;

    // Top border 
    for (i = a->offset[0]; i < a->offset[1]; i++, s++)
      b_add(c, a, s, i, merge_case);
    c->offset[1] = s;

    // Right border 
    if (a->offset[2] != a->offset[1])
    {
      for (j = 0; j < a->dims[2]; j++)
      {
        for (k = 0; k < a->dims[1] + b->dims[1]; k++, s++)
        {
          if (k < a->dims[1])
          {
            i = j * a->dims[1] + k;
            b_add(c, a, s, a->offset[1] + i, merge_case);
          }
          else
          {
            i = j * b->dims[1] + k - a->dims[1];
            b_add(c, b, s, b->offset[1] + i, merge_case);
          }
        }
      }
    }
    c->offset[2] = s;

    // Botoom border 
    for (i = b->offset[2]; i < b->offset[3]; i++, s++)
      b_add(c, b, s, i, merge_case);

    c->offset[3] = s;

    // Left border 
    if (a->offset[4] != a->offset[3])
    {
      for (j = 0; j < a->dims[2]; j++)
      {
        for (k = 0; k < a->dims[1] + b->dims[1]; k++, s++)
        {
          if (k < a->dims[1])
          {
            i = j * a->dims[1] + k;
            b_add(c, a, s, a->offset[3] + i, merge_case);
          }
          else
          {
            i = j * b->dims[1] + k - a->dims[1];
            b_add(c, b, s, b->offset[3] + i, merge_case);
          }
        }
      }
    }
    c->offset[4] = s;

    // Front border 
    if (a->offset[5] != a->offset[4])
    {
      for (j = 0; j < a->dims[1] + b->dims[1]; j++)
      {
        for (k = 0; k < a->dims[0]; k++, s++)
        {
          if (j < a->dims[1])
          {
            i = j * a->dims[0] + k;
            b_add(c, a, s, a->offset[4] + i, merge_case);
          }
          else
          {
            i = (j - a->dims[1]) * (b->dims[0]) + k;
            b_add(c, b, s, b->offset[4] + i, merge_case);
          }
        }
      }
    }
    c->offset[5] = s;

    // Back border 
    if (a->offset[6] != a->offset[5])
    {
      for (j = 0; j < a->dims[1] + b->dims[1]; j++)
      {
        for (k = 0; k < a->dims[0]; k++, s++)
        {
          if (j < a->dims[1])
          {
            i = j * a->dims[0] + k;
            b_add(c, a, s, a->offset[5] + i, merge_case);
          }
          else
          {
            i = (j - a->dims[1]) * b->dims[0] + k;
            b_add(c, b, s, b->offset[5] + i, merge_case);
          }
        }
      }
    }

    c->offset[6] = s;
    c->dims[0] = a->dims[0];
    c->dims[1] = a->dims[1] + b->dims[1];
    c->dims[2] = a->dims[2];
  }
  else
  {

    ulong size_merge_idx = a->offset[6] + b->offset[6] - 2 * (a->offset[5] - a->offset[4]);
    c->merge_idx = malloc(size_merge_idx * sizeof(ulong));
    check_alloc(c->merge_idx, 211);
    c->offset[0] = 0;

    // Top border 
    for (i = a->offset[0]; i < a->offset[1]; i++, s++)
      b_add(c, a, s, i, merge_case);

    for (i = b->offset[0]; i < b->offset[1]; i++, s++)
      b_add(c, b, s, i, merge_case);
    c->offset[1] = s;

    // Right border 
    for (i = a->offset[1]; i < a->offset[2]; i++, s++)
      b_add(c, a, s, i, merge_case);

    for (i = b->offset[1]; i < b->offset[2]; i++, s++)
      b_add(c, b, s, i, merge_case);
    c->offset[2] = s;

    // Bottom border 
    for (i = a->offset[2]; i < a->offset[3]; i++, s++)
      b_add(c, a, s, i, merge_case);

    for (i = b->offset[2]; i < b->offset[3]; i++, s++)
      b_add(c, b, s, i, merge_case);
    c->offset[3] = s;

    // Left border 
    for (i = a->offset[3]; i < a->offset[4]; i++, s++)
      b_add(c, a, s, i, merge_case);

    for (i = b->offset[3]; i < b->offset[4]; i++, s++)
      b_add(c, b, s, i, merge_case);
    c->offset[4] = s;

    // Front border 
    for (i = a->offset[4]; i < a->offset[5]; i++, s++)
      b_add(c, a, s, i, merge_case);
    c->offset[5] = s;

    // Back border 
    for (i = b->offset[5]; i < b->offset[6]; i++, s++)
      b_add(c, b, s, i, merge_case);

    c->offset[6] = s;
    c->dims[0] = a->dims[0];
    c->dims[1] = a->dims[1];
    c->dims[2] = a->dims[2] + b->dims[2];
  }

  if (c->size_curr == 0) // No more merge steps to be done 
    return NULL;

  s = c->size_curr;

  BorderIndex origin, parent;
  ulong curr, ex = s;
  idx c_idx;

  for (i = 0; i < ex; i++)
  {
    curr = i;
    origin = c->border_ori[curr];
    parent = b_levelroot(b_parent(origin));

    while (true)
    {
      if (is_bottom(parent))
      {
        // case 1: parent is bottom 
        c->border_par[curr] = (BorderIndex){.b = c, .i = BOTTOM};
        break;
      }
      c_idx = b_node(parent).border_idx; // index in c, if not BOTTOM 
      if (c_idx != BOTTOM)
      {
        // case 2: parent is in c, avoid duplicates 
        c->border_par[curr] = (BorderIndex){.b = c, .i = c_idx};
        break;
      }
      else
      {
        // case 3: add parent to c 
        (parent.b)->array[parent.i].border_idx = s; /* to avoid duplicates in case 2 */
        c->border_ori[s] = parent;
        c->array[s] = b_node(parent);
        if(merge_case){
          void *b_attr = parent.b->attribute + parent.i * parent.b->size_attr;
          void *c_attr = c->attribute + s * c->size_attr;
          clone_aux_data(NULL, &c_attr, b_attr);
        }
        c->border_par[curr] = (BorderIndex){.b = c, .i = s};

        curr = s;
        s++;
      }
      parent = b_parent_lr(parent);
    }
  }

  c->size_curr = c->size_alloc = c->size_init = s;
  c->array = realloc(c->array, s * sizeof(BoundaryNode));
  check_alloc(c->array, 212);
  if(merge_case){
    c->attribute = realloc(c->attribute, s * c->size_attr);
    check_alloc(c->attribute, 213);
  }
  c->border_par = realloc(c->border_par, s * sizeof(BorderIndex));
  check_alloc(c->border_par, 214);
  c->border_ori = realloc(c->border_ori, s * sizeof(BorderIndex));
  check_alloc(c->border_ori, 215);
  if(merge_case != 1){
    c->border_lr = malloc(s * sizeof(BorderIndex));
    check_alloc(c->border_lr, 216);
    #pragma omp parallel for
    for (i = 0; i < s; ++i)
      c->border_lr[i] = (BorderIndex){.b = c, .i = BOTTOM};
  }


  return c;
}


void merge_branch(BorderIndex x, BorderIndex y, int merge_case)
{
  BorderIndex z, h;
  void *x_attr;
  void *y_attr;
  void *attr = NULL;
  void *attr1 = NULL;
  ulong size_attr = merge_case != 0? x.b->size_attr : 0;

  x = b_levelroot(x);
  y = b_levelroot(y);

  if (merge_case == ATTRIBUTES)
  {
    while (!is_bottom(x) && !bi_equal(x, y))
    {

      z = b_parent_lr(x);
      x_attr = x.b->attribute + x.i * size_attr;
      y_attr = y.b->attribute + y.i * size_attr;

      merge_aux_data(y_attr, x_attr);

      x.b->border_par[x.i] = y;
      x = z;
      y = b_parent_lr(y);
    }
  }


  else if (merge_case == PARENTS_AND_ATTRIBUTES)
  {
    while (!bi_equal(x, y) && !is_bottom(y))
    {

      x_attr = x.b->attribute + x.i * size_attr;
      y_attr = y.b->attribute + y.i * size_attr;

      z = b_parent_lr(x);

      if (!is_bottom(z) && (b_gval(z) >= b_gval(y)))
      {

        if (attr)
          merge_aux_data(x_attr, attr);

        x = z;
      }
      else
      {
        if (b_gval(x) == b_gval(y))
        {
          if (!is_bottom(y.b->border_lr[y.i]))
          {
            x.b->border_par[x.i] = x.b == y.b ? y : y.b->border_lr[y.i];
            x.b->border_lr[x.i] = x.b == y.b ? y.b->border_lr[y.i] : y;
          }
          else if (!is_bottom(x.b->border_lr[x.i]))
          {
            if (y_attr != NULL)
              clone_aux_data(NULL, &attr, y_attr);
            h = b_parent_lr(y);
            y.b->border_par[y.i] = x.b == y.b ? x : x.b->border_lr[x.i];
            y.b->border_lr[y.i] = x.b == y.b ? x.b->border_lr[x.i] : x;
            y = h;
            continue;
          }
          else
          {
            if (x.b != y.b)
            {
              y.b->border_lr[y.i] = x;
              x.b->border_lr[x.i] = y;
            }
            x.b->border_par[x.i] = y;
          }
        }
        else
          x.b->border_par[x.i] = y;

        attr ? merge_to_aux_data(NULL, &attr1, x_attr, attr) : clone_aux_data(NULL, &attr1, x_attr);
        clone_aux_data(NULL, &attr, x_attr);
        if (attr1)
          clone_aux_data(NULL, &x_attr, attr1);

        x = y;
        y = z;
      }
    }

    if (is_bottom(y) && attr)
    {
      while (!is_bottom(x))
      {
        x_attr = x.b->attribute + x.i * size_attr;
        merge_aux_data(x_attr, attr);
        x = b_parent_lr(x);
      }
    }

    delete_aux_data(attr);
    delete_aux_data(attr1);
  }

  else // ONLY PARENTS
  {
    while (!bi_equal(x, y) && !is_bottom(y))
    {

      z = b_parent_lr(x);

      if (!is_bottom(z) && (b_gval(z) >= b_gval(y)))
      {
        x = z;
      }
      else
      {
        if (b_gval(x) == b_gval(y))
        {
          if (!is_bottom(y.b->border_lr[y.i]))
          {
            x.b->border_par[x.i] = x.b == y.b ? y : y.b->border_lr[y.i];
            if (!is_bottom(x.b->border_lr[x.i]))
            {
              (x.b->border_lr[x.i].b)->border_par[x.b->border_lr[x.i].i] = y.b == x.b ? y.b->border_lr[y.i] : y;
              (x.b->border_lr[x.i].b)->border_lr[x.b->border_lr[x.i].i].i = BOTTOM;
            }
            x.b->border_lr[x.i].i = BOTTOM;
          }
          else if (!is_bottom(x.b->border_lr[x.i]))
          {
            h = b_parent_lr(y);
            y.b->border_par[y.i] = x.b == y.b ? x : x.b->border_lr[x.i];
            y.b->border_lr[y.i] = x.b == y.b ? x.b->border_lr[x.i] : x;
            y = h;
            continue;
          }
          else
          {
            if (x.b != y.b)
            {
              y.b->border_lr[y.i] = x;
              x.b->border_lr[x.i] = y;
            }
            x.b->border_par[x.i] = y;
          }
        }
        else
        {
          x.b->border_par[x.i] = y;
        }
        x = y;
        y = z;
      }
    }
  }
}

void merge_borders(Boundary *a, Boundary *b, Direction d, int merge_case) {
  // There are three cases:
  //  HORIZONTALLY: merge the right side of a with the left side of b 
  //  VERTICALLY  : merge the bottom of a with top of b, 
  //  DEPTH       : nerge the front border of b with the back border of a
  
  ulong length, offset;
  
  if (d == HORIZONTAL) {
    length = b->dims[1]*b->dims[2];
    offset = b->dims[1];
  } else if (d == VERTICAL) {
    length = b->dims[0]*b->dims[2];
    offset = b->dims[0];
  } else {
    length = b->dims[0]*b->dims[1];
    offset = b->dims[0];
  }

  BorderIndex 	x, y;
  // traverse border 
  value 	min_prev = 0;
  value 	min_curr = 1;
  bool		test_min = false;

  for (ulong c = 0; c < length; c++) {
    if(c %offset == 0) test_min = false;
    x = idx_i(a, c, d);  
    y = idx_j(b, c, d); 
    min_curr = b_gval(x);
    assert(b_gval(x) == b_gval(y));
    if (!test_min || min_curr > min_prev )
      merge_branch(x, y, merge_case);
    min_prev = min_curr;
    test_min = true;
  }
} /* merge */

void merge_boundary_trees(int grid[3], int grid_cur[3], int base[3], int col, int row, int slice, Boundary **bound_tree, int *n_merge, ulong store_item, int *merged, int merge_case)
{
  int myrank = rank();
  int neighrank;

  while (grid_cur[0] > 1 || grid_cur[1] > 1 || grid_cur[2] > 1)
  {
    if (grid_cur[0] > 1)
    {
      if (!(*merged) && col % (2 * base[0]) == 0)
      {
        neighrank = myrank + base[0];
        if (col < grid[0] - 1 && neighrank < np())
        {
          bound_tree[*n_merge + 1] = receive_boundary(neighrank, store_item, merge_case);
          merge_borders(bound_tree[*n_merge], bound_tree[*n_merge + 1], HORIZONTAL, merge_case);
          debug("Merging horizontally: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          bound_tree[*n_merge + 2] = combine(bound_tree[*n_merge], bound_tree[*n_merge + 1], HORIZONTAL, merge_case);
          *n_merge += 2;
        }
      }
      else if (!(*merged))
      {
        neighrank = myrank - base[0];
        send_boundary(bound_tree[*n_merge], neighrank, merge_case);
        (*merged)++;
      }
      else
      {
        (*merged)++;
      }
      base[0] *= 2;
      grid_cur[0] = (grid_cur[0] + 1) / 2;
    }

    if (grid_cur[1] > 1)
    {
      // Similar to horizontal logic for vertical merge
      if (!(*merged) && row % (2 * base[1]) == 0)
      {
        neighrank = myrank + base[1] * grid[0];
        if (row < grid[1] - 1 && neighrank < np())
        {
          bound_tree[*n_merge + 1] = receive_boundary(neighrank, store_item, merge_case);
          merge_borders(bound_tree[*n_merge], bound_tree[*n_merge + 1], VERTICAL, merge_case);
          debug("Merging vertically: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          bound_tree[*n_merge + 2] = combine(bound_tree[*n_merge], bound_tree[*n_merge + 1], VERTICAL, merge_case);
          debug("Combining vertically: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          *n_merge += 2;
        }
      }
      else if (!(*merged))
      {
        neighrank = myrank - base[1] * grid[0];
        send_boundary(bound_tree[*n_merge], neighrank, merge_case);
        (*merged)++;
      }
      else
      {
        (*merged)++;
      }

      base[1] *= 2;
      grid_cur[1] = (grid_cur[1] + 1) / 2;
    }

    if (grid_cur[2] > 1)
    {
      // Similar to horizontal logic for depth merge
      if (!(*merged) && slice % (2 * base[2]) == 0)
      {
        neighrank = myrank + base[2] * grid[0] * grid[1];
        if (slice < grid[2] - 1 && neighrank < np())
        {
          bound_tree[*n_merge + 1] = receive_boundary(neighrank, store_item, merge_case);
          merge_borders(bound_tree[*n_merge], bound_tree[*n_merge + 1], DEPTH, merge_case);
          debug("Merging in depth: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          bound_tree[*n_merge + 2] = combine(bound_tree[*n_merge], bound_tree[*n_merge + 1], DEPTH, merge_case);
          debug("Combining in depth: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          *n_merge += 2;
        }
      }
      else if (!(*merged))
      {
        neighrank = myrank - base[2] * grid[0] * grid[1];
        send_boundary(bound_tree[*n_merge], neighrank, merge_case);
        (*merged)++;
      }
      else
      {
        (*merged)++;
      }
      base[2] *= 2;
      grid_cur[2] = (grid_cur[2] + 1) / 2;
    }
  }
}

/* Update functions */

Boundary *realloc_b(Boundary *b, ulong size_tree_new, int border_case){

  // debug("REALLOC, oldsize %d size_tree_new %d", b->size_curr, size_tree_new);
 
  b->array         = realloc(b->array,          size_tree_new * sizeof(BoundaryNode));
  check_alloc(b->array, 227);

  b->border_par    = realloc(b->border_par,     size_tree_new * sizeof(BorderIndex));
  check_alloc(b->border_par, 229);

  b->border_lr     = realloc(b->border_lr,      size_tree_new * sizeof(BorderIndex));
  check_alloc(b->border_lr, 230);

  b->border_ori    = realloc(b->border_ori,     size_tree_new * sizeof(BorderIndex));
  check_alloc(b->border_ori, 231);

  if(border_case){
    b->attribute     = realloc(b->attribute,  size_tree_new * b->size_attr);
    check_alloc(b->attribute, 228);
  }

  b->reached       = realloc(b->reached,        size_tree_new * sizeof(bool));
  check_alloc(b->reached, 232);



  if (size_tree_new > b->size_curr) {
    #pragma omp parallel for 
    for (ulong i = b->size_curr; i < size_tree_new; i++)
      b->border_lr[i] = b->border_ori[i] = (BorderIndex) {.b = b, .i = BOTTOM};    
  }
 
  b->size_alloc = size_tree_new;
  return b;
}

void update_nested_relations(Boundary *c)
{
  if (c == NULL)
    return;
  for (size_t i = 0; i < c->size_curr; i++)
  {
    BorderIndex x = (BorderIndex){.b = c, .i = i};
    BorderIndex z = b_levelroot(x);

    if (bi_equal(x, z))
      continue;

    BorderIndex origin_a, origin_c, origin_lr;
    origin_a = c->border_ori[x.i];
    if (is_bottom(origin_a))
      continue;

    origin_c = c->border_ori[z.i];
    origin_lr = b_levelroot(origin_a);
    if (bi_equal(origin_lr, origin_c))
      continue;

    BorderIndex lr_a = origin_a.b->border_lr[origin_a.i], lr_c = origin_c.b->border_lr[origin_c.i];
    if (origin_a.b == origin_c.b)
    {
      origin_a.b->border_par[origin_a.i] = origin_c;
      if (!is_bottom(lr_a) && !is_bottom(lr_c))
      {
        lr_a.b->border_par[lr_a.i] = lr_c;
      }
      else if (!is_bottom(lr_a))
      {
        origin_c.b->border_lr[origin_c.i] = lr_a;
      }
      else if (!is_bottom(lr_c))
      {
        origin_a.b->border_lr[origin_a.i] = lr_c;
      }
    }
    else if (!is_bottom(lr_c) && origin_a.b == lr_c.b)
    {
      origin_a.b->border_par[origin_a.i] = lr_c;
      if (!is_bottom(lr_a))
      {
        lr_a.b->border_par[lr_a.i] = origin_c;
      }
    }
    else if (is_bottom(lr_c) && !is_bottom(lr_a))
    {
      lr_a.b->border_par[lr_a.i] = origin_c;
      origin_c.b->border_lr[origin_c.i] = origin_a;
      origin_a.b->border_lr[origin_a.i] = origin_c;
    }
    else if (is_bottom(lr_c) && is_bottom(lr_a))
    {
      origin_c.b->border_lr[origin_c.i] = origin_a;
      origin_a.b->border_lr[origin_a.i] = origin_c;
    }
  }
}

Boundary *adding_node_to_btree(BorderIndex x, BorderIndex s, int merge_case) {
  // If x does not have a parent with the same level of s, s is added in the boundary tree of x 
  if (x.b->size_curr == x.b->size_alloc)
    x.b = realloc_b(x.b, x.b->size_curr != 1 ? 1.5*x.b->size_curr: 100*x.b->size_curr, merge_case);

  if(merge_case){
    void *s_attr = s.b->attribute + s.i * s.b->size_attr;
    void *x_attr = x.b->attribute + x.b->size_curr * x.b->size_attr;
    clone_aux_data(NULL, &x_attr, s_attr);
  }
  x.b->array[x.b->size_curr] = s.b->array[s.i];
  x.b->reached[x.b->size_curr] =  false;
  x.b->border_par[x.i].b = x.b;
  x.b->border_par[x.i].i = x.b->size_curr;
  x.b->size_curr++;
  return x.b;
}


void update_node_attribute(BorderIndex x, BorderIndex s) {
  // Updating attribute of node x using the more recent node s //
  void *x_attr = x.b->attribute + x.i * x.b->size_attr;;
  void *s_attr = s.b->attribute + s.i * s.b->size_attr;
  clone_aux_data(NULL, &x_attr, s_attr);
}


Boundary *update_branch_with_attr(BorderIndex x, BorderIndex z, BorderIndex s)
{
  //Update the branch of the node x using its levelroot z in the merged tree, and the corresponding updated node s in the combined tree if it exists. 
  // If s does not exist (node has not been added in the combined tree) 

  while (is_bottom(s) && !(x.b->reached[x.i]))
  {
    if (!bi_equal(x, z))
    {
      update_node_attribute(x, z);
      if (x.b == z.b && b_gval(x) == b_gval(z))
      {
        x.b->border_par[x.i] = z;
        x.b->reached[x.i] = true;
        x = x.b->border_par[x.i];
        continue;
      }
      else if (!is_bottom(z.b->border_lr[z.i]) && z.b->border_lr[z.i].b == x.b && z.b->border_lr[z.i].i != x.i && b_gval(x) == b_gval(z))
      {
        x.b->border_par[x.i] = z.b->border_lr[z.i];

        x.b->reached[x.i] = true;
        x = x.b->border_par[x.i];
        continue;
      }
    }

    z = b_parent_lr(z);
  
    if (is_bottom(z))
    {
      x.b->border_par[x.i] = z;
      break;
    }
    if (x.b == z.b)
    {
      x.b->border_par[x.i] = z;
    }
    else
    {
      if (!is_bottom(z.b->border_lr[z.i]))
      {
        x.b->border_par[x.i] = z.b->border_lr[z.i];
      }
      else
      {
        x.b = adding_node_to_btree(x, z, PARENTS_AND_ATTRIBUTES);
        x.b->array[x.b->size_curr - 1].border_idx = z.b->array[z.i].border_idx;
        x.b->border_lr[x.b->border_par[x.i].i] = z;
        z.b->border_lr[z.i] = x.b->border_par[x.i];
      }
    }

    x.b->reached[x.i] = true;
    x = x.b->border_par[x.i];
    if (s.b != NULL){
      s.i = z.b->array[z.i].border_idx;
      if(is_bottom(s) && !is_bottom(z.b->border_lr[z.i]))
         s.i = z.b->border_lr[z.i].b->array[z.b->border_lr[z.i].i].border_idx;
    }
  }

  // If s exists (node has  been added in the combined tree and may have been updated in later steps) 
  if (s.b != NULL)
  {

    BorderIndex origin;
    s = b_levelroot(s);

    while (!(x.b->reached[x.i]) && !is_bottom(s))
    {
      update_node_attribute(x, s);
      origin = s.b->border_ori[s.i];

      if (!bi_equal(x, origin) && !is_bottom(origin))
      {
        if (x.b == origin.b && b_gval(x) == b_gval(origin))
        {
          x.b->border_par[x.i] = origin;
          x.b->reached[x.i] = true;
          x = x.b->border_par[x.i];
          continue;
        }
        // This might solve the Big Bang but I have honestly no idea
        else if (!is_bottom(origin.b->border_lr[origin.i]) && origin.b->border_lr[origin.i].b == x.b && origin.b->border_lr[origin.i].i != x.i && b_gval(x) == b_gval(origin))
        {
          x.b->border_par[x.i] = origin.b->border_lr[origin.i];
          x.b->reached[x.i] = true;
          x = x.b->border_par[x.i];
          continue;
        }
      }

      s = b_parent_lr(s);

      if (!is_bottom(s))
      {
        origin = s.b->border_ori[s.i];
        if (x.b == origin.b)
        {
          x.b->border_par[x.i] = origin;
        }
        else if (!is_bottom(origin) && !is_bottom(origin.b->border_lr[origin.i]))
        {
          x.b->border_par[x.i] = origin.b->border_lr[origin.i];
        }
        else
        {
          x.b = adding_node_to_btree(x, s, PARENTS_AND_ATTRIBUTES);
          x.b->array[x.b->size_curr - 1].border_idx = s.i;
          if (!is_bottom(origin))
          {
            x.b->border_lr[x.b->size_curr - 1] = origin;
            origin.b->border_lr[origin.i] = (BorderIndex){.b = x.b, .i = x.b->size_curr - 1};
          }
          else
            s.b->border_ori[s.i] = (BorderIndex){.b = x.b, .i = x.b->size_curr - 1};
        }
        x.b->reached[x.i] = true;
        x = x.b->border_par[x.i];
      }
      else
      {
        x.b->border_par[x.i].i = BOTTOM;
        x.b->reached[x.i] = true;
      }
    }
  }
  return x.b;
}

Boundary *update_branch_without_attr(BorderIndex x, BorderIndex z, BorderIndex s)
{
  //Update the branch of the node x using its levelroot z in the merged tree, and the corresponding updated node s in the combined tree if it exists. 
  // If s does not exist (node has not been added in the combined tree) 

  while (is_bottom(s) && !(x.b->reached[x.i]))
  {
    if (!bi_equal(x, z))
    {
      if (x.b == z.b && b_gval(x) == b_gval(z))
      {
        x.b->border_par[x.i] = z;
        x.b->reached[x.i] = true;
        x = x.b->border_par[x.i];
        continue;
      }
      else if (!is_bottom(z.b->border_lr[z.i]) && z.b->border_lr[z.i].b == x.b && z.b->border_lr[z.i].i != x.i && b_gval(x) == b_gval(z))
      {
        x.b->border_par[x.i] = z.b->border_lr[z.i];

        x.b->reached[x.i] = true;
        x = x.b->border_par[x.i];
        continue;
      }
    }

    z = b_parent_lr(z);
  
    if (is_bottom(z))
    {
      x.b->border_par[x.i] = z;
      break;
    }
    if (x.b == z.b)
    {
      x.b->border_par[x.i] = z;
    }
    else
    {
      if (!is_bottom(z.b->border_lr[z.i]))
      {
        x.b->border_par[x.i] = z.b->border_lr[z.i];
      }
      else
      {
        x.b = adding_node_to_btree(x, z, PARENTS);
        x.b->array[x.b->size_curr - 1].border_idx = z.b->array[z.i].border_idx;
        x.b->border_lr[x.b->border_par[x.i].i] = z;
        z.b->border_lr[z.i] = x.b->border_par[x.i];
      }
    }

    x.b->reached[x.i] = true;
    x = x.b->border_par[x.i];
    if (s.b != NULL){
      s.i = z.b->array[z.i].border_idx;
      if(is_bottom(s) && !is_bottom(z.b->border_lr[z.i]))
         s.i = z.b->border_lr[z.i].b->array[z.b->border_lr[z.i].i].border_idx;
    }
  }

  // If s exists (node has  been added in the combined tree and may have been updated in later steps) 
  if (s.b != NULL)
  {

    BorderIndex origin;
    s = b_levelroot(s);

    while (!(x.b->reached[x.i]) && !is_bottom(s))
    {
      origin = s.b->border_ori[s.i];

      if (!bi_equal(x, origin) && !is_bottom(origin))
      {
        if (x.b == origin.b && b_gval(x) == b_gval(origin))
        {
          x.b->border_par[x.i] = origin;
          x.b->reached[x.i] = true;
          x = x.b->border_par[x.i];
          continue;
        }
        // This might solve the Big Bang but I have honestly no idea
        else if (!is_bottom(origin.b->border_lr[origin.i]) && origin.b->border_lr[origin.i].b == x.b && origin.b->border_lr[origin.i].i != x.i && b_gval(x) == b_gval(origin))
        {
          x.b->border_par[x.i] = origin.b->border_lr[origin.i];
          x.b->reached[x.i] = true;
          x = x.b->border_par[x.i];
          continue;
        }
      }

      s = b_parent_lr(s);

      if (!is_bottom(s))
      {
        origin = s.b->border_ori[s.i];
        if (x.b == origin.b)
        {
          x.b->border_par[x.i] = origin;
        }
        else if (!is_bottom(origin) && !is_bottom(origin.b->border_lr[origin.i]))
        {
          x.b->border_par[x.i] = origin.b->border_lr[origin.i];
        }
        else
        {
          x.b = adding_node_to_btree(x, s, PARENTS);
          x.b->array[x.b->size_curr - 1].border_idx = s.i;
          if (!is_bottom(origin))
          {
            x.b->border_lr[x.b->size_curr - 1] = origin;
            origin.b->border_lr[origin.i] = (BorderIndex){.b = x.b, .i = x.b->size_curr - 1};
          }
          else
            s.b->border_ori[s.i] = (BorderIndex){.b = x.b, .i = x.b->size_curr - 1};
        }
        x.b->reached[x.i] = true;
        x = x.b->border_par[x.i];
      }
      else
      {
        x.b->border_par[x.i].i = BOTTOM;
        x.b->reached[x.i] = true;
      }
    }
  }
  return x.b;
}

Boundary *update_merged_boundary(Boundary *b, Boundary *c, int merge_case) {
  // Update of the boundary tree b using the combined, already updated, boundary tree c if it exists 
  
  BorderIndex x, z, s;
  b->reached = calloc(b->size_alloc, sizeof(bool)); 
  check_alloc(b->reached, 202);
  Boundary *(*update_branch_func)(BorderIndex, BorderIndex, BorderIndex);

  // Set the function pointer based on the merge_case
  if (merge_case) {
      update_branch_func = &update_branch_with_attr;  // Function with attribute update
  } else {
      update_branch_func = &update_branch_without_attr;    // Function without attribute update
  }
  for (size_t i = 0; i< b->size_curr; i++) {
    x = (BorderIndex) {.b = b, .i = i};
    x = b_id_levelroot(x);
    if (x.b->reached[x.i]) continue;

    z = b_levelroot(x);    
    s = (BorderIndex) {.b = c, .i = z.b->array[z.i].border_idx };  
    if(is_bottom(s) && !is_bottom(z.b->border_lr[z.i])){
      s.i = z.b->border_lr[z.i].b->array[z.b->border_lr[z.i].i].border_idx;
    }

    if (!bi_equal(x, z)){
      if(x.b == z.b){
         x.b->border_par[x.i] = z;
         x.b->reached[x.i] = true;
         x = x.b->border_par[x.i];
        if (x.b->reached[x.i]) continue;
      }
      else if (!is_bottom(z.b->border_lr[z.i]) && z.b->border_lr[z.i].b == x.b && z.b->border_lr[z.i].i != x.i){
          x.b->border_par[x.i] = z.b->border_lr[z.i];
          x.b->reached[x.i] = true;
          x = x.b->border_par[x.i];
           if (x.b->reached[x.i]) continue;
      }
   }
    
    b = update_branch_func(x,z,s);
  }    
  free(b->reached);
  b->reached = NULL;
  return b;
}

void update_attributes(Boundary *a, Boundary *b, Boundary *c) {
  /* Update of the boundary tree b using the combined, already updated, boundary tree c if it exists */
  
  BorderIndex x, z, s;
  a->reached = calloc(a->size_alloc, sizeof(bool)); check_alloc(a->reached, 202);
  b->reached = calloc(b->size_alloc, sizeof(bool)); check_alloc(b->reached, 202);

  for (size_t i = 0; i< a->size_curr; i++) {
    x = (BorderIndex) {.b = a, .i = i};
    if ((x.b->reached[x.i])) continue;

    z = b_levelroot(x);
    s = (BorderIndex) {.b = c, .i = z.b->array[z.i].border_idx };
    if (is_bottom(s))
      update_node_attribute(x,z);
    else
      update_node_attribute(x, s);
    x.b->reached[x.i] = true;
  }

  for (size_t i = 0; i < b->size_curr; i++)
  {
    x = (BorderIndex) {.b = b, .i = i};
    if ((x.b->reached[x.i])) continue;

    z = b_levelroot(x);
    s = (BorderIndex) {.b = c, .i = z.b->array[z.i].border_idx };
    
    if (is_bottom(s))
      update_node_attribute(x,z);
    else
      update_node_attribute(x,s);

    x.b->reached[x.i] = true;
  }

  free(a->reached);
  a->reached = NULL;
  free(b->reached);
  b->reached = NULL;
}

void update_boundary_trees(int grid[3], int grid_cur[3], int base[3], int col, int row, int slice, Boundary **bound_tree, int *n_merge, int *merged, ulong store_item, int merge_case)
{
  int myrank = rank();
  int neighrank;
  while (grid_cur[0] != 1 || grid_cur[1] != 1 || grid_cur[2] != 1)
  {
    if (grid_cur[2] >= grid_cur[1] && grid_cur[2] >= grid_cur[0] && grid_cur[2] > 1)
    {
      if (*merged <= 0 && slice % (2 * base[2]) == 0)
      {
        neighrank = myrank + base[2] * grid[0] * grid[1];
        if (slice < grid[2] - 1 && neighrank < np())
        {
          if(merge_case != 1){
            update_nested_relations(bound_tree[*n_merge]);
            update_merged_boundary(bound_tree[*n_merge - 2], bound_tree[*n_merge], merge_case);
            update_merged_boundary(bound_tree[*n_merge - 1], bound_tree[*n_merge],merge_case);
          } else {
            update_attributes(bound_tree[*n_merge-2],bound_tree[*n_merge-1], bound_tree[*n_merge]);
          }
          send_updated_boundary(bound_tree[*n_merge - 1], neighrank, merge_case);
          debug("Updating in depth: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          free_boundary(bound_tree[*n_merge]);
          free_boundary(bound_tree[*n_merge - 1]);
          *n_merge -= 2;
        }
      }
      else if (*merged <= 0)
      {
        neighrank = myrank - base[2] * grid[0] * grid[1];
        bound_tree[*n_merge] = receive_updated_boundary(bound_tree[*n_merge], neighrank, merge_case);
        (*merged)--;
      }
      else
      {
        (*merged)--;
      }
      base[2] /= 2;
      grid_cur[2] = grid_cur[2] % 2 ? (grid_cur[2] + 1) / 2 : grid_cur[2] / 2;
    }

    if (grid_cur[1] >= grid_cur[0] && grid_cur[1] > grid_cur[2] && grid_cur[1] > 1)
    {
      if (*merged <= 0 && row % (2 * base[1]) == 0)
      {
        neighrank = myrank + base[1] * grid[0];
        if (row < grid[1] - 1 && neighrank < np())
        {
          if(merge_case != 1){
            update_nested_relations(bound_tree[*n_merge]);
            update_merged_boundary(bound_tree[*n_merge - 2], bound_tree[*n_merge], merge_case);
            update_merged_boundary(bound_tree[*n_merge - 1], bound_tree[*n_merge], merge_case);
          } else {
     	      update_attributes(bound_tree[*n_merge-2], bound_tree[*n_merge-1], bound_tree[*n_merge]);       
          }
          send_updated_boundary(bound_tree[*n_merge - 1], neighrank, merge_case);
          debug("Updating vertically: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          free_boundary(bound_tree[*n_merge]);
          free_boundary(bound_tree[*n_merge - 1]);
          *n_merge -= 2;
        }
      }
      else if (*merged <= 0)
      {
        neighrank = myrank - base[1] * grid[0];
        bound_tree[*n_merge] = receive_updated_boundary(bound_tree[*n_merge], neighrank, merge_case);
        (*merged)--;
      }

      else
      {
        (*merged)--;
      }
      base[1] /= 2;
      grid_cur[1] = (grid_cur[1] + 1) / 2;
    }

    if (grid_cur[0] > grid_cur[1] && grid_cur[0] > grid_cur[2] && grid_cur[0] > 1)
    {
      if (*merged <= 0 && col % (2 * base[0]) == 0)
      {
        neighrank = myrank + base[0];
        if (col < grid[0] - 1 && neighrank < np())
        {
          if(merge_case != 1){
            update_nested_relations(bound_tree[*n_merge]);
            update_merged_boundary(bound_tree[*n_merge - 2], bound_tree[*n_merge], merge_case);
            update_merged_boundary(bound_tree[*n_merge - 1], bound_tree[*n_merge], merge_case);
          } else {
            update_attributes(bound_tree[*n_merge-2], bound_tree[*n_merge-1],bound_tree[*n_merge]);
          }
          send_updated_boundary(bound_tree[*n_merge - 1], neighrank, merge_case);
          debug("Updating horizontally: wallclock time = %0.2f",
                (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));
          free_boundary(bound_tree[*n_merge]);
          free_boundary(bound_tree[*n_merge - 1]);
          *n_merge -= 2;
        }
      }
      else if (*merged <= 0)
      {
        neighrank = myrank - base[0];
        bound_tree[*n_merge] = receive_updated_boundary(bound_tree[*n_merge], neighrank, merge_case);
        (*merged)--;
      }

      else
      {
        (*merged)--;
      }
      base[0] /= 2;
      grid_cur[0] = (grid_cur[0] + 1) / 2;
    }
  }
}

/* Local tree update functions */

void bound_to_tree(Node *local_tree, Boundary *b, ulong b_idx, ulong tree_idx){
  local_tree->gval[tree_idx]        = b->array[b_idx].gval;
  local_tree->parent[tree_idx]      = BOTTOM;
} /* bound_to_tree */

Node *correct_local_tree(Node *local_tree, Boundary *b, int merge_case)
{
  // Adding and updating the nodes from the local tree

  idx b_parent;
  ulong size_tree_init = local_tree->size_curr;
  ulong size_tree_new = (local_tree->size_curr + b->size_curr - b->size_init);

  local_tree->parent = realloc(local_tree->parent, size_tree_new * sizeof(idx));
  local_tree->gval = realloc(local_tree->gval, size_tree_new * sizeof(value));

  if (merge_case == PARENTS)
  {
    for (ulong i = 0; i < b->size_curr; i++)
    {
      if (i >= b->size_init)
      {
        b->array[i].index = local_tree->size_curr;
        bound_to_tree(local_tree, b, i, local_tree->size_curr);
        (local_tree->size_curr)++;
      }

      b_parent = b->border_par[i].i;

      if (b_parent == BOTTOM)
        local_tree->parent[b->array[i].index] = BOTTOM;
      else if ((ulong)b_parent < b->size_init)
        local_tree->parent[b->array[i].index] = b->array[b_parent].index;
      else
        local_tree->parent[b->array[i].index] = size_tree_init + b_parent - b->size_init;
    }
  }

  else if (merge_case == ATTRIBUTES)
  {
    void *b_attr;
    void *t_attr;
    for (ulong i = 0; i < b->size_curr; i++)
    {
      b_attr = (char *)b->attribute + i * b->size_attr;
      t_attr = local_tree->attribute + b->array[i].index * local_tree->size_attr;
      clone_aux_data(NULL, &t_attr, b_attr);
    }
  }
  else // PARENTS & ATTRIBUTES
  {
    void *b_attr;
    void *t_attr;
    local_tree->attribute = realloc(local_tree->attribute, size_tree_new * local_tree->size_attr);
    for (ulong i = 0; i < b->size_curr; i++)
    {
      if (i >= b->size_init)
      {
        b->array[i].index = local_tree->size_curr;
        bound_to_tree(local_tree, b, i, local_tree->size_curr);
        (local_tree->size_curr)++;
      }

      b_attr = b->attribute + i * b->size_attr;
      t_attr = local_tree->attribute + b->array[i].index * local_tree->size_attr;
      clone_aux_data(NULL, &t_attr, b_attr);

      b_parent = b->border_par[i].i;

      if (b_parent == BOTTOM)
        local_tree->parent[b->array[i].index] = BOTTOM;
      else if ((ulong)b_parent < b->size_init)
        local_tree->parent[b->array[i].index] = b->array[b_parent].index;
      else
        local_tree->parent[b->array[i].index] = size_tree_init + b_parent - b->size_init;
    }
  }
  return local_tree;
}


/* BOUNDARY CREATION FUNCTIONS */

int calculate_merges(int grid[3], int grid_cur[3], int base[3], int col, int row, int slice, int myrank, int nprocesses) {
  int n_merge = 0;

  while (grid_cur[0] > 1 || grid_cur[1] > 1 || grid_cur[2] > 1) {
      if (grid_cur[0] > 1) {
          if (col % (2 * base[0]) == 0 && col < grid[0] - 1 && myrank + base[0] < nprocesses) {
              n_merge++;
          }
          base[0] *= 2;
          grid_cur[0] = (grid_cur[0] + 1) / 2;
      }
      if (grid_cur[1] > 1) {
          if (row % (2 * base[1]) == 0 && row < grid[1] - 1 && myrank + base[1] * grid[0] < nprocesses) {
              n_merge++;
          }
          base[1] *= 2;
          grid_cur[1] = (grid_cur[1] + 1) / 2;
      }
      if (grid_cur[2] > 1) {
          if (slice % (2 * base[2]) == 0 && slice < grid[2] - 1 && myrank + base[2] * grid[0] * grid[1] < nprocesses) {
              n_merge++;
          }
          base[2] *= 2;
          grid_cur[2] = (grid_cur[2] + 1) / 2;
      }
  }

  return n_merge;
}

void node_to_bound(Node *local_tree, Boundary *b, ulong b_idx, ulong tree_idx, int border_case){
  b->array[b_idx].index         = tree_idx;
  b->array[b_idx].gval          = local_tree->gval[tree_idx];
  b->array[b_idx].border_idx    = local_tree->border_idx[tree_idx];
  if(border_case){
    void *b_attr = b->attribute + b_idx * b->size_attr;
    void *l_attr = local_tree->attribute + tree_idx * local_tree->size_attr;
    clone_aux_data(NULL, &b_attr, l_attr);
  }
} /* node_to_bound */

void add_side(Node *local_tree, Boundary *b, int side, ulong *length, ulong offset, ulong *offset_bd, long *increment, int connectivity, int border_case)
{
  // Adding nodes from local_tree to the boundary tree
  b->offset[side] = *offset_bd;

  ulong s = 0, k, l, indx;

  for (ulong j = 0; j < length[1]; j++)
  {
    for (ulong i = 0; i < length[0]; ++i, s++)
    {
      k = offset + (j * increment[1]) + (i * increment[0]);
      l = k + increment[2];

      if (connectivity >= 8 && local_tree->gval[l] > local_tree->gval[k])
      {
        k = l;
      }
      else if (connectivity < 8 && local_tree->gval[l] < local_tree->gval[k])
      {
        k = l;
      }

      k = get_levelroot(local_tree, k);

      if (local_tree->border_idx[k] != BOTTOM)
        b->merge_idx[*offset_bd + s] = local_tree->border_idx[k];
      else
      {
        indx = b->size_curr++;
        local_tree->border_idx[k] = indx;
        b->merge_idx[*offset_bd + s] = indx;
        node_to_bound(local_tree, b, indx, k, border_case);
      }
    }
  }
  *offset_bd += s;

} /* add_side */

Boundary *add_ancestors(Node *local_tree, Boundary *b, int border_case) {
  // We have to keep track of all nodes that are ancestors which are not in the border 
  // Note that in order to visit every node only once, we need to iterate over the border, and only add ancestors that are not in the border nor already in the map. Administration is done here. 
  
  //The theoretical upper limit of the boundary including ancestors is   (G-1) * (L/2)   where G is the number of gray levels, and L is the length of the border   We do not allocate for this upper limit, as the number of nodes added to a boundary tree was found to be much smaller in practice.   As the size of the ancestors map is dynamic and depending on the input image, we then have to make an assumption for a good initial size of the map. We here chose to make it 5 times bigger, and reallocate 1.5 its current size (1.5 is the most efficient factor)  when needed, which is O(n).
  
  ulong origsize  = b->size_curr;
  idx   parent, bx, bx_par;
  ulong curr;

  for (ulong i = 0; i < origsize; ++i)
  {

    curr = b->array[i].index; // index in tree

    while (true)
    {

      parent = get_parent(local_tree, curr);  

      if (parent == BOTTOM)
      {
        b->border_par[local_tree->border_idx[curr]] = (BorderIndex){.b = b, .i = BOTTOM};
        break; 
      }

      bx_par = local_tree->border_idx[parent]; 
      bx = local_tree->border_idx[curr];      

      if (bx_par != BOTTOM)
      {
        // parent is already in the tree
        b->border_par[bx] = (BorderIndex){.b = b, .i = bx_par};
        break;
      }

      else
      {
        if (b->size_curr == b->size_alloc)
        {
          b->size_alloc = b->size_alloc > 1 ? 1.5 * b->size_alloc : 200;
          b->array = realloc(b->array, b->size_alloc * sizeof(BoundaryNode));
          check_alloc(b->array, 224);

          b->border_par = realloc(b->border_par, b->size_alloc * sizeof(BorderIndex));
          check_alloc(b->border_par, 225);

          if(border_case){
            b->attribute = realloc(b->attribute, b->size_alloc * b->size_attr);
            check_alloc(b->attribute, 226);
          }
        }

        // add parent to border 
        local_tree->border_idx[parent] = b->size_curr;
        node_to_bound(local_tree, b, b->size_curr, parent, border_case);

        // set border_parents
        b->border_par[bx] = (BorderIndex){.b = b, .i = local_tree->border_idx[parent]};
        b->size_curr++; /* move to next empty index */
        curr = parent;  /* next! */
      }
    }
  }
  return b;
} /* add_ancestors */


Boundary *create_boundary(Node* local_tree,  int connectivity, int border_case){
  
  //  The boundary tree is an array that includes the levelroots of the nodes located in the tile border, and their parents.  In case of 4 or 6 connectivity, we merge along the overlapping nodes with the smallest intensity. With 8 or 26 connectivity, we need to merge along the overlapping nodes with th highest intensity, as the case is more complex. 

  Boundary *b = calloc(1, sizeof(Boundary));   check_alloc(b, 217);
  ulong     length[2]    = {0};
  long      increment[3] = {0};
  long      offset[3]    = {0};
  ulong     offset_bd 	 = 0;
  ulong	    size_border	 = 0;
  bool      *border 	   = data_properties.border;
  ulong *dims = data_properties.dims_process;
  local_tree->border_idx = malloc(local_tree->size_curr*sizeof(idx));
  memset(local_tree->border_idx, -1, local_tree->size_curr*sizeof(idx));

  b->dims[0] = dims[0];
  b->dims[1] = dims[1];
  b->dims[2] = dims[2];

  if(border[0]){
    b->dims[0]--;
    offset[1]++;
    offset[2]++;
  }
  if(border[1])
    b->dims[0]--;  
  if(border[2]){
    b->dims[1]--;
    offset[0] +=  dims[0];
    offset[2] +=  dims[0];
  }
  if(border[3])
    b->dims[1]--; 
  if(border[4]){
    b->dims[2]--;
    offset[0] +=  dims[0]*dims[1];
    offset[1] +=  dims[0]*dims[1];
  }
  if(border[5]) 
    b->dims[2]--;

  size_border       = (border[0]+border[1])*b->dims[2]*b->dims[1]
    + (border[2]+border[3]) * b->dims[2]*b->dims[0] + (border[4]+border[5]) * b->dims[0]*b->dims[1];
  b->size_alloc     = 2*size_border;

  b->array          = malloc(b->size_alloc  * sizeof(BoundaryNode)); check_alloc(b->array, 218);
  if(border_case){
    b->size_attr 	    = local_tree->size_attr;
    b->attribute      = malloc(b->size_alloc  * b->size_attr);
  }
  b->border_par     = malloc(b->size_alloc  * sizeof(BorderIndex));  check_alloc(b->border_par, 221);
  b->merge_idx      = malloc(size_border    * sizeof(BorderIndex));  check_alloc(b->merge_idx, 219);

  //     0: Adding the bottom of the tile    //
  
  if(border[2]){
    increment[0] = 1;
    increment[1] = dims[0] * dims[1];
    increment[2] = dims[0];
    length[0]    = b->dims[0];
    length[1]    = b->dims[2];
    add_side(local_tree, b, 0, length, 0+offset[1], &offset_bd, increment, connectivity, border_case);
  } else
    b->offset[0] = offset_bd;

  //    1: Adding the right side of the tile   //

  if(border[1]){
    increment[0] = dims[0];
    increment[1] = dims[0] * dims[1];
    increment[2] = -1;
    length[0]    = b->dims[1];
    length[1]    = b->dims[2];
    add_side(local_tree,    b, 1, length, dims[0]-1+offset[0], &offset_bd, increment, connectivity, border_case);
  } else
    b->offset[1] = offset_bd;
  
  //     2: Adding the top side of the tile  //

  if(border[3]){
    increment[0] = 1;
    increment[1] = dims[0] * dims[1];
    increment[2] = -dims[0];
    length[0]    = b->dims[0];
    length[1]    = b->dims[2];
    add_side(local_tree,   b, 2, length, dims[0] * (dims[1]-1)+offset[1], &offset_bd, increment, connectivity, border_case);
  } else
    b->offset[2] = offset_bd;

  //     2: Adding the left side of the tile  //
  
  if(border[0]){
    increment[0] = dims[0];
    increment[1] = dims[0] * dims[1];
    increment[2] = 1;
    length[0]    = b->dims[1];
    length[1]    = b->dims[2];
    add_side(local_tree,   b, 3, length,  0+offset[0], &offset_bd, increment, connectivity, border_case);
  } else
    b->offset[3] = offset_bd;

  //     4: Adding the front side   //
  
  if(border[4]){
    increment[0] = 1;
    increment[1] = dims[0];
    increment[2] = dims[0]*dims[1];

    length[0]    = b->dims[0];
    length[1]    = b->dims[1];
    add_side(local_tree,   b, 4, length, 0+offset[2], &offset_bd, increment, connectivity, border_case);
  } else
    b->offset[4] = offset_bd;


  //     5: Adding the back side    //

  if(border[5]){
    increment[0] = 1;
    increment[1] = dims[0];
    increment[2] = -dims[0]*dims[1];
    length[0]    = b->dims[0];
    length[1]    = b->dims[1];
    add_side(local_tree, b, 5, length, dims[0]*dims[1]*(dims[2] - 1)+offset[2], &offset_bd, increment, connectivity, border_case);//tochange
  } else {
    b->offset[5] = offset_bd;
  }

  b->offset[6]    = offset_bd ;

  //     6: Parents      //
  b = add_ancestors(local_tree, b, border_case); //tochange

  //       Shrinking     //
  b->array           = realloc(b->array,           b->size_curr * sizeof(BoundaryNode));
  if(border_case)
    b->attribute       = realloc(b->attribute,       b->size_curr * b->size_attr);
  b->border_par      = realloc(b->border_par,      b->size_curr * sizeof(BorderIndex));

  // Allocating the remaining variables //
  if(border_case != 1){
    b->border_lr 	     = malloc(b->size_curr * sizeof(BorderIndex)); check_alloc(b->border_lr, 223);

    #pragma omp parallel for 
    for (ulong i = 0; i < b->size_curr; ++i) 
      b->border_lr[i] = (BorderIndex) {.b = b, .i = BOTTOM};
  } else {
    b->border_lr 	     = NULL;
  }
  b->size_init = b->size_alloc = b->size_curr;
  free(local_tree->border_idx);
  local_tree->border_idx =NULL;

  return b;
}


/* MAIN FUNCTION */

Node *correct_borders(Arguments *args, Node *local_tree, int merge_case) {
  // Initialize variables 
  
  int merged = 0;
  int connectivity = args->pixel_connectivity;
  int grid[3] = {args->mpi_grid[0], args->mpi_grid[1], args->mpi_grid[2]};
  int base[3] = {1, 1, 1};
  int grid_cur[3] = {grid[0], grid[1], grid[2]};
  int myrank = rank();
  int nprocesses = np();
  int n_merge = 0;
  ulong store_item = local_tree->size_attr;
  
  int slice = myrank / (grid[0] * grid[1]);
  int rank2D = myrank % (grid[0] * grid[1]);
  int col = rank2D % grid[0];
  int row = rank2D / grid[0];

  info("Starting MPI correction of the local trees");

  // Calculate the number of merges
  n_merge = calculate_merges(grid, grid_cur, base, col, row, slice, myrank, nprocesses);

  // Allocate boundary tree array
  Boundary **bound_tree = malloc((1 + 2 * n_merge) * sizeof(Boundary *));
  check_alloc(bound_tree, 200);
  bound_tree[0] = create_boundary(local_tree, connectivity, merge_case);

  debug("Initial boundary tree construction done: wallclock time = %0.2f",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

  //Resetting changed variables
  base[0] =  base[1] =  base[2] = 1;
  grid_cur[0] = grid[0];
  grid_cur[1] = grid[1];
  grid_cur[2] = grid[2];
  n_merge = 0;


  // Perform the merging process
  merge_boundary_trees(grid, grid_cur, base, col, row, slice, bound_tree, &n_merge, store_item, &merged, merge_case);

  timing("Boundary trees merged: wallclock time = %0.2f",
        (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));

  //Resetting changed variables
  base[0] /= 2;
  base[1] /= 2;
  base[2] /= 2;
  grid_cur[0] = grid[0];
  grid_cur[1] = grid[1];
  grid_cur[2] = grid[2];
  merged--;

  // Update boundary trees recursively
  update_boundary_trees(grid, grid_cur, base, col, row, slice, bound_tree, &n_merge, &merged, store_item, merge_case);

  debug("Boundary tree size: Init %ld, Final %ld", bound_tree[0]->size_init, bound_tree[0]->size_curr);

  // Final tree correction and cleanup
  local_tree = correct_local_tree(local_tree, bound_tree[0], merge_case);

  free_boundary(bound_tree[0]);
  timing("Boundary trees and local trees updated: wallclock time = %0.2f",
         (float)(times(&tstruct) - start) / (float)sysconf(_SC_CLK_TCK));

  return local_tree;
} /* correct_borders */

/* HELPER FUNCTIONS */
BorderIndex idx_i(Boundary *a, ulong c, Direction d)
{
    /* Return node from a to be merged */
    BorderIndex bi;
    bi.b = a;

    ulong offset_index;

    // Assuming HORIZONTAL is the most common case
    if (d == HORIZONTAL)
        offset_index = a->offset[1];
    else if (d == VERTICAL)
        offset_index = a->offset[2];
    else
        offset_index = a->offset[5];

    bi.i = a->merge_idx[offset_index + c];

    return bi;
}

BorderIndex idx_j(Boundary *b, ulong c, Direction d)
{
    /* Return node from a to be merged */
    BorderIndex bi;
    bi.b = b;

    ulong offset_index;

    // Assuming HORIZONTAL is the most common case
    if (d == HORIZONTAL)
        offset_index = b->offset[3];
    else if (d == VERTICAL)
        offset_index = b->offset[0];
    else
        offset_index = b->offset[4];

    bi.i = b->merge_idx[offset_index + c];

    return bi;
}

void reset_border_idx(Boundary *b)
{
#pragma omp parallel for
  for (ulong i = 0; i < b->size_curr; i++)
    b->array[i].border_idx = BOTTOM;
} /* reset_border_idx */

BorderIndex b_levelroot(BorderIndex bi)
{
  BorderIndex ri = bi;
  if (is_bottom(ri))
    return bi; /* still BOTTOM */
  else
  {
    value gv = b_gval(bi);

    while (!is_bottom(b_parent(ri)) && (gv == b_gval(b_parent(ri))))
    {
      ri = b_parent(ri);
      // if (bi_equal(ri, b_parent(ri)))
      // {
      //   error("[CONNEXION ERROR] Next levelroot %ld in %d");
      //   MPI_Abort(MPI_COMM_WORLD, 240);
      // }
    }
    return ri;
  }
} /* b_levelroot */

BorderIndex b_id_levelroot(BorderIndex bi)
{
  if (is_bottom(bi))
    return bi; /* still BOTTOM */

  BorderIndex ri = bi;
  BorderIndex pi = bi;
  value gv = b_gval(bi);

  while (!is_bottom(b_parent(pi)) && gv == b_gval(b_parent(pi)) && b_parent(pi).b == bi.b)
  {
    pi = b_parent(pi);
    ri = pi;
  }

  return ri;
} /* b_id_levelroot */

bool is_bottom(BorderIndex bi)
{
  return bi.i == BOTTOM;
} /* is_bottom */

bool bi_equal(BorderIndex ai, BorderIndex bi)
{
  return (ai.b == bi.b) && (ai.i == bi.i);
} /* bi_equal */

BorderIndex b_parent_lr(BorderIndex bi)
{
  /* get parent index */
  BorderIndex pi = b_parent(bi);
  /* get levelroot of parent */
  BorderIndex ri = b_levelroot(pi);
  // trace("b_parent_lr pi: %p %ld ri: %p %ld", pi.i, (void *)pi.b, ri.i, (void *)ri.b);
  return ri;
} /* b_parent_lr */

BorderIndex b_parent(BorderIndex bi)
{
  BorderIndex pi = bi.b->border_par[bi.i];
  return pi;
} /* b_parent */

BoundaryNode b_node(BorderIndex bi)
{
  // if (bi.b == NULL)
  // {
  //   error("asking node of non-initilized border %ld", bi.i);
  // }
  return (bi.b)->array[bi.i];
}

value b_gval(BorderIndex bi)
{
  // if (bi.b == NULL)
  // {
  //   error("asking gval of non-initilized border");
  //   MPI_Abort(MPI_COMM_WORLD, 241);
  // }
  // if ((ulong)bi.i > ((bi.b)->size_curr))
  // {
  //   warn("asking gval of index (%ld) higher than size (%ld) of border", bi.i, (bi.b)->size_curr);
  //   return BOTTOM;
  // }

  return ((bi.b)->array[bi.i]).gval;
} /* b_gval */

void free_boundary(Boundary *b)
{
  if (b != NULL)
  {
    free(b->array);
    free(b->attribute);
    free(b->merge_idx);
    free(b->border_par);
    free(b->border_lr);
    free(b->border_ori);
  }
  free(b);
} /* free_boundary */
