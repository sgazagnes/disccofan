
#include "types.h"
#include "attributes.h"
#include "flood.h"
#include "queue.h"
#include "mpihelper.h"
#include "workspace.h"


struct tms 	tstruct;				
clock_t 	start;
ulong 		g_max_levels;
float 		g_max_greyval;

void tree_flood(idx *parent,  void *attribute,  BitArray *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims,  ulong lwb, ulong upb,  int connectivity, int size_att);
void tree_flood_par(idx *parents,  BitArray *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims, ulong lwb, ulong upb, int connectivity);

void create_ranks_inv(value *gvals, ulong **ranks_inv, ulong histos[][MXT_HISTO_SZ], ulong lwb, ulong upb, value_t (*functor)(const value));

BitArray *create_bit_array(ulong size, Allocator *work_alloc);
void bit_array_set(ulong *data, ulong index);
bool bit_array_get(ulong *data, ulong index);
void bit_array_free(BitArray *bit_array);

void fuse_sections(Node *tree,   ulong *dims, uint id, uint i, uint n_t, int connectivity);
void fuse_parents(Node *tree, ulong *dims, uint id, uint i, uint n_t, int connectivity);
void merge_nodes(Node *tree, idx x, idx y);
void merge_nodes_par(Node *tree,idx x, idx y);

int get_neighbors(ulong *dims, ulong lwb, ulong upb, ulong *neighbors, ulong p, ulong x, ulong y, ulong z,  int connectivity);
bool check_neighbor(BitArray *visited, PrioQueue *q, ulong *ranks, ulong* index, ulong* x, ulong* y, ulong *z, ulong* rank, long offset, ulong n_x, ulong n_y, ulong n_z, ulong lwb);
void remaining(BitArray *visited, PrioQueue *q, ulong *ranks, ulong* index, ulong* rank, ulong* dims, ulong lwb, ulong upb, int connectivity);

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Tree Building         */
/*				   */
/* +++++++++++++++++++++++++++++++ */


void build_local_tree(Arguments *args, Node *local_tree, ulong *dims){
  int connectivity =  args->connectivity_arg;
  int nthreads     =  args->threads_arg;
  int *saval 	   = calloc(nthreads, sizeof(int));

  omp_lock_t lock[nthreads];
  for (int i=0; i<nthreads; i++)
    omp_init_lock(&(lock[i]));

  local_tree->parent = NULL;
  local_tree->attribute = NULL;
  ulong *ranks 	= NULL;

  if(nthreads > 1 && args->bpp_arg > 16)
    warn("Using multitle threads with high dynamic range dataset might not work well");
  
  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			/* Thread number */
    ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads) 
      : dims[0]*dims[1]*((id*dims[2])/nthreads) ;              /* Lower bound for current thread */
    ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/nthreads)
      : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);       /* Upper bound for current thread */

    // Create mappings
    
    ulong *ranks_inv  	= NULL;
    create_mappings(local_tree->gval, &ranks_inv, lwb, upb);

    #pragma omp single
    {
      ranks 	       		= malloc(local_tree->size_curr * sizeof(ulong));
      local_tree->parent 	= (idx *) ranks;
      local_tree->attribute 	= malloc(local_tree->size_curr * local_tree->size_attr);
    }
    
    for (ulong i = 0; i != upb-lwb; ++i){
      ranks[ranks_inv[i]+lwb] = i;
    }

    /*    for (ulong i = 0; i != 100; ++i){
      info("%d, %1.50f, %ld", i, local_tree->gval[i], ranks[i]);
      }*/
    //   info("%ld, %1.15f",ranks_inv[upb-lwb-1], local_tree->gval[ranks_inv[upb-lwb-1]]);
    info("Mappings: wallclock time = %0.2f ",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

    // Start construction
    
    Allocator work_alloc;
    allocator_init(&work_alloc, workspace_construct(upb-lwb));
    PrioQueue *q        = create_prio_queue(upb-lwb-1, &work_alloc);
    BitArray *visited   = create_bit_array(upb-lwb, &work_alloc);
    
    init_attrib_array(local_tree->attribute, local_tree->gval, local_tree->border, dims, local_tree->offsets, lwb, upb, local_tree->size_attr);
    
    tree_flood(local_tree->parent, local_tree->attribute, visited, q, ranks+lwb, ranks_inv, dims, lwb, upb, connectivity, local_tree->size_attr);
    
    free(ranks_inv);
    allocator_free(&work_alloc);

    // Tree compression
    
    for (ulong x =lwb ; x < upb; x++){
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
    }

    // Merging threads if needed
    
    if (nthreads > 1) {
      int i = 1;
      int qq = id;
      while (id+i < nthreads && qq%2 == 0){
		
	while (saval[id+i] <= 0){
	  #pragma omp flush // wait
	}
	omp_set_lock(&lock[id+i]);
	saval[id+i]--;
	omp_unset_lock(&lock[id+i]);
	fuse_sections(local_tree, dims, id, i, nthreads, connectivity);
	i *= 2;
	qq /= 2;      
      }
      if(id>0){
	omp_set_lock(&lock[id]);	      
	saval[id]++;
	omp_unset_lock(&lock[id]);
      }

    }
    
  }
    
  free(saval);
  for (int i = 0; i < nthreads; i++)
    omp_destroy_lock(&(lock[i]));

} /* build_local_tree */




void build_local_tree_par(Arguments *args, Node *local_tree, ulong *dims){
  int connectivity =  args->connectivity_arg;
  int nthreads     =  args->threads_arg;
  int *saval 	   = calloc(nthreads, sizeof(int));
  
  omp_lock_t lock[nthreads];
  for (int i=0; i<nthreads; i++)
    omp_init_lock(&(lock[i]));
  
  local_tree->parent = NULL;
  ulong *ranks 	= NULL;

  if(nthreads > 1 && args->bpp_arg > 16)
    warn("Using multitle threads with high dynamic range dataset might not work well");


  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			/* Thread number */
    ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads)
      : dims[0]*dims[1]*((id*dims[2])/nthreads);              /* Lower bound for current thread */
    ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/nthreads)
      : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);;       /* Upper bound for current thread */

    ulong *ranks_inv  	= NULL;
    create_mappings(local_tree->gval, &ranks_inv, lwb, upb);

    #pragma omp single
    {
      ranks 	       		= malloc(local_tree->size_curr * sizeof(ulong));
      local_tree->parent 	= (idx *) ranks;
    }
    
    for (ulong i = 0; i != upb-lwb; ++i){
      ranks[ranks_inv[i]+lwb] = i;
    }

       
    info("Mappings: wallclock time = %0.2f ",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

    // Start construction
    
    Allocator work_alloc;
    allocator_init(&work_alloc, workspace_construct(upb-lwb));
    PrioQueue *q        = create_prio_queue(upb-lwb-1, &work_alloc);
    BitArray *visited   = create_bit_array(upb-lwb, &work_alloc);
   
    
    tree_flood_par(local_tree->parent,  visited, q, ranks+lwb, ranks_inv, dims, lwb, upb, connectivity);
    
    free(ranks_inv);
    allocator_free(&work_alloc);

    // Tree compression
    
    for (ulong x =lwb ; x < upb; x++){
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
    }

    if (nthreads > 1) {
      /* Merge threads sections */
      int i = 1;
      int qq = id;
      while (id+i < nthreads && qq%2 == 0){
		
	while (saval[id+i] <= 0){
	  #pragma omp flush // wait
	}
	omp_set_lock(&lock[id+i]);
	saval[id+i]--;
	omp_unset_lock(&lock[id+i]);

	fuse_parents(local_tree, dims, id, i, nthreads, connectivity);
	
	i *= 2;
	qq /= 2;      
      }
      if(id>0){
	omp_set_lock(&lock[id]);	      
	saval[id]++;
	omp_unset_lock(&lock[id]);
      }
     
    }
  }
  free(saval);
  for (int i = 0; i < nthreads; i++)
    omp_destroy_lock(&(lock[i]));

} /* build_local_tree_par */

void build_local_tree_att(Node *local_tree, ulong *ranks,  ulong *dims){

  bool *visited     = calloc(local_tree->size_curr, sizeof(bool));  check_alloc(visited, 603);
  local_tree->attribute 	= malloc(local_tree->size_curr * local_tree->size_attr);

  init_attrib_array(local_tree->attribute, local_tree->gval, local_tree->border, dims, local_tree->offsets, 0, local_tree->size_curr, local_tree->size_attr);
  for( long i = local_tree->size_curr-1; i >= 0; i--) {
    idx p = ranks[i];
    visited[p] = true;
    idx parent = get_parent(local_tree, p);
    if(parent == BOTTOM) continue;

    merge_aux_data(local_tree->attribute + parent*local_tree->size_attr, local_tree->attribute + p*local_tree->size_attr);

    if(local_tree->gval[p] == local_tree->gval[parent] && visited[parent]){
      parent = get_parent(local_tree, parent);
      if(parent != BOTTOM){
	    
	merge_aux_data(local_tree->attribute + parent*local_tree->size_attr, local_tree->attribute + p*local_tree->size_attr);
	  
      }
    }
  }
  free(visited);

  
}
    // Start construction

  
  


void tree_flood(idx *parents, void *attribute,  BitArray *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims, ulong lwb, ulong upb, int connectivity, int size_att){

  ulong index = lwb;
  ulong rank = ranks[index-lwb];
  bit_array_set(visited->data, index-lwb);
  while(true){
    remaining(visited, q, ranks, &index, &rank, dims, lwb, upb, connectivity);

    if (q->m_levels[0][0] == 0) break;
    rank = q->m_top;
    ulong parent = ranks_inv[rank]+lwb;
    prio_queue_remove(q);

    merge_aux_data((attribute + parent*size_att), (attribute + index*size_att));
    parents[index] = parent;
    index = parent;

  }
  parents[index] = BOTTOM;
} /* tree_flood */

void tree_flood_par(idx *parents,  BitArray *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims, ulong lwb, ulong upb, int connectivity){
  
  ulong index = lwb; 
  ulong rank = ranks[index-lwb];
  bit_array_set(visited->data, index-lwb);
  while(true){
    remaining(visited, q, ranks, &index, &rank, dims, lwb, upb, connectivity);

    if (q->m_levels[0][0] == 0) break;
    rank = q->m_top;
    prio_queue_remove(q);
    ulong parent = ranks_inv[rank]+lwb;

    parents[index] = parent;
    index = parent;
  }
  parents[index] = BOTTOM;
} /* tree_flood_tee_par */

/* ++++++++++++++++++++++++++ */
/*      Threads functions     */
/* ++++++++++++++++++++++++++ */

void fuse_sections(Node *tree,  ulong *dims, uint id, uint i, uint n_t, int connectivity){
  ulong inc  = dims[2] == 1 ? dims[0]: dims[0]*dims[1];
  ulong mdb  = dims[2] == 1 ? dims[0]*(((id+i) * dims[1]) / n_t) : dims[0]*dims[1]*(((id+i) * dims[2]) / n_t);

  ulong *ranks_inv  	= NULL;
  create_mappings(tree->gval, &ranks_inv, mdb - inc, mdb + inc);
  ulong *ranks 	       	= malloc(2*inc * sizeof(ulong));
  
  for (ulong ii = 0; ii != 2*inc; ++ii){
    ranks[ranks_inv[ii]] = ii;
  }

  Allocator work_alloc;
  allocator_init(&work_alloc, workspace_construct(2*inc));
  PrioQueue *q        = create_prio_queue(2*inc-1, &work_alloc);
  BitArray *visited   = create_bit_array(2*inc, &work_alloc);

  ulong index = mdb-inc;
  ulong rank = ranks[index-mdb+inc];
  bit_array_set(visited->data, index-mdb+inc);
  while(true){
    remaining(visited, q, ranks, &index, &rank, dims, mdb-inc, mdb+inc, connectivity);
    if (q->m_levels[0][0] == 0) break;
    rank = q->m_top;
    ulong parent = ranks_inv[rank]+mdb-inc;
    prio_queue_remove(q);
    merge_nodes(tree, index, parent);    
    index = parent;
  }

  allocator_free(&work_alloc);
  free(ranks);
  free(ranks_inv);
}

/*
void fuse_sections(Node *tree,  ulong *dims, uint id, uint i, uint n_t, int connectivity){
  ulong p, u, v, u_x, u_y;
  value min_curr, min_prev;
  bool  test_min = 0;
  ulong inc  = dims[2] == 1 ? dims[0]: dims[0]*dims[1];
  ulong offb  = 2*inc*(id + i);
  ulong offa  = 2*inc*(id + i-1);
  ulong mdb  = dims[2] == 1 ? dims[0]*(((id+i) * dims[1]) / n_t) : dims[0]*dims[1]*(((id+i) * dims[2]) / n_t);

  for (p=0, u = mdb; p < inc; p++, u++){
    u_x = u % dims[0];
    v = u-inc;
    if(u_x == 0) test_min = false;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (!test_min || min_curr > min_prev)
      merge_nodes(tree, u, v);
    min_prev = min_curr;
    test_min = true;
    if(connectivity >= 8){
      if(u_x > 0){
	min_curr = MIN(tree->gval[u], tree->gval[v-1]);
	if (min_curr > min_prev)
	  merge_nodes(tree, u, v-1);
      }     
      if(u_x < dims[0]-1){
	min_curr = MIN(tree->gval[u], tree->gval[v+1]);
	if (!test_min || min_curr > min_prev)
	  merge_nodes(tree, u, v+1);
      }
    }
    if(connectivity >= 26){
      u_y = (u % inc) / dims[0];

      if(u_x > 0 && u_y > 0){
	min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]-1]);
	if (min_curr > min_prev )
	  merge_nodes(tree, u, v-dims[0]-1);
      }
      if(u_y > 0){
	min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]]);
	if (min_curr > min_prev )
	  merge_nodes(tree, u, v-dims[0]);
      }
      if(u_y > 0 && u_x<dims[0]-1){
	min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]+1]);
	if (min_curr > min_prev )
	  merge_nodes(tree, u, v-dims[0]+1);
      }
      if(u_y < dims[1]-1 && u_x>0){
	min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]-1]);
	if (!test_min || min_curr > min_prev )
	  merge_nodes(tree, u, v+dims[0]-1);
      }
      if(u_y < dims[1]-1){
	min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]]);
	if (!test_min || min_curr > min_prev )
	  merge_nodes(tree, u, v+dims[0]);
      }
      if(u_y < dims[1]-1 && u_x<dims[0]-1){
	min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]+1]);
	if (!test_min || min_curr > min_prev )
	  merge_nodes(tree, u, v+dims[0]+1);
      }
    }
    
  }
}
*/


void fuse_parents(Node *tree, ulong *dims, uint id, uint i, uint n_t, int connectivity){
  ulong inc  = dims[2] == 1 ? dims[0]: dims[0]*dims[1];
  ulong mdb  = dims[2] == 1 ? dims[0]*(((id+i) * dims[1]) / n_t) : dims[0]*dims[1]*(((id+i) * dims[2]) / n_t);

  ulong *ranks_inv  	= NULL;
  create_mappings(tree->gval, &ranks_inv, mdb - inc, mdb + inc);
  ulong *ranks 	       	= malloc(2*inc * sizeof(ulong));
  
  for (ulong ii = 0; ii != 2*inc; ++ii){
    ranks[ranks_inv[ii]] = ii;
  }

  Allocator work_alloc;
  allocator_init(&work_alloc, workspace_construct(2*inc));
  PrioQueue *q        = create_prio_queue(2*inc-1, &work_alloc);
  BitArray *visited   = create_bit_array(2*inc, &work_alloc);

  ulong index = mdb-inc;
  ulong rank = ranks[index-mdb+inc];
  bit_array_set(visited->data, index-mdb+inc);
  while(true){
    remaining(visited, q, ranks, &index, &rank, dims, mdb-inc, mdb+inc, connectivity);
    if (q->m_levels[0][0] == 0) break;
    rank = q->m_top;
    ulong parent = ranks_inv[rank]+mdb-inc;
    prio_queue_remove(q);
    merge_nodes_par(tree, index, parent);    
    index = parent;
  }

  allocator_free(&work_alloc);
  free(ranks);
  free(ranks_inv);
}




void merge_nodes(Node *tree, idx x, idx y) {
  void *cor  = NULL;
  void *copa = NULL;
  void *attr = NULL;
  idx h, z;

  if ( tree->gval[y] > tree->gval[x]){
    h=x; x=y; y=h;
  }
  x =  get_levelroot(tree,x);
  y =  get_levelroot(tree,y);
  
  while ((x != y) && (y != BOTTOM)) {
    z = get_parent(tree, x);
    if ((z != BOTTOM) && (tree->gval[z]>=tree->gval[y])) {
      if (cor) {
	merge_aux_data(tree->attribute + x*tree->size_attr, cor);
      }
      x = z;
    } else {
      if(cor)
	merge_to_aux_data(NULL, &copa,tree->attribute + x*tree->size_attr, cor);
      else
	clone_aux_data(NULL, &copa, tree->attribute + x*tree->size_attr);
      clone_aux_data(NULL, &cor, tree->attribute + x*tree->size_attr);
      
      if(copa){
	attr = tree->attribute + x*tree->size_attr;
      	clone_aux_data(NULL, &attr, copa);
      }
      tree->parent[x] = y ;
      x = y;
      y = z;
    }
  }

  if (y == BOTTOM && cor) {
    while(x != BOTTOM) {      
      merge_aux_data(tree->attribute + x*tree->size_attr, cor);

      x = get_parent(tree, x);
    }
  }
  if (cor)  delete_aux_data(cor);
  if (copa) delete_aux_data(copa);
}/* merge_nodes */



void merge_nodes_par(Node *tree,idx x, idx y) {
  idx h, z;

  x =  get_levelroot(tree, x);
  y =  get_levelroot(tree, y);

  if (tree->gval[x] < tree->gval[y]) {
    h=x; x=y; y=h;
  }
  while ((x != y) && (y != BOTTOM)) {
    z = get_parent(tree, x);

    if ((z != BOTTOM) && (tree->gval[z]>=tree->gval[y])) {
      x = z;
    } else {
      tree->parent[x] = y;
      x = y;
      y = z;
    }
  }
}/* merge_nodes_par */


/* +++++++++++++++++++++++++++++++ */
/*				   */
/*         Annexes Functions       */
/*				   */
/* +++++++++++++++++++++++++++++++ */


inline value_t transform_float(const float val) {    
  value_t valu = *((const value_t *)(&val));

  if (valu &    0x80000000)
    return      0xFFFFFFFF - valu;
        
  return valu | 0x80000000;
}

inline value_t transform_dum(const value val)
{
  return val;
}


void gen_histogram(value *gvals, ulong histos[][MXT_HISTO_SZ],  ulong lwb, ulong upb, value_t (*functor)(const value)) {
  memset(histos, 0, NUM_DIGITS * MXT_HISTO_SZ * sizeof(histos[0][0]));
  
  for (ulong i = lwb; i != upb; ++i) {
    int shift = 0;
    value_t gval = functor(gvals[i]);
    for (int digit_nr = 0; digit_nr != NUM_DIGITS; ++digit_nr) {
      int digit = (gval >> shift) & MXT_HISTO_MASK;
      ++histos[digit_nr][digit];
      shift += MXT_HISTO_SZ_LOG2;
    }
  }
}


void exclusive_sum(ulong *it, ulong *it_end){
  ulong sum = 0;
  while (it != it_end) {
    ulong next = *it;
    *it++ = sum;
    sum += next;
  }
}


NO_INLINE void create_mappings(value *gvals, ulong **ranks_inv, ulong lwb, ulong upb) {
  const int   num_digits = NUM_DIGITS;
  value_t (*foo)(const value);

  if(FLOAT_TYPE == 1)
    foo = &transform_float;
  else
    foo = &transform_dum;

  
  ulong histos[num_digits][MXT_HISTO_SZ];
  gen_histogram(gvals, histos, lwb, upb, foo);

  for (int digit_nr = 0; digit_nr != num_digits; digit_nr++) 
    exclusive_sum(histos[digit_nr], histos[digit_nr] + MXT_HISTO_SZ);
  
  create_ranks_inv(gvals, ranks_inv, histos, lwb, upb, foo);
  /* *ranks = (ulong *)  allocator_allocate(work_alloc, upb-lwb, sizeof(ulong));
  
  for (ulong i = 0; i != upb-lwb; ++i)
  (*ranks)[(*ranks_inv)[i]] = i;*/
  
}



void scatter_first_digit(value *gvals, SortItem *pair_it, ulong* histo, ulong lwb, ulong upb, value_t (*functor)(const value)  ){
  for (ulong i = lwb; i != upb; ++i) {
    value_t gval; 
    gval = functor(gvals[i]);
    int digit = gval & MXT_HISTO_MASK;
    pair_it[histo[digit]++] = (SortItem) {gval, i-lwb};
  }
}

void scatter_digit_bit1(value *gvals, ulong *out, ulong* histo, ulong lwb, ulong upb, value_t (*functor)(const value))
{   
  for (ulong i = lwb; i != upb; ++i)
  {
    int digit = functor(gvals[i]) & MXT_HISTO_MASK;
    out[histo[digit]++] = i-lwb;
  }
}

void scatter_digit(SortItem *in, SortItem *out, int digit_nr, ulong* histo, ulong size) {  
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i) {
    const SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair;
  }
}

void scatter_last_digit(SortItem *in, ulong *out, int digit_nr, ulong* histo, ulong size){
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i)  {
    const SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair.rank;
  }
}


void create_ranks_inv(value *gvals, ulong **ranks_inv, ulong histos[][MXT_HISTO_SZ], ulong lwb, ulong upb, value_t (*functor)(const value)) {

  const int num_digits = NUM_DIGITS;
  //  size_t alloc_pos =allocator_pos(work_alloc);
  // *ranks_inv = (ulong *)  allocator_allocate(work_alloc, upb-lwb, sizeof(ulong));
  //  size_t alloc_pos_rank_inv = allocator_pos(work_alloc);

  if (num_digits == 1) {
    *ranks_inv = calloc(upb-lwb, sizeof(ulong));
    scatter_digit_bit1(gvals, *ranks_inv, histos[0], lwb, upb, functor);
    return;
  }

  //allocator_set_pos(work_alloc, alloc_pos);
  SortItem *pairs1 =  malloc((upb-lwb)*sizeof(SortItem));//*ranks_inv;
  //allocator_allocate(work_alloc, upb-lwb, sizeof(SortItem));
  SortItem *pairs2 =  malloc((upb-lwb)*sizeof(SortItem));//*ranks_inv + (upb-lwb)*sizeof(SortItem);
  //allocator_allocate(work_alloc, upb-lwb, sizeof(SortItem));

  SortItem *pairswap = NULL;
  scatter_first_digit(gvals, pairs2, histos[0], lwb, upb, functor);

  int digit_nr = 1;
  for (; digit_nr != NUM_DIGITS - 1; ++digit_nr)
  {
    pairswap = pairs1;
    pairs1 = pairs2;
    pairs2 = pairswap;
    
    scatter_digit(pairs1, pairs2, digit_nr, histos[digit_nr], (upb-lwb));
  }
  
  free(pairs1);
  *ranks_inv = calloc(upb-lwb, sizeof(ulong));
  scatter_last_digit(pairs2, (ulong *) *ranks_inv, digit_nr, histos[digit_nr], (upb-lwb));
  free(pairs2);
  // allocator_set_pos(work_alloc, alloc_pos_rank_inv);
}



/*          Bit Array         */


BitArray *create_bit_array(ulong size, Allocator *work_alloc) {
  BitArray *bit_array = malloc(sizeof(BitArray)); check_alloc(bit_array, 507);
  bit_array->num_words = (size + bits_per_word_log2() - 1) / bits_per_word_log2();
  bit_array->data = allocator_allocate(work_alloc, bit_array->num_words, sizeof(ulong));
  memset(bit_array->data, 0, bit_array->num_words * sizeof(ulong));
    //calloc(bit_array->num_words, sizeof(ulong)); check_alloc( bit_array->data, 508);
  return bit_array;
}

void bit_array_set(ulong *data, ulong index) {
  ulong word_idx = index >> bits_per_word_log2();
  int bit_idx = index & (bits_per_word() - 1);
  data[word_idx] |= (ulong)(1) << bit_idx; 
}

bool bit_array_get(ulong *data, ulong index){
  ulong word_idx = index >> bits_per_word_log2();
  int bit_idx = index & (bits_per_word() - 1);
  return !!(data[word_idx] & ((ulong)(1) << bit_idx));
}

void bit_array_free(BitArray *bit_array) {
  free(bit_array->data);
  free(bit_array);
}

int bit_scan_reverse(ulong val) {
  return sizeof(val) * CHAR_BIT - 1 - __builtin_clzl(val);
}
int bits_per_word_log2(void) {
  return sizeof(unsigned) * CHAR_BIT - 1 - __builtin_clz((sizeof(ulong) * CHAR_BIT));
}
int bits_per_word(void) { 
  return sizeof(ulong) * CHAR_BIT;
}


/* Misc */


bool check_neighbor(BitArray *visited, PrioQueue *q,  ulong *ranks, ulong* index, ulong* x, ulong* y, ulong *z, ulong* rank, long offset, ulong n_x, ulong n_y, ulong n_z, ulong lwb){

  ulong n = *index + offset;
    
  if(bit_array_get(visited->data,  n-lwb))
    return false;

  bit_array_set(visited->data,  n-lwb);


  ulong rank_n = ranks[n-lwb];
  if (*rank > rank_n) {
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



int get_neighbors(ulong *dims, ulong lwb, ulong upb, ulong *neighbors, ulong p, ulong x, ulong y, ulong z, int connectivity){
  int n = 0;
    
  if ((x < dims[0]-1) && (p+1 < upb))        neighbors[n++] = p+1;
  if ((y > 0) && (p - dims[0] >= lwb))       neighbors[n++] = p-dims[0];
  if ((x>0) && (p-1 >= lwb))                 neighbors[n++] = p-1;
  if ((y < dims[1]-1) && (p+dims[0] < upb))  neighbors[n++] = p+dims[0];

  if(connectivity == 8 || connectivity == 26){
    if ((x < dims[0]-1) && (y > 0) && (p+1-dims[0] >=lwb))          neighbors[n++] = p+1-dims[0];
    if ((y > 0) && (x>0) && (p-dims[0]-1 >= lwb))       	    neighbors[n++] = p-dims[0]-1;
    if ((x > 0) && (y < dims[1]-1) && (p-1+dims[0] < upb))          neighbors[n++] = p-1+dims[0];
    if ((y < dims[1]-1) && (x < dims[0]-1) && (p+dims[0]+1 < upb))  neighbors[n++] = p+dims[0]+1;
  }
  
  if (dims[2] > 1 && (connectivity == 6 || connectivity ==26)) {
    if ((z > 0) && (p >= lwb + dims[0]*dims[1]))       neighbors[n++] = p-dims[0]*dims[1];
    if ((z < dims[2]-1) && (p+dims[0]*dims[1] < upb))  neighbors[n++] = p+dims[0]*dims[1];
     if(connectivity == 26) {
 
      if ((z>0) && (y>0) && (x>0) && (p - dims[0]*dims[1]-dims[0]-1>= lwb))
	neighbors[n++] = p-dims[0]*dims[1]-dims[0]-1;
      if ((z>0) && (y>0) && (p - dims[0]*dims[1]-dims[0]>= lwb))
	neighbors[n++] = p-dims[0]*dims[1]-dims[0];
      if ((z>0) && (y>0) && (x < dims[0]-1) && (p - dims[0]*dims[1]-dims[0]+1>= lwb))
	neighbors[n++] = p-dims[0]*dims[1]-dims[0]+1;
      if ((z>0) && (x>0) && (p - dims[0]*dims[1]-1>= lwb))
	neighbors[n++] = p-dims[0]*dims[1]-1;
      if ((z>0) && (x < dims[0]-1) && (p - dims[0]*dims[1]+1 >= lwb))
	neighbors[n++] = p-dims[0]*dims[1]+1;
      if ((z>0) && (y < dims[1]-1) && (x>0) && (p - dims[0]*dims[1]+dims[0]-1>= lwb))
	neighbors[n++] = p-dims[0]*dims[1]+dims[0]-1;
      if ((z>0) && (y < dims[1]-1) && (p - dims[0]*dims[1]+dims[0]>= lwb))
	neighbors[n++] = p-dims[0]*dims[1]+dims[0];
      if ((z>0) && (y < dims[1]-1) && (x < dims[0]-1) && (p - dims[0]*dims[1]+dims[0]+1>= lwb))
	neighbors[n++] = p-dims[0]*dims[1]+dims[0]+1;
      
      if ((z < dims[2]-1) && (y>0) && (x>0) && (p + dims[0]*dims[1]-dims[0]-1 < upb))
	neighbors[n++] = p+dims[0]*dims[1]-dims[0]-1;
      if ((z < dims[2]-1) && (y>0) && (p + dims[0]*dims[1]-dims[0]< upb))
	neighbors[n++] = p+dims[0]*dims[1]-dims[0];
      if ((z < dims[2]-1) && (y>0) && (x < dims[0]-1) && (p + dims[0]*dims[1]-dims[0]+1<upb ))
	neighbors[n++] = p+dims[0]*dims[1]-dims[0]+1;
      if ((z < dims[2]-1) && (x>0) && (p + dims[0]*dims[1]-1< upb))
	neighbors[n++] = p+dims[0]*dims[1]-1;
      if ((z < dims[2]-1) && (x < dims[0]-1) && (p + dims[0]*dims[1]+1 < upb))
	neighbors[n++] = p+dims[0]*dims[1]+1;
      if ((z < dims[2]-1) && (y < dims[1]-1) && (x>0) && (p + dims[0]*dims[1]+dims[0]-1< upb))
	neighbors[n++] = p+dims[0]*dims[1]+dims[0]-1;
      if ((z < dims[2]-1) && (y < dims[1]-1) && (p + dims[0]*dims[1]+dims[0]< upb))
	neighbors[n++] = p+dims[0]*dims[1]+dims[0];
      if ((z < dims[2]-1) && (y < dims[1]-1) && (x < dims[0]-1) && (p+dims[0]*dims[1]+dims[0]+1< upb))
	neighbors[n++] = p+dims[0]*dims[1]+dims[0]+1;
     }
  
  }
  
  return(n);
} /* get_neighbors */



void remaining(BitArray *visited, PrioQueue *q,  ulong *ranks, ulong* index, ulong* rank, ulong* dims, ulong lwb, ulong upb, int connectivity){
  ulong z = *index / (dims[0] * dims[1]);
  ulong x = (*index % (dims[0] * dims[1])) % dims[0];
  ulong y = (*index % (dims[0] * dims[1])) / dims[0];

  bool cond = 1;

  while(cond){
    cond = 
      ((x > 0)           && (*index -1 >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -1, x - 1, y, z, lwb)) ||
      ((x < dims[0] - 1) && (*index + 1 < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) 1, x + 1, y, z, lwb)) ||
      ((y > 0)           && (*index - dims[0] >= lwb)       && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, -dims[0], x, y - 1, z, lwb)) ||
      ((y < dims[1] - 1)  && (*index + dims[0] < upb) && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, dims[0], x, y + 1, z, lwb)) ||
      ((z > 0)  && (*index >= lwb + dims[0]*dims[1])  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, -dims[1]*dims[0], x, y, z - 1, lwb)) ||
      ((z < dims[2] - 1)  && (*index+dims[0]*dims[1] < upb) && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, +dims[1]*dims[0], x, y, z + 1, lwb));
    
    if(connectivity >= 8){
      cond = cond ||
	((x > 0) && (y > 0) && (*index-1-dims[0] >=lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -1-dims[0], x - 1, y - 1, z, lwb)) ||
	((x < dims[0]-1) && (y > 0) && (*index+1-dims[0] >=lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) 1-dims[0], x + 1, y - 1, z, lwb)) ||
	((x > 0) && (y < dims[1]-1) && (*index-1+dims[0] < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -1+dims[0], x-1, y + 1, z, lwb)) ||
	((y < dims[1]-1) && (x < dims[0]-1) && (*index+dims[0]+1 < upb) && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) dims[0]+1, x+1, y + 1, z, lwb)) ;   
    }
      if(connectivity == 26){
      cond = cond || (z > 0 &&
		      (((y>0) && (x>0) && (*index -dims[0]*dims[1]-dims[0]-1 >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -1-dims[0]-dims[0]*dims[1], x-1, y-1, z-1, lwb)) ||
	((y>0) && (*index -dims[0]*dims[1]-dims[0] >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -dims[0]-dims[0]*dims[1], x, y-1, z-1, lwb)) ||
	((y>0) && (x<dims[0]-1) && (*index-dims[0]*dims[1]-dims[0]+1 >= lwb) && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) +1-dims[0]-dims[0]*dims[1], x+1, y-1, z-1, lwb)) ||
	((x>0) && (*index -dims[0]*dims[1]-1 >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -1-dims[0]*dims[1], x-1, y, z-1, lwb)) ||
	((x<dims[0]-1) && (*index -dims[0]*dims[1]+1 >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) +1-dims[0]*dims[1], x+1, y, z-1, lwb)) ||
	((y<dims[1]-1) && (x>0) && (*index -dims[0]*dims[1]+dims[0]-1 >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long)-1+dims[0]-dims[0]*dims[1], x-1, y+1, z-1, lwb)) ||
      	((y<dims[1]-1)  && (*index -dims[0]*dims[1]+dims[0] >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) dims[0]-dims[0]*dims[1], x, y+1, z-1, lwb)) ||
		       ((y<dims[1]-1) && (x<dims[0]-1) && (*index -dims[0]*dims[1]+dims[0]+1 >= lwb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) +1+dims[0]-dims[0]*dims[1], x+1, y+1, z-1, lwb))))
	|| (z<dims[2]-1 &&
	    (((y>0) && (x>0) && (*index +dims[0]*dims[1]-dims[0]-1 < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -1-dims[0]+dims[0]*dims[1], x-1, y-1, z+1, lwb)) ||
	((y>0) && (*index+dims[0]*dims[1]-dims[0] < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -dims[0]+dims[0]*dims[1], x, y-1, z+1, lwb)) ||
	((y>0) && (x<dims[0]-1) && (*index+dims[0]*dims[1]-dims[0]+1  < upb) && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) +1-dims[0]+dims[0]*dims[1], x+1, y-1, z+1, lwb)) ||
	((x>0) && (*index+dims[0]*dims[1]-1  < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) -1+dims[0]*dims[1], x-1, y, z+1, lwb)) ||
	((x<dims[0]-1) && (*index +dims[0]*dims[1]+1  < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) +1+dims[0]*dims[1], x+1, y, z+1, lwb)) ||
	((y<dims[1]-1) && (x>0) && (*index +dims[0]*dims[1]+dims[0]-1  < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long)-1+dims[0]+dims[0]*dims[1], x-1, y+1, z+1, lwb)) ||
      	((y<dims[1]-1)  && (*index +dims[0]*dims[1]+dims[0]  < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) dims[0]+dims[0]*dims[1], x, y+1, z+1, lwb)) ||
	     ((y<dims[1]-1) && (x<dims[0]-1) && (*index +dims[0]*dims[1]+dims[0]+1  < upb)  && check_neighbor(visited, q,  ranks, index, &x, &y, &z, rank, (long) +1+dims[0]+dims[0]*dims[1], x+1, y+1, z+1, lwb))));
      }
  }
}

inline bool is_border(bool border[6], ulong *dims, ulong p) {
  ulong x, y, z;
  x = p % (dims[0] * dims[1]) % dims[0];
  y = p % (dims[0] * dims[1]) / dims[0];
  z = p / (dims[0] * dims[1]);
  return ((x == 0 && border[0]) || (x == dims[0]-1 && border[1]) ||
	  (y == 0 && border[2]) || (y == dims[1]-1 && border[3]) ||
	  (z == 0 && border[4]) || (z == dims[2]-1 && border[5]));
}

bool is_levelroot(Node *tree, idx x) {
  return ((tree->parent[x] == BOTTOM) || (tree->gval[x] != tree->gval[tree->parent[x]]));
} /* is_levelroot */

idx get_levelroot(Node *tree, idx x) {
  /* index based, check whether this node is bottom */
  idx r = x;
  if (r == BOTTOM)
    return BOTTOM;
  value gv = tree->gval[x]; 
  while ((tree->parent[r] != BOTTOM) && (gv == tree->gval[tree->parent[r]])) { 
    r = tree->parent[r];
  }
  /* tree compression */
  while (x != r) {
    idx y = tree->parent[x];
    tree->parent[x] = r;
    x = y;
  }
  return r;
} /* get_levelroot */


idx levelroot(Node *tree, idx index) {
  return get_levelroot(tree, index);
} /* levelroot */

idx get_parent(Node *tree, idx x) {
  return get_levelroot(tree, tree->parent[x]);
} /* get_parent */

void free_tree(Node *tree){
  free(tree->parent);
  free(tree->gval);
  free(tree->attribute);
  free(tree);
}








/******************************/

void local_hist_rs(value *gvals, ulong lwb, ulong upb, ulong *hist, ulong *sorted, ulong step){
  ushort radix;
  memset(hist, 0, sizeof(ulong)*NUMBUCKETS);
    
  if(step == 0) {
    for (ulong i=lwb; i<upb; ++i){
      radix = *( (ushort*) &(gvals[i]));			
      hist[radix]++;		
    }
  } else {
    for (ulong i=lwb; i<upb; ++i){			
      radix = *( (ushort*) (&(gvals[sorted[i]]))+step);			
      hist[radix]++;
    }
  }	     
}


void create_sorted_ars(value *gvals, ulong lwb, ulong upb, ulong *hist, ulong *sorted_new, ulong *sorted_old, ulong step){
  ushort radix;

  if(step == 0) {
    for (ulong i=lwb; i<upb; ++i) {	
      radix = *( (ushort*) &(gvals[i]));
      sorted_new[hist[radix]++] = i; // sort in ascending order
    }
  } else {
    for (ulong i=lwb; i<upb; ++i) {		
      radix = *( (ushort*) (&(gvals[sorted_old[i]]))+step);
      sorted_new[hist[radix]++] = sorted_old[i]; // sort in ascending order
    }		
  }
}


ulong *sort_image_pixels(Arguments *args, Node *tree){
 
  int nthreads = args->threads_arg;
  int numstep = (int) abs(ceil((double)args->bpp_arg /(16)));
  ulong *sortedrs[2];
  ulong *histogramrs[2];
  ulong *histogram;
  histogramrs[0] = (ulong *) malloc(NUMBUCKETS*sizeof(ulong));
  histogramrs[1] = (ulong *) malloc(NUMBUCKETS*sizeof(ulong));
  sortedrs[0] = (ulong *) malloc(tree->size_curr * sizeof(ulong));
  sortedrs[1] = (ulong *) malloc(tree->size_curr * sizeof(ulong));

  ulong *loc_hist = calloc(nthreads*NUMBUCKETS, sizeof(ulong)); check_alloc(loc_hist, 610);

  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();		    /* Thread number */
    ulong 	lwb     = (id*tree->size_curr)/nthreads;    /* Lower bound for current thread */
    ulong 	upb     = ((id+1)*tree->size_curr)/nthreads;     /* Upper bound for current thread */
    for(int step=0; step<numstep; step++) {
	 
      local_hist_rs(tree->gval, lwb, upb, &loc_hist[NUMBUCKETS*id], sortedrs[step%2], step);
      #pragma omp barrier

      if(id == 0){
	histogram = histogramrs[step%2];
	memset(histogram, 0, sizeof(ulong)*NUMBUCKETS);
	ulong prevhist, total= 0;
	for(int j=0; j<NUMBUCKETS; j++){
	  for(int i=0; i<nthreads; i++){	  
	    histogram[j] += loc_hist[NUMBUCKETS*i+j];	  
	  }
	  prevhist = histogram[j];
	  histogram[j] = total;
	  total += prevhist;
	}
	
	for(int j=0; j<NUMBUCKETS; j++){
	  ulong sum = loc_hist[j], curr;
	  for(int i=1; i<nthreads; i++){	  
	    curr = loc_hist[NUMBUCKETS*i+j];
	    loc_hist[NUMBUCKETS*i+j] = sum + histogram[j];
	    sum += curr;
	  }
	  loc_hist[j] = histogram[j];
	}
      }
      #pragma omp barrier
      create_sorted_ars(tree->gval, lwb, upb, &loc_hist[NUMBUCKETS*id], sortedrs[((step+1)%2)], sortedrs[step%2], step);
    }
  }
  // for( ulong i = 0; i < tree->size; i++) 
  //   info("gval%d  sorted %d",tree->gval[i], sortedrs[numstep%2][i]);
  free(sortedrs[(numstep+1)%2]);
  free(histogramrs[(numstep)%2]);
  free(loc_hist);
  return sortedrs[numstep%2];
  //return ranks_inv;
}



/*void build_local_tree2(Arguments *args, Node *local_tree, ulong *dims){
  int connectivity =  args->connectivity_arg;
  int nthreads     =  args->threads_arg;
  int attrib 	   =  args->attribute_arg;
  ulong inc        =  dims[2] == 1 ? dims[0]: dims[0]*dims[1];

  omp_lock_t lock[nthreads];
  for (int i=0; i<nthreads; i++)
    omp_init_lock(&(lock[i]));
  int *saval = calloc(nthreads, sizeof(int));

  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			
    ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads) 
      : dims[0]*dims[1]*((id*dims[2])/nthreads) ;              
    ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/nthreads)
      : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);  

    Allocator work_alloc;
    ulong size_alloc = workspace_size(args, upb-lwb, local_tree->size_att);
    allocator_init(&work_alloc, size_alloc);
    ulong *ranks = NULL;
    ulong *ranks_inv  = NULL;
    create_mappings(local_tree->gval, &ranks, &ranks_inv, lwb, upb);

    PrioQueue *q        = create_prio_queue(upb-lwb-1, &work_alloc);
    BitArray *visited   = create_bit_array(upb-lwb, &work_alloc);
    info("Mappings: wallclock time = %0.2f ",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
    
    init_attrib_array(local_tree->attribute, local_tree->border, dims, local_tree->attr_off, lwb, upb, local_tree->size_att);


    tree_flood_tee_base(local_tree->parent, local_tree->attribute, visited, q, ranks, ranks_inv, dims, lwb, upb, connectivity, local_tree->size_att);

    printf("Build: wallclock time = %0.2f \n",
	   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));

    for (ulong x =lwb ; x < upb; x++){
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
      }
      #pragma omp barrier


    if (nthreads > 1) {
      int i = 1;
      int q = id;
      while (id+i < nthreads && q%2 == 0){
		
	while (saval[id+i] <= 0){
	  #pragma omp flush // wait
	}
	omp_set_lock(&lock[id+i]);
	saval[id+i]--;
	omp_unset_lock(&lock[id+i]);
	fuse_sections(local_tree, dims, id, i, nthreads, connectivity);
	i *= 2;
	q /= 2;      
      }
      if(id>0){
	omp_set_lock(&lock[id]);	      
	saval[id]++;
	omp_unset_lock(&lock[id]);
      }

    }
    
  }

  free(saval);
  for (int i = 0; i < nthreads; i++)
    omp_destroy_lock(&(lock[i]));

} /* build_local_tree2 */
