#include "types.h"
#include "attributes.h"
#include "tree_flood.h"
#include "queue.h"
#include "mpihelper.h"

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Tree Building         */
/*				   */
/* +++++++++++++++++++++++++++++++ */

idx *build_local_tree(Arguments *args, Node *tree, ulong *dims, ulong *attr_off){
  int flood_algo   =  args->flood_arg;
  int connectivity =  args->connectivity_arg;
  int nthreads     =  args->threads_arg;
  int attrib 	   =  args->attribute_arg;

  omp_lock_t lock[nthreads];
  for (int i=0; i<nthreads; i++)
    omp_init_lock(&(lock[i]));
  int *saval = calloc(nthreads, sizeof(int));
  
  
  if(args->bpp_arg < 0 || args->bpp_arg >=16)
    nthreads = 1;
  //AuxDataStore **store =  malloc(nthreads * sizeof(AuxDataStore*)); check_alloc(store, 500);
  tree->parent = malloc(tree->size * sizeof(idx)); 

  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			/* Thread number */
    ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads)
      : dims[0]*dims[1]*((id*dims[2])/nthreads);              /* Lower bound for current thread */
    ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/nthreads)
      : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);;       /* Upper bound for current thread */

    
    tree->store[id] = malloc(sizeof(AuxDataStore)); //check_alloc(store_th, 500);
    AuxDataStore *store_th =  tree->store[id];
    init_aux_data_store(store_th, AttribsArray[attrib].size, upb-lwb);
    //void *store_th = NULL;
    if(flood_algo == 0){
      /* Salembier algorithm */
      
      ulong 	min_idx     = lwb;
      void 	*curr_attr  = NULL;	
      ulong 	*histogram  = calloc(g_max_levels, sizeof(ulong)); check_alloc(histogram,  501);
      idx	*levelroots = calloc(g_max_levels, sizeof(idx));   check_alloc(levelroots, 502);
      bool 	*reached    = calloc(upb - lwb,    sizeof(bool));  check_alloc(reached, 503);
      Queue 	*queue      = create_queue(upb - lwb, g_max_levels);

      for (ulong i = 0; i < g_max_levels; i++) 
	levelroots[i] = BOTTOM;

      for (ulong i = lwb; i < upb; i++) {
	histogram[(value_t) tree->gval[i]]++; 
	if (tree->gval[min_idx] > tree->gval[i])
	  min_idx = i;
      }
      
      set_queue_offsets(queue, histogram, g_max_levels);
      free(histogram);

      queue_add(queue, (value_t) tree->gval[min_idx], min_idx);
      levelroots[(value_t) tree->gval[min_idx]] = min_idx;
      reached[min_idx-lwb]   = true;
      tree_flood_sal(tree, store_th, queue, levelroots, reached, dims, attr_off, lwb, upb, (value_t) tree->gval[min_idx], connectivity, &curr_attr);
      free(queue);
      free(levelroots);
      free(reached);
      
    } else if (flood_algo == 2) {
      /* Wilkinson algorithm */

      ulong 	min_idx = lwb;
      
      for (ulong i = lwb; i < upb; i++) {
	tree->parent[i] = BOTTOM;
	if (tree->gval[min_idx] > tree->gval[i])
	  min_idx = i;
      }
    
      pQueue 	*queue = pQueueCreate(upb-lwb);
      pStack 	*stack = pStackCreate(upb-lwb, queue->array);
    
      tree_flood_wil(tree, store_th, queue, stack, dims, attr_off, lwb, upb, min_idx, connectivity);
    
      pQueueDelete(queue);
      pStackDelete(stack);
      
    } else if (flood_algo == 1) {
	/* Wilkinson improved algorithm */ 
      ulong *ranks        = calloc(upb-lwb, sizeof(ulong));  check_alloc(ranks, 504);
      ulong *ranks_inv    = calloc(upb-lwb, sizeof(ulong));  check_alloc(ranks_inv, 055);
      bool *visited   	  = calloc(upb-lwb, sizeof(bool));  check_alloc(visited, 506);

      create_mappings(tree->gval, ranks, ranks_inv, upb-lwb, lwb, upb);

      PrioQueue *q        = create_prio_queue(upb-lwb);	    
      tree_flood_tee(tree,  store_th, visited, q, ranks, ranks_inv, dims, attr_off, lwb, upb, connectivity, 0);
      free_prio_queue(q);
     
      free(visited);
      free(ranks);
      free(ranks_inv);
      
    }
    if (nthreads > 1) {
      /* Merge threads sections */
      int i = 1;
      int q = id;
      while (id+i < nthreads && q%2 == 0){
		
	while (saval[id+i] <= 0){
	  #pragma omp flush // wait
	}
	omp_set_lock(&lock[id+i]);
	saval[id+i]--;
	omp_unset_lock(&lock[id+i]);

	fuse_sections(tree, store_th, dims, id, i, nthreads, connectivity);
	
	i *= 2;
	q /= 2;      
      }
      if(id>0){
	omp_set_lock(&lock[id]);	      
	saval[id]++;
	omp_unset_lock(&lock[id]);
      }

    }
    realloc_store(tree->store[id]);
  }
  free(saval);
  for (int i = 0; i < nthreads; i++)
    omp_destroy_lock(&(lock[i]));

  return tree->parent;
} /* build_local_tree */


  /* +++++++++++++++++++++++++++++++ */
  /*           Tree Salembier        */
  /* +++++++++++++++++++++++++++++++ */

long tree_flood_sal(Node *tree, AuxDataStore *store, Queue *q,  idx *levelroot, bool *reached, ulong *dims, ulong *attr_off, ulong lwb, ulong upb, long level, int connectivity, void **thisattr) {

  void  *curr_attr  = NULL;
  void  *attr_child = NULL;
  int 	n_neighbors;
  long  fc;
  ulong neighbors[connectivity]; 
  ulong c, p, x, y, z;
  double *init_attr = calloc(4, sizeof(double));
  
  while (queue_is_not_empty(q, level)) {

    p = queue_first(q, level); 
    x = p % (dims[0] * dims[1]) % dims[0] ;
    y = p % (dims[0] * dims[1]) / dims[0];
    z = p / (dims[0] * dims[1]) ;

    init_attr[0] = (double) (x + attr_off[0]);
    init_attr[1] = (double) (y + attr_off[1]);
    init_attr[2] = (double) (z + attr_off[2]);
    init_attr[3] = (double) 1;

    n_neighbors = get_neighbors(dims, lwb, upb, neighbors, p, x, y, z, connectivity);
    
    if (curr_attr){
      if(!is_border(tree->border, dims, p))
	add_to_aux_data(curr_attr, init_attr);    
    } else{
      if(!is_border(tree->border, dims, p))
	curr_attr = new_aux_data(store, init_attr);    
      if (*thisattr){
	if(curr_attr){
	  merge_aux_data(curr_attr, *thisattr);
	}					       
	else{
	  clone_aux_data(store, &curr_attr, *thisattr);
	}
      }
    }
         
    for (int i = 0; i < n_neighbors; i++) {
      c = neighbors[i]; 
      if (!reached[c-lwb]) {
        reached[c-lwb] = true;
        fc = (long) tree->gval[c]; 
        if (levelroot[fc] == BOTTOM) { 
          levelroot[fc] = c; 
        } else {
          tree->parent[c] = levelroot[fc];
        }

        queue_add(q, fc, c);

        if (fc > level) { 
	  attr_child = NULL;
          do {
            fc = tree_flood_sal(tree, store, q, levelroot, reached, dims, attr_off, lwb, upb, fc, connectivity, &attr_child);
            if ((ulong) fc >= g_max_levels) { 
	      if(curr_attr)
		delete_aux_data(curr_attr);
              return fc;
            }
          } while (fc != level);
	  if(attr_child != NULL){
	    if(curr_attr) merge_aux_data(curr_attr, attr_child);
	    else {
	      //    #pragma omp critical
	      clone_aux_data(store, &curr_attr, attr_child);
	    }
	  }	
	}
      }
    }
  }

  long m = (long) level - 1;
  while (m > 0 && levelroot[m] == BOTTOM) m--;
  if (m >= 0) tree->parent[levelroot[level]] = levelroot[m];
  else tree->parent[levelroot[level]] = BOTTOM;
  if(curr_attr != NULL) tree->attribute[levelroot[level]] = curr_attr;
  levelroot[level] = BOTTOM;
  *thisattr = curr_attr;
  free(init_attr);
  return m;
} /* tree_flood_sal */


/* +++++++++++++++++++++++++++++++ */
/*     	     Tree Wilkinson        */
/* +++++++++++++++++++++++++++++++ */

void tree_flood_wil(Node *tree, AuxDataStore *store,  pQueue *queue, pStack *stack, ulong *dims, ulong *attr_off, ulong lwb, ulong upb, ulong min_idx, int connectivity){
  void *curr_attr = NULL;
  int n_neighbors;
  ulong neighbors[connectivity];
  ulong nextpix, p, q, x, y, z, oldtop;
  double *init_attr = calloc(4, sizeof(double));
  pstack_push(stack, min_idx);


  if(!is_border(tree->border, dims, min_idx)){
    x = min_idx % (dims[0] * dims[1]) % dims[0];
    y = min_idx % (dims[0] * dims[1]) / dims[0];
    z = min_idx / (dims[0] * dims[1]);
    init_attr[0] = (double) (x + attr_off[0]);
    init_attr[1] = (double) (y + attr_off[1]);
    init_attr[2] = (double) (z + attr_off[2]);
    init_attr[3] = (double) tree->gval[min_idx];
    curr_attr = new_aux_data(store, init_attr);
    tree->attribute[min_idx] = curr_attr;
  }
  
  tree->parent[min_idx]    = min_idx;
  nextpix = min_idx;

  do {
    p 		= nextpix;
    x 		= p % (dims[0] * dims[1]) % dims[0];
    y 		= p % (dims[0] * dims[1]) / dims[0];
    z		= p / (dims[0] * dims[1]);
    n_neighbors = get_neighbors(dims, lwb, upb, neighbors, p, x, y, z, connectivity);
    
    for (int i=0; i < n_neighbors; i++) {
      q = neighbors[i];
      if (tree->parent[q] == BOTTOM) {
	tree->parent[q] = q;

	if(!is_border(tree->border, dims, q)){
	  x = q % (dims[0] * dims[1]) % dims[0];
	  y = q % (dims[0] * dims[1]) / dims[0];
	  z = q / (dims[0] * dims[1]);
	  init_attr[0] = (double) (x + attr_off[0]);
	  init_attr[1] = (double) (y + attr_off[1]);
	  init_attr[2] = (double) (z + attr_off[2]);
	  init_attr[3] = (double) 1;
	  curr_attr = new_aux_data(store, init_attr);
	  tree->attribute[q] = curr_attr;
	}

        if (tree->gval[q] > tree->gval[p]) {
	  pstack_push(stack,q);
	  nextpix = q;
	  break;
	}
	pQueuePush(tree->gval, queue,q);
      }
    }

    if (nextpix == p) {     
      if (p != pstack_top(stack)) {    

	p =  pQueuePop(tree->gval, queue); 
	tree->parent[p] = pstack_top(stack);


	if(!is_border(tree->border, dims, p)){
	  x = p % (dims[0] * dims[1]) % dims[0];
	  y = p % (dims[0] * dims[1]) / dims[0];
	  z = p / (dims[0] * dims[1]);
	  init_attr[0] = (double) (x + attr_off[0]);
	  init_attr[1] = (double) (y + attr_off[1]);
	  init_attr[2] = (double) (z + attr_off[2]);
	  init_attr[3] = (double) 1;
	  if (tree->attribute[pstack_top(stack)])
	    add_to_aux_data(tree->attribute[pstack_top(stack)], init_attr);
	  else {
	    //    #pragma omp critical
	    curr_attr = new_aux_data(store, init_attr);
	    tree->attribute[pstack_top(stack)] = curr_attr;
	  }
	}

	if (!is_empty_pqueue(queue)){
	  nextpix = pqueue_front(queue);

	  if (tree->gval[nextpix] < tree->gval[p])  
	    nextpix = pstack_top(stack);
	  
	} else {
	    nextpix = pstack_top(stack);
	}
   
      } else {                        
	if (!is_empty_pqueue(queue)) {	  
	  nextpix = pqueue_front(queue);    
	  if (tree->gval[nextpix] < tree->gval[p]) { 
	    oldtop = pstack_pop(stack);
	    p = pstack_top(stack);
	    if (tree->gval[nextpix] > tree->gval[p]) {
	      nextpix = pQueuePop(tree->gval, queue); 
	      pstack_push(stack, nextpix);
              p = nextpix;         
	    } else if (tree->gval[nextpix] < tree->gval[p]) {
              nextpix = p;              
	    }
	    
	    tree->parent[oldtop] = p;        
	    
	    if(!tree->attribute[p] && !is_border(tree->border, dims, p)){
	      x = p % (dims[0] * dims[1]) % dims[0];
	      y = p % (dims[0] * dims[1]) / dims[0];
	      z = p / (dims[0] * dims[1]);
	      init_attr[0] = (double) (x + attr_off[0]);
	      init_attr[1] = (double) (y + attr_off[1]);
	      init_attr[2] = (double) (z + attr_off[2]);
	      init_attr[3] = (double) 1;
	      //   #pragma omp critical
	      curr_attr = new_aux_data(store, init_attr);
	      tree->attribute[p] = curr_attr;
	    }
	      
	    if(tree->attribute[oldtop]) {             
	       if(!is_border(tree->border, dims, p))
		 merge_aux_data(tree->attribute[p], tree->attribute[oldtop]);
	       else if(tree->attribute[p])
		 merge_aux_data(tree->attribute[p], tree->attribute[oldtop]);
	       else{
		 //	 #pragma omp critical
		 clone_aux_data(store, &tree->attribute[p], tree->attribute[oldtop]);
	       }
	    }
	    
	  }
	} else {
	  oldtop = pstack_pop(stack);
	  if (!is_empty_pstack(stack)) {
	    p = pstack_top(stack);
	    tree->parent[oldtop] = p;      	    
	    if(!tree->attribute[p] && !is_border(tree->border, dims, p)){
	      x = p % (dims[0] * dims[1]) % dims[0];
	      y = p % (dims[0] * dims[1]) / dims[0];
	      z = p / (dims[0] * dims[1]);
	      init_attr[0] = (double) (x + attr_off[0]);
	      init_attr[1] = (double) (y + attr_off[1]);
	      init_attr[2] = (double) (z + attr_off[2]);
	      init_attr[3] = (double) 1;
	      curr_attr = new_aux_data(store, init_attr);
	      tree->attribute[p] = curr_attr;
	    }
 
	      
	    if(tree->attribute[oldtop]) {             
	      if(!is_border(tree->border, dims, p))
		merge_aux_data(tree->attribute[p], tree->attribute[oldtop]);
	      else if(tree->attribute[p])
		merge_aux_data(tree->attribute[p], tree->attribute[oldtop]);
	      else{
		//	 #pragma omp critical
		clone_aux_data(store, &tree->attribute[p], tree->attribute[oldtop]);
	      }
	    }
	    
	    nextpix = p;
	  }
	}	
      }
    }
  } while (!is_empty_pqueue(queue) || (!is_empty_pstack(stack)));   
  tree->parent[min_idx] = BOTTOM;
} /* tree_flood_wil */


/* +++++++++++++++++++++++++++++++ */
/*     	     Tree Teeninga         */
/* +++++++++++++++++++++++++++++++ */

void tree_flood_tee(Node *tree, AuxDataStore *store,  bool *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims, ulong *attr_off, ulong lwb, ulong upb, int connectivity, int periodic){

  ulong index = lwb;

  ulong rank_cur = ranks[index-lwb];
  visited[index-lwb] = true;
  double *init_attr = calloc(4, sizeof(double));
  while(true){

    remaining(visited, q, ranks, &index, &rank_cur, dims, lwb, upb, connectivity, periodic);
    if (q->m_levels[0][0] == 0) break;
    rank_cur = q->m_top;
    prio_queue_remove(q);
    ulong parent = ranks_inv[rank_cur]+lwb;

    tree->parent[index] = parent;

    if(!tree->attribute[index] && !is_border(tree->border, dims, index)){

      init_attr[0] = (double) ((index % (dims[0] * dims[1])) % dims[0]  + attr_off[0]);
      init_attr[1] = (double) ((index % (dims[0] * dims[1])) / dims[0]  + attr_off[1]);
      init_attr[2] = (double) (index / (dims[0] * dims[1]) + attr_off[2]);
      init_attr[3] = (double) 1;//tree->gval[index];
      // init_attr[4] = tree->dens ? (double) tree->gval_dens[index]: 0;

      tree->attribute[index] = new_aux_data (store, init_attr);
    }
    
    if(!tree->attribute[parent] && !is_border(tree->border, dims, parent)){
      init_attr[0] = (double) ((parent % (dims[0] * dims[1])) % dims[0]  + attr_off[0]);
      init_attr[1] = (double) ((parent % (dims[0] * dims[1])) / dims[0]  + attr_off[1]);
      init_attr[2] = (double) (parent / (dims[0] * dims[1]) + attr_off[2]);
      init_attr[3] = (double) 1;//tree->gval[parent];
      //init_attr[4] = tree->dens ? (double) tree->gval_dens[parent]: 0;

      tree->attribute[parent] = new_aux_data (store, init_attr);
    }

  
    if(tree->attribute[index]){
      if(!is_border(tree->border, dims, parent) || tree->attribute[parent])
	merge_aux_data(tree->attribute[parent], tree->attribute[index]);
      else{
	//	#pragma omp critical
	clone_aux_data(store, &tree->attribute[parent], tree->attribute[index]);
      }      
    }
    index = parent;
  }
  tree->parent[index] = BOTTOM;
  free(init_attr);
} /* tree_flood_tee */


/* ++++++++++++++++++++++++++ */
/*      Threads functions     */
/* ++++++++++++++++++++++++++ */

void fuse_periodic_3D(Node *tree, AuxDataStore *store, ulong *dims, int connectivity){
  
  value min_curr, min_prev;
  bool test_min = 0;
  ulong u_x, u_y, u_z, v, vx, vy, vz;
  ulong p ;

  
  for (ulong u = 0, p = 0; p < dims[1]*dims[2]; u+=dims[0], p++){
    u_x = (u % (dims[0]*dims[1])) % dims[0];
    u_y = (u % (dims[0]*dims[1])) / dims[0];
    u_z = (u / (dims[0]*dims[1]));
    v   = u+dims[0]-1;
    vx  = (v % (dims[0]*dims[1])) % dims[0];
    vy  = (v % (dims[0]*dims[1])) / dims[0];
    vz  = (v / (dims[0]*dims[1])); 
    ulong vybef = vy > 0? vy-1:  dims[1]-1;
    ulong vyaf  = vy < dims[1]-1? vy+1: 0;
    ulong vzbef = vz > 0? vz-1:  dims[2]-1;
    ulong vzaf  = vz < dims[2]-1? vz+1: 0;


    
    if(u == 0) test_min = false;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (!test_min || min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);
    min_prev = min_curr;
    test_min = true;

    v = vzaf*(dims[0]*dims[1]) + vy*dims[0] + vx;
    //  attr[2] = (long)(u_z + 1) <  dims[2]? 0 : dims[2];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzaf*(dims[0]*dims[1]) + vyaf*dims[0] + vx;
    //  attr[1] =(long) (u_y + 1) <  dims[1]? 0 : dims[1];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzaf*(dims[0]*dims[1]) + vybef*dims[0] + vx;
    // attr[1] = (long)(u_y - 1) >= 0? 0 : dims[1];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vyaf*dims[0] + vx;
    //  attr[2]=0;
    //  attr[1] = (long)(u_y + 1) <  dims[1]? 0 : dims[1];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vybef*dims[0] + vx;
    //  attr[1] =(long) (u_y - 1) >= 0? 0 : dims[1];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzbef*(dims[0]*dims[1]) + vy*dims[0] + vx;
    // attr[2] = (long) (u_z - 1) >= 0 ? 0 : dims[2];
    //  attr[1] = 0;
    
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzbef*(dims[0]*dims[1]) + vyaf*dims[0] + vx;
    //  attr[1] = (long) (u_y + 1) <  dims[1]? 0 : dims[1];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzbef*(dims[0]*dims[1]) + vybef*dims[0] + vx;
    //  attr[1] = (long) (u_y - 1) >= 0? 0 : dims[1];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

  }

  //  attr[1] = dims[1];
  // attr[0] = 0;
  // attr[2] = 0;
  // Top/Bottom
  
   for (ulong u = 0, p =0; p < dims[0]*dims[2]; p++, u++){
    if(p % dims[0] == 0 && p != 0)
      u += dims[1]*(dims[0]-1);
    u_x = (u % (dims[0]*dims[1])) % dims[0];
    u_y = (u % (dims[0]*dims[1])) / dims[0];
    u_z = (u / (dims[0]*dims[1]));
    v   = u+(dims[0]-1)*dims[1];
    vx  = (v % (dims[0]*dims[1])) % dims[0];
    vy  = (v % (dims[0]*dims[1])) / dims[0];
    vz  = (v / (dims[0]*dims[1])); 
    ulong vxbef = vx > 0? vx-1:  dims[0]-1;
    ulong vxaf  = vx <  dims[0]-1? vx+1: 0;
    ulong vzbef = vz > 0? vz-1:  dims[2]-1;
    ulong vzaf  = vz < dims[2]-1? vz+1: 0;

    if(u == 0) test_min = false;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (!test_min || min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);
    min_prev = min_curr;
    test_min = true;

    v = vzaf*(dims[0]*dims[1]) + vy*dims[0] + vx;
    //  attr[2] = (long) (u_z + 1) <  dims[2]? 0 : dims[2];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzaf*(dims[0]*dims[1]) + vy*dims[0] + vxaf;
    //  attr[0] = (long) (u_x + 1) <  dims[0]? 0 : dims[0];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzaf*(dims[0]*dims[1]) + vy*dims[0] + vxbef;
    // attr[0] = (long) (u_x - 1) >= 0? 0 : dims[0];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vy*dims[0] + vxaf;
    //   attr[2] = 0;
    //   attr[0] =(long)  (u_x + 1) <  dims[0]? 0 : dims[0];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vy*dims[0] + vxbef;
    //  attr[0] = (long) (u_x - 1) >= 0? 0 : dims[0];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzbef*(dims[0]*dims[1]) + vy*dims[0] + vx;
    //  attr[2] = (long) (u_z - 1) >= 0 ? 0 : dims[2];
    //  
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzbef*(dims[0]*dims[1]) + vy*dims[0] + vxaf;
    //  attr[0] = (long) (u_x + 1) <  dims[0]? 0 : dims[0];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vzbef*(dims[0]*dims[1]) + vy*dims[0] + vxbef;
    //  attr[0] = (long) (u_x - 1) >= 0? 0 : dims[0];

    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

   }

  // Front/Back
    for (ulong u = 0; u < dims[0]*dims[1]; u++){
    u_x = (u % (dims[0]*dims[1])) % dims[0];
    u_y = (u % (dims[0]*dims[1])) / dims[0];
    u_z = (u / (dims[0]*dims[1]));
    v   = u+dims[0]*dims[1]*(dims[2]-1);
    vx  = (v % (dims[0]*dims[1])) % dims[0];
    vy  = (v % (dims[0]*dims[1])) / dims[0];
    vz  = (v / (dims[0]*dims[1])); 
    ulong vxbef = vx > 0? vx-1:  dims[0]-1;
    ulong vxaf  = vx <  dims[0]-1? vx+1: 0;
    ulong vybef = vy > 0? vy-1:  dims[1]-1;
    ulong vyaf  = vy < dims[1]-1? vy+1: 0;

    //   info("\n%ld %ld, %ld",u_x, u_y, u_z);

    if(u == 0) test_min = false;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (!test_min || min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);
    min_prev = min_curr;
    test_min = true;

    v = vz*(dims[0]*dims[1]) + vyaf*dims[0] + vx;
    // attr[1] = (long) (u_y + 1) <  dims[1]? 0 : dims[1];
    //   info("%ld %ld, %ld",vx, vy, vz);
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);
 
    v = vz*(dims[0]*dims[1]) + vyaf*dims[0] + vxaf;
    // attr[0] = (long) (u_x + 1) <  dims[0]? 0 : dims[0];
    //   info("%ld %ld, %ld",u_x, u_y, u_z);
    //   info("%ld %ld, %ld",vx, vy, vz);
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vyaf*dims[0] + vxbef;
    // attr[0] = (long) (u_x - 1) >= 0? 0 : dims[0];
    //   info("%ld %ld, %ld",u_x, u_y, u_z);
    //  info("%ld %ld, %ld",vx, vy, vz);
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vy*dims[0] + vxaf;
    // attr[0] = (long) (u_x + 1) <  dims[0]? 0 : dims[0];
    // attr[1] = 0;
    //   info("%ld %ld, %ld",u_x, u_y, u_z);
    //  info("%ld %ld, %ld",vx, vy, vz);
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vy*dims[0] + vxbef;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vybef*dims[0] + vx;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vybef*dims[0] + vxaf;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);

    v = vz*(dims[0]*dims[1]) + vybef*dims[0] + vxbef;
    min_curr = MIN(tree->gval[u], tree->gval[v]);
    if (min_curr > min_prev)
      merge_nodes_perio(tree, store,dims, u, v);
      }
  }

void merge_nodes_perio(Node *tree, AuxDataStore *store, ulong *dims, idx x, idx y) {
  void *cor  = NULL;
  void *copa = NULL;
  idx h, z;
  ulong count = 0;


  x =  get_levelroot(tree, x);
  y =  get_levelroot(tree, y);

  if (tree->gval[x] < tree->gval[y]) {
    h=x; x=y; y=h;
  }


  
  while ((x != y) && (y != BOTTOM)) {
    z = get_parent(tree, x);

    if ((z != BOTTOM) && (tree->gval[z]>=tree->gval[y])) {
      if (cor) {
	if(tree->attribute[x])
	  correct_inertiaeor_data(tree->attribute[x], cor);
	else {
	  //	  #pragma omp critical
	  clone_aux_data(store, &(tree->attribute[x]), cor);
	}
      }
      x = z;
    } else {
      if (cor && tree->attribute[x]) merge_to_aux_data(NULL, &copa, tree->attribute[x], cor);
      else if(tree->attribute[x]){
	clone_aux_data(NULL, &copa, tree->attribute[x]);
      }
      if(tree->attribute[x]){
	clone_aux_data(NULL, &cor, tree->attribute[x]);
      } else {
      	delete_aux_data(cor);
	cor= NULL;
      }
      
      if(copa){
	//	#pragma omp critical
      	clone_aux_data(store, &(tree->attribute[x]), copa);
      }
      tree->parent[x] = y;
      x = y;
      y = z;
    }
  }
  if (y == BOTTOM && cor) {
    while(x != BOTTOM) {
      if (tree->attribute[x])
	correct_inertiaeor_data(tree->attribute[x], cor);
      else{
	//	#pragma omp critical
	clone_aux_data(store, &tree->attribute[x], cor);
      }
      x = get_parent(tree, x);
    }
  }
  if (cor)  delete_aux_data(cor);
  if (copa) delete_aux_data(copa);
  
  
}/* merge_nodes */

void fuse_sections(Node *tree, AuxDataStore *store, ulong *dims, uint id, uint i, uint n_t, int connectivity){
  ulong  mdb;
  ulong p, u, v, u_x, u_y;
  value min_curr, min_prev;
  bool test_min = 0;

  mdb = dims[2] == 1 ? dims[0]*(((id+i) * dims[1]) / n_t) : dims[0]*dims[1]*(((id+i) * dims[2]) / n_t);

  if(dims[2] == 1){
    for (p=0, u = mdb; p < dims[0]; p++, u++){
      u_x = u % dims[0];
      v = u-dims[0];
     if(u_x == 0) test_min = false;
     min_curr = MIN(tree->gval[u], tree->gval[v]);
     if (!test_min || min_curr > min_prev)
       merge_nodes(tree, store, u, v);
     min_prev = min_curr;
     test_min = true;
     if(connectivity >= 8){
       if(u_x > 0){
	 min_curr = MIN(tree->gval[u], tree->gval[v-1]);
	 if (min_curr > min_prev)
	   merge_nodes(tree, store, u, v-1);
       }     
       if(u_x < dims[0]-1){
	 min_curr = MIN(tree->gval[u], tree->gval[v+1]);
	 if (!test_min || min_curr > min_prev)
	   merge_nodes(tree, store, u, v+1);
       }
     }
    }
  } else {
    for (p=0, u = mdb; p < (dims[0] * dims[1]); p++, u++){
      u_x = (u % (dims[0]*dims[1])) % dims[0];
      u_y = (u % (dims[0]*dims[1])) / dims[0];
      v   = u - dims[0]*dims[1];

      if(u_x == 0) test_min = false;
      min_curr = MIN(tree->gval[u], tree->gval[v]);
      if (!test_min || min_curr > min_prev )
	merge_nodes(tree, store, u, v);
      min_prev = min_curr;
      test_min = true;
		
      if(connectivity == 26){
	if(u_x > 0 && u_y > 0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]-1]);
	  if (min_curr > min_prev )
	    merge_nodes(tree, store, u, v-dims[0]-1);
	}
	if(u_y > 0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]]);
	  if (min_curr > min_prev )
	    merge_nodes(tree, store, u, v-dims[0]);
	}
	if(u_y > 0 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]+1]);
	  if (min_curr > min_prev )
	    merge_nodes(tree, store, u, v-dims[0]+1);
	}
	if(u_x>0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-1]);
	  if (min_curr > min_prev )
	    merge_nodes(tree, store, u, v-1);	      
	}
	if(u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes(tree, store, u, v+1);
	}
	if(u_y < dims[1]-1 && u_x>0){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]-1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes(tree, store, u, v+dims[0]-1);
	}
	if(u_y < dims[1]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes(tree, store, u, v+dims[0]);
	}
	if(u_y < dims[1]-1 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes(tree, store, u, v+dims[0]+1);
	}
      }
    }
  }
}



void merge_nodes(Node *tree, AuxDataStore *store, idx x, idx y) {
  void *cor  = NULL;
  void *copa = NULL;
  idx h, z;

  x =  get_levelroot(tree, x);
  y =  get_levelroot(tree, y);

  if (tree->gval[x] < tree->gval[y]) {
    h=x; x=y; y=h;
  }
  while ((x != y) && (y != BOTTOM)) {
    z = get_parent(tree, x);

    if ((z != BOTTOM) && (tree->gval[z]>=tree->gval[y])) {
      if (cor) {
	if(tree->attribute[x])
	  merge_aux_data(tree->attribute[x], cor);
	else {
	  //	  #pragma omp critical
	  clone_aux_data(store, &(tree->attribute[x]), cor);
	}
      }
      x = z;
    } else {
      if (cor && tree->attribute[x]) merge_to_aux_data(NULL, &copa, tree->attribute[x], cor);
      else if(tree->attribute[x]){
	clone_aux_data(NULL, &copa, tree->attribute[x]);
      }
      if(tree->attribute[x]){
	clone_aux_data(NULL, &cor, tree->attribute[x]);
      } else {
      	delete_aux_data(cor);
	cor= NULL;
      }
      
      if(copa){
	//	#pragma omp critical
      	clone_aux_data(store, &(tree->attribute[x]), copa);
      }
      tree->parent[x] = y;
      x = y;
      y = z;
    }
  }
  if (y == BOTTOM && cor) {
    while(x != BOTTOM) {
      if (tree->attribute[x])
	merge_aux_data(tree->attribute[x], cor);
      else{
//	#pragma omp critical
	clone_aux_data(store, &tree->attribute[x], cor);
      }
      x = get_parent(tree, x);
    }
  }
  if (cor)  delete_aux_data(cor);
  if (copa) delete_aux_data(copa);
}/* merge_nodes */


/* +++++++++++++++++++++++++++++++ */
/*				   */
/*         Annexes Functions       */
/*				   */
/* +++++++++++++++++++++++++++++++ */



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
  
  if (dims[2] > 1) {
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


int bit_scan_reverse(ulong val) {
  return sizeof(val) * CHAR_BIT - 1 - __builtin_clzl(val);
}
int bits_per_word_log2(void) {
  return sizeof(unsigned) * CHAR_BIT - 1 - __builtin_clz((sizeof(ulong) * CHAR_BIT));
}
int bits_per_word(void) { 
  return sizeof(ulong) * CHAR_BIT;
}


void create_mappings(value *gvals, ulong *ranks, ulong *ranks_inv, ulong size, ulong lwb, ulong upb) {
  ulong histos[NUM_DIGITS][MXT_HISTO_SZ];
  gen_histogram(gvals, histos, lwb, upb);

  for (int digit_nr = 0; digit_nr != NUM_DIGITS; digit_nr++) 
    exclusive_sum(histos[digit_nr], histos[digit_nr] + MXT_HISTO_SZ);
  create_ranks_inv(gvals, ranks_inv, histos, size, lwb, upb);

  if(ranks != NULL)
    for (ulong i = 0; i != upb-lwb; ++i){
      ranks[ranks_inv[i]] = i;
    }
}



void gen_histogram(value *gvals, ulong histos[NUM_DIGITS][MXT_HISTO_SZ],  ulong lwb, ulong upb) {
  memset(histos, 0, NUM_DIGITS * MXT_HISTO_SZ * sizeof(histos[0][0]));

  for (ulong i = lwb; i != upb; ++i) {
    int shift = 0;
    value_t gval;
    if(FLOAT_TYPE == 1)
      gval = transform((gvals[i])) ;
    else
      gval = gvals[i];

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


void create_ranks_inv(value *gvals, ulong *ranks_inv, ulong histos[NUM_DIGITS][MXT_HISTO_SZ], ulong size, ulong lwb, ulong upb) {    
  SortItem **pairs = malloc(2 * sizeof(SortItem*));
  pairs[0] = calloc(size, sizeof(SortItem));
  pairs[1] = calloc(size, sizeof(SortItem));

  scatter_first_digit(gvals, pairs[1], histos[0], lwb, upb);
 
  int digit_nr = 1;
  for (; digit_nr != NUM_DIGITS - 1; ++digit_nr)
  {
    int mod1 = (digit_nr) % 2;
    int mod2 = (digit_nr+1) % 2;
    scatter_digit(pairs[mod1], pairs[mod2], digit_nr, histos[digit_nr], size);
  }

  scatter_last_digit(pairs[1], ranks_inv, digit_nr, histos[digit_nr], size);
  free(pairs[0]);
  free(pairs[1]);
  free(pairs);
}

void scatter_first_digit(value *gvals, SortItem *pair_it, ulong* histo, ulong lwb, ulong upb){
  for (ulong i = lwb; i != upb; ++i) {
    value_t gval;
    if(FLOAT_TYPE == 1)
      gval = transform((gvals[i])) ;
    else
      gval = gvals[i];
 
    int digit = gval & MXT_HISTO_MASK;
    pair_it[histo[digit]++] = (SortItem) {gval, i-lwb};
  }
}

void scatter_digit(SortItem *in, SortItem *out, int digit_nr, ulong* histo, ulong size) {  
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i) {
    SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair;
  }
}

void scatter_last_digit(SortItem *in, ulong *out, int digit_nr, ulong* histo, ulong size){
  int shift = digit_nr * MXT_HISTO_SZ_LOG2;
  for (ulong i = 0; i != size; ++i)  {
    SortItem in_pair = in[i];
    int digit = (in_pair.val >> shift) & MXT_HISTO_MASK;
    out[histo[digit]++] = in_pair.rank;
  }
}

/*          Bit Array         */

BitArray *create_bit_array(ulong size) {
  BitArray *bit_array = malloc(sizeof(BitArray)); check_alloc(bit_array, 507);
  bit_array->num_words = (size + bits_per_word_log2() - 1) / bits_per_word_log2();
  bit_array->data = calloc(bit_array->num_words, sizeof(ulong)); check_alloc( bit_array->data, 508);
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

value_t transform(const float val) {    
  value_t valu = *((const value_t *)(&val));

  if (valu &    0x80000000)
    return      0xFFFFFFFF - valu;
        
  return valu | 0x80000000;
}


bool check_neighbor(bool *visited, PrioQueue *q, ulong *ranks, ulong* index, ulong* x, ulong* y, ulong *z, ulong* rank, ulong *dims, ulong n_x, ulong n_y, ulong n_z, ulong lwb){

  ulong n = n_z*(dims[0]*dims[1])+n_y*dims[0]+n_x;
  if(visited[n-lwb])
    return false;

 visited[n-lwb]=true;
  

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

void remaining(bool *visited, PrioQueue *q, ulong *ranks, ulong* index, ulong* rank, ulong* dims, ulong lwb, ulong upb, int connectivity, int periodic){
  ulong z = *index / (dims[0] * dims[1]);
  ulong x = (*index % (dims[0] * dims[1])) % dims[0];
  ulong y = (*index % (dims[0] * dims[1])) / dims[0];
  bool cond = 1;
  if(!periodic){
    while(cond){
      cond = 
	((x > 0)           && (*index -1 >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x - 1, y, z, lwb)) ||
	((x < dims[0] - 1) && (*index + 1 < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x + 1, y, z, lwb)) ||
	((y > 0)           && (*index - dims[0] >= lwb)       && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y - 1, z, lwb)) ||
	((y < dims[1] - 1)  && (*index + dims[0] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y + 1, z, lwb)) ||
	((z > 0)  && (*index >= lwb + dims[0]*dims[1])  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y, z - 1, lwb)) ||
	((z < dims[2] - 1)  && (*index+dims[0]*dims[1] < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y, z + 1, lwb));
    
      if(connectivity >= 8){
	cond = cond ||
	  ((x > 0) && (y > 0) && (*index-1-dims[0] >=lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x - 1, y - 1, z, lwb)) ||
	  ((x < dims[0]-1) && (y > 0) && (*index+1-dims[0] >=lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x + 1, y - 1, z, lwb)) ||
	  ((x > 0) && (y < dims[1]-1) && (*index-1+dims[0] < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x-1, y + 1, z, lwb)) ||
	  ((y < dims[1]-1) && (x < dims[0]-1) && (*index+dims[0]+1 < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x+1, y + 1, z, lwb)) ;   
      }
      if(connectivity == 26){
	cond = cond || (z > 0 &&
			(((y>0) && (x>0) && (*index -dims[0]*dims[1]-dims[0]-1 >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x-1, y-1, z-1, lwb)) ||
			 ((y>0) && (*index -dims[0]*dims[1]-dims[0] >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y-1, z-1, lwb)) ||
			 ((y>0) && (x<dims[0]-1) && (*index-dims[0]*dims[1]-dims[0]+1 >= lwb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x+1, y-1, z-1, lwb)) ||
			 ((x>0) && (*index -dims[0]*dims[1]-1 >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x-1, y, z-1, lwb)) ||
			 ((x<dims[0]-1) && (*index -dims[0]*dims[1]+1 >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x+1, y, z-1, lwb)) ||
			 ((y<dims[1]-1) && (x>0) && (*index -dims[0]*dims[1]+dims[0]-1 >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x-1, y+1, z-1, lwb)) ||
			 ((y<dims[1]-1)  && (*index -dims[0]*dims[1]+dims[0] >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y+1, z-1, lwb)) ||
			 ((y<dims[1]-1) && (x<dims[0]-1) && (*index -dims[0]*dims[1]+dims[0]+1 >= lwb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x+1, y+1, z-1, lwb))))
	  || (z<dims[2]-1 &&
	      (((y>0) && (x>0) && (*index +dims[0]*dims[1]-dims[0]-1 < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x-1, y-1, z+1, lwb)) ||
	       ((y>0) && (*index+dims[0]*dims[1]-dims[0] < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y-1, z+1, lwb)) ||
	       ((y>0) && (x<dims[0]-1) && (*index+dims[0]*dims[1]-dims[0]+1  < upb) && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x+1, y-1, z+1, lwb)) ||
	       ((x>0) && (*index+dims[0]*dims[1]-1  < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x-1, y, z+1, lwb)) ||
	       ((x<dims[0]-1) && (*index +dims[0]*dims[1]+1  < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x+1, y, z+1, lwb)) ||
	       ((y<dims[1]-1) && (x>0) && (*index +dims[0]*dims[1]+dims[0]-1  < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x-1, y+1, z+1, lwb)) ||
	       ((y<dims[1]-1)  && (*index +dims[0]*dims[1]+dims[0]  < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y+1, z+1, lwb)) ||
	       ((y<dims[1]-1) && (x<dims[0]-1) && (*index +dims[0]*dims[1]+dims[0]+1  < upb)  && check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x+1, y+1, z+1, lwb))));
      }
    }
  } else{
    ulong xbef = (long) (x - 1) >= 0? x-1:  dims[0]-1;
    ulong xaf  = (long) (x + 1) <  dims[0]? x+1: 0;
    ulong ybef = (long) (y - 1) >= 0? y-1:  dims[1]-1;
    ulong yaf  = (long) (y + 1) <  dims[1]?  y+1: 0;
    ulong zbef = (long) (z - 1) >= 0? z-1:  dims[2]-1;
    ulong zaf  = (long) (z + 1) <  dims[2]?  z+1: 0;
     while(cond){
      cond = 
	 check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, y, z, lwb) ||
	 check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, y, z, lwb) ||
	 check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, ybef, z, lwb) ||
	 check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, yaf, z, lwb) ||
	 check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y, zbef, lwb) ||
	 check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, y, zaf, lwb);
    
      if(connectivity >= 8){
	cond = cond ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, ybef, z, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, ybef, z, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, yaf, z, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, yaf, z, lwb) ;   
      }
      if(connectivity == 26){
	cond = cond ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, ybef, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, ybef, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, ybef, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, y, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, y, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, yaf, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, yaf, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, yaf, zbef, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, ybef, zaf, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, ybef, zaf, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, ybef, zaf, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, y, zaf, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, y, zaf, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xbef, yaf, zaf, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, x, yaf, zaf, lwb) ||
	  check_neighbor(visited, q, ranks, index, &x, &y, &z, rank, dims, xaf, yaf, zaf, lwb);
      }
    }
  }
  
}

bool is_border(bool border[6], ulong *dims, ulong p) {
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

void free_tree(Node *tree, ulong size_old){
  /* #pragma omp parallel for
  for (ulong i = 0; i < tree->size; i++ ) {
    if (tree->attribute[i]) {
      delete_aux_data(tree->attribute[i]);
    }
    }*/
  free(tree->parent);
  free(tree->gval);
  free(tree->attribute);
  //  if(tree->parent_qu != NULL) free(tree->parent_qu);
  // if(tree->gval_qu != NULL) free(tree->gval_qu);
  free(tree);
}



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


ulong *sort_image_pixels(Arguments *args, Node *tree, ulong *dims){
 
  int nthreads = args->threads_arg;
  int numstep = (int) abs(ceil((double)args->bpp_arg /(16)));
  ulong *sortedrs[2];
  ulong *histogramrs[2];
  ulong *histogram;
  histogramrs[0] = (ulong *) malloc(NUMBUCKETS*sizeof(ulong));
  histogramrs[1] = (ulong *) malloc(NUMBUCKETS*sizeof(ulong));
  sortedrs[0] = (ulong *) malloc(tree->size * sizeof(ulong));
  sortedrs[1] = (ulong *) malloc(tree->size * sizeof(ulong));

  ulong *loc_hist = calloc(nthreads*NUMBUCKETS, sizeof(ulong)); check_alloc(loc_hist, 610);

  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			/* Thread number */
    //  ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads)
    //     : dims[0]*dims[1]*((id*dims[2])/nthreads);              /* Lower bound for current thread */
    //   ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/ nthreads)
    //     : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);;       /* Upper bound for current thread */

    ulong 	lwb     = (id*tree->size)/nthreads;              /* Lower bound for current thread */
    ulong 	upb     = ((id+1)*tree->size)/nthreads;     /* Upper bound for current thread */
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

