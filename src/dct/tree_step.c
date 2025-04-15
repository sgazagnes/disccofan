#include "types.h"
#include "attributes.h"
#include "tree_step.h"
#include "tree_flood.h"
#include "queue.h"


void build_tree_parents(Arguments *args, Node *tree, ulong *dims){
  int flood_algo   =  args->flood_arg;
  int connectivity =  args->connectivity_arg;
  int nthreads     =  args->threads_arg;


  omp_lock_t lock[nthreads];
  for (int i=0; i<nthreads; i++)
    omp_init_lock(&(lock[i]));
  int *saval = calloc(nthreads, sizeof(int));

  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			/* Thread number */
    ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads)
      : dims[0]*dims[1]*((id*dims[2])/nthreads);              /* Lower bound for current thread */
    ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/nthreads)
      : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);;       /* Upper bound for current thread */

    
    if(flood_algo == 0){
      /* Salembier algorithm */
      
      ulong 	min_idx     = lwb;
      ulong 	*histogram  = calloc(g_max_levels, sizeof(ulong)); check_alloc(histogram,  601);
      idx	*levelroots = calloc(g_max_levels, sizeof(idx));   check_alloc(levelroots, 602);
      bool 	*reached    = calloc(upb - lwb,    sizeof(bool));  check_alloc(reached, 603);
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
      tree_par_sal(tree, queue, levelroots, reached, dims, lwb, upb, (value_t) tree->gval[min_idx], connectivity);

      free(queue);
      free(levelroots);
      free(reached);
      
    } else if (flood_algo == 1) {
	/* Wilkinson improved algorithm */ 
      PrioQueue *q        = create_prio_queue(upb-lwb);	    
      ulong *ranks        = calloc(upb-lwb, sizeof(ulong));  check_alloc(ranks, 604);
      ulong *ranks_inv    = calloc(upb-lwb, sizeof(ulong));  check_alloc(ranks_inv, 604);
      bool *visited   	  = calloc(upb-lwb, sizeof(bool));  check_alloc(visited, 603);

      create_mappings(tree->gval, ranks, ranks_inv, upb-lwb, lwb, upb);
      tree_flood_tee_par(tree, visited, q, ranks, ranks_inv, dims, lwb, upb, connectivity);
      
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

	fuse_parents(tree, dims, id, i, nthreads, connectivity);
	
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

  
  // return &(tree->parent);
} /* build_tree_step */



long tree_par_sal(Node *tree,  Queue *q,  idx *levelroot, bool *reached, ulong *dims, ulong lwb, ulong upb, long level, int connectivity) {

  int 	n_neighbors;
  long  fc;
  ulong neighbors[connectivity]; 
  ulong c, p, x, y, z;
  
  while (queue_is_not_empty(q, level)) {

    p = queue_first(q, level); 
    x = p % (dims[0] * dims[1]) % dims[0];
    y = p % (dims[0] * dims[1]) / dims[0];
    z = p / (dims[0] * dims[1]);

    n_neighbors = get_neighbors(dims, lwb, upb, neighbors, p, x, y, z, connectivity);
             
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
          do {
            fc = tree_par_sal(tree, q, levelroot, reached, dims, lwb, upb, fc, connectivity);
            if ((ulong) fc >= g_max_levels) { 
              return fc;
            }
          } while (fc != level);
	}
      }
    }
  }

  long m = (long) level - 1;
  while (m > 0 && levelroot[m] == BOTTOM) m--;
  if (m >= 0) tree->parent[levelroot[level]] = levelroot[m];
  else tree->parent[levelroot[level]] = BOTTOM;
  levelroot[level] = BOTTOM;
  return m;
} /* tree_flood_sal_par */


void tree_flood_tee_par(Node *tree, bool *visited, PrioQueue *q, ulong *ranks, ulong* ranks_inv, ulong *dims,  ulong lwb, ulong upb, int connectivity){
  
  ulong index = lwb; 
  ulong rank = ranks[index-lwb];
  visited[index-lwb] = true;
  while(true){
    remaining(visited, q, ranks, &index, &rank, dims, lwb, upb, connectivity,0);

    if (q->m_levels[0][0] == 0) break;
    rank = q->m_top;
    prio_queue_remove(q);
    ulong parent = ranks_inv[rank]+lwb;

    tree->parent[index] = parent;
    index = parent;
  }
  tree->parent[index] = BOTTOM;
} /* tree_flood_tee_par */


/*void tree_flood_tee_att(Node *tree, AuxDataStore *store,  bool *visited,  PrioQueue *q, ulong *ranks, ulong* ranks_inv,  ulong *dims,  ulong *attr_off, ulong lwb, ulong upb, int connectivity){

  ulong index = lwb;
  idx parent_l;
  ulong rank = ranks[index-lwb];
  visited[index-lwb] = true;
  while(true){
    remaining(visited, q, ranks, &index, &rank, dims, lwb, upb, connectivity);
    info("Index %d, gval %d", index,tree->gval[index]);
    parent_l = get_parent(tree,index);
    info("Parent_L %d, gval %d", parent_l, tree->gval[parent_l]);

    if(parent_l == BOTTOM) break;
    if(!tree->attribute[index] && !is_border(tree->border, dims, index) && index < upb){
      tree->attribute[index] = new_aux_data(store, (index % (dims[0] * dims[1])) % dims[0]+ attr_off[0] , (index % (dims[0] * dims[1])) / dims[0]+ attr_off[1], index / (dims[0] * dims[1])+ attr_off[2]);
    }
    while(parent_l >= upb && parent_l != BOTTOM){
      if(tree->attribute[index] && tree->attribute[parent_l])
	merge_aux_data(tree->attribute[parent_l], tree->attribute[index]);
      else if( tree->attribute[index])
	clone_aux_data(store, &tree->attribute[parent_l], tree->attribute[index]);
      parent_l = get_parent(tree,parent_l);
    }
    
    info("Parent_L %d, gval %d", parent_l, tree->gval[parent_l]);

    if(parent_l == BOTTOM) break;

    if(!tree->attribute[parent_l] && !is_border(tree->border, dims, parent_l) && parent_l < upb){
      tree->attribute[parent_l] = new_aux_data(store, (parent_l % (dims[0] * dims[1])) % dims[0]+ attr_off[0], (parent_l % (dims[0] * dims[1])) / dims[0]+ attr_off[1], parent_l / (dims[0] * dims[1])+ attr_off[2]);
    }

    if(tree->attribute[index]){
      if((!is_border(tree->border, dims, parent_l) && parent_l < upb) || tree->attribute[parent_l])
	merge_aux_data(tree->attribute[parent_l], tree->attribute[index]);
      else{
	clone_aux_data(store, &tree->attribute[parent_l], tree->attribute[index]);
      }
    }
  
   
    if (q->m_levels[0][0] == 0) break;
    rank = q->m_top;
    if(ranks[parent_l] > rank)
      index = parent_l;
    else{
      index = ranks_inv[rank]+lwb;
      prio_queue_remove(q);
    }
  }
  while(parent_l != -1){
    clone_aux_data(store, &tree->attribute[parent_l], tree->attribute[index]);
    index = parent_l;
    parent_l = tree->parent[parent_l];
  }
} /* tree_flood_tee_att */


void fuse_parents(Node *tree, ulong *dims, uint id, uint i, uint n_t, int connectivity){
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
       merge_nodes_par(tree, u, v);
     min_prev = min_curr;
     test_min = true;
     if(connectivity >= 8){
       if(u_x > 0){
	 min_curr = MIN(tree->gval[u], tree->gval[v-1]);
	 if (min_curr > min_prev)
	   merge_nodes_par(tree, u, v-1);
       }     
       if(u_x < dims[0]-1){
	 min_curr = MIN(tree->gval[u], tree->gval[v+1]);
	 if (!test_min || min_curr > min_prev)
	   merge_nodes_par(tree, u, v+1);
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
	merge_nodes_par(tree, u, v);
      min_prev = min_curr;
      test_min = true;
		
      if(connectivity == 26){
	if(u_x > 0 && u_y > 0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]-1]);
	  if (min_curr > min_prev )
	    merge_nodes_par(tree, u, v-dims[0]-1);
	}
	if(u_y > 0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]]);
	  if (min_curr > min_prev )
	    merge_nodes_par(tree, u, v-dims[0]);
	}
	if(u_y > 0 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]+1]);
	  if (min_curr > min_prev )
	    merge_nodes_par(tree, u, v-dims[0]+1);
	}
	if(u_x>0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-1]);
	  if (min_curr > min_prev )
	    merge_nodes_par(tree, u, v-1);	      
	}
	if(u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_par(tree, u, v+1);
	}
	if(u_y < dims[1]-1 && u_x>0){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]-1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_par(tree, u, v+dims[0]-1);
	}
	if(u_y < dims[1]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_par(tree, u, v+dims[0]);
	}
	if(u_y < dims[1]-1 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_par(tree, u, v+dims[0]+1);
	}
      }
    }
  }
}

void fuse_attributes(Node *tree, AuxDataStore *store, ulong *dims, uint id, uint i, uint n_t, int connectivity){
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
       merge_nodes_att(tree, store, u, v);
     min_prev = min_curr;
     test_min = true;
     if(connectivity >= 8){
       if(u_x > 0){
	 min_curr = MIN(tree->gval[u], tree->gval[v-1]);


	 if (min_curr > min_prev)
	   merge_nodes_att(tree, store, u, v-1);
       }     
       if(u_x < dims[0]-1){
	 min_curr = MIN(tree->gval[u], tree->gval[v+1]);
	 if (!test_min || min_curr > min_prev)
	   merge_nodes_att(tree, store, u, v+1);
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
	merge_nodes_att(tree, store, u, v);
      min_prev = min_curr;
      test_min = true;
		
      if(connectivity == 26){
	if(u_x > 0 && u_y > 0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]-1]);
	  if (min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v-dims[0]-1);
	}
	if(u_y > 0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]]);
	  if (min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v-dims[0]);
	}
	if(u_y > 0 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v-dims[0]+1]);
	  if (min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v-dims[0]+1);
	}
	if(u_x>0){
	  min_curr = MIN(tree->gval[u], tree->gval[v-1]);
	  if (min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v-1);	      
	}
	if(u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v+1);
	}
	if(u_y < dims[1]-1 && u_x>0){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]-1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v+dims[0]-1);
	}
	if(u_y < dims[1]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v+dims[0]);
	}
	if(u_y < dims[1]-1 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval[u], tree->gval[v+dims[0]+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_att(tree, store, u, v+dims[0]+1);
	}
      }
    }
  }
}


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


void merge_nodes_att(Node *tree, AuxDataStore *store, idx x, idx y) {
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
      	clone_aux_data(store, &(tree->attribute[x]), copa);
      }
      x = y;
      y = z;
    }
  }
  if (y == BOTTOM && cor) {
    while(x != BOTTOM) {
      if (tree->attribute[x])
	merge_aux_data(tree->attribute[x], cor);
      else{
	clone_aux_data(store, &tree->attribute[x], cor);
      }
      x = get_parent(tree, x);
    }
  }
  if (cor)  delete_aux_data(cor);
  if (copa) delete_aux_data(copa);
}/* merge_nodes */

