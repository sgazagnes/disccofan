#include "types.h"
#include "moschini.h"
#include "tree_flood.h"
#include "queue.h"
#include "attributes.h"


void local_hist(value *gvals, ulong lwb, ulong upb, ulong *hist){
  ushort radix;
  for (ulong i=lwb; i<upb; ++i){
    radix = *( (ushort*) &(gvals[i]));
    hist[radix]++; // sort in ascending order
  }
  return;
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


/* QUANTIZED */


int boundaries(value *gvals,  ulong *px_start, ulong *px_end, ulong *sorted, ulong size, int nlevel){
  ulong num_pix_level = size/nlevel;
  int qtz_lv_curr = 0;
  value gval_prev = 0;
  ulong px = 0;
  int th_idx = 0;

  while(px<size){
    px_start[th_idx] = px;
    px += (num_pix_level - 1); // point to the last gvalue of that thread.
		
    if(px<size)	{
      gval_prev = gvals[sorted[px]];
      px++;
      
      while((px<size) && ((gval_prev == gvals[sorted[px]]) || (qtz_lv_curr == nlevel-1) ) ){
	px++;
      }
		
      px_end[th_idx] = px - 1;
    } else 
      px_end[th_idx] = size-1;
                       
    qtz_lv_curr++;        
    th_idx++;
  }
     
  return qtz_lv_curr; // return the number of quantization levels used
}

int create_quantized_image(value *gvals, value_t *gvals_qu, ulong *px_start, ulong *px_end, ulong *sorted, ulong size, int nlevel){

 ulong num_rem_pixels = size;
 ulong num_pix_level = size/nlevel;

 int qtz_lv_curr = 0;
 value gval_prev = 0;
 ulong px = 0, count = 0;
 int th_idx = 0;

 while(px<size){

   px_start[th_idx] = px;
   while(px<size && (count < num_pix_level || qtz_lv_curr == nlevel-1)) {
     gvals_qu[sorted[px]] = qtz_lv_curr;
     count++;
     gval_prev = gvals[sorted[px]];
     px++;
   }

   while(px<size && gval_prev == gvals[sorted[px]]) {
     gvals_qu[sorted[px]] = qtz_lv_curr ;
     count++;
     px++;
   }
   px_end[th_idx] = px-1;
                
   num_rem_pixels -= count;
   num_pix_level =  num_rem_pixels/(nlevel-qtz_lv_curr);
        
   qtz_lv_curr++;        
   th_idx++;
   count=0;
 }
 
 return qtz_lv_curr ; // return the number of quantization levels used
}


idx get_levelroot_qu(Node *tree, idx x) {
  /* index based, check whether this node is bottom */
  idx r = x;
  if (r == BOTTOM)
    return BOTTOM;
  value_t gv = tree->gval_qu[x]; 
  while ((tree->parent_qu[r] != BOTTOM) && (gv == tree->gval_qu[tree->parent_qu[r]])) { 
    r = tree->parent_qu[r];
  }
  while (x != r) {
    idx y = tree->parent_qu[x];
    tree->parent_qu[x] = r;
    x = y;
  }
  return r;
} /* get_levelroot */


idx get_parent_qu(Node *tree, idx x) {
  return get_levelroot_qu(tree, tree->parent_qu[x]);
} /* get_parent */


int tree_flood_qu(Node *tree, AuxDataStore *store, Queue *q,  idx *levelroot, bool *reached, ulong *dims, ulong lwb, ulong upb, int level, int connectivity, int qtz_lvl, void **thisattr) {

  void  *curr_attr  = NULL;
  void  *attr_child = NULL;
  int 	n_neighbors;
  int  fc;
  ulong neighbors[connectivity]; 
  ulong c, p, x, y, z;
  
  while (queue_is_not_empty(q, level)) {

    p = queue_first(q, level);
    x = p % (dims[0] * dims[1]) % dims[0];
    y = p % (dims[0] * dims[1]) / dims[0];
    z = p / (dims[0] * dims[1]);

    n_neighbors = get_neighbors(dims, lwb, upb, neighbors, p, x, y, z, connectivity);
    
    if (curr_attr){
      if(!is_border(tree->border, dims, p))
	add_to_aux_data(curr_attr, x, y, z);     
    } else{
      if(!is_border(tree->border, dims, p))
	curr_attr = new_aux_data(store, x, y, z);    
      if (*thisattr){
	if(curr_attr){
	  merge_aux_data(curr_attr, *thisattr);
	} else{
	  clone_aux_data(store, &curr_attr, *thisattr);
	}
      }
    }
    for (int i = 0; i < n_neighbors; i++) {
      c = neighbors[i];
      if (!reached[c-lwb]) {
        reached[c-lwb] = true;

        fc =  tree->gval_qu[c];
        if (levelroot[fc] == BOTTOM) {

          levelroot[fc] = c; 
        } else {
	  if( ((tree->gval[c] < tree->gval[levelroot[fc]]) || ((tree->gval[c] == tree->gval[levelroot[fc]]) && (c < (ulong) levelroot[fc]) )) ){
	    tree->parent_qu[levelroot[fc]] = c; //new lero
	    levelroot[fc] = c; //new lero						
	  }					
	  tree->parent_qu[c] = levelroot[fc];	
        }
        queue_add(q, fc, c);

        if (fc > level) { 
	  attr_child = NULL;
          do {
            fc = tree_flood_qu(tree, store, q, levelroot, reached, dims, lwb, upb, fc, connectivity, qtz_lvl, &attr_child);

            if ((int) fc >= qtz_lvl ) { 
	      if(curr_attr)
		delete_aux_data(curr_attr);
              return fc;
            }
          } while (fc != level);
	  if(attr_child != NULL){
	    if(curr_attr) merge_aux_data(curr_attr, attr_child);
	    else {
	      clone_aux_data(store, &curr_attr, attr_child);
	    }
	  }	
	}
      }
    }
  }

  int m =  level - 1;
  while (m > 0 && levelroot[m] == BOTTOM) m--;
  if (m >= 0) tree->parent_qu[levelroot[level]] = levelroot[m];
  else tree->parent_qu[levelroot[level]] = BOTTOM;
   if(curr_attr != NULL) tree->attribute[levelroot[level]] = curr_attr;
  
  levelroot[level] = BOTTOM;
  *thisattr = curr_attr;
  return m;
} /* tree_flood_qu */


void merge_nodes_qu(Node *tree, AuxDataStore *store, idx x, idx y) {
  void *cor  = NULL;
  void *copa = NULL;
  idx h, z;


  x =  get_levelroot_qu(tree, x);
  y =  get_levelroot_qu(tree, y);

  if ((tree->gval[x] < tree->gval[y]) || ((tree->gval[x] == tree->gval[y]) && (x < y))) {
    h=x; x=y; y=h;
  }
  while ((x != y) && (y != BOTTOM)) {
    z = get_parent_qu(tree, x);

    if ((z != BOTTOM) && ( (tree->gval[z]>tree->gval[y])|| ((tree->gval[z]==tree->gval[y]) && (z>=y))) ) {
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
      tree->parent_qu[x] = y;
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
      x = get_parent_qu(tree, x);
    }
  }
  if (cor)  delete_aux_data(cor);
  if (copa) delete_aux_data(copa);
}/* merge_nodes_qu */


void fuse_sections_qu(Node *tree, AuxDataStore *store, ulong *dims, uint id, uint i, uint n_t, int connectivity){
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
     min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v]);
     if (!test_min || min_curr > min_prev)
       merge_nodes_qu(tree, store, u, v);
     min_prev = min_curr;
     test_min = true;
     if(connectivity >= 8){
       if(u_x > 0){
	 min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v-1]);
	 if (min_curr > min_prev)
	   merge_nodes_qu(tree, store, u, v-1);
       }     
       if(u_x < dims[0]-1){
	 min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v+1]);
	 if (!test_min || min_curr > min_prev)
	   merge_nodes_qu(tree, store, u, v+1);
       }
     }
    }
  } else {
    for (p=0, u = mdb; p < (dims[0] * dims[1]); p++, u++){
      u_x = (u % (dims[0]*dims[1])) % dims[0];
      u_y = (u % (dims[0]*dims[1])) / dims[0];
      v   = u - dims[0]*dims[1];

      if(u_x == 0) test_min = false;
      min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v]);
      if (!test_min || min_curr > min_prev )
	merge_nodes_qu(tree, store, u, v);
      min_prev = min_curr;
      test_min = true;
		
      if(connectivity == 26){
	if(u_x > 0 && u_y > 0){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v-dims[0]-1]);
	  if (min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v-dims[0]-1);
	}
	if(u_y > 0){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v-dims[0]]);
	  if (min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v-dims[0]);
	}
	if(u_y > 0 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v-dims[0]+1]);
	  if (min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v-dims[0]+1);
	}
	if(u_x>0){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v-1]);
	  if (min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v-1);	      
	}
	if(u_x<dims[0]-1){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v+1);
	}
	if(u_y < dims[1]-1 && u_x>0){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v+dims[0]-1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v+dims[0]-1);
	}
	if(u_y < dims[1]-1){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v+dims[0]]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v+dims[0]);
	}
	if(u_y < dims[1]-1 && u_x<dims[0]-1){
	  min_curr = MIN(tree->gval_qu[u], tree->gval_qu[v+dims[0]+1]);
	  if (!test_min || min_curr > min_prev )
	    merge_nodes_qu(tree, store, u, v+dims[0]+1);
	}
      }
    }
  }
}



idx descend_root(value_t *gval, idx *parent, long q, value_t myLev){
  idx curr = q;
  while(gval[ parent[curr]] > myLev) {
    curr = parent[curr];		
  }
  return curr;
}


idx find_root(idx *zpar, idx p) {
  idx first = p;
  while(p != zpar[p])
    {
      p = zpar[p];		
    }
	
   idx nextEl;
  idx next = first;
  while(next != zpar[p])
    {
      nextEl = zpar[next];
      zpar[next] = zpar[p];
      next = nextEl;	
      }
  return zpar[p];		
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



value_t *calculate_quantized_image(Node *tree, ulong *sorted, ulong *px_start, ulong *px_end, int nthreads){
  
  value_t *gval = malloc(tree->size*sizeof(value_t)); check_alloc(gval, 500);
  int qtz_lvl_num = create_quantized_image(tree->gval, gval, px_start, px_end, sorted, tree->size, nthreads);
  if(nthreads != qtz_lvl_num)
    error("Wrong number of quantize levels");
  return gval;
}



idx *build_quant_tree(Arguments *args, Node *tree,  ulong *dims){
  int connectivity =  args->connectivity_arg;
  int nthreads     =  args->threads_arg;
  // int attrib 	   =  args->attribute_arg;
 
  omp_lock_t lock[nthreads];
  for (int i=0; i<nthreads; i++)
    omp_init_lock(&(lock[i]));
  int *saval = calloc(nthreads, sizeof(int));
  
  tree->parent_qu = malloc(tree->size * sizeof(idx)); check_alloc(tree->parent_qu, 400);
  // AuxDataStore **store_th =  malloc(nthreads * sizeof(AuxDataStore*)); check_alloc(store_th, 0);
  //init_aux_data_store(store_th, AttribsArray[attrib].size, upb-lwb);
  
  #pragma omp parallel num_threads(nthreads)
  {
    int 	id      = omp_get_thread_num();			/* Thread number */
    ulong 	lwb     = dims[2] == 1 ? dims[0]*((id*dims[1])/nthreads)
      : dims[0]*dims[1]*((id*dims[2])/nthreads);              /* Lower bound for current thread */
    ulong 	upb     = dims[2] == 1 ? dims[0]*(((id+1)*dims[1])/ nthreads)
      : dims[0]*dims[1]*(((id+1)*dims[2]) / nthreads);;       /* Upper bound for current thread */

    ulong 	min_idx     = lwb;
    void 	*curr_attr  = NULL;	
    ulong 	*histogram  = calloc(nthreads, sizeof(ulong)); check_alloc(histogram,  601);
    idx		*levelroots = calloc(nthreads, sizeof(idx));   check_alloc(levelroots, 602);
    bool 	*reached    = calloc(upb - lwb,    sizeof(bool));  check_alloc(reached, 603);
    Queue 	*queue      = create_queue(upb - lwb, nthreads);

    for (int i = 0; i < nthreads; i++) 
      levelroots[i] = BOTTOM;

    for (ulong i = lwb; i < upb; i++) {
      histogram[tree->gval_qu[i]]++;
      tree->parent_qu[i] = -1;
      if (tree->gval_qu[min_idx] > tree->gval_qu[i])
	min_idx = i;
    }

    //  store_th[id] =  malloc(sizeof(AuxDataStore)); check_alloc(store_th[id], 0);
    //  init_aux_data_store(store_th[id], AttribsArray[attrib].size, upb-lwb);
  
    set_queue_offsets(queue, histogram, nthreads);
    free(histogram);

    queue_add(queue, tree->gval_qu[min_idx], min_idx);
    levelroots[tree->gval_qu[min_idx]] = min_idx;
    reached[min_idx-lwb]   = true;

    tree_flood_qu(tree, NULL, queue, levelroots, reached,  dims, lwb, upb, (value_t) tree->gval_qu[min_idx], connectivity, nthreads, &curr_attr);

    free(queue);
    free(levelroots);
    free(reached);

    int i = 1;
    int qq = id;
    while (id+i < nthreads && qq%2 == 0){
		
      while (saval[id+i] <= 0){
	#pragma omp flush // wait
      }
      omp_set_lock(&lock[id+i]);
      saval[id+i]--;
      omp_unset_lock(&lock[id+i]);
      fuse_sections_qu(tree, NULL, dims, id, i, nthreads, connectivity);
      i *= 2;
      qq /= 2;      
    }
    if(id>0){
      omp_set_lock(&lock[id]);	      
      saval[id]++;
      omp_unset_lock(&lock[id]);
    }
  }

  free(saval);
  for (int i = 0; i < nthreads; i++)
    omp_destroy_lock(&(lock[i]));

  return tree->parent_qu;
}

idx *refine_tree(Arguments *args, Node *tree, ulong *dims, ulong *px_start, ulong *px_end, ulong *sorted){

  int connectivity =  args->connectivity_arg;
  int nthreads     =  args->threads_arg;
  // int attrib 	   =  args->attribute_arg;
  tree->parent = malloc(tree->size*sizeof(idx));
  idx *zpar = malloc(tree->size*sizeof(idx));
  void **attribute = calloc(tree->size ,sizeof(void*));
  memset(zpar, -1, tree->size*sizeof(idx));
  memset(tree->parent, -1, tree->size*sizeof(idx));

  #pragma omp parallel num_threads(nthreads) shared(zpar)
  {
    uint 	id      = omp_get_thread_num();			/* Thread number */
    idx 	lwb = (idx) px_start[id];			/* Lower bound for current thread */
    idx 	upb = (idx)  px_end[id];	
    int numneighbors;
    idx neighbors[connectivity];
    idx  i, p, q, r, ancest, tobezipped;
    ulong x, y, z;


    // AuxDataStore *store_th =  malloc(sizeof(AuxDataStore)); check_alloc(store_th, 0);
    //  init_aux_data_store(store_th, AttribsArray[attrib].size, upb-lwb);
    //
    void *store_th = NULL;
    // go through the sorted pixels of your partition from the highest to the lowest intensity
    for( i = upb; i >= lwb; i--) {
      p = sorted[i];

      zpar[p] = p;
      x = p % (dims[0] * dims[1]) % dims[0];
      y = p % (dims[0] * dims[1]) / dims[0];
      z = p / (dims[0] * dims[1]);

      if(!is_border(tree->border, dims, p))
	attribute[p] = new_aux_data(store_th, x, y, z);
      
      numneighbors = get_neighbors(dims, 0, tree->size, neighbors, p, x, y, z, connectivity);
		
      for (int j=0; j<numneighbors; j++){	
	q = neighbors[j];

	if(tree->gval_qu[q] > id) {
	  ancest = descend_root(tree->gval_qu, tree->parent_qu, q, id);
	  if (tree->parent[ancest] == BOTTOM) {
	    tree->parent[ancest] = p;

	    if(tree->attribute[ancest]){
	      if(attribute[p]){
		merge_aux_data(attribute[p], tree->attribute[ancest]);
	      } else {
		clone_aux_data(store_th, &attribute[p], tree->attribute[ancest]);
	      }
	    }
	  } else {
	    tobezipped = find_root(zpar, tree->parent[ancest]);	       
	       	
	    if(tobezipped != p)  {
	      tree->parent[tobezipped] = p;
	      zpar[tobezipped] = p;
	      if(attribute[tobezipped]){
		if( attribute[p]){
		  merge_aux_data(attribute[p], attribute[tobezipped]);
		} else {
		  clone_aux_data(store_th, &attribute[p], attribute[tobezipped]);
		}
	      }
	    }			
	  }											
	} else if(tree->gval_qu[q] == id)  {	 		
	  if(zpar[q] != BOTTOM)  {
	    r = find_root(zpar, q);
	    if(r != p) {
	      tree->parent[r] = p;
	      zpar[r] = p;
	      if(attribute[r]){
		if(attribute[p]){
		  merge_aux_data(attribute[p], attribute[r]);
		} else {
		  clone_aux_data(store_th, &attribute[p], attribute[r]);
		}
	      }
	    }				
	  }

	}			
      }		
    }
  }

  #pragma omp parallel for
  for(ulong i = 0; i < tree->size; i++)
    if(tree->attribute[i]) free(tree->attribute[i]);
  
  tree->attribute = attribute;
  return tree->parent;
}

