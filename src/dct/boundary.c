#include "types.h"
#include "boundary.h"
#include "communication.h"
#include "attributes.h"
#include "tree_flood.h"

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Main Functions        */
/*				   */
/* +++++++++++++++++++++++++++++++ */

Node *correct_borders(Arguments *args, Node *local_tree, ulong *dims) {
  /* Initialize variables */
  
  Boundary 	**bound_tree;            
  int 		merged      = 0;
  int 		connectivity= args->connectivity_arg;	
  int 		attribute   = args->attribute_arg;
  int 		grid[3]     = {args->grid_arg[0], args->grid_arg[1], args->grid_arg[2]};
  int 		base[3]     = {1, 1, 1};
  int 		grid_cur[3] = {grid[0], grid[1], grid[2]};
  int 		myrank      = rank();
  int 		myslice     = myrank / (grid[0] * grid[1]);
  int 		myrank_2D   = myrank % (grid[0] * grid[1]);
  int 		mycol_2D    = myrank_2D % grid[0];
  int 		myrow_2D    = myrank_2D / grid[0];
  int		n_merge     = 0;
  int 		neighrank;
  ulong 	store_item;


  /* Find the size of the boundary tree array for this process */
  
  if (mycol_2D % 2)
    n_merge = 0;
  else if (myrow_2D % 2){
    if(mycol_2D < grid[0]-1) n_merge++;
  } else {
    base[0] = grid[0] - mycol_2D +1 ;
    base[1] = grid[1] - myrow_2D + 1;
    base[2] = grid[2] - myslice + 1;
    if(mycol_2D < grid[0]-1){
      while (base[0] / 2) {
	n_merge++;  base[0] /= 2;
      }
    }
    if(myrow_2D < grid[1]-1){
      while (base[1] / 2) {
	n_merge++;  base[1] /= 2;
      }
    }
    if(myslice < grid[2]-1){
      while (base[2] / 2) {
	n_merge++;  base[2] /= 2;
      }
    }
    base[0] =  base[1] =  base[2] = 1;
  }
  bound_tree = malloc((1 + 2*n_merge)*sizeof(Boundary*));  check_alloc(bound_tree, 200);
  n_merge = 0;

  /* Initialize the boundary tree of the local component tree */

  bound_tree[0] = create_boundary(local_tree, dims, attribute, connectivity);
  store_item 	= bound_tree[0]->store->size_item;

  timing("First boundary tree construction: wallclock time = %0.2f",
	 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	  
  /* Merging */
  while ((grid_cur[0] > 1) || (grid_cur[1] > 1) || (grid_cur[2] > 1)) {

    /* Horizontally */
    if ((grid_cur[0] > 1)) {
      if (!merged) {
	if (mycol_2D % (2*base[0]) == 0) {
	  if(mycol_2D  < grid[0]-1){
	    neighrank = myrank + base[0];
	    bound_tree[n_merge+1] = receive_boundary(neighrank, store_item);
	    merge(bound_tree[n_merge], bound_tree[n_merge+1], HORIZONTAL);
	    timing("Merging horizontally: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	    bound_tree[n_merge+2] = combine(bound_tree[n_merge], bound_tree[n_merge+1], HORIZONTAL);
	    timing("Combining horizontally: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	    n_merge += 2;
	  }
	} else {
	  neighrank = myrank - base[0];
	  send_boundary(bound_tree[n_merge], neighrank);
	  merged++;
	}
      } else {
	merged++;
      }
      base[0]     *= 2;
      grid_cur[0]++;
      grid_cur[0] /= 2;
    }

    /* Vertically */
    if ((grid_cur[1] > 1)) {
      if (!merged) {
	if (myrow_2D % (2*base[1]) == 0) {
	  if(myrow_2D < grid[1] - 1) { 
	    neighrank = myrank + base[1]*grid[0];
	    bound_tree[n_merge+1] = receive_boundary(neighrank, store_item);
	    merge(bound_tree[n_merge], bound_tree[n_merge+1], VERTICAL);
	    timing("Merging vertically: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));	 
	    bound_tree[n_merge+2] = combine(bound_tree[n_merge], bound_tree[n_merge+1], VERTICAL);
	    timing("Combining vertically: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	    n_merge += 2;
	  }
	} else {
	  neighrank = myrank - base[1]*grid[0];
	  send_boundary(bound_tree[n_merge], neighrank);
	  merged++;
	}
      } else {
	merged++;
      }
      base[1]     *= 2;
      grid_cur[1]++;
      grid_cur[1] /= 2;
    }

    /* In depth */
    if ((grid_cur[2] > 1)) { 
      if (!merged) {
	if (myslice % (2*base[2]) == 0) {
	  if(myslice < grid[2] - 1){
	    neighrank = myrank + base[2]*grid[0]*grid[1];
	    bound_tree[n_merge+1] = receive_boundary(neighrank, store_item);
	    merge(bound_tree[n_merge], bound_tree[n_merge+1], DEPTH);
	    timing("Merging in depth: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));	
	    bound_tree[n_merge+2] = combine(bound_tree[n_merge], bound_tree[n_merge+1], DEPTH);
	    timing("Combining in depth: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));
	    n_merge += 2;
	  }
	} else {
	  neighrank = myrank - base[2]*grid[0]*grid[1];
	  send_boundary(bound_tree[n_merge], neighrank);
	  merged++;
	}
      } else {
	merged++;
      }
      base[2]     *= 2;
      grid_cur[2] /= 2;
    }
  }

  base[0] /= 2;
  base[1] /= 2;
  base[2] /= 2;
  grid_cur[0] = grid[0];
  grid_cur[1] = grid[1];
  grid_cur[2] = grid[2];
  merged--;

  /*       Updating the boundary trees    */
  
  while ((grid_cur[0] != 1) || (grid_cur[1] != 1) || (grid_cur[2] != 1)) {

    /* In depth */
    if (grid_cur[2] >= grid_cur[1] && grid_cur[2] >= grid_cur[0] && grid_cur[2] > 1) {
      if (merged<=0) {
	if (myslice % (2*base[2]) == 0) {
	  if (myslice < grid[2] - 1) {
	    neighrank = myrank + base[2]*grid[0]*grid[1];
	    update_par_step( bound_tree[n_merge]);
	    update(bound_tree[n_merge-2], bound_tree[n_merge]);
	    update(bound_tree[n_merge-1], bound_tree[n_merge]);
	    send_updated_boundary(bound_tree[n_merge-1], neighrank);
	    timing("Updating in depth: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));	
	    free_boundary(bound_tree[n_merge]);
	    free_boundary(bound_tree[n_merge-1]);
	    n_merge -= 2;
	  }
	} else {
	  neighrank = myrank - base[2]*grid[0]*grid[1];
	  bound_tree[n_merge] = receive_updated_boundary(bound_tree[n_merge], neighrank);
	  merged--;
	}
      }else {
	merged--;
      }
      base[2] 	  /= 2;
      grid_cur[2] /= 2;
    }    
  
     /* Vertically */
    if (grid_cur[1] >= grid_cur[0] && grid_cur[1] > grid_cur[2] && grid_cur[1] > 1) {	    
      if (merged<=0) {
	if (myrow_2D % (2*base[1]) == 0) {
	  if(myrow_2D  < grid[1]-1){
	    neighrank = myrank + base[1]*grid[0];
	    update_par_step( bound_tree[n_merge]);

	    update(bound_tree[n_merge-2], bound_tree[n_merge]);
	    update(bound_tree[n_merge-1],  bound_tree[n_merge]);
	    send_updated_boundary(bound_tree[n_merge-1], neighrank);
	    timing("Updating vertically: wallclock time = %0.2f",
		   (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));	
	    free_boundary(bound_tree[n_merge]);
	    free_boundary(bound_tree[n_merge-1]);
	    n_merge -= 2;
	  }
	} else {
	  neighrank = myrank - base[1]*grid[0];      	  
	  bound_tree[n_merge] = receive_updated_boundary(bound_tree[n_merge], neighrank);
	  merged--;
	}
      }else {
	merged--;
      }
      base[1] 	  /= 2;
      grid_cur[1]++;
      grid_cur[1] /= 2;
    }

    /* Horizontally */
    if (grid_cur[0] > grid_cur[1] && grid_cur[0] > grid_cur[2] && grid_cur[0] > 1) {
      if (merged<=0) {
	if (mycol_2D % (2*base[0]) == 0){
	  if (mycol_2D < grid[0]-1) {
	  neighrank = myrank + base[0];
	  update_par_step( bound_tree[n_merge]);	    
	  update(bound_tree[n_merge-2], bound_tree[n_merge]);
	  update(bound_tree[n_merge-1], bound_tree[n_merge]);
	  send_updated_boundary(bound_tree[n_merge-1], neighrank);
	  timing("Updating horizontally: wallclock time = %0.2f",
		 (float)(times(&tstruct) - start)/(float)sysconf(_SC_CLK_TCK));	
	  free_boundary(bound_tree[n_merge]);
	  free_boundary(bound_tree[n_merge-1]);
	  n_merge -= 2;
	  }
	} else {
	  neighrank = myrank - base[0];
	  bound_tree[n_merge] = receive_updated_boundary(bound_tree[n_merge], neighrank);
	  merged--;
	}
      } else {
	merged--;
      } 
      base[0] 	  /= 2;
      grid_cur[0]++;
      grid_cur[0] /= 2;
    }
  }
  local_tree = correct_local_tree(local_tree, bound_tree[0]);
  debug("Boundary tree size: Init %ld, Final %ld", bound_tree[0]->size_init, bound_tree[0]->size_curr);
  free_boundary(bound_tree[0]);
  return local_tree;
} /* correct_borders */


Node *correct_local_tree(Node *local_tree, Boundary *b) {
  /* Adding and updating the nodes from the local tree */

  void 		*b_attr;
  idx		b_parent;
  ulong 	size_tree_init = local_tree->size;
  ulong 	size_tree_new  = (local_tree->size + b->size_curr - b->size_init);

  local_tree->parent = realloc(local_tree->parent, size_tree_new * sizeof(idx));
  local_tree->gval = realloc(local_tree->gval, size_tree_new * sizeof(value));
  local_tree->attribute = realloc(local_tree->attribute, size_tree_new * sizeof(void*)); 

  for (ulong i = 0; i < b->size_curr; i++) {
    if (i >= b->size_init) {
      b->array[i].index = local_tree->size;
      bound_to_tree(local_tree, b, i, local_tree->size);
      (local_tree->size)++;
    }
    
    if (b->attribute_idx[i] != BOTTOM) {
      b_attr = (char *) b->store->data + b->attribute_idx[i] * b->store->size_item;      
      clone_aux_data(NULL, &local_tree->attribute[b->array[i].index], b_attr);
    }
    
    b_parent = b->border_par[i].i;

    if (b_parent == BOTTOM)
      local_tree->parent[b->array[i].index] = BOTTOM;    
    else if ((ulong) b_parent < b->size_init)
      local_tree->parent[b->array[i].index] = b->array[b_parent].index;
    else
      local_tree->parent[b->array[i].index] = size_tree_init + b_parent - b->size_init;   
  }
  
  return local_tree;
} /* apply_changes */

void update_par_step(Boundary *c) {

  if(c == NULL) return;
  for (size_t i = 0; i< c->size_curr; i++) {
    BorderIndex x = (BorderIndex) {.b = c, .i = i};
    BorderIndex z = b_levelroot(x);

    if(bi_equal(x,z)) continue;

    BorderIndex origin_a, origin_c, origin_lr;
    origin_a = c->border_ori[x.i];
    if(is_bottom(origin_a)) continue;
    
    origin_c = c->border_ori[z.i];
    origin_lr = b_levelroot(origin_a);
    if(bi_equal(origin_lr, origin_c)) continue;
    
    BorderIndex lr_a = origin_a.b->border_lr[origin_a.i], lr_c=  origin_c.b->border_lr[origin_c.i];
    if(origin_a.b == origin_c.b){
      origin_a.b->border_par[origin_a.i] = origin_c;
      if(!is_bottom(lr_a) && !is_bottom(lr_c)){
	lr_a.b->border_par[lr_a.i] = lr_c;
      } else if(!is_bottom(lr_a)){
	origin_c.b->border_lr[origin_c.i] = lr_a;
      }
      else if(!is_bottom(lr_c)) {
	origin_a.b->border_lr[origin_a.i] = lr_c;
      }
    } else if (!is_bottom(lr_c) && origin_a.b == lr_c.b){
      origin_a.b->border_par[origin_a.i] = lr_c;
      if(!is_bottom(lr_a)){
	lr_a.b->border_par[lr_a.i] = origin_c;
      } 
    } else if (is_bottom(lr_c) && !is_bottom(lr_a)){
      lr_a.b->border_par[lr_a.i] = origin_c;
      origin_c.b->border_lr[origin_c.i] = origin_a;
      origin_a.b->border_lr[origin_a.i] = origin_c;
    } else if(is_bottom(lr_c) && is_bottom(lr_a)) {
      origin_c.b->border_lr[origin_c.i] = origin_a;
      origin_a.b->border_lr[origin_a.i] = origin_c;
    } 
  }

}

Boundary *update(Boundary *b, Boundary *c) {
  /* Update of the boundary tree b using the combined, already updated, boundary tree c if it exists */
  
  BorderIndex x, z, s;
  b->reached = calloc(b->size_alloc, sizeof(bool)); check_alloc(b->reached, 202);
  for (size_t i = 0; i< b->size_curr; i++) {
    x = (BorderIndex) {.b = b, .i = i};
    x = b_id_levelroot(x);
    if (x.b->reached[x.i]) continue;
    z = b_levelroot(x);    
    s = (BorderIndex) {.b = c, .i = z.b->array[z.i].border_idx };  
    b = update_branch(x,z,s);
  }    

  free(b->reached);
  return b;
}

Boundary *update_branch(BorderIndex x, BorderIndex z, BorderIndex s){
  /* Update the branch of the node x using its levelroot z in the merged tree, and the corresponding updated node s in the combined tree if it exists. */

  /* If s does not exist (node has not been added in the combined tree) */
  while (is_bottom(s)  && !(x.b->reached[x.i])) {
    if(!bi_equal(x,z))
      update_node_attribute(x,z);
    z = b_parent_lr(z);
    if(is_bottom(z)){
      x.b->border_par[x.i] = z;
      break;
    }
    if (x.b == z.b)
      x.b->border_par[x.i] = z;
    else {
      if (!is_bottom(z.b->border_lr[z.i]))
	x.b->border_par[x.i] = z.b->border_lr[z.i];
      else {
	x.b = adding_node(x,z);
	x.b->array[x.b->size_curr-1].border_idx = z.b->array[z.i].border_idx;
	x.b->border_lr[x.b->border_par[x.i].i] = z;
	z.b->border_lr[z.i] = x.b->border_par[x.i];
      }
    }
    x.b->reached[x.i] = true;
    x = x.b->border_par[x.i];
    if(s.b != NULL) s.i = z.b->array[z.i].border_idx;
  }

  /* If s exists (node has  been added in the combined tree and may have been updated in later steps) */
  if(s.b != NULL){

    BorderIndex origin;
    s = b_levelroot(s);

  while (!(x.b->reached[x.i]) && !is_bottom(s)) {
    update_node_attribute(x,s);
    s = b_parent_lr(s);
    
    if (!is_bottom(s)) {
      origin = s.b->border_ori[s.i];
      if (x.b == origin.b) {
	x.b->border_par[x.i] = origin;
      } else if (!is_bottom(origin) && !is_bottom(origin.b->border_lr[origin.i])) {
	x.b->border_par[x.i] = origin.b->border_lr[origin.i];
      } else {      
	x.b = adding_node(x,s);
	x.b->array[x.b->size_curr-1].border_idx = s.i;	 
	if (!is_bottom(origin)) {
	  x.b->border_lr[x.b->size_curr-1] = origin;
	  origin.b->border_lr[origin.i] = (BorderIndex) {.b = x.b, .i = x.b->size_curr-1};
	}
	else
	  s.b->border_ori[s.i] =  (BorderIndex) {.b = x.b, .i = x.b->size_curr-1};
      }
      x.b->reached[x.i] = true;
      x = x.b->border_par[x.i];
    }
    else {
      x.b->border_par[x.i].i = BOTTOM;
      x.b->reached[x.i] = true;
    }   	  
   }
  }
  return x.b; 
}


void update_node_attribute(BorderIndex x, BorderIndex s) {
  /* Updating attribute of node x using the more recent node s */
  void *x_attr = NULL;
  if(x.b->attribute_idx[x.i] != BOTTOM )
    x_attr = (char *)  x.b->store->data + x.b->attribute_idx[x.i] * x.b->store->size_item;
  else
    x.b->attribute_idx[x.i] = x.b->store->item_curr;
  void *s_attr = (char *) s.b->store->data + s.b->attribute_idx[s.i] * s.b->store->size_item;
  clone_aux_data(x.b->store, &x_attr, s_attr);
}

Boundary *adding_node(BorderIndex x, BorderIndex s) {
  /* If x does not have a parent with the same level of s, s is added in the boundary tree of x */
  if (x.b->size_curr == x.b->size_alloc)
    x.b = realloc_b(x.b, x.b->size_curr != 1 ? 1.5*x.b->size_curr: 100*x.b->size_curr);
  void *attr = NULL;
  void *s_attr = (char *) s.b->store->data + s.b->attribute_idx[s.i] * s.b->store->size_item;
  x.b->array[x.b->size_curr] = s.b->array[s.i];
  clone_aux_data(x.b->store, &attr, s_attr);
  x.b->attribute_idx[x.b->size_curr] = x.b->store->item_curr - 1;
  x.b->reached[x.b->size_curr] =  false;
  x.b->border_par[x.i].b = x.b;
  x.b->border_par[x.i].i = x.b->size_curr;
  x.b->size_curr++;
  return x.b;
}
  
Boundary *combine(Boundary *a, Boundary *b, Direction d) {
  /* Given two boundary trees a and b, return the combined tree c */
  
  ulong  size_upb      	= a->size_curr + b->size_curr;
  ulong  size_attr_upb	= a->store->item_curr + b->store->item_curr;
  ulong  size_attrib	= b->store->size_item;  
  Boundary 	*c      = calloc(1, sizeof(Boundary));  	     check_alloc(c, 203);  
  c->array              = calloc(size_upb, sizeof(BoundaryNode));    check_alloc(c->array, 204);
  c->attribute_idx      = calloc(size_upb, sizeof(idx));	     check_alloc(c->attribute_idx, 205);
  c->border_par         = calloc(size_upb, sizeof(BorderIndex));     check_alloc(c->border_par, 206);
  c->border_ori         = calloc(size_upb, sizeof(BorderIndex));     check_alloc(c->border_ori, 207);
  c->store              = malloc(sizeof(AuxDataStore));              check_alloc(c->store, 208);
  
  init_aux_data_store(c->store, size_attrib, size_attr_upb);
  
  reset_border_idx(a);
  reset_border_idx(b);

  ulong s = 0;   /* points to entry point in c */
  ulong i, j, k; /* index local to either a or b */
  
  if (d == HORIZONTAL) { 

    ulong size_merge_idx   = a->offset[6] + b->offset[6] - 2*(a->offset[2]-a->offset[1]);
    c->merge_idx    	   = malloc(size_merge_idx * sizeof(ulong)); check_alloc(c->merge_idx, 209);
    c->offset[0]    	   = 0;
    
    /* Top border */
    if(a->offset[1] != a->offset[0]){
      for (j = 0; j < a->dims[2]; j++) {
	for (k = 0; k < a->dims[0] +  b->dims[0]; k++, s++) {
	  if (k < a->dims[0]) {	  	  
	    i = j * a->dims[0] + k;
	    b_add(c, a, s, a->offset[0] + i);
	  } else {
	    i = j * b->dims[0] + k - a->dims[0];
	    b_add(c, b, s, b->offset[0] + i);
	  } 
	}
      }
    } 
    c->offset[1] = s;

     /* Right border */
    for(i = b->offset[1]; i < b->offset[2]; i++, s++)
      b_add(c, b, s, i);
    c->offset[2] = s;

    /* Bottom border */
    if(a->offset[3] != a->offset[2]){
      for (j = 0; j < a->dims[2]; j++) {
	for (k = 0; k < a->dims[0] + b->dims[0]; k++, s++) {
	  if (k < a->dims[0]) {	  	  
	    i = j*a->dims[0] + k;
	    b_add(c, a, s, a->offset[2] + i);
	  } else {
	    i = j*b->dims[0] + k - a->dims[0];
	    b_add(c, b, s, b->offset[2] + i);
	  } 
	}
      }
    }
    c->offset[3] = s;

    /* Left border */
    for(i = a->offset[3]; i < a->offset[4]; i++, s++)
      b_add(c, a, s, i);
    c->offset[4] = s;

    /* Front border */
    if(a->offset[5] != a->offset[4]){
      for (j = 0; j < a->dims[1]; j++) {
	for (k = 0; k < a->dims[0] + b->dims[0]; k++, s++) {
	  if (k < a->dims[0]) {	  	  
	    i = j * a->dims[0] + k;
	    b_add(c, a, s, a->offset[4] + i);
	  } else {
	    i = j*b->dims[0] + k - a->dims[0];
	    b_add(c, b, s, b->offset[4] + i);
	  } 
	}
      }
    }
    c->offset[5] = s;

    /* Back border */
    if(a->offset[6] != a->offset[5]){
	for (j = 0; j < a->dims[1]; j++) {
	  for (k = 0; k < a->dims[0] + b->dims[0]; k++, s++) {
	    if (k < a->dims[0]) {	  	  
	      i = j*a->dims[0] + k;
	      b_add(c, a, s, a->offset[5] + i);
	    } else {
	      i = j*b->dims[0] + k - a->dims[0];
	      b_add(c, b, s, b->offset[5] + i);
	    } 
	  }
	}
    }
    
    c->offset[6] = s;
    c->dims[0]   = a->dims[0] + b->dims[0];
    c->dims[1]   = a->dims[1];
    c->dims[2]   = a->dims[2];


  } else if (d == VERTICAL) {
    
    ulong size_merge_idx  = a->offset[6] + b->offset[6] - 2*(a->offset[1]-a->offset[0]);
    c->merge_idx          = malloc(size_merge_idx * sizeof(ulong)); check_alloc(c->merge_idx, 210);
    c->offset[0]          = 0;
    
    /* Top border */
    for (i = a->offset[0]; i < a->offset[1]; i++, s++) 
      b_add(c, a, s, i);
    c->offset[1] = s;

    /* Right border */
    if(a->offset[2] != a->offset[1]){
      for (j = 0; j < a->dims[2]; j++) {
	for (k = 0; k < a->dims[1] + b->dims[1]; k++, s++) {
	  if (k < a->dims[1]) {	  	  
	    i = j*a->dims[1] + k;
	    b_add(c, a, s, a->offset[1] + i);
	  } else {
	    i = j*b->dims[1] + k - a->dims[1];
	    b_add(c, b, s, b->offset[1] + i);
	  } 
	}
      }
    }
    c->offset[2] = s; 

    /* Botoom border */
    for (i = b->offset[2]; i < b->offset[3]; i++, s++)  
      b_add(c, b, s, i);
    
    c->offset[3] = s;


    /* Left border */
    if(a->offset[4] != a->offset[3]){
      for (j = 0; j < a->dims[2]; j++) {
	for (k = 0; k < a->dims[1] + b->dims[1]; k++, s++) {
	  if (k < a->dims[1]) {	  	  
	    i = j*a->dims[1] + k;
	    b_add(c, a, s, a->offset[3] + i);
	  } else {
	    i = j * b->dims[1] + k - a->dims[1];
	    b_add(c, b, s, b->offset[3] + i);
	  } 
	}
      }
    }
    c->offset[4] = s;
    
    /* Front border */
    if(a->offset[5] != a->offset[4]){
      for (j = 0; j < a->dims[1] + b->dims[1]; j++) {
	for (k = 0; k < a->dims[0]; k++, s++) {
	  if (j < a->dims[1]) {	  	  
	    i = j*a->dims[0] + k;
	    b_add(c, a, s, a->offset[4] + i);
	  } else {
	    i = (j-a->dims[1])*(b->dims[0]) + k;
	    b_add(c, b, s, b->offset[4] + i);
	  } 
	}
      }
    }    
    c->offset[5] = s;

    /* Back border */
    if(a->offset[6] != a->offset[5]){
      for (j = 0; j < a->dims[1] + b->dims[1]; j++) {
	for (k = 0; k < a->dims[0]; k++, s++) {
	  if (j < a->dims[1]) {	  	  
	    i = j * a->dims[0] + k;
	    b_add(c, a, s, a->offset[5] + i);
	  } else {
	    i = (j - a->dims[1])*b->dims[0] + k;
	    b_add(c, b, s, b->offset[5] + i);
	  } 
	}
      }
    }
    
    c->offset[6] = s;
    c->dims[0] = a->dims[0];
    c->dims[1] = a->dims[1] + b->dims[1];
    c->dims[2] = a->dims[2];

  } else {
    
    ulong size_merge_idx  = a->offset[6] + b->offset[6] - 2*(a->offset[5]-a->offset[4]);
    c->merge_idx   = malloc(size_merge_idx * sizeof(ulong)); check_alloc(c->merge_idx, 211);
    c->offset[0]   = 0;

    /* Top border */
    for(i = a->offset[0]; i<a->offset[1]; i++, s++) 
      b_add(c, a, s, i);

    for(i = b->offset[0]; i<b->offset[1]; i++, s++) 
      b_add(c, b, s, i);    
    c->offset[1] = s;

    /* Right border */
    for(i = a->offset[1]; i<a->offset[2]; i++, s++) 
      b_add(c, a, s, i);
    
    for(i = b->offset[1]; i<b->offset[2]; i++, s++) 
      b_add(c, b, s, i); 
    c->offset[2] = s;

    /* Bottom border */
    for(i = a->offset[2]; i<a->offset[3]; i++, s++) 
      b_add(c, a, s, i);
    
    for(i = b->offset[2]; i<b->offset[3]; i++, s++) 
      b_add(c, b, s, i);   
    c->offset[3] = s;

    /* Left border */
    for(i = a->offset[3]; i<a->offset[4]; i++, s++) 
      b_add(c, a, s, i);
       
    for(i = b->offset[3]; i<b->offset[4]; i++, s++) 
      b_add(c, b, s, i);    
    c->offset[4] = s;

    /* Front border */
    for(i = a->offset[4]; i<a->offset[5]; i++, s++) 
      b_add(c, a, s, i);    
    c->offset[5] = s;

    /* Back border */
    for(i = b->offset[5]; i<b->offset[6]; i++, s++) 
      b_add(c, b, s, i);
    
    c->offset[6] = s;
    c->dims[0] = a->dims[0];
    c->dims[1] = a->dims[1];
    c->dims[2] = a->dims[2] + b->dims[2];
  }
  
  if(c->size_curr == 0) /* No more merge steps to be done */
    return NULL;
  
  s = c->size_curr;

  BorderIndex 	origin,	parent;
  ulong         curr, ex = s;
  idx 		c_idx;
  

  for (i = 0; i < ex; i++) {
    curr = i;
    origin = c->border_ori[curr];
    parent = b_levelroot(b_parent(origin));

    while (true) {
      if (is_bottom(parent)) {
        /* case 1: parent is bottom */
        c->border_par[curr] = (BorderIndex) {.b = c, .i = BOTTOM};
        break;
      }
      c_idx = b_node(parent).border_idx; /* index in c, if not BOTTOM */
      if (c_idx != BOTTOM) {
        /* case 2: parent is in c, avoid duplicates */
        c->border_par[curr] = (BorderIndex) {.b = c, .i = c_idx};	
	break;	     	
      } else {         
	/* case 3: add parent to c */
	(parent.b)->array[parent.i].border_idx = s; /* to avoid duplicates in case 2 */
	c->border_ori[s] = parent;
	c->array[s] = b_node(parent);
	void *attr = NULL;
	if((parent.b)->attribute_idx[parent.i] != BOTTOM){
	  void *b_attr = (char *) (parent.b)->store->data + (parent.b)->attribute_idx[parent.i]*size_attrib;
	  clone_aux_data(c->store, &attr, b_attr);
	  c->attribute_idx[s] = c->store->item_curr - 1;
	} else 
	  c->attribute_idx[s] = BOTTOM;
	c->border_par[curr] = (BorderIndex) {.b = c, .i = s};
	
	curr = s;
	s++;
      }
      parent = b_parent_lr(parent);
    }
  }

  c->size_curr     = c->size_alloc = c->size_init = s;
  c->array         = realloc(c->array,        s*sizeof(BoundaryNode)); check_alloc(c->array, 212);
  c->attribute_idx = realloc(c->attribute_idx,s*sizeof(idx));      check_alloc(c->attribute_idx, 213);
  c->border_par    = realloc(c->border_par,   s*sizeof(BorderIndex)); check_alloc(c->border_par, 214);
  c->border_ori    = realloc(c->border_ori,   s*sizeof(BorderIndex)); check_alloc(c->border_ori, 215);
  c->border_lr     = malloc(s * sizeof(BorderIndex));		      check_alloc(c->border_lr, 216);
  
  #pragma omp parallel for
  for (i = 0; i < s; ++i) 
    c->border_lr[i] = (BorderIndex)  {.b = c, .i = BOTTOM};

  return c;
}


void b_add(Boundary *c, Boundary *b, ulong s, ulong i) {
  /* Adding nodes from b in combined tree c */
  
  void *attr = NULL;
  ulong idx_b = b->merge_idx[i];
  BorderIndex toadd = (BorderIndex) {.b = b, .i = idx_b};
  toadd = b_levelroot(toadd);
  idx_b = toadd.i;
  b = toadd.b;
  
  if(b->array[idx_b].border_idx != BOTTOM)
    c->merge_idx[s] = b->array[idx_b].border_idx;
  else {
    ulong idx_c = c->size_curr++;
    b->array[idx_b].border_idx = idx_c;
    c->merge_idx[s]		= idx_c;
    c->array[idx_c] 		= b->array[idx_b];
    c->border_ori[idx_c]  	= (BorderIndex) {.b = b, .i = idx_b};
    
    if (b->attribute_idx[idx_b] != BOTTOM) {
      void *b_attr = (char *) b->store->data + b->attribute_idx[idx_b] * b->store->size_item;
      clone_aux_data(c->store, &attr, b_attr);
      c->attribute_idx[idx_c] = c->store->item_curr - 1;
    } else
      c->attribute_idx[idx_c] = BOTTOM;
  }
}


void merge(Boundary *a, Boundary *b, Direction d) {
  /* there are three cases:
     HORIZONTALLY: merge the right side of a with the left side of b 
     VERTICALLY  : merge the bottom of a with top of b, 
     DEPTH       : nerge the front border of b with the back border of a */

  ulong length;
  ulong offset;
  
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
  
  /* traverse border */
  BorderIndex 	x, y;
  value 	min_prev = 0;
  value 	min_curr = 1;
  bool		test_min = false;
  for (ulong c = 0; c < length; c++) {
    if(c %offset == 0) test_min = false;
    x = idx_i(a, c, d);  
    y = idx_j(b, c, d);
    min_curr = b_gval(x);
    if (!test_min || min_curr > min_prev )
      merge_b_nodes(x,y);
    min_prev = min_curr;
    test_min = true;
  }
} /* merge */


bool AreSame(value a, value b)
{
    return fabs(a - b) < 0.0000001;
}
bool AreGt(value a, value b)
{
    return a - b > 0.0000001;
}

void merge_b_nodes(BorderIndex x, BorderIndex y){
  BorderIndex   z, h;
  void 		*x_attr;
  void 		*y_attr;
  void 		*attr        = NULL;
  void 		*attr1       = NULL;
  ulong 	size_attrib  = x.b->store->size_item;
  value         xg,yg,zg;
  
  if (x.b->attribute_idx[x.i] == BOTTOM){ /* Switch x and y */
    h = x; x = y; y = h;
  }
  
  x = b_levelroot(x);
  y = b_levelroot(y);

  assert(bi_gval_equal(b_gval(x),b_gval(y))==0);

  while (!bi_equal(x, y) && !is_bottom(y) ) {

    x_attr = x.b->attribute_idx[x.i] != -1 ?
      (char *) x.b->store->data + x.b->attribute_idx[x.i] * size_attrib : NULL;
    y_attr = y.b->attribute_idx[y.i] != -1 ?
      (char *) y.b->store->data + y.b->attribute_idx[y.i] * size_attrib : NULL;

    xg = b_gval(x);
    yg = b_gval(y);

    z = b_parent_lr(x);
    zg = b_gval(z);

    if (!is_bottom(z) && zg >= yg) {

      if(attr){
	if(x_attr != NULL)
	  merge_aux_data(x_attr, attr);
	else{
	  x.b->attribute_idx[x.i] = x.b->store->item_curr;
	  clone_aux_data(x.b->store, &x_attr, attr);
	}
      }
      x = z;
    }  else {

      if (xg == yg) {
	if (!is_bottom(y.b->border_lr[y.i])) {

	  x.b->border_par[x.i] = x.b == y.b ? y : y.b->border_lr[y.i];
	  x.b->border_lr[x.i] = x.b == y.b ? y.b->border_lr[y.i] : y ;
	} else if (!is_bottom(x.b->border_lr[x.i])) {
	  if(y_attr != NULL) clone_aux_data(NULL, &attr, y_attr);
	  h = b_parent_lr(y);
	  y.b->border_par[y.i] = x.b == y.b ? x : x.b->border_lr[x.i];
	  y.b->border_lr[y.i] =  x.b == y.b ? x.b->border_lr[x.i] : x;
	  y = h;
	  continue;
	} else {
	  if(x.b != y.b) {
	    y.b->border_lr[y.i] = x;
	    x.b->border_lr[x.i] = y;
	  } 
	    
	  x.b->border_par[x.i] = y;
	}
      } else
	x.b->border_par[x.i] = y;

      if(x_attr){
	attr ? merge_to_aux_data(NULL, &attr1, x_attr, attr) : clone_aux_data(NULL, &attr1, x_attr);
	clone_aux_data(NULL, &attr, x_attr);
      } else{
      	delete_aux_data(attr);
	attr= NULL;
	x.b->attribute_idx[x.i]= x.b->store->item_curr;
      }
      if(attr1) clone_aux_data(x.b->store, &x_attr, attr1);
      x = y;
      y = z;
    }	
  }

  if (is_bottom(y) && attr) {
    while (!is_bottom(x)) {
      if(x.b->attribute_idx[x.i] != BOTTOM){
	x_attr = (char *) x.b->store->data + x.b->attribute_idx[x.i] * size_attrib;
	merge_aux_data(x_attr, attr);
      }
      else{
	x.b->attribute_idx[x.i] = x.b->store->item_curr;
	clone_aux_data(x.b->store, &x_attr, attr);
      }
      x = b_parent_lr(x);
    }
  }

  delete_aux_data(attr);
  delete_aux_data(attr1);
}

BorderIndex idx_i(Boundary *a, ulong c, Direction d) {
  /* Return node from a to be merged */
  
  BorderIndex bi;
  bi.b = a;

  if (d == HORIZONTAL) 
    bi.i = a->merge_idx[a->offset[1]+c];
  else if (d == VERTICAL) 
    bi.i = a->merge_idx[a->offset[2]+c];
  else 
    bi.i = a->merge_idx[a->offset[5]+c];
        
  return bi;
}

BorderIndex idx_j(Boundary *b, ulong c, Direction d) {
  /* Return node from b to be merged */

  BorderIndex bi;
  bi.b = b;
 
  if (d == HORIZONTAL) 
    bi.i = b->merge_idx[b->offset[3]+c];
  else if (d == VERTICAL) 
    bi.i = b->merge_idx[b->offset[0]+c];
  else 
    bi.i = b->merge_idx[b->offset[4]+c];
    
  return bi;
}
   
   
Boundary *create_boundary(Node* local_tree, ulong *dims, int attrib_choice, int connectivity){
  
  /*  The boundary tree is an array that includes the levelroots of the nodes located in the tile border, and their parents.  In case of 4 or 6 connectivity, we merge along the overlapping nodes with the smallest intensity. With 8 or 26 connectivity, we need to merge along the overlapping nodes with th highest intensity, as the case is more complex. */

  Boundary *b 		 = calloc(1, sizeof(Boundary));   check_alloc(b, 217);
  ulong     length[2]    = {0};
  long      increment[3] = {0};
  long      offset[3]    = {0};
  ulong     offset_bd 	 = 0;
  ulong	    size_border	 = 0;
  bool      *border 	 = local_tree->border;
  
  local_tree->border_idx = malloc(dims[0]*dims[1]*dims[2]*sizeof(idx));
  memset(local_tree->border_idx, -1, dims[0]*dims[1]*dims[2]*sizeof(idx));

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
  b->size_alloc     = 3 * size_border; /* Allocation is bigger to add node's parents */
  b->array          = malloc(b->size_alloc  * sizeof(BoundaryNode)); check_alloc(b->array, 218);
  b->merge_idx      = malloc(size_border    * sizeof(BorderIndex));  check_alloc(b->merge_idx, 219);
  b->attribute_idx  = malloc(b->size_alloc  * sizeof(idx));	     check_alloc(b->attribute_idx, 220);
  b->border_par     = malloc(b->size_alloc  * sizeof(BorderIndex));  check_alloc(b->border_par, 221);
  b->store          = calloc(1, sizeof(AuxDataStore));               check_alloc(b->store, 222);
    
  init_aux_data_store(b->store, AttribsArray[attrib_choice].size, size_border);

  //     0: Face top    //
  
  if(border[2]){
    increment[0] = 1;
    increment[1] = dims[0] * dims[1];
    increment[2] = dims[0];
    length[0]    = b->dims[0];
    length[1]    = b->dims[2];
    add_side(local_tree, b, 0, length, 0+offset[1], &offset_bd, increment, connectivity);
  } else
    b->offset[0] = offset_bd;

  //    1: Face right   //

  if(border[1]){
    increment[0] = dims[0];
    increment[1] = dims[0] * dims[1];
    increment[2] = -1;
    length[0]    = b->dims[1];
    length[1]    = b->dims[2];
    add_side(local_tree,    b, 1, length, dims[0]-1+offset[0], &offset_bd, increment, connectivity);
  } else
    b->offset[1] = offset_bd;
  
  //     2: Face bottom  //

  if(border[3]){
    increment[0] = 1;
    increment[1] = dims[0] * dims[1];
    increment[2] = -dims[0];
    length[0]    = b->dims[0];
    length[1]    = b->dims[2];
    add_side(local_tree,   b, 2, length, dims[0] * (dims[1]-1)+offset[1], &offset_bd, increment, connectivity);
  } else
    b->offset[2] = offset_bd;

  //    3: Face left    //
  
  if(border[0]){
    increment[0] = dims[0];
    increment[1] = dims[0] * dims[1];
    increment[2] = 1;
    length[0]    = b->dims[1];
    length[1]    = b->dims[2];
    add_side(local_tree,   b, 3, length,  0+offset[0], &offset_bd, increment, connectivity);
  } else
    b->offset[3] = offset_bd;

  //     4: Faxe front   //
  
  if(border[4]){
    increment[0] = 1;
    increment[1] = dims[0];
    increment[2] = dims[0]*dims[1];

    length[0]    = b->dims[0];
    length[1]    = b->dims[1];
    add_side(local_tree,   b, 4, length, 0+offset[2], &offset_bd, increment, connectivity);
  } else
    b->offset[4] = offset_bd;


  //     5: Face back    //

  if(border[5]){
    increment[0] = 1;
    increment[1] = dims[0];
    increment[2] = -dims[0]*dims[1];
    length[0]    = b->dims[0];
    length[1]    = b->dims[1];
    add_side(local_tree, b, 5, length, dims[0]*dims[1]*(dims[2] - 1)+offset[2], &offset_bd, increment, connectivity);
  } else {
    b->offset[5] = offset_bd;
  }

  b->offset[6]    = offset_bd ;

  // Parent QU //


  //     6: Parents      //

  b = add_ancestors(local_tree, b);

  // b->parent_qu =  add_parent_qu(local_tree, b);

  /*       Shrinking     */
  
  b->array           = realloc(b->array,           b->size_curr * sizeof(BoundaryNode));
  b->attribute_idx   = realloc(b->attribute_idx,   b->size_curr * sizeof(idx));
  b->border_par      = realloc(b->border_par,    b->size_curr * sizeof(BorderIndex));

  /* Allocating the remaining variables */
  
  b->border_lr 	     = malloc(b->size_curr * sizeof(BorderIndex)); check_alloc(b->border_lr, 223);

  #pragma omp parallel for 
  for (ulong i = 0; i < b->size_curr; ++i) 
    b->border_lr[i] = (BorderIndex) {.b = b, .i = BOTTOM};
  
  b->size_init = b->size_alloc = b->size_curr;
  free(local_tree->border_idx);
  return b;
}


void add_side(Node *local_tree, Boundary *b, int side, ulong *length, ulong offset, ulong *offset_bd, long *increment, int connectivity) {
  /* Adding nodes from local_tree to the boundary tree */
  b->offset[side] = *offset_bd;

  ulong s = 0, k,l,  indx;
  
  for(ulong j = 0; j < length[1]; j++) {
    for (ulong i = 0; i < length[0]; ++i, s++) {
      k = offset + (j * increment[1]) + (i * increment[0]);
      l = k + increment[2];

      if(connectivity >= 8 && local_tree->gval[l] > local_tree->gval[k]){	
	k = l;
      }  else if (connectivity < 8 && local_tree->gval[l] < local_tree->gval[k]){
	k = l;
      }
      
      k = get_levelroot(local_tree, k);

      if(local_tree->border_idx[k] != BOTTOM)
	b->merge_idx[*offset_bd+s] = local_tree->border_idx[k];
      else {
	indx 			   = b->size_curr++;
	local_tree->border_idx[k]  = indx;
	b->merge_idx[*offset_bd+s] = indx;
	node_to_bound(local_tree, b, indx, k);
	if (local_tree->attribute[k]) {
	  void *attr = NULL;
	  clone_aux_data(b->store, &attr, local_tree->attribute[k]);
	  b->attribute_idx[indx] = b->store->item_curr - 1;
	} else{
	  b->attribute_idx[indx] = BOTTOM;
	}
      }
    }
  }
  *offset_bd += s;    

} /* add_side */

/*idx *add_parent_qu(Node *tree, Boundary *b) {
  b->parent_qu = malloc(2*b->offset[6]*sizeof(idx));
  memset(b->parent_qu, -1, 2*b->offset[6]*sizeof(idx));
  // info("%d", b->size_curr);
  for (ulong i = 0; i < b->offset[6]; ++i) {
    ulong curr  = b->array[i].index;
    idx par_qu =  tree->parent_qu[curr];
    info("CUrr %d, parent_qu %d", curr, tree->parent_qu[curr]);

    while(par_qu != BOTTOM && b->parent_qu[tree->border_idx[curr]] == BOTTOM){
      ulong bpar = tree->border_idx[par_qu];
      info("Bpar %d", bpar);

      b->parent_qu[tree->border_idx[curr]] = bpar;
      curr = par_qu;
      par_qu = tree->parent_qu[par_qu];
      info("CUrr %d, parent_qu %d", curr, tree->parent_qu[curr]);

    }
    b->parent_qu[tree->border_idx[curr]] = BOTTOM;
  }
  return b->parent_qu;
  }*/

Boundary *add_ancestors(Node *local_tree, Boundary *b) {
  /* We have to keep track of all nodes that are ancestors which are not in the border */
  /* Note that in order to visit every node only once, we need to iterate over the border, and only add ancestors that are not in the border nor already in the map. Administration is done here. 
  
  The theoretical upper limit of the boundary including ancestors is:
  (G-1) * (L/2)
  where G is the number of gray levels, and L is the length of the border
  We do not allocate for this upper limit, as the number of nodes added to a boundary tree was found to be much smaller in practice 

  As the size of the ancestors map is dynamic and depending on the input image, we then have to make an assumption for a good initial size of the map. We here chose to make it 5 times bigger, and reallocate 1.5 its current size (1.5 is the most efficient factor)  when needed, which is O(n).
  */
  
  ulong origsize  = b->size_curr;
  idx   parent, bx, bx_par;
  ulong curr;
   
  for (ulong i = 0; i < origsize; ++i) {

    curr  = b->array[i].index; /* index in maxtree */
    while (true) {

      parent = get_parent(local_tree, curr); /* index in maxtree */

      if (parent == BOTTOM) {
	b->border_par[local_tree->border_idx[curr]] = (BorderIndex) {.b=b, .i = BOTTOM};
        break; /* next! */
      }

      bx_par = local_tree->border_idx[parent]; /* parent index in border */
      bx     = local_tree->border_idx[curr]; /* index in border */

      if (bx_par != BOTTOM) {
        /* parent is already in the border */
        b->border_par[bx] = (BorderIndex) {.b=b, .i = bx_par};	
	  break;
	}
	        
      else {
	if (b->size_curr == b->size_alloc) {
	  b->size_alloc      = b->size_alloc > 1 ? 1.5*b->size_alloc : 100;
	  b->array           = realloc(b->array,           b->size_alloc * sizeof(BoundaryNode));
	  b->attribute_idx   = realloc(b->attribute_idx,   b->size_alloc * sizeof(idx));
	  b->border_par    = realloc(b->border_par,    b->size_alloc * sizeof(BorderIndex));
	  check_alloc(b->array,         224);
	  check_alloc(b->attribute_idx, 225);
	  check_alloc(b->border_par,  226);
	}

        /* add parent to border */
	void *attr = NULL;
        local_tree->border_idx[parent] = b->size_curr;
	node_to_bound(local_tree, b, b->size_curr, parent);
	if(local_tree->attribute[parent]){
	  clone_aux_data(b->store, &attr, local_tree->attribute[parent]);
	  b->attribute_idx[b->size_curr] = b->store->item_curr - 1;
	} else
	  b->attribute_idx[b->size_curr] = BOTTOM;
	
        /* set border_par */
        b->border_par[bx] = (BorderIndex) {.b = b, .i = local_tree->border_idx[parent]};
        b->size_curr++; /* move to next empty index */
        curr = parent; /* next! */

      }
    }
  }
  return b;
} /* add_ancestors */

void node_to_bound(Node *local_tree, Boundary *b, ulong b_idx, ulong tree_idx){
  b->array[b_idx].index         = tree_idx;
  b->array[b_idx].gval          = local_tree->gval[tree_idx];
  b->array[b_idx].border_idx    = local_tree->border_idx[tree_idx];
  b->attribute_idx[b_idx]       = BOTTOM;
} /* node_to_bound */

void bound_to_tree(Node *local_tree, Boundary *b, ulong b_idx, ulong tree_idx){
  local_tree->gval[tree_idx]    = b->array[b_idx].gval;
  local_tree->parent[tree_idx]     = BOTTOM;
  local_tree->attribute[tree_idx]  = 0;
} /* bound_to_tree */



Boundary *realloc_b(Boundary *b, ulong size_tree_new){

  //  debug("REALLOC, oldsize %d size_tree_new %d", b->size_curr, size_tree_new);
 
  b->array         = realloc(b->array,          size_tree_new * sizeof(BoundaryNode));
  b->attribute_idx = realloc(b->attribute_idx,  size_tree_new * sizeof(idx));
  b->border_par    = realloc(b->border_par,     size_tree_new * sizeof(BorderIndex));
  b->border_lr     = realloc(b->border_lr,      size_tree_new * sizeof(BorderIndex));
  b->border_ori    = realloc(b->border_ori,     size_tree_new * sizeof(BorderIndex));
  b->reached       = realloc(b->reached,        size_tree_new * sizeof(bool));

  check_alloc(b->array, 227);
  check_alloc(b->attribute_idx, 228);
  check_alloc(b->border_par, 229);
  check_alloc(b->border_lr, 230);
  check_alloc(b->border_ori, 231);
  check_alloc(b->reached, 232);

  if (size_tree_new > b->size_curr) {
    #pragma omp parallel for 
    for (ulong i = b->size_curr; i < size_tree_new; i++)
      b->border_lr[i] = b->border_ori[i] = (BorderIndex) {.b = b, .i = BOTTOM};    
  }
 
  b->size_alloc = size_tree_new;
  return b;
}


void reset_border_idx(Boundary *b) {
  #pragma omp parallel for 
  for (ulong i = 0; i < b->size_curr; i++) 
    b->array[i].border_idx = BOTTOM;
} /* reset_border_idx */


BorderIndex b_levelroot(BorderIndex bi) {
  BorderIndex ri = bi;
  if (is_bottom(ri)) 
    return bi; /* still BOTTOM */
  else{
    value gv = b_gval(bi);

    while (!is_bottom(b_parent(ri)) && (gv == b_gval(b_parent(ri)))) {
      ri = b_parent(ri);
      if (bi_equal(ri, b_parent(ri))){
	error("[CONNEXION ERROR] Next levelroot %ld in %d");
	MPI_Abort(MPI_COMM_WORLD, 240);
      }
    }
    return ri;
  }
} /* b_levelroot */


BorderIndex b_id_levelroot(BorderIndex bi) {
  if (is_bottom(bi)) 
    return bi; /* still BOTTOM */
    
  BorderIndex ri = bi;
  BorderIndex pi = bi;  
  value       gv = b_gval(bi);

  while ( !is_bottom(b_parent(pi)) && gv == b_gval(b_parent(pi)) && b_parent(pi).b == bi.b ) {
    pi = b_parent(pi);
    ri = pi;
  }

  return ri;
} /* b_id_levelroot */


bool is_bottom(BorderIndex bi) {
  return bi.i == BOTTOM;
} /* is_bottom */
 
bool bi_equal(BorderIndex ai, BorderIndex bi) {
  return (ai.b == bi.b) && (ai.i == bi.i);
} /* bi_equal */

int compare_floats (const void *a, const void *b)
{
  const float *da = (const float *) a;
  const float *db = (const float *) b;

  return (*da > *db) - (*da < *db);
}

int bi_gval_equal(value a, value b) {
  //info("%f, %f",a,b);
  if(FLOAT_TYPE){
    value_t at = transform(a);
    value_t bt = transform(b);
    if(at < bt) return -1;
    else if(at == bt) return 0;
    else  return 1;
  } 
  if(a < b) return -1;
  else if(a == b) return 0;
  else  return 1;
} /* bi_equal */

BorderIndex b_parent_lr(BorderIndex bi) {
  /* get parent index */
  BorderIndex pi = b_parent(bi);
  /* get levelroot of parent */
  BorderIndex ri = b_levelroot(pi);
  // trace("b_parent_lr pi: %p %ld ri: %p %ld", pi.i, (void *)pi.b, ri.i, (void *)ri.b);
  return ri;
} /* b_parent_lr */

BorderIndex b_parent(BorderIndex bi) {
  BorderIndex pi = bi.b->border_par[bi.i];
  return pi;
} /* b_parent */

BoundaryNode b_node(BorderIndex bi) {
  if (bi.b == NULL) {
    error("asking node of non-initilized border %ld", bi.i);
  }
  return (bi.b)->array[bi.i];
}


value b_gval(BorderIndex bi) {
  if (bi.b == NULL) {
    error("asking gval of non-initilized border");
    MPI_Abort(MPI_COMM_WORLD, 241);
  }
  if ((ulong)bi.i > ((bi.b)->size_curr)) {
    warn("asking gval of index (%ld) higher than size (%ld) of border", bi.i, (bi.b)->size_curr);
    return BOTTOM;
  }

 return ((bi.b)->array[bi.i]).gval;
} /* b_gval */


void free_boundary(Boundary *b) {
  if(b != NULL){
    free(b->array);
    free(b->attribute_idx);
    free(b->merge_idx);
    free(b->border_par);
    free(b->border_lr);
    free(b->border_ori);
    clear_aux_data_store(b->store);    
  }
  free(b);
} /* free_boundary */

