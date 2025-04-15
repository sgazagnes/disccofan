#include "types.h"
#include "attributes.h"
#include "lambdavec.h"
#include "flood.h"
#include "filtering.h"

void tree_filter_min(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);
void tree_filter_direct(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);
void tree_filter_max(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);
void tree_filter_subtractive(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Tree Filtering        */
/*				   */
/* +++++++++++++++++++++++++++++++ */

void tree_filtering(Node *tree, value *out, ulong size_tile,  int decision_choice, int attrib_choice, double lambda){
  
  bool   *reached    = calloc(tree->size_curr, sizeof(bool)); check_alloc(reached, 600);
  
  debug("Applying decision pruning %s", Decisions[decision_choice].name);
  #pragma omp parallel 
  {
    int np_threads = omp_get_num_threads();

    int id	   = omp_get_thread_num();
    ulong lwb 	   = id*size_tile/np_threads;
    ulong upb	   = (id+1)*size_tile/np_threads;
   
    Decisions[decision_choice].filter(tree, out, reached, lwb, upb, AttribsArray[attrib_choice].attribute, lambda);
  }
  if (decision_choice == 2) {
    for (ulong i = 0; i < size_tile; i++) {    
      if (!reached[i]) {
	idx im = (idx) i;
	while ((im != BOTTOM) && (out[im] != tree->gval[im])) {
	  out[im]  	  = tree->gval[im];
	  im              = tree->parent[im];
	}
      } 
    }
  }
  free(reached);
}

void tree_attribute_check(Node *tree, value *out, ulong lwb, ulong upb, double (*attribute)(void *)) {
  for (ulong v = lwb; v < upb; v++) {
    ulong v_lr = get_levelroot(tree, v);
    out[v] = (value) (*attribute)(tree->attribute + v_lr*tree->size_attr);
  }
} /* tree_attribute_check */


void tree_filter_min(Node *tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda) {
  value  val;
  idx	 parent;
  ulong  v, u, w, r;

  for (v = lwb; v < upb; v++) {
    if (!reached[v]) { /* not filtered yet */
      w = get_levelroot(tree, v);
      r = w;
      parent = get_parent(tree, w);
      /* repeat while we're not at the bottom, th */
      while ((parent != BOTTOM) && (!reached[w])) {
	if ((*attribute)(tree->attribute + w*tree->size_attr) < lambda)  r=w;
        w = parent;
        parent = get_parent(tree, w);
      }
      
      if (reached[w] && (out[w] != tree->gval[w])) {
        /* criterion satisfied at level tree[w].filter */
        val = out[w];
	r = w;
      } else if ((*attribute)(tree->attribute +r*tree->size_attr) >= lambda) {
        /* w satisfies criterion */
        val = tree->gval[r]; 
      } else if ((((*attribute)(tree->attribute + r*tree->size_attr) < lambda)) && (tree->parent[r] != BOTTOM)) {
	val = tree->gval[tree->parent[r]]; 
      }else {
        /* criterion cannot be satisfied */
        val = 0;
      }
      /* set filt along par-path from v to w */
      u = v;
      while (u != r) {
	if ((lwb <= u) && (upb > u)){
	  reached[u] = true;
	  out[u] = val;
	}
        u = tree->parent[u];
      }
      if ((lwb <= r) && (upb > r)){
	reached[r] = true;
	out[r] = val;
      }
    }
  }
} /* tree_filter_min */

void tree_filter_direct(Node *tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda)
{
  value val;
  idx parent;
  ulong v, u, w;
  for (v = lwb; v < upb; v++)
  {
    if (!reached[v])
    { /* not filtered yet */
      w = v;
      parent = tree->parent[w];
      /* repeat while we're not at the bottom, th */
      while ((parent != BOTTOM) && (!reached[w]) && ((tree->gval[w] == tree->gval[parent]) || ((*attribute)(tree->attribute + w * tree->size_attr) < lambda)))
      {
        w = parent;
        parent = tree->parent[w];
      }
      //info("%lf", (*attribute)(tree->attribute + w * tree->size_attr) );
      if (reached[w])
      {
        /* criterion satisfied at level tree[w].filter */
        val = out[w];
      }
      else if ((*attribute)(tree->attribute + w * tree->size_attr) >= lambda)
      {
        /* w satisfies criterion */
        val = tree->gval[w];
      }
      else
      {
        /* criterion cannot be satisfied */
        val = 0;
      }

      /* set filt along par-path from v to w */
      u = v;
      while (u != w)
      {
        if ((lwb <= u) && (upb > u))
        {
          out[u] = val;
          reached[u] = true;
        }
        u = tree->parent[u];
      }
      if ((lwb <= w) && (upb > w))
      {
        out[w] = val;
        reached[w] = true;
      }
    }
  }
} /* tree_filter_direct */

void tree_filter_max(Node *tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda) {
  value  val;
  idx	 u, parent;
  ulong  v,  w;

  for (v = lwb; v < upb; v++) {
    if (!reached[v]) { /* not filtered yet */
      w = v;
      parent = tree->parent[w];
      /* repeat while we're not at the bottom, th */
      while ((parent != BOTTOM) && (!reached[w]) && ((tree->gval[w] == tree->gval[parent] ) ||  ((*attribute)(tree->attribute + w*tree->size_attr) < lambda))) {
        w = parent;
        parent = tree->parent[w];	      
      }
      if (reached[w]) {
        /* criterion satisfied at level tree[w].filter */
	val = out[w];
      } else if ((*attribute)(tree->attribute + w*tree->size_attr) >= lambda) {
        /* w satisfies criterion */
        val = tree->gval[w]; 
	u = parent;
	while ((u != BOTTOM) && (out[u] != tree->gval[u])) {
	  if ((lwb <= (ulong) u) && (upb > (ulong) u)){
	    out[u] = tree->gval[u];
	    reached[u] = true;
	  } else {   /* else let thread 0 do the remainning work */
	    reached[u] = false;
	    break;
	  }      
	  u = tree->parent[u];
	}
      } else {
        /* criterion cannot be satisfied */
        val = 0;
      }

      /* set filt along par-path from v to w */
      u = v;
      while ((ulong) u != w) {
	if ((lwb <= (ulong) u) && (upb > (ulong) u)){
	  out[u] = val;
	  reached[u] = true;
	}
        u = tree->parent[u];

      }
      if ((lwb <=  w) && (upb >  w)){
	out[w] = val;
	reached[w] = true;
      }
    }
  }
} /* tree_filter_max */


void tree_filter_subtractive(Node *tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda) {

  ulong	v;
  idx u;
  idx w;
  ulong levpath[g_max_levels];
  value	filter[g_max_levels];
  idx 	i, j;
  
  for (v = lwb; v <  upb; v++) {
    if (!reached[v]) { /* not filtered yet */
      w = get_levelroot(tree, v);
      i = 0;
      while ((w != BOTTOM) && (!reached[w])) {
	levpath[i++] = w; 
	w = get_parent(tree, w);
      } 
      levpath[i] = w;
      if (w != BOTTOM)
	filter[i] = out[w];
      for (j=i-1; j>=0; j--) {
	u = levpath[j];
	if (w != BOTTOM) {
	  if (((*attribute)(tree->attribute + u*tree->size_attr)) < lambda)
	    filter[j] = filter[j+1];
	  else
	    filter[j] = tree->gval[u] + filter[j+1] - tree->gval[w];
	} else if (((*attribute)(tree->attribute + u*tree->size_attr)) >= lambda)
	  filter[j] = tree->gval[u];
	else
	  filter[j]=0; 
	w = (idx) u;
      }
      u = v;
      for (j=0; j<=i; j++) {
	while ((ulong) u != levpath[j]) {
	  if ((lwb <= (ulong) u) && (upb > (ulong) u)){
	    out[u] = filter[j];
	    reached[u] = true;
	  }
	  u = tree->parent[u];
	}
	if ((lwb <= (ulong) u) && (upb > (ulong) u)){
	  out[u] = filter[j];
	  reached[u] = true;
	}
	if (u != BOTTOM) u = tree->parent[u];
      }
    }
  }
} /* tree_filter_subtractive */


/* +++++++++++++++++++++++++++++++ */
/*				   */
/*         Annexes Functions       */
/*				   */
/* +++++++++++++++++++++++++++++++ */


DecisionStruct Decisions[NUMDECISIONS] =
{
  {"Direct", tree_filter_direct},
  {"Min", tree_filter_min},  
  {"Max", tree_filter_max},
  {"Subtractive", tree_filter_subtractive},
};
