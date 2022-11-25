#include "types.h"
#include "attributes.h"
#include "lambdavec.h"
#include "flood.h"
#include "filter.h"

void tree_filter_min(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);
void tree_filter_direct(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);
void tree_filter_max(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);
void tree_filter_subtractive(Node* tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda);
void node_differential_seg(Node *tree, LambdaVec *lvec, value *out_dh, value *temp_dh, value *out_orig, value *out_scale, value *temp_scale, bool *temp_valid, ulong lwb, ulong upb, ulong current,  value *maxDH, value *curDH, value *maxOrig,  int *maxScale, int *curScale,  double (*attribute)(void *));
int find_scale(LambdaVec *lvec, double attribute);
int find_scale_csl(LambdaVec *lvec, double attribute);
/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	     Tree Filtering        */
/*				   */
/* +++++++++++++++++++++++++++++++ */

void tree_filtering(Node *tree, value *out, ulong size_tile,  int decision_choice, int attrib_choice, double lambda){
  
  bool   *reached    = calloc(tree->size_curr, sizeof(bool)); check_alloc(reached, 600);
  
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

void tree_filter_test(Node *tree, float *out, ulong lwb, ulong upb, double (*attribute)(void *)) {
  for (ulong v = lwb; v < upb; v++) {
    ulong v_lr = get_levelroot(tree, v);
    //info("%ld",v );
    out[v] = (float) (*attribute)(tree->attribute + v_lr*tree->size_attr);
    //  info("%f",out[v]);

  }
} /* tree_filter_min */


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


void tree_filter_direct(Node *tree, value *out, bool *reached, ulong lwb, ulong upb, double (*attribute)(void *), double lambda) {
  value  val;
  idx	 parent;
  ulong  v, u, w;

  for (v = lwb; v < upb; v++) {
    if (!reached[v]) { /* not filtered yet */
      w = v;
      parent = tree->parent[w];
      /* repeat while we're not at the bottom, th */
      while ((parent != BOTTOM) && (!reached[w]) && ((tree->gval[w] == tree->gval[parent] )  ||  ( (*attribute)(tree->attribute + w*tree->size_attr) < lambda))) {
        w = parent;
        parent = tree->parent[w];	      
      }
     
      if (reached[w]) {
        /* criterion satisfied at level tree[w].filter */
	val = out[w];
      } else if ( (*attribute)(tree->attribute + w*tree->size_attr) >= lambda) {
        /* w satisfies criterion */
        val = tree->gval[w]; 
      } else {
        /* criterion cannot be satisfied */
        val = 0;
      }

      /* set filt along par-path from v to w */
      u = v;
      while (u != w) {
	if ((lwb <= u) && (upb > u)){
	  out[u] 	= val;
	  reached[u] 	= true;
	}
        u = tree->parent[u];
      }
      if ( (lwb <= w) && (upb > w)){
	out[w] 		= val;
	reached[w]      = true;
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
/*           Diff Profil           */
/*				   */
/* +++++++++++++++++++++++++++++++ */


void tree_differential(Node *tree, ulong size_tile, LambdaVec *lvec, value *out_dh, value *temp_dh, value *out_orig, value *out_scale, value *temp_scale, bool *temp_valid,  double (*attribute)(void *)){
;
#pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int id	   = omp_get_thread_num();
    ulong lwb 	   = id*size_tile/np_threads;
    ulong upb	   = (id+1)*size_tile/np_threads;
    ulong v;
    value curDH, maxDH, maxOrig;
    int maxScale, curScale;

    for (v=lwb; v<upb; v++) {
      if (!temp_valid[v]) {
	node_differential_seg(tree, lvec, out_dh, temp_dh, out_orig, out_scale, temp_scale, temp_valid, lwb, upb, v, &maxDH, &curDH, &maxOrig, &maxScale, &curScale, attribute);
      }
    }
  }
}

void combine_results(Node *tree, ulong size, LambdaVec *lvec, value *out_dh, value *out_dh2, value *out_orig, value *out_orig2, value *out_scale, value *out_scale2) {

  #pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int   id	    = omp_get_thread_num();
    ulong lwb 	    = id*size/np_threads;
    ulong upb	    = (id+1)*size/np_threads;
    int   numscales = lvec->num_lambdas;
    for (ulong i = lwb; i < upb; i++){

      if (out_dh2[i] > out_dh[i]) { /* maxDH of opening > maxDH of closing instance */
	out_scale[i]   = (value) numscales + 1 + out_scale2[i];
	out_dh[i]      = out_dh2[i];
	out_orig[i]    = out_orig2[i];
      } else if (out_dh2[i] < out_dh[i]) {
	out_orig[i]    = g_max_greyval - out_orig[i];
	out_scale[i]++;          
      } else { /* equality */
	out_scale[i]   = 0;
	out_orig[i]    = (out_dh[i] != 0) * (g_max_greyval - tree->gval[i]);
      }
    
    }
  }
}

void node_differential_seg(Node *tree, LambdaVec *lvec, value *out_dh, value *temp_dh, value *out_orig, value *out_scale, value *temp_scale, bool *temp_valid, ulong lwb, ulong upb, ulong current,  value *maxDH, value *curDH, value *maxOrig,  int *maxScale, int *curScale,  double (*attribute)(void *)){
   
  int scale = -1, numscales = lvec->num_lambdas; 
  idx parent;
  value DH = 0;

  if (is_levelroot(tree, current)){

    scale = find_scale_csl(lvec, (*attribute)(tree->attribute + current*tree->size_attr));
    DH = ((tree->parent[current] != BOTTOM) && (scale < numscales)) ?
      tree->gval[current] - tree->gval[tree->parent[current]] : 0;

  }
  if (is_levelroot(tree, current) && ((scale == numscales) || (tree->parent[current] == BOTTOM))){

    if (tree->parent[current] != BOTTOM) {
      *maxScale = numscales;
      *maxOrig  = 0;
    } else {
      *maxScale = scale;
      *maxOrig  = tree->gval[current];
    }

    *maxDH    = 0;
    *curDH    = 0;
    *curScale = *maxScale;

  } else {

    parent = tree->parent[current];

    if (!temp_valid[parent]) { 
      /* go into recursion to set parent values correctly */

      node_differential_seg(tree, lvec, out_dh, temp_dh, out_orig, out_scale, temp_scale, temp_valid, lwb, upb, parent, maxDH, curDH, maxOrig, maxScale, curScale, attribute);

    } else { /* if the parent is valid, copy relevant values */

      *maxScale = (int) out_scale[parent];
      *maxDH    = out_dh[parent];
      *maxOrig  = out_orig[parent];
      *curScale = (int) temp_scale[parent];
      *curDH    = temp_dh[parent];
    } 
    
    if (is_levelroot(tree, current)) {  
      /* if I have a level root, some things might change */

      if (scale == *curScale) {
        /* parent's area is in same scale class,  add current pixel's curDH */
	*curDH += DH;
      } else {
	/* at scale class change, update current scale and DH */
	*curDH = DH;
	*curScale = scale;
      }
      
      if (*curDH >= *maxDH) {
	/* If updated curDH is higher than or equal to the maximum DH found
	   update maxDH, maxScale, and outOrig */
	*maxDH    = *curDH;
	*maxScale = *curScale;
	*maxOrig  = tree->gval[current];
      }       
    }  
  }

  if ((current >= lwb) && (current < upb)) {
    out_scale[current]  = (value) *maxScale;
    out_dh[current]     = *maxDH;
    out_orig[current]   = *maxOrig;      
    temp_scale[current] = (value) *curScale;
    temp_dh[current]    = *curDH;
    temp_valid[current] = true;
    //  info("%d, %d, %d, %d", current,  out_scale[current],out_dh[current] ,out_orig[current]);
  } 
  return;
} /* node_differential_seg */


/* +++++++++++++++++++++++++++++++ */
/*				   */
/*        Pattern Spectra          */
/*				   */
/* +++++++++++++++++++++++++++++++ */



void tree_pattern_spectrum(Node *tree,  ulong size, LambdaVec *lvec, double *copy_attr, value *gvals_par, double* spectrum, int background, double (*area)(void *), double (*attribute)(void *)){

  int numscales = lvec->num_lambdas;
  double *spectrum_th;

  #pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int   id	   = omp_get_thread_num();
    int scale; 
    ulong u,v;
    idx parent;
    value dh;
    double private;

    #pragma omp single
    spectrum_th = calloc(numscales*np_threads, sizeof(double)); check_alloc(spectrum_th, 601);
    if(np() > 1){
      #pragma omp for
      for (v = 0; v < size; v++) {
	if ((copy_attr[v] != -DBL_MAX)  && (tree->parent[v] != BOTTOM || background)) {
	  scale = find_scale(lvec, (*attribute)(tree->attribute + get_levelroot(tree, v)*tree->size_attr));
	  if (scale < numscales) {
	    parent = get_levelroot(tree, tree->parent[v]);
	    private = copy_attr[v];
	    dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[parent]): tree->gval[v];
	    spectrum_th[(numscales)*id + scale] += dh * private;	
	    u = v;
	    if(tree->parent[v] != BOTTOM){
	      while ((tree->parent[parent] != BOTTOM || background) && (gvals_par[v] < tree->gval[parent])) {    
		u = parent;
		scale = find_scale(lvec, (*attribute)(tree->attribute + u*tree->size_attr));
		if (scale>=numscales)
		  break;
		parent = get_levelroot(tree, tree->parent[u]);
		dh = parent != BOTTOM ? (tree->gval[u] - tree->gval[parent]): tree->gval[u];
		spectrum_th[(numscales)*id + scale] += dh * private;
		if(parent == BOTTOM)
		  break;
	      }
	    }
	  }
	}
      }
    } else {
      #pragma omp for
      for (v = 0; v < size; v++) {
	if(is_levelroot(tree, v) && (tree->parent[v] != BOTTOM || background)){
	  dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[tree->parent[v]]): tree->gval[v];
	  if (dh) {
	    scale = find_scale(lvec, (*attribute)(tree->attribute + v*tree->size_attr));
	    if (scale < numscales) 
	      spectrum_th[(numscales)*id + scale] += dh * (*area)(tree->attribute + v*tree->size_attr);
	  }
	}
      }
    }

    #pragma omp for
    for(int i=0; i < numscales; i++) {
      for(int t=0; t < np_threads; t++) {
	spectrum[i] += spectrum_th[t*numscales + i];
      }
    }

  }
  free(spectrum_th);
}/* tree_pattern_spectrum */

void tree_pattern_spectrum2d(Node *tree,  ulong size, LambdaVec *lvec_attr1,LambdaVec *lvec_attr2, double *copy_attr, value *gvals_par, double* spectrum, int background,  double (*area)(void *), double (*attribute1)(void *), double (*attribute2)(void *)){

  int numscales_attr1 = lvec_attr1->num_lambdas;
  int numscales_attr2 = lvec_attr2->num_lambdas;

  double *spectrum_th;

  #pragma omp parallel
  {
    int np_threads = omp_get_num_threads();
    int   id	   = omp_get_thread_num();
    int scale_attr1, scale_attr2; 
    ulong u,v;
    idx parent;
    value dh;
    double private;


    #pragma omp single
    spectrum_th = (double **) calloc(numscales_attr1*numscales_attr2*np_threads, sizeof(double));
    check_alloc(spectrum_th, 601);

    if(np() > 1){
      // To implement correctly with the good copy_attr
      /*
      #pragma omp for
      for (v = 0; v < size; v++) {
	if ((copy_attr[v] != -DBL_MAX)  && (tree->parent[v] != BOTTOM || background)) {
	  scale = find_scale(lvec, (*attribute)(tree->attribute + get_levelroot(tree, v)*tree->size_attr));
	  if (scale < numscales) {
	    parent = get_levelroot(tree, tree->parent[v]);
	    private = copy_attr[v];
	    dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[parent]): tree->gval[v];
	    spectrum_th[(numscales)*id + scale] += dh * private;	
	    u = v;
	    if(tree->parent[v] != BOTTOM){
	      while ((tree->parent[parent] != BOTTOM || background) && (gvals_par[v] < tree->gval[parent])) {    
		u = parent;
		scale = find_scale(lvec, (*attribute)(tree->attribute + u*tree->size_attr));
		if (scale>=numscales)
		  break;
		parent = get_levelroot(tree, tree->parent[u]);
		dh = parent != BOTTOM ? (tree->gval[u] - tree->gval[parent]): tree->gval[u];
		spectrum_th[(numscales)*id + scale] += dh * private;
		if(parent == BOTTOM)
		  break;
	      }
	    }
	  }
	}
      }*/
      error("This is not implemented yet!");
    } else {
      #pragma omp for
      for (v = 0; v < size; v++) {
	if(is_levelroot(tree, v) && (tree->parent[v] != BOTTOM || background)){
	  dh = tree->parent[v] != BOTTOM ? (tree->gval[v] - tree->gval[tree->parent[v]]): tree->gval[v];
	  if (dh) {
	    scale_attr1 = find_scale(lvec_attr1, (*attribute1)(tree->attribute + v*tree->size_attr));
	    scale_attr2 = find_scale(lvec_attr2, (*attribute2)(tree->attribute + v*tree->size_attr));
	    if (scale_attr1 < numscales_attr1 & scale_attr2 < numscales_attr2) 
	      spectrum_th[(numscales_attr1*numscales_attr2)*id + scale_attr1 + scale_attr2*numscales_attr1] +=
		dh * (*area)(tree->attribute + v*tree->size_attr);
	  }
	}
      }
    }

    #pragma omp for
    for(int i=0; i < numscales_attr1*numscales_attr2; i++) {
      for(int t=0; t < np_threads; t++) {
	spectrum[i] += spectrum_th[t*numscales_attr1*numscales_attr2 + i];
      }
    }

  }
  free(spectrum_th);
}/* tree_pattern_spectrum */
/* +++++++++++++++++++++++++++++++ */
/*				   */
/*         Annexes Functions       */
/*				   */
/* +++++++++++++++++++++++++++++++ */


int find_scale(LambdaVec *lvec, double attribute){
  int upper = lvec->num_lambdas-1, lower = 0, mid;

  if (attribute >=  (double) lvec->lambdas[upper])
    return upper;


  mid = (upper + lower) / 2;
  while (mid!=lower) {
    if(attribute >= (double) lvec->lambdas[mid])
      lower = mid;
    else
      upper = mid;

    mid = (upper + lower) / 2;
  }
  return lower;
} /* find_scale */

int find_scale_csl(LambdaVec *lvec, double attribute){
  int upper = lvec->num_lambdas-1, lower = 0, mid;

  if (attribute >=  (double) lvec->lambdas[upper])
    return upper+1;


  mid = (upper + lower) / 2;
  while (mid!=lower) {
    if(attribute >= (double) lvec->lambdas[mid])
      lower = mid;
    else
      upper = mid;

    mid = (upper + lower) / 2;
  }
  return lower;
} /* find_scale */


DecisionStruct Decisions[NUMDECISIONS] =
{
  {"Direct", tree_filter_direct},
  {"Min", tree_filter_min},  
  {"Max", tree_filter_max},
  {"Subtractive", tree_filter_subtractive},
};
