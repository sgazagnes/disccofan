
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