#include "types.h"
#include "flood.h"
#include "workspace.h"
#include "queue.h"

/* +++++++++++++++++++++++++++++++ */
/*				   */
/*     	 Workspace Function        */
/*				   */
/* +++++++++++++++++++++++++++++++ */


size_t workspace_create_mappings(size_t len){
  const int num_digits = (sizeof(value) * CHAR_BIT + MXT_HISTO_SZ_LOG2 - 1) / MXT_HISTO_SZ_LOG2;
    
  if (num_digits == 1) return len * 2 * sizeof(ulong);
    
  return len * 2 * sizeof(SortItem);
}

size_t workspace_construct(ulong max_rank){


  ulong num_words = (max_rank + bits_per_word_log2() - 1) / bits_per_word_log2();
  size_t size_bitarr = num_words * sizeof(ulong);

  int nlevels = num_levels(max_rank);
  ulong count = nlevels;
      
  for (int i = nlevels; i--;)
    {
      max_rank >>= bits_per_word_log2();        
      count += max_rank + 1;
    }
      
  size_t size_queue = count * sizeof(ulong);

  return size_bitarr + size_queue;
}



size_t workspace_size(ulong size_base) 
{
  //size_t size_base     =  dims[0]*dims[1]*dims[2];
  //size_t inc           =  dims[2] == 1 ? dims[0]: dims[0]*dims[1];
  //size_t addsize       =  2*(args->threads_arg-1)*inc;
  //  size_t size_max      =  size_base;// + addsize;
  // size_t space_values  =  sizeof(value) * size_max;
  size_t space_mapping =  workspace_create_mappings(size_base);  
  size_t space_construction = 2 * sizeof(ulong) * size_base + workspace_construct(size_base);
  size_t max_space = MAX(space_construction, space_mapping);

  return max_space;  
}


void allocator_init(Allocator *workalloc, size_t size) {
  workalloc->d_size = size;
  workalloc->d_pos  = 0;
  workalloc->data = (char*)(malloc(size));
  memset(workalloc->data, 0, size); // touch pages
}

void allocator_free(Allocator *workalloc){
  free(workalloc->data);
  // free(workalloc);
}
    
char *allocator_allocate(Allocator *workalloc, size_t count, size_t size_elt)
{
  assert(workalloc->d_pos + count <= workalloc->d_size);
      
  char* ret = (workalloc->data + workalloc->d_pos);
      
  size_t to_alloc = count * size_elt;
      
  workalloc->d_pos += to_alloc;
  return ret;
}
        
char *allocator_current_ptr(Allocator *workalloc) {

  return (workalloc->data + workalloc->d_pos);
}

size_t allocator_pos(Allocator *workalloc) {
  return workalloc->d_pos;
}

void allocator_set_pos(Allocator *workalloc, size_t pos) {workalloc->d_pos = pos; }
    
  
