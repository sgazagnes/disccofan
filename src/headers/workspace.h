#ifndef WORKSPACE_H
#define WORKSPACE_H


size_t workspace_create_mappings(size_t len);
size_t workspace_construct(ulong max_rank);
size_t workspace_size(ulong size_base);
void allocator_init(Allocator *workalloc, size_t size);
void allocator_free(Allocator *workalloc);
char *allocator_allocate(Allocator *workalloc, size_t count, size_t size_elt);
char *allocator_current_ptr(Allocator *workalloc);
size_t pos(Allocator *workalloc);
void set_pos(Allocator *workalloc, size_t pos);
    
#endif
