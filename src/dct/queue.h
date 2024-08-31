#ifndef QUEUE_H
#define QUEUE_H

#define queue_first(hq,h)    		(hq[h].pixels[hq[h].head++]) 
#define queue_add(hq,h,p)    		(hq[h].pixels[hq[h].tail++] = p)
#define queue_is_not_empty(hq,h) 	(hq[h].head != hq[h].tail)
#define is_empty_pstack(stack)        ((stack->pos_cur)==stack->size_max)
#define is_full_pstack(stack)         ( stack->pos_cur == 0 )
#define is_empty_pqueue(queue)        ((queue->size_cur)==0)
#define is_full_pqueue(queue)         ( queue->size_cur == queue->size_max )
#define pstack_top(stack)      	( stack->array[stack->pos_cur + 1])
#define pstack_pop(stack)       	(stack->array[++stack->pos_cur])
#define pstack_push(stack,elem)      	(stack->array[stack->pos_cur--]=elem)
#define pqueue_front(queue)       	(queue->array[1])

Queue *create_queue(ulong imgsize, ulong maxvalue);
void set_queue_offsets(Queue *hq, ulong *numpixelsperlevel, ulong maxvalue); 
void free_queue(Queue *hq);

pStack *pStackCreate(ulong maxsize, ulong *stackArray);
void pStackDelete(pStack *oldstack);
pQueue *pQueueCreate(ulong maxsize);
void pQueueDelete(pQueue *oldqueue);
ulong pQueuePop(value *gvals, pQueue *queue);
void pQueuePush(value *gvals, pQueue *queue, ulong pixpos);
PrioQueue *create_prio_queue(ulong max_rank, Allocator *work_alloc);
void insert_prio_queue(PrioQueue *q, ulong prefix);
void prio_queue_remove(PrioQueue *q);
void free_prio_queue(PrioQueue *q);
ulong num_levels(ulong size);
#endif
