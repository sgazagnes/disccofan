
#include "types.h"
#include "queue.h"
#include "flood.h"
#include "workspace.h"

/* ++++++++++++++++++++++++++ */
/*                            */
/*     Priority Q Salembier   */
/*                            */
/* ++++++++++++++++++++++++++ */

Queue *create_queue(ulong imgsize, ulong maxvalue) {
  Queue *hq;
  hq 		= calloc(maxvalue, sizeof(Queue));      check_alloc(hq, 400);  
  hq->pixels 	= calloc(imgsize, sizeof(ulong));       check_alloc(hq->pixels, 401);  
  return hq;
}

void set_queue_offsets(Queue *hq, ulong *numpixelsperlevel, ulong maxvalue) {
  ulong i;
  hq->head = hq->tail = 0;
  for (i = 1; i < maxvalue; i++) {
    hq[i].pixels = hq[i - 1].pixels + numpixelsperlevel[i - 1];
    hq[i].head = hq[i].tail = 0;
  }
}

void free_queue(Queue *hq) {
  free(hq->pixels);
  free(hq);
}

/* ++++++++++++++++++++++++++ */
/*                            */
/*     Priority Q Wilkinson   */
/*                            */
/* ++++++++++++++++++++++++++ */

pStack *pStackCreate(ulong maxsize, ulong *stackArray){
  pStack *newStack    = malloc(sizeof(pStack));    check_alloc(newStack, 402);
  newStack->array     = stackArray;
  newStack->size_max  = maxsize;
  newStack->pos_cur   = maxsize;
  return newStack;
}

void pStackDelete(pStack *oldstack){
  free(oldstack);
}


pQueue *pQueueCreate(ulong maxsize){
  pQueue *newQueue 	= (pQueue *) malloc(sizeof(pQueue));  check_alloc(newQueue, 403);
  newQueue->array 	= malloc((maxsize+1)*sizeof(ulong));  check_alloc(newQueue->array, 404);
  newQueue->size_cur 	= 0;
  newQueue->size_max 	= maxsize;
  return newQueue;
}


void pQueueDelete(pQueue *oldqueue){
  free(oldqueue->array);
  free(oldqueue);
}

ulong pQueuePop(value* gvals, pQueue *queue){
  ulong val_out = queue->array[1];
  ulong current = 1, moved;
  value val_cur;

  moved = queue->array[queue->size_cur];
  queue->size_cur--;
  val_cur = gvals[moved];

  while ( ((current*2<=queue->size_cur) &&
	   (val_cur< gvals[queue->array[current*2]]))
	  ||
	  ((current*2+1<=queue->size_cur) &&
	   (val_cur< gvals[queue->array[current*2+1]]))
	  ){
    if ((current*2+1<=queue->size_cur) && 
	(gvals[queue->array[current*2]]< 
	 gvals[queue->array[current*2+1]])){
      queue->array[current]= queue->array[current*2+1];
      current+=current+1;
    } else {
      queue->array[current]= queue->array[current*2];
      current+=current;
    }
  }
  queue->array[current]=moved;
  return val_out;
}

void pQueuePush(value *gvals, pQueue *queue, ulong pixpos){
  long current;
  value val_cur = gvals[pixpos];
  queue->size_cur++;
  current=queue->size_cur;
  
  while ((current/2 !=0) && (gvals[queue->array[current/2]]<val_cur)){
    queue->array[current]= queue->array[current/2];
    current=current/2;
  }

  queue->array[current]=pixpos;
}

/* ++++++++++++++++++++++++++ */
/*                            */
/*     Priority Q Teeninga    */
/*                            */
/* ++++++++++++++++++++++++++ */


ulong num_levels(ulong size){
  int bits_per_rank = bit_scan_reverse(size) + 1;
  // number of required levels
  return (bits_per_rank + bits_per_word_log2() - 1) / bits_per_word_log2();
}

PrioQueue *create_prio_queue(ulong max_rank, Allocator *work_alloc) {
  PrioQueue *queue = malloc(sizeof(PrioQueue));
  queue->m_num_levels = num_levels(max_rank);
  queue->m_levels = (ulong **) allocator_allocate(work_alloc, queue->m_num_levels, sizeof(ulong*));
  for (int i = queue->m_num_levels; i--;){
    max_rank >>= bits_per_word_log2();
    queue->m_levels[i] = (ulong *) allocator_allocate(work_alloc, max_rank+1, sizeof(ulong));
    memset(queue->m_levels[i], 0,  (max_rank+1)*sizeof(ulong));
  }
  queue->m_top = 0;
  return queue;
}


void insert_prio_queue(PrioQueue *q, ulong prefix) {

  if (prefix > q->m_top) {
    q->m_top = prefix;

  }
   for (int i = q->m_num_levels; i--;) {     
     ulong word_idx = prefix >> bits_per_word_log2();
     int bit_idx = prefix & (bits_per_word() - 1);
     q->m_levels[i][word_idx] |= (ulong)(1) << bit_idx;
     prefix = word_idx;
   }
}                              


void prio_queue_remove(PrioQueue *q) {
  ulong prefix = q->m_top;
  for (int i = q->m_num_levels; i--;){
    ulong word_idx = prefix >> bits_per_word_log2();
    q->m_levels[i][word_idx] &= ~((ulong)(1) << (prefix & (bits_per_word() - 1)));
    if (q->m_levels[i][word_idx] != 0) {
      for (int j = i; j != q->m_num_levels; ++j)	
	  word_idx = (word_idx << bits_per_word_log2()) ^ bit_scan_reverse(q->m_levels[j][word_idx]);  
	
      q->m_top = word_idx;

      return;
    }
    prefix = word_idx;
  }
}

void free_prio_queue(PrioQueue *q) {
  for(int i = 0; i != q->m_num_levels; i++)
    free(q->m_levels[i]);
  free(q->m_levels);
  free(q);
}
