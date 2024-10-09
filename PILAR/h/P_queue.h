/*
 *  P_QUEUE.H
 *
 *	Include file with the definitions needed for module
 *  P_Queue.C
 *
 *	Requires PORTABLE.H and P_TYPES.H
 *
 *	Designed by:
 *	    J. R. Valverde	    8  - apr - 1990
 */

extern queue_t new_queue();

extern boolean qstore(info_t *, queue_t);      

extern info_t *qretrieve(queue_t);

extern fifo_t new_fifo();

extern boolean enqueue(info_t *, fifo_t);

extern info_t *dequeue(fifo_t);

