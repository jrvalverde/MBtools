/*
 *
 *
 *	MODULE
 *	    P_list.c -- subroutines for list manipulation
 *
 *	SYNOPSIS
 *	(I) Basic routines:
 *		new_list()
 *		new_node()
 *		list_put()
 *		list_insert()
 *		list_append()
 *		list_del_node()
 *	(II) Convenience routines:
 *		list_del_next()
 *		list_length()
 *		list_first()
 *  	    	list_current()
 *		list_last()
 *		list_next()
 *		list_reset()
 *		list_find()					
 *		list_traverse()					
 *		list_del_first()   				
 *		list_delete()					
 *								
 *	DESCRIPTION						
 *	    --	new_list allows the creation of a new list. A 	
 *	list record contains info about the list first and last	
 *	nodes, as well as the number of nodes in the list.	
 *	    --	new_node builds a new node, fills it with the	
 *	value to be stored and adds it to a specified list.	
 *	    --	list_put adds a node to the beginning of a	
 *	list.							
 *	    --	list_insert inserts a node in the middle of a	
 *	list AFTER a specified node. It is useful for building	
 *	topologically ordered lists.				
 *	    --	list_append adds a node to the end of a list.	
 *	    --	list_del_node erases a node from the beginning	
 *	of a list.						
 *								
 *	    --	list_del_next deletes a node anywhere.		
 *	    --	list_length returns the number of nodes.	
 *	    --	list_first returns the data stored in the first	
 *	node of the list.					
 *  	    --  list_current returns the data stored in the current
 *  	active node.
 *	    --	list_last returns the data stored in the last	
 *	node of the list.					
 *	    --	list_next returns the data stored in the next	
 *	node (the one following last consult).			
 *	    --	list_reset restores pointer to current node to	
 *	the first position in the list.				
 *	    --	list_find looks for a given item in a list and	
 *	returns TRUE if it can find it.				
 *	    --	list_traverse traverses a list processing all	
 *	its nodes with a given function.			
 *	    --	list_del_first removes first item in a list.	
 *	    --	list_delete removes any arbitrary item.		
 *								
 *	RETURNS							
 *								
 *	CAVEATS							
 *		The user should check the length of the list	
 *	before calling list_del_node.				
 *								
 *	FILES							
 *	    This module depends on the following files:		
 *		    STDIO.H					
 *		    PORTABLE.H					
 *								
 *	NOTES							
 *		You should use the second set of functions	
 *	and the type list_t whenever possible.			
 *		The first set of functions may be total or	
 *	partially hidden or modified in future releases.			
 *								
 *	SEE ALSO						
 *								
 *	DESIGNED BY						
 *	    José Ramón Valverde Carrillo.			
 *								
 *	HISTORY							
 *	    15 - jul - 1989 (J. R. Valverde)			
 *	    17 - jul - 1989 (J. R. Valverde) Debugged.		
 *	    21 - oct - 1989 (J. R. Valverde) Added second set	
 *		of functions and defined type list_t.		
 *	     1 - aug - 1991 (J. R. Valverde) Ported to Ultrix
 *  	    18 - may - 1997 (J. R. Valverde) Ported to Linux
 *  	    24 - may - 1997 (J. R. Valverde) Ported to FreeBSD
 *								
 *	COPYRIGHT						
 *		Up to now this is public domain. Please, do not	
 *	use it for commercial purposes without the author's	
 *	permission.						
 *
 * $Log: P_list.c,v $
 * Revision 1.2  1997/05/24 17:11:04  jr
 * Ported to FreeBSD
 *
 * Revision 1.1  1997/05/18 18:44:17  jr
 * Initial revision
 *
 */


#include <stdio.h>

#include "portable.h"

#define MOD_DEBUG   0
#define	info_t	void

struct l_node {
    info_t	    *item;
    struct l_node   *next;
};

typedef struct l_node l_node_t;

struct l_head {
    int		l_length;
    l_node_t	*l_first;
    l_node_t	*l_last;
    l_node_t	*l_cur;
};

typedef struct l_head l_head_t;

typedef l_head_t    *list_t;
/*
 * This allows for the external definition of type list_t as
 *
 *  typedef void *list_t;
 *
 * hidding in this way all internal details about lists from the
 * applications programmer.
 */

#define NEW(type)   (type *) malloc(sizeof(type))
                    
/*-------------------------------------------------------------------
	NEW_LIST
	    Creates a new list.

	OUTPUT: -- address of the list newly created, or NULL if
		   there was not enough memory available.

	J. R. VALVERDE
-------------------------------------------------------------------*/

list_t new_list()
{
    list_t new;

    /* allocate memory */
    if (new = NEW(l_head_t)) {
	/* Initialize all components */
	new->l_length = 0;
	new->l_first = new->l_last = new->l_cur = NULL;
    }
    /* return address of new list */
    return (new);
}

/*-------------------------------------------------------------------
	NEW_NODE
		Creates and initializes a new node, inserting it in
	a list.

	INPUT:	VAL -- value to put in the node
		PTR -- pointer to last node in the list

	J. R. VALVERDE
-------------------------------------------------------------------*/

static l_node_t *new_node(val, ptr)
info_t *val;
l_node_t *ptr;
{
    l_node_t *new;

    /* allocate memory and init record */
    if (new = NEW(l_node_t)) {
	new->item = val;
	new->next = ptr;
    }
    /* return the node */
    return (new);
}

/*-------------------------------------------------------------------
	LIST_PUT
	    Inserts a value at the beginning of a list

	INPUT:	DATA -- the value to be added
		LIST -- the list to which it must be added

	OUTPUT:	TRUE -- if success
		FALSE -- if failure

	J. R. VALVERDE
-------------------------------------------------------------------*/

boolean list_put(list, data)
list_t list;
info_t *data;
{
    l_node_t *new, *new_node();

    /* check there is enough memory */
    if ((new = NEW(l_node_t)) == NULL)
    	return (FALSE);
    else {
    	new->item = data;
    	new->next = list->l_first;
	/* Ok. Assign a pointer to the new header */
	list->l_first = new;	    /* link to the new node */
	if (list->l_length == 0)    /* if this is the first node */
            list->l_last = new;	    /* then assign pointer to the last node */
				    /* to this new one */
	list->l_length++;	    /* actualize list length */

   	return (TRUE);
    }
}

/*-------------------------------------------------------------------
	LIST_INSERT
	    Inserts a record inside a list after a specified
	node.

	INPUT:	DATA -- Data to fit in the new record
		NODE -- Node after which the record will be inserted
		LIST -- List to receive the new value.
	OUTPUT: TRUE if success, otherwise FAIL

	J. R. VALVERDE
-------------------------------------------------------------------*/

boolean list_insert(list, node, data)
list_t list;
l_node_t *node;
info_t *data;
{         
    l_node_t *new;
    boolean list_append();

    /* first check for the two special cases */
    if (node == NULL)		    /* It is the first */
	return ( list_put(list, data) );

    if (node == list->l_last)	    /* it is the last */
	return ( list_append(list, data) );

    /* insert new data between this node and the next one */
    if (new = new_node(data, node->next)) {
	node->next = new;	    /* update link and rebuild the list */
       	list->l_length++;	    /* update list length */
	return ( TRUE );
    }
    /* We were not able to introduce the record */
    return ( FALSE );
}

/*-------------------------------------------------------------------
	LIST_APPEND
		Append a value to the end of the list

	INPUT:	DATA -- The value to add
		LIST -- The list.

	J. R. VALVERDE
-------------------------------------------------------------------*/
boolean list_append(list, data)
list_t list;
info_t *data;
{
    l_node_t *new, *new_node();

    /* append data to the end of the list */
    if (new =  new_node(data, NULL)) {
	/* Initialization OK. Add data to the record */
	if (list->l_length)		/* if list is not empty */
            list->l_last->next = new;	/* link the new node */
	else				/* otherwise */
            list->l_first = new;	/* assign ptr to first to the new one */
	list->l_last = new;		/* actualize pointer to last */
	list->l_length++;		/* update list length */
	return (TRUE);
    }
    else
	/* initialization error */
	return (FALSE);
}
                
/*-------------------------------------------------------------------
	LIST_DEL_NODE
		Delete the first node in the list.

	INPUT:	LIST -- List whose first node is to be deleted
	OUTPUT:	-- pointer to the value stored in the deleted record

 	CAUTION: Check value of list->l_length before calling this function
	    (it is better if you want to avoid infinite loops)

	J. R. VALVERDE
-------------------------------------------------------------------*/

info_t *list_del_node(list)
list_t list;
{
    l_node_t *temp;
    info_t *data;

    if (list->l_length == 0)		/* if the list is empty */
	return NULL;
    /* Delete first node */
    temp = list->l_first;		/* save pointer to node */
    data = temp->item;          	/* save pointer to record */
    list->l_first = temp->next;		/* reassign pointer to first */
    list->l_length--;			/* adjust new length */
    if (list->l_length == 0)		/* there are no more records */
	list->l_last = NULL;		/* assign NULL to the pointer to last */
    free( (char *) temp);		/* free the memory used */

    /* return address of unlinked record */
    return (data);
}

/*
    Calla, calla el pastor, con tu dulce tañer
    has tornado a encender mis tristezas de amor
    y a quien llora un placer recordarle un dolor
    es hacerle tener una pena mayor.

		Miguel de Cervantes Saavedra.
*/

/*--------------------------------------------------------------------
	LIST_DEL_NEXT
	    Delete the node in the list that is after this one. (We
	need to know this to be able to set correctly the pointers)

	INPUT:	PREV -- Node previous to the one to be deleted
		LIST -- The list to be decreased.

	OUTPUT:	Address of the item stored in the deleted node.

	J. R. VALVERDE
--------------------------------------------------------------------*/

info_t *list_del_next(list, node)
list_t list;
l_node_t *node;
{
    l_node_t *victim;
    info_t *data;

    /* first item in the list */
    if (node == NULL)
	return (list_del_node(list));

    victim = node->next;

    /* have we gone past the end of the list ? */
    if (victim == NULL)
	return NULL;

    data = victim->item;	    /* save stored data */
    node->next = victim->next;	    /* maintain list continuity */
    list->l_length--;		    /* decrement list size */

    /* was this the last item ? */
    if (victim == list->l_last) {
	list->l_last = node;	    /* we have a new last item */
    }  
    free( victim );		/* kill the node */
    return ( data );
}

/*-------------------------------------------------------------------
	LIST_LENGTH
	    Return number of items stored in a list

	INPUT:	LIST -- the list
	OUTPUT:	Size of the list

	J. R. VALVERDE
-------------------------------------------------------------------*/

int list_length(list)
list_t list;
{
    return (list->l_length);
}

/*-------------------------------------------------------------------
	LIST_FIRST
	    Return first node in the list

	INPUT:	LIST -- the list
	OUTPUT:	Address of the first node in the list

	J. R. VALVERDE
--------------------------------------------------------------------*/

info_t *list_first(list)
list_t list;
{
    list->l_cur = list->l_first;
    return (list->l_first->item);
}

/*--------------------------------------------------------------------
	LIST_CURRENT
	    Return current node in the list

	INPUT:	LIST -- the list
	OUTPUT:	Address of current item in the list

	J. R. VALVERDE
--------------------------------------------------------------------*/

info_t *list_current(list)
list_t list;
{
    return (list->l_cur->item);
}
/*--------------------------------------------------------------------
	LIST_LAST
	    Return last node in the list

	INPUT:	LIST -- the list
	OUTPUT:	Address of last item in the list

	J. R. VALVERDE
--------------------------------------------------------------------*/

info_t *list_last(list)
list_t list;
{
    list->l_cur = list->l_last;
    return (list->l_last->item);
}

/*--------------------------------------------------------------------
	LIST_NEXT
	    Return next item in the list

	INPUT:	LIST -- the list
	OUTPUT:	Address of next item in the list

	J. R. VALVERDE
--------------------------------------------------------------------*/

info_t *list_next(list)
list_t list;
{
    /* if the list is empty or we are already at the end */
    if (list->l_cur == list->l_last)
	return NULL;

    /* else, move one step forward */
    list->l_cur = list->l_cur->next;
    return (list->l_cur->item);
}

/*--------------------------------------------------------------------
	LIST_RESET
	    Restore pointer to current node to the beginning of the
	list.
    	    UNNEEDED: to be removed.

	INPUT:	LIST -- the list

	J. R. VALVERDE
--------------------------------------------------------------------*/

void list_reset(list)
list_t list;
{
    list->l_cur = list->l_first;
}

/*--------------------------------------------------------------------
	LIST_FIND
	    Look for a given value in a list. Since you can store
	any kind of item in a list, you must also provide a comparison
	function to know the properties of each item. It sets the
	current node pointer to the node containing that item.

	INPUT:	DATA -- The value to be looked up.
		COMP -- A comparison function between items.
		    It must receive two pointer arguments and return
	      		> 0	if arg1 > arg2
			= 0	if arg1 = arg2
			< 0	if arg1 < arg2
		LIST -- The list to be searched.
	OUTPUT:	TRUE if the value searched was found.
		FALSE if not.

	J. R. VALVERDE
--------------------------------------------------------------------*/

bool list_find(list, data, comp)
list_t list;
info_t *data;
int (*comp)(info_t *, info_t *);
{
    l_node_t *tmp;

    tmp = list->l_first;
    while (tmp) {
	if ( (*comp)(data, tmp->item) == 0) {
	    list->l_cur = tmp;
	    return TRUE;
	}
	tmp = tmp->next;
    }
    return FALSE;
}

/*--------------------------------------------------------------------
	LIST_TRAVERSE
	    Traverse a list and process all the items in the list.
	You must provide a function to process the items.

	INPUT:	LIST -- The list to be traversed.
		PROCESS -- What we want to do with each item in the list
		    It must be a void function that receives a pointer
		    to an item as argument.

	J. R. VALVERDE
--------------------------------------------------------------------*/

void list_traverse(list, process)
list_t list;
void (*process)(info_t *);
{
    l_node_t *tmp;

    tmp = list->l_first;
    while (tmp) {
	(*process)(tmp->item);
	tmp = tmp->next;
    }
}

/*--------------------------------------------------------------------
	LIST_DEL_FIRST
	    Remove first item in the list and get its address

	INPUT:	LIST -- the list
	OUTPUT:	Address of the first item in the list

	J. R. VALVERDE 
--------------------------------------------------------------------*/ 

info_t *list_del_first(list)
list_t list;
{
    return (list_del_node(list));
}

/*--------------------------------------------------------------------
	LIST_DELETE
	    Remove an item in the list. You must provide a comparison
	function to know when the item has been found.

	INPUT:	LIST -- the list
		DATA -- The value to be deleted
		COMP -- A comparison function to know wether we have
	    found the item to delete. It must receive two pointers and
	    return 0 if the two items are the same.

	OUTPUT: Address of the item removed from the list. The reason
	    is that we don't know if the space occupied by it must be
	    freed or not, only the caller can know.

	J. R. VALVERDE
--------------------------------------------------------------------*/

info_t *list_delete(list, data, comp)
list_t list;
info_t *data;
int (*comp)(info_t*, info_t*);
{
    l_node_t *prev, *tmp;

    /* test if list is empty */
    if (list->l_length == 0)
	return NULL;

    /* Run over the list carrying a pointer to the previous item until
    we find the node to be deleted, and then, remove that node */

    prev = NULL;
    tmp = list->l_first;

    while (tmp) {
	if ( (*comp)(data, tmp->item) == 0) {
	    /* found, delete it. */
	    if (tmp == list->l_first)
	    	return(list_del_first(list));
    	    else
	    	return (list_del_next(list, prev));
	}
	prev = tmp;
	tmp = tmp->next;
    }

    /* the item has not been found */
    return NULL;
}

                    
/*-------------------------------------------------------------------
    DEBUGGING ROUTINES
-------------------------------------------------------------------*/

#if MOD_DEBUG != 0

#include <stdio.h>

int iary[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 100 };

int icmp(pi1, pi2)
int *pi1, *pi2;
{
    return ((*pi1 > *pi2)? 1 : ((*pi1 == *pi2) ? 0 : -1 ));
}

void putpint(pi)
int *pi;
{
    if (pi != NULL) 
    	printf("%d ", *pi);
    else
    	puts("Error: null value got during traversal");
}

main()
{
    list_t mylist;
    l_node_t *mynode;
    int i, *p;
    
    mylist = new_list();
    
    puts("Testing list_put(): Is this 4 3 2 1 0 ?");
    for (i = 0; i < 5; i++)
        list_put(mylist, &iary[i]);
    for (mynode = mylist->l_first; mynode != NULL; mynode = mynode->next)
    {
    	p = (int *) mynode->item;
    	printf("%d ", *p);
    }
    puts("\n");
    puts("Testing list_append(): is this 4 3 2 1 0 5 6 7 8 9 ?");
    for (i = 5; i < 10; i++)
        list_append(mylist, &iary[i]);
    for (mynode = mylist->l_first; mynode != NULL; mynode = mynode->next)
    {
    	p = (int *) mynode->item;
    	printf("%d ", *p);
    }
    puts("\n");
    puts("Testing list_insert(): is this 4 3 2 1 0 100 5 6 7 8 9 ?");
    mynode = mylist->l_first;
    for (i = 0; i < 4; i++) mynode = mynode->next;
    list_insert(mylist, mynode, &iary[10]);
    for (mynode = mylist->l_first; mynode != NULL; mynode = mynode->next)
    {
    	p = (int *) mynode->item;
    	printf("%d ", *p);
    }
    puts("\n");
    puts("Testing list_del_node(): is this 5 6 7 8 9 ?");
    for (i = 0; i < 6; i++)
    	list_del_node(mylist);
    for (mynode = mylist->l_first; mynode != NULL; mynode = mynode->next)
    {
    	p = (int *) mynode->item;
    	printf("%d ", *p);
    }
    puts("\n");
    puts("Testing list_del_next(): is this 5 6 8 9 ?");
    mynode = mylist->l_first;
    mynode = mynode->next;
    list_del_next(mylist, mynode);
    for (mynode = mylist->l_first; mynode != NULL; mynode = mynode->next)
    {
    	p = (int *) mynode->item;
    	printf("%d ", *p);
    }
    puts("\n");
    puts("Testing list_length(): is it 4?");
    printf("%d\n\n", list_length(mylist));
    puts("Testing list_first(): is it 5 ?");
    p = (int *) list_first(mylist);
    printf("%d\n\n", *p);
    puts("Testing list_last(): is it 9 ?");
    p = (int *) list_last(mylist);
    printf("%d\n\n", *p);
    puts("Testing list_next(): is this 5 6 8 9 ?");
    for (p = (int *) list_first(mylist);
    	p != NULL;
    	p = list_next(mylist))
    	printf("%d ", *p);
    puts("\n");
    puts("Testing list_traverse(): is this 5 6 8 9 ?");
    list_traverse(mylist, putpint);
    puts("\n");
    puts("Testing list_find(): is it OK ?");
    p = &iary[8];
    if (! list_find(mylist, p, icmp)) {
    	puts("Not OK");
    } else {
    	p = &iary[2];
    	if (list_find(mylist, p, icmp)) {
    	    puts("Not OK");
    	} else {
    	    puts("OK\n");
    	}
    }
    puts("Testing list_delete(): is this 5 6 9?");
    p = &iary[8];
    p = list_delete(mylist, p, icmp);
    list_traverse(mylist, putpint);
    if (p != NULL)
        printf(" (removed value %d)\n", *p);
}



#endif
