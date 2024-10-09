/**
 *  \file P_dlist.h
 *  \brief Definitions for the module P_dlist.c
 *
 *  \pre portable.h
 *  \pre P_types.h
 *
 *  \see portable.h
 *  \see p_types.h
 *
 *  \author   J. R. Valverde
 *  \date   8-apr-1990
 */

/** dlist_t new_dlist()
 *  Create a doubly linked list.
 *
 *  \return A pointer to the newly creted doubly linked list.
 *
 *  \author J. R. VALVERDE
 */
extern dlist_t new_dlist();

/** dl_node_t *new_dnode(info_t *, dl_node_t *, dl_node_t *)
 *  Initialize a node for a doubly linked list.
 *
 *  \param DATA -- information to store in this node.
 *  \param  PREV -- previous node in the list.
 *  \param  NEXT -- next node in the list.
 *  \return A pointer to the newly created node.
 *
 *  \author J. R. VALVERDE
 *
 */
extern dl_node_t *new_dnode(info_t *, dl_node_t *, dl_node_t *);

/** boolean dlist_put(info_t *, dlist_t )
    Insert an item at the beginning of a doubly linked list.

    \param  DATA -- The information to add to the front of the list.
    \param  LIST -- The list that will receive the information.
    \return TRUE if it is able to add the new item, FALSE if not.

    \author J. R. VALVERDE
 */
extern boolean dlist_put(info_t *, dlist_t );

/** boolean dlist_append(info_t *, dlist_t )
	Append data to the end of a doubly linked list.

    \param  DATA -- Information to append to the list
    \param  LIST -- List that will receive the item
    \return TRUE if all goes OK, FALSE if not.

    \author J. R. VALVERDE
 */
extern boolean dlist_append(info_t *, dlist_t );

/** info_t *dlist_del_first(dlist_t)
	Erase the first node from a doubly linked list.

    \warning	The value of list->dl_length should be verified before
	calling this routine.
    \param  LIST -- The list that will be deprived of a node
    \return A pointer to the information stored in the deleted node

    \author J. R. VALVERDE
 */
extern info_t *dlist_del_first(dlist_t);

/** bool dlist_insert(info_t *, dl_node_t *, dlist_t)
	    Insert a node after a given one in a doubly linked list

	\param	DATA -- The data to store.
	\param	NODE -- The node after which we are to insert the new one.
	\param	LIST -- The list
	\return	TRUE if the node could succesfully be inserted, 
	    	FALSE if an error occurred.

	\author J. R. VALVERDE
 */
extern bool dlist_insert(info_t *, dl_node_t *, dlist_t);

/** info_t *dlist_del_last(dlist_t)
	    Delete last node in a doubly linked list

	\param	LIST -- The list
	\return	Address of the item stored in the last node of the list

	\author J. R. VALVERDE
 */
extern info_t *dlist_del_last(dlist_t);

/** info_t *dlist_del_next(info_t *, dlist_t)
	    Delete next node in the list.

	\param	NODE -- The node to be deleted.
	\param LIST -- The list to be decreased.

	\return	Address of the item stored in the deleted node.

	\author J. R. VALVERDE
 */
extern info_t *dlist_del_next(info_t *, dlist_t);

/** info_t *dlist_del_node(dl_node_t *, dlist_t)
	    Delete this node in a doubly linked list.

	\param	VICTIM -- The node to be deleted
	\param	LIST -- The list.
	\return The address of the item stored in that node

	\author J. R. VALVERDE
 */
extern info_t *dlist_del_node(dl_node_t *, dlist_t);

/** int dlist_length( dlist_t )
	    Return number of items stored in a list

	\param	LIST -- the list
	\return	Size of the list

	\author J. R. VALVERDE
 */
extern int dlist_length( dlist_t );

/** info_t *dlist_first( dlist_t )
	    Return first node in the list

	\param	LIST -- the list
	\return	Address of the first node in the list

	\author J. R. VALVERDE
 */
extern info_t *dlist_first( dlist_t );

/** info_t *dlist_last( dlist_t )
	    Return last node in the list

	\param	LIST -- the list
	\return	Address of last item in the list

	\author J. R. VALVERDE
 */
extern info_t *dlist_last( dlist_t list);

/** info_t *dlist_current(list)
	    Return current node in the list

	\param	LIST -- the list
	\return	Address of current item in the list

	\author J. R. VALVERDE
 */
extern info_t *dlist_current(dlist_t list);

/** info_t *dlist_next( dlist_t )
	    Return next item in the list

	\param	LIST -- the list
	\return	Address of next item in the list

	\author J. R. VALVERDE
 */
extern info_t *dlist_next( dlist_t );

/** info_t *dlist_prev( dlist_t )
	    Return previous item in the list

	\param	LIST -- The list
	\return Address of the item.

	\author J. R. VALVERDE
 */
extern info_t *dlist_prev( dlist_t );

/** info_t *dlist_begset( dlist_t )
	    Restore pointer to current node to the beginning of the list.

	\param	LIST -- the list

	\author J. R. VALVERDE
 */
extern info_t *dlist_begset( dlist_t );

/** info_t *dlist_endset( dlist_t )
	    Restore pointer to current node to the ending of the list.

	\param	LIST -- the list

	\author J. R. VALVERDE
 */
extern info_t *dlist_endset( dlist_t );

/** bool dlist_ffind(dlist_t, info_t *, int (*)())
	    Look for a given value in a list. 
	
	Since you can store
	any kind of item in a list, you must also provide a comparison
	function to know the properties of each item. It sets the
	current node pointer to the node containing that item.
	Lookes for in the forward direction.

	\param	DATA -- The value to be looked up.
	\param	COMP -- A comparison function between items,
		    it must receive two pointer arguments and return
			\li > 0	if arg1 > arg2
			\li = 0	if arg1 = arg2
			\li < 0	if arg1 < arg2
	\param	LIST -- The list to be searched.
	\return	TRUE if the value searched was found,
		FALSE if not.

	\author J. R. VALVERDE
 */
extern bool dlist_ffind(dlist_t, info_t *, int (*)());

/** bool dlist_bkfind(dlist_t, info_t *, int (*)())
	    Look for a given value in a list. 
	
	Since you can store
	any kind of item in a list, you must also provide a comparison
	function to know the properties of each item. It sets the
	current node pointer to the node containing that item.
	Lookes for in the backward direction.

	\param	DATA -- The value to be looked up.
	\param	COMP -- A comparison function between items;
		    it must receive two pointer arguments and return
			\li > 0	if arg1 > arg2
			\li = 0	if arg1 = arg2
			\li < 0	if arg1 < arg2
	\param	LIST -- The list to be searched.
	\return	TRUE if the value searched was found,
		FALSE if not.

	\author J. R. VALVERDE
 */
extern bool dlist_bkfind(dlist_t, info_t *, int (*)());

/** void dlist_ftraverse(dlist_t, void (*)())
    Traverse a list forwards.
    
	    Traverse a list and process all the items in the list going
	in forward direction.
	You must provide a function to process the items.

	\param	LIST -- The list to be traversed.
	\param	PROCESS -- What we want to do with each item in the list;
		    it must be a void function that receives a pointer
		    to an item as argument.

	\author J. R. VALVERDE
 */
extern void dlist_ftraverse(dlist_t, void (*)());

/** void dlist_bktraverse(dlist_t, void (*)())
    Traverse a list backwards.
    
	    Traverse a list and process all the items in the list going
	in backward direction.
	You must provide a function to process the items.

	\param	LIST -- The list to be traversed.
	\param	PROCESS -- What we want to do with each item in the list
		    It must be a void function that receives a pointer
		    to an item as argument.

	\author J. R. VALVERDE
 */
extern void dlist_bktraverse(dlist_t, void (*)());

/** info_t *dlist_fdelete(dlist_t, info_t *, int (*)())
	    Remove an item in the list. 
	    
	    You must provide a comparison
	function to know when the item has been found. This function
	scans the list in the forward direction.

	\param	LIST -- the list
	\param	DATA -- The value to be deleted
	\param	COMP -- A comparison function to know wether we have
	    found the item to delete. It must receive two pointers and
	    return 0 if the two items are the same.

	\return Address of the item removed from the list. The reason
	    is that we don't know if the space occupied by it must be
	    freed or not, only the caller can know.

	\author J. R. VALVERDE
 */
extern info_t *dlist_fdelete(dlist_t, info_t *, int (*)());

/** info_t *dlist_bkdelete(dlist_t, info_t *, int (*)())
	    Remove an item in the list. 
	    
	    You must provide a comparison
	function to know when the item has been found. This function
	scans the list in the backward direction.

	\param	LIST -- the list
	\param	DATA -- The value to be deleted
	\param	COMP -- A comparison function to know wether we have
	    found the item to delete. It must receive two pointers and
	    return 0 if the two items are the same.

	\return Address of the item removed from the list. The reason
	    is that we don't know if the space occupied by it must be
	    freed or not, only the caller can know.

	\author J. R. VALVERDE
 */
extern info_t *dlist_bkdelete(dlist_t, info_t *, int (*)());

/*
 *  Nothing compares 2 u!
 */

