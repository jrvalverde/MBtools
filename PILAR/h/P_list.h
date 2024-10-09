/*
 *  P_LIST.H
 *
 *	Include file containing the definitions for the module
 *  P_List.C
 *
 *	Requires PORTABLE.H and P_TYPES.H
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */

extern list_t new_list();

extern l_node_t *new_node(info_t *, l_node_t *);

extern boolean list_put(info_t *, list_t );

extern boolean list_insert(info_t *, l_node_t *, list_t );

extern boolean list_append(info_t *, list_t );

extern info_t *list_del_node(list_t );

extern info_t *list_del_next(l_node_t *, list_t);

extern int list_length(list_t);

extern info_t *list_first(list_t);

extern info_t *list_last(list_t);

extern info_t *list_next(list_t);

extern void list_reset(list_t);

extern bool list_find(list_t, info_t *, int (*)());

extern void list_traverse(list_t, void (*)());

extern info_t *list_del_first(list_t);

extern info_t *list_delete(list_t, info_t *, int (*)());

/*
 * See you
 */
