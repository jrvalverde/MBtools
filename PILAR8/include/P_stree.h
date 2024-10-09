/*
 *  P_STREE.H
 *
 *	Include file to provide access to the routines for
 *  manipulation of a generic binary Search Tree, from module
 *  P_Stree.C
 *	These routines assume a given minimum tree structure,
 *  but you should be able to use them also with enlarged
 *  structures as far as the original/initial information,
 *  remains unchanged.
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */



struct stnode_t {
    info_t	    *item;
    struct stnode_t *left;
    struct stnode_t *right;
    };

typedef struct stnode_t *stree_t;	/* search tree type */

extern void st_preorder(stree_t, void (*)());

extern void st_inorder(stree_t, void (*)());

extern void st_postorder(stree_t, void (*)());

extern stree_t st_search(stree_t, info_t *, int (*)());

extern bool st_insert(stree_t, info_t *, int (*)());

extern stree_t st_delete(stree_t, info_t *, int (*)());

/*
 *  End of P_STREE.H
 */

