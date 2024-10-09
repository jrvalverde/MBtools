/**
 *  \file   P_AVL.h
 *  \brief  A class to manipulate AVL trees.
 *
 *	File to include to provide access to the routines for
 *  special manipulation of AVL trees in module P_AVL.C. The
 *  remaining functions must be provided by the module P_Stree.C
 *  (binary search tree subroutines).
 *
 *  Notes:
 *	AVL trees are balanced trees: a tree is balanced if for
 *  each node the depths of its two subtrees differ in at most 1.
 *	Use only if information recovery is considerably
 *  more frequent than insertion.
 *
 *  \see    	    P_AVL.c
 *  \see    	    P_stree.h
 *  \see    	    P_stree.c
 *
 *  \author	    J. R. Valverde
 *  \date   	    8-Apr-1990
 */

/** \typedef	void *item_t;
 *  \brief  	A type definition for a generic object store.
 */
typedef void *item_t;

/** \struct avl_node
 *  \brief  A node in an AVL tree.
 *
 *  The AVL tree node is defined to be compatible with generic
 *  binary tree nodes, so that generic binary search tree routines can
 *  be used on AVL trees.
 */
struct avl_node {
    item_t	    key;    /**< The object stored in this node. */
    struct avl_node *left;  /**< Left child of the node. */
    struct avl_node *right; /**< Right child of the node. */
    flag	    equi;   /**< \internal Whether this node's subtree is equilibrated or not. */
};

/** \typedef    struct avl_node *avl_tree;
 *  \brief  	A type definition for an AVL tree
 */
typedef struct avl_node *avl_tree;

/** \fn     void avl_insert(avl_tree *tree, item_t *item, int (*compare)(), bool *h);
 *  \brief  Search+Insert a value in an AVL tree.
 *
 *  If the item is already stored in the tree, sets h to TRUE.
 *  If the item is not in the tree, store it.
 *  To actually search without inserting the key, use the generic
 *  binary tree functions.
 *
 *  \param  tree    The AVL tree to operate upon.
 *  \param  item    The item to search/store.
 *  \param  compare A function to compare two items in the tree.
 *  \param  h	    A flag set to FALSE to find out if the key was in the tree.
 *
 *  \warning
 *	The root of the tree might change (?) in this process.
 * P will keep pointing to it in any case.
 * We should return a pointer to the deleted key so it can
 * be freed by the caller. This is horriby inefficient since
 * the key must be passed through all the recursive calls.
 */
extern void avl_insert(avl_tree *tree, item_t *item, int (*compare)(), bool *h);

/** \fn     avl_delete(avl_tree *tree, item_t *item, int(*compare)(), bool *h);
 *  \brief  Delete a key from an AVL tree.
 *
 *  If the requested key was not found in the tree, h will be set to FALSE.
 *
 *  \param  tree    The AVL tree to operate upon.
 *  \param  item    The item to remove.
 *  \param  compare A function to compare two items in the tree.
 *  \param  h	    A flag set to TRUE to find out if the key was in the tree.
 */
extern void avl_delete(avl_tree *tree, item_t *item, int(*compare)(), bool *h);

/*
 *
 */
