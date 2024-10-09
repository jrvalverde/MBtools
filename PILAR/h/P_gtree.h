/*
 *  P_GTREE.H
 *
 *	Include file with the definitions needed to access the
 *  generic tree subroutines module P_GTree.C
 *
 *	Requires PORTABLE.H, P_LIST.H and, hence, P_TYPES.H
 *
 *	Notes: NOW it contains the type definitions as well
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */

#ifndef info_t
#define info_t	void
#endif

#define MAX	10	    /* Maximum number of nodes in an N-ary tree */


/*----------------------------------------------------------------*/

/* node structure for a binary tree */

struct binary_node {
    info_t		*element;	/* node's data */
    struct binary_node	*left;		/* left child node */
    struct binary_node	*right;		/* right child node */
    struct binary_node	*parent;	/* parent node */
};

typedef struct binary_node *b_tree_t;

/*----------------------------------------------------------------*/

/* node structure for three children by node */
          
struct triple_node {
    info_t		*element;
    struct triple_node	*left;
    struct triple_node	*middle;
    struct triple_node	*right;
    struct triple_node	*parent;
};

typedef struct triple_node *t_tree_t;

/*-----------------------------------------------------------------*/

/* for more than three nodes it is better to use an array of nodes.
   This array may be built as an array of pointers to child nodes
    as follows:
*/

/* node structure for a fixed number of children nodes (MAX) */

struct n_ary_node {
    info_t	    *element;	    /* node's data */
    struct n_ary_node *child[MAX];  /* array of children nodes */
    struct n_ary_node *parent;	    /* parent node */
};

typedef struct n_ary_node *n_tree_t;

/*---------------------------------------------------------------*/

/* node structure for an arbitrary number of children */

struct dynamic_node {
    info_t		*element;	/* node's data */
    struct list_t    	*children;	/* header of the list of childs */
    struct dynamic_node	*parent;	/* parent node */
};

typedef struct dynamic_node *g_tree_t;	/* generic tree */

/*----------------------------------------------------------------*/

/* structure of a node for a linked list of trees (a forest) */

struct forest_node {
    g_tree_t		*fst_item;
    struct forest_node	*fst_next;
};

typedef struct forest_node forest_n_t;

/* header for a linked list of trees (forest) */

struct forest_head {
    int	  	fst_len;
    forest_n_t	*fst_first, *fst_last;
};

typedef struct forest_head *forest_t;

/*
 * Generic tree manipulation subroutines.
 */

extern g_tree_t init_gtree(info_t *, g_tree_t);

extern g_tree_t new_gtree(info_t *);

extern bool add_child(info_t *, g_tree_t);

extern bool add_sibling(info_t *, g_tree_t);

extern l_node_t *leftmost(g_tree_t);

extern g_tree_t root(g_tree_t);

extern is_leaf(g_tree_t);

extern void pre_order(g_tree_t, void (*)());

extern void in_order(g_tree_t, void (*)());

extern void post_order(g_tree_t, void (*)());

/*
 * And this is it!
 */

