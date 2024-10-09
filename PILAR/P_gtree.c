/*
 *  P_GTREE.C
 *
 *	Generic module containing soubroutines to manipulate
 *	a "generic" tree type (i. e. a tree with an undefined
 *	number of subtrees at each node) for the project PILAR
 *	(Protein Identification and Library Access Resource)
 *
 *  Subroutines
 *
 *	g_tree_t gt_init(parent, data)	-- create a new tree
 *	g_tree_t parent;
 *	info_t *data;
 *
 *	g_tree_t gt_new(data)		-- begin a new generic tree
 *	info_t *data;
 *
 *	bool gt_add_child(t_node, data)	    -- add a new child to a node
 *	g_tree_t t_node;
 *	info_t *data;
 *
 *	bool gt_add_sibling(t_node, data)   -- add a brother to a node
 *	g_tree_t t_node;
 *	info_t *data;
 *
 *	l_node_t *gt_leftmost(t_node)	    -- find gt_leftmost child of a node
 *	g_tree_t t_node;
 *
 *	g_tree_t gt_root(t_node)	    -- find gt_root of a tree given
 *	g_tree_t t_node;		    -- any node in it.
 *
 *	bool gt_is_leaf(t_node)		    -- has branches this node?
 *	g_tree_t t_node;
 *
 *	void gt_pre_order(t_node, visit)    -- make a pre order traversal
 *	g_tree_t t_node;		    -- of a tree
 *	void (*visit)(g_tree_t);
 *
 *	void gt_in_order(t_node, visit)	    -- make an in-order traversal
 *	g_tree_t t_node;
 *	void (*visit)(g_tree_t);
 *
 *	void gt_post_order(t_node, visit)   -- make a post order traversal
 *	g_tree_t t_node;
 *	void (*visit)(g_tree_t);
 *
 *  Files:
 *	This module is dependent on the following files:
 *	    STDIO.H
 *	    PORTABLE.H
 *	    P_LIST.C
 *
 *  Caveats:
 *
 *  Designed by:
 *	José Ramón Valverde Carrillo
 *
 *  History:
 *	5 - jan - 1990	(Compiles OK)
 *     22 - jul - 1990	(Debugged)
 *	1 - aug - 1991	(Ported to Ultrix)
 *     18 - may - 1997	(Ported to Linux)
 *     24 - may - 1997	(Ported to FreeBSD)
 *
 *  Copyright:
 *	This is Public Domain software. You can freely use it
 *  as far as you do not use this code for commercial purposes.
 *
 *	J. R. Valverde	    jrvalverde@es.embnet.org	    	OK
 *  	    	    	    jrvalverde@ebi.ac.uk    	    	OUTDATED
 *			    JRamon@uamed.uam.es	(Internet)  	OUTDATED
 *  	    	    	    JRAMON@EMDUAM51	(BitNet, EARN)	OUTDATED
 *
 * $Log: P_gtree.c,v $
 * Revision 1.2  1997/05/24 17:09:58  jr
 * Ported to FreeBSD
 *
 * Revision 1.1  1997/05/18 20:33:56  jr
 * Initial revision
 *
 */

#include <stdio.h>

#include "portable.h"

#define MOD_DEBUG   0

/*
	Oui, je vois l'heure; il est l'Eternité !
			    Baudelaire. Le spleen de Paris.
*/

#define info_t	void
#define MAX	10

/* ----------- List data structures ------------ */

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
    list_t   		children;	/* header of the list of childs */
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
    In fame mea taberna,
    In nocte mea lucerna,
    Recte me semper guberna.
	Charles Baudelaire.
*/

/*------------------------------------------------------------------*/
/*			GENERIC SUBROUTINES			    */
/*------------------------------------------------------------------*/

/*  GT_INIT
	Create a new tree.

    INPUT:  DATA -- Data to store in order to begin the tree
	    PARENT -- Supposed parent of the new tree. Null if root.

    RETURN: A pointer to the new tree or subtree.

    J. R. Valverde
*/

g_tree_t gt_init(parent, data)
info_t *data;
g_tree_t parent;
{
    g_tree_t new;
    list_t new_list();

    /* check memory allocation and initialize */
    if (new = (g_tree_t) malloc( sizeof(struct dynamic_node) ) ) {
	new->element = data;
	/* check creation of a new list of childs */
	if ((new->children = new_list()) == NULL) {
	    free( new );
	    return NULL;
	}
	new->parent = parent;
    }

    return (new);
}

/*  GT_NEW
	begin a new generic tree.

    INPUT: DATA -- Information to begin the tree with.
    RETURNS: the address of the gt_root node.

    J. R. VALVERDE
*/

g_tree_t gt_new(data)
info_t *data;
{
    g_tree_t	gt_init();

    return (gt_init(NULL, data));
}

/*  GT_ADD_CHILD
	Add a child to any node in a tree
    Assumes that the argument t_node is not NULL

    INPUT:  DATA -- Information to be stored in the node
	    T_NODE -- Parent of the new node
    RETURN: TRUE if success.
*/

bool gt_add_child(t_node, data)
info_t *data;
g_tree_t t_node;
{
    g_tree_t new, gt_init();
    bool list_append();

    /* verify initialization */
    if (new = gt_init(t_node, data))
	/* link new node in the list of children */
	return (list_append(t_node->children, new));

    /* bad initialization */
    return FALSE;
}

/*  GT_ADD_SIBLING
	add a new brother an any tree node (but gt_root)
    assume that t_node->parent is not NULL

    INPUT:  DATA -- The information to be stored in the new brother
	    T_NODE -- The father (or mother) of the new node
    RETURN: TRUE if success.
        
    J. R. Valverde
*/

bool gt_add_sibling(t_node, data)
info_t *data;
g_tree_t t_node;
{
    bool gt_add_child();

    /* link new sibling in the list of childs of the parent node */
    return (gt_add_child(t_node->parent, data));
}
        
/*  GT_LEFTMOST
	Find gt_leftmost child of a node in a tree
    Assumes that t_node is not NULL

    It should be called as follows:

    l_node_t	*ptr, *gt_leftmost();
    g_tree_t	tree;
    ...
    ptr = gt_leftmost(tree);

    INPUT: T_NODE -- The subtree whose childs are to be explored
    RETURN: a pointer of type l_node_t to the first node in the list
  	of children.

    J. R. Valverde
*/

l_node_t *gt_leftmost(t_node)
g_tree_t t_node;
{
    return (t_node->children->l_first);
}

/*  GT_ROOT
	Find gt_root of a tree given any node
   Assumes that received node is not NULL

    INPUT:  T_NODE -- A pointer to any node in a tree.
    RETURN: A pointer to the gt_root node of the tree.

    J. R. Valverde
*/

g_tree_t gt_root(t_node)
g_tree_t t_node;
{
    /* look for a node with parent equal to null */
    while (t_node->parent)
	t_node = t_node->parent;

    /* return found gt_root */
    return (t_node);
}

/*  GT_IS_LEAF
	Has branches this node?
    Assumes that t_node is not NULL

    INPUT:  T_NODE -- The node to be explored
    RETURN: TRUE if it is a leaf node, i. e. it has no branches.

    J. R. Valverde
*/

bool gt_is_leaf(t_node)
g_tree_t t_node;
{
    /* has this node any child ? */ 
    if (t_node->children->l_length)	/* YES ==> it is not a leaf node */
	return (FALSE);
    else				/* NO ==> it is a leaf */
	return (TRUE);
}

/*  GT_PRE_ORDER
	Make a preorder traversal of the tree
    Assumes that function visit exists
    Assumes that t_node is not NULL

    INPUT:  T_NODE -- A pointer to the tree or subtree to be traversed
	    VISIT -- A pointer to the procedure to be applied to each
		visited node.

    J. R. Valverde
*/

void gt_pre_order(t_node, visit)
g_tree_t t_node;
void (*visit)( g_tree_t );
{
    l_node_t *temp, *gt_leftmost();

    if (t_node == NULL)
	return;
    /* first process gt_root node */
    (*visit)(t_node);

    /* now take every child (if any) from left to right */
    if (t_node->children->l_length) {
	temp = gt_leftmost(t_node);
	while (temp) {
	    gt_pre_order(temp->item, visit);
	    temp = temp->next;
	}
    }
}

/*  GT_IN_ORDER
	Make a in-order traversal of a tree
    Assumes that visit() exists
    Assumes that t_node is not NULL

    INPUT:  T_NODE -- A pointer to the tree or subtree to be traversed
	    VISIT -- A pointer to the procedure to apply when visiting
		a node

    J. R. Valverde
*/

void gt_in_order(t_node, visit)
g_tree_t t_node;
void (*visit)(g_tree_t);
{
    l_node_t *temp;

    if (t_node == NULL) return;

    if (t_node->children->l_length == 0)	/* no more childs */
	(*visit)(t_node);
    else {
	temp = gt_leftmost(t_node);
	gt_in_order(temp->item, visit);
	(*visit)(t_node);
	temp = temp->next;
	while (temp) {
	    gt_in_order(temp->item, visit);
	    temp = temp->next;
	}
    }
}

/*  GT_POST_ORDER
	Make a post-order traversal of a tree
    Assumes visit() exists
    Assumes t_node is not NULL

    INPUT:  T_NODE -- Pointer to the tree or subtree to be traversed
	    VISIT -- Pointer to the function to do when visiting a node

    J. R. Valverde
*/

void gt_post_order(t_node, visit)
g_tree_t t_node;
void (*visit)(g_tree_t);
{
    l_node_t *temp;

    if (t_node == NULL) return;

    /* first process children (if any) from left to right */
    if (t_node->children->l_length) {
	temp = gt_leftmost(t_node);
	while (temp) {
	    gt_post_order(temp->item, visit);
	    temp = temp->next;
	}
    }
    /* and then process gt_root node */
    (*visit)(t_node);
}

/*--------------------------------------------------------------------
		    DEBUGGING ROUTINES
--------------------------------------------------------------------*/
#if MOD_DEBUG == 1

int data[31] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
	      	11,12,13,14,15,16,17,18,19, 20,
	   	21,22,23,24,25,26,27,28,29, 30, 100};

void show(t)
g_tree_t t;
{
    int *pi;

    pi = t->element;
    printf("\t%d", *pi);
}

void show_leaf(t)
g_tree_t t;
{
    int *pi;

    pi = t->element;
    if (gt_is_leaf(t))
	printf("\t%d", *pi);
}

main()
{
    g_tree_t tree;
    l_node_t *node;
    int i;

    tree = gt_new(&data[0]);

    gt_add_child(tree, &data[1]);

    node = gt_leftmost(tree);
    tree = node->item;
    gt_add_sibling(tree, &data[11]);
    gt_add_sibling(tree, &data[21]);
    gt_add_child(tree, &data[2]);
    gt_add_child(tree, &data[12]);
    gt_add_child(tree, &data[22]);

    node = node->next;
    tree = node->item;
    gt_add_child(tree, &data[3]);
    gt_add_child(tree, &data[13]);
    gt_add_child(tree, &data[23]);


    node = node->next;
    tree = node->item;
    gt_add_child(tree, &data[4]);
    gt_add_child(tree, &data[14]);
    gt_add_child(tree, &data[24]);

    node = gt_leftmost(tree);
    tree = node->item;
    gt_add_child(tree, &data[5]);
    gt_add_child(tree, &data[15]);
    gt_add_child(tree, &data[25]);

    node = node->next;
    tree = node->item;
    gt_add_child(tree, &data[6]);
    gt_add_child(tree, &data[16]);
    gt_add_child(tree, &data[26]);

    node = node->next;
    tree = node->item;          
    gt_add_child(tree, &data[7]);
    gt_add_child(tree, &data[17]);
    gt_add_child(tree, &data[27]);

    tree = gt_root(tree);
    node = gt_leftmost(tree);
    node = gt_leftmost(node->item);
    tree = node->item;
    gt_add_child(tree, &data[8]);
    node = gt_leftmost(tree);
    gt_add_sibling(tree, &data[18]);
    gt_add_sibling(tree, &data[28]);

    tree = gt_root(tree);
    printf("\n--------\n");
    gt_pre_order(tree, show);
    printf("\n--------\n");
    gt_in_order(tree, show);
    printf("\n--------\n");
    gt_post_order(tree, show);

    printf("\n--------\n");
    gt_in_order(tree, show_leaf);
    printf("\n--------\n");
}

/*
    The tree produced is as follows:

				    1
				    |
		+-------------------+-------------------+
		|		    |			|
		2		    12			22
	+---+---+---+---+	+---+---+   +-----------+-----------+
	|   |	|   |	|	|   |	|   |		|	    |
	3   13	23  19	29	4   14	24  5		15	    25
    +---+				    |		|	    |
    |					+---+---+   +---+---+	+---+---+
    9					|   |	|   |	|   |	|   |	|
					6   16	26  7	17  27	8   18	28

*/

#endif
