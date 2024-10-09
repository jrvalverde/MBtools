/**
 *  \file   P_AVL.c
 *  \brief  Module to operate on AVL trees.
 *
 *	Subroutines to manipulate generic AVL-trees
 *  (Avelson, Velskii & Landis) for PILAR (Protein
 *  Identification and Library Access Resource).
 *
 *  Subroutines:
 *	avl_insert(p, x, compare, h)	-- Search & insertion in an AVL tree
 *	avl_tree_t *p;
 *	item_t x;
 *	int (*compare)();
 *	bool *h;
 *
 *	avl_delete(p, x, compare, h)	-- Delete a key in an AVL tree
 *	avl_tree_t *p;
 *	item_t x;
 *	int (*compare)();
 *	bool *h;
 *
 *  Notes:
 *	AVL trees are balanced trees: a tree is balanced if for
 *  each node the depths of its two subtrees differ in at most 1.
 *	Use only if information recovery is considerably
 *  more frequent than insertion.
 *
 *  Caveats:
 *
 *  Designed by:
 *	José Ramón Valverde Carrillo.
 *
 *  History:
 *	6-jan-1990    First implementation
 *	6-jul-1990    Frist debugging session
 *	1-aug-1991    Bug in AVL_delete removed
 *	1-aug-1991    Ported to ULTRIX
 *	sometime in 94	Ported to OSF/1 and fixed again bug an avl_delete 
 *		    that slipped in the latest freeze (tape backup used an 
 *		    old version).
 *	19-apr-1997   Ported to Linux and fixed again same bug that had
 *		    slipped through again in last freeze (CD-ROM master
 *		    used old wrong tape backup).
 *  	24-may-1997   Ported to FreeBSD
 *
 *  Copyright:
 *	This is Public Domain software. DO NOT use it for
 *  commercial purposes without prior written permission
 *  from the author.
 *	J. R. Valverde	    jr@es.embnet.org			OK
 *			    jose@ebi.ac.uk			OUTDATED
 *  	    	    	    JRValverde@biomed.iib.uam.es    	OUTDATED
 *			    JRamon@uamed.uam.es			OUTDATED
 *			    JRAMON@EMDUAM51	(BITNET, EARN)	OUTDATED
 *
 *  	    $Log: P_AVL.c,v $
 *  	    Revision 2.1  2003/06/06 09:17:01  jr
 *  	    Añadido nombre de versión y quitado código de depuración
 *
 *  	    Revision 1.2  2003/06/06 09:01:28  jr 
 * Activado el código de validación
 *
 * Revision 1.1  2003/06/06 08:52:39  jr
 * Initial revision
 *
 * Revision 1.2  1997/05/24 17:08:54  jr
 * Ported to FreeBSD
 *
 * Revision 1.1  1997/05/15 15:31:33  jr
 * Initial revision
 *
 */

char *version="$Revision";

#include <stdio.h>
#include <stdlib.h>

#include "portable.h"
/*#include "P_utils.h"		/ * For error routines */
/*#include "P_stree.h"		/ * For debug traversal routines */

#define MOD_DEBUG FALSE

extern void horror(char *routine, char *message);


typedef void *item_t;

struct avl_node {
    item_t	    key;
    struct avl_node *left;
    struct avl_node *right;
    flag	    equi;	/* -1, 0, +1 (left-, no-, right- skewed) */
};

typedef struct avl_node *avl_tree_t;

/*
 *  AVL_INSERT
 *	Search an AVL tree
 *
 *  J. R. Valverde
 */

void avl_insert(p, x, compare, h)
avl_tree_t *p;
item_t x;
int (*compare)();
bool *h;		/* FALSE at beginning of search */
{
    avl_tree_t p1, p2;
    int cmp;
                                                 
    if (*p == NULL) {
	/* the key is not in the tree, insert it */
	p1 = (avl_tree_t) malloc( sizeof(struct avl_node) );
	if (! p1) {
	    horror("AVL_INSERT", "Not enough memory");
	    return;
	}
	*h = TRUE;
	p1->key = x;
	p1->left = p1->right = NULL;
	p1->equi = 0;
	(*p) = p1;
    }
    else {
	/* look for the key in the tree */
   	cmp = (*compare)(x, (*p)->key);
	if (cmp < 0) {
	    avl_insert(&((*p)->left), x, compare, h);
	    if (*h) {
		/* inserted: left branch has grown */
		switch ((*p)->equi) {
		    case 1: /* the tree was unbalanced to the right,
			    now it is balanced */
			    (*p)->equi = 0;
			    *h = FALSE;
			    break;
		    case 0: /* the tree was balanced, now it is
			    unbalanced to the left */
			    (*p)->equi = -1;
			    break;
		    case -1:/* the tree was already unbalanced to
			    the left: restore equilibrium */
			    p1 = (*p)->left;
			    if (p1->equi == -1) {
				/* single left-left rotation */
				(*p)->left = p1->right;
				p1->right = (*p);
				(*p)->equi = 0;
		 		(*p) = p1;
			    }
			    else {
				/* double left-right rotation */
				p2 = p1->right;
				p1->right = p2->left;
				p2->left = p1;
				(*p)->left = p2->right;
				p2->right = (*p);
				if (p2->equi == -1)
				    (*p)->equi = 1;
				else
				    (*p)->equi = 0;
				if (p2->equi == 1)
				    p1->equi = -1;
				else
				    p1->equi = 0;
				(*p) = p2;
			    }
			    (*p)->equi = 0;
			    *h = FALSE;
			    break;
		}   /* end switch */
	    }	    /* end if (*h) */
	}	    /* end if (minor) */
	else if (cmp > 0) {
	    avl_insert(&((*p)->right), x, compare, h);
	    if (*h) {
		/* inserted: right branch has grown */
		switch ((*p)->equi) {
		    case -1: /* the tree was unbalanced to the
			    left, now it is balanced */
			    (*p)->equi = 0;
			    *h = FALSE;
			    break;
		    case  0: /* the tree was balanced, so it is
			    now unbalanced to the right */
			    (*p)->equi = 1;
			    break;
		    case 1: /* the tree was already unbalanced to
			    the right, restore equilibrium */
			    p1 = (*p)->right;
			    if (p1->equi == 1) {
				/* single right-right rotation */
				(*p)->right = p1->left;
				p1->left = (*p);
				(*p)->equi = 0;
				(*p) = p1;
			    }
			    else {
				/* double right-left rotation */
				p2 = p1->left;
				p1->left = p2->right;
				p2->right = p1;
				(*p)->right = p2->left;
				p2->left = (*p);
				if (p2->equi == 1)
				    (*p)->equi = -1;
				else
				    (*p)->equi = 0;
				if (p2->equi == -1)
				    p1->equi = 1;
				else
				    p1->equi = 0;
				(*p) = p2;
			    }
			    (*p)->equi = 0;
			    *h = FALSE;
			    break;
		}   /* end switch */
	    }	    /* end if (*h) */
	}	    /* end if (GREATER) */
	else
	    /* already in the tree */
 	    *h = FALSE;
    }		    /* end of look for the key in the tree */
}

/*
 *  AVL_DELETE
 *	Delete a key in an AVL tree
 *
 *  J. R. Valverde
 */

static void equileft(p, h)
avl_tree_t *p;
bool *h;
{
    avl_tree_t p1, p2;
    int	e1, e2;

    /* *h == TRUE ==> left branch has decreased */
    switch ((*p)->equi) {
	case -1: (*p)->equi = 0;	/* we're now balanced */
		break;
	case 0:	(*p)->equi = 1;		/* ie. we're skewed to the right */
		*h = FALSE;
		break;
	case 1: /* restore equilibrium */
		p1 = (*p)->right;
		e1 = p1->equi;
		if (e1 >= 0) {
		    /* single right-right rotation */
		    /*
		     * p->right = p->right-left
		     * p->right->left = p
		     * p = p->right
		     *
		     *      \                \
		     *       P                R
		     *      / \              / \
		     *     L   R    --->    P   RR
		     *        / \          / \
		     *       RL  RR       L   RL
		     *
		     * i.e. P is substitued by its right son
		     * and moved to its R-son's left branch which we know
		     * to be equal or shorter than the right one.
		     * The left branch of R-son is now the
		     * R-son of P
		     *
		     */
		    (*p)->right = p1->left;
		    p1->left = (*p);
		    if (e1 == 0) {
			(*p)->equi = 1;
			p1->equi = -1;
			*h = FALSE;
		    }
		    else {
			(*p)->equi = 0;
			p1->equi = 0;
		    }
		    (*p) = p1;
		}
		else {
		    /* double right-left rotation */
		    /*
		     * p->right->left = p->right->left->right
		     * p->right->left->right = p->right
		     * p->right = p->right->left->left
		     * p->right->left->left = p
		     * p = p->right-left
		     *
		     *       \                     \
		     *        P                    RL
		     *       / \                 /    \
		     *      L   R      --->     P      R
		     *         / \             / \    / \
		     *        RL  RR          L RLL  RLR RR
		     *       / \
		     *     RLL  RLR
		     *
		     * R is skewed to the left. So we move its left node
		     * up (in place of P) and share RL's childs among P
		     * and R
		     */
		    p2 = p1->left;
		    e2 = p2->equi;
		    p1->left = p2->right;
		    p2->right = p1;
		    (*p)->right = p2->left;
		    p2->left = (*p);
		    if (e2 == 1)
			(*p)->equi = -1;
		    else
			(*p)->equi = 0;
		    if (e2 == -1)
			p1->equi = 1;
		    else
			p1->equi = 0;
		    (*p) = p2;
		    p2->equi = 0;
		}
		break;
    }	    /* end switch */
}

static void equiright(p, h)
avl_tree_t *p;
bool *h;
{
    avl_tree_t p1, p2;
    int e1, e2;

    /* Converse of equileft */
    /* *h == TRUE ==> right branch has decreased */
    switch ( (*p)->equi ) {
	case 1: (*p)->equi = 0;
		break;
	case 0: (*p)->equi = -1;
		*h = FALSE;
		break;
	case -1: /* restore equilibrium */
		p1 = (*p)->left;
		e1 = p1->equi;
		if (e1 <= 0) {
		    /* single left-left rotation */
		    (*p)->left = p1->right;
		    p1->right = (*p);
		    if (e1 == 0) {
			(*p)->equi = -1;
			p1->equi = 1;
			*h = FALSE;
		    }
		    else {
			(*p)->equi = 0;
			p1->equi = 0;
		    }
		    (*p) = p1;
		}
		else {
		    /* double left-right rotation */
		    p2 = p1->right;
		    e2 = p2->equi;
		    p1->right = p2->left;
		    p2->left = p1;
		    (*p)->left = p2->right;
		    p2->right = (*p);
		    if (e2 == -1)
			(*p)->equi = 1;
		    else
			(*p)->equi = 0;
		    if (e2 == 1)
			p1->equi = -1;
		    else
			p1->equi = 0;
		    (*p) = p2;
		    p2->equi = 0;
		}
		break;
    }
}

/* *** Ugly, ugly, ugly... *** */
static avl_tree_t q = NULL;

static void del(r, h)
avl_tree_t *r;
bool *h;	/* *h == FALSE */
{
    extern avl_tree_t q;	/* has been assigned by caller avl_delete */

    if ((*r)->right != NULL) {
	del( &((*r)->right), h);
	if (*h)
	    equiright(r, h);
    }
    else {
	/* we move r's key to p	instead of moving the node r
	upwards recursively which is more costly; we sacrifice
	r instead of p (that's why we make q = r, so when at
	the end of recursion q is freed, it is r that is released */
	q->key = (*r)->key;
	(*r) = (*r)->left;
	q = (*r);
	*h = TRUE;		
    }
}

void avl_delete(p, x, compare, h)
avl_tree_t *p;
item_t x;
int (*compare)();
bool *h;
{
    /*
     * NOTE !!!
     *	The root of the tree might change (?) in this process.
     * P will keep pointing to it in any case.
     * We should return a pointer to the deleted key so it can
     * be freed by the caller. This is horriby inefficient since
     * the key must be passed through all the recursive calls.
     */
    extern avl_tree_t q;
    int cmp;

    if ((*p) == NULL) {
	/* the key is not in the tree */
	*h = FALSE;
	return;
    }
    else {
	cmp = (*compare)(x, (*p)->key);
	if (cmp < 0) {
	    avl_delete( &((*p)->left), x, compare, h);
	    if (*h) equileft(p, h);
	}
	else if (cmp > 0) {
	    avl_delete( &((*p)->right), x, compare, h);
	    if (*h) equiright(p, h);
	}
	else {
	    /* delete *p */
	    q = (*p);
	    /* The trivial cases */
	    if ((*p)->right == NULL) {
		(*p) = (*p)->left;
		*h = TRUE;
	    }
	    else if (q->left == NULL) {
		(*p) = (*p)->right;
		*h = TRUE;
	    }
	    else {
/*		del(&(q->left), h); */
		del(&((*p)->left), h);
		if (*h) equileft(p, h);
	    }
	    free(q);
	    q = NULL;
	}
    }
}

/*--------------------------------------------------------------------
		    DEBUGGING ROUTINES
--------------------------------------------------------------------*/


#if MOD_DEBUG

#define MINOR	    -1
#define EQUAL	    0
#define GREATER	    1

extern void st_preorder(avl_tree_t, void (*)());
extern void st_inorder(avl_tree_t, void (*)());
extern void st_postorder(avl_tree_t, void (*)());
extern void *st_search(avl_tree_t, void *, int (*)());
extern bool st_insert(avl_tree_t, void *, int (*)());
extern void *st_delete(avl_tree_t, void *, int(*)());
extern void *st_item(avl_tree_t);

int data[30] = { 0, 1, 2, 4, 3, 9, 7, 8, 6,10,
		20,15,17,13,18,16,14,12,19,11,
		32,40,50,27,68,41,58,46,5,100 };

int icmp( i1, i2 )
int *i1, *i2;
{
    return ((*i1 > *i2)? GREATER : (*i1 == *i2)? EQUAL : MINOR);
}

void show( i )
int *i;
{
    fprintf(stderr, "%4d", *i);
}

void horror(s1, s2)
char *s1, *s2;
{
    fprintf(stderr, "%s: %s\n", s1, s2);
}

void main()
{
    static int i;
    static int *pi;
    static avl_tree_t tree = NULL, chk;
    static bool is_in;

    is_in = FALSE;
    avl_insert(&tree, &data[0], icmp, &is_in);
    for (i = 1; i < 30; i++) {
	is_in = FALSE;
	avl_insert(&tree, &data[i], icmp, &is_in);
    }

    fprintf(stderr, "\nThis is preorder...\n");	
    st_preorder(tree, show);

    fprintf(stderr, "\nThis is inorder...\n");
    st_inorder(tree, show);

    fprintf(stderr, "\nAnd this postorder...\n");
    st_postorder(tree, show);

    i = 20; pi = &i;
    chk =  st_search(tree, pi, icmp);
    pi = chk->key;
    fprintf(stderr, "\n\nIs this 20?... %d\n", *pi);
    i = 5; pi = &i;
    chk = st_search(tree, pi, icmp);
    pi = chk->key;
    fprintf(stderr, "\nIs this 5?... %d\n", *pi);
    i = 1000; pi = &i;
    chk = st_search(tree, pi, icmp);
    if (chk == NULL)
	fprintf(stderr, "\n1000 not found... OK\n");
    else
	fprintf(stderr, "\nError, horror, where does 1000 come from?\n");

    fprintf(stderr, "\nLet us delete 20 values...\n");
    pi = &i;
    for (i = 1; i <= 20; ++i) {
	is_in = FALSE;
	avl_delete(&tree, pi, icmp, &is_in);
    }

    st_inorder(tree, show);
    fprintf(stderr, "\n Was this 0, 27, 32, 40, 41, 46, 50, 58, 68, 100 ?\n");
}

#endif


/*
 *  The end.
 */
