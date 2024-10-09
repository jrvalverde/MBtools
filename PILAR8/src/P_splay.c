struct splay_node {
	void *key;
	struct splay_node *left;
	struct splay_node *right;
	struct splay_node *parent;
};

typedef struct splay_node splay_node_t;
typedef struct splay_node *splay_tree;

static void zig_left(tree)
splay_node_t *tree;
{
	splay_node_t *parent,*right;
	
	parent = tree->parent;
	right = tree->right;
	tree->parent = NULL;
	if (right != NULL)
		right->parent = parent;
	parent->left = right;
	parent->parent = tree;
}

static void zig_right(tree)
splay_node_t *tree;
{
	splay_node_t *parent, *left;
	
	parent = tree->parent;
	left = tree->left;
	tree->parent = NULL;
	if (left != NULL)
		left->parent = parent;
	parent->right = left;
	parent->parent = tree;
}

static void zig_zig_left(t)
splay_node_t *t;
{
	splay_node_t *p, 		/* parent */
				 *pr, 		/* parent's right */
				 *g, 		/* grand parent */
				 *r, 		/* right */
				 *gp;		/* grand-grand parent */
	
	p = t->parent;
	g = p->parent;
	gp = g->parent;
	r = t->right;
	pr = p->right;
	p->parent = t;
	g->left = pr;
	g->parent = p;
	t->parent = gp;
	t->right = p;
	p->left = r;
	if (r != NULL)
		r->parent = p;
	if (pr != NULL)
		pr->parent = g;
	if (gp != NULL) {
		if (gp->left == g)
			gp->left = t;
		else
			gp->right = t;
	}
}

static void zig_zig_right(t)
splay_node_t *t;
{
	splay_node_t *p, *pl, *g, *l, *gp;
	
	p = t->parent;
	g = p->parent;
	gp = g->parent;
	l = t->left;
	pl = p->left;
	p->left = g;
	p->parent = t;
	g->right = pl;
	g->parent = p;
	t->parent = gp;
	t->left = p;
	p->right = l;
	if (l != NULL)
		l->parent = p;
	if pl != NULL)
		pl->parent = g;
	if (gp != NULL) {
		if (gp->left == g)
			gp->left = t;
		else
			gp->right = t;
	}
}

static void zig_zag_left(t)
splay_node_t *t;
{
	splay_node_t *p, *gp, *ggp, *l, *r;
	
	p = t->parent;
	gp = p->parent;
	ggp = gp->parent;
	l = t->left;
	r = t->right;
	gp->parent = t;
	p->parent = t;
	p->right = l;
	gp->left = r;
	if (l != NULL)
		l->parent = p;
	if (r != NULL)
		r->parent = gp;
	if (ggp != NULL) {
		if (ggp->left == gp)
			ggp->left = t;
		else
			ggp->right = t;
	}
}

static void zig_zag_right(t)
splay_node_t *t;
{
	splay_node_t *p, *gp, *ggp, *l, *r;
	
	p = t->parent;
	gp = p->parent;
	ggp = gp->parent;
	l = t->left;
	r = t->right;
	t->left = gp;
	t->right = p;
	t->parent = ggp;
	p->parent = t;
	gp->parent = t;
	gp->right = l;
	p->left = r;
	if (l != NULL)
		l->parent = gp;
	if (r != NULL)
		r->parent = p;
	if (ggp != NULL) {
		if (ggp->left == gp) 
			ggp->left = t;
		else
			ggp->right = t;
	}
}

/* Splay functions
 * Assumptions: Tree is pointed to by global variable T, type TREENODE *
 *	External function Compare_key(KEYTYPE k1, KEYTYPE k2) returns
 *		0 if the two are equal, -1 if first key is less than
 *		second, 1 if first is greater than second
 *	External rotation functions
 */
 
static void Splay(tree, k, compare)
splay_node_t *tree;
info_t *k;
int (*compare)(info_t *, info_t *);
{
	int gle;
	/* Assume global tree T. Traverse T looking for key. */
	
	while (gle = (*compare)(T->key, k) != 0) {
		if (gle > 0) {
			if (T->left)
				T = T->left;
			else
				break;
		}
		else {
			if (T->right)
				T = T->right;
			else
				break;
		}
	}
	
	/* T now points to the node containing the key k, or to the inorder
	 * sucessor or predecessor to k. We don't really care which at this
	 * point. T will be the root when its parent pointer is null.
	 */
	while (T->parent) {
		if (T->parent->parent) {
			/* zig-zig or zig-zag */
			if (T->parent->parent->left == T->parent) {
				/* parent is a left child */
				if (T->parent->left == T) {
					/* left child of left child */
					zig_zig_left(T);
				}
				else {
					/* right child of left child */
					zig_zag_left(T);
				}
			}
			else {
				/* parent is a right child */
				if (T->parent->right == T) {
					/* right child of right child */
					zig_zig_right(T);
				}
				else {
					/* left child of right child */
					zig_zag_right(T);
				}
			}
		}
		else {
			/* zig */
			if (T->parent->left == T) {
				zig_left(T);
			}
			else {
				zig_right(T);
			}
		}
	}
}

/*
 * Splay tree standard operations -- Find, Delete and Insert.
 * Assumptions
 *	Tree is pointed to by global variable T, type TREENODE *
 *	External function Destroy_key recovers memory used by a key.
 */

info_t *Find(k, compare)
info_t *k;
int (*compare)(info_t *, info_t *);
{
	Splay(k, compare);
	/* if k was in the tree, it's now the root node */
	if ((*compare)(T->key, k) == 0)
		return T->key;
	else
		return NULL;
}

void Delete(k, compare)
info_t *k;
int (*compare)(info_t *, info_t *);
{
	splay_node_t *l, *r;
	
	/* bring target node to the root */
	Splay(k, compare);
	/* make sure the key was in the tree... */
	if ((*compare)(T->key, k) == 0) {
		/* Detach the node and dispose of it */
		l = T->left;
		r = T->right;
		Destroy_node(T);
		l->parent = NULL;
		r->parent = NULL;
		/* Splay the left subtree to find new root of the tree:
		   the key won't be found and we get the inorder sucessor
		   or predecessor of the node which becomes the new root. */
		T = l;
		Splay(k, compare);
		/* Root and left subtree are fine now, just attach the right
		subtree */
		T->right = r;
	}
}

void Insert(k, compare)
info_t *k;
int (*compare)(info_t *, info_t *);
{
	splay_node_t *x;
	
	Splay(k, compare);	/* find its place */
	/* insert the new node into the tree */
	Attach_Node(k);		/* could be all done with st_insert() alone??? */
	
	/* Now splay to bring the new node to the root. This is
	somewhat inefficient because Splay searches the tree all over again */
	Splay(k, compare);
}
