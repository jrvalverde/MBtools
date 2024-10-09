/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*			      					*
*	MODULE							*
*	    P_TYPES.H -- Type definitions for the Protein	*
*	Identification and Library Access Resource.		*
*								*
*	SYNOPSIS						*
*								*
*	DESCRIPTION						*
*								*
*	RETURNS							*
*								*
*	CAVEATS							*
*	    Not finished yet !					*
*								*
*	FILES							*
*								*
*	NOTES							*
*								*
*	SEE ALSO						*
*								*
*	DESIGNED BY						*
*	    José Ramón Valverde Carrillo (FM / UAM)		*
*								*
*	HISTORY							*
*	     8 - jul - 1989 (J R Valverde) Version #1		*
*	    26 - sep - 1989 (J R Valverde) Version #2		*
*	     8 - oct - 1989 (J R Valverde) Version #3		*
*	     8 - apr - 1990 (J R Valverde) Version #4		*
*								*
*	COPYRIGHT						*
*						    YoEgo.	*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

	    /* -- oOo -- */

#define generic	void *	 /* -- type for objects from generic packages -- */

#define info_t	void	/* -- info type to use with generic structures -- */

typedef int counter;		    /* -- counter and index for arrays -- */

typedef float distance_t;		/* -- distance between regions -- */

typedef char key_t[KEYLEN + 1];	    /* -- keyword for store with index -- */
                                                     

/*-------------------------------------------------------------------
	    -- QUEUE STRUCTURES --
-------------------------------------------------------------------*/

struct st_queue {
    counter q_spos;		    /* start position	*/
    counter q_rpos;		    /* rear position	*/
    counter q_size;		    /* queue size	*/
    info_t  **queue;		    /* the actual queue	*/
    };

typedef struct st_queue	*queue_t;   /* -- fixed-length circular queue -- */

typedef struct q_node {
    info_t	    *item;	    /* pointer to data stored in the node   */
    struct q_node   *next;	    /* pointer to next node in the queue    */
    } q_node_t;		    /* -- node of a dynamic queue -- */

typedef struct {
    int		q_size;		    /* size of the queue    */
    q_node_t	*q_front;	    /* front element	    */
    q_node_t	*q_rear;	    /* last (growing) end   */
    } *fifo_t;		    /* -- dynamic queue -- */


/*-------------------------------------------------------------------
	    -- LIST STRUCTURES --
-------------------------------------------------------------------*/
                                     
struct l_node {
    info_t	    *item;   	    /* info record stored in the node */
    struct l_node   *next;	    /* link to next node in the list  */
    };

typedef struct l_node l_node_t;	/* -- single linked list node -- */

struct l_head {
    int		l_length;	    /* number of nodes in the list    */
    l_node_t	*l_first;	    /* link to first node in the list */
    l_node_t	*l_last;	    /* link to last node 	      */
    };

typedef struct l_head l_head_t;	/* -- single linked list header info -- */
typedef l_head_t *list_t;	/* -- singly linked list	     -- */

/*--------------------------------------------------------------------
	    -- DOUBLY LINKED LIST STRUCTURES --
--------------------------------------------------------------------*/

struct dbly_node {		    /* structure of a 2-linked node	*/
    info_t	    *item;	    /* pointer to record stored in it	*/
    struct dbl_node *prev;	    /* pointer to previous node		*/
    struct dbl_node *next;	    /* pointer to successor node	*/
    };

typedef struct dbly_node dl_node_t; /* -- doubly linked list node    -- */

typedef struct {
    int		dl_length;	    /* number of nodes in the list	*/
    dl_node_t	*dl_first;	    /* pointer to the list's first node	*/
    dl_node_t	*dl_last;	    /* pointer to the list's last node	*/
    } dl_head_t;		/* -- doubly linked list header info -- */
typedef dl_head_t *dlist_t;	/* -- doubly linked list	     -- */

/*-------------------------------------------------------------------
	    -- STACK STRUCTURES --
-------------------------------------------------------------------*/

typedef struct s_node {
    info_t	    *item;	    /* information stored in this node	*/
    struct s_node   *next;	    /* underlying node in the stack    	*/
    }	s_node_t;	    /* -- node in the stack 		     -- */

typedef struct {
    int		s_length;	    /* size of the stack		*/
    s_node_t	*s_top;		    /* upper item stored in the stack	*/
    } *lifo_t;		    /* -- stack (LIFO) header information    -- */


/*-------------------------------------------------------------------
	    -- PROTEIN STRUCTURES --
-------------------------------------------------------------------*/

typedef struct {
    long    i_pos;			/* item position in DB	    */
    key_t   i_key;			/* keyword for item	    */
    } index_t;		/* -- a simple index into database	 -- */

typedef struct {
    counter	r_begin;		/* position of first aa	    */
    counter	r_length;		/* length of the region	    */
    counter 	r_offset;		/* offset to first region   */
    flag	r_status;		/* used in backtracking	    */
    } region_t;		/* -- region characteristics 		 -- */

typedef struct {
    index_t	p_index;		/* index into prot. library */
    char	p_name[PNAME_LEN + 1];	/* name of this protein	    */
    char	p_seq[MAX_SEQ_LEN + 1];	/* aminoacid sequence	    */
    int		p_seqlen;		/* sequence length	    */
    int		p_nofeats;		/* no of features measured  */
    float	**p_aaval;		/* values by feature / aa   */
    list_t 	*p_regions;	/* lists of significative regions   */
    } protein_t;	/* -- basic protein characteristics	 -- */


/*-------------------------------------------------------------------
	    -- DESCRIPTOR STRUCTURES --
-------------------------------------------------------------------*/


typedef float feature_table_t[AA_NUMBER];	/* -- ça  même -- */

typedef struct {
    string	    f_name;	    /* feature name		    */
    feature_table_t f_values;	    /* values for every aa	    */
    float	    f_maxvalue;	    /* max. amd min. values in the  */
    float	    f_minvalue;	    /* table (used after normaliz.) */
    flag	    f_inverted;	    /* invert values.		    */
    int		    f_window;	    /* window size to apply	    */
    float	    f_threshold;    /* threshold to define a region */
    int		    f_minrglen;	    /* minimum region length	    */
    int		    f_maxrglen;	    /* maximum region length	    */
    int		    f_no_regions;   /* minimum number of regions    */
    float	    f_moment;	    /* moment angle in degrees	    */
    } feature_t;

typedef struct {
    string	d_name;		    /* name of descriptor	    */
    int		d_minplen;	    /* minimum protein length	    */
    l_head_t	*d_features;	    /* the features		    */
    l_head_t	*d_rules;	    /* the rules		    */
    region_t	***d_regions;	    /* map of regions to cover	    */
    int		d_maxregno;	    /* maximum number of regions    */
    } descriptor_t;	/* -- descriptor to check		 -- */

/* The rule will be one of the following: */

/* Rule # 1: Physical separation between two regions	*/
/*  dist (F1, R1, p1) - (F2, R2, p2) = <min .. max>	*/
typedef struct {
	/* position in the first region */
	int feat1;
	int reg1;
	int pos1;
	/* position in the second region */
	int feat2;
	int reg2;
	int pos2;
	/* allowed distance limits */
	int dmin;
	int dmax;
} distance_rule_t;

/* Rule # 2: Is there a given aminoacid sequence    */
/* beginning at a specific position ?		    */
/*	Seq "sq" in (P1, R1, min .. max)	    */
typedef struct {
	string	sr_seq;
	int sr_feat;
	int sr_reg;
	int sr_minpos;
	int sr_maxpos;
} sequence_rule_t;

/* Rule # 3: How similar is its profile to the pattern region ?	*/
/*	Prof (P, R) (values,,,) < min		   		*/
typedef struct {
	int	pr_feat;
	int	pr_reg;
	int	pr_len;
	float	*pr_values;
	float	pr_mindiff;
} profile_rule_t;

typedef struct {
    int	    type;			    /* type of rule */
    flag    status;			    /* status (visited or not) */
    union {
	distance_rule_t	    distance_rule;
	sequence_rule_t	    sequence_rule;
	profile_rule_t	    profile_rule;
    } data;				    /* The actual rule */
}   rule_t;		/* -- General rule type for evaluation -- */


/*-------------------------------------------------------------------
	    -- RULE PARSER STRUCTURES --
-------------------------------------------------------------------*/


typedef struct {
	flag	tk_type;	    /* type of token */
	char	*tk_data;	    /* data of the token */
    } token_t;		/* -- token used by the parser -- */


			/* -- oOo -- */

