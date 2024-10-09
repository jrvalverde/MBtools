/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE:							*
*	    P_THINK.C : Intelligence routines for use in PILAR	*
*	(Protein Identification and Library Access Resource)	*
*								*
*	SYNOPSIS:						*
*	    This module contains the routines used to check the	*
*	proteins in the database and assert the rules specified	*
*	by the user and contained in a descriptor.		*
*								*
*	DESCRIPTION:						*
*	    This module contains the following subroutine:	*
*		assert_rules()					*
*								*
*	CAVEATS:						*
*	    Not fully debugged. Distance rules seem to work Ok.	*
*	Other rules are still unchecked.			*
*								*
*	NOTES:							*
*	    Right now the subroutine only follows all the rules	*
*	and returns TRUE iff all the rules are fulfilled. It is	*
*	my wish to modify the algorithm used for processing of	*
*	AND/OR graphs.						*
*	    To add a new rule you need to create a function	*
*	to evaluate it and insert it in the CASE construct of	*
*	function assert_rule().					*
*								*
*	    IN ORDER TO KEEP TRACK OF THE RECURSION WE NEED	*
*	TO KNOW OF THE INTERNALS OF A LIST STRUCTURE IN THIS	*
*	PACKAGE. USING LIST CALLS WOULD RESULT IN NOT BEING	*
*	ABLE TO	BACKTRACK CORRECTLY.				*
*								*
*	FILES:							*
*	    STDIO.H						*
*	    PORTABLE.H						*
*	    P_DEFS.H						*
*	    P_TYPES.H						*
*	    P_EXTERN.H						*
*								*
*	SEE ALSO:						*
*								*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo.			*
*								*
*	HISTORY:						*
*	    jun/sept - 89 (J R VALVERDE) Various approaches	*
*	    10 - sep - 89 (J R VALVERDE) Written		*
*	    16 - sep - 89 (J R VALVERDE) Implemented		*
*	     6 - oct - 89 (J R VALVERDE) Correct assert_dist	*
*	     5 - nov - 89 (J R VALVERDE) Support for seq & prof	*
*	     1 - aug - 91 (J R VALVERDE) Ported to Ultrix	*
*	    30 - sep - 91 (J R VALVERDE) Moved init_rules()	*
*				here from P_init.c		*
*								*
*	COPYRIGHT:						*
*	    This is Public Domain software. Please note that	*
*	you cannot use this software for commercial purposes	*
*	without prior written consent from the author.		*
*							YoEgo	*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>
#include <ctype.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"

#include "P_list.h"

#define _V_5_

#define DISTANCE    1
#define SEQUENCE    2
#define PROFILE	    3

static protein_t *prot;
static descriptor_t *dscr;
static l_head_t **regions;
static float **values;
static region_t ***rmap;

extern void init_rules(list_t);

static int assert_rule(l_node_t *p_rule);
static int assert_distance(l_node_t *p_rule);
static int check_dist(region_t *reg_1, l_node_t *p_rule);
static int assert_sequence(l_node_t *p_rule);
static int assert_profile(l_node_t *p_rule);

/*
    ... et comme il n'est aucune science qui requiert une application
    plus grande et sincère que celle cy, soit fort soigneux de lire et
    relire souvent les enseignments que je te donne, àfin qu'êtant
    gravés en ta mémoire, tu ne sois pas capable de faire des manquements
    qui te pouvoient être prejudiciables, au lieu que prenant tous les
    soins nècessaires, la science te deviendra autant facile quelle te
    paroitra difficille au premier abord, et pour çe moyen de terrestre
    et opaque que tu es, tu deviendras subtil et clarivident par la
    communication des susdits esprits.

		--- Les clavicules de Salomon ---
*/

/*--------------------------------------------------------------------

    ASSERT_RULES

	This function checks that a protein adapts to the characteristics
    specified by a given descriptor.
	It makes a depth first traversal of a graph containing all the
    constraints imposed on the protein structure. The constraints are
    expressed as a list of rules ordered in the same order that the
    graph should be traversed.
	It will traverse unconnected graphs composed of several chunks
    or bursts of rules independent one group from the other.
	Upon exit the regions that have been considered will be tagged
    as ASSIGNED in their r_status field. On entry these tags MUST BE all
    set to ~ASSIGNED, since the routine does NOT CHECK that.

    INPUT:  P -- The protein to test.
	    D -- The descriptor to check against.

    OUTPUT: TRUE if the protein can be accepted as modelled by the
	    description specified.
	    FALSE if not.

    J. R. VALVERDE

---------------------------------------------------------------------*/

boolean assert_rules(p, d)
protein_t *p;
descriptor_t *d;
{
    l_node_t *p_feat, *p_rule;
    feature_t *feat;
    rule_t *the_rule;
    int i, j, k, m;

    extern protein_t *prot;
    extern descriptor_t *dscr;
    region_t ***RegMap;
    extern region_t ***rmap;	/* rmap is a bidimensional array of pointers
				to regions: [feat# x maxreg#]. It is used
				to keep track of assigned and unassigned
				regions during the checkout process. */

    prot = p;
    dscr = d;
    regions = prot->p_regions;
    values = prot->p_aaval;
    rmap = dscr->d_regions;
    /* initialize rmap to be empty */
#ifdef _V_4_
    for (i = 0, p_feat = dscr->d_features->l_first;
	p_feat != NULL;
	p_feat = p_feat->next, i++) {
	/* for every feature */
	feat = (feature_t *) p_feat->item;
	for (j = 0; j < feat->f_no_regions; ++j)
	    /* set all needed regions to NULL */
	    rmap[i][j] = NULL;
    }
#else
    k = dscr->d_features->l_length;	    /* Number of features */
    m = dscr->d_maxregno;		    /* maximum number of regions */
#ifdef _V_5_ 
   RegMap = rmap;

    for (i = 0; i < k; ++i)
	memset(rmap[i], 0, m * sizeof(region_t *));
#else
    for (i = 0; i < k; ++i)		    /* for all the rmap */
	for (j = 0; j < m; ++j)
	    rmap[i][j] = NULL;		    /* set all rmap to NULL */
#endif
#endif

    /* initialize rules graph */                                 
    init_rules(dscr->d_rules);

    /* For every rule */
    for (p_rule = dscr->d_rules->l_first, i = dscr->d_rules->l_length;
	 (p_rule != NULL) && i;
	 p_rule = p_rule->next, --i) {
	the_rule = (rule_t *) p_rule->item;
	if (the_rule->status == EXPLORED)
	    continue;
	/* this avoids problems with unconnected graphs */
	if (assert_rule(p_rule) == ~OK)
	    return FALSE;
    }
    return TRUE;	    /* iff every rule is OK */
}

/*------------------------------------------------------------------*/

static int assert_rule(p_rule)
l_node_t *p_rule;
/*
 *  ASSERT RULE
 *	This functions selects the appropriate one for testing
 *  the actual rule according to its type. It first checks the
 *  type and then calls the actually checking function.
 *
 *  INPUT: The rule to check
 *  OUTPUT: The result of the check.
 *
 *  J. R. Valverde
 */
{
    extern int assert_disance();
    extern int assert_sequence();
    extern int assert_profile();
    extern void horror(char *, char *);

    /* If there are no more rules, then we are done */
    if (p_rule == NULL)
	return (OK);

    /* If this path has already been explored... */
    if ( ((rule_t *)p_rule->item)->status == EXPLORED )
	return (OK);

    switch ( ((rule_t *)p_rule->item)->type ) {

	case DISTANCE:
	    return (assert_distance(p_rule));

	case SEQUENCE:
	    return (assert_sequence(p_rule));

	case PROFILE:
	    return (assert_profile(p_rule));
	default:
	    /* Let the user know something is going wrong */
	    horror("ASSERT_RULE:", "Unknown rule");
	    /* and skip this anomalous rule */
	    return ( assert_rule(p_rule->next) );
    }
}

/*----------------------------------------------------------------*/
/*
    "It seems very pretty," she said when she had finished it,
    "but it's RATHER hard to understand!"

    Lewis Carroll. Through the looking glass and what Alice found there.
*/
/*----------------------------------------------------------------*/

extern int check_dist(region_t *, l_node_t *);

static int assert_distance(p_rule)
l_node_t *p_rule;
/*
 *  ASSERT_DISTANCE
 *
 *	This function tries to find two regions of defined characteristics
 *  that fall within the physical separation restraints imposed by the
 *  rule. First it tries to find the first region and then it calls
 *  check_distance() to find the second.
 *
 *  INPUT:  The rule to test
 *  OUTPUT: TRUE if this rule and its childs are all verified.
 *	    FALSE if any of them fails.
 *
 *  J. R. Valverde
 */
{
    rule_t *the_rule;
    distance_rule_t *r;
    region_t *reg_1;
    l_node_t *p_reg;
    extern region_t ***rmap;

    the_rule = (rule_t *) p_rule->item;
    r = &(the_rule->data.distance_rule);

    /* Has this region been assigned before ? */
    if (rmap[r->feat1][r->reg1] != NULL) {
	/* If reg_1 has already been assigned */
	reg_1 = rmap[r->feat1][r->reg1];

	/* look for a suitable reg_2 */
	if (check_dist(reg_1, p_rule) == OK)
	    return (OK);
	else
	    return(~OK);
    }

    else {
	/* Fi Rj is not already known	*/
	/*  Find one:			*/
	/*	for every Rj in Fi	*/
	for (p_reg = regions[r->feat1]->l_first;
	     p_reg != NULL;
	     p_reg = p_reg->next) {

	    reg_1 = (region_t *) p_reg->item;

	    if (reg_1->r_status == ASSIGNED)
		/* This region is assigned to any other in the pattern */
		continue;

	    rmap[r->feat1][r->reg1] = reg_1;
	    reg_1->r_status = ASSIGNED;

	    /* look for a suitable reg_2 */
	    if (check_dist(reg_1, p_rule) == OK)
		return (OK);
	    else {
		/* try another reg_1 */
		rmap[r->feat1][r->reg1] = NULL;
		reg_1->r_status = ~ASSIGNED;
		continue;
	    }
	}	/* end of for loop */

	/* if we arrive here then reg_1 == NULL */
	/* there is no such a region number 1: the rule has failed */
	rmap[r->feat1][r->reg1] = NULL;
	the_rule->status = ~EXPLORED;
	return (~OK);
    }
}

/*-----------------------------------------------------------------*/

static int check_dist(reg_1, p_rule)
region_t *reg_1;
l_node_t *p_rule;
/*
 *  CHECK_DIST
 *
 *	Given an initial region this function tries to find a
 * second one of the given characteristics within the limits of
 * separation imposed by the rule.
 *	If it can find the second one, then it recursively checks
 *  the following child rules before definitely accepting the rule.
 *
 *  INPUT:  The first region and the rule describing the constraints
 *	to find the second region.
 *
 *  J. R. Valverde
 */
{
    int distance;
    int p1, p2;
    rule_t *the_rule;
    distance_rule_t *r;
    region_t *reg_2;
    l_node_t *pointer;
    extern region_t ***rmap;
    extern protein_t *prot;

    /* we already have a possible reg_1 */
    /* now look for reg_2 */
    the_rule = (rule_t *) p_rule->item;
    r = &(the_rule->data.distance_rule);

    if ((reg_2 = rmap[r->feat2][r->reg2]) != NULL) {

	/* if reg_2 has already been assigned then
	use it to check the distance */
        if (r->pos1 == -1)  p1 = reg_1->r_begin + reg_1->r_length;
	else		    p1 = reg_1->r_begin + r->pos1;
	if (r->pos2 == -1)  p2 = reg_2->r_begin + reg_2->r_length;
	else		    p2 = reg_2->r_begin + r->pos2;
	distance = p2 - p1;
	if ((distance < r->dmin) || (r->dmax < distance))
	    /* if out of limits we have to backtrack and try
	    another set of regions */
	    return (~OK);
	else {
	    /* check following rules */
	    the_rule->status = EXPLORED;
	    if (assert_rule(p_rule->next) == OK) {
		/* following rules are OK (no conflicts) */
		return (OK);
	    }
	    else {
		/* following rules do conflict with this choice */
		the_rule->status = ~EXPLORED;
		return (~OK);
	    }
	}   
    } /* end IF REG ALREADY ASSIGNED */

    else {
	/* Search for a possible reg2 */
	for (pointer = regions[r->feat2]->l_first;
	     pointer != NULL;
	     pointer = pointer->next) {

	    reg_2 = (region_t *) pointer->item;
	    if (reg_2->r_status == ASSIGNED)
		/* this one is not available */
		continue;

	    /* IF REG i,j IS NOT IN USE CHECK RESTRICTION */
	    if (r->pos1 == -1)	p1 = reg_1->r_begin + reg_1->r_length;
	    else		p1 = reg_1->r_begin + r->pos1;
	    if (r->pos2 == -1)	p2 = reg_2->r_begin + reg_2->r_length;
	    else		p2 = reg_2->r_begin + r->pos2;
	    distance = p2 - p1;
	    if ((r->dmin <= distance) && (distance <= r->dmax)) {
		/* Bingo !! : set rmap and status */
		rmap[r->feat2][r->reg2] = reg_2;
		reg_2->r_status = ASSIGNED;
		the_rule->status = EXPLORED;
		/* check following rules */
		if (assert_rule(p_rule->next) == OK) {
		    return (OK);
		}
		else {
		    /* This region does not fit since the following
		    rules have found conflicts. Try another one. */
		    reg_2->r_status = ~ASSIGNED;
		    rmap[r->feat2][r->reg2] = NULL;
		    the_rule->status = ~EXPLORED;
		    continue;
		}
	    }
	    else {
		/* this region is not acceptable, try another one */
		continue;
	    }
	} /* end FOR EVERY REGION j IN FEATURE i */

	/* No more remaining regions: */
	/* There is no such reg_2 under these constraints for that first one */
	rmap[r->feat2][r->reg2] = NULL;
	the_rule->status = ~EXPLORED;	    
	return (~OK);

    }  /* end IF REG_2 NOT ALREADY ASSIGNED */

}


/*-----------------------------------------------------------------*/

static int assert_sequence(p_rule)
l_node_t *p_rule;
/*
 *  Stub for the function assert_sequence.
 *  It only ignores this rule.
 */
{
    int pos, end_pos, len, result;
    rule_t *the_rule;
    sequence_rule_t *r;
    region_t *reg;
    l_node_t *pointer;
    extern region_t ***rmap;
    extern protein_t *prot;

    /* Look for reg */
    the_rule = (rule_t *) p_rule->item;
    r = &(the_rule->data.sequence_rule);

    if ((reg = rmap[r->sr_feat][r->sr_reg]) != NULL) {

	/* if reg has already been assigned then
	use it to check the rule:
	for each possible position in the region
	compare with pattern secuence: */

	if (r->sr_maxpos == -1)
	/* we must finish at the END of the region */
	    end_pos = reg->r_begin + reg->r_length;
	else
	    end_pos = reg->r_begin + r->sr_maxpos;

	pos = reg->r_begin + r->sr_minpos;
	if (pos < 0) pos = 0;		/* test for extreme minpos & maxpos */
	if (end_pos > prot->p_seqlen) end_pos = prot->p_seqlen;
	len = strlen(r->sr_seq);
	for (;
	     pos < end_pos;
	     pos++) {
	    result = strncmp(prot->p_seq + pos, r->sr_seq, len);
            if (result == 0) break;
	}

	if (result != 0)
	    /* if the sequence was not found then we have to
	    backtrack and try another region */
	    return (~OK);
	else {
	    /* check following rules */
	    the_rule->status = EXPLORED;
	    if (assert_rule(p_rule->next) == OK) {
		/* following rules are OK (no conflicts) */
		return (OK);
	    }
	    else {
		/* following rules do conflict with this choice */
		the_rule->status = ~EXPLORED;
		return (~OK);
	    }
	}   
    } /* end IF REG ALREADY ASSIGNED */

    else {
	/* Search for a possible reg */
	for (pointer = regions[r->sr_feat]->l_first;
	     pointer != NULL;
	     pointer = pointer->next) {

	    reg = (region_t *) pointer->item;
	    if (reg->r_status == ASSIGNED)
		/* this one is not available */
		continue;

	    /* IF REG i,j IS NOT IN USE CHECK RESTRICTION */
	    if (r->sr_maxpos == -1)
	    /* we must finish at the END of the region */
		end_pos = reg->r_begin + reg->r_length;
	    else
		end_pos = reg->r_begin + r->sr_maxpos;

	    len = strlen(r->sr_seq);
	    pos = reg->r_begin + r->sr_minpos;
	    if (pos < 0) pos = 0;		/* test for extreme positions */
	    if (end_pos > prot->p_seqlen) end_pos = prot->p_seqlen;
	    for (;
		 pos < end_pos;
	 	 pos++) {
		result = strncmp(prot->p_seq + pos, r->sr_seq, len);
            	if (result == 0) break;
	    }

	    if (result == 0) {
		/* Bingo !! : set rmap and status */
		rmap[r->sr_feat][r->sr_reg] = reg;
		reg->r_status = ASSIGNED;
		the_rule->status = EXPLORED;
		/* check following rules */
		if (assert_rule(p_rule->next) == OK) {
		    return (OK);
		}
		else {
		    /* This region does not fit since the following
		    rules have found conflicts. Try another one. */
		    reg->r_status = ~ASSIGNED;
		    rmap[r->sr_feat][r->sr_reg] = NULL;
		    the_rule->status = ~EXPLORED;
		    continue;
		}
	    }
	    else {
		/* this region is not acceptable, try another one */
		continue;
	    }
	} /* end FOR EVERY REGION j IN FEATURE i */

	/* No more remaining regions: */
	/* There is no such reg under these constraints */
	rmap[r->sr_feat][r->sr_reg] = NULL;
	the_rule->status = ~EXPLORED;	    
	return (~OK);

    }  /* end IF REG NOT ALREADY ASSIGNED */
}

/*----------------------------------------------------------------*/

static int assert_profile(p_rule)
l_node_t *p_rule;
/*
 *  Stub for the function assert_profile.
 *  It only ignores this rule.
 */
{
    return assert_rule(p_rule->next);

#ifdef LO_QUE_YO_TE_DIGA
    int pos, dist, min_dist;     
    rule_t *the_rule;
    sequence_rule_t *r;
    region_t *reg_2;
    l_node_t *pointer;
    extern region_t ***rmap;
    extern protein_t *prot;

    /* Look for reg */
    the_rule = (rule_t *) p_rule->item;
    r = &(the_rule->data.profile_rule);

    if ((reg = rmap[r->pr_feat][r->pr_reg]) != NULL) {

	/* if reg has already been assigned then
	use it to check the rule: */

	min_dist = F_INFINITY;
	for (pos = reg->r_begin;
	     pos < reg->r_begin + reg->r_length - r->pr_len;
	     pos++) {
	    dist = cmp_profiles(prot->p_values[r->pr_feat]+pos, r->pr_values);
            if (dist < min_dist) min_dist = dist;
	}

	if (min_dist >= r->pr_mindiff)
	    /* if the difference was too high then we have to
	    backtrack and try another region */
	    return (~OK);
	else {
	    /* check following rules */
	    the_rule->status = EXPLORED;
	    if (assert_rule(p_rule->next) == OK) {
		/* following rules are OK (no conflicts) */
		return (OK);
	    }
	    else {
		/* following rules do conflict with this choice */
		the_rule->status = ~EXPLORED;
		return (~OK);
	    }
	}   
    } /* end IF REG ALREADY ASSIGNED */

    else {                                   
	/* Search for a possible reg */
	for (pointer = regions[r->sr_feat]->l_first;
	     pointer != NULL;
	     pointer = pointer->next) {

	    reg = (region_t *) pointer->item;
	    if (reg->r_status == ASSIGNED)
		/* this one is not available */
		continue;

	    /* IF REG i,j IS NOT IN USE CHECK RESTRICTION */

	    min_dist = F_INFINITY;
	    for (pos = reg->r_begin;
		 pos < reg->r_begin + reg->r_length - r->pr_len;
	 	 pos++) {
	    dist = cmp_profiles(prot->p_values[r->pr_feat]+pos, r->pr_values);
            if (dist < min_dist) min_dist = dist;
	    }

	    if (min_dist < r->pr_mindiff)
		/* Bingo !! : set rmap and status */
		rmap[r->sr_feat][r->sr_reg] = reg;
		reg->r_status = ASSIGNED;
		the_rule->status = EXPLORED;
		/* check following rules */
		if (assert_rule(p_rule->next) == OK) {
		    return (OK);
		}
		else {
		    /* This region does not fit since the following
		    rules have found conflicts. Try another one. */
		    reg->r_status = ~ASSIGNED;
		    rmap[r->sr_feat][r->sr_reg] = NULL;
		    the_rule->status = ~EXPLORED;
		    continue;
		}
	    }
	    else {
		/* this region is not acceptable, try another one */
		continue;
	    }
	} /* end FOR EVERY REGION j IN FEATURE i */

	/* No more remaining regions: */
	/* There is no such reg under these constraints */
	rmap[r->sr_feat][r->sr_reg] = NULL;
	the_rule->status = ~EXPLORED;	    
	return (~OK);

    }  /* end IF REG NOT ALREADY ASSIGNED */
#endif /* LO_QUE_YO_TE_DIGA */

}

/*--------------------------------------------------------------------
    INIT_RULES
	Initialize status of every rule in the rules graph. Since now
    we play with a plain set of rules with no defined topography nor
    relationships, we only have to deal with a single rules list.
	This function simply sets the status of every rule to unexplored
    in order for the graph traversal algorithm to begin with a clean
    set of rules.

    INPUT:  RULES -- A list of rules
    OUTPUT: The list has been initialized.

    J. R. VALVERDE
--------------------------------------------------------------------*/

static void init_rule_status( r )
rule_t *r;
{
    r->status = ~EXPLORED;
}

void init_rules(rules)
list_t rules;
{
    extern void list_traverse(list_t, void (*)());

#ifdef V_3
    rule_t *r;		    /* the actual rule */
    l_node_t *node;	    /* the node in the list */

    for (node = rules->l_first;
	node != NULL;
	node = node->next) {
	r = (rule_t *) node->item;
	r->status = ~EXPLORED;
    }
#else
    list_traverse( rules, init_rule_status );
#endif
}

/*------------------------------------------------------------------*/


/*
    Aieth Gadol Leolaus Adonai

    "Adonai (the Lord) will be Great in the Ethernity".
*/
