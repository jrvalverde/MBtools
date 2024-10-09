/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE							*
*	    P_REGION -- region manipulation routines 		*
*								*
*	SYNOPSIS						*
*	    find_regions()                                  	*
*	    align_region()                              	*
*	    aa_dist()						*
*								*
*	DESCRIPTION						*
*	    Routines for manipulating protein regions in 	*
*	the PILAR. Contains functions for finding all regions   *
*	a given sequence of a protein ( find_regions() ),     	*
*	aligning a region with a specified pattern and          *
*	calculate its minimum distance (align_region()).    	*
*								*
*	RETURNS							*
*	    find_regions adds to the protein a list with	*
*	all the regions found.					*
*	    align_region computes the distance and position	*
*	of alignment between pattern and region.		*
*	    aa_dist returns the distance between the given	*
*	values of two aminoacids as a float.			*
*								*
*	CAVEATS			       				*
*								*
*	FILES							*
*	    Uses files:						*
*		    P_DEFS.H					*
*		    P_TYPES.H					*
*		    P_LIST.C					*
*								*
*	CAVEATS:						*
*								*
*	NOTES							*
*	    Needs information about the protein structure	*
*	as potential values per aa. in p_aaval.			*
*	It also needs to know the limit and threshold values	*
*	for defining a region.					*
*	    Debugged through other programs.			*
*	   							*
*	SEE ALSO						*
*	    P_TYPES.H -- for the definition of types used	*
*	    P_LIST.C  -- for the list manipulation funcs.	*
*								*
*	DESIGNED BY						*
*	    José Ramón Valverde Carrillo			*
*  								*
*	HISTORY							*
*	    ( 1 - jul - 1989) Implementation (IBM-PC AT).	*
*	    (17 - jul - 1989) Ported (uVAX) and debugged.	*
*	    (20 - aug - 1989) Corrected for use of feature_t.	*
*	    (30 - sep - 1989) Modified to use descriptor_t	*
*		instead of feature_t for multiple features.	*
*	    ( 1 - aug - 1991) Ported to Ultrix.			*
*	 							*
*	COPYRIGHT						*
*	    Public domain. Not allowed to get money from it	*
*	without prior written consent of the author. Otherwise	*
*	    Non eres palma, eres retama,			*
*	    Eres ciprés de triste rama,				*
*	    Eres desdicha, desdicha malal.			*
*	    Tus campos rompan tristes volcanes			*
*	    Non veas placeres, sinon pesares			*
*	    Cubran tus flores los arenales.			*
*			    (Endecha a Guillén Peraza)		*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"

#include "P_list.h"
#include "P_utils.h"

static void skip_gap(float *p_aaval, int p_seqlen);
static void calc_length(float *p_aaval, int p_seqlen);
static void annotate_region(list_t p_regions);


/*-------------------------------------------------------------------

    FIND_REGIONS

	Finds all regions of interest in the protein and stores
    them in a linked list (for now).

    INPUT:  PROT: -- The protein to analyze
	    DESC: -- The descriptor to apply
    OUTPUT: TRUE  -- Found enough regions for every feature
	    FALSE -- Not enough regions: reject protein.

    J. R. VALVERDE
-------------------------------------------------------------------*/

/* -- private variables -- */

static counter	pos = 0;	    /* current position in the sequence */
static counter	length = 0;	    /* current region length 		*/
static counter	beginning = 0;	    /* beginning of current region 	*/
static float	threshold = 0.0;    /* minimum value to consider a region */

static char find[] = "FIND REGIONS";
static char annotate[] = "ANNOTATE REGION";
static char nem[] = "NOT ENOUGH MEMORY";

/* -- the actual function -- */

status find_regions(prot, d)
protein_t *prot;
descriptor_t *d;
{
    counter i;
    l_node_t *f_node;
    feature_t *feat;
    list_t *regions;
    float **values;

    extern counter pos, length;
    extern list_t new_list();

    values = prot->p_aaval;
    regions = (list_t *) calloc(list_length(d->d_features), sizeof(list_t));
    if (regions == NULL)
	letal(find, nem);
    prot->p_regions = regions;

#ifdef V_3
    for (f_node = d->d_features->l_first, i = 0;
	f_node != NULL;
	f_node = f_node->next, ++i) {

	/* get the feature */
	feat = f_node->item;
#else
    list_reset(d->d_features);
    for (feat = (feature_t *) list_first(d->d_features), i = 0;
	feat != NULL;
	feat = (feature_t *) list_next(d->d_features), ++i) {
#endif

	pos = 0;
	threshold = feat->f_threshold;
     
	/* allocate space for regions list */
	if ( (regions[i] = new_list()) == NULL ) {
	    letal(find, nem);
    	}

	/* Find regions for this feature */

    	do {
	    /* skip non significative values */
	    skip_gap(values[i], prot->p_seqlen);

	    /* compute current region length */
	    calc_length(values[i], prot->p_seqlen);

	    /* if lengthy enough and not too short:  */
	    if ((length >= feat->f_minrglen) && (length <= feat->f_maxrglen))
		/* add it to the list.	     */
        	annotate_region(regions[i]);
    	}
    	while (pos < prot->p_seqlen);

	/* check minimum number of regions */
#ifdef V_3
	if (regions[i]->l_length < feat->f_no_regions)
#else
	if (list_length(regions[i]) < feat->f_no_regions)
#endif
	    /* There is no sense in continuing.
	    This protein is not worth enough */
	    return FALSE;
    }

    /* We got all regions and enough of them all */
    return TRUE;
}

static void skip_gap(p_aaval, p_seqlen)
float *p_aaval;
int p_seqlen;
/*
 *   Skips non significative values and sets counter into protein to the
 *   first position whose value is equal or greater than the threshold.
 */
{
    extern counter pos;
    extern float threshold;

    while ((p_aaval[pos] < threshold) && (pos < p_seqlen) )
	pos++;
}

static void calc_length(p_aaval, p_seqlen)
float *p_aaval;
int p_seqlen;
/*
 *   runs along the sequence until it finds the end of the region
 *   and sets the values of its beginning, end and length.
 */
{
    extern counter pos, beginning, end, length;
    extern float threshold;

    beginning = pos;
    while ( (p_aaval[pos] >= threshold) && (pos < p_seqlen) )
	pos++;
    length = pos - beginning;
}

static void annotate_region(p_regions)
list_t p_regions;
/*
 *   Add region to the list of found regions.
 */
{
    region_t *prgn;	    /* pointer to new region */
    region_t *pfirst;	    /* pointer to the first region */

    /* allocate space for a region record */
    if ((prgn = NEW(region_t)) == NULL) {
	letal(annotate, nem);
      }

    /* fill in the form */              
    prgn->r_begin = beginning;
    prgn->r_length = length;
    prgn->r_status = ~ASSIGNED;

    /* offset is displacement of the beginning of the current region
    with respect to the beginning of the first one. */
#ifdef V_3
    if (p_regions->l_length == 0)
	prgn->r_offset = 0;
    else
	prgn->r_offset = beginning
	    - ((region_t *) (p_regions->l_first->item))->r_begin;
#else
    if (list_length(p_regions) == 0)
	prgn->r_offset = 0;
    else {
	pfirst = (region_t *) list_first(p_regions);
	prgn->r_offset = beginning - pfirst->r_begin;
    }
#endif

    /* store it at the end of the list */
    if ( list_append(prgn, p_regions) ) {
	return;
    }
    else {
	letal(annotate, nem);
    }
}


/*-------------------------------------------------------------------
    ALIGN_REGION
    Finds the minimum distance between a given region and the
    specified pattern.

    INPUT:  PROTEIN -- the protein in which the region is
	    REGION -- the region whose distance we need
	    PATTERN -- the pattern to compare with
    OUTPUT: DIST -- minimum distance computed (PLUS_INFINITY if
		no alignment is possible).
	    OFFSET -- offset of the region with respect to the
		specified pattern (NOT_ALIGNED if no align-
		ment is possible).

    J. R. VALVERDE
-------------------------------------------------------------------*/

#ifdef ALIGN_REGION
/* Cambiar definicion a float align_region(reg, pat) */

void align_region(protein, region, pattern, dist, offset)
protein_t *protein;
region_t  *region;
pattern_t *pattern;
float     *dist;
int       *offset;
{
    counter pat_beg, aa_offset;
    float tot_dist;
    extern float aa_dist();

    *dist = F_INFINITY;	/* maximum float */
    *offset = NOT_ALIGNED;
    /* check to see if region length is out of bounds */
    if ( (region->r_length < pattern->pt_lowlen) ||
	 (region->r_length > pattern->pt_upplen) ) {
	return;
    }

    /* if not, then region length is less or equal than pattern length;
       offset aa by aa the region over the pattern and compute distance
       storing the alignment with lesser distance.	    */
    pat_beg = 0;
    do {
	tot_dist = 0.0;
	for (aa_offset = 0; aa_offset < region->r_length; aa_offset++) {
            tot_dist += aa_dist(pattern->pt_val[pat_beg + aa_offset],
		protein->p_aaval[region->r_begin + aa_offset]);
	}
	tot_dist = sqrt(tot_dist / region->r_length);
	if (tot_dist < *dist) {
      	    *offset = aa_offset;
      	    *dist = tot_dist;
	}
	pat_beg++;
    }
    while (pat_beg + region->r_length < pattern->pt_length);
}


/*-------------------------------------------------------------------
    AA_DIST

    Returns the distance between two single amino acid associated
    values.

    INPUT:  VAL1, VAL2: the two values
    OUTPUT: The computed distance.

    J. R. VALVERDE
-------------------------------------------------------------------*/

float aa_dist(val1, val2)
float val1, val2;
{
    float v1, v2;

    /* Check for zero values to avoid division by zero
    and replace with small constants.	*/
    v1 = (val1 == 0.0)? 1e-5 : val1;
    v2 = (val2 == 0.0)? 1e-5 : val2;

    return ((v1 - v2) * (v1 - v2)) / (v1 * v2);
}

#endif
