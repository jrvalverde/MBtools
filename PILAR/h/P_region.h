/*
 *  P_REGION.H
 *
 *	Include file containing the definitions for the module
 *  P_Region.C
 *
 *	Requires PORTABLE.H for definition of type status
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */


extern status find_regions(protein_t *, descriptor_t *);

/* MUCHO OJO A ESTO * /
extern void align_region(protein_t *, region_t *, pattern_t *, float *, int *);

extern float aa_dist(float, float);
/ * */
/*
 * Don't worry for the former comment in spanish.
 */

