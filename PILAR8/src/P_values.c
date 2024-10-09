/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE							*
*	    P_VALUES -- assign values to a sequence and		*
*	find regions of interest.				*
*								*
*	SYNOPSIS						*
*	    assign_values()					*
*								*
*	DESCRIPTION						*
*	    This module contains the functions used for		*
*	assigning the appropriate values to a sequence accor-	*
*	ding to a feature values table (assign_values()) and	*
*	for averaging the values so obtained using a window of	*
*	predefined length without considering (window_soften())	*
*	and considering a moment (moment_soften(),		*
*	AVmoment_soften).					*
*								*
*	RETURNS							*
*								*
*	CAVEATS							*
*								*
*	FILES							*
*	    stdio.h						*
*	    math.h						*
*	    float.h						*
*	    portable.h						*
*	    p_defs.h						*
*	    p_types.h						*
*	    p_extern.h						*
*								*
*	NOTES							*
*								*
*	SEE ALSO						*
*								*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo			*
*								*
*	HISTORY							*
*	    José R. Valverde (30 - jun - 1989)			*
*		-- initial implementation.			*
*	    J. R. Valverde (12-jul-89) -- Debugged.		*
*	    J. R. Valverde (16-jul-89) -- Ported to Vax.	*
*	    J. R. Valverde (20-aug-89) -- Use of feature_t.	*
*	    J. R. Valverde (30-sep-89) -- Use descriptor_t	*
*	    instead of feature_t to get multiple region lists.	*
*	    J. R. Valverde (07-oct-89) -- Window softening	*
*	    subroutine rewritten unsing queue techniques to	*
*	    decrease the number of floating point operations.	*
*	    J. R. Valverde (08-apr-90) -- Added moment routines	*
*	    J. R. Valverde (01-jul-91) -- Removed AV's routines	*
*	    J. R. Valverde (01-aug-91) -- Ported to Ultrix	*
*								*
*	COPYRIGHT						*
*	    This is initially public domain software.		*
*	Please note that if you want to earn any profit from 	*
*	it you should contact the designer.			*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#ifndef DEBUG
/* #define DEBUG    3 */
#define DEBUG	    -1
#endif

#include <stdio.h>
#include <ctype.h>
#include <math.h>

/* #include <float.h> */

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"

#include "P_list.h"
#include "P_utils.h"

#define PI	3.1415926538979323846
#define SQR(x)  ((x) * (x))

extern int aa_table[];
/*
void moment_soften(float *, int , int , float);
*/
void moment_soften();

/*-------------------------------------------------------------------
    ASSIGN_VALUES

    	Assign values to a protein's sequence based on a feature.

    INPUT:  PROT -- The protein to wich we wish to assign values
	    D -- The descriptor to measure in the protein
    OUTPUT: Protein values field filled with the corresponding
	values for every aminoacid and every feature to consider.

    J. R. VALVERDE
---------------------------------------------------------------------*/

static float *aa_values;

void assign_values(prot, d)
protein_t *prot;
descriptor_t *d;
{
    counter i, f_number;
    l_node_t *f_node;
    feature_t *feat;
    float **values;
    extern void window_soften();

    /* allocate auxiliary array for temporal storage of
    intermediate not softened values */
    aa_values = (float *) calloc(prot->p_seqlen, sizeof(float));
    if (aa_values == NULL) {
	letal("ASSIGN VALUES", "Not enough memory");
    }

    /* allocate array of arrays of values */
#ifdef V_3
    values = (float **) calloc(i=d->d_features->l_length,
				    sizeof(float *));
#else
    values = (float **) calloc(list_length(d->d_features),
				    sizeof(float *));
#endif
    if (values == NULL) {
	letal("ASSIGN VALUES", "Not enough memory");
    }

#ifdef V_3
    prot->p_nofeats = d->d_features->l_length;
#else
    prot->p_nofeats = list_length(d->d_features);
#endif

    /* Run over every feature obtaining its values */
#ifdef V_3
    for (f_node = d->d_features->l_first, f_number = 0;
	f_node != NULL;
	f_node = f_node->next, f_number++) {

	/* get the feature */
	feat = f_node->item;
#else
    list_reset(d->d_features);
    for (feat = (feature_t *) list_first(d->d_features), f_number = 0;
	feat != NULL;
	feat = (feature_t *) list_next(d->d_features), f_number++) {
#endif
	/* assign the values */
	for (i = 0; i < prot->p_seqlen; i++)
       	    aa_values[i] = feat->f_values[ aa_table[ prot->p_seq[i] ] - 1 ];

	/* and soften them */
	values[f_number] = calloc(prot->p_seqlen, sizeof(float));
	if (values[f_number] == NULL)
	    letal("ASSIGN VALUES", "Not enough memory");
	if (feat->f_moment < 0.0)	    /* moment must be > 0.0 */
	    window_soften(  values[f_number],
			    prot->p_seqlen,
			    feat->f_window );
	else
	    moment_soften(  values[f_number],
			    prot->p_seqlen,
			    feat->f_window,
			    feat->f_moment );
    }
    prot->p_aaval = values;

    free(aa_values);
}

/*---------------------------------------------------------------------
    WINDOW_SOFTEN

	Promediates the values of every aa with its neighbours.

    INPUT:  P_AAVAL -- the array that will receive the values.
	    P_SEQLEN -- length of the sequence whose values are to be softened.
	    WSIZE --  the size of window to use.
    OUTPUT: The values for every aa. have been promediated.

    J. R. VALVERDE
---------------------------------------------------------------------*/

void window_soften(p_aaval, p_seqlen, wsize)
float *p_aaval;
int p_seqlen;
int wsize;
{    
    int whalf_1, whalf_2;
    float sum;
    int rear, front, i;
    extern float *aa_values;

    whalf_2 = whalf_1 = wsize >> 1;
    if ((wsize & 1) == 0) --whalf_2;

    rear = front = i = 0;
    sum = 0.0;

    /* first aa: setup initial sum */
    /* first value plus a half window of values */

    for (front = 0; front <= whalf_2; front++) {
	sum += aa_values[front];
    }
    p_aaval[0] = sum / front;
    /* front equals now whalf_2 + 1 which is also the number of
    values considered */

    /* first half window: rear is still 0 */

    for (i = 1; i <= whalf_1; ++i) {
	sum += aa_values[front++];
	p_aaval[i] = sum / front;
    }
    /* We have now one half window on the left, with
    front pointng one value beyond the right half window
    (and, since we are 0-offset, equal to the number of
    values considered, which should be one full window) */

    /* The bulk of values */

    for (; i < (p_seqlen - whalf_2); ++i) {
	sum += aa_values[front++];	/* add next value */
	sum -= aa_values[rear++];	/* and subtract the leftish one */
	p_aaval[i] = sum / wsize;
    }
    /* At the end of the for block front always points one
    value beyond the rightest one, wich is the next value to
    be added and rear to the first one in the just considered
    window, which is the next to be eliminated */

    /* And the last half window: front won't increment */

    for (; i < p_seqlen; ++i) {
	sum -= aa_values[rear++];
	p_aaval[i] = sum / (p_seqlen - rear);
    }
    /* We keep subtracting the tracking value, but now we cannot
    add more values. After the subtraction rear points to the first
    value in the window. Since we are 0-offset, subtracting it from
    the number of values yields the window size	*/

}

/*---------------------------------------------------------------------
    MOMENT_SOFTEN

	Promediates the values of every aa with its neighbours taking
    into account the angle of each residue with respect to the peptide
    backbone (computing hence its moment, as described by Eisenberg et
    al. J.Mol.Biol. 179, 125-142; Finer-Moore & Stroud, P.N.A.S. 81,
    155-159).

    INPUT:  P_AAVAL -- the array that will receive the values.
	    P_SEQLEN -- length of the sequence whose values are to be softened.
	    WSIZE --  the size of window to use.
    OUTPUT: The values for every aa. have been promediated.

    J. R. VALVERDE
---------------------------------------------------------------------*/

void moment_soften(p_aaval, p_seqlen, wsize, degrees)
float *p_aaval;
int p_seqlen;
int wsize;
float degrees;
{    
    int whalf_1, whalf_2;
    float sumcos, sumsin, 
	delta, 
	frontdelta, reardelta;
    int rear, front, i;
    extern float *aa_values;

    whalf_2 = whalf_1 = wsize >> 1;
    if ((wsize & 1) == 0) --whalf_2;

    rear = front = i = 0;

    delta = degrees * 2.0 * PI / 360.0;
    frontdelta = reardelta = sumcos = sumsin = 0.0;

    /* first aa: setup initial sum */
    /* first value plus a half window of values */

    for (front = 0; front <= whalf_2; front++) {
	frontdelta += delta;
	sumcos += aa_values[front] * cos(frontdelta);
	sumsin += aa_values[front] * sin(frontdelta);
    }
    p_aaval[0] = sqrt( SQR(sumcos) + SQR(sumsin) ) / front;
    /* front equals now whalf_2 + 1 which is also the number of
    values considered */

    /* first half window: rear is still 0 */

    for (i = 1; i <= whalf_1; ++i) {
	frontdelta += delta;
	sumcos += aa_values[front] * cos(frontdelta);
	sumsin += aa_values[front] * sin(frontdelta);
	p_aaval[i] = sqrt( SQR(sumcos) + SQR(sumsin) ) / ++front;
    }
    /* We have now one half window on the left, with
    front pointing one value beyond the right half window
    (and, since we are 0-offset, equal to the number of
    values considered, which should be one full window) */

    /* The bulk of values */
    for (; i < (p_seqlen - whalf_2); ++i) {
	/* compute deltas */
	frontdelta += delta; reardelta += delta;
	/* add next value */
	sumcos += aa_values[front] * cos(frontdelta);
	sumsin += aa_values[front] * sin(frontdelta);
	/* and subtract the leftish one */
	sumcos -= aa_values[rear] * cos(reardelta);
	sumsin -= aa_values[rear] * sin(reardelta);
	front++; rear++;
	p_aaval[i] = sqrt( SQR(sumcos) + SQR(sumsin) ) / wsize;
    }
    /* At the end of the for block front always points one
    value beyond the rightest one, wich is the next value to
    be added and rear to the first one in the just considered
    window, which is the next to be eliminated */

    /* And the last half window: front won't increment */

    for (; i < p_seqlen; ++i) {
	reardelta += delta;
	sumcos -= aa_values[rear] * cos(reardelta);
	sumsin -= aa_values[rear] * sin(reardelta);
	rear++;
	p_aaval[i] = sqrt(SQR(sumcos) + SQR(sumsin)) / (p_seqlen - rear);
    }
    /* We keep subtracting the tracking value, but now we cannot
    add more values. After the subtraction rear is set to the first
    value in the window. Since we are 0-offset, subtracting it from
    the number of values yields the window size	*/

}

#ifdef OLD_VERSION

/*---------------------------------------------------------------------
    WINDOW_SOFTEN

	Promediates the values of every aa with its neighbours.
    This is the old version, far more slower than the former one.

    INPUT:  P_AAVAL -- the array that will receive the values.
	    P_SEQLEN -- length of the sequence whose values are to be softened.
	    WSIZE --  the size of window to use.
    OUTPUT: The values for every aa. have been promediated.

    J. R. VALVERDE
---------------------------------------------------------------------*/

void window_soften(p_aaval, p_seqlen, wsize)
float *p_aaval;
int p_seqlen;
int wsize;
{    
    register counter pos, wpos;
    int whalf_1, whalf_2;
    extern float *aa_values;

    whalf_1 = wsize >> 1;
    whalf_2 = whalf_1 + (wsize & 1);

    /* first half window aminoacids have not full left end */

    for (pos = 0; pos < whalf_1; pos++) {
	p_aaval[pos] = 0.0;
	for (wpos = 0; wpos < (pos + whalf_2); wpos++)
            p_aaval[pos] += aa_values[wpos];
	p_aaval[pos] /= (pos + whalf_2);
    }

    /* most of the sequence can completely be computed */

    for (pos = whalf_1; pos < (p_seqlen - whalf_2); pos++) {
	p_aaval[pos] = 0.0;
	for (wpos = (pos - whalf_1);wpos < (pos + whalf_2); wpos++)
            p_aaval[pos] += aa_values[wpos];
	p_aaval[pos] /= wsize;
    }

    /* last half window is again a special case */

    for (pos = p_seqlen - whalf_2;
     	pos < p_seqlen;
     	pos++)
    {
	p_aaval[pos] = 0.0;
	for (wpos = pos - whalf_1; wpos < p_seqlen; wpos++)
            p_aaval[pos] += aa_values[wpos];
	p_aaval[pos] /= (p_seqlen - (pos - whalf_1));
    }

}


/*---------------------------------------------------------------------
    MOMENT_SOFTEN

	Promediates the values of every aa with its neighbours.
    Similar to window_soften, but this routine takes into account
    the moment of every aminoacid.
    This is the old version, far more slower than the former one.

    INPUT:  P_AAVAL -- the array that will receive the values.
	    P_SEQLEN -- length of the sequence whose values are to be softened.
	    WSIZE --  the size of window to use.
    OUTPUT: The values for every aa. have been promediated.

    J. R. VALVERDE
---------------------------------------------------------------------*/

void moment_soften(p_aaval, p_seqlen, wsize, degrees)
float *p_aaval;
int p_seqlen;
int wsize;
float degrees;		/* turn in degrees per residue */
{    
    register counter pos, wpos;
    float sumcos, sumsin, delta;
    int whalf_1, whalf_2;
    extern float *aa_values;

    whalf_1 = wsize >> 1;
    whalf_2 = whalf_1 + (wsize & 1);

    delta = degrees * 2.0 * PI / 360.0;	    /* convert to radians */

    /* first half window aminoacids have not full left end */

    for (pos = 0; pos < whalf_1; pos++) {
	p_aaval[pos] = sumcos = sumsin = 0.0;
	for (wpos = 0; wpos < (pos + whalf_2); wpos++) {
            sumcos += aa_values[wpos] * cos(delta * wpos);
	    sumsin += aa_values[wpos] * sin(delta * wpos);
	}
	p_aaval[pos] = sqrt( SQR(sumcos) + SQR(sumsin) );
	p_aaval[pos] /= (pos + whalf_2);
    }

    /* most of the sequence can completely be computed */

    for (pos = whalf_1; pos < (p_seqlen - whalf_2); pos++) {
	p_aaval[pos] = sumcos = sumsin = 0.0;
	for (wpos = (pos - whalf_1);wpos < (pos + whalf_2); wpos++) {
	    sumcos += aa_values[wpos] * cos(delta * wpos);
	    sumsin += aa_values[wpos] * sin(delta * wpos);
	}
	p_aaval[pos] = sqrt( SQR(sumcos) + SQR(sumsin) );
	p_aaval[pos] /= wsize;
    }

    /* last half window is again a special case */

    for (pos = p_seqlen - whalf_2;
     	pos < p_seqlen;
     	pos++)
    {
	p_aaval[pos] = sumcos = sumsin = 0.0;
	for (wpos = pos - whalf_1; wpos < p_seqlen; wpos++) {
	    sumcos += aa_values[wpos] * cos(delta * wpos);
	    sumsin += aa_values[wpos] * sin(delta * wpos);
	}
    	p_aaval[pos] = sqrt( SQR(sumcos) + SQR(sumsin) );
	p_aaval[pos] /= (p_seqlen - (pos - whalf_1));
    }

}

#endif

/*
	Este amor que quiere ser
	acaso pronto será;
	pero ¿cuándo ha de volver
	lo que acaba de pasar?
	Hoy dista mucho de ayer.
	¡Ayer es nunca jamás!
		Antonio Machado. Humorismos, fantasías, apuntes.
*/

/*------------------------------------------------------------------*/
/*	DEBUGGING ROUTINE					    */
/*------------------------------------------------------------------*/

#if DEBUG == 3

#include "P_global.h"

char seqfile[80];
char parfile[80];
char outfile[80];
protein_t prot;
descriptor_t desc;
int i, j;

FILE *outfp;

horror(s1, s2)
char *s1, *s2;
{
    printf("ERROR. HORROR -- %s: %s\n", s1, s2);
}

extern int read_protein();
extern int read_descriptor();

main()
{
    puts("ASSIGN VALUES TO A SEQUENCE");

    /* Get sequence file */
    printf("Test sequence file name : ? ");
    scanf("%s", seqfile);
    if (read_protein(seqfile, &prot) == FAIL)
    exit(-1);

    /* Get descriptor to apply */
    printf("\nDescriptor to apply: ? ");
    scanf("%s", parfile);
    if ((i = read_descriptor(parfile, &desc)) == FAIL)
    exit(-2);


    /* Get output file */
    printf("\nFile name for output results: ? ");
    scanf("%s", outfile);
    if ((outfp = fopen(outfile, "w+")) == NULL)
    exit(-3);

    /* assign values according to parameter */
    assign_values(&prot, &desc);

    /* save results */
    fprintf(outfp, "Protein has %d aminoacids\n", prot.p_seqlen);
    for (i = 0; i < desc.d_features->l_length; ++i) {
	fprintf("\nFeature # %d:\n", i+1);
        for (j = 0; j < prot.p_seqlen; j++) {
	    fprintf(outfp,"%2.3f ", aaval[i][j]);
	    if ((j+1 % 10) == 0) fprintf(outfp, "\n");
	}
    }

    /* see you later, alligator ! */
    fclose(outfp);

    exit(EXIT_OK);
}

#endif

/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*		    END OF MODULE P_VALUES.C			*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/
