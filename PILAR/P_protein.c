/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE							*
*	    P_PROTEIN.C -- Protein manipulation module		*
*								*
*	SYNOPSIS						*
*	    init_protein()					*
*	    read_protein() -- Reads a protein from a file.	*
*	    get_protein()  -- Gets information for a protein.	*
*	    protein_fname() -- Returns name of current file.	*
*	    protein_name() -- Return name of a protein.		*
*	    protein_seq_length() -- Length of the sequence of	*
*		a given protein.				*
*	    protein_sequence() -- Sequence of a protein.	*
*	    protein_feat_values() -- Values for a given feature	*
*	    protein_feat_regions() -- Regions found for a feat.	*
*	    init_protein() -- Initialize a protein structure	*
*		(freeing all dynamic memory it may occupy)	*
*								*
*	DESCRIPTION						*
*	Protein Identification and Library Access Resource	*
*	(PILAR).						*
*								*
*	RETURNS							*
*								*
*	CAVEATS							*
*								*
*	FILES							*
*	    STDIO.H						*
*	    PORTABLE.H						*
*	    P_DEFS.H						*
*	    P_TYPES.H						*
*	    P_EXTERN.H						*
*								*
*	NOTES							*
*								*
*	SEE ALSO						*
*								*
*	DESIGNED BY						*
*	    José Ramón Valverde Carrillo			*
*								*
*	HISTORY							*
*	    21 - jul - 90 -- Built after joining several spread	*
*			     soubroutines.			*
*	     1 - aug - 91 -- Ported to Ultrix.			*
*	    30 - sep - 91 -- init_protein() moved here from	*
*			     module P_init.c			*
*	    19 - jan - 91 -- init_protein() tuned for GREATER	*
*			     performance.			*
*								*
*	COPYRIGHT						*
*	    All rights reserved for the author. This module is	*
*	hereby placed in the public domain. You CAN NOT sell or	*
*	otherwise make money from it without prior written	*
*	consent of the author.					*
*							YoEgo	*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#define DEBUG 0

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"

#include "P_list.h"
#include "P_db.h"

/* Uncomment the following line in System III, System V or Xenix 3.0 */
/* #include "P_utils.h"		/ * for rindex */

#define _P_INIT_  FALSE
#if _P_INIT_   == TRUE

/*
	Y en la cosa nunca vista
	de tus ojos me he buscado:
	en el ver con que me miras.
			Abel Martín
*/

/*-------------------------------------------------------------------
    INIT_PROTEIN
	Initialize all fields of a protein structure.

    INPUT:  PROTEIN -- The structure to be initialized.
    OUTPUT: -- The protein has been initialized.

    J. R. VALVERDE
-------------------------------------------------------------------*/

void init_protein(the_protein)
protein_t *the_protein;
{
    register counter i;
    region_t *pregion;
    list_t *regions;
    float **values;

    /* first clear index */
    the_protein->p_index.i_pos = 0L;
    the_protein->p_index.i_key[0] = '\0';

/*	RESERVED FOR FUTURE USE
*
*    for (i = 0; i < KEYWORDS_NO; i++) {
*	if (the_protein->p_index.i_key[i] != NULL)
*           free(the_protein->p_index.i_key[i]);
*	the_protein->p_index.i_key[i] = NULL;
*    }
*/
    /* clear protein name */
    the_protein->p_name[0] = '\0';

    /* clear sequence */
    for (i = 0; i < MAX_SEQ_LEN; i++) {
	the_protein->p_seq[i] = '\0';
    }
    the_protein->p_seqlen = 0;

    /* free arrays with the values for the previous protein */
    if (the_protein->p_aaval != NULL) {
	/* free each feature values */
	values = the_protein->p_aaval;
	for (i = 0; i < the_protein->p_nofeats; ++i)
	    free(values[i]);
	/* free global feat. values array */
        free(values);
	the_protein->p_aaval = NULL;
    }

    /* clear regions list (if it exists) */
    if (the_protein->p_regions != NULL) {
	/* first free each list independently */
	regions = the_protein->p_regions;
	for (i = 0; i < the_protein->p_nofeats; ++i) {
	    if (regions[i] == NULL)
		continue;
    	    /* delete every node (and its data !) in the list */
	    while (list_length(regions[i]) != 0) {
		pregion = (region_t *) list_del_first(regions[i]);
		free( pregion );
	    }            
	    /* dispose of the list */
	    free( regions[i] );
    	    regions[i] = NULL;
	}
	/* and then free the array of list headers */
	free(regions);
	the_protein->p_regions = NULL;
    }

    /* reset counter of measured features in this protein */
    the_protein->p_nofeats = 0;

}

#endif

/*---------------------------------------------------------------------
       	READ_PROTEIN

	Reads a protein from the specified file.

	INPUT:  P_FILNAM: the name of the file containing the protein info.
	OUTPUT: THE_PROT filled with the info found in the file.
	RETURNS: SUCCES if no problem, FAIL otherwise.

	J. R. Valverde
----------------------------------------------------------------------*/

static char p_filnam[MAXPATH + 1] = {'\0'};

status read_protein(p_filename, the_prot)
char *p_filename;
protein_t *the_prot;
{
    db_t p_file;
    int  result;

    /* Open specified file containing the information */
    if ((p_file = db_open(p_filename, DB_UNKNOWN)) == NULL) {
	horror("Read Protein", "Cannot open specified protein file");
	return FALSE;
    }

    result = db_next_protein(p_file, the_prot);

    db_close(p_file);

    /* save protein file name */
    strcpy(p_filnam, p_filename);

    return result;
}

/*--------------------------------------------------------------

    GET_PROTEIN

    Asks the user for the file containing the protein sequence
    and reads it.

    INPUT: THE_PROTEIN  structure to fill with protein info.
    OUTPUT: SUCCESS if read's been OK. FAIL otherwise.

    J. R. Valverde
---------------------------------------------------------------*/

status get_protein(the_protein)
protein_t *the_protein;
{
    char *p;
    extern char p_filnam[MAXPATH + 1];

    printf("\nName of file containing protein sequence: ? ");
    scanf("%s", p_filnam);

    if ( (p = rindex(p_filnam, ENDDIRCHAR)) == NULL)
	p = p_filnam;
    if ( rindex(p_filnam, FILEXTCHAR) == NULL)
	strcat(p_filnam, P_EXT);

    return ( read_protein(p_filnam, the_protein) );
}                     

/*-----------------------------------------------------------
    PROTEIN_FNAME

    Returns the name of the last protein file read.

    INPUT:  P_NAME -- String to store the name.
    OUTPUT: SUCCESS if so, FAIL if there is no current
	protein name.

    CAUTION: This function is SYSTEM DEPENDENT. Supposes
	that p_name is long enough.

    J. R. Valverde
------------------------------------------------------------*/

status protein_fname(p_name)
char *p_name;
{
    char *nam;
    extern char p_filnam[MAXPATH + 1];

    if (p_filnam[0] == '\0')
	return FAIL;

    nam = rindex(p_filnam, ENDDIRCHAR);
    if (nam == NULL)
	nam = p_filnam;
    else
	nam++;

    strcpy(p_name, nam);
}

/*--------------------------------------------------------------------
    PROTEIN_NAME
	Return the name of the specified protein read.

    INPUT:  PROT   -- protein.
    OUTPUT: A pointer to the name of the protein.

    CAUTION: This function returns a pointer to the content of
	the name field in the given protein. Should you want to
	keep the name after re-filling this field, you must
	save it elsewhere with STRCPY.

    J. R. VALVERDE

--------------------------------------------------------------------*/

char *protein_name(prot)
protein_t *prot;
{
    return prot->p_name;
}

/*--------------------------------------------------------------------
    PROTEIN_SEQ_LENGTH
	Return length of sequence of a given protein

    INPUT:  PROT -- the protein
    OUTPUT: length of the sequence

    CAUTION: The length will vary if you re-fill the protein with
	    new data.

    J. R. VALVERDE
--------------------------------------------------------------------*/

int protein_seq_length( prot )
protein_t *prot;
{
    return prot->p_seqlen;
}

/*--------------------------------------------------------------------
    PROTEIN_SEQUENCE
	Return a pointer to the sequence of a protein.

    INPUT:  PROT    -- the protein
    OUTPUT: A pointer to the sequence of the protein

    CAUTION:	The sequence will vary if you re-fill the protein
	data by re-reading it. Hence, you must save it elsewhere
	if you'll need this data thereafter.

    J. R. VALVERDE

--------------------------------------------------------------------*/

char *protein_sequence( prot )
protein_t *prot;
{
    return prot->p_seq;
}

/*--------------------------------------------------------------------
    PROTEIN_FEAT_VALUES
	Return an array with the values obtained for the given
    parameter in the specified protein.

    INPUT:  PROT    -- the protein
	    FEAT    -- the number of the feature
    OUTPUT: an array of values.

    CAUTION: The values may vary if you re-read the protein and re-assign
	them. You should hence save them if you'll need thereafter.
	YOU SHOULD HAVE CALLED ASSIGN_VALUES BEFORE CALLING THIS FUNCTION

    J. R. VALVERDE

--------------------------------------------------------------------*/

float *protein_feat_values( prot, f_no )
protein_t *prot;
int f_no;
{
    return prot->p_aaval[f_no];
}

/*--------------------------------------------------------------------
    PROTEIN_FEAT_REGIONS
	Return a list of significative regions found for the
    specified protein and a given parameter.

    INPUT:  PROT    -- the protein
	    F_NO    -- the parameter number

    OUTPUT: a list of regions.

    CAUTION:	You must have called ASSIGN_VALUES and FIND_REGIONS before
	calling this function.

    J. R. VALVERDE
--------------------------------------------------------------------*/

list_t protein_feat_regions( prot, f_no )
protein_t *prot;
int f_no;
{
    return prot->p_regions[f_no];
}

/*-------------------------------------------------------------------
    INIT_PROTEIN
	Initialize all fields of a protein structure.

    INPUT:  PROTEIN -- The structure to be initialized.
    OUTPUT: -- The protein has been initialized.

    J. R. VALVERDE
-------------------------------------------------------------------*/

void init_protein(the_protein)
protein_t *the_protein;
{
    register counter i;
    region_t *pregion;
    list_t *regions;
    float **values;

    /* first clear index */
    the_protein->p_index.i_pos = 0L;
/*
    for (i = 0; i < KEYWORDS_NO; i++) {
	if (the_protein->p_index.key[i] != NULL)
            free(the_protein->p_index.key[i]);
	the_protein->p_index.key[i] = NULL;
    }
*/
    /* clear protein name */
    the_protein->p_name[0] = '\0';

    /* clear sequence: */
    /* NOTE: profiling has shown that a database run may spend about
    44% (or more) of processing time (!) in the three lines below: * /
	for (i = 0; i < MAX_SEQ_LEN; i++) {
	    the_protein->p_seq[i] = '\0';
	}
    / * Thus, though it is less secure, I've decided to change to this
    other mechanism. It means that database reading routines MUST be
    sure that they add a '\0' char after the sequence. */

    the_protein->p_seq[0] = '\0';

    the_protein->p_seqlen = 0;

    /* free arrays with the values for the previous protein */
    if (the_protein->p_aaval != NULL) {
	/* free each feature values */
	values = the_protein->p_aaval;
	for (i = 0; i < the_protein->p_nofeats; ++i)
	    free(values[i]);
	/* free global feat. values array */
        free(values);
	the_protein->p_aaval = NULL;
    }

    /* clear regions list (if it exists) */
    if (the_protein->p_regions != NULL) {
	/* first free each list independently */
	regions = the_protein->p_regions;
	for (i = 0; i < the_protein->p_nofeats; ++i) {
	    if (regions[i] == NULL)
		continue;
    	    /* delete every node (and its data !) in the list */
#ifdef V_3
    	    while (regions[i]->l_length) {
        	pregion = (region_t *) delete_node(regions[i]);
        	free(pregion);
        	pregion = NULL;
    	    }
    	    /* free regions list header */
    	    free(regions[i]);
#else
	    while (list_length(regions[i]) != 0) {
		pregion = (region_t *) list_del_first(regions[i]);
		free( pregion );
	    }
	    /* dispose of the list */
	    free( regions[i] );
#endif
    	    regions[i] = NULL;
	}
	/* and then free the array of list headers */
	free(regions);
	the_protein->p_regions = NULL;
    }

    /* reset counter of measured features in this protein */
    the_protein->p_nofeats = 0;

}


/*--------------------------------------------------------------------
    Produces and does not own
    Produces and does not appropiate
    Directs without looking for dominion
    Once achieved does not claim
    Not claiming does not lose.

			    Lao Tse. Tao Te King.
--------------------------------------------------------------------*/


/*------------------------------------------------------------------*/
/*		    D E B U G G I N G 				    */
/*------------------------------------------------------------------*/

#if DEBUG == 2

#include "P_global.h"

protein_t prot;
char prfile[80];
char infile[80];
char outfile[80];
FILE *fp1;
FILE *fp2;
int no_aa, count;                       

main()
{
    static int i, j, end;

    puts("TESTING MODULE PROTEIN");
    puts("======================\n");
                           
    puts("READ_PROTEIN");
    printf("Name of protein file: ? ");
    scanf("%s", prfile);
    if (read_protein(prfile, &prot) != SUCCESS)
	exit(-1);

}
#endif
