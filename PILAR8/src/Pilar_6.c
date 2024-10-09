/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*									*
*	MODULE:								*
*	    PILAR_6.C							*
*									*
*	SYNOPSIS:							*
*	    Body of the main program for PILAR_6 (6-th program of	*
*	the Protein Identification and Library Access Resource).	*
*									*
*	DESCRIPTION:							*
*	    The program is invoked by using the command line		*
*									*
*	PILAR_6 -l library {-i index } -o output -d descriptor		*
*									*
*	If an index file is not specified then library access is	*
*	sequential.							*
*									*
*	RETURNS:							*
*	    On return the program has checked the library or		*
*	the proteins specified in the index file if success, otherwise	*
*	it returns an error code.		   			*
*	    								*
*	CAVEATS:							*
*									*
*	FILES:				   				*
*									*
*	NOTES:								*
*									*
*	SEE ALSO:							*
*									*
*	DESIGNED BY:							*
*	    José Ramón Valverde Carrillo.				*
*									*
*	HISTORY:							*
*	    20 - aug - 1989 (J.R.Valverde - Initiation)			*
*									*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

/*
    ... so exceedingly beautiful that the Sun himself. although he saw
    her very often, was enchanted every time she came out into the
    sunshine.

		The brothers Grimm. The Frog Prince.
*/

/*----------------------------------------------------------------------
			INCLUDE FILES
----------------------------------------------------------------------*/

#include <stdio.h>
#include <errno.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_extern.h"
#include "P_global.h"

/*----------------------------------------------------------------------
			PRIVATE DEFINITIONS
----------------------------------------------------------------------*/

#define MIN_ARG_COUNT	7
#define DEBUG		0
/* #define DEBUG	    -6 */

/*----------------------------------------------------------------------
			MODULE'S VARIABLES
----------------------------------------------------------------------*/

char lib_file[MAXPATH + 1];	    /* library file name    */
char ind_file[MAXPATH + 1];	    /* index file name	*/
char out_file[MAXPATH + 1];	    /* output file name	*/
char dsc_file[MAXPATH + 1];	    /* pattern file name    */

FILE *out_fp;

flag ind_flag = NO;	    /* use an index file ?  */

static char program[MAXPATH + 1];
                                 
static db_t database;	    /* database in search */

/*----------------------------------------------------------------------
			FORWARD DEFINITIONS
----------------------------------------------------------------------*/

void save_name(char *);
void show_usage();
void parse_arguments(int, char *[]);
int  new_protein(protein_t *);
int check_protein(protein_t *, descriptor_t *);
void save_protein(protein_t *);

/*------------------------------------------------------------------------
			MAIN PROGRAM
------------------------------------------------------------------------*/

main(argc, argv)
int argc;
char *argv[];
{
    static int i;
    static protein_t p;	    /* This initializes p to zeroes */
    static descriptor_t d;

    save_name("PILAR_6");
    if ( isatty(1) ) {	    /* Is stdout a terminal? */
	printf("Homo sum, nihil humanum a me alienum puto...\n");
	et_voila();
	printf("PILAR_6 - J. R. Valverde - %s\n", __DATE__);
    }

    parse_arguments(argc, argv);

    /* Read descriptor */
    if (! read_descriptor(dsc_file, &d)) {
	letalf("%s -- When reading descriptor file %s", program, dsc_file);
    }

    /* open files */
    if ( (database = db_open(lib_file, DB_UNKNOWN)) == NULL ) {
	letalf("%s -- Cannot open database %s", program, lib_file);
    }
    if (ind_flag == YES)
	/* open index file for reading */
	if (! idx_open(ind_file)) {
	    letalf("%s -- Cannot open index file %s", program, ind_file);
     	}
    if ((out_fp = fopen(out_file, "w+")) == NULL) {
	letalf("%s -- Cannot open output file %s", program, out_file);
    }

    /* Process database */
    i = 0;
    init_protein( &p );
    while ( new_protein( &p ) ) {
	i++;
      	if (check_protein( &p, &d ) == TRUE)
	    save_protein( &p );
	init_protein( &p );
    }

    /* Close all open files */
    db_close( database );
    idx_close();
    fclose(out_fp);

#ifdef MEM_DEBUG
/*    Mem_Display(stdout); */
#endif

   exit(EXIT_OK);

}                                                             

/*----------------------------------------------------------------------
    PARSE_ARGUMENTS
	Parse command line arguments. Checks arguments passed
    to the program and sets the global variables LIB_FILE, OUT_FILE,
    DSC_FILE and optionally IND_FILE to the appropriate values.
                         
    INPUT:  ARGC -- Count of command line arguments
	ARGV -- Array containing the arguments.
    OUTPUT: If a syntax error is found outputs a warning message
	telling the correct usage and aborts the program.

    J. R. VALVERDE
-----------------------------------------------------------------------*/

void parse_arguments(argc, argv)
int argc;
char *argv[];  
{
    counter i;
    flag    opts = 0;

    if (argc < MIN_ARG_COUNT) {
	show_usage();
	exit(BAD_ARGS);
    }

    for (i = 1; i < argc; i++) {
	if (argv[i][0] == OPTION_DELIMITER) {
	    /* it's a flag */
	    switch (argv[i][1]) {
		case 'l':
		case 'L':
		    BitSet(opts, BIT_0);      
		    strcpy(lib_file, argv[++i]);
		    break;
		case 'i':
		case 'I':
		    ind_flag = YES;
		    strcpy(ind_file, argv[++i]);
		    break;
		case 'o':
		case 'O':
		    BitSet(opts, BIT_1);
		    strcpy(out_file, argv[++i]);
		    break;
		case 'd':
		case 'D':
		    BitSet(opts, BIT_2);
		    strcpy(dsc_file, argv[++i]);
		    break;
		default:
		    show_usage();
		    exit(BAD_ARGS);
		    break;
	    }
	}
	else {
	    show_usage();
	    exit(BAD_ARGS);
	}
    }
    if (!BitTest(opts, BIT_0) ||
	!BitTest(opts, BIT_1) ||
	!BitTest(opts, BIT_2)  )
    {
	show_usage();
	exit(BAD_ARGS);
    }
}

/*--------------------------------------------------------------------
    SAVE_NAME
	Saves current program name. If there is no name, uses PILAR by
    default.

    INPUT:  NAME -- Current program name

    J. R. Valverde
--------------------------------------------------------------------*/

void save_name(name)
char *name;
{
    if (*name == '\0')
	strcpy(program, "PILAR_6");
    else
	strncpy(program, name, MAXPATH);
}

/*--------------------------------------------------------------------
    SHOW_USAGE
	Display a message indicating how to use the program.

    J. R. Valverde
--------------------------------------------------------------------*/

void show_usage()
{
    fprintf(stderr, "%s:\n", program);
    fprintf(stderr,
"Scan for a description of protein structure in a sequence database\n");
    fprintf(stderr, 
     	"\nUsage: %s %cl library [ %ci index ] %co output %cd descriptor\n",
	program, OPTION_DELIMITER, OPTION_DELIMITER,
	OPTION_DELIMITER, OPTION_DELIMITER);
}

/*--------------------------------------------------------------------
    NEW_PROTEIN
	Get next protein to analyze.
    If we are using an index file, get index and look for the correct
    sequence. Else access sequentially the database.

    INPUT:  PROTEIN -- Space to store the info read.

    J. R. Valverde
---------------------------------------------------------------------*/
         
int new_protein( protein )
protein_t *protein;
{
    index_t next;
    int res;
    extern flag ind_flag;

    if (ind_flag == YES) {
	/* If we are using an index file */
	res = get_index(&next);
	if (res != EOF && res != 0) {
	    res = db_seek_protein(database, &next, protein);
	    if (res == FALSE) {
#ifdef __VAXC__
		printf("ERROR: errno = %d, vaxc$errno = %d (%s)\n",
		    errno, errno, strerror(errno));
#else
	printf("ERROR: errno = %d (%s)\n", errno, strerror(errno));
#endif
		return FALSE;
	    }
	    return res;
	}
	else
	    return FALSE;
    }
    else {
	res = db_next_protein(database, protein);
	return res;
    }
}

/*--------------------------------------------------------------------
    SAVE_PROTEIN
	Save the relevant info about this protein.

    INPUT:  PROTEIN -- The protein whose information is to be saved

    J. R. Valverde
--------------------------------------------------------------------*/

void save_protein(protein)
protein_t *protein;
{
    counter i;

    fprintf(out_fp, "%ld\t%s", protein->p_index.i_pos, protein->p_index.i_key);
#if (DEBUG == -6)
    fprintf(out_fp, " : Protein: %s", protein->p_name);
    fprintf(out_fp, "\nHas %u aminoacids", protein->p_seqlen);
    fprintf(out_fp, "\nIt's sequence is:");
    for (i = 0; i < protein->p_seqlen; ++i) {
	if ( (i % 80) == 0 )
	    fprintf(out_fp, "\n");
	fprintf(out_fp, "%c", protein->p_seq[i]);
    }
#endif
    fprintf(out_fp, "\n");
}
                                 
/*-------------------------------------------------------------------
    CHECK_PROTEIN
	Check that the current protein fulfills all imposed
    restrictions.

    INPUT:  P -- The protein to check.
	    D -- The descriptor that it has to match
    OUTPUT: TRUE -- Iff the protein passes ALL the tests
	    FALSE -- If the protein should be rejected.

    J. R. VALVERDE
--------------------------------------------------------------------*/

boolean check_protein(p, d)
protein_t *p;
descriptor_t *d;
{
    status result;
    feature_t *feat;
    l_node_t *node;
    region_t *reg;

    /* first discard too short proteins */
    if (p->p_seqlen < d->d_minplen)
	return FALSE;

    /* check for the existence of enough regions */
    assign_values(p, d);
    if (! find_regions(p, d))
	return FALSE;
    /* and then check out the rules */
    if ( assert_rules(p, d) )
	return TRUE;	/* It (hopefully) passed all the checks */
    else
	return FALSE;
}

/*
    hi o kesu ya
    kokoro gake nasu
    tsuki no mae

    Dull light,
    Cliff of the soul
    facing the moon.

		Japanese haiku.
*/

