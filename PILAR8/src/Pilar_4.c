/*                                                         
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*									*
*	MODULE:								*
*	    PILAR_4.C							*
*									*
*	SYNOPSIS:							*
*	    Body of the main program for PILAR_4 (4-th program of	*
*	the Protein Identification and Library Access Resource).	*
*									*
*	DESCRIPTION:							*
*	    The program is invoked by using the command line		*
*									*
*	PILAR_4 -l library { -i index } {-o output} {-v}		*
*									*
*	If an index file is not specified then library access is	*
*	sequential. This program extracts names of proteins from a	*
*	library, either using an index file, or accessing it in a	*
*	sequential form. If -v is expressed, more information is	*
*	printed about each protein.					*
*									*
*	RETURNS:							*
*	    On return the program has checked the library or		*
*	the proteins specified in the index file if success, otherwise	*
*	it returns an error code.		   			*
*	    								*
*	CAVEATS:							*
*									*
*	FILES:				   				*
*	    This program uses the following files and modules:		*
*		STDIO.H							*
*		ERRNO.H							*
*		PORTABLE.H						*
*		P_DEFS.H						*
*		P_TYPES.H						*
*		P_EXTERN.H						*
*		P_GLOBAL.H						*
*									*
*	NOTES:								*
*	    To be enhanced for searching an index file by keywords.	*
*									*
*	SEE ALSO:							*
*									*
*	DESIGNED BY:							*
*	    José Ramón Valverde Carrillo.				*
*									*
*	HISTORY:							*
*	    13 - aug - 1989 (J.R. Valverde - Initiation)		*
*	    15 - aug - 1989 (J.R. Valverde - Depuration)		*
*	    17 - aug - 1989 (J.R. Valverde - Purged exceding code)	*
*									*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/
                                                                         
/*----------------------------------------------------------------------
			INCLUDE FILES
----------------------------------------------------------------------*/

#include <stdio.h>
#include <errno.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_global.h"

#include "P_db.h"
#include "P_init.h"
#include "P_utils.h"

#include "mshell.h"

/*----------------------------------------------------------------------
			PRIVATE DEFINITIONS
----------------------------------------------------------------------*/

#define MIN_ARG_COUNT	3

#define DEBUG	0
/* #define DEBUG    -4 */

/*----------------------------------------------------------------------
			MODULE'S VARIABLES
----------------------------------------------------------------------*/

char lib_file[MAXPATH + 1];	    /* library file name    */
char ind_file[MAXPATH + 1];	    /* index file name	*/
char out_file[MAXPATH + 1];	    /* output file name	*/
char pat_file[MAXPATH + 1];	    /* pattern file name    */

FILE *out_fp;

flag ind_flag = NO;	    /* use an index file ?  */
flag out_flag = NO;	    /* use an output file ? */
flag ver_flag = NO;	    /* verbose mode? */

static char program[MAXPATH + 1];

static db_t db;


/*----------------------------------------------------------------------
			FORWARD DEFINITIONS
----------------------------------------------------------------------*/

void save_name(char *);
void show_usage();
void parse_arguments(int, char *[]);
int  new_protein(protein_t *);
void save_protein(protein_t *);

/*------------------------------------------------------------------------
			MAIN PROGRAM
------------------------------------------------------------------------*/

main(argc, argv)
int argc;
char *argv[];
{
    static protein_t p;	    /* This initializes p to zeroes */
    static descriptor_t d;

    if ( isatty(1) ) {
	printf("PILAR_4: Extract protein name from a library.\n\t%s\n",
		__DATE__);
	et_voila();
    }
    save_name("PILAR_4");
    parse_arguments(argc, argv);

    /* open files */
    if ((db = db_open(lib_file, DB_UNKNOWN)) == NULL) {
	letalf("%s -- Cannot open database %s", program, lib_file);
    }
    if (ind_flag == YES)
	/* open index file for reading */
	if (! idx_open(ind_file)) {
	    letalf("%s -- Cannot open index file %s", program, ind_file);
     	}
    if (out_flag == YES) {
	if ((out_fp = fopen(out_file, "w+")) == NULL) {
	    letalf("%s -- Cannot open output file %s", program, out_file);
	}
    }
    else
	out_fp = stdout;

    /* Process database */
    while ( new_protein( &p ) ) {
	save_protein( &p );
	init_protein( &p );
    }

    /* Close all open files */
    db_close(db);
    idx_close();
    if (out_fp != stderr)
	fclose(out_fp);

    Mem_Display(stdout);

    exit(EXIT_OK);

}                                                             

/*----------------------------------------------------------------------
    PARSE_ARGUMENTS
	Parse command line arguments. Checks arguments passed
    to the program and sets the global variables LIB_FILE, OUT_FILE,
    and optionally IND_FILE to the appropriate values.
                         
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
		case 'v':
		case 'V':
		    ver_flag = YES;
		    break;
		case 'o':
		case 'O':
		    out_flag = YES;
		    strcpy(out_file, argv[++i]);
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
    if (!BitTest(opts, BIT_0) )	{
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
	strcpy(program, "PILAR_4");
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
    fprintf(stderr, 
     	"\nUsage: %s %cl library { %ci index } { %co output } { %cv }\n",
	program, OPTION_DELIMITER, OPTION_DELIMITER, OPTION_DELIMITER,
		OPTION_DELIMITER);
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

    static int i = 0;

    i++;
    if (ind_flag == YES) {
	if (get_index(&next)) {
	    res = db_seek_protein(db, &next, protein);
	    if (res == FALSE) {
		/***** NOTE: SYSTEM DEPENDENT CODE *****/
#ifdef __VAXC__
		fprintf(stderr, "ERROR: errno = %d, vaxc$errno = %d\n%s\n",
		    errno, vaxc$errno, strerror(errno));
#else
		fprintf(stderr, "ERROR: errno = %d\n%s\n", errno, strerror(errno));
#endif
		fprintf(stderr, "\ni = %d\nprotein: %s\nindex = %ld\n", i,
			protein->p_name, protein->p_index.i_pos);
	    }
	    return res;
	}
	else
	    return FALSE;
    }
    else {
	res = db_next_protein(db, protein);
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

    fprintf(out_fp, "\nProtein: %s", protein->p_name);
    if (ver_flag) {	    /* Verbose mode */
	fprintf(out_fp, "\nHas %u aminoacids", protein->p_seqlen);
	fprintf(out_fp, "\nIt's sequence is:");
	for (i = 0; i < protein->p_seqlen; ++i) {
	    if ( (i % 80) == 0 )
		fprintf(out_fp, "\n");
	    fprintf(out_fp, "%c", protein->p_seq[i]);
	}
	fprintf(out_fp, "\n");
    }
}

/*-------------------------------------------------------------------
|								    |
|			That's all folks !			    |
|								    |
-------------------------------------------------------------------*/

