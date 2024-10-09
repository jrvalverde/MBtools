/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE	     						*
*	    PILAR_7.C						*
*								*
*	SYNOPSIS						*
*	    PILAR_7 -- Seventh PILAR program. Find regions in a	*
*	protein according to a description.			*
*								*
*	DESCRIPTION						*
*	    This program is invoked with the following command:	*
*	PILAR_7 -p protein -o output -d descriptor {-s}		*
*								*
*	    It reads the protein file specified and then it	*
*	reads the descriptor and applies all the features	*
*	contained in it to the protein. With the data found it	*
*	writes an output file with the specified name. If no	*
*	-o option is used the it writes to standard output.	*
*	The -s flag is used to specify a shorter output, i. e.	*
*	the output omits the list of values per a.a.		*
*								*
*	RETURNS							*
*	    On return it sends back a completion status to the	*
*	shell indicating the success or not of the process.	*
*								*
*	CAVEATS			       				*
*								*
*	FILES							*
*	    Uses files:						*
*		    P_DEFS.H					*
*		    P_TYPES.H					*
*		    P_EXTERN.H					*
*		    P_GLOBAL.H					*
*								*
*	CAVEATS:						*
*								*
*	NOTES							*
*	   							*
*	SEE ALSO						*
*								*
*	DESIGNED BY						*
*	    José Ramón Valverde Carrillo			*
*  								*
*	HISTORY							*
*	     8 - oct - 89 (J.R.Valverde) Written and debugged	*
*		(Happy Birthday to you...)			*
*	     1 - aug - 91 (J.R.Valverde) Ported to Ultrix.	*
*	 							*
*	COPYRIGHT						*
*	    Public domain. Not allowed to get money from it	*
*	without prior written consent of the author.		*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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

#define MIN_ARG_COUNT	5

#define DEBUG	0
/* #define DEBUG    -7 */

/*----------------------------------------------------------------------
			MODULE'S VARIABLES
----------------------------------------------------------------------*/

char prot_file[MAXPATH + 1];	    /* protein file name    */
char out_file[MAXPATH + 1];	    /* output file name	    */
char pat_file[MAXPATH + 1];	    /* pattern file name    */

FILE *out_fp;
flag out_flag = FALSE;
flag short_flag = FALSE;
protein_t prot;
descriptor_t desc;

static char program[MAXPATH + 1];
                                 

/*----------------------------------------------------------------------
			FORWARD DEFINITIONS
----------------------------------------------------------------------*/

static void save_name(char *);
static void show_usage();
static void parse_arguments(int, char *[]);
static void save_values();
static void save_regions();

/*------------------------------------------------------------------------
			MAIN PROGRAM
------------------------------------------------------------------------*/

main(argc, argv)
int argc;
char *argv[];
{
    status result;

    /* Si vales, bene est. Ego valeo */
    if ( isatty(1) ) {		    /* Is stdout a terminal ? */
	printf("PILAR_7: Apply a descriptor to a protein.\n");
	printf("\tJ.R.Valverde - %s\n", __DATE__);
	printf("=========================================\n\n");
	et_voila();
    }

    save_name("PILAR_7");
    parse_arguments(argc, argv);

    if (read_protein(prot_file, &prot) == FAIL)
	letalf("%s -- Cannot read protein file %s", program, prot_file);

    if (! read_descriptor(pat_file, &desc))
	letalf("%s -- Something is wrong with the descriptor %s",
	    program, pat_file);

    if (! out_flag)
	out_fp = stdout;
    else {
	out_fp = fopen(out_file, "w+");
	if (out_fp == NULL) {
	    letalf("%s -- Couldn't open %s to save data", program, out_file);
	}
    }

    fprintf(out_fp, "Protein: %s is %d aa long\n",
	prot.p_name, prot.p_seqlen);

    if (prot.p_seqlen < desc.d_minplen) {
	fprintf(out_fp, "%s -- The protein is too short", program);
	exit( EXIT_OK );
    }

    assign_values(&prot, &desc);
    if (short_flag == YES)
	save_values();

    result = find_regions(&prot, &desc);
    save_regions();
    if (result == FAIL) {
	fprintf(out_fp, "\nProtein can be rejected without testing rules.\n");
	fclose(out_fp);
    }
    else if (assert_rules(&prot, &desc)) {
	fprintf(out_fp, "The protein fulfills all the rules.\n");
    }
    else {
	fprintf(out_fp, "The protein does not fit the required rules.\n");
    }

    fclose(out_fp);
#ifdef MEM_DEBUG
    Mem_Display(stdout);
#endif
    exit( EXIT_OK );
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

static void parse_arguments(argc, argv)
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
		case 'p':
		case 'P':
		    BitSet(opts, BIT_0);      
		    strcpy(prot_file, argv[++i]);
		    break;
		case 's':
		case 'S':
		    short_flag = YES;
		    break;
		case 'o':
		case 'O':
		    out_flag = YES;
		    strcpy(out_file, argv[++i]);
		    break;
		case 'd':
		case 'D':
		    BitSet(opts, BIT_1);
		    strcpy(pat_file, argv[++i]);
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
    if (!BitTest(opts, BIT_0) || !BitTest(opts, BIT_1) )    {
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

static void save_name(name)
char *name;
{
    if (*name == '\0')
	strcpy(program, "PILAR_7");
    else
	strncpy(program, name, MAXPATH);
}

/*--------------------------------------------------------------------
    SHOW_USAGE
	Display a message indicating how to use the program.

    J. R. Valverde
--------------------------------------------------------------------*/

static void show_usage()
{
    fprintf(stderr, 
     	"\nUsage: %s %cp protein { %co output } %cd descriptor\n",
	program, OPTION_DELIMITER, OPTION_DELIMITER, OPTION_DELIMITER);
}

/*--------------------------------------------------------------------
    SAVE_VALUES
	Save the values per aa and window softened obtained by the
    protein.

    J. R. Valverde
--------------------------------------------------------------------*/

static void save_values()
{
    float **values;		/* table of protein aa. values		*/
    l_node_t *f_node;		/* current node in the list of feats.	*/
    feature_t *feat;		/* current feature being processed	*/
    counter fcount,		/* number of the current feature	*/
	    aacount;		/* count into the array of aa. values	*/

    extern protein_t prot;
    extern descriptor_t desc;
    extern char pat_file[MAXPATH + 1];

    fprintf(out_fp, "\nDescriptor used is: %s\n", pat_file);

    values = prot.p_aaval;
    /* Loop over every feature in the descriptor */
    for (f_node = desc.d_features->l_first, fcount = 0;
	 f_node != NULL;
	 f_node = f_node->next, fcount++) {
	/* get feature */
	feat = (feature_t *) f_node->item;
	fprintf(out_fp, "\n%d:\tValues for parameter %s are:\n",
	    (fcount + 1), feat->f_name);

	/* output the values */
	for (aacount = 0; aacount < prot.p_seqlen; aacount++) {
	    if ((aacount % 10) == 0)
		fprintf(out_fp, "\n%d:\t", aacount + 1);
	    fprintf(out_fp, "%2.3f ", values[fcount][aacount]);
     	}
	fprintf(out_fp, "\n");
    }
}

/*--------------------------------------------------------------------
    SAVE_REGIONS
	Write to the output file a list of the regions found for each
    parameter considered in the descriptor.

    J. R. Valverde
--------------------------------------------------------------------*/

static void save_regions()
{
    l_node_t	*fnode;		    /* current node in the feats. list	*/
    feature_t	*feat;		    /* current feature being considered	*/
    l_head_t	**regions;	    /* table with the lists of regions	*/
    l_node_t	*rnode;		    /* current node in the regions list	*/
    region_t	*reg;		    /* current region being output	*/
    counter	fcount,		    /* count of features		*/
		rcount;		    /* count of regions			*/

    fprintf(out_fp, "\nList of regions found:\n");
    regions = prot.p_regions;

    for (fnode = desc.d_features->l_first, fcount = 0;
	 fnode != NULL;
	 fnode = fnode->next, fcount++) {

	/* get feature */
	feat = (feature_t *) fnode->item;
	fprintf(out_fp, "\n%d:\tRegions for parameter %s are:\n",
	    (fcount + 1), feat->f_name);

	/* output regions */
	for (rnode = regions[fcount]->l_first, rcount = 0;
	     rnode != NULL;
	     rnode = rnode->next, rcount++) {

	    reg = (region_t *) rnode->item;
	    fprintf(out_fp, "\t%d.%d:\t- begin %d\t- length %d\n",
		fcount+1, rcount+1, reg->r_begin, reg->r_length);
	    }
	fprintf(out_fp, "\n");

	/* check for premature end of regions evaluation */
	if (regions[fcount]->l_length < feat->f_no_regions) {
	    fprintf(out_fp, "Too few regions\n");
	    /* No more features have been evaluated */
	    break;
	}
    }
}
