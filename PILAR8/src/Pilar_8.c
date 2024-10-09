/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE	     						*
*	    PILAR_8.C						*
*								*
*	SYNOPSIS						*
*	    PILAR_8 -- Eigth PILAR program. Find regions in a	*
*	protein according to a description. SAME AS PILAR_7.	*
*	THIS IS AN VT240 VERSION WITH ReGis GRAPHICS SUPPORT.	*
*								*
*	DESCRIPTION						*
*	    This program is invoked with the following command:	*
*	PILAR_8 -p protein -o output -d descriptor -g -s -n	*
*								*
*	    It reads the protein file specified and then it	*
*	reads the descriptor and applies all the features	*
*	contained in it to the protein. With the data found it	*
*	writes an output file with the specified name. If no	*
*	-o option is used then it writes to standard output.	*
*	It also draws several graphics in the screen showing	*
*	the corresponding information if -g is used.		*
*	If the -s flag is used then the report will be a short	*
*	one, i. e. without aminoacid values.			*
*	If the -n flag is given then no report will be shown.	*
*								*
*	RETURNS							*
*	    On return it sends back a completion status to the	*
*	shell indicating the success or not of the process.	*
*								*
*	CAVEATS			       				*
*	    I have used this program only on a VT240 terminal.	*
*	It uses ReGis graphics mode, so it may run on other	*
*	ReGis supporting devices, BUT I'M NOT SURE AT ALL.	*
*								*
*	FILES							*
*	    Uses files:						*
*		    P_DEFS.H					*
*		    P_TYPES.H					*
*		    P_EXTERN.H					*
*		    P_GLOBAL.H					*
*								*
*		    CURSES.H					*
*		    P_GRAPHIC.X					*
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
*	    25 - jun - 90 (J.R.Valverde) Added -s flag and	*
*		enhanced usage()				*
*	     5 - jan - 91 (J.R.Valverde) Added -n flag.		*
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
#include <curses.h>
#include <math.h>
#include <float.h>


#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_extern.h"
#include "P_global.h"
#include "P_graph.x"

/*----------------------------------------------------------------------
			PRIVATE DEFINITIONS
----------------------------------------------------------------------*/

#define MIN_ARG_COUNT	5

#define DEBUG	0
/* #define DEBUG    -8 */

/*----------------------------------------------------------------------
			MODULE'S VARIABLES
----------------------------------------------------------------------*/

char prot_file[MAXPATH + 1];	    /* protein file name    */
char out_file[MAXPATH + 1];	    /* output file name	    */
char pat_file[MAXPATH + 1];	    /* pattern file name    */

FILE *out_fp;
flag out_flag = NO;
flag grf_flag = NO;
flag short_flag = NO;
flag silent_flag = NO;
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
static void draw_data();
static void init_map();
static void draw_axis();
static void draw_values();
static void draw_regions(int);

/*
	A quoi pourroit-il servir de prier d'amour une grande dame
	qui ne manqueroit pas d'en rire à nos dépens?

		Le Roman de Renart.
*/

/*------------------------------------------------------------------------
			MAIN PROGRAM
------------------------------------------------------------------------*/

main(argc, argv)
int argc;
char *argv[];
{
    status result;

    /* Si vales, bene est. Ego valeo */
    if ( isatty(1) ) {		    /* Is stdout a terminal? */
	printf("PILAR_8: Apply a descriptor to a protein.\n");
	printf("    J.R.Valverde - (VT 240) - %s\n", __DATE__);
	printf("=========================================\n\n");
	et_voila();
    }

    save_name("PILAR_8");
    parse_arguments(argc, argv);

    if (! read_protein(prot_file, &prot))
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

    result = find_regions(&prot, &desc);

    if (! silent_flag) {
	if (! short_flag)
	    save_values();
	save_regions();
    }

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

    if (out_flag)
	fclose(out_fp);
    if (grf_flag == YES) {
	printf("\nPress any key to see the graphs...\n");
	draw_data();
    }

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
		case 'o':
		case 'O':
		    out_flag = YES;
		    strcpy(out_file, argv[++i]);
		    break;
		case 'g':
		case 'G':
		    grf_flag = YES;
		    break;
		case 's':
		case 'S':
		    short_flag = YES;
		    break;
		case 'n':
		case 'N':
		    silent_flag = YES;
		    short_flag = YES;
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
	strcpy(program, "PILAR_8");
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
    "\nUsage: %s %cp protein %cd descriptor {%co output} {%cg} {%cs} {%cn}\n",
	program, OPTION_DELIMITER, OPTION_DELIMITER, OPTION_DELIMITER,
	OPTION_DELIMITER, OPTION_DELIMITER, OPTION_DELIMITER);
    fprintf(stderr,
	"\tg\tUse this flag to get a graphic output\n%s%s",
	"\ts\tUse this flag to get a short report\n",
	"\tn\tUse this flag to inhibit reports.\n");
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
    region_t	*rgn;		    /* current region being output	*/
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
	    int rgnlen;

	    rgn = (region_t *) rnode->item;
	    rgnlen = rgn->r_begin + rgn->r_length;
	    fprintf(out_fp, "\t%d.%d:\t- begin %d\t- end %d\t- length %d\n",
		fcount+1, rcount+1, rgn->r_begin, rgnlen, rgn->r_length);
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

/*--------------------------------------------------------------------
    DRAW_DATA
	This subroutine displays all he data available about the
    protein in a graphical form to the screen using ReGis graphics
    mode.

    J. R. Valverde
--------------------------------------------------------------------*/

#define BLACK	0
#define BLUE	1
#define RED	2
#define GREEN	3

static int  x_min,
	    x_max,
	    x_len,
	    y_min,
	    y_max,
	    y_len;

static float x_step;

static void draw_data()
{
    l_node_t	*fnode;
    feature_t	*feat;
    counter	fcount;
    char 	ch;
    float	true_threshold;

    ch = getchar();
    initscr();
    clear();
    init_map();

    for (fnode = desc.d_features->l_first, fcount = 0;
	 fnode != NULL;
	 fnode = fnode->next, fcount++) {
	feat = (feature_t *) fnode->item;
	initscr();
	clear(); refresh();
	printf("\tPlot for feature %s\t%s\n",
	    feat->f_name, (feat->f_inverted)? "(inverted)" : "");

	/* Threshold was normalized. Un-normalize it */
	true_threshold = feat->f_threshold;
	if (feat->f_inverted) true_threshold = 1.0 - true_threshold;
	true_threshold =
	    (true_threshold * (feat->f_maxvalue - feat->f_minvalue)) +
		feat->f_minvalue;

	printf("\n\tMax. val. = %+2.3f, min. val. = %+2.3f, threshold = %+2.3f\n",
	    feat->f_maxvalue, feat->f_minvalue, true_threshold);

	printf("\n\tWindow = %d, region range = [%d-%d]",
	    feat->f_window, feat->f_minrglen, feat->f_maxrglen);
	if (feat->f_moment >= 0.0)
	    printf(", moment = %+3.3f", feat->f_moment);

	enterReGis(1);
	init_map();
	draw_axis();
	draw_values( fcount, feat->f_threshold );
	draw_regions( fcount );
	ch = getchar();
	EraseScreen();
	exitReGis();
    }
    clear(); refresh();
}

/*--------------------------------------------------------------------
	INIT_MAP
	    Initialize colour map to the desired colours.

	J. R. Valverde
--------------------------------------------------------------------*/

static void init_map()
{
#ifdef FULL_GRAPHICS
    enterReGis(1);
    SelMapHLS(0, 0, 0, 0);	    /* Map location 0 is black	*/
    SelMapHLS(1, 0, 50, 100);	    /* Map location 1 is blue	*/
    SelMapHLS(2, 120, 50, 100);	    /* Map location 2 is red	*/
    SelMapHLS(3, 240, 50, 100);	    /* Map location 3 is green	*/
    exitReGis();
#else
    return;
#endif
}

/*--------------------------------------------------------------------
	DRAW_AXIS
	    This procedure draws the set of axis to use in the display

	J. R. Valverde
--------------------------------------------------------------------*/

static float itof(x) int x; { return (x); }

static void draw_axis()
{
    x_min = 25;
    x_max = 775;
    y_min = 120;
    y_max = 360;
    x_len = 750;
    y_len = 240;

    x_step = itof(x_len) / itof(prot.p_seqlen);

    SelFgIntensity( BLACK );
    MoveTo(x_min, y_min);
    SelFgIntensity( GREEN );
    LineTo(x_min, y_max);
    LineTo(x_max, y_max);
    LineTo(x_max, y_min);
    LineTo(x_min, y_min);
    SelFgIntensity( BLACK );
}

/*--------------------------------------------------------------------
	DRAW_VALUES
	    This functions displays the values for the specified
	parameter.

	INPUT:	FCOUNT -- The feature whose values are to be displayed
		THRESHOLD -- The feature's threshold

	J. R. Valverde
--------------------------------------------------------------------*/

static void draw_values(fcount, threshold)
counter fcount;
float threshold;
{
    float	val;
    float	**all_values, *values;
    counter	aacount;
    float	x_coord, y_coord;
    int		x, y;

    all_values = prot.p_aaval;
    values = all_values[fcount];

    /* threshold line */
    SelFgIntensity( BLUE );
    y_coord = itof(y_max) - (threshold * itof(y_len));
    y = ceil(y_coord);
    MoveTo(x_min, y);
    LineTo(x_max, y);
    SelFgIntensity( BLACK );

    /* First point */
    MoveTo(x_min, y_max );  
    SelFgIntensity( GREEN );
    val = values[0];
    x_coord = itof(x_min);
    y_coord = itof(y_max) - (val * itof(y_len));
    x = floor(x_coord);
    y = ceil(y_coord);
    LineTo(x, y);

    for (aacount = 1; aacount < prot.p_seqlen; aacount++) {
	val = values[aacount];	    /* values between 0 and 1 */
	x_coord += x_step;
	y_coord = itof(y_max);
	y_coord -= (val * itof(y_len));	/* scale to 0 .. 100 */
	x = floor(x_coord);
	y = ceil(y_coord);
	if (val < threshold)
	    SelFgIntensity( GREEN );
	else
	    SelFgIntensity( RED );
	LineTo(x, y);
    }
    SelFgIntensity( BLACK );
}

/*--------------------------------------------------------------------
	DRAW_REGIONS
	    Display the regions found for this feature (if any) under
	the plot of values.

	INPUT:	FCOUNT	-- The feature whose regions are to be shown.

	J. R. Valverde
--------------------------------------------------------------------*/

extern region_t ***RegMap;

static void draw_regions(fcount)
counter fcount;
{
    l_head_t	**all_regions, *regions;
    l_node_t	*rnode;
    region_t	*rgn;
    float	x_coord;
    int		x, y;

    all_regions = prot.p_regions;
    if (all_regions[fcount] == NULL)
	return;
    regions = all_regions[fcount];

    SelShadeYRef(440);
    y = 460;

    for (rnode = regions->l_first;
	 rnode != NULL;
	 rnode = rnode->next) {
	rgn = (region_t *) rnode->item;

	SelFgIntensity( BLACK );
	x_coord = (itof(rgn->r_begin) * x_step);
	x_coord += itof(x_min);
	x = floor(x_coord);
	MoveTo(x, y);
	if (rgn->r_status == ASSIGNED)
	    /* We test this because the regions may not
	    have been tested for rule compliance (if there
	    weren't the required regions for each param.).
	    Thus, then they would be ~ASSIGNED */
	    SelFgIntensity( BLUE );
	else
	    SelFgIntensity( RED );
	x_coord = itof(rgn->r_begin + rgn->r_length);
	x_coord *= x_step;
	x_coord += itof(x_min);
	x = floor(x_coord);
	LineTo(x, y);
    }
    SelFgIntensity( BLACK );
}

/*
	THE END
*/

