/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE	     						*
*	    PILAR_9.C						*
*								*
*	SYNOPSIS						*
*	    PILAR_9 -- ninth PILAR program. Counts the electric	*
*	interactions between two sequences treating them as two	*
*	linear arrays of +/- 1 charges (i.e. just computes  the	*
*	net number of interactions).				*
*	THIS IS AN VT240 VERSION WITH ReGis GRAPHICS SUPPORT.	*
*								*
*	DESCRIPTION						*
*	    This program is invoked with the following command:	*
*	PILAR_9 -p protein -o output -g -n -t			*
*								*
*	    It reads the protein(s) file(s) specified and then	*
*	computes the electrostatic interactions considering	*
*	the(m) as a linear array of amino acids with net charge	*
*	of + - 1.						*
*	-o option not used, then it writes to standard output.	*
*	It also draws several graphics in the screen showing	*
*	the corresponding information if -g is used.		*
*	If the -n flag is given then no report will be shown.	*
*	The -t flag truncates the values at 0-0 offset in self	*
*	interactions.						*
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
*	    15 - feb - 91 (JRValverde) First version.		*
*	     1 - aug - 91 (JRValverde) Ported to Ultrix.	*
*	 							*
*	COPYRIGHT						*
*	    Public domain. Not allowed to earn money from it	*
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
#include <string.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_extern.h"
#include "P_global.h"
#include "P_graphic.h"

/*----------------------------------------------------------------------
			PRIVATE DEFINITIONS
----------------------------------------------------------------------*/

#define MIN_ARG_COUNT	2

#define DEBUG	0
/* #define DEBUG    -8 */

/*----------------------------------------------------------------------
			MODULE'S VARIABLES
----------------------------------------------------------------------*/

char prot1_file[MAXPATH + 1];	    /* protein file name    */
char prot2_file[MAXPATH + 1];	    /* second protein file  */
char out_file[MAXPATH + 1];	    /* output file name	    */

FILE *out_fp;
flag self_flag = NO;
flag trunc_flag = NO;
flag out_flag = NO;
flag grf_flag = NO;
flag silent_flag = NO;
protein_t *prot1, *prot2;
int *values, maxval, minval, val_no;

static char program[MAXPATH + 1];

/*----------------------------------------------------------------------
			FORWARD DEFINITIONS
----------------------------------------------------------------------*/

static void save_name(char *);
static void show_usage();
static void parse_arguments(int, char *[]);
static void assign_charges(protein_t *);
static void save_values();
static void draw_data();
static void init_map();
static void draw_axis();
static void draw_values();

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
	printf("PILAR_9: Count interactions between two sequences.\n");
	printf("    J.R.Valverde - (VT 240) - %s\n", __DATE__);
	printf("====================================================\n\n");
	et_voila();
    }

    save_name("PILAR_9");
    parse_arguments(argc, argv);

    if ( (prot1 = getmem(sizeof(protein_t))) == NULL)
	letalf("%s -- Run out of memory", program);
    if ( (prot2 = getmem(sizeof(protein_t))) == NULL)
	letalf("%s -- Run out of memory", program);

    if (! read_protein(prot1_file, prot1))
	letalf("%s -- Cannot read protein file %s", program, prot1_file);

    if (prot2_file[0]) {
	if (! read_protein(prot2_file, prot2))
	    letalf("%s -- Cannot read protein file %s", program, prot2_file);
    }
    else {
	memcpy(prot2, prot1, sizeof(protein_t) );
	self_flag = YES;
    }

    if (! out_flag)
	out_fp = stdout;
    else {
	out_fp = fopen(out_file, "w+");
	if (out_fp == NULL) {
	    letalf("%s -- Couldn't open %s to save data", program, out_file);
	}
    }

    fprintf(out_fp, "Protein: %s is %d aa long\n",
	prot1->p_name, prot1->p_seqlen);
    fprintf(out_fp, "Protein: %s is %d aa long\n",
	prot2->p_name, prot2->p_seqlen);

    assign_charges(prot1);
    assign_charges(prot2);

    count_interactions(prot1, prot2);

    if (! silent_flag) {
	save_values();
    }

    if (out_flag)
	fclose(out_fp);
    if (grf_flag == YES) {
	printf("\nPress any key to see the graphs...\n");
	draw_data();
    }
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
		    if (! BitTest(opts, BIT_0)) {
			BitSet(opts, BIT_0);
			strcpy(prot1_file, argv[++i]);
		    }
		    else
			strcpy(prot2_file, argv[++i]);
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
		case 'n':
		case 'N':
		    silent_flag = YES;
		    break;
		case 't':
		case 'T':
		    trunc_flag = YES;
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
    if ( ! BitTest(opts, BIT_0)) {
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
	strcpy(program, "PILAR_9");
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
    "\nUsage: %s %cp protein [%co output] [%cg] [%ct] [%cn]\n",
	program, OPTION_DELIMITER, OPTION_DELIMITER,
	OPTION_DELIMITER, OPTION_DELIMITER, OPTION_DELIMITER);
    fprintf(stderr,
	"\tg\tUse this flag to get a graphic output\n%s%s",
	"\tt\tUse this flag to skip 0-0 offset count in self-interactions\n",
	"\tn\tUse this flag to inhibit reports.\n");
}

/*--------------------------------------------------------------------
    ASSIGN_CHARGES
	Assigns net charges of +/- 1 to the charged residues in the
    provided sequence. Considered residues are:
	Basic:	K, R, H
	Acid:	D, E

    J. R. Valverde
--------------------------------------------------------------------*/

void assign_charges(prot)
protein_t *prot;
{
    float *charges;
    int i, len;

    /* allocate array of arrays of values in the protein */
    /* we are using only one array of values */
    prot->p_aaval = (float **) malloc( sizeof(float *) );
    if (prot->p_aaval == NULL)
	letalf("%s -- Not enough memory (ugh!)", program);

    /* allocate array of values */
    len = prot->p_seqlen;
    charges = (float *) calloc(len, sizeof(float));
    prot->p_aaval[0] = charges;

    /* assign charges */
    for (i = 0; i < len; ++i)
	switch (prot->p_seq[i]) {
	    case 'K':
	    case 'k':
	    case 'R':
	    case 'r':
	    case 'H':
	    case 'h':
		charges[i] = 1.0;
		break;
	    case 'D':
	    case 'd':
	    case 'E':
	    case 'e':
		charges[i] = -1.0;
		break;
	    default:
		charges[i] = 0.0;
		break;
	}
}

/*--------------------------------------------------------------------
    COUNT_INTERACTIONS
	Count number of interactions between two proteins considered as
    linear arrays of aminoacids.

    J. R. Valverde
--------------------------------------------------------------------*/

count_interactions(prot1, prot2)
protein_t *prot1;
protein_t *prot2;
{
    int init1, init2,
	end1, end2,
	p1, p2,
	count, pos;
    float *vals1, *vals2;
    extern int *values, maxval, minval;

    /* unneeded */
    if (prot1->p_seqlen < prot2->p_seqlen) {
	protein_t *prot3;	/* swap proteins (locally to the function) */

	prot3 = prot1;
	prot1 = prot2;
	prot2 = prot3;
    }

    val_no = count = prot1->p_seqlen + prot2->p_seqlen - 1;
    if ( (values = (int *) calloc( count, sizeof(int) )) == NULL)
	letalf("\n%s: Not enough memory to hold count of interactions\n",
	    program);
    while (count--)
	values[count] = 0;
    maxval = minval = 0;
    end1 = prot1->p_seqlen - 1;
    end2 = prot2->p_seqlen - 1;
    vals1 = prot1->p_aaval[0];
    vals2 = prot2->p_aaval[0];

    for (init1 = 0, pos = 0; pos <= end2; pos++) {

	/* ****** NOTE ******* */
	/* If self interaction is computed, at offsets 0 - 0
	we will obtain either a huge maximum or a huge minimum.
	Therefore we will ARTIFICIALLY obscure them if the user
	desires to do so.
	    For this we require BOTH, self-interaction and
	user confirmation to be true.			    */

	if (self_flag && trunc_flag)
	    if (pos == end2)
		break;		    /* skip offset 0-0 */
	/* ******************* */

	init2 = end2 - pos;
	count = 0;
	for (p1 = init1, p2 = init2;
	    p1 <= end1 && p2 <= end2;
	    p1++, p2++) {
	    if ((vals1[p1] != 0.0) && (vals2[p2] != 0.0))
		if (vals1[ p1 ] == vals2[ p2 ])
		    count--;
		else
		    count++;
	}
	values[pos] = count;
	if (count > maxval)
	    maxval = count;
	else
	    if (count < minval)
		minval = count;
    }

    /* ***** NOTE **** */
    /* Restore status and store a dummy value for
    offset 0-0 if the user requested it to be ignored
    in a self-interaction computation		    */

    if (self_flag && trunc_flag)
	values[ end2 ] = 0;

    for (init2 = 0, pos = 1; pos <= end1; pos++) {
	/* initial pos = 1 'cos init1=0 and init2=0 has already been
	evaluated in the preceding loop */
	init1 = pos;
	count = 0;
	for (p1 = init1, p2 = init2;
	    p1 <= end1 && p2 <= end2;
	    p1++, p2++) {
	    if ((vals1[p1] != 0.0) && (vals2[p2] != 0.0))
		if (vals1[ p1 ] == vals2[ p2 ])
		    count--;
		else
		    count++;
	}
	values[ end2 + pos ] = count;
	if (maxval < count)
	    maxval = count;
	else
	    if (minval > count)
		minval = count;
    }
}

/*--------------------------------------------------------------------
    SAVE_VALUES
	Save the values of electrostatic interactions count for each
    displacement of one protein over the other.

    J. R. Valverde
--------------------------------------------------------------------*/

static void save_values()
{
    counter i;
    extern FILE *out_fp;
    extern int *values;
    extern int val_no;

    for (i = 0; i < val_no; i++)
	fprintf(out_fp, "%d\n", values[i]);
}

/*--------------------------------------------------------------------
    DRAW_DATA
	This subroutine displays all the data available about the
    interaction in a graphical form to the screen using ReGis graphics
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
    char 	ch;

    ch = getchar();
    initscr();
    clear();
    init_map();

    initscr();
    clear();
    enterReGis(1);
    EraseScreen();
    draw_axis();
    draw_values();
    ch = getchar();
    EraseScreen();
    exitReGis();
}

/*--------------------------------------------------------------------
	INIT_MAP
	    Initialize colour map to the desired colours.

	J. R. Valverde
--------------------------------------------------------------------*/

static void init_map()
{
    return;
#ifdef NOTDEF
    enterReGis(1);
    SelMapHLS(0, 0, 0, 0);	    /* Map location 0 is black	*/
    SelMapHLS(1, 0, 50, 100);	    /* Map location 1 is blue	*/
    SelMapHLS(2, 120, 50, 100);	    /* Map location 2 is red	*/
    SelMapHLS(3, 240, 50, 100);	    /* Map location 3 is green	*/
    exitReGis();
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

    SelFgIntensity( BLACK );
    LineTo(x_min, y_min);
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

static void draw_values()
{
    extern int	*values, maxval, minval, val_no;
    counter	count;
    float	val;
    float	x_coord, y_coord, y_step;
    int		x, y, threshold;

    /* +/- line */
    SelFgIntensity( BLUE );
    y = threshold = ((y_max - y_min) / 2) + y_min;
    MoveTo(x_min, y);
    LineTo(x_max, y);
    SelFgIntensity( BLACK );

    /* Compute limits */
    /* The algorithm used to draw the values will plot them as the net
    sum of attractive minus repulsive interactions with attractive
    interactions being plot upside the zero line, and repulsive downside.
    To achieve a proporcional look, the limits of the plot are adjusted
    according to max ( | maxval |, | minval | ) so that the zero line
    falls just in the middle of the plot.
    Therefore, + limit == - limit == max( maxval, |minval| )	    */

    if (maxval < abs(minval))
	maxval = abs(minval);	    /* we'll use maxval as the limit */

    x_step = itof(x_len) / itof(val_no);
    y_step = itof(y_len / 2) / itof(maxval);

    /* First point */
    MoveTo(x_min, y_max );  
    SelFgIntensity( GREEN );
    val = itof( values[0] );
    x_coord = itof(x_min);
    y_coord = val * y_step;
    y = threshold - floor(y_coord);
    if (val > 0.0) {
	SelFgIntensity( GREEN );
    }
    else {
	SelFgIntensity( RED );
    }
    x = floor(x_coord);
    /* LineTo(x, y) */
    MoveTo(x, y);

    for (count = 1; count < val_no; count++) {
	val = itof( values[count] );
	x_coord += x_step;
	y_coord = val * y_step;
	y = threshold - floor(y_coord);
	if (val > 0.0) {
	    SelFgIntensity( GREEN );
	}
	else {
	    SelFgIntensity( RED );
	}
	x = floor(x_coord);
	LineTo(x, y);
    }
    SelFgIntensity( BLACK );
}

/*
	THE END
*/

