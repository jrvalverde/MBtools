/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*									*
*	MODULE:								*
*	    P_DESCRIPTOR.C						*
*									*
*	DESCRIPTION:							*
*									*
*	SYNOPSIS:							*
*	    read_descriptor() -- Reads descriptor relevant info.	*
*	    descriptor_name() -- Return name of descriptor file.	*
*									*
*	CAVEATS:							*
*	    Debugged?							*
*									*
*	NOTES:								*
*									*
*	DESIGNED BY:							*
*	    José Ramón Valverde Carrillo.				*
*									*
*	HISTORY:							*
*	    4 - oct - 91 -- JRValverde. Added support for 'aliases' in	*
*	file names through global environment variables.		*
*									*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#ifndef DEBUG
/* #define DEBUG    2 */
#define DEBUG	   -1
#endif

/*--------------------------------------
	GENERAL DEFINITIONS
--------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

extern char *getenv(char *name);

/*--------------------------------------
	PARTICULAR DEFINITIONS
--------------------------------------*/

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"

#include "P_list.h"
#include "P_parser.h"
#include "P_utils.h"

extern int aa_table[];

/*---------------------------------------
	FORWARD DEFINITIONS
---------------------------------------*/

extern int read_feature(char *, feature_t *);

extern status get_feature(feature_t *);

extern status feature_name(char *);

static bool get_dfeatures(FILE *, descriptor_t *);

static bool get_drules(FILE *, descriptor_t *);

static bool get_new_rule( descriptor_t * );

static bool get_distance_rule( distance_rule_t * );

static bool get_sequence_rule( sequence_rule_t * );

static bool get_profile_rule( profile_rule_t * );

/*
	Gracias, Petenera mía;
	por tus ojos me he perdido;
	era lo que yo quería.
			Abel Martín
*/

/*---------------------------------------------------------------

    READ_DESCRIPTOR

    Reads descriptor relevant info.

----------------------------------------------------------------*/

/*	PSEUDOCODE:

read_descriptor(file, descr)
char *file;
descriptor_t *descr;
{

    Read descriptor:

    1 -- Features to evaluate:
        1.1 -- Number of features
        1.2 -- Name of files containing the values
        1.3 -- For every feature:
	    Read feature
    2 -- Limit values
        2.1 -- For every feat.:
	    threshold
	    window size
	    size limits
        2.2 -- error limits.
	2.3 -- minimum protein length
 	2.4 -- moment angle ( < 0 degrees ==> NO MOMENT TO CONSIDER)
    3 -- Rules to evaluate.
        3.1 -- While there exist tokens:
	    new rule
	    allocate space for it
	    read rule values
	3.2 -- Setup descriptor data structures
	    regions map

}
*/

static char d_filnam[MAXPATH + 1] = {'\0'};
static int max_reg = 0;

status read_descriptor(filename, descr)
char *filename;
descriptor_t *descr;
{
    FILE *d_file;
    region_t ***map;
    int i, no_features;
    char *auxname;

    /* open file for read */
    if ((d_file = fopen(filename, "r")) == NULL) {
	if ((auxname = getenv(filename)) == NULL) {
	    return 0;
	}
	else {
	    if ((d_file = fopen(auxname, "r")) == NULL) {
		return 0;
	    }
	} 
    }
    
    /* save descriptor name */
    strncpy(d_filnam, filename, MAXPATH);

    if (! get_dfeatures(d_file, descr)) return FALSE;

    if (! get_drules(d_file, descr)) return FALSE;

    fclose(d_file);

    /* allocate map of regions to be covered */

    /*	first, the array of feature regions arrays */
    no_features = list_length(descr->d_features);
    map = (region_t ***) calloc( no_features, sizeof(region_t **));
    if (map == NULL)
	letal("READ DESCRIPTOR", "not enough memory");

    /*	second, the feature regions arrays: these are
    arrays of pointers to the actual regions */
    for (i = 0; i < no_features; ++i) {
	map[i] = (region_t **) calloc( max_reg, sizeof(region_t *));
	if ( map[i] == NULL)
	    letal("READ DESCRIPTOR", "not enough memory");
    }

    descr->d_regions = map;
    descr->d_maxregno = max_reg;

    return TRUE;
}

/*------------------------------------------------------------------*/

static bool get_dfeatures(d_file, descr)
FILE *d_file;
descriptor_t *descr;
/*
*   get features for the descriptor specified from the file indicated
*
*   INPUT:  D_FILE  -- File to read descriptor features from
*	    DESCR   -- Structure to store the information obtained
*   OUTPUT: TRUE if all has gone OK, FALSE otherwise
*
*   J. R. Valverde
*/
{
    counter count, i;		/* loop counters			    */
    int minplen = 0;		/* minimum protein length for this feature  */
    int no_data = 0;		/* number of features to scan		    */
    feature_t *pfeat;		/* the feature being read		    */
    char f_name[MAXPATH + 1];	/* name of the file with feat's values  */
    char mode = '+';		/* should we invert feature's values ?	    */
                                          
    extern int max_reg;

    f_name[0] = '\0';
    /* get features to evaluate */
    fscanf(d_file, "%*s%d\n", &no_data);
    if (no_data == 0) {
	return FALSE;
    }
    descr->d_features = new_list();
    descr->d_minplen = 0;
    /* read each feature. If not, continue */
    for (count = 0; count < no_data; count++) {

	fscanf(d_file, "%s %c", f_name, &mode);
	pfeat = (feature_t *) NEW( feature_t );
	if (! read_feature(f_name, pfeat) ) {
	    free(pfeat);
	    for (i = 0; i < 5; ++i)
		/* skip its arguments */
		fscanf(d_file, "%*s%*s");
            continue;
	}

	if (mode == '-') {
	    antinormalize(pfeat->f_values, AA_NUMBER);
	    pfeat->f_inverted = TRUE;
	}
	else {
	    normalize(pfeat->f_values, AA_NUMBER);
	    pfeat->f_inverted = FALSE;
	}

	pfeat->f_name = (char *) calloc((strlen(f_name) + 1), sizeof(char));
	strcpy(pfeat->f_name, f_name);
	fscanf(d_file, "%*s%d", &(pfeat->f_window));
	fscanf(d_file, "%*s%f", &(pfeat->f_threshold));
	fscanf(d_file, "%*s%d", &(pfeat->f_minrglen));
	fscanf(d_file, "%*s%d", &(pfeat->f_maxrglen));
	fscanf(d_file, "%*s%d", &(pfeat->f_no_regions));
	fscanf(d_file, "%*s%f",	&(pfeat->f_moment));

	/* state general information:
	    normalize threshold,
	    max number of regions
	    and minimum protein length. */
	pfeat->f_threshold = (pfeat->f_threshold - pfeat->f_minvalue)
			    / (pfeat->f_maxvalue - pfeat->f_minvalue);
	if (pfeat->f_inverted == TRUE)
	    pfeat->f_threshold = 1.0 - pfeat->f_threshold;
	if (pfeat->f_no_regions > max_reg)
	    max_reg = pfeat->f_no_regions;
	minplen = pfeat->f_no_regions * pfeat->f_minrglen;
	if (minplen < pfeat->f_window)
	    minplen = pfeat->f_window;
	if (minplen > descr->d_minplen)
	    descr->d_minplen = minplen;
	list_append(pfeat, descr->d_features);
    }
    return TRUE;
}

/*-----------------------------------------------------------
    READ_FEATURE

    Reads data values from the especified feature
    values file.

    INPUT: F_FILNAM: Name of file containing values.
	Values must be stored in the following format:
	"AminoAcid (1 letter code) HTAB Value CR/LF"

    OUTPUT: FEATURE: Feature table filled with the values.

    RETURNS: Number of amino acid values read.

    J. R. Valverde
-----------------------------------------------------------*/

static char f_filnam[MAXPATH + 1] = {'\0'};

status read_feature(f_filename, pfeat)
char *f_filename;
feature_t *pfeat;
{                           
    FILE *f_file;
    char aa, *auxname;
    float *feature;
    int no_aa_read;
    float aa_value;
    float min, max;

    /* first try to open the file */
    if ((f_file = fopen(f_filename, "r")) == NULL) {
	/* Maybe it was an alias */
	if ( (auxname = getenv(f_filename)) == NULL ) {
	    /* No, it wasn't */
	    horror("Read feature", "Cannot open specified file");
	    fprintf(stderr, " *** %s ***\n", f_filename);
	    return FALSE;
	}
	else {
	    /* Yes, try with it */
	    if ((f_file = fopen(auxname, "r")) == NULL) {
		horror("Read feature", "Cannot open specified file");
		fprintf(stderr, " *** %s ***\n", f_filename);
		return FALSE;
	    }
	}
    }


    /* save feature file name */
    strcpy(f_filnam, f_filename);

    /* initialize feature table to zeroes */
    feature = pfeat->f_values;
    for (no_aa_read = 0; no_aa_read < AA_NUMBER; no_aa_read++)
	feature[no_aa_read] = 0.0;

    /* Loop until all values have been read */
    no_aa_read = 0;
    min = F_INFINITY;
    max = - F_INFINITY;
    do {

      if (feof(f_file)) break;

      /* get next value and check that reading has gone OK */
      if (fscanf(f_file, "%c %f", &aa, &aa_value) != 2) {
	/* if not, then retry loop */
	continue;
      }

      /* update table: the value is stored in the position
      specified by the aminoacid letter (both, in upper
      and lower case).
      We force the value to fit into the table by making
      an AND of the position with 0x7F, which sets it
      between 0 and 127 */
      aa = ASCII(aa);
      if (aa_table[aa] > 0) {
        /* set the value in position */
        feature[aa_table[aa] - 1] = aa_value;
        no_aa_read++;
	/* save max or min */
	if (aa_value > max)
	    max = aa_value;
	if (aa_value < min)
	    min = aa_value;
      }
    }
    while (no_aa_read < AA_NUMBER);
    pfeat->f_maxvalue = max;
    pfeat->f_minvalue = min;

    fclose(f_file);

    return (no_aa_read);
}


/*--------------------------------------------------------------

    GET_FEATURE

    Asks the user for the feature to be used in the search
    and reads it from the file.

    INPUT: THE_FEATURE  table to read data values into.
    OUTPUT: SUCCESS if read's been OK. FAIL otherwise.

    J. R. Valverde
---------------------------------------------------------------*/

status get_feature(the_feature)
feature_t *the_feature;
{
    char *f;
    extern char f_filnam[MAXPATH + 1];

    printf("\nName of file containing feature values: ? ");
    scanf("%s", f_filnam);
    if ((f = rindex(f_filnam, ENDDIRCHAR)) == NULL)
	f = f_filnam;
    if (rindex(f, FILEXTCHAR) == NULL)
	strcpy(f_filnam, F_EXT);

    return ( read_feature(f_filnam, the_feature) );
}                     

/*---------------------------------------------------------------
    FEATURE_NAME
	Returns the name of the last feature read.

    INPUT:  F_NAME -- String to store the name.
    OUTPUT: SUCCESS if so, FAIL if there is not current
	feature name.

    CAUTION: This function is HIGHLY SYSTEM DEPENDENT.
            Supposes that f_name is long enough to hold
	the name.

    J. R. Valverde
---------------------------------------------------------------*/

status feature_name(f_name)
char *f_name;
{
    char *nam;
    extern char f_filnam[MAXPATH + 1];

    if (f_filnam[0] == '\0')
	return FAIL;

    nam = rindex(f_filnam, ENDDIRCHAR);
    if (nam == NULL)
	nam = f_filnam;
    else
	nam++;

    strcpy(f_name, nam);
}


/*-----------------------------------------------------------------*/

static bool get_drules(d_file, descr)
FILE *d_file;
descriptor_t *descr;
/*
 * get_drules()
 *  get from a description file the rules for a descriptor.
 *
 *  INPUT:  D_FILE  -- Description file with position pointer set to
 *		beginning of rules section
 *	    DESCR   -- Descriptor structure to hold the data
 *  OUTPUT: TRUE if success, FALSE if not.
 *
 */
{
    char *data;
    flag type;

    parse(d_file);

    descr->d_rules = new_list();
    if (descr->d_rules == NULL) return FALSE;
    data = retrieve_token( &type );
    if (type == NO_MORE_TOKENS) {
	/* There are no rules ! */ 
	return TRUE;
    }
    while (get_new_rule(descr));
    return TRUE;
}

/*-----------------------------------------------------------------*/

static bool get_new_rule(descr)
descriptor_t *descr;
/*
 * get_new_rule()
 *  retrieve another rule from the set of tokens read
 *
 *  INPUT:  Descriptor into which we are to store the rule
 *  OUTPUT: TRUE if rule could be successfully read, FALSE if not
 */
{
    rule_t *the_rule;
    char *data = NULL;
    flag type;
    bool more = TRUE;	    

    the_rule = NEW( rule_t );
    /* get rule type */
    data = retrieve_token( &type );
    /* data should be NULL */
    if (type == NO_MORE_TOKENS) {
	free(the_rule);
	return FALSE;	    /* We have finished */
    }
    the_rule->type = type;
    switch (type) {
	case DISTANCE:
	    more = get_distance_rule( &(the_rule->data.distance_rule) );
	    break;
	case SEQUENCE:
	    more = get_sequence_rule( &(the_rule->data.sequence_rule) );
	    break;
	case PROFILE:
	    more = get_profile_rule( &(the_rule->data.profile_rule) );
	    break;
	default:
	    free(the_rule);
     	    return FALSE;   /* and finish */
    }
    if (! list_append(the_rule, descr->d_rules)) {
	letal("GET RULE", "cannot add rule to list");
    }
    return (more);
}

/*-----------------------------------------------------------------*/

static bool get_distance_rule( r )
distance_rule_t *r;
/*
 * get_distance_rule()
 *  retrieve a rule of separation between regions from the set of
 *  tokens read.
 *
 *  INPUT:  R	-- Structure to store the data about the rule.
 *  OUTPUT: TRUE if everything has gone OK, FALSE otherwise.
 */
{
    char *data;
    flag type;
    extern int max_reg;

    data = retrieve_token( &type );	/* parameter 1 */
    r->feat1 = atoi(data) -1;		/* internally we use 0-offset */
    free(data); 

    data = retrieve_token( &type );	/* region 1 */
    r->reg1 = atoi(data) -1;
    if (r->reg1 > max_reg)
	max_reg = r->reg1;
    free(data);

    data = retrieve_token( &type );	/* position 1 */
    r->pos1 = atoi(data) -1;
    free(data);

    data = retrieve_token( &type );	/* parameter 2 */
    r->feat2 = atoi(data) -1;
    free(data);

    data = retrieve_token( &type );	/* region 2 */
    r->reg2 = atoi(data) -1;
    if (r->reg2 > max_reg)
	max_reg = r->reg2;
    free(data);

    data = retrieve_token( &type );	/* position 2 */
    r->pos2 = atoi(data) -1;
    free(data);

    data = retrieve_token( &type );	/* minimum distance */
    r->dmin = atoi(data);
    free(data);

    data = retrieve_token( &type );	/* maximum distance */
    r->dmax = atoi(data);
    free(data);

    /* check for more rules */
    data = retrieve_token( &type );
    if (type == NO_MORE_TOKENS)
	return FALSE;	    /* finished */
    else
	return TRUE;
}

/*-----------------------------------------------------------------*/

static bool get_sequence_rule( r )
sequence_rule_t *r;
/*
 * get_sequence_rule()
 *  Get a rule specifying a sequence to be searched in a region
 *  from the stream of tokens read.
 *
 *  INPUT:  R -- structure to store the data got
 *  OUTPUT: TRUE if no problem has been found.
 *	    FALSE otherwise
 */
{
    char *data;
    flag type;
    extern int max_reg;

    /* the sequence is a special case: we should NOT free data */
    data = retrieve_token( &type );	/* sequence */
    r->sr_seq = data;

    data = retrieve_token( &type );	/* parameter */
    r->sr_feat = atoi(data) -1;		/* internally we use 0-offset */
    free(data);

    data = retrieve_token( &type );	/* region */
    r->sr_reg = atoi(data) -1;
    if (r->sr_reg > max_reg)
	max_reg = r->sr_reg;
    free(data);

    data = retrieve_token( &type );	/* minimum position */
    r->sr_minpos = atoi(data) -1;
    free(data);

    data = retrieve_token( &type );	/* maximum position withn region */
    r->sr_maxpos = atoi(data) -1;
    free(data);

    /* check if we have more rules to add */
    data = retrieve_token( &type );
    if (type == NO_MORE_TOKENS)
	return FALSE;
    else
	return TRUE;
}

/*-----------------------------------------------------------------*/

static bool get_profile_rule( r )
profile_rule_t *r;
/*
 * get_profile_rule()
 *  Get a rule that describes a profile to adjust a given region to
 *  from the set of tokens found in the description file
 *
 *  INPUT:  R -- Structure to store the information retrieved
 *  OUTPUT: TRUE if the information has been succesfully obtained
 *	    FALSE if any error is detected.
 */
{
    char *data;
    flag type;
    int i;
    float *val_arr;
    extern int max_reg;
    

    data = retrieve_token( &type );	/* parameter */
    r->pr_feat = atoi(data) -1;
    free(data);

    data = retrieve_token( &type );	/* region */
    r->pr_reg = atoi(data) -1;
    if (r->pr_reg > max_reg)
	max_reg = r->pr_reg;
    free(data);

    data = retrieve_token( &type );
    val_arr = calloc((MAX_SEQ_LEN + 1), sizeof(float));
    i = 0;
    while (type == VALUE) {
	sscanf(data, "%f", &val_arr[i]);
	i++;
	free(data);
	data = retrieve_token( &type );
    }
    /* i points now one element beyond the last one, which is mindiff */
    r->pr_mindiff = val_arr[--i];
    /* all previous values are those for the region */
    /* i now points to mindiff as a subindex, but, since indexes are
    zero offset, it is also just the length of the array */
    r->pr_values = realloc(val_arr, i * sizeof(float));
    r->pr_len = i;	/* the array begins at zero */

    /* check if there remains any rule */
    /* type != VALUE... Is it NO_MORE_TOKENS or NEW_RULE ??? */
    if (type == NO_MORE_TOKENS)
	return FALSE;
    else
	return TRUE;
}

/*------------------------------------------------------------------*/
/*		    D E B U G G I N G 				    */
/*------------------------------------------------------------------*/

#if DEBUG == 2

#include "P_global.h"

char infile[80];
char outfile[80];
FILE *fp1;
FILE *fp2;
int no_aa, count;                       
feature_table_t feat1, feat2;
descriptor_t d;
l_head_t *the_rules;
rule_t *my_rule;
distance_rule_t *dr;
sequence_rule_t *sr;
profile_rule_t  *pr;

main()
{
    static int i, j, end;

    puts("TESTING MODULE P_DESCRIPTOR");
    puts("===========================\n");

    puts("READ_FEATURE");
    printf("Name of feature file: ? ");
    scanf("%s", infile);
    if ((no_aa = read_feature(infile, feat1)) == FAIL)
	exit(-2);
    printf("\nFile contained values for %d aa.\n", no_aa);
    for (count = 0; count < AA_NUMBER; count++)
	printf("%c --> %2.3f\n", aa_1[count], feat1[count]);

    puts("GET_FEATURE");
    if ((no_aa = get_feature(infile, &feat2[20])) == FAIL)
	exit(-3);
    printf("\nFile contained values for %d aa.\n", no_aa);
    for (count = 0; count < AA_NUMBER; count++)
	printf("%s --> %2.3f\n", aa_3[count], feat2[count]);

    puts("READ_DESCRIPTOR");
    read_descriptor("DESCR.TEST", &d);

    printf("\n\nDescriptor: DESCR.TEST\n\n");
    printf("Number of features: %d\n", list_length(d.d_features));
    printf("Number of rules: %d\n", list_length(d.d_rules));

    the_rules = d.d_rules;
    end = list_length(d.d_rules);
    for (i = 0; i < end; ++i) {
	my_rule = (rule_t *) list_del_first(the_rules);
	printf("Rule no. %d: \n", i + 1);
	switch (my_rule->type) {
	    case DISTANCE:
		dr = &my_rule->data;
		printf("distance (%d, %d, %d) -- (%d, %d, %d) = %d .. %d\n",
		    dr->feat1, dr->reg1, dr->pos1,
		    dr->feat2, dr->reg2, dr->pos2,
		    dr->dmin, dr->dmax);
		break;
	    case SEQUENCE:
		sr = &my_rule->data;
		printf("sequence \"%s\" in (%d, %d, %d .. %d)\n",
		    sr->sr_seq, sr->sr_feat, sr->sr_reg,
		    sr->sr_minpos, sr->sr_maxpos);
		break;
	    case PROFILE:
		pr = &my_rule->data;
		printf("profile (%d, %d) (", pr->pr_feat, pr->pr_reg);
		for (j = 0; j < pr->pr_len; ++j)
		    printf("%g, ", pr->pr_values[j]);
		printf(") < %g\n", pr->pr_mindiff);
	 	break;
	    default:
		printf("\n*** UNKNOWN RULE ***\n\n");
		break;
	}
    }
    exit(EXIT_OK);
}

#endif
