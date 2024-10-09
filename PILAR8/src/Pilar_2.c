/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*  	MODULE:							*
*	    PILAR_2.C						*
*								*
*	DESCRIPTION:						*
*	    This is the main module for the program PILAR_2	*
*	It searches index files and compares them obtaining	*
*	another wich contains the intersection of the specified	*
*	ones.							*
*								*
*	SYNOPSIS:						*
*	    Main module for PILAR_2				*
*								*
*	RETURNS:						*
*	    The program returns with a failure status if	*
*	it can not open any of the files. The value returned	*
*	is representative of the failed file. Otherwise it	*
*	returns with a value of zero.				*
*								*
*	CAVEATS:						*
*								*
* 	FILES:							*
*								*
*	NOTES:							*
*								*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo.			*
*								*
*	HISTORY:						*
*	    21 - jul - 1989 (J.R.Valverde - Written)		*
*	    21 - jul - 1989 (J.R.Valverde - Debugged)		*
*	     1 - aug - 1991 (J.R.Valverde - Ultrix port)	*
*								*
*	COPYRIGHT:						*
* 								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#define DEBUG	0
/* #define DEBUG    -2	*/

#include <stdio.h>
#include <curses.h>

/* Pseudocode:
    Get files to scan (up to 8).
    Since the library was searched sequentially,
    they are already sorted in increasing order.
    Read indexes from every file and compare with the rest.
    Annotate only common indexes:
	Carry a count of indexes read on every file
	Carry a count of common indexes
	Write common indexes.
    Output counts to the user.
*/

#define MAX_INFILES 8
#define MAXPATH	    128
#define EXIT_OK	    1

main()
{                          
    int i;			/* loop counter */
    int max_val_file;		/* file holding the current maximum
				    value */
    int no_files;		/* number of files to scan */
    int match;			/* current number of matches */

    int icount[MAX_INFILES];	/* number of input values per
				    file */
    int ocount;			/* number of output values */

    long max_val;		/* current maximum value */
    long cur_val;		/* current value */
    
    FILE *outfp;		/* pointer to output file */
    FILE *infp[MAX_INFILES];	/* input files */

    char outfilnam[MAXPATH];	/* output file name */
    char infilnam[MAX_INFILES][MAXPATH];    /* input files names */ 
                                               
    /* Begin with a clean sheet (not with a shit) */
    initscr();
    clear();
                  
    /* Hello, world, it's me ! */
    printf("\t\tP I L A R -- (V. 1.0)\n");
    printf("\t\t=====================\n");
    printf("\n\n");
    printf(
"This program scans index files and finds the common proteins\n");

    /* Prompt for the number of files to scan */
    for (;;) {
	printf("\nNumber of files to be scanned: ? ");
	scanf("%d", &no_files);
	if (no_files > MAX_INFILES) {
	    printf("\nNo more than %d input files allowed !",
		MAX_INFILES);
	    continue;
	}
	else
	    break;
    }

    /* Open the input files */
    printf("\nPlease, enter name of files to be scanned :");
    for (i = 0; i < no_files; i++) {
	printf("\n - Name of file # %d : ? ", i + 1);
	scanf("%s", infilnam[i]);
	infp[i] = fopen(infilnam[i], "r");
	if (infp[i] == NULL) {
 	    exit(-(i + 1));
	}
    }

    /* open output file */
    printf("\nPlease, enter name of file to store the intersection:");
    printf("\n - Name of output file: ? ");
    scanf("%s", outfilnam);
    printf("\n");
    outfp = fopen(outfilnam, "w+");
    if (outfp == NULL) {
	exit( 0 );
    }

    /* Now process the input files until feof */
    /* initialize */

    for (i = 0; i < no_files; i++) {
	rewind(infp[i]);
	icount[i] = 0;
    }

    max_val_file = 0;	/* for instance */
    max_val = cur_val = 0L;
    ocount = 0;
    fscanf(infp[max_val_file], "%ld", &max_val);
    match = 1;
    icount[max_val_file]++;

    /* proceed */
    for (;;) {
	for (i = 0; i < no_files; i++) {
	    if (i == max_val_file)
		/* this is the one we are checking against */
		continue;

	    do {
		if (feof(infp[i]))
		    goto Out_Of_Endless_Loop;

		fscanf(infp[i], "%ld", &cur_val);
		icount[i]++;
		if (cur_val == max_val) {
		    /* other match found */
		    match++;
		    if (match == no_files) {
			/* all values match: save */
			fprintf(outfp, "%ld\n", max_val);
			ocount++;

			/* get new max_val */
			max_val_file = 0;	/* I say */
			fscanf(infp[max_val_file], "%ld", &max_val);
			icount[max_val_file]++;
		    }
		}
		else if (cur_val > max_val) {
		    /* we have another current maximum value */
		    max_val_file = i;
		    max_val = cur_val;
		    match = 1;
		}
	    }
	    while (cur_val < max_val);
	}
    }
Out_Of_Endless_Loop:

    fclose(outfp);
    for (i = 0; i < no_files; i++)
	fclose(infp[i]);

    printf("\n-------------------------------------------------\n");
    printf("\t\tSummary");
    printf("\n-------------------------------------------------\n");

    printf("\nNumber of proteins found in each file:\n");
    for (i = 0; i < no_files; i++)
	printf(" -- File %s: %d proteins\n", infilnam[i], icount[i]);

    printf("\nRemaining common proteins written to output file:\n");
    printf(" -- File %s: %d proteins\n", outfilnam, ocount);

    printf("\n\nBye\n");

    exit(EXIT_OK);
}
