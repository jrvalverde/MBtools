/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE:	P_UTILS.C					*
*								*
*	    Contains several utility routines for PILAR.	*
*								*
*	SYNOPSIS:						*
*		iMaxMin()					*
*		fMaxMin()					*
*		dMaxMin()					*
*		rMaxMin()					*
*								*
*		moment()					*
*		KSone()						*
*		KStwo()						*
*								*
*		normalize()					*
*		antinormalize()					*
*								*
*		et_voila()					*
*		horror()					*
*		letal()						*
*		letalf()					*
*								*
*		index()						*
*		rindex()					*
*		tmpname()					*
*		getmem()					*
*		strindex()					*
*		is_empty()					*
*								*
*		getch()						*
*		getche()					*
*		ungetch()					*
*								*
*	CAVEATS:						*
*		rMaxMin, KSone and KStwo are still untested.	*
*								*
*	NOTES:							*
*		use i, f, dMAXMIN for integer, float or double	*
*	arrays. For more complex data use rMaxMin. This is a	*
*	recursive function which must be used when cost of	*
*	comparing items is elevated.				*
*		KSone and KStwo are for CUMULATIVE distri-	*
*	butions and MUST receive parameter data arrays SORTED	*
*	IN ASCENDING ORDER.					*
*		Functions INDEX and RINDEX can be eliminated in	*
*	environments compatible with UNIX - V7, BSD-UNIX or	*
*	XENIX 2.3. Or they can be substituted by strchr() and	*
*	strrchr() on System III or V and XENIX 3.0 compatible	*
*	environments.						*
*		Functions GETCH, GETCHE and UNGETCH are used	*
*	for buffered input from keyboard. They REQUIRE CURSES	*
*	to be initialized through INITSCR() before using them.	*
*	These functions must be commented under MS-DOS.		*
*								*
*	SEE ALSO:						*
*	    Numerical Recipes. The Art of Scientific		*
*	    Computing.	Cambridge University Press.		*
*								*
*	HISTORY:						*
*	    20 - jun - 1989 (First written on an IBM-AT)	*
*	    10 - jul - 1989 (normalize/anti.. added)		*
*	    12 - jul - 1989 (Partially debugged)		*
*	    17 - jul - 1989 (Ported to VAX and debugged)	*
*	    20 - jul - 1989 (index / rindex added)		*
*	    16 - sep - 1989 (Buffered keyboard input added)	*
*	    08 - apr - 1990 (Added strindex and is_empty)	*
*	    24 - may - 1990 (Added et_voila)			*
*	     1 - aug - 1991 (Ported to ULTRIX)			*
*								*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo.			*
*								*
*	Parts have been taken from "Press, Flannery, Teukolsky	*
*	and Vetterling. Numerical Recipes. Cambridge University	*
*	Press." and are therefore property of the authors.	*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#ifndef DEBUG
#define DEBUG    -1
/* #define DEBUG	    4  */
#endif


#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <malloc.h>
#include <stdarg.h>
#include <curses.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"


#define MAX(x, y)   ((x) > (y))? (x) : (y)

/* extern declarations */
extern void horror(const char *, const char *);	/* forward definition */

/* forward declarations */
static _MaxMin(counter i, counter j, void **fmax, void **fmin);
static void *_max(void *item1, void *item2);
static void *_min(void *item1, void *item2);
static float probKS(float alam);


/*---------------------------------------------------------------
    IMAXMIN:
	Find maximum and minimum of an integer array.

    INPUT:  ARR -- array of data to search
	    N -- number of data in the array
    OUTPUT: MAX -- maximum value found
	    MIN -- minimum value found
    J. R. Valverde
---------------------------------------------------------------*/

void iMaxMin(arr, n, imax, imin)
int arr[];	/* array[1..n] of integer values */
counter n;
int *imax, *imin;
/*
    Set MAX to the maximum and MIN to the minimum of ARR[1..N]
*/                                                              
{
    register counter i;

    *imax = *imin = arr[0];
    for (i = 1; i < n; i++)
	if (arr[i] > *imax)
            *imax = arr[i];
    else
        if (arr[i] < *imin)
	    *imin = arr[i];
}

/*---------------------------------------------------------------
	FMAXMIN:
	    Find maximum and minimum of a float array.

	INPUT:	ARR -- array of data to search
		N -- number of data in the array
	OUTPUT:	MAX -- maximum value found
		MIN -- minimum value found
	J. R. Valverde
---------------------------------------------------------------*/

void fMaxMin(arr, n, fmax, fmin)
float arr[];	    /* array[1..n] of integer values */
counter n;
float *fmax, *fmin;
/*
    Set MAX to the maximum and MIN to the minimum of ARR[1..N]
*/
{
    register counter i;

    *fmax = *fmin = arr[0];
    for (i = 1; i < n; i++)
	if (arr[i] > *fmax)
            *fmax = arr[i];
    else
        if (arr[i] < *fmin)
	    *fmin = arr[i];
}

/*---------------------------------------------------------------
    DMAXMIN:
	    Find maximum and minimum of a double array.

    INPUT:  ARR -- array of data to search
	    N -- number of data in the array
    OUTPUT: MAX -- maximum value found
	    MIN -- minimum value found
    J. R. Valverde
---------------------------------------------------------------*/

void dMaxMin(arr, n, dmax, dmin)
double arr[];	    /* array[1..n] of integer values */
counter n;
double *dmax, *dmin;
/*
    Set MAX to the maximum and MIN to the minimum of ARR[1..N]
*/
{
    register counter i;

    *dmax = *dmin = arr[0];
    for (i = 1; i < n; i++)
	if (arr[i] > *dmax)
            *dmax = arr[i];
    else
        if (arr[i] < *dmin)
	    *dmin = arr[i];
}

/*---------------------------------------------------------------
    RMAXMIN:
    Recursively find maximum and minimum of an item array.
    This function requires on the average approximately 25 %
    less comparisons than the previous ones.
    However, due to the overhead from recursion, it may be
    slower unless comparisons between items are much more
    expensive than comparisons between indexes.

    INPUT:  ARR -- array of data to search
	    N -- number of data in the array
    OUTPUT: MAX -- maximum value found
	    MIN -- minimum value found
    J. R. Valverde
---------------------------------------------------------------*/

static void **_mmarr;	/* auxiliary array of pointers to items */

static int (*_comp)();	/* auxiliary pointer to comparation function */

void rMaxMin(arr, n, rmax, rmin, compare)
void *arr[];		/* array[1..n] of pointers to items */
counter n;
void **rmax, **rmin;
int (*compare)();	/* a comparation function returning
			    -1 if first arg < second arg
			     0 if both equal
			    +1 if first arg > second arg */
{
    _mmarr = arr;
    _comp = compare;

    _MaxMin(0, n, rmax, rmin);
}

static _MaxMin(i, j, fmax, fmin)
counter i, j;
void **fmax, **fmin;
/*
    _mmarr is a global array of N items.
    _comp is a global comparation function.
    Parameters I and J are indexes such that 1 <= i <= j <= n.
    The result is asigning to FMAX and FMIN the upper and lower
    values in _MMARR[I..J] respectively.
*/
{
    counter mid;
    void *gmax, *gmin, *hmax, *hmin;

    if (i == j) {
	*fmax = *fmin = _mmarr[i];
    } else

    if (i == j - 1) {
	if ( (*_comp)(_mmarr[i], _mmarr[j]) == -1) {
	    /* arr[i] < arr[j] */
      	    *fmax = _mmarr[j];
	    *fmin = _mmarr[i];
	}
	else {
      	    *fmax = _mmarr[i];
      	    *fmax = _mmarr[j];
	}
    }

    else {
	mid = (i + j) >> 1;
	_MaxMin(i, mid, &gmax, &gmin);
	_MaxMin(mid + 1, j, &hmax, &hmin);
	*fmax = _max(gmax, hmax);
	*fmin = _min(gmin, hmin);
    }
}

static void *_max(item1, item2)
void *item1, *item2;
{
    if ((*_comp)(item1, item2) == 1)
      return item1;
    else
      return item2;
}

static void *_min(item1, item2)
void *item1, *item2;
{
    if ((*_comp)(item1, item2) == -1)
      return item1;
    else
      return item2;
}

/*
    We advise you, Loddfafnir, to follow the counsel,
    It will do you well, if you take it,
    It will do you good, if you follow it,
    
	Snorri Sturluson, Havamal.
*/

/*---------------------------------------------------------------
    MOMENT:
	Subroutine for calculating moments of a distribution.
	Taken from:
	    Numerical Recipes. The Art of Scientific Computing.

---------------------------------------------------------------*/

#define CUBE(x)	    ((x) * (x) * (x))

#define SQUARE(x)   ((x) * (x))

void moment(data, n, ave, adev, sdev, var, skew, kurt)
float *data;		/* array[1..N] of data */
counter n;		/* number of data */
float *ave, 		/* average of distribution */
      *adev,		/* average deviation */
      *sdev,		/* standard deviation */
      *var,		/* variance */
      *skew,		/* skewness */
      *kurt;		/* kurtosis */
/*
    Given an array of DATA of length N, this routine returns its mean
    AVE, average deviation ADEV, standard deviation SDEV, variance VAR,
    skewness SKEW, and kurtosis KURT.
*/
{
    float s, p;
    counter j;

    if (n <= 1) {
	horror("MOMENT", "N must be at least 2");
	return;
    }

    /* first compute mean */
    s = 0.0;
    for (j = 0; j < n; j++)
	s += data[j];
    *ave = s / n;

    /* then initialize all the rest */
    *adev = *var = *skew = *kurt = 0.0;

    for (j = 0; j < n; j++) {
	s = data[j] - (*ave);
	p = s * s;
	*var += p;
	p *= s;
	*skew += p;
	p *= s;
	*kurt += p;
    }

    *adev /= n;
    *var  /= (n - 1);
    *sdev = sqrt(*var);

    if (*var != 0.0) {
	*skew /= (n + CUBE(*sdev));
	*kurt = *kurt / (n * SQUARE(*var)) - 3.0;
    }
    else {
	horror("MOMENT", "No skew or kurtosis when zero variance");
	*skew = 0.0;
	*kurt = 0.0;
    }
}

#undef SQUARE

#undef CUBE

/*---------------------------------------------------------------
    KOLMOGOROV-SMIRNOV STIMATION:
	Subroutines for comparing two distributions of data.
	Taken from:
	    Numerical Recipes. The Art of Scientific Computing.

---------------------------------------------------------------*/

void KSone(data, n, func, d, prob)
float data[];	    /* array [1 .. n] */
counter n;
float (*func)();
float *d,
      *prob;
/*
    Given an array of N values, DATA, and given an user-supplied
    function of a single variable FUNC which is a cumulative
    distribution function ranging from 0 (for smallest values of
    its arguments) to 1 (for largest values of its arguments),
    this routine returns the K-S statistic D, and the significance
    level PROB. Small values of PROB show that the cumulative
    distribution function of DATA is significantly different from
    FUNC.
    The array DATA MUST BE sorted in ascending order prior to calling
    this routine.
*/
{
    counter en, j;
    float f0, fn, ff, dt;

    en = n;
    *d = 0.0;
    f0 = 0.0;		/* Data's c.d.f. before the next step	*/

    for (j = 0; j < n; j++) {	/* loop over the sorted data points */
	fn = j / en;		/* Data's c.d.f. after this step    */
	ff = (*func)(data[j]);	/* Compare to user supplied function	*/
	dt = MAX(fabs(f0 - ff), fabs(fn - ff));	/* maximum distance */
	if (dt > *d) *d = dt;
	    f0 = fn;
    }

    *prob = probKS( sqrt(en) * (*d));	/* compute significance	    */
}

public KStwo(data1, n1, data2, n2, d, prob)
float data1[];
counter n1;
float data2[];
counter n2;
float *d;
float *prob;
/*
    Given an array DATA1 of N1 values, and an array DATA2 of N2
    values, this routine returns the K-S statistic D, and the
    significance level PROB for the null hypothesis that the data
    sets are drawn from the same distribution. Small values of
    PROB show that the cumulative distribution function of DATA1 is
    significantly different from that of DATA2. The arrays DATA1
    and DATA2 MUST BE sorted before calling this routine.
*/
{
    counter en1, en2, j1, j2;
    float fn1, fn2, d1, d2, dt;

    en1 = n1;
    en2 = n2;
    j1 = j2 = 0;    	    /* next value of DATA1 and DATA2 to be  */
    fn1 = fn2 = 0.0;	    /*  processed.		*/
    *d = 0.0;

    while ( (j1 < n1) && (j2 < n2) ) {	/* if we are not done ...   */
	d1 = data1[j1];
	d2 = data2[j2];
	if (d1 <= d2) {	    /* next step is in DATA1	    */
            fn1 = j1 / en1;
            j1++;
	}
	if (d2 <= d1) {	    /* next step is in DATA2	    */
            fn2 = j2 / en2;
            j2++;
	}
	dt = fabs(fn2 - fn1);
	if (dt > *d) *d = dt;
    }

    *prob = probKS(sqrt(en1 * en2 / en1 + en2) * (*d));
}

static float probKS(alam)
float alam;
{
#define	EPS1	0.001
#define	EPS2	1.0e-8

    float a2, fac, tmp, term, termbf;
    counter j;

    a2 = -2.0 * alam * alam;
    fac = 2.0;
    tmp = 0.0;
    termbf = 0.0;	/* previous term in sum */

    for (j = 1; j <= 100; j++) {
	term = fac * exp(a2 *j * j);
	tmp += term;
	if ((fabs(term) <= EPS1 * termbf) || (fabs(term) <= EPS2 * tmp))
            return tmp;
	fac = -fac;	/* alternating signs in sum */
	termbf = fabs(term);
    }

    return 1.0;		/* get here only by failing to converge */
}


/*-------------------------------------------------------------------
    NORMALIZE
	Given an array of float data this procedure
    finds the maximum and minimum and normalizes all the
    values to fall between 0.0 and 1.0

    INPUT:  DATA    -- Array of values to normalize
	    N	-- The dimension of the array
    OUTPUT: The array DATA has all its values normalized.

    J. R. Valverde
-------------------------------------------------------------------*/

void normalize(data, n)
float *data;
counter n;
{
    float max, min, tmp;
    register counter i;

    fMaxMin(data, n, &max, &min);
    tmp = (max - min);
    for (i = 0; i < n; i++) {
	data[i] -= min;
	data[i] /= tmp;
    }
}

/*-------------------------------------------------------------------
    ANTINORMALIZE
	Anti-normalizes an array of data (i. e. instead of
    assigning 0 to the lowest and 1 to the highest value it does
    the inverse, assigning 1 to the lowest and 0 to the highest
    value).

    INPUT:  DATA    -- The array of values to be anti-normalized
	    N	-- The number of values in DATA
    OUTPUT: DATA is anti-normalized.

    J. R. Valverde
-------------------------------------------------------------------*/

void antinormalize(data, n)
float *data;
counter n;
{
    counter i;
    float max, min, tmp;

    fMaxMin(data, n, &max, &min);
    tmp = max - min;
    for (i = 0; i < n; i++) {
	data[i] -= min;
	data[i] /= tmp;
	data[i] = 1.0 - data[i];
    }
}

/*
    Sobre la pena duermo, sólo y uno,
    pena es mi pan y pena mi batalla.

  	    Miguel Hernández.
*/

/*-------------------------------------------------------------------
    ET_VOILA
	This routine is designed to be called once by the main
    program as a presentation card.

    J. R. Valverde
-------------------------------------------------------------------*/

void et_voila()
{
    printf("\n");
    printf("                            ------\n");
    printf("                            --  --\n");
    printf("                 / / /      ------      \\ \\ \\\n");
    printf("                / /\n");
    printf("               / / /                      \\ \\ \\\n");
    printf("\n");
    printf("              | | |                        | | |\n");
    printf("                  |        P I L A R       | |\n");
    printf("              | | |                        | | |\n");
    printf("\n");
    printf("               \\ \\ \\                      / / /\n");
    printf("                \\                        / / /\n");
    printf("                 \\ \\ \\      --  --      / / /\n");
    printf("                            ------\n");
    printf("                            --  --\n");
    printf("\n");
}

/*---------------------------------------------------------------
    HORROR
	Subroutine to print horror messages.

    INPUT:  S1: string containing calling routine name
	    S2: string containing the horrible error message
    OUTPUT: None.

    J. R. Valverde
---------------------------------------------------------------*/

void horror(s1, s2)
const char *s1, *s2;
{
  fprintf(stderr, "ERROR. HORROR: %s -- %s\n", s1, s2);
}

/*---------------------------------------------------------------
    LETAL
	Print a letal error message and die.

    INPUT:  NAM -- name of the caller
	    MSG -- message

    J. R. Valverde
---------------------------------------------------------------*/

void letal(nam, msg)
char *nam, *msg;
{                             
    printf("\nLETAL ERROR: %s -- %s", nam, msg);
    printf("\nSic transit gloria mundi...\n");

    exit(LETAL_ERROR);
}

/*---------------------------------------------------------------
    LETALF
	Print a formated fatal error message and die.

    INPUT:  FORMAT -- Format to use in the error message
	    VARS   -- Optional list of arguments.

    J. R. Valverde
---------------------------------------------------------------*/

static char letal_buffer[128];

void letalf(char *format, ...)
{
    va_list argptr;

/*  fprintf(stderr, "\nLETAL ERROR: "); */
    va_start(argptr, format);
    vsprintf(letal_buffer, format, argptr);
    va_end(argptr);
    fprintf(stderr, "LETAL ERROR: %s\nSic transit gloria mundi...\n",
	letal_buffer);
/*  fprintf(stderr, "\nSic transit gloria mundi...\n"); */

    exit(LETAL_ERROR);
}

/*-------------------------------------------------------------------
    IS_EMPTY
	Return wether a given line is empty (i. e. only contains
    blanks) or not

    J. R. Valverde
--------------------------------------------------------------------*/

bool is_empty(line)
char *line;
{
    while (line && isspace(*line)) line++;
    return (line)? FALSE : TRUE;
}


#if defined(SYSTEM_III) || defined(SYSTEM_V) || defined(XENIX_3)

/*---------------------------------------------------------------
    INDEX
	returns pointer to first occurrence of a
    character in a string.

    INPUT:  S -- String to search
	    C -- Character to find
    OUTPUT: Pointer to first occurrence of C in S
	    NULL if C is not in S

    J. R. Valverde
---------------------------------------------------------------*/

char *index(s, c)
register char *s, c;
{
    do {
    if (*s == c)
        return (s);
    } while (*s++);

    return NULL;
}

/*----------------------------------------------------------------
    RINDEX
     	Find last occurrence of a character in a string.

    INPUT:  S -- String to search
	    C -- Character to find
    OUTPUT: Pointer to last occurrence of C in S
	    NULL if C is not in S

    J. R. Valverde
----------------------------------------------------------------*/

char *rindex(s, c)
register char *s, c;
{
    register char *r;

    r = NULL;
    do {
    if (*s == c)
        r = s;
    } while (*s++);

    return (r);
}

#endif

/*------------------------------------------------------------------
    TMPNAME
	Return a temporal name (new with each different call)

    INPUT:
	NAM:	string to store the name
    OUTPUT:
	A pointer to the string. Or a NULL ptr else.

    J. R. Valverde
-------------------------------------------------------------------*/

char *tmpname(nam)
char *nam;
{
    static char count = 0;

    if (nam == NULL)
	if ( (nam = (char *) malloc(MAXPATHLEN)) == NULL)
	    return NULL;

    strcpy(nam, "PPA");
    sprintf(nam+3, "%x%x.TMP%c", getpid(), count++, '\0');
    return nam;
}

/*----------------------------------------------------------------
    STRINDEX
	Locate first occurrence of a substring inside a bigger
    string and return its position.

    INPUT:  STR	-- String to be searched
	    SUB	-- Substring to be found
    OUTPUT: Position if Sub is in Str,
	    < 0 if sub is not in str.

    J. R. Valverde
--------------------------------------------------------------------*/

int strindex(str, sub)
char *str, *sub;
{
    int lensub, top, i;

    lensub = strlen(sub);
    top = strlen(str) - lensub;

    if (top < 0)
	return -2;	/* substring longer than string */

    for (i = 0; i <= top; i++, str++)
	if (strncmp(str, sub, lensub) == 0)
	    return i;
    return -1;
}

/*----------------------------------------------------------------
    GETMEM
	Allocate memory and test the process. If it cannot allocate
    the memory, the program is aborted with a letal error message.

    INPUT:  SIZE    -- Size of memory to allocate (in bytes)
    OUTPUT: A pointer to the allocated memory block.

    J. R. VALVERDE
---------------------------------------------------------------*/

void *getmem(size)
unsigned int size;
{
    void *ptr;

    ptr = malloc( size );	/* Get cleared memory */
    if (ptr == NULL)
	letal( "GET MEMORY", "not enough memory");
    return (ptr);
}

/*----------------------------------------------------------------
    GETCH
	Get a character from keyboard using buffered input
    without echoing.

    OUTPUT: character obtained or -1 if EOF.

    J. R. VALVERDE
----------------------------------------------------------------*/

#define BUFFERSIZE	100                 

static char buf[BUFFERSIZE];	    /* buffer for ungetch */
static int bufp = 0;		    /* next free position in buffer */

#ifdef getch
#undef getch
#endif

int getch()
{
    return ((bufp > 0) ? buf[--bufp] : wgetch(stdscr));
}

/*---------------------------------------------------------------
    GETCHE
	Get next character from keyboard echoing to screen
    usig buffered input.

    OUTPUT: character read or -1 if EOF

    J. R. VALVERDE
----------------------------------------------------------------*/

int getche()
{
    extern bool initialized;
    int ch;

    ch = getch();
    putchar(ch);
    return (ch);
}

/*----------------------------------------------------------------
    UNGETCH
	Return a character to the keyboard buffer

    INPUT:  C -- Character to restore
    OUTPUT: C -- If success
	    EOF -- if failure

    J. R. VALVERDE
------------------------------------------------------------------*/

int ungetch(c)		    /* restore a character */
int c;
{
    if (bufp > BUFFERSIZE) {
	fprintf(stderr, "\007");
	return EOF;
    }
    else {
	buf[bufp++] = c;
	return c;
    }
}
                     
/*--------------------------------------------------------------*/
/*	   	     	    DEBUGGING				*/
/*--------------------------------------------------------------*/

#if DEBUG == 4

int iarr[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
float farr[10] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
double darr[10] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };

int imax, imin;
float fmax, fmin;
double dmax, dmin;

counter n;

float ave, adev, sdev, var, skew, kurt;

main()
  {
    int i;

    et_voila();

    n = 10;
    iMaxMin(iarr, n, &imax, &imin);
    fMaxMin(farr, n, &fmax, &fmin);
    dMaxMin(darr, n, &dmax, &dmin);
    printf("\nRESULTS:\n");
    printf("INT    -- MAX = %4d - MIN = %3d\n", imax, imin);
    printf("FLOAT  -- MAX = %2.1f - MIN = %2.1f\n", fmax, fmin);
    printf("DOUBLE -- MAX = %2.1f - MIN = %2.1f\n", dmax, dmin);
    moment(farr, n, &ave, &adev, &sdev, &var, &skew, &kurt);
    printf("\n\nAVERAGE = %g\nAVERAGE DEVIATION = %g\n", ave, adev);
    printf("\nSTANDARD DEVIATION = %g\nVARIANCE = %g\n", sdev, var);
    printf("\nSKEWNESS = %g\nKURTOSIS = %g\n", skew, kurt);
    normalize(farr, n);
    for (i = 0; i < 10; i++)
      printf("%g\t", farr[i]);
    printf("\n\n");
    antinormalize(farr, n);
    for (i = 0; i < 10; i++)
      printf("%g\t", farr[i]);

    letalf("\nGoodbye %s world\n", "cruel");

  }

#endif
