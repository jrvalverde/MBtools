/*
 *  P_UTILS.H
 *
 *	Include file containing the definitions for the module
 *  P_Utils.C
 *
 *	Requires PORTABLE.H
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */


extern void iMaxMin(int *, counter , int *, int *);

extern void fMaxMin(float *, counter , float *, float *);

extern void dMaxMin(double *, counter , double *, double *);

extern void rMaxMin(void **, counter , void **, void **,
	    int (*)());

extern void moment(float *, counter , float *, float *,
	    float *, float *, float *, float *);

extern void KSone(float *, counter, float (*)(),
	    float *, float *);

extern void KStwo(float *, counter,
	  float *, counter,
	  float *, float *);

extern void normalize(float *, counter);

extern void antinormalize(float *, counter);

extern void et_voila();

extern void horror(const char *, const char *);

extern void letal(char *, char *);

extern void letalf(char *, ...);

extern bool is_empty(char *);

#if defined(SYSTEM_III) || defined(SYSTEM_V) || defined(XENIX_3)

extern char *index(char *, char);

extern char *rindex(char *, char);

#endif

extern int strindex(char *, char *);

extern void *getmem(unsigned int);

/* extern int getch(); */

extern int getche();

extern int ungetch(int);

/*
 *  END OF P_UTILS.H
 */
