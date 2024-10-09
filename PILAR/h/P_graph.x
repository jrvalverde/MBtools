/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE							*
*	    P_GRAPHIC.X -- external definitions for PILAR	*
*	graphic routines.					*
*								*
*	SYNOPSIS						*
*								*
*	DESCRIPTION						*
*	    Contains all the external definitions for graphics	*
*	(VT240) modules of the Protein Identification and	*
*	Library Access Resource. Include it for making all	*
*	the needed functions globally available.		*
*								*
*	CAVEATS							*
*	    This is dynamically modified as the project is	*
*	being built up.						*
*								*
*	FILES							*
*	    P_GRAPHIC.H						*
*			       					*
*	NOTES							*
*	    Some definitions may be never be used. They are 	*
*	there just for historical reasons.			*
*								*
*	SEE ALSO						*
*	    It may be worth having a look to the mentioned	*
*	files to get more information about these routines.	*
*								*
*	DESIGNED BY						*
*	    José Ramón Valverde Carrillo.			*
*								*
*	HISTORY							*
*	    7 - oct - 1989  J. R. Valverde			*
*								*
*	COPYRIGHT						*
*	    As far as I know this is supposed to be public	*
*	domain software. Please, do not use it for commercial	*
*	purposes without the written permission of the author.	*
*	(I could become very angry if I became very hungry just	*
*	because of you).					*
*						    J. R.	*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#ifndef _P_GRAPHIC_H_
#include "P_graphic.h"
#endif

/*
 *	P_GRAPHIC1.C
 */

#define InReGis(mode)	printf("%c%dp", DCS, mode)

#define OutReGis    putchar(ST);

extern void enterReGis(ReGisMode_t);

extern void exitReGis();

/*
 *	P_GRAPHIC2.C
 */

extern void SetPVFactor(int);

extern void SelStdPattern(int);

extern void SelBinPattern(char *);

extern void SelPatternFactor(int);

extern void EnableRevPattern(int);

extern void EnableBitMap(int);

extern void SelFgIntensity(int);

extern void SelFgHLS(int, int, int);

extern void SelMapHLS(int, int, int, int);

extern void SelWriteMode(char);

extern void EnableShading(int);

extern void SelShadeYRef(int);

extern void SelShadeXRef(int);

extern void SelChShade(char);

extern void EraseScreen();

/*
 *	P_GRAPHIC3.C
 */

extern void Point();

extern void MoveTo(int, int);

extern void LineTo(int, int);

/*
 *	P_GRAPHIC4.C
 */

extern void Circle(int, int);

extern void RCircle(int, int);

extern void Arc(int, int, int);

extern void RArc(int, int, int);

/*
 *	    The End...
 */
