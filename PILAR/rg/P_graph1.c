/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE:							*
*	    P_GRAPHIC1.C					*
*								*
*  	SYNOPSIS:						*
*	    This module contains the routines for entering	*
*	and exiting ReGis graphics mode.			*
*      								*
*	DESCRIPTION:						*
*	    InReGis()						*
*	    enterReGis()					*
*	    OutReGis()						*
*	    exitReGis()						*
*								*
*	CAVEATS:						*
*								*
*	NOTES:							*
*	    InReGis and OutReGis are implemented as macros	*
*								*
*	FILES:							*
*								*
*	SEE ALSO:						*
*	    VT240 Series. Programmer Pocket Guide.		*
*								*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo.			*
*								*
*	HISTORY:						*
*	    22 - jul - 1989 J.R.Valverde (implementation)	*
*								*
*	COPYRIGHT:						*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>

#include "P_graphic.h"

/*-------------------------------------------------------------------
    INREGIS
    ENTERREGIS
	Function and macro to enter ReGis mode. This mode is
    available only through VT200 and VT100 modes only.

    INPUT:  MODE    RG_PREV_COMMAND Enter ReGis at previous command
	    level
	    RG_HIGH_COMMAND Causes vt240 to enter ReGis at
	    highest command level
	    RG_PREV_DISPLAY Enter at previous command level
	    with commands displayed on the screen's bottom
	    line
	    RG_HIGH_DISPLAY Enter at highest command level
	    with commands displayed on the screen's bottom
	    line.

    J. R. Valverde
--------------------------------------------------------------------*/

#define InReGis(mode)	printf("%c%dp", DCS, mode)

void enterReGis(mode)
ReGisMode_t mode;
{
    printf("%c%dp", DCS, mode);                                   
}

/*--------------------------------------------------------------------
    OUTREGIS
    EXITREGIS
	Function and macro that terminate a ReGis device
    control string. Exits ReGis mode and returns to text mode.

    J. R. Valverde
--------------------------------------------------------------------*/

#define OutReGis    putchar(ST);

void exitReGis()
{
    putchar(ST);
}
