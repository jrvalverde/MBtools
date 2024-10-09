/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE:							*
*	    P_GRAPHIC4.C					*
*								*
*	SYNOPSIS:						*
*	    Curve commands.					*
*								*   
*	DESCRIPTION:						*
*	    Circle()						*
*	    RCircle()						*
*	    Arc()						*
*	    RArc()						*
*								*
*  	CAVEATS:						*
*								*
*	NOTES:							*
*								*
*	FILES:							*
*	    STDIO.H						*
*	    P_GRAPHIC.H						*
*	    P_GRAPHIC.X						*
*								*
*	SEE ALSO:						*
*	    VT240 programmers pocket guide. (Digital Equipment	*
*	Corporation).						*
*								*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo.			*
*								*
*  	HISTORY:						*
*								*
*  	COPYRIGHT:						*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>

#include "P_graphic.h"

/*--------------------------------------------------------------------
    Circle
	Circle with center at current position. X, Y defines
    a point on the circumference of the circle. The X, Y value
    can be absolute, relative, or absolute-relative.

    INPUT:  X, Y Coordinates of the point in the circumference

    J. R. Valverde
--------------------------------------------------------------------*/

void Circle(x, y)
int x, y;
{
    printf("C[%d,%d]", x, y);
}

/*--------------------------------------------------------------------
    RCircle
	Circle with center at specified position. X, Y defines
    a point to serve as the circle center, while current active
    position defines a point on the circumference. The X, Y values
    can be absolute, relative, or absolute/relative

    J. R. Valverde
--------------------------------------------------------------------*/

void RCircle(x, y)
int x, y;
{
    printf("C(C)[%d,%d]", x, y);
}

/*--------------------------------------------------------------------
    Arc
	Arc with center at current position. X, Y defines the
    starting point for drawing the arc. The signed value of the
    degrees (+ or -) determines wich direction the arc is drawn from
    that point: + for counterclockwise, and - for clockwise.
    The X, Y value can be absolute, relative or absolute/relative

    J. R. Valverde
--------------------------------------------------------------------*/

void Arc(x, y, degrees)
int x, y, degrees;
{
    printf("C(A%+d)[%d,%d]", degrees, x, y);
}

/*--------------------------------------------------------------------
    RArc
	Arc with center at specified position. X, Y defines
    the center, while the current active position is the point from
    which the arc is drawn. The signed value of degrees (+ or -)
    determines which direction the arc is drawn: + for counter-
    clockwise and - for clockwise. The X, Y value can be absolute,
    relative, or absolute/relative.

    J. R. Valverde
--------------------------------------------------------------------*/

void RArc(x, y, degrees)
int x, y, degrees;
{
    printf("C(A%+dC)[%d,%d]", degrees, x, y);
}

/*--------------------------------------------------------------------
    ClosedCurve
	Closed curve sequence. Defines a closed curve graphic
    image built from interpolation of X, Y positions specified 
    within the option. The X, Y values can be absolute, relative
    or absolute/relative

    J. R. Valverde
--------------------------------------------------------------------*/

