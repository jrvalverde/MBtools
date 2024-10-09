/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	MODULE:							*
*	    P_GRAPHIC3.C					*
*								*
*	SYNOPSIS:						*
*	    Routines for vector commands.			*
*								*
*  	DESCRIPTION:						*
*	    Point()						*
*	    MoveTo()						*
*	    LineTo()						*
*	    _begBoundSeq()					*
*	    _begUnboundSeq()					*
*	    _endSequence()					*
*								*
*	CAVEATS:						*
*								*
*	NOTES:							*
*								*
*	FILES:							*
*								*
*	SEE ALSO:						*
*	    VT240 series Programmer Pocket Guide		*
*								*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo			*
*								*
*  	HISTORY:						*
*								*
*	COPYRIGHT:						*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#include <stdio.h>

/*--------------------------------------------------------------------
    MoveTo
	Position using X, Y values. X and Y can be absolute,
    relative or absolute/relative

    INPUT:  X, Y -- coordinates

    J. R. Valverde
--------------------------------------------------------------------*/

void MoveTo(x, y)
int x, y;
{
    printf("P[%d,%d]", x, y);
}

/*--------------------------------------------------------------------
    Point
	Draw dot. Used to write to a single pixel defined
    by the current active position. No cursor movement occurs.

    J. R. Valverde
--------------------------------------------------------------------*/

void Point()
{               
    printf("V[]");
}

/*--------------------------------------------------------------------
    LineTo
	Draw line. (X,Y) value defines an end point for a line
    to be drawn from the current active position. The (X, Y) value
    can be absolute, relative, or absolute/relative.

    INPUT:  X	Coord. x of the destination point
	Y   Coord. y of the destination point

    J. R. Valverde
--------------------------------------------------------------------*/

void LineTo(x, y)
int x, y;
{
    printf("V[%d,%d]", x, y);
}

/*--------------------------------------------------------------------
    _begBoundSeq
	begin a bounded sequence. Stores current active position
    for reference at the end of the sequence.

    J. R. Valverde
--------------------------------------------------------------------*/

void _begBoundSeq()
{
    printf("V(B)");
}

/*--------------------------------------------------------------------
    _begUnboundSeq
	Start an unbounded sequence. Stores a dummy position
    for reference at the end of the sequence

    J. R. Valverde
--------------------------------------------------------------------*/

void _begUnboundSeq()
{
    printf("V(S)");
}

/*--------------------------------------------------------------------
    _endSequence
	End of sequence. References last stored (Bounded) or 
    (unbounded) value. If value referenced was stored by Bound, a
    line is drawn from the active position where End is sensed,
    to the location stored by (bound). If value referenced was
    stored by (Unbound), no line is drawn, and active position
    remains at current position.

    J. R. Valverde
--------------------------------------------------------------------*/

void _endSequence()
{
    printf("V(E)");
}
                                                                    
