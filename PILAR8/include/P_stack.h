/*
 *  P_STACK.H
 *
 *	Include file with the definitions needed for module
 *  P_Stack.C
 *
 *	Requires PORTABLE.H and P_TYPES.H
 *
 *	Designed by:
 *	    J. R. Valverde	    8  - apr - 1990
 */


extern lifo_t new_lifo();

extern boolean push(info_t *, lifo_t);

extern info_t *pop(lifo_t);

