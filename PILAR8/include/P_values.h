/*
 *  P_VALUES.H
 *
 *	Include file containing the definitions for the module
 *  P_Values.C
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */


extern void assign_values(protein_t *, descriptor_t *);

extern void window_soften(float *, int, int);

extern void moment_soften(float *, int, int, float);

extern void AVmoment_soften(float *, int, int, float);

