/*
 *  P_PARSER.H
 *
 *	Header file with the definitions pertaining to module
 *  P_Parser.C. It is usually #included via P_Extern.H
 *
 *	Designed by:
 *	    J. R. Valverde	8  - apr - 1990
 */


extern void parse(FILE *);

extern void store_token(flag, char *);

extern char *retrieve_token(flag *);

/*
 *  What else did you expect ?
 */

