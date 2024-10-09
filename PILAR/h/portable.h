/*---------------------------------------------------------*
 * PORTABLE.H -- type and macro definitions for portable C *
 *                                                         *
 * PORTABLE C AND UNIX SYSTEM PROGRAMMING.		   *
 * J E Lapin                                               *
 * Prentice-Hall, Englewood Hills, 1987			   *
 * 							   *
 * Enhanced august the 24th, 1988 by J R Valverde	   *
 *---------------------------------------------------------*/

#ifdef MEM_DEBUG
#include "mshell.h"
#endif

/*----------------------------------------------------------------
    Standard types
----------------------------------------------------------------*/

#ifdef bool
#undef bool
#endif

typedef int	bool;		/* >= 16 bits used as boolean	*/
typedef int	boolean;
typedef char	flag;		/* >= 8 bits used as boolean	*/

typedef unsigned short word;	/* >= 16 bit  unsigned type	*/

typedef int	rchar;		/* I/O function return values	*/
typedef char   *string;		/* for readability		*/

typedef int	status;		/* function status on return	*/

/* If void is not supported, uncomment the following line.	*/
/* typedef int	void;	*/
/* ------>    #define	DEFAULT	int        Default size		*/


/*--------------------------------------------------------------
   Use this length and type for filename variables -- use
   256 on 4.2BSD since it is defined as such in dir.h	
--------------------------------------------------------------*/

#ifndef MAXPATHLEN
#define MAXPATHLEN  64
#endif

typedef char	path[MAXPATHLEN];


/*--------------------------------------------------------------
   Define byte and BYTE to portably support an unsigned
   8 bit data type (see Section 2.9.1 in text).		
--------------------------------------------------------------*/

/*	If char is signed by default and unsigned keyword is
	not allowed with char, byte and BYTE are:
typedef char	  byte;				* 8 bit unsigned type	* /
#define BYTE(x)	  (byte) ((x) & 0xff)		* BYTE truncates data	*/

/*	else they are:	*/
typedef unsigned char byte;
#define BYTE(x)	  (byte) ((x) & 0xff)


/*---------------------------------------------------------------
   Standard macros
---------------------------------------------------------------*/

#define HIBYTE(x)	(((x) >> 8) & 0xff)	/* hi byte of short	*/
#define LOBYTE(x)	((x) & 0xff)		/* lo byte of short	*/
#define HIWORD(x)	(((x) >> 16) & 0xffffL)	/* upper half of long	*/
#define LOWORD(x)	((x) & 0xffffL)		/* lower half of long	*/
#define CHAR(x)		((x) & 0x7F)		/* Truncate to 7 bits	*/
#define ASCII(x)	((x) & 0x7F)		/* Convert to ascii	*/
#define WORD(x)		((x) & 0xffffL)		/* Truncate to 16 bits	*/
#define DECODE(x)	((int) ((x) - '0'))	/* int value of a digit	*/

/* long unsigned right shift */
#define LURSHIFT(n, b) (((long) (n) >> (b)) & (0x7fffffffL >> ((b) - 1)))

/* number of elements in array a */
#define DIM(a)		(sizeof(a)/sizeof(*(a)))

/*----------------------------------------------------------------
   Macros for bit manipulation
----------------------------------------------------------------*/

#define BIT_0	0x01
#define BIT_1	0x02
#define BIT_2	0x04
#define BIT_3	0x08
#define BIT_4	0x10
#define BIT_5	0x20
#define BIT_6	0x40
#define BIT_7	0x80
#define BIT_8	0x100
#define BIT_9	0x200
#define BIT_10	0x400
#define BIT_11	0x800
#define BIT_12	0x1000
#define BIT_13	0x2000
#define BIT_14	0x4000
#define BIT_15	0x8000
#define BIT_16	0x10000
#define BIT_17	0x20000
#define BIT_18	0x40000
#define BIT_19	0x80000
#define BIT_20	0x100000
#define BIT_21	0x200000
#define BIT_22	0x400000
#define BIT_23	0x800000
#define BIT_24	0x1000000
#define BIT_25	0x2000000
#define BIT_26	0x4000000
#define BIT_27	0x8000000
#define BIT_28	0x10000000
#define BIT_29	0x20000000
#define BIT_30	0x40000000
#define BIT_31	0x80000000

#define BitTest(x, bit)	((x) & (bit))

#define BitSet(x, bit)  ((x) |= (bit))

#define BitClr(x, bit)	((x) &= ~(bit))

#define BitAnd(x, y)	((x) & (y))

#define BitOr(x, y)	((x) | (y))

#define BitXor(x, y)	((x) ^ (y))

#define BitNot(x)	(~(x))

#define BitShift(x, count)	( ((count) > 0)? (x) << (count) : (x) >> (count) )


/*----------------------------------------------------------------
   Scope control pseudo-keywords
----------------------------------------------------------------*/

#define public			/* public is C default scope	*/
#define private	static		/* static really means private	*/
#define	LOCAL	auto		/* Local var on 68000	   	*/
#define	EXTERN	extern		/* External variable	   	*/
#define	MLOCAL	static		/* Local to module	   	*/
#define	GLOBAL			/* Global variable	   	*/


/*----------------------------------------------------------------
   Standard constants
----------------------------------------------------------------*/

#ifndef TRUE
#define TRUE	(1)		/* for use with booleans	*/
#endif

#ifndef FALSE
#define FALSE	(0)
#endif

#define	YES	(1)		/* "TRUE"			*/
#define	NO	(0)		/* "FALSE"			*/

#define SUCCEED	( 0)		/* for use in exit()		*/
#define FAILURE	(-1)		/* for exit() & error returns	*/
#define	FAIL	(-1)		/* Function failure return val  */
#define SUCCESS	( 0)		/* Function success return val  */

#define	EOL	'\n'		/* end of line char		*/
#define EOS	'\0'		/* end of string char		*/
#define EOP	'\14'		/* end of page (FF) char	*/

#ifndef EOF
#define	EOF	(-1)		/* EOF Value			*/
#endif

#define	FOREVER	for(;;)		/* Infinite loop declaration    */
/*   #define	NULL	0	   Null pointer value	    	*/
#define NULLPTR (char *) 0	/* Pointer to NULL		*/


/*----------------------------------------------------------------
   portable.h ends here
----------------------------------------------------------------*/


