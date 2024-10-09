/**
 * \file P_graphic.h
 * \brief Auxiliary definitions for graphics
 *
 *  \author J. R. Valverde
 */
#define _P_GRAPHIC_H_

/*----------------------------------------------------------------*/

#define ESC	'\033'
/** \def ESC
    escape sequences --> ESC seq ASCII
  	    only uses 7 bit characters
	    can be used in 7 or 8 bit environments. */
#define _ESC	0x1B

#define CSI	'\233'
/** \def CSI
    Control sequences --> CSI graphic ASCII seq
	    It can be also represented in 7 bits with code
	    extension: ESC [
	    It can only be used in an 8 bit environment, but
	    is preferable whenever possible.	*/
#define _CSI	0x9B


/* \def DCS__ST	
    Device control strings --> DCS ctrl str ST
	Only for 8 bit environments.
	In 7 bit environments it should be substituted by
	    DCS --> ESC P
	    ST  --> ESC /	    */
#define DCS	'\220'
#define _DCS	0x90
#define ST	'\234'
#define _ST	0x9C
                                          
/*-----------------------------------------------------------------*/

#define NP_ALLOW    1
#define NP_DISABLE  0

#define MAX_X	    799
#define MAX_Y	    479

typedef struct {
    int x;
    int y;
    } Point_t;

/*------------------------------------------------------------------

	    ReGis Graphics Protocol
	       Controls Mode

------------------------------------------------------------------*/

typedef enum {
    RG_PREV_COMMAND,
    RG_HIGH_COMMAND,
    RG_PREV_DISPLAY,
    RG_HIGH_DISPLAY,
    } ReGisMode_t;

