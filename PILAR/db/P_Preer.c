/*
 *  P_PREER.C
 *	This module should contain the function needed to get
 *  a protein from a Preer-format protein database. It is NOT
 *  currently implemented, and the function just returns with
 *  a status of FALSE (i. e. no protein read).
 *
 *  Synopsis:
 *	This module provides the function
 *	    status Preer_next_protein(db, prot)
 *	    db_t *db;
 *	    protein_t *prot;
 *
 *	which is to be installed by the database initialization
 *  process in module P_DBANK.C in the database description
 *  sructure. The function will be called by reference thereafter
 *  as a field in the db structure.
 *
 *  Notes:
 *	The function is not implemented, it ignores its arguments,
 *  does nothing, and returns FALSE.
 *
 *  Files:
 *	Needs the following files:
 *	    PORTABLE.H
 *	    P_TYPES.H
 *	    P_DEFS.H
 *	    P_EXTERN.H
 *	    P_DBANK.H
 *
 *  History:
 *	08 - apr - 1990	Stub version.
 *	01 - aug - 1991 Ultrix port.
 *
 *  Designed by:
 *	José Ramón Valverde Carrillo.
 *	JRamon@emduam51.bitnet
 *
 *  Copyright:
 *	I hereby place this source code in the Public Domain.
 *  Code is provided "as is", with no warranties of any type.
 *  Commercial use without prior consent of the author is
 *  absolutely prohibited.
 */

#define _P_DBANK_
#define _P_PREER_

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_extern.h"
#include "P_dblow.h"


/*-------------------------------------------------------------------
    PREER_NEXT_PROTEIN
	Gets next protein from the database. It sequentially access the
    database and reads the protein name and sequence. Also it gets
    the index of the protein into the database and obtains the access
    keywords.

    INPUT:  PROT -- Structure of type protein to store the read sequence
	    DB	-- Database to get the protein from
    OUTPUT: TRUE -- If it could read a new protein
	    FALSE -- otherwise

    J. R. Valverde
 --------------------------------------------------------------------*/

status Preer_next_protein(db, prot)
db_t db;
protein_t *prot;
{
    return FALSE;
}

/*-----------------------------------------------------------------*/
