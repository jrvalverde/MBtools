/*
 *  P_UWGCG.C
 *	This module contains the function needed to get
 *  a protein from an University of Wisconsin / Genetics
 *  Computer Group formatted protein file.
 *
 *  Synopsis:
 *	This module provides the function
 *	    status UWGCG_next_protein(db, prot)
 *	    db_t *db;
 *	    protein_t *prot;
 *
 *	which is to be installed by the database initialization
 *  process in module P_DBANK.C in the database description
 *  sructure. The function will be called by reference thereafter
 *  as a field in the db structure.
 *
 *  Notes:
 *	UW/GCG format files contain only one sequence per file.
 *	The routine expects that db->db_buffer is already filled
 *	with the last comment line of the file, which is:
 *
 *	SEQ_NAME    Length: #	Date	Time	Check: # ..
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
 *	08 - apr - 1990	First version. Routine is extracted from
 *		its former location at P_DBANK.C and renamed.
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
#define _P_UWGCG_

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_extern.h"
#include "P_dblow.h"


/*-------------------------------------------------------------------
    UWGCG_NEXT_PROTEIN
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

status UWGCG_next_protein(db, prot)
db_t db;
protein_t *prot;
{
    int seqlen;
    char *aa;
    flag done = FALSE;

    /* On entry we should be either at EOF if the routine has already been
    called once, or at the last comment line, which contains seq.
    name, length, date/time, CRC and the last ".." */

    /* save current position */
    prot->p_index.i_pos = db->db_pos;

    memset(prot->p_name, '\0', PNAME_LEN);
    strncpy(prot->p_name, &db->db_buffer[1], PNAME_LEN - 1);

    /* get sequence until we find end of file */
    seqlen = 0;
    db->db_pos = ftell(db->db_fp);
    while ( fgets(db->db_buffer, MAXLINE, db->db_fp) != NULL) {
      /* while not eof */

      /* Just a precaution. The user may have concatenated several
      output files. This could lead to an error. Alternatively we
      could process this event. As for now this feature is disabled. */
      if (strindex(db->db_buffer, "ID  ") == 0)
	break;

      for (aa = db->db_buffer; *aa != '\0'; aa++) {
	if ((aa_table[ ASCII(*aa) ]) > 0) {
          /* if it's an aminoacid add it to the sequence */
          prot->p_seq[seqlen] = ASCII(*aa);
          seqlen++;
      	}
	if (seqlen > MAX_SEQ_LEN) {
	    done = TRUE;
	    break;
	}
      }
      db->db_pos = ftell(db->db_fp);

    }

    prot->p_seq[seqlen] = '\0';
    prot->p_seqlen = seqlen;
    if (prot->p_seqlen == 0)	    /* if no sequence has been obtained */
      return FALSE;

    return TRUE;
}                                  

/*
    Elle dormait: son doigt tremblait, sans améthyste
    Et nu, sous sa chemise: après un soupir triste,
    Il s'arrêta, levant au nombril la batiste.

    Stéphane Mallarmé. "...Mysticis umbraculis".
*/
