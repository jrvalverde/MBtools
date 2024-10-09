/*
 *  P_PIR.C
 *	This module contains the function needed to get
 *  a protein from a PIR-format protein database.
 *
 *  Synopsis:
 *	This module provides the function
 *	    status PIR_next_protein(db, prot)
 *	    db_t *db;
 *	    protein_t *prot;
 *
 *	which is to be installed by the database initialization
 *  process in module P_DBANK.C in the database description
 *  sructure. The function will be called by reference thereafter
 *  as a field in the db structure.
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
#define _P_PIR_

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_extern.h"
#include "P_dblow.h"

static bool next_entry(db_t);

/*-------------------------------------------------------------------
    PIR_NEXT_PROTEIN
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

status PIR_next_protein(db, prot)
db_t db;
protein_t *prot;
{
    int seqlen;
    char *aa;

    if (! next_entry(db))
      return FALSE;

    prot->p_index.i_pos = db->db_pos;

    /* db_buffer contains the entry line */
    memset(prot->p_name, '\0', PNAME_LEN);
    strncpy(prot->p_name, &db->db_buffer[1], PNAME_LEN - 1);

    /* skip junk line */
    db->db_pos = ftell(db->db_fp);
    if (fgets(db->db_buffer, MAXLINE, db->db_fp) == NULL)
	return FALSE;	    /* EOF */

    /* get sequence until we find end of file or another new entry,
    ignoring comment lines beginning with a ';' */
    seqlen = 0;
    db->db_pos = ftell(db->db_fp);
    while ( fgets(db->db_buffer, MAXLINE, db->db_fp) != NULL) {
      /* while not eof */

      if (db->db_buffer[0] == '>')
	/* new entry: exit loop */
	break;

      for (aa = db->db_buffer; *aa != '\0'; aa++)
	if ((aa_table[ ASCII(*aa) ]) > 0) {
          /* if it's an aminoacid add it to the sequence */
          prot->p_seq[seqlen] = ASCII(*aa);
          seqlen++;                  
      	}

      db->db_pos = ftell(db->db_fp);
    }

    prot->p_seq[seqlen] = '\0';
    prot->p_seqlen = seqlen;
    if (prot->p_seqlen == 0)	    /* if no sequence has been obtained */
      return FALSE;

    return TRUE;
}                                  

/*-----------------------------------------------------------------*/

static bool next_entry(db)
db_t db;
/*
 * Search for next entry into the database.
 *    It can be recognized because the line must begin with a '>'
 *   character.
 *    Also actualize file pointer as needed.
 */
  {
    while (db->db_buffer[0] != '>') {
      db->db_pos = ftell(db->db_fp);
      if (fgets(db->db_buffer, MAXLINE, db->db_fp) == NULL) /* eof */
	return FALSE;
    }
    return TRUE;
  }
