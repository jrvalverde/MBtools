/*          
 *  P_GCGSW.C
 *	This module contains the function needed to get
 *  a protein from a GCG/SwissProt-format protein database.
 *
 *  Synopsis:
 *	This module provides the function
 *	    status GCGSW_next_protein(db, prot)
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
 *	28 - jun - 1991	First version. Routine is extracted from
 *		its former location at MPQ.C, adapted, and renamed.
 *	01 - aug - 1991 Ultrix port.
 *
 *  Designed by:
 *	José Ramón Valverde Carrillo.
 *	JRamon@uamed.uam.es
 *
 *  Copyright:
 *	I hereby place this source code in the Public Domain.
 *  Code is provided "as is", with no warranties of any type.
 *  Commercial use without prior consent of the author is
 *  absolutely forbidden.
 */

#define _P_DBANK_
#define _P_EMBL_

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"
#include "P_extern.h"
#include "P_dblow.h"

status GCGSW_next_protein(db, prot)
db_t db;
protein_t *prot;
{
    int seqlen;
    char *aa;
    static bool next_entry();

    if (! next_entry())
      return FALSE;

    prot->p_index.i_pos = db->db_pos;

    /* db_buffer contains now first line: ">>>> access code" */
    strncpy(prot->p_index.i_key, &db->db_buffer[4], KEYLEN);
    /* get protein name explicited in next line */
    memset(prot->p_name, '\0', PNAME_LEN);

    db->db_pos = ftell(db->db_fp);
    if (fgets(db->db_buffer, MAXLINE, db->db_fp) == NULL)
	return FALSE;
    strncpy(prot->p_name, &db->db_buffer[6], PNAME_LEN - 1);

    /* Get sequence */
    seqlen = 0;
    db->db_pos = ftell(db->db_fp);
    if (fgets(db->db_buffer, MAXLINE, db->db_fp) == NULL)
	return FALSE;	    /* EOF */

    do {
	if (strindex(db->db_buffer, ">>>>") == 0)
	    break;	    /* end of sequence */

	for (aa = db->db_buffer; *aa != '\0'; aa++)
	    if ((aa_table[ ASCII(*aa) ]) > 0) {
            /* if it's an aminoacid add it to the sequence */
          	prot->p_seq[seqlen] = ASCII(*aa);
          	seqlen++;                  
      	    }
	db->db_pos = ftell(db->db_fp);
    }	    /* while not eof */
    while ( fgets(db->db_buffer, MAXLINE, db->db_fp) != NULL );

    prot->p_seq[seqlen] = '\0';
    prot->p_seqlen = seqlen;
    if (prot->p_seqlen == 0)
	return FALSE;

    return TRUE;
}                                  

/*-----------------------------------------------------------------*/

static bool next_entry(db)
db_t db;
/*
 * Search for next entry into the database.
 *    It can be recognized because the line must begin with a '>>>>'
 *   sequence.
 *    Also actualize file pointer as needed.
 */
  {
    while (strindex(db->db_buffer, ">>>>") != 0) {
      db->db_pos = ftell(db->db_fp);
      if (fgets(db->db_buffer, MAXLINE, db->db_fp) == NULL) /* eof */
	return FALSE;
    }
    return TRUE;
  }

