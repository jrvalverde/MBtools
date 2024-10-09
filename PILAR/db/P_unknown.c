/*
 *  P_UNKNOWN.C
 *	This module contains the function needed to get
 *  a protein from an Unknown-format protein database.
 *
 *  Synopsis:
 *	This module provides the function
 *	    status Unknown_next_protein(db, prot)
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
#define _P_UNKNOWN_

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
    UNKNOWN_NEXT_PROTEIN
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

status Unknown_next_protein(db, prot)
db_t db;
protein_t *prot;
{
    int seqlen, i;
    char *aa;
    flag nonseq;
    static char my_buff[MAXLINE];

    prot->p_index.i_pos = db->db_pos;

    memset(prot->p_name, '\0', PNAME_LEN);
    strncpy(prot->p_name, db->db_buffer, PNAME_LEN - 1);

    /* skip anomalous lines until we find a valid one */
    seqlen = 0, nonseq = FALSE;
    db->db_pos = ftell(db->db_fp);
    while ( fgets(db->db_buffer, MAXLINE, db->db_fp) != NULL) {
      /* while not eof */

      /* if there is a nonvalid char, suppose this is not
      a sequence line: first save growing seq in a temp
      buffer. If it is not correct, then discard it, else
      add the buffer to the sequence */
      for (aa = db->db_buffer, i=0; *aa != '\0'; aa++)
	if ((aa_table[ ASCII(*aa) ]) > 0) {
            /* if it's an aminoacid add it to the sequence */
	    my_buff[i] = ASCII(*aa);
	    i++;
	}
	else {
	    if ((isspace(*aa)) || (isdigit(*aa)))
		continue;
	    else {
		nonseq = TRUE;
		break;
	    }
	}
      
      if (!nonseq) {		/* save this first piece of sequence */
	strcpy(prot->p_seq, my_buff);
	seqlen = i;
	continue;		/* jump to sequence read loop */
      }
    }


    /* get sequence until we find end of file or another new entry */
    db->db_pos = ftell(db->db_fp);
    while ( fgets(db->db_buffer, MAXLINE, db->db_fp) != NULL) {

      for (aa = db->db_buffer, i=0; *aa != '\0'; aa++)
	if ((aa_table[ ASCII(*aa) ]) > 0) {
            /* if it's an aminoacid add it to the sequence */
	    my_buff[i] = ASCII(*aa);
	    i++;
	}
	else {	/* ignore spaces and digits */
	    if ((isspace(*aa)) || (isdigit(*aa)))
		continue;
	    else {  /* The line is not a seq line. */
		nonseq = TRUE;
		break;
	    }
	}

      if (nonseq) break;
      strcpy(&(prot->p_seq[seqlen]), my_buff);
      seqlen += i;

      db->db_pos = ftell(db->db_fp);
    }

    prot->p_seq[seqlen] = '\0';
    prot->p_seqlen = seqlen;
    if (prot->p_seqlen == 0)	    /* if no sequence has been obtained */
      return FALSE;

    return TRUE;
}                                  

/*-----------------------------------------------------------------*/

