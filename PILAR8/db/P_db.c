/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*								*
*	P_DB.C	       						*
*								*
*	DESCRIPTION:						*
*	    Get sequences from a protein database.		*
*								*
*	SYNOPSIS:						*
*	This module provides functions for database access:	*
*		db_open()					*
*		db_next_protein()				*
*		db_seek_protein()				*
*		db_close()					*
*		db_name()					*
*								*
*		idx_open()					*
*		get_index()					*
*		idx_close()					*
*		idx_name()					*
*								*
*	CAVEATS:			 			*
*								*
*	NOTES:							*
*	    Unfinished yet.					*
*								*
*	FILES:							*
*	    This module uses the following files:		*
*		PORTABLE.H					*
*		P_DEFS.H					*
*		P_TYPES.H					*
*		P_EXTERN.H					*
*								*
*	HISTORY:						*
*	    16 - jun - 1989 First version. (IBM-AT).		*
*	    20 - jun - 1989 Ported to VAX-VMS.			*
*	    12 - jul - 1989 Debugged.				*
*	     4 - aug - 1989 Use aa_table to get AAs.		*
*	    16 - sep - 1989 Functions renamed.			*
*	     1 - jul - 1991 Added support for GCG-SwissProt	*
*	     1 - aug - 1991 Ported to ULTRIX			*
*			          				*
*	DESIGNED BY:						*
*	    José Ramón Valverde Carrillo.			*
*								*
*	COPYRIGHT:						*
*	    (C) YoEgo.	Please, would you be so kind		*
*	not to use this code without express permission ?	*
*		Thanks.        					*
*						YoEgo.		*
*								*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#define _P_DBANK_

/*-----------------------------------------------------------------
	MACRO FOR DEBUGGING
-----------------------------------------------------------------*/

#ifndef DEBUG                                                   
/* #define DEBUG    1 */
#define DEBUG	    -1
#endif

/*-------------------------------------------------------------------
	    INCLUDE FILES
-------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#ifdef vms
#include <file.h>
#endif

#if DEBUG == 1
#include <time.h>
#include "P_global.h"	    /* only #included in MAIN modules */
#ifdef __TURBOC__
#include <graphics.h>
#endif
#endif                              

#include "portable.h"
#include "P_defs.h"
#include "P_types.h"

#include "P_utils.h"

#ifndef SEEK_SET
#define SEEK_SET	0
#endif

extern char *getenv(char *name);

/*--------------------------------------------------------------------
	    LOCAL TYPES, EXTERN FUNCTIONS AND DEFINITIONS.
--------------------------------------------------------------------*/

typedef struct db_st{
    FILE    *db_fp;
    int	    db_fd;
    char    db_nam[MAXPATH];
    char    db_buffer[MAXLINE];
    int	    db_linelen;
    long    db_pos;
    int	    db_type;
    status  (*db_next_protein)(); /* (struct db_st *, protein_t *) */
    }	*db_t;

/* these extern definitions are pertinent only for the database
    management module, so they are not put on the global file
    P_EXTERN.H, but only here */

extern status Unknown_next_protein(db_t, protein_t *);
extern status IG_next_protein(db_t, protein_t *);
extern status Strider_next_protein(db_t, protein_t *);
extern status GenBank_next_protein(db_t, protein_t *);
extern status PIR_next_protein(db_t, protein_t *);
extern status Pearson_next_protein(db_t, protein_t *);
extern status EMBL_next_protein(db_t, protein_t *);
extern status UWGCG_next_protein(db_t, protein_t *);
extern status Zuker_next_protein(db_t, protein_t *);
extern status Fitch_next_protein(db_t, protein_t *);
extern status Preer_next_protein(db_t, protein_t *);
extern status GCGSW_next_protein(db_t, protein_t *);

#define DB_UNKNOWN	    0
#define DB_STANFORD	    1
#define DB_INTELLIGENETICS  2
#define DB_STRIDER	    3
#define DB_GENBANK	    4
#define DB_NBRF		    5
#define DB_PIR		    6
#define DB_PEARSON	    7
#define DB_EMBL		    8
#define DB_SWISSPROT	    9
#define DB_UWGCG	    10
#define DB_ZUKER	    11
#define DB_FITCH	    12
#define DB_PREER	    13
#define DB_GCGSW	    14


/*--------------------------------------------------------------------
	DATA BASE MANIPULATION VARIABLES
--------------------------------------------------------------------*/

static const char DB_OPEN[] = "OPEN DATABASE";
static const char CANT_OPEN[] = "CANNOT OPEN PROTEIN SEQUENCES FILE";

/*-------------------------------------------------------------------
	INDEX FILE MANIPULATION VARIABLES
-------------------------------------------------------------------*/

static FILE *idx_fp = NULL;
static int  idx_fd  = 0;
static char idx_nam[MAXPATH + 1] = "";
static long idx_pos = 0;
static int  idx_size;

/*-------------------------------------------------------------------
    OPEN DATABASE
	Opens a protein database.
	It first tries to open the file. If succeeded then it reads
    the first line in the database and sets the pointer into
    the file appropriately. Also saves the current database name
    for possible use later.

    INPUT:  DB_FILNAM -- the database file name (db_filnam).
	    DB_TYPE -- the type of database (if known).
			0 i unknown.
    	    OUTPUT:  -- pointer to a database structure.

    J. R. Valverde
--------------------------------------------------------------------*/

db_t db_open(db_filnam, db_type)
char *db_filnam;
int db_type;
  {
#if (DEBUG == 1)
    static db_t bank;
    static int i;
    static bool is_UWGCG, done;
#else
    db_t bank;
    int i; char *auxnam;
    bool is_UWGCG, done;
#endif

    /* Allocate space for database struct */
    bank = (db_t) NEW( struct db_st );
    if (bank == NULL) return NULL;
    bank->db_fp = NULL;
    bank->db_fd = -1;			/* invalid descriptor */
    bank->db_nam[0] = '\0';
    bank->db_buffer[0] = '\0';
    bank->db_linelen = 80;		/* Default value */
    bank->db_pos = 0L;
    bank->db_type = db_type;

    /* Try to open data base file */
    if ((bank->db_fp = fopen(db_filnam, "rw")) == NULL) {
	/* Is db_filnam an alias?			 */
	/* CAUTION: does not cover case where db_filnam  */
	/* is an alias of an alias of... as is common in */
	/* UW-GCG. To solve this a for loop is needed.   */
	if ( (auxnam = getenv(db_filnam)) == NULL) {
	    horror(DB_OPEN, CANT_OPEN);
	    return NULL;
	}
	else {
	    if ((bank->db_fp = fopen(auxnam, "rw")) == NULL) {
 		horror(DB_OPEN, CANT_OPEN);
		return NULL;
	    }
	}
    }
/*	ALTERNATIVE OPEN ALGORITHM (to be checked):
    while ( (auxnam = getenv(db_filnam)) != NULL)
	dn_filnam = auxnam;
    if ((bank->db_fp = fopen(db_filnam, "rw")) == NULL) {
	horror(DB_OPEN, CANT_OPEN);
	return NULL;
    }
*/

    /* Save protein data base file name */
    strncpy(bank->db_nam, db_filnam, MAXPATH);

    /* set pointer to the beginning of the data base */    
    rewind(bank->db_fp);

    if (db_type == DB_UNKNOWN) {
	/* Check for UW-GCG format since it can masquerade as any
	other format */
	i = 0;
	do {
	    /* Get current position of file pointer */
	    /* First sequence begins always at offset 0L */
	    bank->db_pos = ftell(bank->db_fp);

	    /* We can have anything as a comment until a line
	    that ends as "Check: # .." */
	    fgets(bank->db_buffer, MAXLINE, bank->db_fp);
	    i++;		/* Count # of lines read */
	    if ((strindex(bank->db_buffer, "..") > 0) &&
		(strindex(bank->db_buffer, "Check:") >= 0))
		is_UWGCG = TRUE;
	    else is_UWGCG = FALSE;
	
	    if (( is_UWGCG)		/* decided */
		|| feof(bank->db_fp)	/* file finished */
		|| (i > 500))		/* ECOLAC UW/GenBank doc header is 300 lines */
		done = TRUE;
	    else done = FALSE;

	    if ((i < 5) && (bank->db_buffer[0] == ';')) {
		/* fix for ToIG of UWGCG, also NBRF/EMBL */
		is_UWGCG = FALSE;
		done = TRUE;
	    }
	    if ((i < 5) && (strindex(bank->db_buffer, ">>>>") == 0)) {
		/* test for GCG/Swissprot database:
		It begins with a line with several fields:
	NAME: .. LN: .. SN: .. SW REL: .. RELDATE: .. FORDATE: .. TYPE: ..
 		But I check here for the first data line. This lets the
		file ptr pointing to it, as expected. */
		bank->db_type = DB_GCGSW;
		is_UWGCG = TRUE;
		done = TRUE;
	    }
	}
	while (! done);
	if (is_UWGCG) {
	    /* and it is not the Sw protein database... */
	    if (bank->db_type != DB_GCGSW)
		bank->db_type = DB_UWGCG;
	}
	else {
	    rewind(bank->db_fp);
	    /* skip empty lines */
	    do {
		bank->db_pos = ftell(bank->db_fp);
		fgets(bank->db_buffer, MAXLINE, bank->db_fp);
		if (feof(bank->db_fp)) return NULL;
	    } while ( (bank->db_buffer[0] == '\n') ||
		      (bank->db_buffer[0] == '\0') );
	}
	
	if (! is_UWGCG)

	if (bank->db_buffer[0] == ';') {
	    if (strindex(bank->db_buffer, "Strider") > 0) {
		/* DNA Strider format */
		bank->db_type = DB_STRIDER;
	    }
	    else {
		/* Stanford/IG format */
		bank->db_type = DB_INTELLIGENETICS;
	    }
	}
	else if (strindex(bank->db_buffer, "LOCUS") == 0) {
	    /* GenBank format */
	    bank->db_type = DB_GENBANK;
	}
	else if ((strindex(bank->db_buffer, ">DL;") == 0)   /* NBRF */
	      || (strindex(bank->db_buffer, "> ") == 0) ) { /* PIR */
	    /* NBRF/PIR format */
	    bank->db_type = DB_PIR;
	}
	else if (bank->db_buffer[0] == '>') {
	    /* Pearson format */
	    bank->db_type = DB_PEARSON;
	}
	else if (strindex(bank->db_buffer, "ID   ") == 0) {
	    /* EMBL format */
	    bank->db_type = DB_EMBL;
	}
	else if (bank->db_buffer[0] == '(') {
	    /* Zuker format */
	    bank->db_type = DB_ZUKER;
	}
	else {
	    /* test for Fitch format */
	    int bufflen, is_Fitch;

	    fgets(bank->db_buffer, MAXLINE, bank->db_fp);
	    /* Do NOT update file pointer, to keep it pointing to the
	    beginning of the entry.
		This should be a seq line in 3 letter blocks separated
	    by spaces */
	    bufflen = strlen(bank->db_buffer);
	    i = 1;
	    do {
		if ( ((i % 4 == 0) && (bank->db_buffer[i] = ' '))
 		  || ((i % 4 != 0) && (bank->db_buffer[i] != ' ')) )
		    is_Fitch = TRUE;
		else
		    is_Fitch = FALSE;
		i++;
 	    }
	    while ( (i < bufflen) && (is_Fitch) );
	    if (is_Fitch) {
		/* Fitch format */
		bank->db_type = DB_FITCH;
		/* unget the test line read */
		fseek(bank->db_fp, bank->db_pos, 0);
		fgets(bank->db_buffer, MAXLINE, bank->db_fp);
	    }
	    else {
		/* Unknown format */
		bank->db_type = DB_UNKNOWN;
		rewind(bank->db_fp);
		bank->db_pos = ftell(bank->db_fp);
		fgets(bank->db_buffer, MAXLINE, bank->db_fp);
	    }
	}
    }
    switch (bank->db_type) {
    case    DB_UWGCG:	bank->db_next_protein = UWGCG_next_protein;
			break;
    case    DB_STRIDER:	bank->db_next_protein = Strider_next_protein;
			break;
    case    DB_INTELLIGENETICS: 
    case    DB_STANFORD:bank->db_next_protein = IG_next_protein;
			break;
    case    DB_GENBANK:	bank->db_next_protein = GenBank_next_protein;
			break;
    case    DB_NBRF:
    case    DB_PIR:	bank->db_next_protein = PIR_next_protein;
			break;
    case    DB_PEARSON:	bank->db_next_protein = Pearson_next_protein;
			break;
    case    DB_SWISSPROT:
    case    DB_EMBL:	bank->db_next_protein = EMBL_next_protein;
			break;
    case    DB_ZUKER:	bank->db_next_protein = Zuker_next_protein;
			break;
    case    DB_FITCH:	bank->db_next_protein = Fitch_next_protein;
			break;
    case    DB_PREER:	bank->db_next_protein = Preer_next_protein;
			break;
    case    DB_GCGSW:	bank->db_next_protein = GCGSW_next_protein;
			break;
    case    DB_UNKNOWN:
    default:		bank->db_next_protein = Unknown_next_protein;
			break;
    }

    return bank;
  }
                                         
/*-------------------------------------------------------------------
    GET PROTEIN
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

status db_next_protein(db, prot)
db_t db;
protein_t *prot;
  {
    return (*(db->db_next_protein))(db, prot);
  }

/*-------------------------------------------------------------------
    SEEK PROTEIN
	Gets the information about a protein from the current database.
    It access the file directly by using the supplied index.
    Obtains all the info about the protein and updates its structure
    appropriately

    INPUT:  P_INDEX -- Index of the sequence to retrieve
	    THE_PROT -- Protein structure to receive the sequence
	    DB -- Database to get the protein from
    OUTPUT: TRUE if success, FALSE if not.

    J. R. Valverde
 -------------------------------------------------------------------*/

status db_seek_protein(db, p_index, the_prot)
db_t db;
index_t *p_index;
protein_t *the_prot;
{
    fseek(db->db_fp, p_index->i_pos, SEEK_SET);
    db->db_pos = ftell(db->db_fp);
    if (fgets(db->db_buffer, MAXLINE - 1, db->db_fp) == NULL) return FALSE;      	    return db_next_protein(db, the_prot);
}

/*-------------------------------------------------------------------
    GET DBNAME
	Returns the current database file name by copying it into
    the supplied string.

    INPUT: DB_FILNAM -- string to which the name must be copied.
	    The routine ASSUMES that this string is long enough.

    J. R. Valverde
--------------------------------------------------------------------*/

void db_name(db, db_filnam)
db_t db;
char *db_filnam;
  {
    strcpy(db_filnam, db->db_nam);
  }

/*-------------------------------------------------------------------
    CLOSE DATABASE
	Closes specified database

    INPUT:  DB -- database to close

    J. R. Valverde
 -------------------------------------------------------------------*/

void db_close(db)
db_t db;
  {
    fclose(db->db_fp);
    free( db );
  }

/*--------------------------------------------------------------------
    DB_TYPE
	Return the format of the specified database

    INPUT: DB -- Database to check
    OUTPUT: the database format.

    J. R. VALVERDE
--------------------------------------------------------------------*/

int db_type(db)
db_t db;
{
    return (db->db_type);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		INDEX MANIPULATION ROUTINES

    But oh, beamish nephew, beware of the day,
    If your Snark be a Boojum! For then
    You will softly and suddenly vanish away,
    And never be met with again!

		Lewis Carroll. The hunting of the Snark.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*-------------------------------------------------------------------
    IDX_OPEN
	Open index file

    INPUT:  IDX_FILE -- Name of index file
	    MODE -- access mode
    OUTPUT: TRUE if success, FALSE if not.

    J. R. Valverde
--------------------------------------------------------------------*/

status idx_open(idx_file)
char *idx_file;
  {
    extern FILE *idx_fp;
    extern char idx_nam[MAXPATH + 1];

    if ((idx_fp = fopen(idx_file, "rw")) == NULL)
	return FALSE;

    rewind(idx_fp);

    strncpy(idx_nam, idx_file, MAXPATH);
                       
    return TRUE;
  }

/*--------------------------------------------------------------------
    IDX_CLOSE
	Close index file.

    J. R. Valverde
--------------------------------------------------------------------*/

status idx_close()
  {
     extern FILE *idx_fp;

     return (fclose(idx_fp));
  }

/*--------------------------------------------------------------------
    IDX_NAME
	Return name of last index file read

    INPUT:  NAME -- string to store the name in. It is assumed that
	it is long enough.

    J. R. Valverde
--------------------------------------------------------------------*/

status idx_name(name)
char *name;
  {
    extern char idx_nam[MAXPATH + 1];

    strcpy(name, idx_nam);
  }
                                
/*--------------------------------------------------------------------
    GET_INDEX
	Read index from index file

    INPUT:  INDEX   -- structure to hold the read index

    J. R. Valverde
--------------------------------------------------------------------*/

status get_index(index)
index_t *index;
  {
    static char idx_buff[MAXLINE];
    extern FILE *idx_fp;

    if (feof(idx_fp)) return FALSE;
    if (fgets(idx_buff, (MAXLINE - 1), idx_fp) == NULL)
	return FALSE;
    if ( sscanf(idx_buff, "%ld", &(index->i_pos)) != 0) {
	idx_buff[0] = '\0';
	return TRUE;
    }
    else
	return FALSE;
  }

/*--------------------------------------------------------------------
    SEARCH_INDEX
	Look for an index in the file given its associated
    keyword.
--------------------------------------------------------------------*/

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		    F O R    D E B U G G I N G

    ¡ Si, para muy poco tiempo !
    Más, como cada minuto
    Puede ser mi eternidad,
    ¡ Qué poco tiempo más único !

	    Juan Ramón Jiménez.	La Obra Desnuda.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#if DEBUG == 1
/*
horror(str1, str2)
char *str1, *str2;
{
  printf("ERROR. HORROR:\n%s: %s\n", str1, str2);
}
*/
protein_t   prot;

int histo[200] = { 0 };

main()
  {
    /* declared static to avoid register optimization on the VAX */
    static char dbname[128];
    static db_t db;
    static unsigned int i;
    static int result, maxlen;
    static index_t maxpidx;
    static time_t begtime, endtime;
    static long scnds;
    static FILE *fp;
#ifdef __TURBOC__
    int g_driver, g_mode;   /* for IBM-PC */
#endif

    i = 0;
    result = SUCCESS;
    maxlen = 0;
    
    printf("Database name: ");
    scanf("%s", dbname);

    time(&begtime);
    printf("BEGIN: %s\n", ctime(&begtime));

    db = db_open(dbname, DB_UNKNOWN);
    while ( ! feof(db->db_fp) )
      {
	if (! db_next_protein(db, &prot)) break;
	i++;
	if (maxlen < prot.p_seqlen) {
      	  maxlen = prot.p_seqlen;
	  maxpidx.i_pos = prot.p_index.i_pos;
	}
	histo[prot.p_seqlen / 40] += 1;
#ifdef __TURBOC__
	printf("%5d\t", i);
	if (kbhit()) {
	  getch();
	  break;
	}
#endif
      }                      

    time(&endtime);
    printf("\n\nEND: %s\n", ctime(&endtime));
    scnds = endtime - begtime;
    printf("\n%d seconds elapsed\n\n", scnds);
    getchar();

    printf("\n\n MAX LENGTH FOUND = %d at offset %ld\n",
	maxlen, maxpidx.i_pos);
    fseek(db->db_fp, maxpidx.i_pos, SEEK_SET);
    fgets(db->db_buffer, MAXLINE, db->db_fp);
    do {
	puts(db->db_buffer);
	fgets(db->db_buffer, MAXLINE, db->db_fp);
    } while (db->db_buffer[0] != '>');
    db->db_buffer[0] = '\0';
    fseek(db->db_fp, maxpidx.i_pos, SEEK_SET);
    fgets(db->db_buffer, MAXLINE, db->db_fp);
    db_next_protein(db, &prot);
    printf("\nNAME = %s\n", prot.p_name);
    printf("LENGTH = %d\n", prot.p_seqlen);
    puts("SEQUENCE: ");
    for (i = 0; i < prot.p_seqlen; i++)
      putchar(prot.p_seq[i]);
    db_seek_protein(db, &maxpidx, &prot);
    db_close(db);

    puts("\n");
    fp = fopen("DATA.OUT", "w+");
    for (i = 0; i < 100; i++) {
      printf("<%4d=\t%4d\t", (i+1) * 40, histo[i]);
      fprintf(fp, "%d\n", histo[i]);
    }
    fclose(fp);

#ifdef __TURBOC__
/*  -- for IBM-PC --	*/

    while (! kbhit()); getch();
    detectgraph(&g_driver, &g_mode);
    initgraph(&g_driver, &g_mode, "C:\\TURBO_C\\BGI");
    setcolor(1);
    rectangle(10, 10, 610, 210);
    outtextxy(10, 225, "0.0");
    settextjustify(RIGHT_TEXT, TOP_TEXT);
    outtextxy(605, 15, "800");
    outtextxy(605, 215, "4000 aa");
    moveto(10, 210);
    for (i = 0; i < 100; i++) {
      lineto( 10 + ((i + 1) * 6), 210 - (histo[i] / 4) );
    }
    while (! kbhit()); getch();
    closegraph();

#endif
}

#endif
