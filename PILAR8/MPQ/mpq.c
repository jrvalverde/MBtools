#include <stdio.h>

#define FALSE 0
#define TRUE 1

#define MAXLINE 	512
#define ASCII(x)	((x) & 0x7F)

int scores[ 20000 ];
char sequence[ 20001 ];
char name [MAXLINE + 1];

char buff[MAXLINE + 1];
long index, idx;

char dbnam[512];
char outnam[512];
FILE *db;
FILE *ifil;
FILE *out;

int aa_table[] = {
    -1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,-1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1, 1, 3, 2, 3, 4, 5, 6, 7, 8, 0, 9,10,11,12, 0,
    13,14,15,16,17, 0,18,19,20,21, 4, 0, 0, 0, 0, 0,
     0, 1, 3, 2, 3, 4, 5, 6, 7, 8, 0, 9,10,11,12, 0,
    13,14,15,16,17, 0,18,19,20,21, 4, 0, 0, 0, 0, 0,
};


main(argc, argv)
int argc;
char *argv[];
{
    int i, j;
    char c;

    for (i = 0; i < 20000; i++)
	scores[i] = 0;
    strncpy(dbnam, getenv("sw"), 511); 
    if ((db = fopen(dbnam, "r")) == NULL) {
	if (argc < 2) exit(4);
	if ((db = fopen(argv[1], "r")) == NULL) {
	    printf("No puedo abrir %s\n", argv[1]);
	    exit(4);
	}
    }
    if ((ifil = fopen("mpq.idx", "w+")) == NULL) {
	printf("No puedo abrir MPQ.IDX\n");
	exit(4);
    }
    if (fgets(buff, MAXLINE, db) == NULL) { /* eof */
	printf("No puedo leer la primera línea\n");
	exit(4);
    }
    i = j = 0;
    while (get_next_protein()) {

	j++;	/* count proteins found */

	/* locate first "PQ" */
	i = strindex(sequence, "PQ");

	/* Save index */
	if ((i <= 10) && (i > 0))
	    fprintf(ifil, "%ld\t%s\n", index, name);

	/* store statictics */
	scores[++i]++;
    }
    fclose(db);
    fclose(ifil);

    if ((out = fopen("mpq.out", "w")) == NULL) {
	printf("No puedo guardar MPQ.OUT\n");
	exit( 4 );
    }

    fprintf(out, "Número total de proteínas:\t%d\n\n", j);
    /* print statistics */
    if (scores[0] != 0)
	fprintf(out, "Proteinas sin el PQ:\t%d\tproteinas\n", scores[0]);
    for (i = 1; i < 20000; i++)
	if (scores[i] != 0)
	    fprintf(out, "Proteinas con PQ en posición\t%d:\t%d proteinas\n",
		    i, scores[i]);

    fclose(out);
    exit( 1 );
}

static int next_entry()
/*
 * Search for next entry into the database.
 *    It can be recognized because the line must begin with a '>>>'
 *   sequence.
 *    Also actualize file pointer as needed.
 */
  {
    while (strindex(buff, ">>>") != 0) {
      idx = ftell(db);
      if (fgets(buff, MAXLINE, db) == NULL)	/* eof */
	return FALSE;
    }
    return TRUE;
  }

static get_next_protein()
{
    int seqlen;
    char *aa;
    extern int next_entry();

    if (! next_entry())
	return FALSE;
    index = idx;

    /* buff contains now first line: ">>>> Access code" */
    strncpy(name, &buff[4], 10);
    /* get protein name */
    idx = ftell(db);
    if (fgets(buff, MAXLINE, db) == NULL)
	return FALSE;

    /* Get sequence */
    seqlen = 0;
    idx = ftell(db);
    if (fgets(buff, MAXLINE, db) == NULL)
	return FALSE;	    /* EOF */

    do {

      if (strindex(buff, ">>>") == 0)
	break;		    /* end of sequence */

      for (aa = buff; *aa != '\0'; aa++)
	if ((aa_table[ ASCII(*aa) ]) > 0) {
	  /* if it's an aminoacid add it to the sequence */
	  sequence[seqlen] = ASCII(*aa);
	  seqlen++;		     
	}
      idx = ftell(db);
    }		    /* while not eof */
    while ( fgets(buff, MAXLINE, db) != NULL );

    sequence[seqlen] = '\0';
    return TRUE;
}                                  

int strindex(str, sub)
char *str, *sub;
{
    int lensub, top, i;

    lensub = strlen(sub);
    top = strlen(str) - lensub;

    if (top < 0)
	return -1;		    /* substring longer than string */

    for (i = 0; i <= top; i++, str++)
	if (strncmp(str, sub, lensub) == 0)
	    return i;
    return -1;
}
