
#include <QuickDraw.h>
#include <MacTypes.h>
#include <FontMgr.h>
#include <WindowMgr.h>
#include <MenuMgr.h>
#include <TextEdit.h>
#include <DialogMgr.h>
#include <EventMgr.h>
#include <DeskMgr.h>
#include <FileMgr.h>
#include <ToolboxUtil.h>
#include <ControlMgr.h>
#include <StdFilePkg.h>

#include "stdio.h"
#include "portable.h"
#include "math.h"
#include "pascal.h"
#include <storage.h>
#include <unix.h>

/*------------------------ DEFINICIONES PROPIAS ------------------------*/

#define KB		0
#define MD		1
#define OTHERS	2

/*------------------------ TIPOS Y VARIABLES --------------------------*/

typedef char str255[256];

struct standard {
		char *nam;			/* nombre del est‡ndar					*/
		int num;			/* nœmero de est‡ndares					*/
		double *siz;		/* matr’z con tama–os de los est‡ndar	*/
		double *mov;		/* matr’z con movilidad de los est‡ndar	*/
		double *mcal;		/* matr’z con movs. calculadas			*/
		double *scal;		/* matr’z con tams. calculados			*/
		flag	modified;
			}
		std = {NULL, 0, NULL, NULL, NULL, NULL, FALSE};

struct standard *pstd = &std;

/* en el programa definitivo deber‡ inicializarse
   al principio										*/

struct problem {
		char *nam;			/* nombre del problema					*/
		int num;			/* nœmero de bandas						*/
		double *siz;		/* matr’z con los tama–os problema		*/
		double *mov;		/* matr’z con las movilidades problema	*/
		struct problem *next, *prev;
		flag modified;
			}
		prob = {NULL, 0, NULL, NULL, NULL, NULL, FALSE};

struct problem *pprob = &prob;
/* Lo mismo */

struct session {
		char *nam;
		struct standard *std;
		struct problem	*prb;	/* matr’z de problemas */
		flag sztyp;			/* tipo del tama–o: Kilobases, MegaDaltons.. */
		flag saved;
				}
		sess = {NULL, NULL, NULL, KB, FALSE};

struct session *psess = &sess;

flag StdFlag = FALSE;

flag SessFlag = FALSE;


main()
  {
    Point UpLtCorner;
    str255 prompt, origName;
    char *cfn;
    char ch;
    int volume, i, j;
    FILE *pf;
    struct problem *prblist;
    double c1, c2, c3, sterror;
    extern struct standard *pstd;
    extern struct problem *pprob;
    extern struct session *psess;
    
    InitGraf(&thePort);
	InitFonts();
	InitWindows();
	InitMenus();
	TEInit();
	InitDialogs(0L);
	InitCursor();
	FlushEvents( everyEvent, 0 );

	/* Apertura inocente de un fichero
	strcpy(prompt, "kk-kk-kk-kk");
	cfn =  CtoPstr(prompt);
	printf("P %s\n%s\n", prompt, cfn);
	while (! kbhit()); getch();
    NewFile(&prompt, &volume);
    cfn = PtoCstr(prompt);
    printf("C %s\n%s\n", prompt, cfn);
	while (! kbhit()); getch();
    if ((pf = fopen(cfn, "a")) == NULL)
      puts("ERROR");
    fclose(pf);
    while (! kbhit());					*/
    
    i = Advise("\pPESOS MOLECULARES");
    if (i == 2)
      exit();
    
    /* Chequeos varios */
    cls();
    do
      {
        i = newStd(pstd);
      }
    while (i != SUCCESS);
    calcCoefs(&c1, &c2, &c3, pstd);
    calcStd(c1, c2, c3, pstd, &sterror);
    showStd(pstd, &sterror);
    
    do
      {
        i = prgPrbl(pprob);
      }
    while (i != SUCCESS);
    calcPrbSiz(c1, c2, c3, pprob);
    showPrbl(pprob, pstd);
    cls();
    printf("ÀOtro problema? (S / N) ");
    ch = getchar();
   	cls();
   	prblist = pprob;
   	/* pprob apunta a la base de la lista de problemas */
   	while ((ch == 's') || (ch == 'S'))
   	  {
   	    /* introducimos un nuevo problema al final de la lista */
   	    addPrbl(prblist);
   	    /* y actualizamos el puntero de la lista al final de
   	    la misma.											*/
   	    prblist = prblist->next;
   	    i = prgPrbl(prblist);
    	calcPrbSiz(c1, c2, c3, prblist);
        showPrbl(prblist, pstd);
        cls();
        printf("ÀOtro problema? (S / N) ");
        ch = getchar();
      }
    /* Preparaci—n de una sesi—n de trabajo para guardarla en fichero */
    psess->std = pstd;
    psess->prb = pprob;		/* apunta a la base de la lista de problemas */
    saveSession(psess);
    SayOK("\p       That's all folks !");
    exit();
   }
    
