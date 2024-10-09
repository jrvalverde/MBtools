/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*																*
*							KERNEL . C							*
*																*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

#define DEBUG	1

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
		char *nam;			/* nombre del estándar					*/
		int num;			/* número de estándares					*/
		double *siz;		/* matríz con tamaños de los estándar	*/
		double *mov;		/* matríz con movilidad de los estándar	*/
		double *mcal;		/* matríz con movs. calculadas			*/
		double *scal;		/* matríz con tams. calculados			*/
		flag	modified;
			};

extern struct standard std;

struct problem {
		char *nam;			/* nombre del problema					*/
		int num;			/* número de bandas						*/
		double *siz;		/* matríz con los tamaños problema		*/
		double *mov;		/* matríz con las movilidades problema	*/
		struct problem *next, *prev;
		flag modified;
			};

extern struct problem *prbl;	/* También prb, pero no hace falta
								que éste módulo sepa que es *prbl	*/

struct session {
		char *nam;
		struct standard *std;
		struct problem	*prb;	/* matríz de problemas */
		flag sztyp;			/* tipo del tamaño: Kilobases, MegaDaltons.. */
		flag saved;
				};

extern struct session sess;

extern flag StdFlag;

extern flag SessFlag;

/**********************************************************************/
/*					BORRADORES PARA EL DESARROLLO					  */

#ifdef VALE   
   /* Imprimir la mov y tamaño calculados, calcular y
    sumar los restos y calcular e imprimir el error estándar
    del ajuste	*/
    
    /* fín lógico del programa */
    puts("\n");
    
    /* preguntar
    	¿Volver a usar el programa?
    	¿Los mismos estándares?
    	¿Corregir estándares?
    	¿Nuevos estándares?
    */
#endif

/*******************************************************************/

/* Pseudocódigo:
	0 ¿Leer Sesión?
		SESSFLAG = TRUE
	1 Leer estándar		(IF ! SESSFLAG)
		1.1 ¿Fichero?
			1.1.1 ¿movilidades?
				(si mov = CR ignorarla)
			1.1.2 ¿corregir?
			1.1.3 ¿guardar?
		1.2 ¿Nuevo?
			1.2.1 ¿nombre?
			1.2.2 ¿número de bandas?
			1.2.3 ¿peso, movilidad?
			1.2.4 ¿corregir?
			1.2.5 ¿guardar como fichero las movilidades?
		1.3 STDFLAG = TRUE
	2 Calcular coeficientes de regresión cuadrática (IF STDFLAG)
	3 Presentar el cálculo aplicado al estándar
		3.1 ¿Corregir estándar? ¿guardar? goto 2
	4 Pedir problema	(IF ! SESSFLAG) && (IF STDFLAG)
		4.1 ¿fichero?
		4.2 ¿nuevo?
			4.2.1 ¿nombre?
			4.2.2 ¿número de bandas?
			4.2.3 ¿peso?
	5 Calcular movils. problema	
	6 presentar movils. problema
	(IF SESSFLAG LOOP 5 UNTIL PROB.NEXT == NULL)
	7 ¿FIN?
		7.1 ¿nuevo problema?
		7.2 ¿nueva sesión (nuevos stds)?
		7.3 ¿salvar sesión?
			7.3.1 ¿Nombre?
			7.3.2 Std: mov, siz
			7.3.3 while sess.prb.next != NULL salvar prob.
		7.3 ¿fín?
*/

/*----------------- Leer estándares de un fichero ---------------*/
  
public int leeStdFile(size)
double *size;
  {
  /* Pide un fichero que contenga los valores estándar
  y los lee en un array intermedia hasta que se conozcan
  las movilidades */
    int result;
    str255 pfn;			/* Pascal file name */
    int volume;			/* Volumen en que está el fichero */
    FILE *pf;			/* el fichero */
    int i;
    
    /* Obtener nombre del fichero */
    strcpy(pfn, CtoPstr(""));
    if (OldFile(pfn, &volume) == 0)
      return FAIL;
    else
      {
        std.nam = PtoCstr(pfn);
        if (strlen(std.nam) == 0)
          return FAIL;
        /* abrir el fichero */
        pf = fopen(std.nam, "r");
        if (pf == NULL)
          return FAIL;
        fscanf(pf, "%d", &std.num);
        /* si no contiene nada */
        if (std.num == 0)
          {
            *std.nam = '\0';
            fclose(pf);
            return FAIL;
          }
        /* reservar espacio para los valores de tamaños */
        size = (double *) calloc(std.num, sizeof(double));
        /* leer los valores */
        for (i = 0; i < std.num; i++)
          fscanf(pf, "%f\t", &size[i]);
        /* dejarlo todo bien atado */
        fclose(pf);
       }
     return SUCCESS;
   }

/*---------------- Preguntar Movilidades Estándar ----------------*/

public int prgStdMov(size)
double *size;
/* 	Dados unos tamaños preexistentes del estándar, pregunta
  las movilidades correspondientes a cada fragmento.
    Si la movilidad introducida es 0 ignora el fragmento en
  cuestión, ajustando después std.num apropiadamente.		*/
  {
    int i, j, nulos;
    char movil[80];
    double *movs;
    
    cls();
    movs = (double *) calloc(std.num, sizeof(double));
    puts("Fragmento  Tamaño  Movilidad");
    puts("============================");
    puts("	  ?			?		0 = ignorar banda");
    for (i = nulos = 0; i < std.num; i++)
      {
        /* mostrar el peso estándar y preguntar su movilidad.
        Si mov = 0, ignorar banda */
        gotoxy(3, i + 5); printf("%d", i);
        gotoxy(13, i + 5); printf("%3.3f", size[i]);
        gotoxy(22, i + 5); gets(&movil);
        if (*movil == '\0')
          movs[i] = 0;
        else
          movs[i] = atoi(movil);
        if (movs[i] <= 0)
          nulos++;
      }
      
    /* Pasar los valores a std */
    if ((std.num - nulos) == 0) /* si no se desea ningún fragmento */
      {
        /* limpiar la memoria, innecesaria ya */
        std.num = 0;
        free(size);
        free(movs);
        return FAIL;
      }
    std.siz = (double *) calloc(std.num - nulos, sizeof(double));
    std.mov = (double *) calloc(std.num - nulos, sizeof(double));
    for (i = j = 0; i < std.num; i++)
      {
        if (movs[i] > 0)
          {
            /* si debemos considerar el fragmento, copiarlo */
            std.siz[j] = size[i];
            std.mov[j] = movs[i];
            j++;
          }
      }
    /* reajustar std.num */
    std.num -= nulos;
    /* liberar memoria innecesaria */
    free(size);
    free(movs);
    return SUCCESS;
  }

/*------------- Emplear fich para obtener estándares ------------*/

public int StdFile()
  {
    double *size;
    
    if (leeStdFile(size) == FAIL)
      return FAIL;
    if (prgStdMov(size) == FAIL)
      return FAIL;
    if (!StdFlag)
      StdFlag = TRUE;
    return SUCCESS;
  }
  
/*--------------------------------------------------------------------*/

public eraseStd(std)
struct standard *std;
  {
    if (std->nam != NULL)
      {
        free(std->nam);
        std->nam = NULL;
      }
    if (std->mov != NULL)
      {
        free(std->mov);
        std->mov = NULL;
      }
    if (std->siz != NULL)
      {
        free(std->siz);
        std->siz = NULL;
      }
    if (std->mcal != NULL)
      {
        free(std->mcal);
        std->mcal = NULL;
      }
    if (std->scal != NULL)
      {
        free(std->scal);
        std->scal = NULL;
      }
    std->num = 0;
  }

/*--------------------------------------------------------------------*/

freeStd(pstd)
struct standard *pstd;
  {
    if (pstd == NULL)
      return;
    eraseStd(pstd);
    free(pstd);
  }

/*--------------------------------------------------------------------*/

public erasePrbl(pprbl)
struct problem *pprbl;
  {
    if (pprbl->nam != NULL)
      {
        free(pprbl->nam);
        pprbl->nam = NULL;
      }
    if (pprbl->mov != NULL)
      {
        free(pprbl->mov);
        pprbl->mov = NULL;
      }
    if (pprbl->siz != NULL)
      {
        free(pprbl->siz);
        pprbl->siz = NULL;
      }
    pprbl->num = 0;
  }

/*---------------------------------------------------------------------*/

freePrbl(pprbl)
struct problem *pprbl;
/* Elimina un problema de la lista de problemas */
  {
    struct problem *prbABorrar;
    
    prbABorrar = pprbl;
    /* si no hay nada que borrar */
    if (pprbl == NULL)
      return;
    /* borramos el contenido */
    erasePrbl(pprbl);
    /* si es el primero de la lista */
    if (pprbl->prev == NULL)
      {
        /* y no es el único */
        if (pprbl->next != NULL)
          {
            pprbl->next->prev = NULL;
            /* readjudicamos el puntero para evitar inconsistencias */
            pprbl = pprbl->next;
            free(prbABorrar);
            return;
          }
        else
          {
            /* si es el único */
            pprbl = NULL;
            free(prbABorrar);
            return;
          }
      }
    else
    /* no es el primero de la lista */
      {
        /* si estamos al final de la lista */
        if (pprbl->next == NULL)
          {
            /* actualizamos el final de la cola */
            pprbl->prev->next == NULL;
            pprbl = pprbl->prev;
            /* liberamos la estructura */
            free(prbABorrar);
            return;
          }
        else
          /* no está al principio ni al final */
          {
            /* actualizamos el orden de la lista y
               dejamos de lado a pprbl				*/
            pprbl->prev->next = pprbl->next;
            pprbl = pprbl->next;
            /* liberamos el espacio ocupado por pprbl */
            free(prbABorrar);
            return;
          }
      }
  }

/*---------------------------------------------------------------------*/

public eraseSess(psess)
/*
		¡¡¡		MUCHO OJO CON ESTA FUNCION		!!!
*/
struct session *psess;
  {
    struct problem *pprb1, *pprb2;
    
    free(psess->nam);
    eraseStd(psess->std);	/* std en session es un puntero a standard */
    free(psess->std);
    psess->std = NULL;
    /* Borramos toda la lista de problemas */
    pprb1 = psess->prb;
    do
      {
        pprb2 = pprb1;
        pprb1 = pprb1->next;
        freePrbl(pprb2);
      }
    while (pprb1->next != NULL);
    psess->saved = FALSE;
    /* Ahora la sesión está vacía */
  }

/*---------------------------------------------------------------------*/

freeSess(psess)
struct session *psess;
  {
    eraseSess(psess);
    free(psess);
  }

/*---------------------------------------------------------------------*/

public int leeSession(sess)
struct session sess;
  {
    int i;
    FILE *pf;
    char *fichnam;
    int vRef;
    flag done;
    struct problem *pprob;
    
    /* Pedimos el nombre del fichero */
    fichnam = malloc(sizeof(str255));
    NewFile(fichnam, vRef);
    if (sess.nam != NULL)
      free(sess.nam);
    sess.nam = PtoCstr(fichnam);
    
    /* Intentamos abrirlo */
    pf = fopen(sess.nam, "r");
    if (pf == NULL)
      return FAIL;
      
    /* Intentamos leer el contenido:
    
    	1. nombre de la sesión:
    	si fallamos es que no hay nada de interés */
    fscanf(pf, "%s\n", &sess.nam);
    if (feof(pf))
      {
        sess.nam = NULL;
        fclose(pf);
        return FAIL;
      }
      
    /*	2. valores estándar:
    	primero el nombre y si lo logramos
    	después todo lo demás.					*/
    sess.std = (struct standard *) malloc(sizeof(struct standard));
    sess.std->nam = malloc(sizeof(str255));
    fscanf(pf, "%s\n", &sess.std->nam);
    if (feof(pf))
      {
        free(sess.std->nam);
        sess.std->nam = NULL;
        free(sess.std);
        sess.std = NULL;
        fclose(pf);
        return FAIL;
      }
    fscanf(pf, "%d\n", &sess.std->num);
    
    /*	para los valores hemos de reservar espacio de memoria */
    sess.std->mov = (double *) calloc(sess.std->num, sizeof(double));
    sess.std->siz = (double *) calloc(sess.std->num, sizeof(double));
    for (i = 0; i < sess.std->num; i++)
      {
        fscanf(pf, "%g\t%g\n", &sess.std->mov[i], &sess.std->siz[i]);
      }

    /* intentamos leer primer problema */
    sess.prb = (struct problem *) malloc(sizeof(struct problem));
    pprob = sess.prb;
    
    /* si lo hay, lo leeremos con los siguientes */
    pprob->nam = malloc(sizeof(str255));
    fscanf(pf, "%s\n", &pprob->nam);
   
    /* si no lo hay hemos terminado */
    while (! feof(pf))
      {
        /* leemos los valores */
        fscanf(pf, "%d\n", &pprob->num);
        pprob->mov = (double *) calloc(pprob->num, sizeof(double));
        for (i = 0; i < pprob->num; i++)
          {
            fscanf(pf, "%g\t%g\n", &pprob->mov[i], &pprob->siz[i]);
          }
        
        /* buscamos el siguiente */
        pprob->next = (struct problem *) malloc(sizeof(struct problem));
        pprob = pprob->next;
        pprob->nam = malloc(sizeof(str255));
        fscanf(pf, "%s\n", &pprob->nam);
      }
    
    /* al salir del bucle tenemos reservado un pprob y
       un pprob->nam que no existen: los eliminamos		*/
    free(pprob->nam);
    free(pprob);
    pprob = NULL;
    /* esto es porque así sess.prb = NULL si ya el primero
       no estaba, y pprob->next = NULL si es que no había
       más datos detrás de éstos.						*/
    
    /* y adios */
    fclose(pf);
    sess.saved = TRUE;
    return SUCCESS;
  }

/*--------------------------------------------------------------------*/

