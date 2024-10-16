/*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*																*
*							KERNEL . C							*
*																*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*/

/* A�adir ordenaci�n de pesos y tama�os */

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

struct standard
{
  char *nam;			/* nombre del est�ndar                                  */
  int num;			/* n�mero de est�ndares                                 */
  double *siz;			/* matr�z con tama�os de los est�ndar   */
  double *mov;			/* matr�z con movilidad de los est�ndar */
  double *mcal;			/* matr�z con movs. calculadas                  */
  double *scal;			/* matr�z con tams. calculados                  */
  flag modified;
};

extern struct standard *pstd;

struct problem
{
  char *nam;			/* nombre del problema                                  */
  int num;			/* n�mero de bandas                                             */
  double *siz;			/* matr�z con los tama�os problema              */
  double *mov;			/* matr�z con las movilidades problema  */
  struct problem *next, *prev;
  flag modified;
};

extern struct problem *prbl;	/* Tambi�n prb, pero no hace falta
				   que �ste m�dulo sepa que es *prbl    */

struct session
{
  char *nam;
  struct standard *std;
  struct problem *prb;		/* matr�z de problemas */
  flag sztyp;			/* tipo del tama�o: Kilobases, MegaDaltons.. */
  flag saved;
};

extern struct session sess;

extern flag StdFlag;

extern flag SessFlag;


#define SCR_WIDTH	80
#define SCR_HEIGHT	25

/*------------------------------ cls -----------------------------*/

public
cls ()
/* Borra la pantalla  */
/* OK. 1 - nov - 1988 */
{
  printf ("\f");
}

/*-------------------------- Men� Principal ----------------------*/

public int
mainMenu ()
/*		OK.		1 - nov - 1988		*/
{
  char opcion[80];
  int result;

  cls ();
  fflush (stdin);
  puts ("\t\t\tBIENVENIDO");
  puts ("\t\t\t==========");
  gotoxy (10, 5);
  puts ("1. Leer Sesi�n anterior");
  gotoxy (10, 7);
  puts ("2. Leer Est�ndar");
  gotoxy (10, 9);
  puts ("3. Introducir Problema");
  gotoxy (10, 11);
  puts ("4. Salvar sesi�n");
  gotoxy (10, 13);
  puts ("5. Salir");
  do
    {
      gotoxy (10, 20);
      puts ("                       ");
      gotoxy (10, 20);
      gets (opcion);
      if (*opcion == '\0')
	result = 0;
      else
	result = atoi (opcion);
    }
  while ((result < 1) || (result > 5));
  cls ();
  fflush (stdin);
  return result;
}

/*------------------------- Men� de est�ndares ----------------------*/

public int
stdMenu ()
/*		OK.		1 - nov - 1988		*/
{
  char opcion[80];
  int result;

  cls ();
  fflush (stdin);
  puts ("\t\t\tLEER ESTANDAR");
  puts ("\t\t\t=============");
  gotoxy (10, 10);
  puts ("1. Leer fichero existente");
  gotoxy (10, 12);
  puts ("2. Introducir nuevos valores");
  gotoxy (10, 14);
  puts ("3. Corregir valores");
  gotoxy (10, 16);
  puts ("4. Salvar valores");
  gotoxy (10, 18);
  puts ("5. Men� Principal");
  do
    {
      gotoxy (10, 20);
      puts ("                        ");
      gotoxy (10, 20);
      gets (opcion);
      if (*opcion == '\0')
	result = 0;
      else
	result = atoi (opcion);
    }
  while ((result < 1) || (result > 5));
  cls ();
  fflush (stdin);
  return result;
}

/*------------------------ Est�ndar es nuevo ----------------------*/

public int
newStd (pstd)
     struct standard *pstd;
/*		OK.	1 - nov - 1988		*/
{
  if (prgStdNam (pstd) == FAIL)
    return FAIL;
  if (prgNumStd (pstd) == FAIL)
    return FAIL;
  if (prgStd (pstd) == FAIL)
    return FAIL;
  return SUCCESS;
}

/*--------------- Preguntar el nombre del est�ndar ----------------*/

public int
prgStdNam (std)
     struct standard *std;
/*		OK.		1 - nov - 1988		*/
{
  cls ();
  fflush (stdin);
  puts ("\t\t\tINTRODUCIR ESTANDAR");
  puts ("\t\t\t===================");
  printf ("\n\tNombre del est�ndar: ");
  std->nam = malloc (sizeof (str255));
  while (!kbhit ());
  gets (std->nam);
  if (strlen (std->nam) == 0)
    {
      free (std->nam);
      return FAIL;
    }
  else
    return SUCCESS;
}

/*---------------- Preguntar n�mero de fragmentos ----------------*/

public int
prgNumStd (std)
     struct standard *std;
/*		OK.		1 - nov - 1988		*/
{
  char num[80];

  cls ();
  fflush (stdin);
  puts ("\t\t\tINTRODUCIR ESTANDAR");
  puts ("\t\t\t===================");
  printf ("\n\tN�mero de fragmentos = ");
  gets (num);
  if (*num == '\0')
    std->num = 0;
  else
    std->num = atoi (num);
  if (std->num <= 3)
    {
      int item;

      item = Advise ("\pEl n�mero de est�ndares es insuficiente.");
      if (item == 2)		/* Cancel */
	return FAIL;
      /* si no, es que ha pulsado OK y desea seguir */
    }
  if (std->num <= 0)
    {
      std->num = 0;
      return FAIL;
    }
  return SUCCESS;
}

/*----------------- Preguntar valores del est�ndar ------------------*/

public int
prgStd (std)
     struct standard *std;
/*
	Pregunta los valores de peso y movilidad del est�ndar.
	Ignora los valores introducidos como cero.
	Ajusta Std.Num apropiadamente.
	
	OK.		1 - nov - 1988
*/
{
  int i, j, nulos, linea;
  double *size, *movs;
  char ssiz[80], smov[80];

  size = (double *) calloc (std->num, sizeof (double));
  movs = (double *) calloc (std->num, sizeof (double));
  cls ();
  fflush (stdin);
  puts ("\t\t\tINTRODUCIR ESTANDAR");
  puts ("\t\t\t===================");
  puts ("\nFragmento	Movilidad	Peso");
  puts ("============================");
  linea = 4;
  for (i = nulos = 0; i < std->num; i++)
    {
      linea++;
      if (linea == 23)
	{
	  for (j = 5; j <= 24; j++)
	    {
	      gotoxy (0, j);
	      printf ("                                          ");
	    }
	  linea = 5;
	}
      gotoxy (0, linea);
      printf ("%3d\t\t\t", i + 1);
      gets (smov);
      gotoxy (24, linea);
      gets (ssiz);
      if (ssiz[0] == '\0')
	size[i] = 0;
      else
	size[i] = atof (ssiz);
      if ((smov[0] == '\0'))
	movs[i] = 0;
      else
	movs[i] = atof (smov);
      if ((size[i] <= 0) || (movs[i] <= 0))
	nulos++;
    }
  if ((std->num - nulos) <= 3)
    {
      int item;

      item = Advise ("\p�El n�mero de est�ndares es insuficiente!");
      if (item == 2)
	nulos = std->num;
      /* para que a continuaci�n termine en el siguiente if */
      /* si no, es que desea continuar */
    }
  if ((std->num - nulos) == 0)
    {
      std->num = 0;
      free (size);
      free (movs);
      return FAIL;
    }

  /* Copiar los valores a std */

  /* 1. Reservar espacio de memoria */
  std->siz = (double *) calloc (std->num - nulos, sizeof (double));
  std->mov = (double *) calloc (std->num - nulos, sizeof (double));

  /* 2. Copiar solamente los no nulos */
  for (i = j = 0; i < std->num; i++)
    {
      if ((movs[i] > 0) && (size[i] > 0))
	{
	  /* hacerlo mediante inserci�n en lista ordenada */
	  std->siz[j] = size[i];
	  std->mov[j] = movs[i];
	  j++;
	}
    }

  /* Reajustar std.num */
  std->num -= nulos;

  /* liberar memoria innecesaria */
  free (size);
  free (movs);
  return SUCCESS;
}

/*------------------------- Preguntar problemas -----------------------*/

public int
prgPrbl (prbl)
     struct problem *prbl;
/*		OK. 1 - nov - 1988		*/
{
  int result;

  result = prgPrNam (prbl);
  if (result == FAIL)
    return FAIL;
  result = prgPrNum (prbl);
  if (result == FAIL)
    return FAIL;
  result = prgPrMov (prbl);
  if (result == FAIL)
    return FAIL;
  return SUCCESS;
}

/*------------------------ Hallar nombre de Prbl ----------------------*/

public int
prgPrNam (prbl)
     struct problem *prbl;
/* OK.	1 - nov - 1988		*/
{
  cls ();
  fflush (stdin);
  prbl->nam = malloc (sizeof (str255));
  puts ("\t\t\tINTRODUCIR PROBLEMA");
  puts ("\t\t\t===================\n");
  puts ("\t\tIntroduzca el nombre de los valores problema\n");
  gets (prbl->nam);
  if (strlen (prbl->nam) == 0)
    {
      free (prbl->nam);
      return FAIL;
    }
  else
    return SUCCESS;
}

/*---------------------- Hallar n�mero de bandas ----------------------*/

public int
prgPrNum (prbl)
     struct problem *prbl;
/*		OK.		1 - nov - 1988		*/
{
  char num[80];

  cls ();
  fflush (stdin);
  puts ("\t\t\tINTRODUCIR PROBLEMA");
  puts ("\t\t\t===================\n");
  puts ("\t\tIndique el n�mero de bandas problema\n");
  gets (num);
  if (*num == '\0')
    prbl->num = 0;
  else
    prbl->num = atoi (num);
  if (prbl->num <= 0)
    {
      prbl->num = 0;
      return FAIL;
    }
  return SUCCESS;
}

/*----------------------- Pedir movilidades prbl ----------------------*/

public int
prgPrMov (prbl)
     struct problem *prbl;
/*		OK.		1 - nov - 1988		*/
{
  int i, j;
  int nulos, linea;
  char stmov[80];
  double *mov;

  fflush (stdin);
  if (prbl->num <= 0)
    {
      SayError ("\pError: �Horror!", " N�mero incorrecto de bandas problema");
      return FAIL;
    }

  mov = (double *) calloc (prbl->num, sizeof (double));
  cls ();
  fflush (stdin);
  puts ("\t\t\tINTRODUCIR PROBLEMA");
  puts ("\t\t\t===================\n");
  puts ("Banda		Movilidad");
  puts ("===================");
  linea = 4;
  for (i = nulos = 0; i < prbl->num; i++)
    {
      linea++;
      if (linea == 23)
	{
	  for (j = 5; j < 24; j++)
	    {
	      gotoxy (0, j);
	      printf ("                               ");
	    }
	  linea = 5;
	}
      gotoxy (0, linea);
      printf ("%d\t\t\t", i + 1);
      gets (stmov);
      if (strlen (stmov) == 0)
	{
	  mov[i] = 0;
	  nulos++;
	}
      else
	{
	  mov[i] = atof (stmov);
	  if (mov[i] <= 0)
	    {
	      mov[i] = 0;
	      nulos++;
	    }
	}
    }
  if ((prbl->num - nulos) == 0)
    {
      prbl->num = 0;
      return FAIL;
    }
  else
    {
      prbl->mov = (double *) calloc (prbl->num - nulos, sizeof (double));
      for (i = j = 0; i < prbl->num; i++)
	{
	  if (mov[i] > 0)
	    {
	      prbl->mov[j] = mov[i];
	      j++;
	    }
	}
      prbl->num -= nulos;
    }
  free (mov);
  return SUCCESS;
}

/*------------------------ Calcular coeficientes ----------------------*/

public int
calcCoefs (coef1, coef2, coef3, pstd)
     double *coef1, *coef2, *coef3;
     struct standard *pstd;
/* Calcular los coeficientes de regresi�n cuadr�tica
	tales que
		 2			1			0
		x  Coef3 + x  Coef2 + x   Coef1				*/
/*		OK.		10 - nov - 1988		*/
{
  double tmp = 0.0,
    tmp_2 = 0.0,
    mvtmp = 0.0,
    coeff1, coeff2, coeff3, arrtmp, a2, a3, a4, p1, p2, r1, r2, d;
  int i;

  coeff1 = *coef1;
  coeff2 = *coef2;
  coeff3 = *coef3;
  if ((pstd->num == 0) || (pstd->mov == NULL) || (pstd->siz == NULL))
    {
      SayError ("Error: �Horror!", " No hay est�ndar disponible.");
      return FAIL;
    }
  tmp = tmp_2 = mvtmp = 0.0;
  for (i = 0; i < pstd->num; i++)	/* n1 = n� de �standares */
    {
      arrtmp = log (pstd->siz[i]);
      tmp += arrtmp;
      tmp_2 += arrtmp * arrtmp;
      mvtmp += pstd->mov[i];
    }
  tmp /= pstd->num;
  tmp_2 /= pstd->num;
  mvtmp /= pstd->num;
  a2 = a3 = a4 = p1 = p2 = 0.0;
  for (i = 0; i < pstd->num; i++)
    {
      arrtmp = log (pstd->siz[i]);
      r1 = arrtmp - tmp;
      r2 = (arrtmp * arrtmp) - tmp_2;
      a2 += r1 * r1;
      a3 += r1 * r2;
      a4 += r2 * r2;
      p1 += r1 * pstd->mov[i];
      p2 += r2 * pstd->mov[i];
    }
  d = ((a2 * a4) - a3 * a3);
  a2 /= d;
  a3 = -a3 / d;
  a4 /= d;
  coeff3 = a3 * p1 + a2 * p2;
  coeff2 = a4 * p1 + a3 * p2;
  coeff1 = mvtmp - coeff2 * tmp - coeff3 * tmp_2;
  *coef1 = coeff1;
  *coef2 = coeff2, *coef3 = coeff3;
  return SUCCESS;
}

/*-----------------------------------------------------------------*/

public int
calcStd (coef1, coef2, coef3, pstd, perrStd)
     double coef1, coef2, coef3;
     struct standard *pstd;
     double *perrStd;
/*		OK.		10 - nov - 1988		*/
{
  int i;
  double z, f, z1, d;

  if ((pstd->num == 0) || (pstd->mov == NULL) || (pstd->siz == NULL))
    {
      SayError ("Error: �Horror!", " No hay est�ndar disponible.");
      return FAIL;
    }
  /* Reservar espacio de memoria para std->mcal y std->scal */
  /* primero borrarlos si existen */
  if (pstd->mcal != NULL)
    free (pstd->mcal);
  if (pstd->scal != NULL)
    free (pstd->scal);
  pstd->mcal = (double *) calloc (pstd->num, sizeof (double));
  pstd->scal = (double *) calloc (pstd->num, sizeof (double));
  z1 = 0;
  for (i = 0; i < pstd->num; i++)
    {
      z = log (pstd->siz[i]);
      f = coef1 + (coef2 * z) + (coef3 * z * z);	/* mov calculada */
      z1 += ((f - pstd->mov[i]) * (f - pstd->mov[i]));
      z = pstd->mov[i];
      d = (coef2 * coef2) - (4.0 * coef3 * (coef1 - z));
      z = -(sqrt (d) + coef2) / (coef3 + coef3);
      pstd->mcal[i] = f;
      pstd->scal[i] = exp (z);
    }
  *perrStd = sqrt (z1 / (pstd->num - 3.0));
  return SUCCESS;
}

/*---------------------------------------------------------------------*/

public int
showStd (pstd, perrStd)
     struct standard *pstd;
     double *perrStd;
/*		OK.		10 - nov - 1988		*/
{
  int i, j, linea;

  if ((pstd->num == 0) || (pstd->mov == NULL) || (pstd->siz == NULL))
    {
      SayError ("\pError: �Horror!", " No hay est�ndar disponible.");
      return FAIL;
    }

  /* si los valores a�n no han sido calculados, calcularlos */
  if ((pstd->mcal == NULL) || (pstd->scal == NULL))
    {
      double coef1, coef2, coef3, errStd;

      calcCoefs (&coef1, &coef2, &coef3, pstd);
      calcStd (coef1, coef2, coef3, pstd, &errStd);
      *perrStd = errStd;
    }

  /* Escribir cabeceras */
  cls ();

  /* Imprimir la mov y tama�o calculados, calcular y
     sumar los restos y calcular e imprimir el error est�ndar
     del ajuste       */
  gotoxy (0, 0);
  printf ("\t			Ajuste del estandard %s\n", pstd->nam);
  puts
    ("\t		Movilidad						Tama�o");
  puts
    ("\tReal		Calculada			Real		Calculado");
  puts ("\t=====================================================\n");
  linea = 4;
  for (i = 0; i < pstd->num; i++)
    {
      linea++;
      if (linea == 23)
	{
	  while (!kbhit ());
	  getch ();
	  for (j = 0; j < 24; j++)
	    {
	      gotoxy (0, j);
	      printf
		("                                                                      ");
	    }
	  linea = 5;
	}
      gotoxy (0, linea);
      printf ("%d\t%3.3f\t\t%3.3f\t\t\t\t%3.3f\t\t%3.3f\n", i + 1,
	      pstd->mov[i], pstd->mcal[i], pstd->siz[i], pstd->scal[i]);
    }
  puts ("\n  ERROR ESTANDAR ESTIMADO");
  puts ("===========================\n");
  printf ("%g\n", *perrStd);

  /* esperar */
  while (!kbhit ());
  getch ();
  cls ();
}

/*------------------------------------------------------------------*/

public int
calcPrbSiz (coef1, coef2, coef3, pprbl)
     double coef1, coef2, coef3;
     struct problem *pprbl;
/*		OK.		10 - nov - 1988		*/
{
  int i;
  double d, z;

  if (pprbl == NULL)
    return FAIL;
  if ((pprbl->num == 0) || (pprbl->mov == NULL))
    {
      SayError ("Error: �Horror!", " No dispongo de valores problema.");
      return FAIL;
    }
  if (pprbl->siz != NULL)
    free (pprbl->siz);
  pprbl->siz = (double *) calloc (pprbl->num, sizeof (double));
  for (i = 0; i < pprbl->num; i++)
    {
      d = (coef2 * coef2) - (4.0 * coef3 * (coef1 - pprbl->mov[i]));
      z = -(sqrt (d) + coef2) / (coef3 + coef3);
      pprbl->siz[i] = exp (z);
    }
  return SUCCESS;
}

/*---------------------------------------------------------------------*/

public int
showPrbl (pprbl, pstd)
     struct problem *pprbl, *pstd;
/*		OK.		10 - nov - 1988		*/
{
  int i, j, linea;

  if (pprbl == NULL)
    return FAIL;
  if (pprbl->mov == NULL)
    {
      SayError ("Error: �Horror!", " No hay valores problema.");
      return FAIL;
    }
  else if (pprbl->siz == NULL)
    {
      double coef1, coef2, coef3, errStd;

      if (calcCoefs (&coef1, &coef2, &coef3, pstd, &errStd) == FAIL)
	return FAIL;
      calcPrbSiz (coef1, coef2, coef3, pprbl);
      /* no hace falta chequear su salida: s�lo falla si
         pprbl->mov == NULL y ya est� chequeado                             */
    }
  cls ();
  gotoxy (0, 0);
  printf ("\t	Ajuste de los pesos de %s\n", pprbl->nam);
  puts ("\tMovilidad						Tama�o");
  puts ("\t======================================\n");
  linea = 3;
  for (i = 0; i < pprbl->num; i++)
    {
      linea++;
      if (linea == 23)
	{
	  while (!kbhit ());
	  getch ();
	  for (j = 5; j < 24; j++)
	    {
	      gotoxy (0, j);
	      printf
		("                                                                 ");
	    }
	  linea = 5;
	}
      gotoxy (0, linea);
      printf ("%d\t%3.3f\t\t\t\t\t\t%3.3f\n", i + 1, pprbl->mov[i],
	      pprbl->siz[i]);
    }

  while (!kbhit ());
  getch ();
  cls ();
}

/*---------------------------------------------------------------------*/

public
addPrbl (prblist)
     struct problem *prblist;
/* a�ade un elemento a la lista solicitada */
/*	OK.		10 - nov - 1988		*/
{
  struct problem *pNuevo;

  pNuevo = (struct problem *) malloc (sizeof (struct problem));

  /* inicializamos el nuevo elemento */
  pNuevo->nam = NULL;
  pNuevo->num = 0;
  pNuevo->siz = NULL;
  pNuevo->mov = NULL;
  pNuevo->next = NULL;
  pNuevo->modified = FALSE;
  pNuevo->prev = prblist;

  /* y lo a�adimos */
  if (prblist == NULL)
    {
      prblist = pNuevo;
    }
  else if (prblist->next == NULL)
    {
      prblist->next = pNuevo;
    }
  else
    {
      /* si no es el �ltimo de la lista
         lo insertamos a continuaci�n de prblist */
      prblist->next->prev = pNuevo;
      pNuevo->next = prblist->next;
      prblist->next = pNuevo;
    }
}

/*---------------------------------------------------------------------*/

public int
saveSession (psess)
     struct session *psess;
/*	OK.		21 - nov - 1988		*/
{
  int i;
  FILE *pf;
  struct problem *pprob;

  if (psess == NULL)		/* No hay sesi�n */
    return SUCCESS;
  if (psess->saved == TRUE)
    return SUCCESS;		/* ya est� salvada */
  if (psess->std == NULL)
    return SUCCESS;		/* est� vac�a   */
  /* si NO hay valores problema no se guardar�n, s�lo los est�ndar. */

  /* Pedimos el nombre del fichero */
  if (psess->nam == NULL)
    {
      int result, vRef;
      str255 filename;

      strcpy (filename, "\pMySession");
      result = NewFile (&filename, &vRef);
      if (result == 0)
	/* si ha elegido cancelar es que ha cambiado
	   de opini�n y no quiere salvar                              */
	return FAIL;
      PtoCstr (filename);
      /* reservamos para (char *nam)
         el tama�o m�ximo que puede ocupar */
      psess->nam = calloc (1, sizeof (str255));
      strcpy (psess->nam, filename);
    }

  /* intentamos abrirlo */
  pf = fopen (psess->nam, "w");
  if (pf == NULL)
    return FAIL;

  /* Empezamos guardando su nombre y el est�ndar */
  fprintf (pf, "%s\n", psess->nam);
  fprintf (pf, "%s\n", psess->std->nam);
  fprintf (pf, "%d\n", psess->std->num);
  for (i = 0; i < psess->std->num; i++)
    {
      fprintf (pf, "%g\t", psess->std->mov[i]);
      fprintf (pf, "%g\n", psess->std->siz[i]);
    }

  /* Si hay datos problema los salvamos tambi�n */
  /* si hay movilidades problema pero no se han calculado los tama�os
     se calculan pero no se ense�an al usuario */
  pprob = psess->prb;
  while (pprob != NULL)
    {
      fprintf (pf, "%s\n", pprob->nam);
      fprintf (pf, "%d\n", pprob->num);
      if (pprob->siz == NULL)
	{
	  double coef1, coef2, coef3, errStd;

	  /* calcular los tama�os problema */
	  if (calcCoefs (&coef1, &coef2, &coef3, psess->std, &errStd) == FAIL)
	    return FAIL;
	  calcPrbSiz (coef1, coef2, coef3, pprob);
	}
      for (i = 0; i < pprob->num; i++)
	{
	  /* pero s�lo guardamos las movilidades */
	  fprintf (pf, "%g\t%g\n", pprob->mov[i], pprob->siz[i]);
	}
      pprob = pprob->next;
    }
  fclose (pf);
  psess->saved = TRUE;
  return SUCCESS;
}

/*---------------------------------------------------------------------*/
