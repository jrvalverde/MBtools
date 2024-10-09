/*------------------------------------------------------------------*
*                                                                   *
*    SIGMA.C                                                        *
*                                                                   *
*   Programa para el cálculo de las funciones sigma que determinan  *
*   distintos ángulos de giro de la hélice B del DNA de acuerdo a   *
*   las reglas de Calladine-Dickerson.                              *
*                                                                   *
*    Caveats: El dibujo de las gráficas se realiza con matemática   *
*        entera y presenta problemas de redondeo.                   *
*        Secuencias muy largas no son bien discernibles con éste    *
*        sistema de ventanas.                                       *
*                                                                   *
*    Utiliza los ficheros de cabecera STDIO.H, CTYPE.H,UNIX.H y     *
*            MATH.H                                                 *
*        y las librerías MATH, PRLINK.LIB, STDIO, UNIX y MACTRAPS   *
*                                                                   *
*    Notas:                                                         *
*        - Algunos identificadores de función pueden presentar      *
*        problemas de portabilidad: Calc_Sigma?(), dibuja_sigma?() y*
*        dibuja_ejes(). En éste caso se sugiere sustituirlos por    *
*        c_sig_?(), dibsig?() y dib_ejes().                         *
*        - No lee secuencias de un fichero.                         *
*        - No permite manipulación detallada de fragmentos de una   *
*        secuencia introducida.                                     *
*                                                                   *
*    Véase:                                                         *
*        Shapiro et al. Nucl Acids Res, 14: 75-86                   *
*        Tung & Harvey, Nucl Acids Res, 12: 3343-3356               *
*        Tung & Harvey, Nucl Acids Res, 14: 381-387                 *
*                                                                   *
*   Versión en C para el MacIntosh ajustada al estándar UNIX.       *
*                                                                   *
*   © José Ramón Valverde Carrillo, 1988.                           *
*                                                                   *
*------------------------------------------------------------------*/
    
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#ifdef __LIGHTSPEED_C__
#include "unix.h"
#include "Lightspeed C:Mac #includes:MacTypes.h"
#include "Lightspeed C:Mac #includes:QuickDraw.h"
#include "Lightspeed C:Mac #includes:WindowMgr.h"
#include "Lightspeed C:Mac #includes:FontMgr.h"
#else
#include <unistd.h>
#include <curses.h>
#include <stdlib.h>

#include <g2.h>
#include <g2_PS.h>
#include <g2_FIG.h>
#include <g2_X11.h>

#endif

/* Dimensión máxima de la secuencia a analizar */

#define NOPBASES    1000
#define NOENLACES    1000

/* #defines para las gráficas */

#define MINX    0
#define MINY    0
#define MAXX    512
#define MAXY    342

#ifdef _G2_H
static int grdev;  /* the virtual device we'll use for drawing */

#define move mymove

void mymove(int x, int y)
{
    extern int grdev;
    
    g2_move(grdev, x, y);
}

void line(int x1, int y1, int x2, int y2)
{
    extern int grdev;
    
    g2_line(grdev, x1, y1, x2, y2);
}

void labelat(int x, int y, char *text)
{
    extern int grdev;
    
    g2_string(grdev, x, y, text);
}

void eraseplot()
{
    extern int grdev;
    
    g2_clear(grdev);
}

int getche()
{
    return getchar();
}

void initin()
{
    int psdev, xdev, pngdev;
    
    initscr();	    /* curses */

    grdev = g2_open_vd();
    psdev = g2_open_PS("output.ps", g2_A5, g2_PS_land);
    g2_attach(grdev, psdev);
    xdev = g2_open_X11(MAXX, MAXY);
    g2_attach(grdev, xdev);
}

#endif

/* estructura de datos para almacén de los valores computados.
        El cómputo se hace siguiendo a Calladine-Dickerson        */

typedef struct
    {
        short sigma1;
        short sigma2;
        short sigma3;
        short sigma4;
    } sigma;


/* variables globales */

float baseh;
float basev;
float topeh; 
float topev;

char secuencia[NOPBASES] = {'\0'};
sigma enlace[NOENLACES] = {0};

/*-----------------------------------------------------------*/

int lee_secuencia()
/*                                                        *
   Devuelve el número de bases leídas (el índice de la
   última será este número menos uno).
   Permite edición sencilla de los datos introducidos.
   Convierte entrada a minúsculas.
*                                                        */
{
    register int i;
    register char ch;

    i = 0;
    while (i < NOPBASES)
        {
            /* Leer siguiente letra */
            ch = getche();
            if (ch == '\n')
            {
                /* Es el fín */
                secuencia[i] = '\0';
                break;
            }
            if (ch == '\b')
            {
                /* Corrección */
                if (i == 0)
                  {
                      break;
                  }
                  --i;
                  putchar(' '); putchar('\b');
                  continue;
             }
             if (isalpha(ch))
               ch = (isupper(ch))? tolower(ch) : ch;
             if ((ch != 'a') && (ch != 't')
                && (ch != 'g') && (ch != 'c'))
                {
                    /* letra errónea */
                    putchar('\b');
                    putchar(' ');
                    putchar('\b');
                    putchar('\7');
                    continue;
                }
             secuencia[i] = ch;
             ++i; 
            
        }
    return(i);
}

/*-------------------------------------------------------------*/

Calc_Sigma1(secuencia, no_bases, enlace)
char *secuencia;
unsigned int no_bases;
sigma *enlace;
/*                                                                     *
   Calcula los valores de sigma 1 (rotación media por par de bases)
   correspondientes a la secuencia en estudio. Los valores calculados
   son enteros que oscilan entre [-4, +3].
   Calcula sigma 1 para no_bases bases a partir de la indicada por
   secuencia y almacena el resultado en una secuencia de valores a
   partir del enlace indicado.
   Cada unidad equivale a 2.1° centrados a 35.6°, que son los valores
   apropiados para una estructura B-DNA.
*                                                                    */
{
    register int i;
    
    /* Primer par de bases es un caso especial */
    if ((secuencia[0] == 'g') || (secuencia[0] == 'a'))
    {    /* x-R-?-x */
        if ((secuencia[1] == 't') || (secuencia[1] == 'c'))
        {    /* x-R-Y-x */
            enlace[0].sigma1 = -2;
            enlace[1].sigma1 = 1;
        }
        /* else x-R-R-x sigma1 = 0 */
    }
    else
    {    /* x-Y-?-x */
        if ((secuencia[1] == 'g') || (secuencia[1] == 'a'))
        {    /* x-Y-R-x */
            enlace[0].sigma1 = -4;
            enlace[1].sigma1 = 2;
        }
        /* else x-Y-Y-x sigma1 = 0 */
    }
    
    /* El grueso del peloton */
    for (i = 1; i < no_bases; ++i)
    {
        if ((secuencia[i] == 'g') || (secuencia[i] == 'a'))
        {    /* x-R-?-x */
            if ((secuencia[i + 1] == 't') || (secuencia[i + 1] == 'c'))
            {    /* x-R-Y-x */
                enlace[i - 1].sigma1 += 1;
                enlace[i].sigma1     -= 2;
                enlace[i + 1].sigma1 += 1;
            }
            /* else x-R-R-x sigma1 += 0 */
        }
        else
        {    /* x-Y-?-x */
            if ((secuencia[i + 1] == 'g') || (secuencia[i + 1] == 'a'))
            {    /* x-Y-R-x */
                enlace[i - 1].sigma1 += 2;
                enlace[i].sigma1     -= 4;
                enlace[i + 1].sigma1 += 2;
            }
            /* else x-Y-Y-x sigma1 += 0 */
        }
    }
    
    /* Ultima base es un caso especial */
    /* Su valor ya está calculado, pero sobra el último valor
        calculado.*/
    enlace[no_bases + 1].sigma1 = 0;
}

/*------------------------------------------------------------*/

Calc_Sigma2(secuencia, no_bases, enlace)
char *secuencia;
int no_bases;
sigma *enlace;
/*                                                                        *
   Calcula los valores de sigma 2 (desviación de un eje normal a
   la base respecto al eje de la hélice) para la secuencia en estudio.
   Los valores obtenidos oscilan entre [-5, +7].
   Cada unidad significa 1.1° centrados a -1.0° que son los parámetros
   apropiados para un dúplex de B-DNA.
*                                                                        */
{
    register int i;
    
    /* Primer par de bases es un caso especial */
    if ((secuencia[0] == 'g') || (secuencia[0] == 'a'))
    {    /* x-R-?-x */
        if ((secuencia[1] == 't') || (secuencia[1] == 'c'))
        {    /* x-R-Y-x */
            enlace[0].sigma2 = -2;
            enlace[1].sigma2 = 1;
        }
        /* else x-R-R-x sigma2 = 0 */
    }
    else
    {    /* x-Y-?-x */
        if ((secuencia[1] == 'g') || (secuencia[1] == 'a'))
        {    /* x-Y-R-x */
            enlace[0].sigma2 = 4;
            enlace[1].sigma2 = -2;
        }
        /* else x-Y-Y-x sigma2 = 0 */
    }
    
    /* El grueso del peloton */
    for (i = 1; i < no_bases; ++i)
    {
        if ((secuencia[i] == 'g') || (secuencia[i] == 'a'))
        {    /* x-R-?-x */
            if ((secuencia[i + 1] == 't') || (secuencia[i + 1] == 'c'))
            {    /* x-R-Y-x */
                enlace[i - 1].sigma2 += 1;
                enlace[i].sigma2     -= 2;
                enlace[i + 1].sigma2 += 1;
            }
            /* else x-R-R-x sigma2 += 0 */
        }
        else
        {    /* x-Y-?-x */
            if ((secuencia[i + 1] == 'g') || (secuencia[i + 1] == 'a'))
            {    /* x-Y-R-x */
                enlace[i - 1].sigma2 -= 2;
                enlace[i].sigma2     += 4;
                enlace[i + 1].sigma2 -= 2;
            }
            /* else x-Y-Y-x sigma2 += 0 */
        }
    }
    
    /* Ultima base es un caso especial */
    /* Su valor queda calculado, pero sobra el
        valor del enlace siguiente. */
    enlace[no_bases + 1].sigma2 = 0;
}

/*----------------------------------------------------------------*/
Calc_Sigma3(secuencia, no_bases, enlace)
char *secuencia;
int no_bases;
sigma *enlace;
/*                                                                    *
   Calcula los valores de sigma 3 (desplazamiento del esqueleto 
   azúcar fosfato hacia las Py) para la secuencia en estudio.
   Los valores obtenidos oscilan entre [-3, +3].
   Cada unidad equivale a 15.6° centrados a 0°
*                                                                    */
{
    register int i;
    int no_pasos;
    
    no_pasos = no_bases - 1;
    for (i = 0; i < no_pasos; ++i)
      {
        if ((secuencia[i] == 'a') || (secuencia[i] == 'g'))
          { /* Si la base en curso es R */
            if ((secuencia[i + 1] == 't') || (secuencia[i + 1] == 'c'))
              {  /* R-Y */
                enlace[i].sigma3     += 1;
                enlace[i + 1].sigma3 -= 1;
              }
            /* else R-R : quedan como estaban */
          }
        else
          { /* la base en curso es Y */
            if (( secuencia[i + 1] == 'a') || (secuencia[i + 1] == 'g'))
              {   /* Y-R */
                enlace[i].sigma3     -= 2;
                enlace[i + 1].sigma3 += 2;
              }
            /* else Y-Y : quedan intactos */
          }
      }
}

/*----------------------------------------------------------------*/
Calc_Sigma4(secuencia, no_bases, enlace)
char *secuencia;
int no_bases;
sigma *enlace;
/*                                                                    *
   Calcula los valores de sigma 4 (giro de propulsión para cada
   par de bases) para la secuencia en estudio.
   Los valores obtenidos oscilan entre [-3, 0].
   Cada unidad equivale a -3.6° centrados a 24.3°
*                                                                    */
{
    register int i, no_pasos;
    
    no_pasos = no_bases - 1;
    for (i = 0; i < no_pasos; ++i)
      {
        if ((secuencia[i] == 'a') || (secuencia[i] == 'g'))
          {  /* Si la base en curso es R */
            if ((secuencia[i + 1] == 't') || (secuencia[i + 1] == 'c'))
              {  /* R-Y */
                enlace[i].sigma4     -= 1;
                enlace[i + 1].sigma4 -= 1;
              }
            /* else R-R : quedan como estaban */
          }
        else
          { /* la base en curso es Y */
            if (( secuencia[i + 1] == 'a') || (secuencia[i + 1] == 'g'))
              {   /* Y-R */
                enlace[i].sigma4     -= 2;
                enlace[i + 1].sigma4 -= 2;
              }
            /* else Y-Y : quedan intactos */
          }
      }
}

/*----------------------------------------------------------------*/

dibuja_ejes()
/* Dibuja unas cuantas líneas que sirven de ejes  */
{
    int i, x1, x2, y1, y2;
    
    /* eje horizontal: las y-es permanecen constantes */
    x1 = baseh; x2 = topeh;
    y1 = y2 = basev;
    line(x1, y1, x2, y2);
    y1 = y2 = topev; 
    line(x1, y1, x2, y2);
    
    /* eje vertical: las x-s permanecen constantes */
    x1 = x2 = baseh;
    y1 = basev; y2 = topev;
    line(x1, y1, x2, y2);
    x1 = x2 = topeh;
    line(x1, y1, x2, y2);
    
    /* marcas, tanto en el eje X como en el Y */
    for (i = 0; i <=10; ++i)
      {
        y1 = y2 = basev - (((basev - topev) / 10) * i);
        x1 = baseh; x2 = x1 + 5;
        line(x1, y1, x2, y2);
        x1 = x2 = baseh + (((topeh - baseh) / 10) * i);
        y1 = basev; y2 = y1 - 5;
        line(x1, y1, x2, y2);
      }
      
}

/*----------------------------------------------------------------*/

dibuja_sigma1(inicio, num)
int inicio,            /* nucleótido de partida */
    num;            /* número de nucleótidos a considerar */
/*                                                                *
  Dibuja una gráfica que representa los valores de sigma 1 para
  los enlaces correspondientes a los nucleótidos indicados 
*                                                                */
{ 
    float inch, incv, xnueva, ynueva, xvieja, yvieja;
    int x1, y1, x2, y2;
    register int i;
    
    dibuja_ejes();
    x1 = baseh;
    y1 = topev;
    move(x1, y1); 
    label("SIGMA 1: cómputo del ángulo de torsión");
    
    inch = (topeh - baseh) / num;
    incv = (basev - topev) / 8.0;
    
    xvieja = baseh;
    yvieja = basev - ((enlace[inicio].sigma1 + 4.0) * incv);
    for (i = 1; i < num; ++i)
      {
        xnueva = xvieja + inch;
        ynueva = basev - ((enlace[inicio + i].sigma1 + 4.0) * incv);
        x1 = floor(xvieja);
        y1 = floor(yvieja);
        x2 = floor(xnueva);
        y2 = floor(ynueva);
        line(x1, y1, x2, y2);
        xvieja = xnueva;
        yvieja = ynueva;
      }
}

/*----------------------------------------------------------------*/

dibuja_sigma2(inicio, num)
int inicio,            /* nucleótido de partida */
    num;            /* número de nucleótidos a considerar */
/*                                                                *
   Dibuja una gráfica que representa los valores de sigma 2 para
   los enlaces correspondientes a los nucleótidos indicados
*                                                                */
{ 
    float inch, incv, xnueva, ynueva, xvieja, yvieja;
    int x1, y1, x2, y2;
    register int i;
    
    dibuja_ejes();
    x1 = baseh;
    y1 = topev;
    move(x1, y1); 
    label("SIGMA 2: cómputo del ángulo de rotación de los pares de bases");
    
    inch = (topeh - baseh) / num;
    incv = (basev - topev) / 12;
    
    xvieja = baseh;
    yvieja = basev - ((enlace[inicio].sigma2 + 6) * incv);
    for (i = 1; i < num; ++i)
      {
        xnueva = xvieja + inch;
        ynueva = basev - ((enlace[inicio + i].sigma2 + 6) * incv);
        x1 = floor(xvieja);
        y1 = floor(yvieja);
        x2 = floor(xnueva);
        y2 = floor(ynueva);
        line(x1, y1, x2, y2);
        xvieja = xnueva;
        yvieja = ynueva;
      }
}

/*----------------------------------------------------------------*/

dibuja_sigma3(inicio, num)
int inicio,            /* nucleótido de partida */
    num;            /* número de nucleótidos a considerar */
/*                                                                     *
   Dibuja una gráfica que representa los valores de sigma 3 para
   los enlaces correspondientes a los nucleótidos indicados 
*                                                                    */
{ 
    float inch, incv, xnueva, ynueva, xvieja, yvieja;
    int x1, y1, x2, y2;
    register int i;
    
    dibuja_ejes();
    x1 = baseh;
    y1 = topev;
    move(x1, y1); 
    label("SIGMA 3: cómputo del ángulo de desplazamiento hacia las Py");
    
    inch = (topeh - baseh) / num;
    incv = (basev - topev) / 6;
    
    xvieja = baseh;
    yvieja = basev - ((enlace[inicio].sigma3 + 3) * incv);
    for (i = 1; i < num; ++i)
      {
        xnueva = xvieja + inch;
        ynueva = basev - ((enlace[inicio + i].sigma3 + 3) * incv);
        x1 = floor(xvieja);
        y1 = floor(yvieja);
        x2 = floor(xnueva);
        y2 = floor(ynueva);
        line(x1, y1, x2, y2);
        xvieja = xnueva;
        yvieja = ynueva;
      }
}

/*----------------------------------------------------------------*/

dibuja_sigma4(inicio, num)
int inicio,            /* nucleótido de partida */
    num;            /* número de nucleótidos a considerar */
/*                                                                    *
   Dibuja una gráfica que representa los valores de sigma 4 para
   los enlaces correspondientes a los nucleótidos indicados 
*                                                                    */
{ 
    float inch, incv, xnueva, ynueva, xvieja, yvieja;
    int x1, y1, x2, y2;
    register int i;
    
    dibuja_ejes();
    x1 = baseh;
    y1 = topev;
    move(x1, y1); 
    label("SIGMA 4: cómputo del ángulo de propulsión para cada par de bases");
    
    inch = (topeh - baseh) / num;
    incv = (basev - topev) / 5;
    
    xvieja = baseh;
    yvieja = basev - ((enlace[inicio].sigma4 + 4) * incv);
    for (i = 1; i < num; ++i)
      {
        xnueva = xvieja + inch;
        ynueva = basev - ((enlace[inicio + i].sigma4 + 4) * incv);
        x1 = floor(xvieja);
        y1 = floor(yvieja);
        x2 = floor(xnueva);
        y2 = floor(ynueva);
        line(x1, y1, x2, y2);
        xvieja = xnueva;
        yvieja = ynueva;
      }
}

/*----------------------------------------------------------------*/

dib_sig1_2(inicio, num)
int inicio,            /* nucleótido de partida */
    num;            /* número de nucleótidos a considerar */
/*                                                                *
  Dibuja una gráfica que representa los valores de sigma 1 al
  cuadrado para los enlaces correspondientes a los nucleótidos
  indicados 
*                                                                */
{ 
    float inch, incv, xnueva, ynueva, xvieja, yvieja;
    int x1, y1, x2, y2;
    register int i;
    
    dibuja_ejes();
    x1 = baseh;
    y1 = topev;
    move(x1, y1); 
    label("SIGMA 1: cómputo del ángulo de torsión");
    
    inch = (topeh - baseh) / num;
    incv = (basev - topev) / 64.0;
    
    xvieja = baseh;
    yvieja = basev - ((enlace[inicio].sigma1 * enlace[inicio].sigma1
            + 16.0) * incv);
    for (i = 1; i < num; ++i)
      {
        xnueva = xvieja + inch;
        ynueva = basev - ((enlace[inicio + i].sigma1 * enlace[inicio].sigma1
                + 16.0) * incv);
        x1 = floor(xvieja);
        y1 = floor(yvieja);
        x2 = floor(xnueva);
        y2 = floor(ynueva);
        line(x1, y1, x2, y2);
        xvieja = xnueva;
        yvieja = ynueva;
      }
}

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/

main()
/* Programa principal. */
{
    char ch;
    int i;
    int no_datos, old_cx, old_cy, cx, cy, inch, incv;
    
    initin();
    
    /* Lee secuencia */
    eraseplot();
    puts("Poh' favó, meta la zecuencia pa 'analizá\n");
    puts("Aquí 'ebajo\n");
    no_datos = lee_secuencia();
    
    /* Calcula valores de sigma para la secuencia introducida */
    Calc_Sigma1(&secuencia[0], no_datos, &enlace[0]);
    Calc_Sigma2(&secuencia[0], no_datos, &enlace[0]);
    Calc_Sigma3(&secuencia[0], no_datos, &enlace[0]);
    Calc_Sigma4(&secuencia[0], no_datos, &enlace[0]);
    
    /* Escribe los valores de sigma */
    puts("\n\n");
    puts("Valores calculados de Sigma 1");
    puts("=============================");
    for (i = 0; i < no_datos; i++)
        printf("%g° ", ((enlace[i].sigma1 * 2.1) + 35.6));
    puts("\n\n");
    puts("valores calculados de Sigma 2");
    puts("=============================");
    for (i = 0; i < no_datos; i++)
        printf("%g° ", ((enlace[i].sigma2 * 1.1) -1.0));
    puts("\n\n");
    puts("Valores calculados de Sigma 3");
    puts("=============================");
    for (i = 0; i < no_datos; i++)
        printf("%g° ", (enlace[i].sigma3 * 15.6));
    puts("\n\n");
    puts("Valores calculados de Sigma 4");
    puts("=============================");
    for (i = 0; i < no_datos; i++)
        printf("%g° ", ((enlace[i].sigma4 * -3.6) + 24.3));
    
    /* Nota de aviso */
    puts("\n");
    puts("Siempre: pulsar [c] para continuar");
    while ((ch = getch()) != 'c');
    eraseplot();
    
    eraseplot();
    /* Gráfica de Sigma 1 */
    baseh = 10.0;  basev = 135.0;
    topeh = 475.0; topev = 10.0;
    
    dibuja_sigma1(0, no_datos);
    while ((ch = getch()) != 'c');
    
    /* Gráfica de sigma 2 */
    baseh = 10.0; basev = 270.0;
    topeh = 475.0; topev = 145.0;
    
    dibuja_sigma2(0, no_datos);
    while ((ch = getch()) != 'c');
    
    eraseplot();
    
    /* Gráfica de sigma 3 */
    baseh = 10.0;  basev = 135.0;
    topeh = 475.0; topev = 10.0;
    
    dibuja_sigma3(0, no_datos);
    while ((ch = getch()) != 'c');
    
    /* Gráfica de sigma 4 */
    baseh = 10.0; basev = 270.0;
    topeh = 475.0; topev = 145.0;
    
    dibuja_sigma4(0, no_datos);
    while ((ch = getch()) != 'c');
    
    eraseplot();
    /* Gráfica de sigma 1 ** 2 */
    baseh = 10.0;  basev = 135.0;
    topeh = 475.0; topev = 10.0;
    
    dib_sig1_2(0, no_datos);
    while ((ch = getch()) != 'c');
    
    /* Adios, muy buenas */
    eraseplot();
    puts("");
    move(200, 200);
    label("Adios, muy buenas");
    exit(1);
}

/*------------------------ Fin del mundo -------------------------*/
