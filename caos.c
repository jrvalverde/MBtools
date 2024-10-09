#ifndef _stdioh_
#include <stdio.h>
#endif

#ifndef _unixh_
#include <unix.h>
#endif

#ifndef _mathh_
#include <math.h>
#endif

#define MAX_X	400
#define MAX_Y	250

mapa_logistico(erase)
int erase;
{
  	double x, a;
  	int i, cx, cy;
  	
  /*
  	Vamos a representar el mapa de un oscilador ca—tico
  	conocido como mapa log’stico.
  	El oscilador es
  		dx / dt = a.x.(1 - x),   0Ê< a <= 4
  	y veremos como van variando los valores de x
  	conforme vamos variando a.
  */
	/* Empezamos con a = 2.95 */
	a = 2.95;
	x = 0.3;
	if (erase == 1)
	  eraseplot();
	do
	 {
	   for (i = 1; i <= 200; i++)
	     x = a * x * (1 - x);
	   for (i = 1; i <= 300; i++)
	     {
	       x = a * x * (1 - x);
	       cx = floor(MAX_X * (a - 2.95));
	       cy = floor(MAX_Y * x);
	       point(cx, cy);
	     }
	   a += 2.5e-3;
	 }
	while ((a <= 3.95) && (! kbhit()));
	putchar('\007');
	while (! kbhit());
  }
  
double f(y, c)
double y, c;
  { 
  	return (c * y + 2 * (1 - c) * (y * y) / (1 + y * y));
  }
  
Gumowski_and_Mira()
  {
    int i, j,
    	cx, cy;
    double c, x, y, xx;
    
    eraseplot();
    gotoxy(30, 2);
    puts("MAPA DE GUMOWSKI & MIRA");
    gotoxy(30, 3);
    puts("=======================");
    gotoxy(20, 5);
    puts("x\' <- y");
    gotoxy(20, 6);
    puts("y\' <- -x + 2 F(y)");
    gotoxy(20, 7);
    puts("F(y) <- c.y + 2.(1 - c).y.y / (1 + y.y)");
    gotoxy(10, 10);
    printf("Valor de c: "); scanf("%f", &c);
    gotoxy(10, 15);
    printf("Valor inicial de x: "); scanf("%f", &x);
    gotoxy(10, 20);
    printf("Valor inicial de y: "); scanf("%f", &y);
    eraseplot();
    do
      {
        xx = y;
        y = -x + 2 * f(y, c);
        x = xx;
        cx = floor(x * 50) + (MAX_X >> 1);
        cy = floor(y * 50) + (MAX_Y >> 1);
        point(cx, cy);
      }
    while (! kbhit());
    putchar('\007');
    while (! kbhit());
  }
  

atractor_de_Henon()
{
  double a, b, c, 
         x, xx, y;
  int cx, cy;
  
  eraseplot();
  a = 3.1678;
  b = 0.3;
  y = x = 1;
  do
    {
      xx = y;
      y = b * x + a * y - (y * y);
      x = xx;
      cx = floor((y * 100) + 100);
      cy = MAX_Y - floor((x * 62.5) + 30);
      point(cx, cy);
    }
  while (! kbhit());
  putchar('\007');
  while(! kbhit());
 }
 
 main()
   {char ch;
     do {
     atractor_de_Henon();
     ch = getche();
     Gumowski_and_Mira();
     ch = getche();
     } while ((ch = getche()) != 'f');
   }
 