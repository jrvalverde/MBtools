#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include <grx20.h>
#include <grxkeys.h>


#define MAX_X	400
#define MAX_Y	250

mapa_logistico()
{
  	double x, a;
  	int i, cx, cy;
	int cnt;
  	
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
	
	cnt = 0;
	do
	 {
	   cnt++;
	   for (i = 1; i <= 200; i++)
	     x = a * x * (1 - x);
	   for (i = 1; i <= 300; i++)
	     {
	       x = a * x * (1 - x);
	       cx = floor(MAX_X * (a - 2.95));
	       cy = floor(MAX_Y * x);
	       GrPlot(cx, cy, GrBlack());
	     }
	   a += 2.5e-3;
	 }
	while ((a <= 3.95) && (cx < MAX_X));
  }
  
double f(y, c)
double y, c;
  { 
  	return (c * y + 2 * (1 - c) * (y * y) / (1 + y * y));
  }

/*
A formula is named after the two Physicists (or mathematicians) that discovered it. Their names are Gumowski and Mira.  They did experiments at the CERN research facility in Geneva Switzerland.  They were trying to calculate (or simulate) the trajectories of elementary particles like protons that move at high speeds in an accelerator, a circular channel with the diameter of a tin can but several meters long. Gumowski and Mira used the formula to simulate the orbits of the particles.

One of the formulas that they used is as follows:

X(n+1) = B*Y + F(X)
Y(n+1) = -X + F(X(n+1))

The function "F" is considered their model.  One of the main functions they used for their model is:

F(X) = A*X + (1-A)*2*X*X/(1+X*X)

in which "A" is a parameter to be chosen freely and is usually anywhere in the range of -1 to 1. "B" is a very sensitive constant and usually stays at a value that is very close to 1.000.  If the constant "B" is slightly increased to a value of 1.001, then the trajectory will usually expand (or spiral outward to infinity). If the constant "B" is slightly decreased to something like a value of "0.999", then the trajectory will contract (or spiral inward) towards the attractor points.

Another very important parameter is the initial starting point of the trajectory. This is sometimes called the "Seed" value.  The seed is represented by the initial values of X and Y.  Typical initial values of X and Y can be somewhere in the range of -20 to 20.

For "A=1", F(x) reduces to X. The model (F) is then completly regular, all is stable, and chaos is absent. As A deviates more from 1, more chaos will occur and the computer pictures get more interesting.

In summary, the formula can be written as follows:

X(n+1) = B*Y + F(n)
Y(n+1) = -X + F(X(n+1))

where: F(n) = A*X + (1 - AA)*2*X*X) / (1 + Y*Y)

For programming purposes, this formula can be written as follows:

T = X
X = B*Y + W
W = A*X + (1 - A)*2*X*X/(1 + X*X)
Y = W - T

The trajectories are then plotted on a 2 dimensional graph. The graph coordinates usually range from about -20 to 20 in the X axis and -20 to 20 in the Y axis. 
*/
Gumowski_and_Mira()
  {
    int i, j,
    	cx, cy;
    double c, x, y, xx;
    int cnt;
    extern int grdev;
    
/*
    g2_string(grdev, 30, 180, "MAPA DE GUMOWSKI & MIRA");
    g2_string(grdev, 30, 160, "=======================");
    g2_string(grdev, 20, 140, "x\' <- y");
    g2_string(grdev, 20, 120, "y\' <- -x + 2 F(y)");
    g2_string(grdev, 20, 100, "F(y) <- c.y + 2.(1 - c).y.y / (1 + y.y)");
    g2_string(grdev, 20, 80, "-1 <= c <= 1, -20 <= x, y <= 20");
    g2_string(grdev, 20, 60, "Valor de c: "); scanf("%f", &c);
    g2_string(grdev, 20, 40, "Valor inicial de x: "); scanf("%f", &x);
    g2_string(grdev, 20, 20, "Valor inicial de y: "); scanf("%f", &y);
*/
/*    c = 0.9; x = 10; y = 10; */
    c = -0.662358470; x = 14.717684295; y = -13.660209039;
    
    cnt = 0;
    do
      {
        cnt++;
        xx = y;
        y = -x + 2 * f(y, c);
        x = xx;
        cx = floor(x * 50) + (MAX_X >> 1);
        cy = floor(y * 50) + (MAX_Y >> 1);
        GrPlot(cx, cy, GrBlack());
      }
    while (cnt < 100000);
  }


atractor_de_Henon()
{
  double a, b, c, 
         x, xx, y;
  int cx, cy;
  int cnt;
  GrColor blue;
  extern int grdev;


    GrSetRGBcolorMode();
    /* RGB colors are in 0xRRGGBB format */
    GrSetColor(blue, 0, 0, 128);
  
  a = 3.1678;
  b = 0.3;
  y = x = 1;
  cnt = 0;
  do
    {
      cnt++;
      xx = y;
      y = b * x + a * y - (y * y);
      x = xx;
      cx = floor((y * 100) +50);
      cy = MAX_Y - floor((x * 62.5) + 30);
      GrPlot(cx, cy, blue);
    }
  while (cnt < 100000);
 }
 
 main()
   {
     char ch;
     extern int grdev;
     int xdev, pdev;
     
    GrSetDriver("xwin");

/*    GrSetMode(GR_default_graphics); */
    GrSetMode(GR_width_height_color_graphics, MAX_X, MAX_Y, (GrColor)256);
    
    GrClearScreen(GrWhite());
     
    mapa_logistico();
    GrKeyRead();

    GrClearScreen(GrWhite());
    atractor_de_Henon();
    GrKeyRead();
     
    GrClearScreen(GrWhite());
    Gumowski_and_Mira();
    puts("Done.");
    GrKeyRead();
     
/*     do {
     mapa_logistico();
     ch = getchar();
     atractor_de_Henon();
     ch = getchar(); 
     Gumowski_and_Mira();
     ch = getchar();
     } while ((ch = getchar()) != 'f'); */
   }
 
