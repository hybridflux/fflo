/*   hutils.c - utilility routines for Hc2-program                         

Copyright (C) by Sherryl Manalo, 

functions inv, corr, Date, Time are (C) of 
Ewald Schachinger <schachinger@itp.tu-graz.ac.at>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Contact the author at: sm@antitachyon.com
	
*/

#include <stdio.h>
#include  "hutils.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>


/********** calculates trapezoidal integration over function ****************/
   double trapint(double beg, double end, double mesh, double (*fun)(double))
   /****************************************************************************/
   { int i;
      double dy,y=beg,res=0.0;  
   
      dy = (end-beg)/mesh;
   
   
      for (i = 0; i < mesh; i++)
      { 
                
         res += 0.5*((*fun)(y)+(*fun)(y+dy))*dy;
/*         printf("%f\t %E\t %E\n",y, (*fun)(y),res); */ 
      	y += dy;
      
      }
   
      return res;
   }



/***************** calculates factorial product of integer ******************/
   double fact(int m)  
   /****************************************************************************/
   { 	int j;
      double fac=1.0;
      for (j=m; j>0; j--)
      {
         fac *= j;
      }
   
      return fac;
   }

/******************** calculates (-1)^m *************************************/
   int signone(int m)
   /****************************************************************************/
   {	
      if (!fmod(m,2))
         return 1; 
   
      return -1;
   }

/*****************************************************************************/
   double corr (double v,double x[],double y[])
   /*****************************************************************************/
   {
      double s,x21,x31,y21,y31,xq21,xq31,c,b;
   
      x21 = x[1]-x[0];
      x31 = x[2]-x[0];
      y21 = y[1]-y[0];
      y31 = y[2]-y[0];
      xq21 = x[1]*x[1]-x[0]*x[0];
      xq31 = x[2]*x[2]-x[0]*x[0];
      c = (y31*x21-y21*x31)/(xq31*x21-xq21*x31);
      b = (y21-c*xq21)/x21;
      s = b+2.0*c*x[1];
      return ((y[1]-v)/s);
   }

/****************************************************************************/
   void inv (double a[],double b[])
   /****************************************************************************/
   {
      double x1,x2,x3,x4,t1,t12,t13,t2,t22,t23,t3,t32,t33,t4,t42,t43,d1,d2;
   
      x1 = b[0];
      x2 = b[1];
      x3 = b[2];
      x4 = b[3];
      t1 = a[0];
      t12 = t1*t1;
      t13 = t12*t1;
      t2 = a[1];
      t22 = t2*t2;
      t23 = t22*t2;
      t3 = a[2];
      t32 = t3*t3;
      t33 = t32*t3;
      t4 = a[3];
      t42 = t4*t4;
      t43 = t42*t4;
      d1 = ((t2-t1)*(t43-t13)-(t4-t1)*(t23-t13))*
         ((t2-t1)*(t32-t12)-(t3-t1)*(t22-t12))-
         ((t2-t1)*(t33-t13)-(t3-t1)*(t23-t13))*
         ((t2-t1)*(t42-t12)-(t4-t1)*(t22-t12));
      d2 = ((x4-x1)*(t2-t1)-(x2-x1)*(t4-t1))*
         ((t2-t1)*(t32-t12)-(t3-t1)*(t22-t12))-
         ((x3-x1)*(t2-t1)-(x2-x1)*(t3-t1))*
         ((t2-t1)*(t42-t12)-(t4-t1)*(t22-t12));
      b[3] = d2/d1;
      b[2] = ((x3-x1)*(t2-t1)-(x2-x1)*(t3-t1)-
         b[3]*((t2-t1)*(t33-t13)-(t3-t1)*(t23-t13)))/
         ((t2-t1)*(t32-t12)-(t3-t1)*(t22-t12));
      b[1] = (x2-x1-b[2]*(t22-t12)-b[3]*(t23-t13))/(t2-t1);
      b[0] = x1-b[1]*t1-b[2]*t12-b[3]*t13;
   }

/*****************************************************************************/
   void Date(char *s)
   /*****************************************************************************/
   {
      time_t clock;
      struct tm *t;
   
      time(&clock);
      t = localtime(&clock);
      sprintf(s,"%02d-%02d-%4d",t->tm_mday,t->tm_mon+1,t->tm_year+1900);
   }
/*****************************************************************************/
   void Time(char *s)
   /*****************************************************************************/
   {
      time_t clock;
      struct tm *t;
   
      time(&clock);
      t = localtime(&clock);
      sprintf(s,"%02d:%02d:%02d",t->tm_hour,t->tm_min,t->tm_sec);
   }

