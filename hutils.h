/* hutils.h  - Header file for hutils.c    


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
   double trapint(double, double, double, double (*)(double));
   double fact(int);
   int signone(int);
   double corr(double,double [],double []);
   void inv(double [],double []);
   void Date(char *);
   void Time(char *);
   double trapint(double, double, double, double (*)(double));

