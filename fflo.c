/*
fflo - calculates the fflo critical field of superconductors 
Copyright (C) by Sherryl Manalo

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
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "hutils.h"
#include "hgauss.h"

#define  DIM        100000
#define  CON	    0.561459295424
#define  SDIM       35
#define  STR32      32

   double fun(double);
   double hcDelt();
   void iterD(int);
   void initV();
   void Page1();
   void calcH(int);
   void calcgap(int);
   void extra();
   void calcs(int,double);

   double sm,hkapt,summ[DIM],TOL,argsin,lntdtc,mesh,argfun,argcos,stepsize,
   hsav,sdamp,mucon,argcos,mufac,mudamp,TOLOR,argqnow,argq,qstart,
   argqdif,qdamp,dfak,TOLq,hsav,d_angle1;

   float s_tc[DIM];
   double hc_2[DIM],d_angle=0.0,qstep;

   int ntemp,maxiter,calcgapv,rangehmax,calcmu,INITIAL=0,MODUS;

   FILE *orig,*fp,*fil,*gaps,*deri;

   char *dummy,hc2orig[STR32],hc2lst[STR32],hc2gap[STR32],mufile[STR32];

/*****************************************************************************/
   int main(int argc,char *argv[])
   /*****************************************************************************/
   { FILE *ftemp;
      float x[SDIM];
      int N=0,i;
   
      if(argc < 3)   
      {
         printf("\n\n%s Copyright (c) 2001 Sherryl Manalo <sm@antitachyon.com>\n\n",argv[0]);
         printf("%s comes with ABSOLUTELY NO WARRANTY; for details read the\nCOPYING file that comes with this package.\n\n",argv[0]);
       	printf("usage: %s configfile temperaturefile savename\n",argv[0]);
         exit(0);
      }
   
   
      fp = fopen(argv[1],"r");
   
      while (!feof(fp))
      {
         fscanf(fp, "%s  %f",&dummy,&x[N]);        /* read config-file */
         N++;
      }
   
      ftemp = fopen(argv[2],"r");
   
   /*         assign parameters       */  
      sm = x[0];
      mucon = (mufac = x[1]);
      calcmu = x[2]; 
      mudamp = x[3];
      TOLOR = x[4];
      hkapt = (hsav = x[5]);
      ntemp = x[6];
      maxiter = x[7];
      TOL = x[8];
      TOLq = x[9];
      sdamp = x[10];
      mesh = x[11];
      calcgapv = x[12];
      stepsize = x[13];
      rangehmax = x[14];
      qstart = (argqdif = x[15]);
      qdamp = x[16];
      MODUS = x[17];
      d_angle = x[18]*M_PI/180;
		d_angle1 = x[19]*M_PI/180;
      qstep = x[20];
		
   
      strcpy(hc2lst,argv[3]);
      strcpy(hc2orig,argv[3]);
      strcpy(hc2gap,argv[3]);
      strcpy(mufile,argv[3]);
   
      strcat(hc2lst,".lst");
      strcat(hc2orig,".orig");
      strcat(hc2gap,".gap");
      strcat(mufile,".mu");
   
      N = 0;
      while (!feof(ftemp))
      {
         fscanf(ftemp,"%f",&s_tc[N]);             /* read temperature-file */
         N++;
      } 
   
      Page1();
   
      if (calcgapv) {
         initV();
         for (i=0;i<ntemp;i++)
         {
            calcgap(i); }
      
      }
      else {
         iterD(ntemp); 
         extra();
      }
   
      fclose(fp);
      return 0; 
   }

/*****************************************************************************/
   void Page1()
   /*****************************************************************************/
   {  char ddate[11],ttime[11];
   
      Time(ttime);
      Date(ddate);
      fil = fopen(hc2lst,"w");
      fprintf(fil,"Day and Time of Run: %s//%s\n\n\n\n",ddate,ttime);
      fputs("Input Data:\n",fil);
      fprintf(fil,"Parameter sm:          %f\n",sm);
      fprintf(fil,"Parameter mucon:       %f\n",mucon);
      fprintf(fil,"Ist value of H*vf^2:   %f\n",hkapt);
      fprintf(fil,"Maximum iterations:    %i\n",maxiter);
      fprintf(fil,"Precision:             %E\n\n",TOL);
      fclose(fil);
   
   }

/*****************************************************************************/
   void initV()      /*  initialization of matrix-and gapvalues */
   /*****************************************************************************/
   {
      int i;
   
      for (i = 0; i<DIM; i++)
      { summ[i] = 0.0; } 
   
   }

/*****************************************************************************/
   void iterD(int ntem)     /*  iteration of temperatures */
   /*****************************************************************************/
   { int inum,j,itemp;
      double x[3],y[3],newdif,olddif,korr;
      initV();
   
      for (itemp = 0; itemp<ntem; itemp++)
      {       
         printf("Num of temp: %i, current tempnum: %i, t = %f\n\n",ntem,itemp+1,s_tc[itemp]);
         /*  fprintf(fil,"Num of temp: %i, current tempnum: %i, t = %f\n\n",ntem,itemp+1,s_tc[itemp]); */
         if (calcmu == 1) {
            INITIAL = 1;
            calcH(itemp);
            hc_2[itemp]=hkapt;
            INITIAL = 0;
         
            olddif = newdif = 1; 
            for (inum = 0; inum<maxiter && (olddif || newdif); inum++) 
            { 
               olddif = newdif;
               x[0] = mufac==0.0 ? -0.001 : mufac*0.99999;
               x[1] = mufac;
               x[2] = mufac==0.0 ? 0.001 : mufac*1.00001;
               for (j=0; j<=2; j++) {
                  mufac=x[j];
                  calcH(itemp);
                  y[j] = hkapt-hc_2[itemp];
               }
               korr = corr(0.0,x,y);
               newdif = fabs(korr) > TOLOR;
               mufac = (mufac-korr+mudamp*mufac)/(1.0+mudamp);
               printf("\nIt. step: %4i  mufac = %E  change in mufac = %E\n\n",inum,mufac,korr); 
               printf("\nIt. step: %4i  hc2ursp = %E  hc2now = %E\n\n",inum,hc_2[itemp],hkapt); }
            printf("Value of mufac at temp %f: %f\n",s_tc[itemp],mufac);
            deri = fopen(mufile,"a");            
            fprintf(deri,"%f %f\n",s_tc[itemp],mufac);
            fclose(deri);
            orig = fopen(hc2orig,"a");
            fprintf(orig,"%f\t%E\n",s_tc[itemp],hkapt);
            fclose(orig); 
         }
         
         else {
            calcH(itemp); } }
   
   }

/*****************************************************************************/
   void calcH(int run)
   /*****************************************************************************/
   {
      int qinum,j;
      double qq[4],qx[4],qy[4],qkorr,qnewdif,qolddif,qkorrt;
   
   
      qolddif = qnewdif = 1;
   
      for (qinum = 0; qinum<maxiter && (qolddif || qnewdif); qinum++) 
      {
         qolddif = qnewdif;
         qx[0] = argqdif==0.0 ? -0.001 : argqdif*0.999;
         qx[1] = argqdif==0.0 ? 0.0 : argqdif;
         qx[2] = argqdif==0.0 ? 0.001 : argqdif*1.001;
         qx[3] = argqdif==0.0 ? 0.002 : argqdif*1.002;
      
         for (j=0; j<=3; j++) 
         { calcs(run,qx[j]);  
            qy[j] = hkapt; 
            if (j ==1) {hsav = hkapt; }
         }
      
         for(j=0;j<=2;j++) {
            qq[j] = (qy[j+1]-qy[j])/0.001;
         
         }
         qkorr = corr(0.0,qx,qq);
         if (fabs(qkorr) < TOLq)
         { qnewdif = 0; 
            qolddif = 0; }
         argqdif = (argqdif-qkorr+qdamp*argqdif)/(1.0+qdamp);
         printf("\nIt. step: %4i  qkappa-t = %E  change in qkappa = %E\n\n",qinum,argqdif,qkorr);
         if (fabs(qkorr+qkorrt) < 1.0E-15)
         { 
            if (qkorr>0) argqdif = argqdif+qkorr/2.0;
            if (qkorr<0) argqdif = argqdif-qkorr/2.0;
            printf("Oscillation of value! qkorr+qkorrt: %E\n",qkorr+qkorrt);
            qnewdif = 0;
            qolddif = 0;
         }
         qkorr = qkorr+qstep; 
         qkorrt = qkorr;
      
      }
   
      if (!calcmu) {
         orig = fopen(hc2orig,"a");
         fprintf(orig,"%f\t%E\n",s_tc[run],hsav);
         fil = fopen(hc2lst,"a");
         fprintf(fil,"%f\t%E\n",s_tc[run],argqdif);
         fclose(orig); 
         fclose(fil); 
         summ[run-1]=hkapt; } 
   }

/***************************************************************************/
   void calcs(int run,double qxi)
   /***************************************************************************/
   { int i,inum;
      double x[3],y[3],korr,korrt,olddif,newdif;
   
      lntdtc = log(s_tc[run]);
   
   
      olddif = newdif = 1;
   
      for (inum = 0; inum<maxiter && (olddif || newdif); inum++) 
      { 
         olddif = newdif;
         x[0] = hkapt==0.0 ? -0.001 : hkapt*0.999;
         x[1] = hkapt;
         x[2] = hkapt==0.0 ? 0.001 : hkapt*1.001;
      
         for (i=0; i<=2; i++) {
            if (INITIAL) { mufac = 1.0; }	
            argcos = x[i]*CON*mufac/s_tc[run];
            argqnow= qxi/s_tc[run];
            argsin = x[i]*CON*sm/s_tc[run];
            y[i] = hcDelt()+lntdtc;   
         }
      
         korr = corr(0.0,x,y);
         newdif = fabs(korr) > TOL;
         hkapt = (hkapt-korr+sdamp*hkapt)/(1.0+sdamp);
         if (fabs(korr+korrt) < 1.0E-15) 
         { 
            if (korr>0) hkapt = hkapt+korr/2.0; 
            if (korr<0) hkapt = hkapt-korr/2.0;
            printf("Oscillation of value! korr+korrt: %E\n",korr+korrt); 
            newdif = 0; 
            olddif = 0; 
         }
         korrt = korr;
         printf("\nIt. step: %4i  kappa-t = %E  change in kappa = %E <---------\n\n",inum,hkapt,korr);
      
      } 
   
   
   }


/*****************************************************************************/
   void extra()
   /*****************************************************************************/
   { double slope,hc_0;
      slope = (summ[ntemp-1]-summ[ntemp-2])/(s_tc[ntemp-1]-s_tc[ntemp-2]);    
      hc_0 = summ[ntemp-1] - slope*s_tc[ntemp-1]; 
      orig = fopen(hc2orig,"a");
      fprintf(orig,"%f\t%E\n",0.0,hc_0);
      fclose(orig); 
   }
/*****************************************************************************/
   double hcDelt()
   /*****************************************************************************/
   {  int i,uppersin,maxit = 1000,n_int;
      double beg = 0.0, end = 30.0,res = 0.0;
      short error;
   
      uppersin = (int)(1/mesh);
   
   /********** calculation of the mesh-points **************/
      for (i = 1; i<=uppersin; i++)
      { 
         argfun = argsin*(cos(2*i*mesh*M_PI)*cos(d_angle1)-sin(2*i*mesh*M_PI)*sin(d_angle1));
         if (!MODUS)
         {  argq = argqnow*cos(2*i*mesh*M_PI); 
            dfak = 1.0; }
         else if (MODUS==1)
         {  argq = argqnow*cos(2*i*mesh*M_PI); 
            dfak = 1.0 + cos(4*2*i*M_PI*mesh); }
         else 
         {  argq = argqnow*(cos(2*i*mesh*M_PI)*cos(d_angle)+sin(2*i*mesh*M_PI)*sin(d_angle)); 
            dfak = 1.0 + cos(4*2*i*M_PI*mesh); }
      
         gauss(beg,end,TOL,maxit,&summ[i],&n_int,&error,fun);
      }
   
   /********** mesh-summation, summ[uppersin+1] = 0.0 , 2*PI/2PI ********/
      for (i=0; i<uppersin; i++)
      { res += 0.5*mesh*(summ[i+1]+summ[i]); } 
   
      res += 0.5*(1-uppersin*mesh)*summ[uppersin];
      return res;
   }

/*****************************************************************************/
   double fun(double x)
   /*****************************************************************************/
   { 
      if (INITIAL || (sm==0.0))
         return dfak*(1-cos(argcos*x)*cos(argq*x)-sin(argcos*x)*sin(argq*x))/sinh(x);
      else
         return dfak*(1-(cos(argcos*x)*cos(argq*x)+sin(argcos*x)*sin(argq*x))*sin(argfun*x)/(argfun*x))/sinh(x); 
   }
/*****************************************************************************/
   void calcgap(int run)
   /*****************************************************************************/
   { int i;
      double qkapt=qstart;
   
   
      for (i = 0; i<=rangehmax; i++)
      { 
         qkapt = qkapt - stepsize;
         calcs(run,qkapt);
      
         gaps = fopen(hc2gap,"a");
         fprintf(gaps,"%E\t %E\n",qkapt,hkapt);
         printf("%E\t %E\n",qkapt,hkapt); 
         fclose(gaps); 
      
      }
   
   }
