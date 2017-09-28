/*   gauss.c   Adaptive Gauss Quadrature    

Copyright (C) of Ewald Schachinger  <schachinger@itp.tu-graz.ac.at> */

#include  <stdio.h>
#include  <math.h>
#include  <string.h>
#include  <stdlib.h>

#define  LEGENDRE    16

struct i_data {
   double lb,vlb,ub,vub,est;
};
struct s_item {
   struct i_data info;
   struct s_item *next;
};
struct t_legend {
   double root,weight;
};
static struct t_legend tnc[] = {{0.0,0.0},
                 {0.0950125098376370440185,0.189450610455068496285},
                 {0.281603550778258913230,0.182603415044923588867},
                 {0.458016777657227386342,0.169156519395002538189},
                 {0.617876244402643748447,0.149595988816576732081},
                 {0.755404408355003033895,0.124628971255533872052},
                 {0.865631202387831743880,0.095158511682492784810},
                 {0.944575023073232576078,0.062253523938647892863},
                 {0.989400934991649932596,0.027152459411754094852},
                 {-0.0950125098376370440185,0.189450610455068496285},
                 {-0.281603550778258913230,0.182603415044923588867},
                 {-0.458016777657227386342,0.169156519395002538189},
                 {-0.617876244402643748447,0.149595988816576732081},
                 {-0.755404408355003033895,0.124628971255533872052},
                 {-0.865631202387831743880,0.095158511682492784810},
                 {-0.944575023073232576078,0.062253523938647892863},
                 {-0.989400934991649932596,0.027152459411754094852}};
static struct s_item *stack;
static double ga_qua(double,double,double (*)(double));
static void push(struct i_data *),pop(struct i_data *),
       te_ini(double,double,double,int,struct i_data *,double *,int *,
              short *,short *,double (*)(double)),
       add_int(double,double *,int *,int,struct i_data *,short *,short *),
       div_int(double,double,double,double,struct i_data *);

/*****************************************************************************/
static double ga_qua(double lo_lim,double up_lim,double (*fun)(double))
/*****************************************************************************/
{
   double slope,intercept,transroot,sum;
   int t;

   slope = (up_lim-lo_lim)/2.0;
   intercept = (up_lim+lo_lim)/2.0;
   sum = 0.0;
   for (t=1; t<=LEGENDRE; t++) {
      transroot = slope*tnc[t].root+intercept;
      sum += tnc[t].weight*(*fun)(transroot);
   }
   return(sum*slope);
}
/*****************************************************************************/
static void push(struct i_data *data)
/*****************************************************************************/
{
   struct s_item *newnode;

   if ((newnode = (struct s_item *)calloc(1,sizeof(struct s_item))) == NULL) {
      fprintf(stderr,"\nNot enough memory!\n");
      exit(0);
   }
   memcpy((char *)&newnode->info,(char *)data,sizeof(struct i_data));
   newnode->next = stack;
   stack = newnode;
}
/*****************************************************************************/
static void pop(struct i_data *data)
/*****************************************************************************/
{
   struct s_item *oldnode;

   oldnode = stack;
   stack = oldnode->next;
   memcpy((char *)data,(char *)&oldnode->info,sizeof(struct i_data));
   free((char *)oldnode);
}
/*****************************************************************************/
static void te_ini(double lo_lim,double up_lim,double tolerance,int maxint,
            struct i_data *data,double *integral,int *numint,short *finished,
            short *error,double (*fun)(double))
/*****************************************************************************/
{
   *error = 0;
   if (tolerance < 0.0) *error = 1;
   if (maxint < 0) *error = 2;
   if (!*error) {
      data->lb = lo_lim;
      data->vlb = (*fun)(lo_lim);
      data->ub = up_lim;
      data->vub = (*fun)(up_lim);
      data->est = ga_qua(lo_lim,up_lim,fun);
   }
   *integral = 0.0;
   *numint = 0;
   *finished = 0;
}
/*****************************************************************************/
static void add_int(double newest,double *integral,int *numint,int maxint,
             struct i_data *data,short *finished,short *error)
/*****************************************************************************/
{
   *integral = *integral+newest;
   if (*numint == maxint) {
      *finished = 1;
      *error = 3;
   }  else  {
      if (stack == (struct s_item *)NULL)
         *finished = 1;
      else
         pop(data);
   }
}
/*****************************************************************************/
static void div_int(double middle,double v_middle,double newest1,double newest2,
             struct i_data *data)
/*****************************************************************************/
{
   struct i_data s;

   s.lb = data->lb;
   s.vlb = data->vlb;
   s.ub = middle;
   s.vub = v_middle;
   s.est = newest1;
   push(&s);
   data->lb = middle;
   data->vlb = v_middle;
   data->est = newest2;
}
/*****************************************************************************/
void gauss(double lo_lim,double up_lim,double tolerance,int maxint,
           double *integral,int *numint,short *error,double (*fun)(double))
/*****************************************************************************/
{
   double spacing,estimate,newest1,newest2,newest,middle,v_middle,ga_qua(),
          x1,x2;
   struct i_data data;
   short finished;

   te_ini(lo_lim,up_lim,tolerance,maxint,&data,integral,numint,
          &finished,error,fun);
   stack = (struct s_item *)NULL;
   if (!*error) {
      do {
         *numint = *numint+1;
         spacing = (data.ub-data.lb)/2.0;
         middle = data.lb+spacing;
         v_middle = (*fun)(middle);
         estimate = data.est;
         newest1 = ga_qua(data.lb,middle,fun);
         newest2 = ga_qua(middle,data.ub,fun);
         newest = newest1+newest2;
         x1 = fabs(estimate-newest);
         x2 = fabs(tolerance*newest);
         if (x1<=x2 || *numint==maxint) {
            add_int(newest,integral,numint,maxint,&data,&finished,error);
         }  else  {
            div_int(middle,v_middle,newest1,newest2,&data);
         }
      } while (!finished);
   }
}
