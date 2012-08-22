/***************************************************************
  Name of routine: HDpoly_fit
  Description: 
         Do a least-square polynomial fit:
            y = a(0) + a(1)*x +a(2)*x^2+...+a(degree)*x^degree
  Parameters:
         x - (I) an n-element array of independent variables
         y - (I) an n-element array of dependent variables
         a - (O) an degree+1 element of coefficient
         n - (I) size of array of x or y 
    degree - (I) the degree if the polynomial to fit
  Modification history: 
         written by Ziqin Pan , Novemeber ,2004
***************************************************************/
        
  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headas_polyfit.h"

void HDfuncs (double x, double * afunc, unsigned int m ) {
   int i;
   afunc[0] =1.;
   for ( i=1; i <m; i++) {
      afunc[i] = afunc[i-1]*x;
   }
}


void HDpoly_fit(double * x, double * y, double * a, int n, int degree) {

   int ma;
   double * sig;
   double ** u;
   double ** v;
   double * w;
   double * chisq;


   int i;

   ma = degree+1;
   sig =(double*) malloc(n*sizeof(double));
   w =(double*) malloc(n*sizeof(double));
   chisq =(double*) malloc(n*sizeof(double));
   u =(double**) malloc(n*sizeof(double*));
   v =(double**) malloc(n*sizeof(double*));

   for (i=0; i<n; i++) {
      sig[i]=1;
      u[i] = (double*) malloc(n*sizeof(double));
      v[i] = (double*) malloc(n*sizeof(double));
   }

   HDsvdfit (x,y,sig,n,a,ma,u,v,w,n,n,chisq,HDfuncs);

   for (i=0; i<n; i++) {
     if (u[i]) free(u[i]);
     if (v[i]) free(v[i]);
   }
   if (u) free(u);
   if (v) free(v);

}
/*
int main(int argc, char ** argv) {
      int n=11;
      double x[11];
      double y[11];
      double a[3];
      int i;

      for(i=0 ; i<n; i++) {
         x[i]= 0.1*i;
      }
 
      y[0]=0.25;
      y[1]=0.16;
      y[2]=0.09;
      y[3]=0.04;
      y[4]=0.01;
      y[5]=0.00;
      y[6]=0.01;
      y[7]=0.04;
      y[8]=0.09;
      y[9]=0.16;
      y[10]=0.25;

      fprintf(stderr,"hello\n");

      poly_fit(x,y,a,n,2);
 
      for (i=0; i<=2; i++) {
	fprintf(stderr,"a[%d]=%f\n",i,a[i]);
      }
}
*/



