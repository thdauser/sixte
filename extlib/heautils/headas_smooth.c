/************************************************************
 Name of routine: HDsmooth
 Description:
    Do a boxcar average on input data.
 Parameters:
    input  - the unsmoothed array
    output - the smoothed array
    num    - the size of array
    width  - the width of the boxcar
 Modification History: 
    Writen by: Ziqin Pan, Novemeber,2004

************************************************************/

#include <stdio.h>
#include <stdlib.h>
void HDsmooth(float * input, float * output, int num, int width) {
        int i,j,k,b;
        float total;
        float * tmp;

        b = width/2;
        tmp = (float *) calloc(num,sizeof(float));

        for (i=0; i< num; i++) {
                k=0;
                total=0.0;
		for ( j=i-b; j <=i+b; j++ ) {
                     if (j >=0 && j <num) {
                        k++;
                        total +=*(input+j);
                     }
                }
                tmp[i] = total/k;
         }
         for (i=0; i<num; i++) {
             output[i]=tmp[i];
         }
/*
         for (i=0; i<num; i++) {
             fprintf(stderr,"%d,tmp[%d]=%f\n",num,i,tmp[i]);
             output[i]=tmp[i];
             fprintf(stderr,"output[%d]=%f\n",i,output[i]);
         }
*/
         free(tmp);
}
