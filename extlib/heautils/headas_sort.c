/***********************************************************
 Routine name: HDsort
 Description:
          Do quick sort on the input array. 
          unlike C qsort which returns sorted data, it returns 
          sorted index instead of data.
 Parameters:
          array - (IN) unsorted data array
          index - (IN/OUT) array index
          n     - (IN) size of array
 Modification History:
          writen by Ziqin Pan, Novemeber 2004
****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "headas_polyfit.h"

typedef struct _data {
   float base;
   int index;
} Data; 

void HDsort(float * array, int * index, int n) {

   Data * data;
   int i;

   data = (Data *) calloc(n, sizeof(Data));

   for (i=0; i<n; i++) {
	data[i].base =array[i];
        data[i].index = index[i];
   }

   qsort(data,n,sizeof(Data),HDcmp);
  

   for (i=0; i<n; i++) {
        index[i] = data[i].index;
   }
   free (data);
}
   



int HDcmp(const void * data1, const void * data2) {

   if ( ((Data *) data1)->base < ((Data *) data2)->base ) {
         return -1;
   } else if ( ((Data *) data1)->base  == ((Data *) data2)->base ) {
         return 0;
   } else {
         return 1;
   }

}
