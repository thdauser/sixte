#ifndef REPIX_H
#define REPIX_H 1

#include "sixt.h"
#include "eventarray.h"
#include "testimg.h"


/////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////

//re-pixels from bigger pixels to smaller ones that fit without reminder into former bigger ones
//only for square pixels!

void repixNoReminder(void* arg_dataBeforeRepix, void* arg_dataAfterRepix, int type, 
		     int Size1, int Size2, double pixelwidth_big, double pixelwidth_small);

  //arg_dataBeforeRepix: pointer to array which shall be re-pixeled
  //arg_dataAfterRepix: pointer to array where re-pixeled version is stored
  //type: data type of both arrays as declared in 'testimg.h'
  //Size1/2: sizes of array which shall be re-pixeled (in pixels)
  //pixelwidth_big: size of pixels of given array before re-pixelization (in meters)
  //pixelwidth_small: size of pixels to which the array shall be re-pixeled to (in meters)


#endif /* REPIX_H */
