#ifndef BALANCING_H
#define BALANCING_H 1


#include "sixt.h"
#include "squarepixels.h"
#include "reconstruction.h"
#include "eventfile.h"


/////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////

typedef struct {
  double** Bmap; //the balancing-array data built from reconstruction-array-data.
  int naxis1, naxis2;    // Width of the image [pixel]
}BalancingArray;



/////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////

BalancingArray* newBalancingArray(int* const status);

BalancingArray* getBalancingArray(ReconArray* recon, SquarePixels* detector_pixels, EventFile* ef, int* const status);

#endif
