#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H 1


#include "codedmask.h"
#include "sixt.h"
#include "squarepixels.h"
#include "repix.h"

/////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////
typedef struct {
  double** Rmap; //the reconstruction-array data built from mask-array-data.
  double open_fraction; //Transparency of the mask.
                        //Ratio #(transparent pixels)/#(all pixels)
  int naxis1, naxis2;   // Width of the image [pixel]
}ReconArray;


/////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////
ReconArray* newReconArray(int* const status);

ReconArray* getReconArray(CodedMask* mask, SquarePixels* detector_pixels, int* const status);

double* SaveReconArray1d(ReconArray* recon, int* status);

void FreeReconArray1d(double* ReconArray1d);

void FreeReconArray(ReconArray** const recon);

//void FreeReconImage(ReconArray* recon);


#endif
