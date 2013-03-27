#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H 1


#include "codedmask.h"
#include "sixt.h"
/////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////
typedef struct {
  double** Rmap; //the reconstruction-array data built from mask-array-data.
  double**RImage; //the reconstructed image which is twice as large as the sky-image
  double open_fraction; //Transparency of the mask.
                        //Ratio #(transparent pixels)/#(all pixels)


  int naxis1, naxis2;    // Width of the image [pixel]
}ReconArray;

typedef struct {
  double** Rmap; //the reconstruction-array data built from mask-array-data.
  double open_fraction; //Transparency of the mask.
                        //Ratio #(transparent pixels)/#(all pixels)
  int naxis1, naxis2;    // Width of the image [pixel]
  double cdelt1, cdelt2; // Width of one pixel [m]
  double crpix1, crpix2; // Reference pixel [pixel]
  double crval1, crval2; // Value at reference pixel [m]
} ReconArrayFits;


/////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////

ReconArray* newReconArray(int* const status);

ReconArrayFits* newReconArrayForFits(int* const status);

ReconArray* getReconArray(const CodedMask* const mask, int* const status);

ReconArrayFits* getReconArrayForFits(const CodedMask* const mask, int* const status);

void FreeReconArray(ReconArray* recon);

void FreeReconArrayFits(ReconArrayFits* recon);

void FreeReconImage(ReconArray* recon);

void FreeReconArray1d(double* ReconArray1d);

double* SaveReconArray1d(ReconArray* recon, int* status);

void SaveReconArrayFits(ReconArrayFits* recon, char* filename, int* status);

#endif
