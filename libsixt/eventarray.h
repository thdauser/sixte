#ifndef EVENTARRAY_H
#define EVENTARRAY_H 1

#include "sixt.h"
#include "squarepixels.h"
#include "comaeventfile.h"
#include "codedmask.h"

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////

typedef struct {
  double** EventArray;
  int rawx, rawy;
  double charge;

  int naxis1, naxis2;
}ReadEvent;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////

ReadEvent* newEventArray(int* const status);

ReadEvent* getEventArray(CodedMask* mask, double Size1, double Size2, int* status);

ReadEvent* getEventArrayReBin(const CodedMask* const mask, SquarePixels* detector_pixels, int* status);

int readEventList_nextRow(CoMaEventFile* ef, ReadEvent* ea);

double* SaveEventArray1d(ReadEvent* ea, int* status);

void FreeEventArray1d(double* EventArray1d);

void FreeEventArray(ReadEvent* ea);

#endif
