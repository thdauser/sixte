/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

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

ReadEvent* getEventArray(int Size1, int Size2, int* status);

ReadEvent* getEventArrayReBin(const CodedMask* const mask, SquarePixels* detector_pixels, int* status);

int readEventList_nextRow(CoMaEventFile* ef, ReadEvent* ea);

double* SaveEventArray1d(ReadEvent* ea, int* status);

void FreeEventArray1d(double* EventArray1d);

void FreeEventArray(ReadEvent* ea);

#endif
