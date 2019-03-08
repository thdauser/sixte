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

#ifndef SKYIMAGE_H
#define SKYIMAGE_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////
typedef struct{
  //all the 2-d arrays below are addressed via the same integer-indizes
  double** pixel;   //actual value in particular pixel
  double** ra;      //corresponding RA-value
  double** dec;     //corresponding dec-value

  int naxis1, naxis2;   // Width of the image [pixel]
  float cdelt1, cdelt2; // Width of one pixel [deg]
  float crpix1, crpix2; // [pixel]
  float crval1, crval2; // [deg]

  float minra, maxra;   // Maximum right ascension covered by the image [deg]
  float mindec, maxdec; // Maximum declination covered by the image [deg]
} SkyImage;

struct SkyImageParameters {
  /** Dimensions of the sky image [pixels]. */
  int naxis1, naxis2;
  /** Dimensions of the sky image [deg]. */
  float minra, maxra;
  float mindec, maxdec;
  /** Width of the image pixels [deg]. */
  float cdelt1, cdelt2;
  /** WCS parameters. */
  float crpix1, crpix2;
  float crval1, crval2;
};

/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////

/** Constructor for SkyImage objects returning the bare data
    structure without any memory allocated for the pixels. */
SkyImage* get_SkyImage();

/** Constructor for SkyImage objects generating an empty sky
    image with allocated memory for the requested image dimensions. */
SkyImage* getEmptySkyImage(struct SkyImageParameters* sip, int* status);

void fillRaDecArrays(SkyImage* si);

/** Destructor for SkyImage objects. */
void free_SkyImage(SkyImage* si);

void saveSkyImage(SkyImage* si, char* filename, int* status);

#endif /* SKYIMAGE_H */
