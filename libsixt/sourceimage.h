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


   Copyright 2007-2014 Christian Schmid, FAU
*/

#ifndef SOURCEIMAGE_H
#define SOURCEIMAGE_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Object containing the extended source image from one FITS image extension. */
typedef struct {
  /** Pixel probability distribution. */
  double** pixel;

  /** If this flag is set (1), the SourceImage pixels contain a
      normalized probability distribution. If not (0), the pixels
      contain a real image. */
  int accumulated;

  /** Total photon rate of the entire SourceImage. This is the sum of
      the pixels read from the FITS file. In the computer memory the
      image is normalized to 1. */
  double total_rate;

  /** Time of the last photon emitted from this SourceImage. */
  double t_last_photon;

  int naxis1, naxis2;   /**< Width of the image [pixel]. */
  float cdelt1, cdelt2; /**< Width of one pixel [rad]. */
  float crpix1, crpix2; /**< [pixel] */
  float crval1, crval2; /**< [rad] */

  float minra, maxra;   /**< Maximum right ascension covered by the image [rad]. */
  float mindec, maxdec; /**< Maximum declination covered by the image [rad]. */

} SourceImage;


/** Catalog of all extended SourceImage Objects. */
typedef struct {
  /** Total number of ClusterImage elements in the catalog. */
  int nimages;
  /** Catalog of ClusterImage elements. */
  SourceImage** images; /* nimages */

} SourceImageCatalog;


struct SourceImageParameters {
  /** Dimensions of the source image. */
  int naxis1, naxis2;
  /** Width of the image pixels [rad]. */
  float cdelt1, cdelt2;
  /** WCS parameters. */
  float crpix1, crpix2;
  float crval1, crval2;
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for SourceImage objects returning the bare data
    structure without any memory allocated for the pixels. */
SourceImage* get_SourceImage();

/** Constructor for SourceImage objects generating an empty source
    image with allocated memory for the requested image dimensions. */
SourceImage* getEmptySourceImage(struct SourceImageParameters* sip, int* status);

/** Constructor for SourceImage objects. Reads a SourceImage from a
    FITS file and stores it in the corresponding data structure. */
SourceImage* get_SourceImage_fromFile(char* filename, int* status);

/** Constructor for SourceImage objects. Reads a SourceImage from a
    FITS file and stores it in the corresponding data structure. */
SourceImage* get_SourceImage_fromHDU(fitsfile* fptr, int* status);

/** Destructor for SourceImage objects. */
void free_SourceImage(SourceImage* si);

/** Store the SourceImage in a FITS file. */
void saveSourceImage(SourceImage* si, char* filename, int* status);

/** Constructor for the SourceImageCatalog. */
SourceImageCatalog* get_SourceImageCatalog();

/** Destructor for the SourceImageCatalog. */
void free_SourceImageCatalog(SourceImageCatalog* sic);

/** Add new SourceImage object to the SourceImageCatalog. The new
    image is loaded from the FITS HDU designated by fptr and added to
    the SourceImageCatalog. The return value is the error status. */
int addSourceImage2Catalog(SourceImageCatalog* sic, fitsfile* fptr);

/** Determine a random pixel according to the probability distribution
    given by the SourceImage. */
void getRandomSourceImagePixel(SourceImage* si, int* x, int* y,
			       int* const status);


#endif /* SOURCEIMAGE_H */

