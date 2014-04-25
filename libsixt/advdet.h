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


   Copyright 2014 Thorsten Brand, FAU
*/

#ifndef ADVDET_H
#define ADVDET_H 1

#include "sixt.h"
#include "impact.h"
#include "point.h"
#include "xmlbuffer.h"

////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////

/** Data structure describing a pixel with arbitrary geometry.
    Current implementation: only rectangulars, parallel to detector
    coordinate system. */
typedef struct{
  
  /** x-shift of pixel in respect to the detector reference point */
  double sx;
  
  /** y-shift of pixel in respect to the detector reference point */
  double sy;  

  /** Width of the pixel [m], centered around pixel reference point. */
  double width;

  /** Height of the pixel [m], centered around pixel reference point. */
  double height;
  
  /** Index of pixel in detector structure. */
  int pindex;

}AdvPix;


/** Data structure describing the geometry of a pixel detector with
    arbitrary pixel geometry. */
typedef struct{
  
  /** x-shift of detector in respect to the focal point */
  double sx;
  
  /** y-shift of detector in respect to the focal point */
  double sy;
 
  /** Number of pixels. */
  int npix;
  
  /** Counter for operations on pixels */
  int cpix;

  /** array of pixels. */
  AdvPix *pix;

  /** File name (without path contributions) of the FITS file
      containing the XML detector definition. */
  char* filename;

  /** Path to the FITS file containing the XML detector definition. */
  char* filepath;

}AdvDet;

/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////

/** Read the advanced detector syntax from the specified XML */
void parseAdvDetXML(AdvDet* const det, 
	       const char* const filename,
	       int* const status);

/** Constructor. Allocates memory for a new AdvDet data structure. */
AdvDet* newAdvDet(int* const status);

/** Loads AdvDet information from XML file */
AdvDet* loadAdvDet(const char* const filename,
		     int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the AdvDet data structure to NULL. */
void destroyAdvDet(AdvDet **det, int* const status);

/** Function testing if an impact lies inside a pixel. */
int CheckAdvPixImpact(AdvPix pix, Impact *imp);

/** Function calculating the exact impact position and time in
    pixel coordinates. */
void CalcAdvPixImpact(AdvPix pix, Impact *imp, Impact *piximp);

/** Function determining the pixel indices which have an impact from one 
    event. Gives the number of pixels that were hit.*/
int AdvImpactList(AdvDet *det, Impact *imp, long **pixindex, Impact **piximp);


#endif /* ADVDET_H */