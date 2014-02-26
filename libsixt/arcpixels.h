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

#ifndef ARCPIXELS_H
#define ARCPIXELS_H 1

#include "sixt.h"
#include "point.h"
#include "genericdetector.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** One single ArcPixel. */
typedef struct {
  /** Charge stored in this detector pixel. The charge from several
      photon impacts is added up until the pixel has to be cleared,
      because its output voltage exceeds the input range of the
      subsequent electronics. Actually the variable contains not the
      electrical charge but the corresponding energy in [keV]. */
  float charge; 

  /** Time of the last photon impact in this pixel. This time is
      required in order to determine for the next photon impact,
      whether the shaping process of the previous pulse is still in
      progress. */
  double last_impact;

  /** Start time of the next reset time interval. If a reset of a
      pixel has to be performed after an event at the time t, this
      variable is set to the value t+shaping_time, such that the
      variable denotes the beginning of the reset interval. The pixel
      is insensitive to incident photons until the end of the reset
      period. Photons falling on the pixel during that time are marked
      with the corresponding event grade. */
  double reset_from;

} ArcPixel;


/** Array of ArcPixels. */
typedef struct {
  /** 2-dimensional pixel array with nrings rows of different lengths
      (each has the length npixels[ring]) .*/
  ArcPixel** array;
  
  /** Number of pixel rings. */
  int nrings; 
  /** Number of pixels in each individual ring. The array has nrings
      elements. */
  int* npixels;

  /** Radii of the individual detector rings. The array has nrings
      elements each specifying the outer radius of the respective
      pixel ring. The inner radius of the ring i is given by the
      radius of the ring i-1 (for i>0). */
  double* radius;
  /** Offset angle of the first pixel in each ring ([rad]). In
      general the pixels may not start the 3 o'clock position in the
      rings, but the border of pixel with the index 0 may have an
      angular offset from the 3 o'clock position. These offsets are
      measured in counter-clock-wise direction. */
  double* offset_angle;

  /** Width of the spokes of the mask on the HTRS detector ([m]). If
      no mask is used, the width has to be set to 0. */
  double mask_spoke_width;

} ArcPixels;


struct ArcPixelsParameters {
  int nrings;
  int* npixels;
  double* radius;
  double* offset_angle;
  double mask_spoke_width;
};


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Initialization routine for the ArcPixels data structure. Sets the
    basic properties and allocates memory for the pixel array. The
    return value is the error status. */
int initArcPixels(ArcPixels* ap, struct ArcPixelsParameters* app);

/** Clean up the ArcPixels data structure. E.g. release allocated
    memory. */
void cleanupArcPixels(ArcPixels* ap);

/** Determine the ArcPixel affected by a photon impact. The function
    return value is a number, either 1 or 2, giving the number of
    affected pixels. The function returns the affected pixel rings and
    the pixel numbers within these rings (parameters 'ring' and
    'number', which are both arrays of length 1 or 2 for single and
    double events respectively). If the HTRS detector is operated with
    a mask on top, this has to be given as a parameter. Otherwise the
    parameter can also be set to NULL. */
int getArcPixelSplits(ArcPixels* ap, GenericDetector* gd,
		      struct Point2d position, 
		      int* ring, int* number);

/** Determine the ArcPixel that contains the specified position given
    in polar coordinates. The function returns the affected pixel ring
    and the pixel number within this ring (parameters 'ring' and
    'number'). */
void getArcPixelFromPolar(ArcPixels* ap,
			  double radius, double angle,
			  int* ring, int* number);

/** Determine the absolute pixel index from a given ring and the pixel
    number within this ring. Currently there are two numbering schemes
    used for the HTRS: The first is 2-dimensional: each pixel is
    specified by its ring and its number within this ring. The second
    is linear: each pixel has a unique number for the whole
    detector. This function reads the 2-dimensional coordinates in the
    former system and returns the corresponding unique pixel number in
    the latter linear system. */
int getArcPixelIndex(ArcPixels* ap, int ring, int number);

/** Determine the polar coordinates 'radius' and 'angle' for a given
    position in 2-dimensional carteesian coordinates. The return
    values are within the following margins: radius [0; infinity),
    angle[0:2pi). */
void getPolarCoordinates(struct Point2d position, 
			 double* radius, double* angle);

/** Check if the specified impact position lies on the mask of the
    HTRS detector. If the photon is absorbed by the mask, the return
    value is 1. Otherwise the function returns 0. */
int HTRSisPositionOnMask(ArcPixels* ap, int pixel, int number,
			 double radius, double angle);


#endif /* ARCPIXELS_H */

