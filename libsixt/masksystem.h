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


   Copyright 2007 - 2018: Christian Schmid, Mirjam Oertel, FAU.
   Manuel Castro, National Institute for Space Research (INPE),
		 Brazil; under grant #2017/00968-6,
		 São Paulo Research Foundation (FAPESP).

   Copyright 2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                  Erlangen-Nuernberg
*/

#ifndef MASKSYSTEM_H
#define MASKSYSTEM_H 1

#include "sixt.h"
#include "xmlbuffer.h"


typedef struct {

//  /** Telescope data structure. */
//  GenTel* tel;
//
//  /** Detector data structure. */
//  GenDet* det;
//
  /** File name (without path contributions) of the FITS file
      containing the XML detector definition. */
  char* filename;

  /** Path to the FITS file containing the XML detector definition. */
  char* filepath;
//
//  /** FITS header keyword TELESCOP. */
//  char* telescop;
//
//  /** FITS header keyword INSTRUME. */
//  char* instrume;

//------Mask parameters-----

//Mask dimension (in m)
float x_mask;
float y_mask;

//Mask-detector distance (in m)
float mask_distance;

//Mask pattern - FITS file
char* mask_pattern;

//Collimator height (in m)
float collimator_height;

//Wall thickness (for protomirax,in m )
float wall_thickness;

//Flag (protomirax=1, mirax=0)
int flag;

//-----Detector parameters-----

//Detector dimensions ( in m)
float x_det;
float y_det;

//detector pixel width (in m)
float det_pixelwidth;

//
float det_width;

//Number of pixels
int width;

//
float DCU_length;

//
float DCU_gap;

//
float DCA_gap;

//
float pixelwidth;

//
float repixsize;

} MaskSystem;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new MaskSystem data structure and
    initializes it with the values from the specified XML definition
    file. The second parameter determines, whether the cosmic ray
    detector background model should be activated. */
MaskSystem* newMaskSystem(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the MaskSystem data structure to NULL. */
void destroyMaskSystem(MaskSystem** const det, int* const status);

/** Parse the MaskSystem definition from an XML file. */
MaskSystem* loadMaskSystem(const char* const filename,
		     const unsigned int seed,
		     int* const status);


#endif /* MASKSYSTEM_H */
