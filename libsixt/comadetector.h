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


   Copyright 2007-2014 Christian Schmid, Mirjam Oertel, FAU
*/

#ifndef COMADETECTOR_H
#define COMADETECTOR_H 1


#include "sixt.h"
#include "impact.h"
#include "eventlist.h"
#include "squarepixels.h"
#include "comaevent.h"
#include "comaeventfile.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {
  SquarePixels* pixels;

  /** Event list FITS file for the eROSITA-specific events. */
  CoMaEventFile* eventfile; 

} CoMaDetector;


struct CoMaDetectorParameters {
  struct SquarePixelsParameters pixels;

  char* eventfile_filename;
  char* eventfile_template;  
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for CoMaDetector. Get memory and set up the
    configuration of a CoMaDetector object. The routine also calls the
    init routines of the underlying data structures. */
CoMaDetector* getCoMaDetector(struct CoMaDetectorParameters* parameters, 
			      int* status);

/** Destructor of the CoMaDetector data structure. Release allocated
    memory and call clean-up routines of underlying structures. */
void freeCoMaDetector(CoMaDetector* det);

/** Add a new photon impact to the CoMaDetector pixels. The
    generated charge is determined according to the detector response.
    If the charge cloud size is set, split events are calculated
    according to a Gaussian charge cloud model.  The new charge is
    added to the charge already contained in the detector pixel, so
    pileup effects are taken into account. */
int addImpact2CoMaDetector(CoMaDetector* det, Impact* impact);

int addImpact2CoMaDetector_protoMirax(CoMaDetector* det, Impact* impact);


#endif /* COMADETECTOR_H */

