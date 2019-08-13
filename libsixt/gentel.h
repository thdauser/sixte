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
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#ifndef GENTEL_H
#define GENTEL_H 1

#include "sixt.h"
#include "arf.h"
#include "phabkg.h"
#include "psf.h"
#include "vignetting.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic X-ray telescope. The characteristic properties for a
    particular telescope are defined in a specific XML file. */
typedef struct {

  /** Detector and telescope ARF containing the effective area. */
  char* arf_filename;
  struct ARF* arf;

  /** Telescope PSF. */
  PSF* psf;

  /** Telescope vignetting function. */
  Vignetting* vignetting;

  /** Focal length of the X-ray telescope [m]. */
  float focal_length;

  /** Diameter of the FoV [rad]. In the XML file the diameter is given
      in [deg], but it is converted to [rad] when parsing the XML
      file. */
  float fov_diameter;

  /** Number of Photons imaged by this telescope */
  long num_imaged;

} GenTel;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenTel data structure and
    initializes it with the values from the specified XML definition
    file. */
GenTel* newGenTel(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenTel data structure to NULL. */
void destroyGenTel(GenTel** const tel);

/** Check if any photons have been imaged by this telescope
    and print a warning if there weren't*/
void check_if_imaged(const GenTel* const tel);



#endif /* GENTEL_H */
