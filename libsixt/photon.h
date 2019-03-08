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

#ifndef PHOTON_H
#define PHOTON_H 1

#include "sixt.h"

#include "rmf.h"
#include "vector.h"


////////////////////////////////////////////////////////////////////////
// Definitions.
////////////////////////////////////////////////////////////////////////


// Counter for the number of entirely generated photons.
// extern long photon_counter;


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Contains all information about a single photon in the sky. */
typedef struct {
  /** Real time, when the photon is falling on the detector (in
      [s]). */
  double time;

  /** Photon energy in [keV]. */
  float energy;

  /** Right ascension and declination of photon position [rad]. */
  double ra, dec;

  /** Unique photon identifier. */
  long ph_id;

  /** Unique source identifier for the originating X-ray source. */
  long src_id;

} Photon;


/** Structure containing a photon and a pointer to the next photon in the
    time-ordered photon list. */
struct PhotonOrderedListEntry {
  Photon photon;
  struct PhotonOrderedListEntry *next;  // pointer to the next entry
};


/** Entry in the binary tree that stores the generated photons. */
struct PhotonBinaryTreeEntry {
  /** Photon data. */
  Photon photon;

  /** Pointer to entry with smaller time value. */
  struct PhotonBinaryTreeEntry* sptr;
  /** Pointer to entry with greater time value. */
  struct PhotonBinaryTreeEntry* gptr;
};


//////////////////////////////////////////////////////////////////////////
//   Function declarations.
//////////////////////////////////////////////////////////////////////////


/** Constructor. */
Photon* newPhoton(int* const status);

/** Copy Photon data structure. */
void copyPhoton(Photon* const dest, const Photon* const source);


#endif  /* PHOTON_H */
