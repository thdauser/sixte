#ifndef PHOTON_H
#define PHOTON_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#include "extendedsources.h"
#include "pointsources.h"
#include "rmf.h"
#include "vector.h"


////////////////////////////////////////////////////////////////////////
// Definitions.
////////////////////////////////////////////////////////////////////////


// Counter for the number of entirely generated photons.
long photon_counter;


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

