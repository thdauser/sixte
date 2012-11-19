#ifndef RMF_H
#define RMF_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Load an RMF/RSP matrix and the corresponding EBOUNDS from a
    response file. If the flag 'NORMALIZE_RMF' is set, the RSP is
    renormalized to an RMF in such a way that the sum of each matrix
    row is 1. */
struct RMF* loadRMF(char* filename, int* const status);

/** Destructor for the RMF data structure. Warning: As there is no
    internal destructor for the RMF data structure in the HEASP
    library, the memory allocated by the function ReadRMFMatrix is not
    realeased. */
void freeRMF(struct RMF* const rmf);

/** Returns a randomly selected channel for a photon of the given
    input energy. The channel number can start at 0, 1, or any other
    positive value defined in the RMF. The return value '-1' means
    that no appropriate channel could be determined. */
void returnRMFChannel(struct RMF *rmf, 
		      const float energy, 
		      long* const channel);

/** Determines the PHA channel corresponding to a given energy
    according to the EBOUNDS table of the detector response. The
    routine performs a binary search to obtain the PHA channel the
    specified energy lies within. The energy has to be given in the
    same units as the EBOUNDS are, which usually is [keV]. Note that
    the routine is NOT doing an RMF randomization of the measured
    channel. In case the channel cannot be determined, because there
    is no RMF available or the desired energy lies outside the covered
    range, a negative value is returned. */
long getEBOUNDSChannel(const float energy, const struct RMF* const rmf);

/** Determine the charge corresponding to a particular PHA channel
    according to the EBOUNDS table. The input channel must have the
    same offset as in the EBOUNDS table. I.e. if the first channel in
    the EBOUNDS has the number 1, the numbering starts at 1. If the
    first channel has the number 0, the numbering starts at 0.  The
    returned energy is given in the same units as the EBOUNDS,
    (usually [keV]). The boundary flag determines, whether the lower
    (-1), randomized mean (0), or upper (else) boundary of the energy
    bin are returned. */
float getEBOUNDSEnergy(long channel,
		       const struct RMF* const rmf, 
		       const int boundary,
		       int* const status);


#endif /* RMF_H */
