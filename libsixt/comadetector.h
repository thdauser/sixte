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

/** Add a new photon impact to the CoMaDetector pixels.  The
    generated charge is determined according to the detector response.
    If the charge cloud size is set, split events are calculated
    according to a Gaussian charge cloud model.  The new charge is
    added to the charge already contained in the detector pixel, so
    pileup effects are taken into account. */
int addImpact2CoMaDetector(CoMaDetector* det, Impact* impact);


#endif /* COMADETECTOR_H */

