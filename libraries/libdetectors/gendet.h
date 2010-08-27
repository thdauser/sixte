#ifndef GENDET_H 
#define GENDET_H 1

#include "sixt.h"
#include "gendetline.h"
#include "genericdetector.h"
#include "impact.h"
#include "clocklist.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/////////////////////////////////////////////////////////////////
// Constants
/////////////////////////////////////////////////////////////////


#define GENDET_TIME_TRIGGERED  1
#define GENDET_EVENT_TRIGGERED 2


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic pixelized X-ray detector. The data structure is designed
    in a generic way. The characteristic properties for a particular
    detector are defined in a detector specific XML file. */
typedef struct {

  /** Detector dimensions. Width and height [pixels]. */
  int xwidth, ywidth;

  /** Reference pixel. */
  float xrpix, yrpix;
  /** Reference value [m]. */
  float xrval, yrval;
  /** Pixel width [m]. */
  float xdelt, ydelt;

  /** Array of pointers to pixel lines. */
  GenDetLine** line;

  /** Detector response matrix. The RSP file that is originally loaded
      may also contain ARF contributions. But as they already have
      been taken into account in the generation of the input spectra
      for the X-ray sources, the ARF contributions have to be removed
      by normalizing the RSP matrix. */
  struct RMF* rmf;

  /** Flag for detector readout trigger. The readout can be triggered
      either by an incoming photon event (GENDET_EVENT_TRIGGERED) or
      by a timing clock (GENDET_TIME_TRIGGERED). */
  int readout_trigger;

  /** List of clock operations for time-triggered detectors. */
  ClockList* clocklist;

} GenDet;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenDet data structure. */
GenDet* newGenDet(const char* const filename, int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenDet data structure to NULL. */
void destroyGenDet(GenDet** det);

/** Add a new photon impact to the detector. */
void addGenDetPhotonImpact(GenDet* const det, const Impact* const impact, int* const status);


#endif /* GENDET_H */
