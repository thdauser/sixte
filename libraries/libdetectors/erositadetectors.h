#ifndef EROSITADETECTORS_H
#define EROSITADETECTORS_H 1


#include "sixt.h"
#include "impact.h"
#include "eventfile.h"
#include "erositaeventfile.h"
#include "erositaevent.h"
#include "genericdetector.h"
#include "squarepixels.h"


/** Number of eROSITA telescopes / framestore pnCCDs. */
#define NeROSITATELESCOPES 7


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {
  GenericDetector generic;
  SquarePixels pixels[NeROSITATELESCOPES];

  /** Integration time of the entire pnCCD (!) detector array. (=
      Span of time between 2 subsequent readouts). */
  double integration_time; 

  /** Current readout time. The end of the integration time /
      beginning of dead time. */
  double readout_time;

  /** Number of the current frame. */
  long frame; 

  /** Event list FITS file for the eROSITA-specific events. */
  eROSITAEventFile eventlist;

} eROSITADetectors;


struct eROSITADetectorsParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double integration_time;
  double t0;
  char* eventlist_filename;
  char* eventlist_template;  
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Set up the configuration of a eROSITADetectors. The routine also
    calls the init routines of the underlying data structures. */
int initeROSITADetectors(eROSITADetectors*, struct eROSITADetectorsParameters*);

/** Clean up the eROSITADetectors data structure. Release allocated
    memory and call clean-up routines of underlying structures. */
int cleanupeROSITADetectors(eROSITADetectors*);

/** This routine is called for readout of the eROSITADetectors. The
    routine itself checks, whether a readout is necessary according to
    the current time and readout time. If it's time to do a readout,
    the routine initiates the readout process by calling
    readouteROSITADetectors(). */
int checkReadouteROSITADetectors(eROSITADetectors*, double time);

/** Read out the eROSITADetectors. The measured events are stored in
    an event list FITS file. */
inline int readouteROSITADetectors(eROSITADetectors*);

/** Add a new photon impact to the eROSITADetectors pixels. The
    generated charge is determined according to the detector response.
    If the charge cloud size is set, split events are calculated
    according to a Gaussian charge cloud model. The new charge is
    added to the charge already contained in the detector pixel, so
    pileup effects are taken into account. */
int addImpact2eROSITADetectors(eROSITADetectors*, Impact*);


#endif /* EROSITADETECTORS_H */

