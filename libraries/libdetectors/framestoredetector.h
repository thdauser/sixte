#ifndef FRAMESTOREDETECTOR_H
#define FRAMESTOREDETECTOR_H 1


#include "sixt.h"
#include "impact.h"
#include "eventfile.h"
#include "erositaeventfile.h"
#include "erositaevent.h"
#include "genericdetector.h"
#include "squarepixels.h"


////////////////////////////////////////////////////////////////////////
// Type declarations.
////////////////////////////////////////////////////////////////////////


typedef struct {
  GenericDetector generic;
  SquarePixels pixels;

  /** Integration time of the entire pnCCD (!) detector array.
   * (= Span of time between 2 subsequent readouts). */
  double integration_time; 

  /** Current readout time. The end of the integration 
   * time / beginning of dead time. */
  double readout_time; 

  /** Number of the current frame. */
  long frame; 

  /** Event list FITS file for the eROSITA-specific events. */
  eROSITAEventFile eventlist; 

} FramestoreDetector;


struct FramestoreDetectorParameters {
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


/** Set up the configuration of a FramestoreDetector. 
 * The routine also calls the init routines of the underlying data structures. */
int initFramestoreDetector(FramestoreDetector*, struct FramestoreDetectorParameters*);

/** Clean up the FramestoreDetector data structure. 
 * Release allocated memory and call clean-up routines of underlying structures. */
int cleanupFramestoreDetector(FramestoreDetector*);

/** This routine is called for readout of the FramestoreDetector.
 * The routine itself checks, whether a readout is necessary according to the current
 * time and readout time. If it's time to do a readout, the routine initiates the 
 * readout process by calling readoutFramestoreDetector(). */
int checkReadoutFramestoreDetector(FramestoreDetector*, double time);

/** Read out the FramestoreDetector.
 * The measured events are stored in an event list FITS file. */
inline int readoutFramestoreDetector(FramestoreDetector*);

/** Add a new photon impact to the FramestoreDetector pixels.
 * The generated charge is determined according to the detector response.
 * If the charge cloud size is set, split events are calculated according to
 * a Gaussian charge cloud model.
 * The new charge is added to the charge already contained in the detector pixel, 
 * so pileup effects are taken into account. */
int addImpact2FramestoreDetector(FramestoreDetector*, Impact*);


#endif /* FRAMESTOREDETECTOR_H */

