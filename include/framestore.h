#ifndef FRAMESTORE_H
#define FRAMESTORE_H 1

#include "sixt.h"
#include "impactlist.h"
#include "eventlist.h"
#include "genericdetector.h"
#include "squarepixels.h"


typedef struct {

  GenericDetector generic;
  SquarePixels pixels;

  /** Integration time of the entire pnCCD (!) detector array.
   * (= Span of time between 2 subsequent readouts). */
  double integration_time; 

  double readout_time; /**< Current readout time. The end of the integration 
			* time / beginning of dead time. */

  long frame; /**< Number of the current frame. */

} FramestoreDetector;


struct FramestoreDetectorParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;
  double integration_time;
  double t0;
};



////////////////////////////////////////////////////////

/** Set up the configuration of a FramestoreDetector. 
 * The routine also calls the init routines of the underlying data structures. */
int initFramestoreDetector(FramestoreDetector*, struct FramestoreDetectorParameters*);

/** Clean up the FramestoreDetector data structure. 
 * Release allocated memory and call clean-up routines of underlying structures. */
void cleanupFramestoreDetector(FramestoreDetector* fd);

/** This routine is called for readout of the FramestoreDetector.
 * The routine itself checks, whether a readout is necessary according to the current
 * time and readout time. If it's time to do a readout, the routine calls the function
 * readoutFramestoreDetector(). */
int checkReadoutFramestoreDetector(FramestoreDetector*, double time, struct Eventlist_File*);

/** Read out the FramestoreDetector.
 * The measured events are stored in an event list FITS file. */
inline int readoutFramestoreDetector(FramestoreDetector*, struct Eventlist_File*);

/** Add a new photon impact to the FramestoreDetector pixels.
 * The generated charge is determined according to the detector response.
 * If the charge cloud size is set, split events are calculated according to
 * a Gaussian charge cloud model.
 * The new charge is added to the charge already contained in the detector pixel, 
 * so pileup effects are taken into account. */
void addImpact2FramestoreDetector(FramestoreDetector*, Impact*);

/** Returns 1 if the detector is sensitive in the pixel (x,y) at the specified time. */
//inline int FramestoreDetectorIsSensitive(int x, int y, FramestoreDetector* fd, double time);


#endif /* FRAMESTORE_H */

