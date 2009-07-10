#ifndef FRAMESTORE_H
#define FRAMESTORE_H 1

#include "sixt.h"
#include "impactlist.h"
#include "genericdetector.h"
#include "squarepixels.h"


typedef struct {

  GenericDetector generic;
  SquarePixels pixels;

  /** Integration time of the entire pnCCD (!) detector array.
   * (= Span of time between 2 subsequent readouts). */
  double integration_time; 

  long frame; /**< Number of the current frame. */

  double readout_time; /**< Current readout time. The end of the integration 
			* time / beginning of dead time. */

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

void readoutFramestoreDetector(FramestoreDetector*);
void addImpact2FramestoreDetector(FramestoreDetector*);



#endif /* FRAMESTORE_H */

