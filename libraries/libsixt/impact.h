#ifndef IMPACT_H
#define IMPACT_H 1

#include "sixt.h"
#include "point.h"


/** Impact of a photon on the detector plane. */
typedef struct {
  /** Arrival time of the photon on the detector [s]. */
  double time;
  
  /** Photon energy [keV]. */
  float energy;
  
  /** Impact position of the photon on the detector [m]. */
  struct Point2d position;

  /** Unique photon identifier. */
  long ph_id;

  /** Unique source identifier for the originating X-ray source. */
  long src_id;
} Impact;


#endif /* IMPACT_H */
