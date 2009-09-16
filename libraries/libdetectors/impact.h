#ifndef IMPACT_H
#define IMPACT_H 1

#include "sixt.h"
#include "point.h"


/** Impact of a photon on the detector plane. */
typedef struct {
  double time;
  float energy;
  struct Point2d position;
} Impact;


#endif /* IMPACT_H */
