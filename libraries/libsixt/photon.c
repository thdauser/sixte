#include "photon.h"


Photon* newPhoton(int* const status)
{
  Photon* ph = (Photon*)malloc(sizeof(Photon));
  CHECK_NULL(ph, *status, "memory allocation for Photon failed");

  // Set initial values.
  ph->ra     = 0.;
  ph->dec    = 0.;
  ph->time   = 0.;
  ph->energy = 0.;
  ph->ph_id  = 0;
  ph->src_id = 0;

  return(ph);
}


void copyPhoton(Photon* const dest, const Photon* const source)
{
  dest->time   = source->time;
  dest->ra     = source->ra;
  dest->dec    = source->dec;
  dest->energy = source->energy;
  dest->ph_id  = source->ph_id;
  dest->src_id = source->src_id;
}

