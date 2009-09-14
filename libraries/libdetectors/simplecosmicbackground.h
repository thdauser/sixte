#ifndef SIMPLECOSMICBACKGROUND_H
#define SIMPLECOSMICBACKGROUND_H 1


#include "sixt.h"
#include "spectrum.h"

typedef struct {
  /** Background spectrum. */
  Spectrum *spectrum;

  /** Rate of detector background events. */
  float rate; 
} SimpleCosmicBackground;


#endif /* SIMPLECOSMICBACKGROUND */
