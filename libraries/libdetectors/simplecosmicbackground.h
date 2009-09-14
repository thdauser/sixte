#ifndef SIMPLECOSMICBACKGROUND_H
#define SIMPLECOSMICBACKGROUND_H 1


#include "sixt.h"
#include "spectrum.h"


/** Simple model for Cosmic X-ray detector background. */
typedef struct {
  /** Background spectrum. */
  Spectrum spectrum;

  /** Rate of detector background events. */
  float rate; 
} SimpleCosmicBackground;


struct SimpleCosmicBackgroundParameters {
  char* spectrum_filename;
  float rate;
};


////////////////////////////////////////////////////////


/** Set up the configuration of the SimpleCosmicBackground. */
int initSimpleCosmicBackground(SimpleCosmicBackground*, 
			       struct SimpleCosmicBackgroundParameters*);

/** Clean up the SimpleCosmicBackground data structure. 
 * Release allocated memory and call clean-up routines of underlying structures. */
int cleanupSimpleCosmicBackground(SimpleCosmicBackground*);


#endif /* SIMPLECOSMICBACKGROUND_H */
