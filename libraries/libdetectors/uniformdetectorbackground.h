#ifndef UNIFORMDETECTORBACKGROUND_H
#define UNIFORMDETECTORBACKGROUND_H 1


#include "sixt.h"
#include "spectrum.h"


/** Simple model for Cosmic X-ray detector background. */
typedef struct {
  /** Background spectrum. */
  Spectrum spectrum;

  /** Rate of detector background events. */
  float rate; 
} UniformDetectorBackground;


struct UniformDetectorBackgroundParameters {
  char* spectrum_filename;
  float rate;
};


////////////////////////////////////////////////////////


/** Set up the configuration of the UniformDetectorBackground. */
int initUniformDetectorBackground(UniformDetectorBackground*, 
				  struct UniformDetectorBackgroundParameters*);

/** Clean up the UniformDetectorBackground data structure. 
 * Release allocated memory and call clean-up routines of underlying structures. */
int cleanupUniformDetectorBackground(UniformDetectorBackground*);


#endif /* UNIFORMDETECTORBACKGROUND_H */
