#ifndef UNIFORMDETECTORBACKGROUND_H
#define UNIFORMDETECTORBACKGROUND_H 1


#include "sixt.h"
#include "sixt_random.h"
#include "impact.h"
#include "spectrum.h"
#include "photon.h"
#include "squarepixels.h"
#include "rmf.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Simple model for Cosmic X-ray detector background. */
typedef struct {
  /** Background spectrum. */
  Spectrum spectrum;
  /** Rate of detector background events. */
  float rate; 

  /** Next background event (Impact). */
  Impact nextImpact;

} UniformDetectorBackground;


struct UniformDetectorBackgroundParameters {
  char* spectrum_filename;
  float rate;
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Set up the configuration of the UniformDetectorBackground. */
int initUniformDetectorBackground(UniformDetectorBackground*, 
				  struct UniformDetectorBackgroundParameters*);

/** Clean up the UniformDetectorBackground data structure. Release
    allocated memory and call clean-up routines of underlying
    structures. */
int cleanupUniformDetectorBackground(UniformDetectorBackground*);

/** Create new UniformDetectorBackground event (Impact object). The
    background event is generated with random position on the
    detector, random energy according to spectrum, and random event
    time according to Poisson statistics and average background
    rate. */
int createUniformDetectorBackgroundImpact(UniformDetectorBackground*, SquarePixels*,
					  const struct ARF* const arf);


#endif /* UNIFORMDETECTORBACKGROUND_H */
