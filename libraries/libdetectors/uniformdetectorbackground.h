#ifndef UNIFORMDETECTORBACKGROUND_H
#define UNIFORMDETECTORBACKGROUND_H 1


#include "sixt.h"
#include "sixt_random.h"
#include "impactlist.h"
#include "spectrum.h"
#include "photon.h"
#include "squarepixels.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


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


////////////////////////////////////////////////////////


/** Set up the configuration of the UniformDetectorBackground. */
int initUniformDetectorBackground(UniformDetectorBackground*, 
				  struct UniformDetectorBackgroundParameters*);

/** Clean up the UniformDetectorBackground data structure. 
 * Release allocated memory and call clean-up routines of underlying structures. */
int cleanupUniformDetectorBackground(UniformDetectorBackground*);


int createUniformDetectorBackgroundImpact(UniformDetectorBackground*, SquarePixels*,
					  struct RMF*);


#endif /* UNIFORMDETECTORBACKGROUND_H */
