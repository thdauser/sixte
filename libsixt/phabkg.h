#ifndef PHABKG_H 
#define PHABKG_H 1

#include "sixt.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic square pixel grid. */
typedef struct {
  /** Number of spectral bins / channels. */
  long nbins;

  /** PHA channel numbers. */
  long* channel;
  
  /** Background event rate distribution (counts/s/bin/cm^2). */
  float* distribution;

  /** GSL random number generator. */
  gsl_rng *randgen;

} PHABkg;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Loads the PHA spectrum required for the model from
    the specified file. The PHA file must contain two columns: CHANNEL
    and RATE. The entries in the CHANNEL column correspond to the
    channels in the EBOUNDS extension of the RMF. The entries in the
    RATE column represent the background event rate in this particular
    energy channel per second and per illuminated detector area
    (counts/s/bin/cm^2). */
PHABkg* newPHABkg(const char* const filename, int* const status);

/** Destructor. */
void destroyPHABkg(PHABkg** const phabkg);

/** Determine a number of background events according to the specified
    model. The number of events is determined by the area of the
    detector, the rates given in the input PHA spectrum, and the
    length of the regarded interval according to Poisson
    statistics. */
long* PHABkgGetEvents(const PHABkg* const phabkg, 
		      /** Regarded time interval in [s]. */
		      const double interval, 
		      /** Detector area in [m^2]. */
		      const float area,
		      unsigned int* const nevts, 
		      int* const status);


#endif /* PHABKG_H */
