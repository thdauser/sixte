#ifndef PHABKG_H 
#define PHABKG_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic square pixel grid. */
typedef struct {
  /** Number of spectral bins / channels. */
  long nbins;

  /** PHA channel numbers. */
  long* channel;

  /** Background event rate distribution. The units are either
      [counts/s/bin/deg^2] or [counts/s/bin/m^2] depending on
      whether the background model is assigned to a telescope or a
      detector. */
  float* distribution;

  /** Time of the next background event. */
  double tnext;

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
    [counts/s/bin/m^2] or per sky angle [counts/s/bin/deg^2],
    depending on whether the model is assigned to a telescope or a
    detector. */
PHABkg* newPHABkg(const char* const filename, int* const status);

/** Destructor. */
void destroyPHABkg(PHABkg** const phabkg);

/** Determine the PHA value and time of a background event according
    to the specified spectral distribution. The function returns an
    individual event. The time differences between the events are
    exponentially distributed. The average rate is determined by the
    distribution given in the PHA data set multiplied with the
    scaling factor. */
int getPHABkgEvent(PHABkg* const phabkg,
		   /** Scaling factor for the count rate
		       distribution. Must be given in [m^2] or
		       [deg^2]. */
		   const float scaling,
		   const double tstart,
		   /** Upper limit for the time of the background
		       event [s]. */
		   const double tstop,
		   double* const t,
		   long* const pha,
		   int* const status);


#endif /* PHABKG_H */
