#ifndef TIMINGCCD_H
#define TIMINGCCD_H 1


#include "sixt.h"
#include "impact.h"
#include "eventfile.h"
#include "genericdetector.h"
#include "squarepixels.h"


/*
// If this flag is set, the charge cloud distribution in split events
// is determined by the exponential model proposed by Konrad Dennerl.
// If not, a Gaussian charge cloud model is assumed.
#define EXPONENTIAL_SPLITS 1

// If this flag is activated, the framestore readout algorithm checks
// for split events and marks them according to the eROSITA pattern
// naming scheme.
#define FD_DETECT_PATTERNS 1

// Define the value that is assigned to PAT_INF of invalid events.
// According to the eROSITA event file specification this should be '0',
// but as '0' is also assigned to singles, we choose '-1' instead.
#define FD_INVALID_PATTERN (-1)

// Maximum number of split pixels per readout frame.
#define MAX_N_SPLIT_LIST (1000)
*/

/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


typedef struct {
  GenericDetector generic;
  SquarePixels pixels;

  /** Integration time for one read-out step (one detector line). */
  double step_time; 

  /** Current time. */
  double time; 

} TimingCCD;


struct TimingCCDParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double step_time;
  double t0;
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Set up the configuration of a TimingCCD. The routine also calls
    the init routines of the underlying data structures. */
int initTimingCCD(TimingCCD* tc, struct TimingCCDParameters* tcp);

/** Clean up the TimingCCD data structure. Release allocated
    memory and call clean-up routines of underlying structures. */
int cleanupTimingCCD(TimingCCD* tc);

/** This routine is called for readout of the TimingCCD. The routine
    itself checks, whether a readout is necessary according to the
    current time and step_time. If it's time to do a readout, the
    routine initiates the readout process by calling
    readoutTimingCCD(). */
int checkReadoutTimingCCD(TimingCCD* tc, double time);

/** Read out the TimingCCD. The measured events are stored
    in an event list FITS file. */
inline int readoutTimingCCD(TimingCCD* tc);

/** Add a new photon impact to the TimingCCD pixels. The generated
    charge is determined according to the detector response. If the
    charge cloud size is set, split events are calculated according to
    a Gaussian charge cloud model. The new charge is added to the
    charge already contained in the detector pixel, so pileup effects
    are taken into account. */
int addImpact2TimingCCD(TimingCCD* tc, Impact* impact);


#endif /* TIMINGCCD_H */

