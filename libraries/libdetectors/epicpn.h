#ifndef EPICPN_H
#define EPICPN_H 1


#include "sixt.h"
#include "impact.h"
#include "eventfile.h"
#include "genericdetector.h"
#include "squarepixels.h"


// In XMM-Newton EPIC-pn timing mode all pixels in RAWY direction are
// collapsed to an linear image in RAWX direction. The time resolution 
// is 0.03 ms (pn, 64x200) or 1.75 ms (MOS, 100x600).


// If this flag is set, the charge cloud distribution in split events
// is determined by the exponential model proposed by Konrad Dennerl 
// for the eROSITA camera.
// If not, a Gaussian charge cloud model is assumed.
#define EPICpn_EXPONENTIAL_SPLITS 1

// If this flag is activated, the readout algorithm checks for split 
// events and marks them according to the XMM-Newton pattern naming 
// scheme.
#define EPICpn_DETECT_PATTERNS 1

// Define the value that is assigned to PAT_INF of invalid events.
// According to the eROSITA event file specification this should be '0',
// but as '0' is also assigned to singles, we choose '-1' instead.
#define EPICpn_INVALID_PATTERN (-1)

/*
// Maximum number of split pixels per readout frame.
#define MAX_N_SPLIT_LIST (1000)
*/


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** This data structure can be used for EPIC-pn like detectors. */
typedef struct {
  GenericDetector generic;
  SquarePixels pixels;

  /** Integration time for one detector line. */
  double step_time; 

  /** Current time. */
  double time; 

} EPICpn;


struct EPICpnParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double step_time;
  double t0;
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Set up the configuration of the EPIC-pn. The routine also calls
    the init routines of the underlying data structures. */
int initEPICpn(EPICpn* ep, struct EPICpnParameters* epp);

/** Clean up the EPICpn data structure. Release allocated memory and
    call clean-up routines of underlying structures. */
int cleanupEPICpn(EPICpn* ep);

/** This routine is called for readout of the EPIC-pn. The routine
    itself checks, whether a readout is necessary according to the
    current time and step_time. If it's time to do a readout, the
    routine initiates the readout process by calling
    readoutEPICpn(). */
int checkReadoutEPICpn(EPICpn* ep, double time);

/** Read out the EPIC-pn. The measured events are stored
    in an event list FITS file. */
inline int readoutEPICpn(EPICpn* ep);

/** Add a new photon impact to the EPIC-pn pixels. The generated
    charge is determined according to the detector response. If the
    charge cloud size is set, split events are calculated according to
    a Gaussian charge cloud model. The new charge is added to the
    charge already contained in the detector pixel, so pileup effects
    are taken into account. */
int addImpact2EPICpn(EPICpn* ep, Impact* impact);


#endif /* EPICPN_H */

