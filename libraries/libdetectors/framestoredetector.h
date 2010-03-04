#ifndef FRAMESTOREDETECTOR_H
#define FRAMESTOREDETECTOR_H 1


#include "sixt.h"
#include "impact.h"
#include "eventfile.h"
#include "erositaeventfile.h"
#include "erositaevent.h"
#include "genericdetector.h"
#include "squarepixels.h"


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


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


typedef struct {
  GenericDetector generic;
  SquarePixels pixels;

  /** Integration time of the entire pnCCD (!) detector array. (=
      Span of time between 2 subsequent readouts). */
  double integration_time; 

  /** Current readout time. The end of the integration time /
      beginning of dead time. */
  double readout_time; 

  /** Number of the current frame. */
  long frame; 

  /** Event list FITS file for the eROSITA-specific events. */
  eROSITAEventFile eventlist; 

} FramestoreDetector;


struct FramestoreDetectorParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double integration_time;
  double t0;
  char* eventlist_filename;
  char* eventlist_template;  
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Set up the configuration of a FramestoreDetector. The routine
    also calls the init routines of the underlying data structures. */
int initFramestoreDetector(FramestoreDetector* fd, 
			   struct FramestoreDetectorParameters* fdp);

/** Clean up the FramestoreDetector data structure. Release allocated
    memory and call clean-up routines of underlying structures. */
int cleanupFramestoreDetector(FramestoreDetector*);

/** This routine is called for readout of the FramestoreDetector. The
    routine itself checks, whether a readout is necessary according to
    the current * time and readout time. If it's time to do a readout,
    the routine initiates the * readout process by calling
    readoutFramestoreDetector(). */
int checkReadoutFramestoreDetector(FramestoreDetector*, double time);

/** Read out the FramestoreDetector. The measured events are stored
    in an event list FITS file. */
inline int readoutFramestoreDetector(FramestoreDetector*);

/** Add a new photon impact to the FramestoreDetector pixels.  The
    generated charge is determined according to the detector response.
    If the charge cloud size is set, split events are calculated
    according to a Gaussian charge cloud model.  The new charge is
    added to the charge already contained in the detector pixel, so
    pileup effects are taken into account. */
int addImpact2FramestoreDetector(FramestoreDetector*, Impact*);


/** Marker routine for split partners. The routine scans the pixels
    around a certain event on the FramestoreDetector in order to find
    neighboring events belonging to the same split pattern. */
void fdMarkEvents(eROSITAEvent* list, int* nlist, 
		  int* maxidx, int* minidx,
		  FramestoreDetector* fd, 
		  int x, int y);


/** Identify split patterns according to the eROSITA definition given
    in the event file specification. The event list has to consist of
    all events belonging to the same pattern. The second parameter
    gives the number of entries in the list, i.e., the number of
    events involved in the split pattern, and the third and forth
    parameter give the indeces of the events with maximum and minimum
    energy value in the list. */
void fdPatternIdentification(eROSITAEvent* list, const int nlist, 
			     const int maxidx, const int minidx);


#endif /* FRAMESTOREDETECTOR_H */

