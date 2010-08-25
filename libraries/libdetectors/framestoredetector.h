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
#define FD_EXPONENTIAL_SPLITS 1

// Define the value that is assigned to PAT_INF of invalid events.
// According to the eROSITA event file specification this should be '0',
// but as '0' is also assigned to singles, we choose '-1' instead.
#define FD_INVALID_PATTERN (-1)

// Maximum number of split pixels per readout frame.
#define FD_MAX_N_SPLIT_LIST (1000)


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


typedef struct {
  GenericDetector generic;
  SquarePixels* pixels;

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

  /** The split threshold gives the fraction of the energy of the
      initially detected event which is used as threshold for the
      search of split partners in the surrounding pixels. According to
      the readout strategy proposed by Konrad Dennerl each frame is
      analysed with a particular fixed event threshold given in keV or
      PHA. If a pixel with a charge above this event threshold is
      found, the neighboring pixels are investigated with a lower
      so-called split threshold, which is a particular percentage of
      the charge in the initially detected pixel. Therefore the value
      of the variable "split_threshold" is a fraction between 0 and
      1. */
  float split_threshold;

  /** Flag, whether split events are simulated or not. */
  int make_splits;

} FramestoreDetector;


struct FramestoreDetectorParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double integration_time;
  double t0;
  float split_threshold;
  int make_splits;

  char* eventlist_filename;
  char* eventlist_template;  
};


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Set up the configuration of a FramestoreDetector. The routine
    calls the init routines of the underlying data structures. */
FramestoreDetector* newFramestoreDetector(struct FramestoreDetectorParameters* fdp,
					  int* status);

/** Destroy the FramestoreDetector data structure. Release allocated
    memory and call destructor routines of underlying structures. */
int destroyFramestoreDetector(FramestoreDetector** fd);

/** This routine is called for readout of the FramestoreDetector. The
    routine itself checks, whether a readout is necessary according to
    the current time and readout time. If it's time to do a readout,
    the routine initiates the readout process by calling
    readoutFramestoreDetector(). */
int checkReadoutFramestoreDetector(FramestoreDetector* fd, double time);

/** Read out the FramestoreDetector. The measured events are stored
    in an event list FITS file. */
inline int readoutFramestoreDetector(FramestoreDetector* fd);

/** Add a new photon impact to the FramestoreDetector pixels.  The
    generated charge is determined according to the detector response.
    If the charge cloud size is set, split events are calculated
    according to a Gaussian charge cloud model.  The new charge is
    added to the charge already contained in the detector pixel, so
    pileup effects are taken into account. */
int addImpact2FramestoreDetector(FramestoreDetector* fd, Impact* impact);

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

