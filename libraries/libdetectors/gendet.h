#ifndef GENDET_H 
#define GENDET_H 1

#include "sixt.h"
#include "gendetline.h"
#include "genpixgrid.h"
#include "gensplit.h"
#include "genevent.h"
#include "geneventfile.h"
#include "genericdetector.h"
#include "impact.h"
#include "clocklist.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/////////////////////////////////////////////////////////////////
// Constants
/////////////////////////////////////////////////////////////////


typedef enum {
  GENDET_TIME_TRIGGERED  = 1,
  GENDET_EVENT_TRIGGERED = 2
} ReadoutTrigger;


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic pixelized X-ray detector. The data structure is designed
    in a generic way. The characteristic properties for a particular
    detector are defined in a detector specific XML file. */
typedef struct {

  /** Detector pixel dimensions. */
  GenPixGrid* pixgrid;

  /** Array of pointers to pixel lines. */
  GenDetLine** line;

  /** Detector response matrix. The RSP file that is originally loaded
      may also contain ARF contributions. But as they already have
      been taken into account in the generation of the input spectra
      for the X-ray sources, the ARF contributions have to be removed
      by normalizing the RSP matrix. */
  struct RMF* rmf;

  /** Lower and upper threshold in units of [keV]. These thresholds
      are applied before converting the pixel charge to a PHA
      value. */
  float lo_keV_threshold, up_keV_threshold;
  /** Lower and upper PHA threshold. These thresholds are applied
      after converting the pixel charge to a PHA value. */
  long lo_PHA_threshold, up_PHA_threshold;

  /** Split model. */
  GenSplit* split;

  /** Flag for detector readout trigger. The readout can be triggered
      either by an incoming photon event (GENDET_EVENT_TRIGGERED) or
      by a timing clock (GENDET_TIME_TRIGGERED). */
  ReadoutTrigger readout_trigger;

  /** List of clock operations for time-triggered detectors. */
  ClockList* clocklist;

  /** Charge transfer efficiency (CTE). In a line shift of the pixel
      array the charges in the shifted pixels are multiplied by this
      value in order to account for losses due to the shift. */
  float cte;

  /** Event file for the output of the detected events. */
  GenEventFile* eventfile;
  /** File name of the template for the event list FITS file. */
  char eventfile_template[MAXMSG];

} GenDet;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenDet data structure. */
GenDet* newGenDet(const char* const filename, int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenDet data structure to NULL. */
void destroyGenDet(GenDet** const det, int* const status);

/** Set the output event list file. */
void GenDetSetEventFile(GenDet* const det, const char* const filename, 
			int* const status);

/** Add a new photon impact to the detector. */
void addGenDetPhotonImpact(GenDet* const det, const Impact* const impact, int* const status);

/** Operate the time-triggered elements of the GenDet detector up to
    the specified point of time. */
void operateGenDetClock(GenDet* const det, const double time, int* const status);

/** Shift the lines of the GenDet detector pixel array by one line
    into the direction of the read-out node in line 0. The charges in
    line 1 are added to the charges in line 0, such that the content
    of line 0 is not lost. The charges of all shifted pixels are
    multiplied with the CTE factor in order to account for possible
    losses. */
void GenDetLineShift(GenDet* const det);

/** Read-out a particular line of the GenDet pixel array and store the
    charges in the output event file. After read-out the charges
    in the pixels are deleted. */
void GenDetReadoutLine(GenDet* const det, const int lineindex, const int readoutindex, 
		       const double time, int* const status);


#endif /* GENDET_H */
