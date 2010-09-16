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

  /** Lower and upper readout threshold in units of [keV]. These
      thresholds are applied in the read-out routine before converting
      the pixel charge to a PHA value. Pixel charges below this
      threshold will be discarded. */
  float threshold_readout_lo_keV, threshold_readout_up_keV;
  /** Lower and upper readout threshold in units of [PHA
      channel]. These thresholds are converted to the corresponding
      charge values [keV] according to the EBOUNDS table in the
      detector response file. */
  long threshold_readout_lo_PHA, threshold_readout_up_PHA;

  /** Lower primary event threshold given in unit of [keV]. This value
      is used in the pattern recognition algorithm to find pixels with
      a sufficient charge for a separate event. */
  float threshold_event_lo_keV;

  /** Lower split threshold given as a fraction of the charge in the main
      pixel. This value is used in the pattern recognition
      algorithm. */
  float threshold_split_lo_fraction;

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

/** Set the output event list file. The output file is newly generated
    from a given FITS template. */
void GenDetNewEventFile(GenDet* const det, const char* const filename, 
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
void GenDetReadoutLine(GenDet* const det, const int lineindex, 
		       const int readoutindex, int* const status);

/** Clear a particular line of the GenDet pixel array. */
void GenDetClearLine(GenDet* const det, const int lineindex);


#endif /* GENDET_H */
