#ifndef GENDET_H 
#define GENDET_H 1

#include "sixt.h"
#include "arf.h"
#include "badpixmap.h"
#include "clocklist.h"
#include "codedmask.h"
#include "erodetbkgrndgen.h"
#include "event.h"
#include "eventlistfile.h"
#include "gendetline.h"
#include "genericdetector.h"
#include "genpixgrid.h"
#include "impact.h"
#include "point.h"
#include "psf.h"
#include "rmf.h"
#include "vignetting.h"
#include "xmlbuffer.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


typedef enum {
  GENDET_TIME_TRIGGERED  = 1,
  GENDET_EVENT_TRIGGERED = 2
} ReadoutTrigger;


typedef enum {
  GS_NONE        = 0,
  GS_GAUSS       = 1,
  GS_EXPONENTIAL = 2
} GenSplitType;


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic split event generator. */
typedef struct {

  /** Type of the GenSplit model. */
  GenSplitType type;

  /** Split parameter 1. For the Gaussian split model this parameter
      represents the charge cloud size in [m] at a photon energy of 1
      keV. For the exponential model this parameter represents the
      denominator in the exponential distribution term. */
  float par1;

  /** Split parameter 2. This value is only defined for the Gaussian
  charge cloud model. It represents the slope of the function giving
  the charge cloud size vs. photon energy. (ccs(E) = par1 + par2
  *sqrt(E)). */
  float par2;

} GenSplit;


/** Generic pixelized X-ray detector. The data structure is designed
    in a generic way. The characteristic properties for a particular
    detector are defined in a detector specific XML file. */
typedef struct {

  /** Telescope name, which is stored in the 'TELESCOP' header
      keyword in FITS files. */
  char* telescope;

  /** Detector pixel dimensions. */
  GenPixGrid* pixgrid;

  /** Array of pointers to pixel lines. */
  GenDetLine** line;

  /** Detector response matrix. The RMF file that is originally loaded
      may also contain ARF contributions. But as these already have
      been taken into account in the generation of the input spectra
      for the X-ray sources, the ARF contributions have to be removed
      by normalizing the RSP matrix. */
  struct RMF* rmf;
  /** Detector and telescope ARF containing the effective area. */
  struct ARF* arf;

  /** Telescope PSF. */
  PSF* psf;
  /** Telescope vignetting function. */
  Vignetting* vignetting;
  /** Focal length of the X-ray telescope [m]. */
  float focal_length;
  /** Diameter of the FoV [rad]. In the XML file the diameter is given
      in [deg], but it is converted to [rad] when parsing the XML
      file. */
  float fov_diameter;
  /** If the telescope is a coded mask telescope and not an imaging
      telescope, we have to specify the coded mask pattern file
      instead of the PSF. */
  CodedMask* coded_mask;

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

  /** Lower split threshold given in units of [keV]. This value is
      used in the pattern recognition algorithm. */
  float threshold_split_lo_keV;
  /** Lower split threshold given as a fraction of the charge in the
      main pixel. This value is used in the pattern recognition
      algorithm. */
  float threshold_split_lo_fraction;

  /** Split model. */
  GenSplit* split;

  /** Flag, whether a eroBackground model is available for the cosmic
      ray detector background. */
  int erobackground;

  /** Flag for detector readout trigger. The readout can be triggered
      either by an incoming photon event (GENDET_EVENT_TRIGGERED) or
      by a timing clock (GENDET_TIME_TRIGGERED). */
  ReadoutTrigger readout_trigger;

  /** List of clock operations for time-triggered detectors. */
  ClockList* clocklist;

  /** Time for one read-out frame [s]. */
  double frametime;

  /** Flag whether there has been any photon interaction since the
      last new frame. */
  int anyphoton;

  /** Charge transfer efficiency (CTE). In a line shift of the pixel
      array the charges in the shifted pixels are multiplied by this
      value in order to account for losses due to the shift. */
  float cte;

  /** Bad pixel map. */
  BadPixMap* badpixmap;

  /** File name (without path contributions) of the FITS file
      containing the XML detector definition. */
  char* filename;

  /** Path to the FITS file containing the XML detector definition. */
  char* filepath;

} GenDet;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new GenDet data structure and
    initializes it with the values from the specified XML definition
    file. */
GenDet* newGenDet(const char* const filename, int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenDet data structure to NULL. */
void destroyGenDet(GenDet** const det, int* const status);

/** Add a new photon impact to the detector. The function return value
    is the number of affected valid detector pixels. */
int addGenDetPhotonImpact(GenDet* const det, const Impact* const impact, 
			  EventListFile* const elf, int* const status);

/** Operate the time-triggered elements of the GenDet detector up to
    the specified point of time. */
void operateGenDetClock(GenDet* const det, EventListFile* const elf,
			const double time, int* const status);

/** Set the current detector time in the clocklist to the specified
    value. The default value for the start time is 0. */
void setGenDetStartTime(GenDet* const det, const double t0);

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
void GenDetReadoutLine(GenDet* const det, 
		       const int lineindex, 
		       const int readoutindex, 
		       EventListFile* const elf,
		       int* const status);

/** Clear a particular line of the GenDet pixel array. */
void GenDetClearLine(GenDet* const det, const int lineindex);

/** This function is called if a bad pixel is encountered and has to
    be applied to the detector pixel array. The parameter 'value' has
    to be added to the bad pixel at 'x' and 'y'. */
void encounterGenDetBadPix(void* const data, 
			   const int x, const int y, 
			   const float value);

/** Constructor for GenSplit data structure. */
GenSplit* newGenSplit(int* const status);

/** Destructor for GenSplit data structure. */
void destroyGenSplit(GenSplit** const split);

/** Determine split events for a particular photon impact and add the
    fractional charges to the affected pixels. The function return
    value is the number of valid affected pixels inside the detector,
    where charge has been added to. */
int makeGenSplitEvents(GenDet* const det,
		       const struct Point2d* const impact,
		       const float charge,
		       const long ph_id, const long src_id,
		       const double time,
		       EventListFile* const elf,
		       int* const status);

/** Parse the GenDet definition from an XML file. */
void parseGenDetXML(GenDet* const det, const char* const filename, 
		    int* const status);


#endif /* GENDET_H */
