#ifndef GENERICDETECTOR_H 
#define GENERICDETECTOR_H 1

#include "sixt.h"
#include "sixt_random.h"
#include "gaussianchargecloud.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

// The GNU Scientific Library Errorfunction is used to calculate charge 
// distribution of split events (assuming a Gaussian shape for the carge cloud).
#include <gsl/gsl_sf_erf.h>


#define INVALID_PIXEL (-1)   // flags an invalid pixel


/** Generic x-ray detector. Contains common data/specifications that
    are defined for all different kinds of detectors.  The
    GenericDetector data structure is usually included as element of
    more specific detector models like, e.g. WFIDetector or
    FramestoreDetector.*/
typedef struct {

  /** Properties of Gaussian shaped charge clouds. */
  GaussianChargeCloud gcc;

  /** Lower detector PHA threshold [PHA channels].  Events with a
      lower PHA value are dismissed.  If the PHA threshold is set to
      "-1" the energy threshold is taken into account. */
  long pha_threshold;
  /** Lower detector energy threshold [keV]. 
   * Events with a lower energy / pixel charge are dismissed. 
   * This value is only regarded if the pha_threshold is set to "-1". */
  float energy_threshold; 

  /** Detector response matrix. Includes the RMF and the detector-specific
   * response elements like filter transmission and quantum efficiency.
   * So the sum of each line of the response matrix HAS TO BE less
   * or equal 1 (sum <= 1) !
   * The RMF can be normalized explicitly to be a real RMF without 
   * photon loss due response effects by setting the compiler flag 
   * "-DNORMALIZE_RMF".
   * The mirror specific elements are treated SEPARATELY in the photon
   * imaging process. */
  struct RMF* rmf;

} GenericDetector;


/** Data required to initialize the GenericDetector data structure.
    The meaning of the inidividual parameters is the same as in the
    GenericDetector data structure. */
struct GenericDetectorParameters {
  double ccsigma;

  long pha_threshold;
  float energy_threshold;

  char* rmf_filename;
};


////////////////////////////////////////////////////////////////////


/** Set up the initial configuraton for the GenericDetector data
    structure.  This routines sets the values of the GenericDetector
    elements according to the given parameters. */
int initGenericDetector(GenericDetector*, struct GenericDetectorParameters*);

/** Load an RMF/RSP matrix and the corresponding EBOUNDS from a
    response file.  If the compiler flag '-DNORMALIZE_RMF' is set, the
    RSP is renormalized to an RMF in such a way that the sum of each
    matrix row/column(?) is 1. */
struct RMF* loadRMF(char* filename, int* status);

/** Determines the PHA channel corresponding to a given energy
    according to the EBOUNDS table of the detector response.  The
    routine performs a binary search to obtain the PHA channel the
    specified energy lies in. The energy has to be given in the same
    unit as the EBOUNDS are.  That is usually [keV]. */
long getChannel(float energy, struct RMF* rmf);

/** Determine the charge corresponding to a particular PHA channel
    according to the EBOUNDS table.  The input channel must have the
    same offset as in the EBOUNDS table. I.e. if the first channel in
    the EBOUNDS has the number 1, the numbering starts at 1. If the
    first channel has the number 0, the numbering starts at 0.  The
    returned energy is given in the same units as the EBOUNDS. That is
    usually [keV]. */
float getEnergy(long channel, struct RMF* rmf);

/** Calculates the Gaussian integral using the GSL complementary error
    function. */
inline double gaussint(double x);


#endif /* GENERICDETECTOR_H */
