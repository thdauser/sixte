#ifndef GENERICDETECTOR_H 
#define GENERICDETECTOR_H 1

#include "sixt.h"
#include "random_sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

// The GNU Scientific Library Errorfunction is used to calculate charge 
// distribution of split events (assuming a Gaussian shape for the carge cloud).
#include <gsl/gsl_sf_erf.h>


#define INVALID_PIXEL (-1)   // flags an invalid pixel


/** Generic x-ray detector. Contains common data/specifications that are defined for
 * all different kinds of detectors. */
typedef struct {
  double ccsigma; /**< Charge cloud sigma [m]. This quantity is used to calculate size of 
		   * the charge cloud. */
  double ccsize; /**< Size of the charge cloud [m]. Defined as three times ccsigma. */

  long pha_threshold; /**< lower detector PHA threshold [PHA channels]. */
  float energy_threshold; /**< Lower detector energy threshold [keV]. */

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


/** Data required to initialize the GenericDetector data structure. */
struct GenericDetectorParameters {
  double ccsigma;

  long pha_threshold;
  float energy_threshold;

  char* rmf_filename;
};


////////////////////////////////////////////////////////////////////


/** Set up the initial configuraton for the GenericDetector data structure. */
int initGenericDetector(GenericDetector*, struct GenericDetectorParameters*);

/** Load an RMF/RSP matrix and the corresponding EBOUNDS from a response file. 
 * If the compiler flag '-DNORMALIZE_RMF' is set, the RSP is renormalized to an RMF in
 * such a way that the sum of each matrix row/column(?) is 1. */
struct RMF* loadRMF(char* filename, int* status);

/** Determines the PHA channel corresponding to a given energy according to the EBOUNDS
 * table of the detector response. */
long getChannel(float energy, struct RMF* rmf);

/** Determine the charge corresponding to a particular PHA channel according to 
 * the EBOUNDS table. */
float getEnergy(long channel, struct RMF* rmf);

/** Calculates the Gaussian integral using the GSL complementary error function. */
inline double gaussint(double x);


#endif /* GENERICDETECTOR_H */
