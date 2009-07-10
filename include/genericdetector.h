#ifndef GENERICDETECTOR_H 
#define GENERICDETECTOR_H 1

#include "sixt.h"
#include "heasp.h"


/** Generic x-ray detector. Contains common data/specifications that are defined for
 * all different kinds of detectors. */
typedef struct {
  double ccsigma; /**< Charge cloud sigma [m]. This quantity is used to calculate size of 
		   * the charge cloud. */
  double ccsize; /**< Size of the charge cloud [m]. Defined as three times ccsigma. */

  long pha_threshold; // lower detector PHA threshold [PHA channels]
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


int initGenericDetector(GenericDetector*, struct GenericDetectorParameters*);


#endif /* GENERICDETECTOR_H */
