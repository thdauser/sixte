#ifndef FRAMESTORE_H
#define FRAMESTORE_H 1

#include "sixt.h"
#include "detectors.h"


/** Properties specific for a Framestore CCD detector. */
struct FramestoreProperties {
  double integration_time; /**< Integration time of the entire pnCCD (!) 
			    * detector array.
			    * (= Span of time between 2 subsequent readouts). */
  //  long frame; /**< Number of the current frame. */
};


struct FramestoreParameters {
  double integration_time;
};

////////////////////////////////////////////////////////


/** Setup the configuration of a framestore detector. */
int init_FramestoreDetector(Detector*, struct DetectorParameters, struct FramestoreParameters);

/** Read out the framestore detector IF necessary. 
 * If the integration interval since the last readout operation is exceeded,
 * the entire pixel array is read out and the events are stored in an event
 * file. 
 * After the readout action the pixel array is cleared. */
void readout_FramestoreDetector(void*, double time, struct Eventlist_File*, 
				int *status);

/** Add a new photon impact to the pixels of the CCD.
 * The resulting charge is determined from the photon energy and the detector response.
 * If the selected charge cloud size is greater than 0, split events are taken into account.
 * The new charge is added to the existing charge in the pixels, i.e. pile up is also taken
 * into account. */
void add_Impact2FramestoreDetector(void* det, struct Impact* impact);


#endif /* FRAMESTORE_H */

