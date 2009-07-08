#ifndef FRAMESTORE_H
#define FRAMESTORE_H 1

#include "sixt.h"
#include "detectors.h"


/** Properties specific for a Framestore CCD detector. */
struct FramestoreProperties {
  double integration_time; /**< Integration time of the entire pnCCD (!) 
			    * detector array.
			    * (= Span of time between 2 subsequent readouts). */
  long frame; /**< Number of the current frame. */

  double ccsigma; /**< Charge cloud sigma [m]. This quantity is used to calculate size of 
		   * the charge cloud. */
  double ccsize; /**< Size of the charge cloud [m]. Defined as three times ccsigma. */

};


struct FramestoreParameters {
  double integration_time;
  double ccsigma;
};

////////////////////////////////////////////////////////


/** Setup the configuration of a framestore detector. */
int init_FramestoreDetector(Detector*, struct FramestoreParameters);

/** Read out the framestore detector IF necessary. 
 * If the integration interval since the last readout operation is exceeded,
 * the entire pixel array is read out and the events are stored in an event
 * file. */
inline void readout_FramestoreDetector(void*, double time, struct Eventlist_File*, 
				       int *status);


#endif /* FRAMESTORE_H */

