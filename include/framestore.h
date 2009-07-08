#ifndef FRAMESTORE_H
#define FRAMESTORE_H 1

#include "sixt.h"
#include "detectors.types.h"


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


#endif /* FRAMESTORE_H */

