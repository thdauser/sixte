#ifndef HTRSDETECTOR_H
#define HTRSDETECTOR_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "eventfile.h"
//#include "htrseventfile.h"
#include "impactlist.h"


/** Model for the HTRS detector on IXO.
 * This data structure contains the data required for the simulation of the HTRS. 
 * It inherits some properties of the GenericDetector and SquarePixels data structures.
 * The HTRSDetector data structure can be initialized by calling the initHTRSDetector() function
 * with a HTRSDetectorParameters data structure containing the desired setup.
 * A new photon Impact can be added to the HTRSDetector array by the function 
 * addImpact2HTRSDetector().
 * Finally after the simulation when the data structure is not required any more, the 
 * cleanupHTRSDetector() routine should be called to release allocated memory and close open
 * file connections.
 */
typedef struct {

  /** Generic Detector properties like, e.g., the detector response. */
  GenericDetector generic;

} HTRSDetector;


/** Parameters of the HTRSDetector model.  
 * This data structure contains the parameters for setting up the HTRSDetector 
 * data structure.  
 * It is used as input for the initHTRSDetector() routine.  
 * For documentation of the inidividual parameters see HTRSDetector. 
 */
struct {
  struct GenericDetectorParameters generic;

  char* eventlist_filename;
  char* eventlist_template;
};


#endif /* HTRSDETECTOR_H */

