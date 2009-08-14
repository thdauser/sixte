#ifndef HTRSDETECTOR_H
#define HTRSDETECTOR_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "hexagonalpixels.h"
#include "eventfile.h"
#include "htrseventfile.h"
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
  /** Array of hexagonal pixels. */
  HexagonalPixels pixels;

  /** Dead time of a pixel after an event detection.
   * If a new photon arrives during the dead time after a previous event in the same pixel,
   * the new photon will not be detected. This is a model for the physical
   * readout and clearing process of the detector pixel. */
  double dead_time;

  /** Output event list. 
   * The events read out from the detector array are written to this event file that must
   * have the HTRS-specific format. */
  HTRSEventFile eventlist;

} HTRSDetector;


/** Parameters of the HTRSDetector model.  
 * This data structure contains the parameters for setting up the HTRSDetector 
 * data structure.  
 * It is used as input for the initHTRSDetector() routine.  
 * For documentation of the inidividual parameters see HTRSDetector. 
 */
struct HTRSDetectorParameters {
  struct GenericDetectorParameters generic;
  struct HexagonalPixelsParameters pixels;

  double dead_time;

  char* eventlist_filename;
  char* eventlist_template;
};


////////////////////////////////////////////////////////


/** Set up the configuration of a HTRSDetector. 
 * The routine is responsible to set up the initial the HTRSDetector configuration which
 * is given in the HTRSDetectorParameters data structure.
 * It has to take care of allocating the required memory for the pixel array and to 
 * create an event file for the output of the measured data.
 * For some of these tasks it simply calls the init routines of the underlying 
 * data structures. 
 */
int initHTRSDetector(HTRSDetector*, struct HTRSDetectorParameters*);

/** Clean up the HTRSDetector data structure. 
 * This routine should be called when the HTRSDetector data structure 
 * is not required any more.
 * It takes care of releasing allocated memory and closes open file connections.
 * If applicable it calls clean-up routines of underlying data structures. 
 */
int cleanupHTRSDetector(HTRSDetector*);

/** Add a photon impact to the HTRSDetector pixel array.
 * This is the standard routine to be called for the simulation of the HTRSDetector.
 * For a new photon incident on the detector this routine determines the resulting
 * generated charge from the detector response and stores the event in the output event file.
 * Split events are taken into account based on a Gaussian charge cloud shape.
 */
int addImpact2HTRSDetector(HTRSDetector*, Impact*);


#endif /* HTRSDETECTOR_H */

