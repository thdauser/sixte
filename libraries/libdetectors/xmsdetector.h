#ifndef XMSDETECTOR_H
#define XMSDETECTOR_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "squarepixels.h"
#include "eventfile.h"
#include "xmseventfile.h"
#include "impactlist.h"


/** Model for the Transition Edge Sensor (TES) Calorimeter /
 * X-ray Microcalorimeter Spectrometer (XMS) detector on IXO.
 * This data structure contains the data required for the simulation of the XMS. 
 * It inherits some properties of the GenericDetector and SquarePixels data structures.
 * The XMSDetector data structure can be initialized by calling the initXMSDetector() function
 * with a XMSDetectorParameters data structure containing the desired setup.
 * A new photon Impact can be added to the XMSDetector array by the function 
 * addImpact2XMSDetector().
 * Finally after the simulation when the data structure is not required any more, the 
 * cleanupXMSDetector() routine should be called to release allocated memory and close open
 * file connections.
 */
typedef struct {

  /** Generic Detector properties like, e.g., the detector response. */
  GenericDetector generic;
  /** Array of square pixels. */
  SquarePixels pixels;


  /** Output event list. 
   * The events read out from the detector array are written to this event file that must
   * have the XMS-specific format. */
  XMSEventFile eventlist;

} XMSDetector;


/** Parameters of the XMSDetector model.  
 * This data structure contains the parameters for setting up the XMSDetector 
 * data structure.  
 * It is used as input for the initXMSDetector() routine.  
 * For documentation of the inidividual parameters see XMSDetector. */
struct XMSDetectorParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  char* eventlist_filename;
  char* eventlist_template;
};


////////////////////////////////////////////////////////


/** Set up the configuration of a XMSDetector. 
 * The routine is responsible to set up the initial the XMSDetector configuration which
 * is given in the XMSDetectorParameters data structure.
 * It has to take care of allocating the required memory for the pixel array and to 
 * create an event file for the output of the measured data.
 * For some of these tasks it simply calls the init routines of the underlying 
 * data structures. 
 */
int initXMSDetector(XMSDetector*, struct XMSDetectorParameters*);

/** Clean up the XMSDetector data structure. 
 * This routine should be called if the XMSDetector data structure is not required any more.
 * It takes care of releasing allocated memory and closes open file connections.
 * If applicable it calls clean-up routines of underlying data structures. 
 */
int cleanupXMSDetector(XMSDetector* wd);

/** Add a photon impact to the XMSDetector pixel array.
 * This is the standard routine to be called for the simulation of the XMSDetector.
 * For a new photon incident on the detector this routine determines the resulting
 * generated charge from the detector response and stores the event in the output event file.
 * Split events are taken into account based on a Gaussian charge cloud shape.
 */
int addImpact2XMSDetector(XMSDetector*, Impact*);


#endif /* XMSDETECTOR_H */
