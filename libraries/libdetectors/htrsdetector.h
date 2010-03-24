#ifndef HTRSDETECTOR_H
#define HTRSDETECTOR_H 1

#include "sixt.h"
#include "genericdetector.h"

#ifdef HTRS_HEXPIXELS
#define HTRS_ANYPIXELS 1
#include "hexagonalpixels.h"
#endif
#ifdef HTRS_ARCPIXELS
#define HTRS_ANYPIXELS 1
#include "arcpixels.h"
#endif

// Check whether any of the possible detector pixel layouts (hexagonal
// or arc) has been selected. If not, break with an error message.
#ifndef HTRS_ANYPIXELS
#error "Error: No pixel type for HTRS detector selected!"
#endif

#include "eventfile.h"
#include "htrseventfile.h"
#include "impact.h"


/** Model for the HTRS detector on IXO. This data structure contains
    the data required for the simulation of the HTRS. It inherits some
    properties of the GenericDetector and SquarePixels data
    structures.  The HTRSDetector data structure can be initialized by
    calling the initHTRSDetector() function with a
    HTRSDetectorParameters data structure containing the desired
    setup. A new photon Impact can be added to the HTRSDetector array
    by the function addImpact2HTRSDetector(). Finally after the
    simulation when the data structure is not required any more, the
    cleanupHTRSDetector() routine should be called to release
    allocated memory and close open file connections.
 */
typedef struct {

  /** Generic Detector properties like, e.g., the detector response. */
  GenericDetector generic;

#ifdef HTRS_HEXPIXELS
  /** Array of hexagonal pixels. */
  HexagonalPixels pixels;
#endif
#ifdef HTRS_ARCPIXELS
  /** Array of ArcPixels. */
  ArcPixels pixels;
#endif

  /** Shaping time for a pulse. The pixel is insensitive to further
      photons during that time span after a photon detection. If a new
      photon arrives during the shaping time after a previous event in
      the same pixel, the new photon will not be detected. This is a
      model for the physical readout and clearing process of the
      detector pixel. */
  double shaping_time;

  /** Time required to reset a HTRS detector pixel. When the output
      voltage of the pixel exceeds the input range of the subsequent
      electronics, the collected charge in the pixel has to be
      deleted, which is called reset. During the time required for the
      reset, the pixel is insensitive for incident photons. */
  double reset_time;

  /** Output event list. The events read out from the detector array
      are written to this event file that must have the HTRS-specific
      format. */
  HTRSEventFile eventlist;

  /** Number of events. */
  long nevents;
  /** Number of single events. */
  long nsingles;
  /** Number of double split events. */
  long ndoubles; 

} HTRSDetector;


/** Parameters of the HTRSDetector model. This data structure
    contains the parameters for setting up the HTRSDetector data
    structure. It is used as input for the initHTRSDetector()
    routine. For documentation of the inidividual parameters see
    HTRSDetector. */
struct HTRSDetectorParameters {
  struct GenericDetectorParameters generic;

#ifdef HTRS_HEXPIXELS
  struct HexagonalPixelsParameters pixels;
#endif
#ifdef HTRS_ARCPIXELS
  struct ArcPixelsParameters pixels;
#endif

  double shaping_time;
  double reset_time;

  char* eventlist_filename;
  char* eventlist_template;
};


////////////////////////////////////////////////////////


/** Set up the configuration of a HTRSDetector.  The routine is
    responsible to set up the initial the HTRSDetector configuration
    which is given in the HTRSDetectorParameters data structure.  It
    has to take care of allocating the required memory for the pixel
    array and to create an event file for the output of the measured
    data.  For some of these tasks it simply calls the init routines
    of the underlying data structures.
 */
int initHTRSDetector(HTRSDetector* hd, struct HTRSDetectorParameters* parameters);

/** Clean up the HTRSDetector data structure.  This routine should be
    called when the HTRSDetector data structure is not required any
    more.  It takes care of releasing allocated memory and closes open
    file connections.  If applicable it calls clean-up routines of
    underlying data structures.
 */
int cleanupHTRSDetector(HTRSDetector* hd);

/** Add a photon impact to the HTRSDetector pixel array.  This is the
    standard routine to be called for the simulation of the
    HTRSDetector.  For a new photon incident on the detector this
    routine determines the resulting generated charge from the
    detector response and stores the event in the output event file.
    Split events are taken into account based on a Gaussian charge
    cloud shape.
 */
int addImpact2HTRSDetector(HTRSDetector* hd, Impact* impact);


#endif /* HTRSDETECTOR_H */

