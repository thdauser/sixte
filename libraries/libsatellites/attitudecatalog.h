#ifndef ATTITUDECATALOG_H
#define ATTITUDECATALOG_H 1

#include "sixt.h"
#include "attitudefile.h"
#include "vector.h"
#include "telescope.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Entry of the AttitudeCatalog. */
typedef struct {
  /** Point of time for which this attitude is valid. */
  double time;
  /** Telescope pointing direction. */
  Vector nz;
  /** Defines the detector x-direction.
   * The x-axis doesn't necessarily have to point in the direction of the telescope
   * motion, but can be distorted by the roll-angle. */
  Vector nx;
} AttitudeEntry;

/** Catalog containing the attitude information for an X-ray telescope. */
typedef struct {
  /** Number of AttituideEntry elements in the AttitudeCatalog. */
  long nentries;
  /** Individual AttitudeEntry elements giving the attitude of the telescope
   * at a particular point of time. */
  AttitudeEntry* entry;
} AttitudeCatalog;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Constructor for the AttitudeCatalog. Load data from the specified
    file for the time interval from t0 to t0+timespan. */
AttitudeCatalog* get_AttitudeCatalog(const char* filename, double t0, double timespan, 
				     int* status);
/** Constructor for the AttitudeCatalog. Load the entire attitude
    catalog from the specified file (not only within a certain
    interval). */
AttitudeCatalog* getEntireAttitudeCatalog(const char* filename, int* status);

/** Destructor for the AttitudeCatalog. */
void free_AttitudeCatalog(AttitudeCatalog* ac);


#endif /* ATTITUDECATALOG_H */

