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

  /** Defines the detector x-direction. The x-axis doesn't necessarily
      have to point in the direction of the telescope motion, but can
      be distorted by the roll-angle. */
  Vector nx;

  /** Roll-angle ([rad]). */
  float roll_angle;

  // TODO Keep either the roll_angle or nx. 

} AttitudeEntry;


/** Catalog containing the attitude information for an X-ray
    telescope. */
typedef struct {
  /** Number of AttituideEntry elements in the AttitudeCatalog. */
  long nentries;

  /** Individual AttitudeEntry elements giving the attitude of the
      telescope at a particular point of time. */
  AttitudeEntry* entry;

  /** Index of the currently selected entry in the catalog. */
  long current_entry;

  /** Alignment flag. If the rollangle should be determined with
      respect to the equatorial plane, the value of the alignment flag
      should be 0. If the alignment is with respect to the motion of
      the telescope axis, the value of the alignment flag should be
      1. */
  int alignment;

} AttitudeCatalog;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Constructor for the AttitudeCatalog. Allocate memory for the
    object. */
AttitudeCatalog* getAttitudeCatalog(int* const status);

/** Get a new AttitudeCatalog object and load the data from the
    specified file. The routine loads the entire attitude data from
    the file. After reading it checks, whether the required time
    interval is a subset of the data provided in the attitude file. */
AttitudeCatalog* loadAttitudeCatalog(const char* filename, 
				     int* const status);

/** Destructor for the AttitudeCatalog. */
void freeAttitudeCatalog(AttitudeCatalog** const ac);

/** Determine the telescope pointing direction at a specific time. */
Vector getTelescopeNz(AttitudeCatalog* const ac, 
		      const double time,
		      int* const status);

/** Determine the 3 axes vectors for the telescope coordinate
    system. */
void getTelescopeAxes(AttitudeCatalog* const ac,
		      Vector* const nx,
		      Vector* const ny,
		      Vector* const nz,
		      const double time, 
		      int* const status);

/** Determine the roll-angle ([rad]) at a specific time. */
float getRollAngle(AttitudeCatalog* const ac, 
		   const double time, 
		   int* const status);

/** Return an empty AttitudeEntry object with default values. */
AttitudeEntry defaultAttitudeEntry();


#endif /* ATTITUDECATALOG_H */

