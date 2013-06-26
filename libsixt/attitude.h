#ifndef ATTITUDE_H
#define ATTITUDE_H 1

#include "sixt.h"
#include "attitudefile.h"
#include "vector.h"
#include "telescope.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Entry of the Attitude collection of pointing directions. */
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


/** Collection containing the temporal evolution of the attitude. */
typedef struct {
  /** Number of AttituideEntry elements in the Attitude. */
  long nentries;

  /** Individual AttitudeEntry elements giving the attitude of the
      telescope at a particular point of time. */
  AttitudeEntry* entry;

  /** Index of the currently selected entry in the collection. */
  long current_entry;

  /** Alignment flag. If the rollangle should be determined with
      respect to the equatorial plane, the value of the alignment flag
      should be 0. If the alignment is with respect to the motion of
      the telescope axis, the value of the alignment flag should be
      1. */
  int alignment;

} Attitude;


/////////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////////


/** Constructor for the Attitude data structure. Allocates memory for
    the object. */
Attitude* getAttitude(int* const status);

/** Get a new Attitude object and load the data from the specified
    file. The routine loads the entire attitude data from the
    file. After reading it checks, whether the required time interval
    is a subset of the data provided in the attitude file. */
Attitude* loadAttitude(const char* filename, int* const status);

/** Destructor for the Attitude data structure. */
void freeAttitude(Attitude** const ac);

/** Determine the telescope pointing direction at a specific time. */
Vector getTelescopeNz(Attitude* const ac, 
		      const double time,
		      int* const status);

/** Determine the 3 axes vectors for the telescope coordinate
    system. */
void getTelescopeAxes(Attitude* const ac,
		      Vector* const nx,
		      Vector* const ny,
		      Vector* const nz,
		      const double time, 
		      int* const status);

/** Determine the roll-angle ([rad]) at a specific time. */
float getRollAngle(Attitude* const ac, 
		   const double time, 
		   int* const status);

/** Return an empty AttitudeEntry object with default values. */
AttitudeEntry defaultAttitudeEntry();


#endif /* ATTITUDE_H */

