#ifndef ATTITUDE_H
#define ATTITUDE_H 1

#include "sixt.h"
#include "attitudefile.h"
#include "vector.h"
#include "telescope.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


typedef enum {
  /** Alignment of the nx vector of the RAWX/RAWY coordinate system
      along the north direction. */
  ATTNX_NORTH =0,

  /** Alignment of the nx vector along the direction of motion of the
      pointing axis. */
  ATTNX_MOTION=1,

} AttNxAlign;


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////%///////////////////////////////////////////


/** Entry of the Attitude collection of pointing directions. */
typedef struct {
  /** Point of time for which this attitude is valid. */
  double time;

  /** Telescope pointing direction. */
  Vector nz;

  /** Roll-angle ([rad]). */
  float roll_angle;

} AttitudeEntry;


/** Collection containing the temporal evolution of the attitude. */
typedef struct {
  /** Number of AttituideEntry elements in the Attitude. */
  long nentries;

  /** Individual AttitudeEntry elements giving the attitude of the
      telescope at a particular point of time. */
  AttitudeEntry* entry;

  /** Index of the currently selected entry in the collection. */
  long currentry;

  /** Alignment flag. Determines the reference direction for the
      alignment of the nx vector. The nx vector is determined by the
      rotation of this reference direction around the roll angle. */
  AttNxAlign align;

  /** MJDREF as specified in the FITS header of the attitude file. */
  double mjdref;

  /** TSTART and TSTOP. */
  double tstart, tstop;

} Attitude;


////////////////////////////////////////////////////////////////
// Function Declarations.
////////////////////////////////////////////////////////////////


/** Constructor for the Attitude data structure. Allocates memory for
    the object. */
Attitude* getAttitude(int* const status);

/** Get a new Attitude object and load the data from the specified
    file. The routine loads the entire attitude data from the
    file. After reading it checks, whether the required time interval
    is a subset of the data provided in the attitude file. */
Attitude* loadAttitude(const char* const filename, int* const status);

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

/** Produce a pointing attitude. */
Attitude* getPointingAttitude(const double mjdref,
			      const double tstart,
			      const double tstop,
			      const double ra,
			      const double dec,
			      int* const status);

/** Check if the interval specified by tstart and tstop with respect
    to mjdred is covered by the attitude. If that is not the case, an
    error is returned. */
void checkAttitudeTimeCoverage(const Attitude* const ac,
			       const double mjdref,
			       const double tstart,
			       const double tstop,
			       int* const status);


#endif /* ATTITUDE_H */

