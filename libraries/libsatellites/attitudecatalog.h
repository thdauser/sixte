#ifndef ATTITUDECATALOG_H
#define ATTITUDECATALOG_H 1

#include "sixt.h"
#include "attitudefile.h"
#include "vector.h"
#include "telescope.h"


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



/** Constructor for the AttitudeCatalog. */
AttitudeCatalog* get_AttitudeCatalog(const char attitude_filename[],
				     double t0, double timespan, int* status);

/** Destructor for the AttitudeCatalog. */
void free_AttitudeCatalog(AttitudeCatalog* ac);


#endif /* ATTITUDECATALOG_H */

