#ifndef LADIMPACT_H 
#define LADIMPACT_H 1

#include "sixt.h"
#include "point.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Generic event on a pixelized X-ray detector. */
typedef struct {
  
  /** Panel ID. */
  long panel;

  /** Module ID. */
  long module;

  /** Element ID. */
  long element;

  /** Impact position of the photon on the element [m]. */
  struct Point2d position;

  /** Photon energy [keV]. */
  float energy;

  /** Arrival time of the photon [s]. */
  double time;

  /** Identifier of the photon. */
  long ph_id;

  /** Identifier of the corresponding source (defined in the SIMPUT
      source catalog). */
  long src_id;

} LADImpact;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor for LADImpact data structure. Initializes pointers with
    NULL and variables with their default values. */
LADImpact* getLADImpact(int* const status);

/** Destructor for LADImpact data structure. */
void freeLADImpact(LADImpact** const impact);


#endif /* LADIMPACT_H */
