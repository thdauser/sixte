#ifndef TELESCOPE_H
#define TELESCOPE_H 1

#include "vector.h"

struct Telescope {
  /* Unit vectors defining the telescope attitude. */

  /** Direction of the telescope axis (in simplest case same direction
      as r). */
  Vector nz;
  /** Usually points in the direction of motion, but is perpendicular
      to nz. */
  Vector nx;
  /** Perpendicular to nx and nz */
  Vector ny; 
};


#endif /* TELESCOPE_H */

