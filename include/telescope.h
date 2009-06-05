#ifndef TELESCOPE_H
#define TELESCOPE_H 1

#include "vector.h"

struct Telescope {
  double fov_diameter;    // diameter of the FOV (in [rad])
  double focal_length;    // focal length of the telescope

  double time;
  Vector r;              // position vector of the satellite
  Vector v;              // velocity vector of the satellite

  /* Unit vectors defining the telescope attitude */
  /** the direction of the telescope axis (in simplest case same direction as r) */
  Vector nz;
  /** usually points in the direction of motion, but is perpendicular to nz */
  Vector nx;
  /** perpendicular to nx and nz */
  Vector ny; 
};

#endif /* TELESCOPE_H */

