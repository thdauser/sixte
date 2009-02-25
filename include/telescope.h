#ifndef TELESCOPE_H
#define TELESCOPE_H 1

#include "vector.h"

struct Telescope {
  double fov_diameter;    // diameter of the FOV (in [rad])
  double focal_length;    // focal length of the telescope

  double time;
  struct vector r;              // position vector of the satellite
  struct vector v;              // velocity vector of the satellite

  // Unit vectors defining the telescope attitude:
  // nz is the irection of the telescope axis (in simplest case same direction as r),
  // nx usually points in the direction of motion, but is perpendicular to nz, and
  // ny is perpendicular to nx and nz
  struct vector nz, nx, ny;     
};

#endif /* TELESCOPE_H */

