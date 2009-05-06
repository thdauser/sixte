#ifndef TELESCOPE_H
#define TELESCOPE_H 1

#include "vector.h"

struct Telescope {
  double fov_diameter;    // diameter of the FOV (in [rad])
  double focal_length;    // focal length of the telescope

  double time;
  struct vector r;              // position vector of the satellite
  struct vector v;              // velocity vector of the satellite

  /* Unit vectors defining the telescope attitude */
  /** the direction of the telescope axis (in simplest case same direction as r) */
  struct vector nz;
  /** usually points in the direction of motion, but is perpendicular to nz */
  struct vector nx;
  /** perpendicular to nx and nz */
  struct vector ny; 
};

#endif /* TELESCOPE_H */

