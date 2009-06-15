#include "create_rnd_source.h"


void create_rnd_source_position(double *rasc, double *dec) {
  double r2d;
  double x,y,z;

  do {
    // create random position inside a cube [-1;1]^3
    x = (get_random_number()*2)-1.0;
    y = (get_random_number()*2)-1.0;
    z = (get_random_number()*2)-1.0;

    // check, if the point is inside the unit sphere
    if ((x*x + y*y + z*z)<=1.0){
      // calculate right ascension and declination
      r2d = sqrt(x*x + y*y);
      if (y < 0.0) {
        *rasc = 2.0*M_PI - acos(x/r2d);
      } else {
        *rasc = acos(x/r2d);
      }
      *dec = atan2(z,r2d);
      break;
    }
  } while (1);
}





// creates a random source intensity according to a powerlaw distribution with index alpha
float get_rnd_intensity(double powerlaw_index) {
  return((float)pow(get_random_number(), 1./powerlaw_index));
}

