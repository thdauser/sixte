#include <stdio.h>
#include <math.h>
#include "random.h"


// creates a source on the surface of the unit sphere
//    output --> right ascension (radians), declination (radians), intensity (I >= 1.0)
void create_rnd_source_position(double *rasc, double *dec);

// creates a random source intensity according to a powerlaw distribution
float get_rnd_intensity(double powerlaw_index);
