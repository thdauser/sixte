#include <stdio.h>
#include <stdlib.h>

#include "check_fov.h"
#include "vector.h"

int main()
{
  // x is the vector, where we check, wheter it is inside the FOV
  struct vector x = unit_vector(0.0, 3.1415);

  // constants
  const struct vector x0 = {1.0, 0.0, 0.0};   // direction of the telescope axis
  const struct vector h1 = {0.0, 0.0, 1.0};   // normal vector 1
  const struct vector h2 = {0.0, 1.0, 0.0};   // normal vector 2
  const double dec_max = 0.5;                 // half-height of the FOV
  const double rasc_max = 0.6;                // half-width of the FOV
  const double sin_dec_max = sin(dec_max);
  const double sin_rasc_max = sin(rasc_max);
  const double min_align = cos(sqrt(dec_max*dec_max + rasc_max*rasc_max));

  switch (check_fov(x, x0, /* h1, h2, sin_dec_max, sin_rasc_max, */ min_align)) {
    case 1:
      printf("Source outside FOV (declination out of range)\n");
      break;
    case 2:
      printf("Source outside FOV (right ascension out of range)\n");
      break;
    case 3:
      printf("Source outside FOV (declination and right ascension out of range)\n");
      break;
    default:
      printf("Source inside FOV\n");
      break;
  }

//  printf("%lf\t%lf\t%lf\t%lf\n", x.x, x.y, x.z, scalar_product(x,x));

  return EXIT_SUCCESS;
}
