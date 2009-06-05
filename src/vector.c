#include "vector.h"


/** Creates a unit vector for specified right ascension and declination.
 * Angles have to be given in [rad]. */
Vector unit_vector(const double rasc, const double dec)
{
  Vector x;
  double cos_dec = cos(dec);

  x.x = cos_dec * cos(rasc);
  x.y = cos_dec * sin(rasc);
  x.z = sin(dec);

  return(x);
}



/** Returns a normalized vector (length 1.0, same direction).*/
Vector normalize_vector(Vector x) {
  double l;         // length of the vector x
  Vector y;  // normalized vector

  l = sqrt(pow(x.x,2.0)+pow(x.y,2.0)+pow(x.z,2.0));

  y.x = x.x/l;
  y.y = x.y/l;
  y.z = x.z/l;

  return(y);
}



/** Calculates the scalar product of two vector structures.*/
inline double scalar_product(Vector* const x, Vector* const y)
{
  return(x->x * y->x + x->y * y->y + x->z * y->z);
}



/** Calculates the vector product of two vectors. */
Vector vector_product(Vector x, Vector y) {
  Vector z;  // return vector

  z.x = x.y*y.z-x.z*y.y;
  z.y = x.z*y.x-x.x*y.z;
  z.z = x.x*y.y-x.y*y.x;

  return(z);
}


////////////////////////////////////////////////////////////////
Vector vector_difference(Vector x2, Vector x1) {
  Vector z;  // return vector

  z.x = x2.x-x1.x;
  z.y = x2.y-x1.y;
  z.z = x2.z-x1.z;

  return(z);
}

 

/////////////////////////////////////////////////////////////////
Vector interpolate_vec(Vector v1, double t1, 
			      Vector v2, double t2, 
			      double time) {
  Vector pos;
  
  pos.x = v1.x + (time-t1)/(t2-t1)*(v2.x-v1.x);
  pos.y = v1.y + (time-t1)/(t2-t1)*(v2.y-v1.y);
  pos.z = v1.z + (time-t1)/(t2-t1)*(v2.z-v1.z);

  return(pos);
}



/////////////////////////////////////////////////////////////////
void calculate_ra_dec(Vector v, double* ra, double* dec)
{
  // Determine the declination:
  *dec = asin(v.z/sqrt(scalar_product(&v, &v)));

  // Determine the right ascension:
  *ra = atan2(v.y, v.x);
}

