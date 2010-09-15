#include "vector.h"


Vector unit_vector(const double ra, const double dec)
{
  Vector x;
  double cos_dec = cos(dec);

  x.x = cos_dec * cos(ra);
  x.y = cos_dec * sin(ra);
  x.z = sin(dec);

  return(x);
}



Vector normalize_vector(Vector x) {
  double l;  // length of the vector x
  Vector y;  // normalized vector

  l = sqrt(pow(x.x,2.0)+pow(x.y,2.0)+pow(x.z,2.0));

  y.x = x.x/l;
  y.y = x.y/l;
  y.z = x.z/l;

  return(y);
}
/////////////////////////////////////////////////////
void normalize_vector_fast(Vector* v) {
  double l;  // length of the vector x
  l = sqrt(pow(v->x,2.0)+pow(v->y,2.0)+pow(v->z,2.0));

  v->x /= l;
  v->y /= l;
  v->z /= l;
}



/** Calculates the scalar product of two vector structures.*/
double scalar_product(const Vector* const x, const Vector* const y)
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



Vector vector_difference(Vector x2, Vector x1) {
  Vector z;  // return vector

  z.x = x2.x-x1.x;
  z.y = x2.y-x1.y;
  z.z = x2.z-x1.z;

  return(z);
}

 

Vector interpolate_vec(Vector v1, double t1, 
		       Vector v2, double t2, 
		       double time) {
  Vector pos;
  
  pos.x = v1.x + (time-t1)/(t2-t1)*(v2.x-v1.x);
  pos.y = v1.y + (time-t1)/(t2-t1)*(v2.y-v1.y);
  pos.z = v1.z + (time-t1)/(t2-t1)*(v2.z-v1.z);

  return(pos);
}



Vector interpolateCircleVector(Vector v1, Vector v2, double phase)
{
  Vector x1 = normalize_vector(v1); // Use as first base vector.
  Vector x2 = normalize_vector(v2);
  Vector r; // Return value. 

  // Calculate cosine of angle between x1 and x2 (v1 and v2) [rad].
  double cosine_value = scalar_product(&x1, &x2);

  if (fabs(cosine_value) < cos(0.1/3600*M_PI/180.)) { 
    // The misalignment between the 2 vectors is more than 1 arcsec.
    // This is important to check for the subsequent algorithm, 
    // because the vectors should not be aligned parallel or 
    // anti-parallel.

    // Angle between x1 and x2:
    double phi = acos(cosine_value); 

    // Calculate the second base vector spanning the plane of 
    // the circle.
    Vector d = { .x=x2.x-cos(phi)*x1.x, 
		 .y=x2.y-cos(phi)*x1.y, 
		 .z=x2.z-cos(phi)*x1.z };
    x2 = normalize_vector(d); 
    
    // Determine the angle corresponding to the phase.
    r.x = cos(phase*phi)*x1.x + sin(phase*phi)*x2.x;
    r.y = cos(phase*phi)*x1.y + sin(phase*phi)*x2.y;
    r.z = cos(phase*phi)*x1.z + sin(phase*phi)*x2.z;

  } else { 
    // There is quasi no motion at all, so perform a linear
    // interpolation.
    r.x = x1.x + phase*(x2.x-x1.x);
    r.y = x1.y + phase*(x2.y-x1.y);
    r.z = x1.z + phase*(x2.z-x1.z);
    r = normalize_vector(r);
  }

  return(r);
}



void calculate_ra_dec(Vector v, double* ra, double* dec)
{
  // Determine the declination:
  *dec = asin(v.z/sqrt(scalar_product(&v, &v)));

  // Determine the right ascension:
  *ra = atan2(v.y, v.x);
}



double getVectorDimensionValue(Vector* vec, int dimension) 
{
  if (0==dimension) {
    return(vec->x);
  } else if (1==dimension) {
    return(vec->y);
  } else {
    return(vec->z);
  }
}

