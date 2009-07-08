#include "vector.h"


/////////////////////////////////////////////////////
Vector unit_vector(const double ra, const double dec)
{
  Vector x;
  double cos_dec = cos(dec);

  x.x = cos_dec * cos(ra);
  x.y = cos_dec * sin(ra);
  x.z = sin(dec);

  return(x);
}



/////////////////////////////////////////////////////
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
inline void normalize_vector_fast(Vector* v) {
  double l;  // length of the vector x
  l = sqrt(pow(v->x,2.0)+pow(v->y,2.0)+pow(v->z,2.0));

  v->x /= l;
  v->y /= l;
  v->z /= l;
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
Vector interpolate_vec2(Vector v1, double t1, Vector v2, double t2, double time)
{
  Vector r; // r(time)
  Vector x1 = normalize_vector(v1); // use as first base vector
  Vector x2 = normalize_vector(v2);

  double alpha = acos(scalar_product(&x1, &x2)); // angle between v1 and v2 [rad]

  if (alpha > 0.1/3600*M_PI/180.) { // Greater than 0.1 arcsec.
    Vector d = {.x=v2.x-cos(alpha)*v1.x, 
		.y=v2.y-cos(alpha)*v1.y, 
		.z=v2.z-cos(alpha)*v1.z};
    x2 = normalize_vector(d); // second base vector
    
    double dalpha = alpha*(time-t1)/(t2-t1); // angle-difference

    r.x = cos(dalpha)*x1.x + sin(dalpha)*x2.x;
    r.y = cos(dalpha)*x1.y + sin(dalpha)*x2.y;
    r.z = cos(dalpha)*x1.z + sin(dalpha)*x2.z;			      
  } else { // Quasi no motion at all.
    r = v1;
  }

  return(r);
}


/////////////////////////////////////////////////////////////////
void calculate_ra_dec(Vector v, double* ra, double* dec)
{
  // Determine the declination:
  *dec = asin(v.z/sqrt(scalar_product(&v, &v)));

  // Determine the right ascension:
  *ra = atan2(v.y, v.x);
}

