#ifndef VECTOR_H
#define VECTOR_H 1

#include <stdlib.h>
#include <math.h>


// data structure with three double-valued components (= 3d vector)
typedef struct {
  double x;
  double y;
  double z;
} Vector;


// Returns a unit vector pointing in the specified direction of 
// right ascension and declination
Vector unit_vector(const double rasc, const double dec);

// returns a normalized vector (length 1.0, same direction)
Vector normalize_vector(Vector);

// calculates the scalar product of two vectors
inline double scalar_product(Vector* const, Vector* const);

// calculates the vector product of two vectors
Vector vector_product(Vector, Vector);

/** Calculates the difference between two vectors. */
Vector vector_difference(Vector x2, Vector x1);

/** Function interpolates between two vectors at time t1 and t2 for the specified time 
 * and returns the interpolated vector. */
Vector interpolate_vec(Vector v1, double t1, Vector v2, 
		       double t2, double time);
Vector interpolate_vec2(Vector v1, double t1, Vector v2, double t2, double time);

/** Function determines the equatorial coordinates of right ascension 
 * and declination for a given vector pointing in a specific direction. 
 * The angles are calculated in [rad].
 * The given vector doesn't have to be normalized. */
void calculate_ra_dec(Vector v, /**< Direction. Does not have to be normalized. */
		      double* ra, /**< Right ascension. Unit: [rad], Interval: [-pi;pi]. */ 
		      double* dec /**< Declination. Unit: [rad], Interval: [-pi/2;pi/2]. */);


#endif /* VECTOR_H */
