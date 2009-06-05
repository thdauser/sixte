#ifndef VECTOR_H
#define VECTOR_H 1

#include <stdlib.h>
#include <math.h>


// data structure with three double-valued components (= 3d vector)
struct vector {
  double x;
  double y;
  double z;
};


// Returns a unit vector pointing in the specified direction of 
// right ascension and declination
struct vector unit_vector(const double rasc, const double dec);

// returns a normalized vector (length 1.0, same direction)
struct vector normalize_vector(struct vector);

// calculates the scalar product of two vectors
inline double scalar_product(struct vector* const, struct vector* const);

// calculates the vector product of two vectors
struct vector vector_product(struct vector, struct vector);

/** Calculates the difference between two vectors. */
struct vector vector_difference(struct vector x2, struct vector x1);

/** Function interpolates between two vectors at time t1 and t2 for the specified time 
 * and returns the interpolated vector. */
struct vector interpolate_vec(struct vector v1, double t1, struct vector v2, 
			      double t2, double time);

/** Function determines the equatorial coordinates of right ascension 
 * and declination for a given vector pointing in a specific direction. 
 * The angles are calculated in [radians].
 * The given vector doesn't have to be normalized. */
void calculate_ra_dec(struct vector v, /**< Direction. Does not have to be normalized. */
		      double* ra, /**< Right ascension. Unit: [rad], Interval: [-pi;pi]. */ 
		      double* dec /**< Declination. Unit: [rad], Interval: [-pi/2;pi/2]. */);


#endif /* VECTOR_H */
