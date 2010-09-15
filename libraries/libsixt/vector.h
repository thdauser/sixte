#ifndef VECTOR_H
#define VECTOR_H 1

#include "sixt.h"

/** 3-dimensional vector. Data structure with three double-valued
    components. */
typedef struct {
  double x;
  double y;
  double z;
} Vector;


/** Creates a unit vector for specified right ascension and
    declination. Angles have to be given in [rad]. */
Vector unit_vector(const double ra, const double dec);

/** Returns a normalized vector of length 1.0 with the same direction
    as the original vector).*/
Vector normalize_vector(Vector);
/** Faster than normalize_vector. Deals with pointer instead of
    handling structures at function call. */
void normalize_vector_fast(Vector* v);

/** Calculates the scalar product of two vectors. */
double scalar_product(const Vector* const v1, const Vector* const v2);

/** Calculates the vector product of two vectors. */
Vector vector_product(Vector, Vector);

/** Calculates the difference between two vectors. */
Vector vector_difference(Vector x2, Vector x1);

/** Function interpolates between two vectors at time t1 and t2 for
    the specified time and returns the interpolated vector. */
Vector interpolate_vec(Vector v1, double t1, Vector v2, 
		       double t2, double time);

/** Interpolate between 2 vectors assuming that they discribe a great
    circle on the unit sphere. The parameter 'phase' should have a
    value in the interval [0,1]. The return value is a normalized
    vector. */
Vector interpolateCircleVector(Vector v1, Vector v2, double phase);

/** Function determines the equatorial coordinates of right ascension
    and declination for a given vector pointing in a specific
    direction. The angles are calculated in [rad]. The given vector
    doesn't have to be normalized. */
void calculate_ra_dec(Vector v, /**< Direction. Does not have to be normalized. */
		      double* ra, /**< Right ascension. Unit: [rad], Interval: [-pi;pi]. */ 
		      double* dec /**< Declination. Unit: [rad], Interval: [-pi/2;pi/2]. */);

/** Returns the value of the k-th dimension of a vector (k=0 ->
    x-value, k=1 -> y-value, k=2 -> z-value). */
double getVectorDimensionValue(Vector* vec, int dimension);


#endif /* VECTOR_H */
