#ifndef VECTOR_H
#define VECTOR_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** 3-dimensional vector. Data structure with three double-valued
    components. */
typedef struct {
  double x;
  double y;
  double z;
} Vector;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Creates a unit vector for specified right ascension and
    declination. Angles have to be given in [rad]. */
Vector unit_vector(const double ra, const double dec);

/** Returns a normalized vector of length 1.0 with the same direction
    as the original vector).*/
Vector normalize_vector(const Vector v);
/** Faster than normalize_vector. Deals with pointer instead of
    handling structures at function call. */
void normalize_vector_fast(Vector* const v);

/** Calculates the scalar product of two vectors. */
double scalar_product(const Vector* const v1, const Vector* const v2);

/** Calculates the vector product of two vectors. */
Vector vector_product(const Vector v1, const Vector v2);

/** Calculates the difference between two vectors. */
Vector vector_difference(const Vector x2, const Vector x1);

/** Function interpolates between two vectors at time t1 and t2 for
    the specified time and returns the interpolated vector. */
Vector interpolate_vec(const Vector v1, const double t1, 
		       const Vector v2, const double t2, 
		       const double time);

/** Interpolate between 2 vectors assuming that they describe a great
    circle on the unit sphere. The parameter 'phase' should have a
    value in the interval [0,1]. The return value is a normalized
    vector. */
Vector interpolateCircleVector(const Vector v1, 
			       const Vector v2, 
			       const double phase);

/** Function determines the equatorial coordinates right ascension and
    declination for a given vector pointing in a specific
    direction. The angles are calculated in [rad]. The given vector
    doesn't have to be normalized. */
void calculate_ra_dec(/** Direction. Does not have to be normalized. */
		      const Vector v, 
		      /** Right ascension. Units: [rad], Interval: [0;2*pi]. */ 
		      double* const ra,
		      /** Declination. Units: [rad], Interval: [-pi/2;pi/2]. */
		      double* const dec); 

/** Returns the value of the k-th dimension of a vector (k=0 ->
    x-value, k=1 -> y-value, k=2 -> z-value). */
double getVectorDimensionValue(const Vector* const vec, const int dimension);


#endif /* VECTOR_H */
