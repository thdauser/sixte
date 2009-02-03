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
struct vector unit_vector(double rasc, double dec);

// returns a normalized vector (length 1.0, same direction)
struct vector normalize_vector(struct vector);

// calculates the scalar product of two vectors
double scalar_product(struct vector, struct vector);

// calculates the vector product of two vectors
struct vector vector_product(struct vector, struct vector);

// interpolates between to orbit positions for the specified time
// and returns a vector
struct vector interpolate_vec(struct vector v1, double t1, struct vector v2, 
			      double t2, double time);


#endif
