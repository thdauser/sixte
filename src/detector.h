#ifndef DETECTOR_H
#define DETECTOR_H 1

#include "sixt.h"
#include "random_sixt.h"
#include "fits_pha.h"

#define GAUSSIAN_CHARGE_CLOUDS 1 // Assume a Gaussian charge cloud shape.

// Reads out the entire detector and creates event list entries for the 
// measured photons.
static inline void readout(Detector*, struct Eventlist_File*, int *status);


// Reads out a particular detector line and creates event list entries 
// for the measured photons.
static inline void readout_line(Detector*, int line, struct Eventlist_File*, 
				int *fitsstatus);


// Determines the index of the minimum value in an array of distances to 
// the pixel borders.
static inline int min_dist(double array[], int directions);


// This routine clears the entire detector pixel array 
// (i.e., all created charges are removed, e.g., after read out).
static inline void clear_detector(Detector*);


// This routine clears a particular line of the detector pixel array 
// (i.e., all created charges are removed, e.g., after read out).
static inline void clear_detector_line(Detector*, int line);


// Add background photons to the detector pixels  according to a given 
// background spectrum and count-rate.
static void insert_background_photons(Detector*, PointSource background, 
				      double integration_time);


// Returns the pixel index that corresponds to the pixel segment which is
// defined by the 3 given line indices.
static inline int htrs_get_lines2pixel(int* l, Detector*);


// Returns the (integer) pixel coordinates of the pixel which is specified 
// by its number.
static inline struct Point2i htrs_get_pixel2icoordinates(int pixel, Detector*);


// This function calculates the value of the linear function f(x) = m*x + t
static inline double linear_function(double x, double m, double t);


// Returns the lower line index of one particular line group with specified m
// and variable t for a given point in the 2D plane.
static inline int htrs_get_line(struct Point2d point, double m, double dt, Detector*);


// Get the line indices of the lines from each of the 3 groups of lines that define
// the hexagonal pixel shape.
static inline void htrs_get_lines(struct Point2d, Detector*, int* l);



#endif /* DETECTOR_H */

