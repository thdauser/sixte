/** 
 * Header file for "measurement_det.c".
 * Contains definitions/function headers needed to handle detector stuff.
 */

#ifndef DETECTOR_H
#define DETECTOR_H 1


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

// The GNU Scientific Library Errorfunction is used to calculate charge 
// distribution of split events (assuming a Gaussian shape for the carge cloud).
#include <gsl/gsl_sf_erf.h>

// HEAdas header files
#include "headas.h"
#include "headas_error.h"
#include "fitsio.h"

#include "event_list.h"
#include "fits_pha.h"
#include "random.h"
#include "photon.h"
#include "global_constants.h"
#include "point.h"

#include "detector_types.h"
#include "sources_types.h"



// This routine is called after each photon event. It takes care of the detector 
// action. I.e. if the integration time is over, it reads out the detector array 
// and clears the pixels. It can also perform the readout of individual detector 
// lines, depending on the detector model (e.g. DEPFET for WFI in contrast to 
// framestore for eROSITA). If the exposure time is not exceeded it simply does 
// nothing.
void framestore_detector_action(struct Detector*, double time, 
				struct source_cat_entry background,
				struct Event_List_File*, int *status);


// This routine is implemented for the DEPFET detector on (WFI on IXO) with 
// the same purpose as the corresponding routine for the eROSITA framestore detector.
void depfet_detector_action(struct Detector*, double time, 
			    struct source_cat_entry background,
			    struct Event_List_File*, int *status); 


// This routine is implemented for the TES microcalorimeter (on IXO).
void tes_detector_action(struct Detector*, double time, 
			 struct source_cat_entry background,
			 struct Event_List_File*, int *status); 


// This routine is implemented for the HTRS (on IXO).
void htrs_detector_action(struct Detector*, double time, 
			  struct source_cat_entry background,
			  struct Event_List_File*, int *status); 


// Reads out the entire detector and creates event list entries for the 
// measured photons.
void readout(struct Detector, struct Event_List_File*, int *status);


// Reads out a particular detector line and creates event list entries 
// for the measured photons.
void readout_line(struct Detector, int line, struct Event_List_File*, 
		  int *fitsstatus);


// Returns a detector PHA channel for the given photon energy according to the RMF.
// Caution: This PHA channel doesn't have to be equivalent to the photon energy. 
// Depending on the detector redistribution matrix the energy can result in one 
// of several possible PHA channels with certain probability.
// If the energy is above the highest available energy bin in the RMF, the 
// return value is "-1".
long detector_rmf(float energy, struct RMF rmf);


// Get the PHA channel that corresponds to a particular charge.
long get_pha(float, struct Detector detector);

// Get the charge that corresponds to a particular PHA channel according to 
// the ebounds table.
float get_charge(long, struct Ebounds);


// Function allocates memory for the detector array.
int get_detector(struct Detector *detector);


// Get memory for detector EBOUNDS matrix and fill it with data from FITS file.
int get_ebounds(struct Ebounds *ebounds, int *Nchannels, const char filename[]);


// Release memory of detector EBOUNDS.
void free_ebounds(struct Ebounds ebounds);


// Load the detector response matrix from the given RMF file.
int get_rmf(struct Detector *detector, char *rmf_name);


// Release memory of detector response matrix.
void free_rmf(struct RMF);


// This routine clears the entire detector pixel array 
// (i.e., all created charges are removed, e.g., after read out).
void clear_detector(struct Detector detector);


// This routine clears a particular line of the detector pixel array 
// (i.e., all created charges are removed, e.g., after read out).
void clear_detector_line(struct Detector detector, int line);


/*
// Creates split events according to the size of the charge cloud.
// Assumption: +/- 3*ccsigma contain 99% of the charge cloud.
void split_events(double det_xa, double det_ya, int det_xi, int det_yi, 
		  struct Detector, double pha[3][3]);
*/


// Add background photons to the detector pixels  according to a given 
// background spectrum and count-rate.
void insert_background_photons(struct Detector detector, 
			       struct source_cat_entry background, 
			       double integration_time);


// This function returns '1', if the specified detector pixel is active at 
// the time 'time'. If the pixel is, e.g., cleared at this time it cannot 
// receive a charge. In that case the function returns '0'.
int detector_active(int x, int y, struct Detector detector, double time);


// This function calculates the value of the linear function f(x) = m*x + t
double linear_function(double x, double m, double t);


// Returns the lower line index of one particular line group with specified m
// and variable t for a given point in the 2D plane.
int htrs_get_line(struct Point2d point, double m, double dt, struct Detector detector);


// Get the line indices of the lines from each of the 3 groups of lines that define
// the hexagonal pixel shape.
void htrs_get_lines(struct Point2d, struct Detector, int* l0, int* l1, int* l2);


// This function determines the pixel index for a given 2D floating point.
// The point lies within the hexagonal HTRS pixel.
int htrs_get_pixel(struct Detector, struct Point2d, int*** pixel_relations);


// This routine performs the initialization of the HTRS detector.
// The return value is the error status.
int htrs_get_detector(struct Detector*);


// Determines the pixel coordinates for a given point on a 2D array of square pixels.
// The function returns all pixels that are affected due to a splitting of the charge
// cloud. The coordinates of the affected pixels are stored in the x[], and y[] arrays,
// the charge fraction in each pixel in fraction[]. The number of totally affected 
// pixels is given by the function's return value. 
int get_pixel_square(struct Detector, struct Point2d, 
		     int* x, int* y, double* fraction);




#endif

