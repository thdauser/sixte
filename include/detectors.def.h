/** 
 * Header file for "measurement_det.c".
 * Contains definitions/function headers needed to handle detector stuff.
 */

#ifndef DETECTORS_DEF_H
#define DETECTORS_DEF_H (1)


#include "fits_pha.h"
#include "random.h"
#include "photon.h"
#include "point.h"

#include "detectors.types.h"
#include "astrosources.types.h"



// This routine is called after each photon event. It takes care of the detector 
// action. I.e. if the integration time is over, it reads out the detector array 
// and clears the pixels. It can also perform the readout of individual detector 
// lines, depending on the detector model (e.g. DEPFET for WFI in contrast to 
// framestore for eROSITA). If the exposure time is not exceeded it simply does 
// nothing.
void detector_action(Detector*, double time, struct Eventlist_File*, int* status);

inline void framestore_detector_action(void*, double time, 
					      struct Eventlist_File*, int *status);
inline void depfet_detector_action(void*, double time, 
					  struct Eventlist_File*, int *status); 
inline void tes_detector_action(void*, double time, struct Eventlist_File*, 
				       int *status); 
inline void htrs_detector_action(void*, double time, 
					struct Eventlist_File*, int *status);


// Returns a detector PHA channel for the given photon energy according to the RMF.
// Caution: This PHA channel doesn't have to be equivalent to the photon energy. 
// Depending on the detector redistribution matrix the energy can result in one 
// of several possible PHA channels with certain probability.
// If the energy is above the highest available energy bin in the RMF, the 
// return value is "-1".
long detector_rmf(float energy, RMF*);


// Get the PHA channel that corresponds to a particular charge.
long get_pha(float, Detector*);

// Get the charge that corresponds to a particular PHA channel according to 
// the ebounds table.
float get_charge(long, Ebounds*);


// Constructor: function allocates memory for the detector array.
Detector* get_Detector(int*);
int get_DetectorPixels(Detector*, int*);
// Destructor: function releases memory of detector.
//void free_Detector(Detector* detector); // TODO


// Get memory for detector EBOUNDS matrix and fill it with data from FITS file.
int get_ebounds(Ebounds*, int *Nchannels, const char filename[]);
// Release memory of detector EBOUNDS.
void free_ebounds(Ebounds *);


// Load the detector response matrix from the given RMF file.
int get_rmf(Detector *, char* rmf_name);


// Release memory of detector response matrix.
void free_rmf(RMF *);


// This function returns '1', if the specified detector pixel is active at 
// the time 'time'. If the pixel is, e.g., cleared at this time it cannot 
// receive a charge. In that case the function returns '0'.
int detector_active(int x, int y, Detector*, double time);
int htrs_detector_active(int x, int y, Detector*, double time);


// This function determines the integer pixel coordinates for a given 
// 2D floating point. The point lies within the hexagonal HTRS pixel.
int htrs_get_pixel(Detector*, struct Point2d, 
		   int* x, int* y, double* fraction);


// Constructor: this routine performs the initialization of the HTRS
// detector. The return value is the error status.
Detector* htrs_get_Detector(int *);
// Destructor: release all dynamically allocated memory in the HTRS
// detector structure.
void htrs_free_Detector(Detector*);


// Determines the pixel coordinates for a given point on a 2D array of 
// square pixels. The function returns all pixels that are affected due 
// to a splitting of the charge cloud. The coordinates of the affected 
// pixels are stored in the x[] and y[] arrays, the charge fraction in 
// each pixel in fraction[]. The number of totally affected 
// pixels is given by the function's return value. 
int get_pixel_square(Detector*, struct Point2d, 
		     int* x, int* y, double* fraction);




#endif  /*  DETECTORS_DEF_H */

