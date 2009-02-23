#ifndef DETECTORS_H
#define DETECTOR_H (1)


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

// The GNU Scientific Library Errorfunction is used to calculate charge 
// distribution of split events (assuming a Gaussian shape for the carge cloud).
#include <gsl/gsl_sf_erf.h>

// HEAdas header files
#include "pil.h"
#include "headas.h"
#include "headas_error.h"
#include "fitsio.h"

#include "fits_pha.h"
#include "random.h"
#include "photon.h"
#include "point.h"

#include "sources_types.h"
#include "event_list_types.h"

#define INVALID_PIXEL (-1)   // flags an invalid pixel
#define HTRS_N_PIXELS (37)  // Total number of pixels in the HTRS array


int n_events, n_dead, n_interframe, n_outside;



// Define the different detector Types.
enum DetectorTypes {
  FRAMESTORE=1, 
  DEPFET    =2, 
  TES       =3,
  HTRS      =4
};


// Data type for the detector pixels.
struct {  // union
  // charge stored in the pixel
  float charge;    

  // list of photon arrival times in the pixel (needed for TES)
  double arrival; 
} Pixel;


// Data structure for storing the detector EBOUNDS 
// (relation PHA channel -> [E_min; E_max]).
struct {
  long channel;
  float E_min, E_max;
} Ebounds_Row;

struct Ebounds {
  struct Ebounds_Row *row;
};



// Data structure for storing one single line of the detector redistribution matrix.
struct {
  float E_low, E_high;
  int N_grp;
  int F_chan[1024];
  int N_chan[1024];
  float *matrix;
} RMF_Row;

// Data structure to store the detector response matrix.
struct {
  int Nrows;               // number of rows in the detector redistribution matrix
  int Ncols;               // number of columns in the detector redistribution matrix
  struct RMF_Row *row;     // data array for the detector RMF data
} RMF;



// Detector data structure.
struct Detector {
  enum DetectorTypes type; // Detector Type (framestore, depfet, ...) 

  // Detector array (contains the charge created by the x-ray 
  // photons or additional data)
  struct Pixel **restrict pixel;

  int width;               // width (and height) of the detector 
                           // (number of [integer pixels])
  int offset;              // offset of the detector array [integer pixels], 
                           // the physical origin of the detector (at the center) 
                           // has the array-index 'offset'
  double pixelwidth;       // width of a single pixel in the detector array [m]

  double integration_time; // Integration time of the entire pnCCD (!) 
                           // detector array
                           // (= span of time between 2 subsequent readouts).
  double dead_time;        // Necessary time to read out one line of the DEPFET 
                           // or one pixel of the HTRS (!) detectors.
  long frame;              // Number of the current frame.

  double ccsigma;          // charge cloud sigma (needed to calculate size of 
                           // the charge cloud) [m]
  double ccsize;           // size of the charge cloud [m]
  long low_threshold;      // lower detector threshold (in PHA)
  
  int Nchannels;           // Number of detector PHA channels
  struct Ebounds ebounds;  // Detector energy bounds (relation PHA channel -> 
                           // [E_min; E_max])
  struct RMF rmf;          // RMF


  // This is a pointer to the routine, which is called after each photon event.
  // Its task is to manage the detector action, i.e., perform the readout process
  // if it necessary.
  void (*detector_action) (struct Detector*, double, 
			   struct source_cat_entry,
			   struct Event_List_File*, int *);


  // DEPFET specific parameters:

  double readout_time;     // Current readout time (i.e., the end of the 
                           // integration time/beginning of dead time).
  double clear_time;       // Time required to clear a row of pixels on the detector.
  int readout_line;        // Current readout line of the DEPFET detector 
                           // (not used for framestore).



  // HTRS specific parameters:

  double a;   // length of one egde of the hexagonal pixels
  double h;   // height of one of the six equilateral triangles that define
              // the hexagonal structure of the pixels.

  // -- 2 different access numbering schemes --
  // 2d array contains linear numbers
  int** htrs_icoordinates2pixel;   
  // linear array contains 2d integer coordinates
  struct Point2i* htrs_pixel2icoordinates; 

  // Data structure to obtain a pixel from given coordinates
  int*** htrs_lines2pixel; 

};



/////////////////////////////////////////////////////////////////////////////////



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
// the same purpose as the corresponding routine for the eROSITA framestore 
// detector.
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


// Determines the index of the minimum value in an array of distances to 
// the pixel borders.
int min_dist(double array[], int directions);


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


// Add background photons to the detector pixels  according to a given 
// background spectrum and count-rate.
void insert_background_photons(struct Detector detector, 
			       struct source_cat_entry background, 
			       double integration_time);


// This function returns '1', if the specified detector pixel is active at 
// the time 'time'. If the pixel is, e.g., cleared at this time it cannot 
// receive a charge. In that case the function returns '0'.
int detector_active(int x, int y, struct Detector detector, double time);
int htrs_detector_active(int x, int y, struct Detector detector, double time);


// This function calculates the value of the linear function f(x) = m*x + t
double linear_function(double x, double m, double t);


// Returns the lower line index of one particular line group with specified m
// and variable t for a given point in the 2D plane.
int htrs_get_line(struct Point2d point, double m, double dt, struct Detector);


// Get the line indices of the lines from each of the 3 groups of lines that define
// the hexagonal pixel shape.
void htrs_get_lines(struct Point2d, struct Detector, int* l);


// This function determines the integer pixel coordinates for a given 
// 2D floating point. The point lies within the hexagonal HTRS pixel.
int htrs_get_pixel(struct Detector, struct Point2d, 
		   int* x, int* y, double* fraction);


// Returns the pixel index that corresponds to the pixel segment which is
// defined by the 3 given line indices.
int htrs_get_lines2pixel(int* l, struct Detector detector);


// Returns the (integer) pixel coordinates of the pixel which is specified 
// by its number.
struct Point2i htrs_get_pixel2icoordinates(int pixel, struct Detector detector);


// This routine performs the initialization of the HTRS detector.
// The return value is the error status.
int htrs_get_detector(struct Detector*);


// Release all dynamically allocated memory in the HTRS detector structure.
void htrs_free_detector(struct Detector* detector);


// Determines the pixel coordinates for a given point on a 2D array of 
// square pixels. The function returns all pixels that are affected due 
// to a splitting of the charge cloud. The coordinates of the affected 
// pixels are stored in the x[] and y[] arrays, the charge fraction in 
// each pixel in fraction[]. The number of totally affected 
// pixels is given by the function's return value. 
int get_pixel_square(struct Detector, struct Point2d, 
		     int* x, int* y, double* fraction);




#endif  /* DETECTORS_H */

