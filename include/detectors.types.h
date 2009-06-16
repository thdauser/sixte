#ifndef DETECTOR_TYPES_H
#define DETECTOR_TYPES_H 1

#include "fitsio.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


#include "eventlist.types.h"
#include "detectors.enum.h"



/** Represents a detector pixel. */
struct Pixel {  // union
  /** Charge stored in the pixel. */
  float charge;    

  /** List of photon arrival times in the pixel (needed for TES). */
  double arrival; 
};



/** Detector data structure. */
typedef struct {
  DetectorTypes type;      /**< Detector Type (framestore, depfet, ...). */

  // Detector array (contains the charge created by the x-ray 
  // photons or additional data)
  struct Pixel ** pixel;

  int width;               /**< Width (and height) of the detector 
			    *(number of [integer pixels]). */
  int offset;              /**< Offset of the detector array [integer pixels], 
			    * the physical origin of the detector (at the center)
			    * has the array-index 'offset'. */
  double pixelwidth;       /**< Width of a single pixel in the detector array [m]. */

  double integration_time; // Integration time of the entire pnCCD (!) 
                           // detector array
                           // (= span of time between 2 subsequent readouts).
  double dead_time;        // Necessary time to read out one line of the DEPFET 
                           // or one pixel of the HTRS (!) detectors.
  long frame;              // Number of the current frame.

  double ccsigma;          // charge cloud sigma (needed to calculate size of 
                           // the charge cloud) [m]
  double ccsize;           // size of the charge cloud [m]
  
  long pha_threshold;      // lower detector PHA threshold [PHA channels]
  float energy_threshold;  // lower detector energy threshold [kev]
  // If the PHA threshold is 0, the energy threshold is used.
  
  //  int Nchannels;           // Number of detector PHA channels
  //  Ebounds ebounds;         // Detector energy bounds (relation PHA channel -> 
                               // [E_min; E_max])
  struct RMF* rmf;


  // This is a pointer to the routine, which is called after each photon event.
  // Its task is to manage the detector action, i.e., perform the readout process
  // if it necessary.
  void (*action) (void*, double, struct Eventlist_File*, int *);


  // DEPFET specific parameters:

  double readout_time;     // Current readout time (i.e., the end of the 
                           // integration time/beginning of dead time).
  double clear_time;       // Time required to clear a row of pixels on the detector.
  int readout_line;        // Current readout line of the DEPFET detector 
                           // (not used for framestore).
  int readout_directions;  // Either 1 or 2.


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

} Detector;



#endif /* DETECTORS_TYPES_H */

