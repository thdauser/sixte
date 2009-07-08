#ifndef DETECTOR_TYPES_H
#define DETECTOR_TYPES_H 1

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


#include "eventlist.types.h"
#include "detectors.enum.h"


/** Represents a detector pixel. */
struct Pixel {
  /** Charge stored in the pixel. */
  float charge;    
  /** List of photon arrival times in the pixel (needed for TES). */
  double arrival; 
};



/** Detector data structure. */
typedef struct {
  DetectorTypes type; /**< Detector Type (FRAMESTORE, WFI, TES, HTRS, ...). */

  // Detector array (contains the charge created by the x-ray 
  // photons or additional data)
  struct Pixel ** pixel;

  int width; /**< Width (and height) of the detector (number of [integer pixels]). */

  int offset; /**< Offset of the detector array [integer pixels]. 
	       * The physical origin of the detector (at the center of the detector) 
	       * has the array-index 'offset'. */

  double pixelwidth; /**< Width of a single pixel in the detector array [m]. */


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
  float energy_threshold; /**< Lower detector energy threshold [keV]. */
  // If the PHA threshold is -1, the energy threshold is used.
  
  struct RMF* rmf; /**< Detector response matrix. Includes the RMF and the detector-specific
		    * response elements like filter transmission and quantum efficiency.
		    * So the sum of each line of the response matrix HAS TO BE less
		    * or equal 1 (sum <= 1) !
		    * The RMF can be normalized explicitly to be a real RMF without 
		    * photon loss due response effects by setting the compiler flag 
		    * "-DNORMALIZE_RMF".
		    * The mirror specific elements are treated SEPARATELY in the photon
		    * imaging process. */


  // This is a pointer to the routine, which is called after each photon event.
  // Its task is to manage the detector action, i.e., perform the readout process
  // if it necessary.
  void (*readout) (void*, double, struct Eventlist_File*, int *);


  /** Detector-specific elements. Pointer to a data structure that contains 
   * detector-specific properties for the different detector type. 
   * The data structure is allocated and assigned to the pointer by the
   * init_{DETECTOR_TYPE}(...) routine.*/
  void* specific;


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

