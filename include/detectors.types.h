#ifndef DETECTOR_TYPES_H
#define DETECTOR_TYPES_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif

#include "impactlist.h"
#include "eventlist.types.h"
#include "detectors.enum.h"


/** Parameters for initialization of the detector data structure. */
struct DetectorParameters {
  int width; /**< Width of the detector in pixels. */
  double pixelwidth;
  double ccsigma;

  char rmf_filename[MAXMSG];
  long pha_threshold;
  float energy_threshold;

  double t0;
};


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

  /** Detector pixel array. Contains the charge created by the x-ray 
   * photons, and additional data if necessary. */
  struct Pixel ** pixel;

  int width; /**< Width (and height) of the detector (number of [integer pixels]). */

  int offset; /**< Offset of the detector array [integer pixels]. 
	       * The physical origin of the detector (at the center of the detector) 
	       * has the array-index 'offset'. */

  double pixelwidth; /**< Width of a single pixel in the detector array [m]. */

  double ccsigma; /**< Charge cloud sigma [m]. This quantity is used to calculate size of 
		   * the charge cloud. */
  double ccsize; /**< Size of the charge cloud [m]. Defined as three times ccsigma. */

  double dead_time; // Necessary time to read out one line of the DEPFET 
                    // or one pixel of the HTRS (!) detectors.

  long frame; /** Number of the current frame. */

  double readout_time; /**< Current readout time. The end of the integration 
			* time / beginning of dead time. */
  
  long pha_threshold; // lower detector PHA threshold [PHA channels]
  float energy_threshold; /**< Lower detector energy threshold [keV]. */
  // If the PHA threshold is -1, the energy threshold is used.
  
  /** Detector response matrix. Includes the RMF and the detector-specific
   * response elements like filter transmission and quantum efficiency.
   * So the sum of each line of the response matrix HAS TO BE less
   * or equal 1 (sum <= 1) !
   * The RMF can be normalized explicitly to be a real RMF without 
   * photon loss due response effects by setting the compiler flag 
   * "-DNORMALIZE_RMF".
   * The mirror specific elements are treated SEPARATELY in the photon
   * imaging process. */
  struct RMF* rmf;


  /** Detector-specific elements. Pointer to a data structure that contains 
   * detector-specific properties for the different detector type. 
   * The data structure is allocated and assigned to the pointer by the
   * init_{DETECTOR_TYPE}(...) routine.*/
  void* specific;


  // DEPFET specific parameters:
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


  /** Pointer to the routine that is called to read out the detector.
   * The routine has to decide on its own, whether a readout is required, as e.g.
   * the integration time is expired.
   * If NO readout is to be performed, the routine simply does nothing.
   * Otherwise it calls the appropriate routines to store the events from the detector
   * in the output event list. 
   * The routine is also responsible for clearing the detector pixels after the readout. */
  void (*readout) (void*, double, struct Eventlist_File*, int*);

  /** Pointer to the routine that is called when a new photon hits the detector.
   * The routine determines the generated charge according to the photon energy and
   * the detector response. Split events are taken into account if 'ccsize>0.'. */
  void (*add_impact) (void*, Impact*);

} Detector;



#endif /* DETECTORS_TYPES_H */

