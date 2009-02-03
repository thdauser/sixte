#ifndef DETECTOR_TYPES_H
#define DETECTOR_TYPES_H 1

#include "fitsio.h"
#include "sources_types.h"
#include "event_list_types.h"


#define HTRS_N_PIXELS 37  // total number of pixels in the HTRS array


// Define the different detector Types.
enum DetectorTypes {
  FRAMESTORE=1, 
  DEPFET    =2, 
  TES       =3,
  HTRS      =4
};


// Data type for the detector pixels.
struct Pixel {  // union
  // charge stored in the pixel
  float charge;    

  // list of photon arrival times in the pixel (needed for TES)
  double arrival; 
};


// Data structure for storing the detector EBOUNDS 
// (relation PHA channel -> [E_min; E_max]).
struct Ebounds_Row {
  long channel;
  float E_min, E_max;
};

struct Ebounds {
  struct Ebounds_Row *row;
};



// Data structure for storing one single line of the detector redistribution matrix.
struct RMF_Row {
  float E_low, E_high;
  int N_grp;
  int F_chan[1024];
  int N_chan[1024];
  float *matrix;
};

// Data structure to store the detector response matrix.
struct RMF {
  int Nrows;               // number of rows in the detector redistribution matrix
  int Ncols;               // number of columns in the detector redistribution matrix
  struct RMF_Row *row;     // data array for the detector RMF data
};



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
  double pixelwidth;       // width of a single pixel in the detector array [mu m]

  double integration_time; // Integration time of the entire framestore(!) CCD 
                           // detector array
                           // (= span of time between 2 subsequent readouts).
  long frame;              // Number of the current frame.

  double ccsigma;          // charge cloud sigma (needed to calculate size of 
                           // the charge cloud)
  double ccsize;           // size of the charge cloud [real pixels]
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
  double dead_time;        // Necessary time to read out one line of the DEPFET (!) 
                           // detector.
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



#endif
