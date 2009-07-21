#ifndef WFIDETECTOR_H
#define WFIDETECTOR_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "squarepixels.h"
#include "eventfile.h"
#include "impactlist.h"
#include "wfieventfile.h"


typedef struct {

  GenericDetector generic;
  SquarePixels pixels;

  /** Time required to read out one line of the WIF detector. */
  double line_readout_time;

  /** Time required to clear a line of pixels on the detector. 
   * Clear time is a part of the readout time. */
  double line_clear_time;

  /** Number of readout directions. Either 1 or 2. */
  int readout_directions;

  /** Current readout time. */
  double readout_time;

  /** Current readout line of the WFI detector. */
  int readout_lines[2];

  /** Number of the current frame. */
  long frame; 

  /** Output event list. */
  WFIEventFile eventlist;

} WFIDetector;


struct WFIDetectorParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double line_readout_time, line_clear_time;
  int readout_directions;
  char* eventlist_filename;
  char* eventlist_template;

  double t0;
};


////////////////////////////////////////////////////////


/** Set up the configuration of a WFIDetector. 
 * The routine also calls the init routines of the underlying data structures. */
int initWFIDetector(WFIDetector*, struct WFIDetectorParameters*);

/** Clean up the WFIDetector data structure. 
 * Release allocated memory and call clean-up routines of underlying structures. */
int cleanupWFIDetector(WFIDetector* wd);

/** Check out whether the WFI detector needs to be read out. */
int checkReadoutWFIDetector(WFIDetector*, double time);

/**  Read out the charge from the currently active readout lines (one or two) 
 * of the WFI detector. */
inline int readoutLinesWFIDetector(WFIDetector*);

/** Add a photon impact to the WFI detector pixel array.
 * The generated charge is determined according to the detector RSP.
 * Split events are taken into account based on a Gaussian-shape charge cloud model.
 * The charge in the in the pixels is summed up, i.e., realistic pile up is implemented. */
int addImpact2WFIDetector(WFIDetector*, Impact*);

/** Returns 1 if the detector is sensitive in the pixel (x,y) at the specified time. */
inline int WFIDetectorIsSensitive(int y, WFIDetector*, double time);


#endif /* WFIDETECTOR */
