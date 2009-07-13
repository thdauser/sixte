#ifndef WFIDETECTOR_H
#define WFIDETECTOR_H 1

#include "sixt.h"
#include "genericdetector.h"
#include "squarepixels.h"


typedef struct {

  GenericDetector generic;
  SquarePixels pixels;

  /** Time required to clear a row of pixels on the detector. */
  double row_clear_time;

  /** Current readout line of the WFI detector. */
  int readout_line;

  /** Number of readout directions. Either 1 or 2. */
  int readout_directions;

} WFIDetector;


struct WFIDetectorParameters {
  struct GenericDetectorParameters generic;
  struct SquarePixelsParameters pixels;

  double t0;
};


#endif WFIDETECTOR
