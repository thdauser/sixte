#ifndef COMAEVENT_H
#define COMAEVENT_H 1


/** Coded Mask specific event. */
typedef struct {
  double time;
  /** Pixel charge represented by the photon energy [keV]. */
  float charge; 
  /** Pixel coordinates starting at 0. */
  int rawx, rawy;
} CoMaEvent;


#endif /* COMAEVENT_H */

