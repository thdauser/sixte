#ifndef COMAEVENT_H
#define COMAEVENT_H 1


/** Coded Mask specific event. */
typedef struct {
  double time;
  float energy; /**< Event energy [keV] */
  int xi, yi; /**< Pixel coordinates starting at 0. */
} CoMaEvent;


#endif /* COMAEVENT_H */

