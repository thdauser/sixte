#ifndef HTRSEVENT_H
#define HTRSEVENT_H 1


/** HTRS-specific event. */
typedef struct {
  /** Detection time. */
  double time; 

  /** PHA channel. */
  long pha; 

  /** Number of the pixel, where the event is detected. Numbering
      starts at 0. */
  int pixel; 

  /** Exact impact position of the photon on the detector. */
  double x, y; 

} HTRSEvent;


#endif /* HTRSEVENT_H */

