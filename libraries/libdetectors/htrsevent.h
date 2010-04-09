#ifndef HTRSEVENT_H
#define HTRSEVENT_H 1


/** HTRS-specific event. */
typedef struct {
  /** Detection time. */
  double time; 

  /** PHA channel. */
  long pha; 

  /** Photon energy [keV]. */
  float energy;

  /** Number of the pixel, where the event is detected. Numbering
      starts at 0. */
  int pixel; 

  /** Event grade. This gives information about the energy and time
      resolution of the event. There are the following event grades:

      - 0: nominal time and energy resolution 

      - 1: nominal time but no energy resolution 

      - 2: event cannot be distinguished from the previous event
  */
  int grade;

  /** Exact impact position of the photon on the detector. */
  double x, y; 

} HTRSEvent;


#endif /* HTRSEVENT_H */

