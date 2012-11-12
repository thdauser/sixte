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

      -  0: event is measured with the slow shaper

      -  1: event is measured with the fast shaper

      -  2: event is mixed up with another event (pile-up) because they
            cannot be distinguished even with the fast shaper

      - -1: lost during reset

  */
  int grade;
  

  /** Exact impact position of the photon on the detector. */
  double x, y; 

} HTRSEvent;


#endif /* HTRSEVENT_H */

