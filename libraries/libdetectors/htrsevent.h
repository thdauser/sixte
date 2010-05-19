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

      - 2: event is measured, but there is at least one subsequent
	   event that cannot be distinguished from this event.

      - 3: event cannot be distinguished from the previous event
      
      - 4: event is lost during a pixel reset.

  */
  int grade1;

  /** Event grade. This gives information about the energy and time
      resolution of the event. There are the following event grades:

      -  0: nominal time and energy resolution 

      -  1: degraded energy information

      - -1: lost during reset

  */
  int grade2;
  
  /** Pile-up flag. If there is another event within the fast shaping
      time, both events cannot be distinguished from each other. In
      that case the pile-up flag is set to 1, otherwise it is 0. */
  int pileup;

  /** Exact impact position of the photon on the detector. */
  double x, y; 

} HTRSEvent;


#endif /* HTRSEVENT_H */

