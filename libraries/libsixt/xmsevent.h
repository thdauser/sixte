#ifndef XMSEVENT_H
#define XMSEVENT_H 1


/** XMS-specific event. */
typedef struct {
  double time;
  long pha;
  int xi, yi;
  
  /** Energy Accuracy of this event. There are several categories
      designating the accuracy of the energy determination of this
      event: 0=unclassified, 1=high, 2=intermediate, ... */
  int grade;

  /** TES microcalorimeter pixel in which the event is measured. If
      the event is detected in the inner TES microcalorimeter array
      with the high energy resolution, array has the value 1. For the
      outer array with lower energy resolution it has the value 2. */
  int array;

} XMSEvent;


#endif /* XMSEVENT_H */

