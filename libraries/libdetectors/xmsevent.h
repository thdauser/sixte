#ifndef XMSEVENT_H
#define XMSEVENT_H 1


/** XMS-specific event. */
typedef struct {
  double time;
  long pha;
  int xi, yi;
  
  /** Energy Accuracy of this event. 
   * There are several categories designating the accuracy of the energy determination
   * of this event: 0=unclassified, 1=high, 2=intermediate, ... */
  int grade;
} XMSEvent;


#endif /* XMSEVENT_H */

