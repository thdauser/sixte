#ifndef XMSEVENT_H
#define XMSEVENT_H 1


/** XMS-specific event. */
typedef struct {
  double time;
  long pha;
  int xi, yi;
} XMSEvent;


#endif /* XMSEVENT_H */

