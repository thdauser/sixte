#ifndef HTRSEVENT_H
#define HTRSEVENT_H 1


/** HTRS-specific event. */
typedef struct {
  double time; /**< Detection time. */
  long pha; /**< PHA channel. */
  int pixel; /**< Number of the pixel, where the event is detected. Numbering starts at 0. */
} HTRSEvent;


#endif /* HTRSEVENT_H */

