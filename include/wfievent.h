#ifndef WFIEVENT_H
#define WFIEVENT_H 1


/** WFI-specific event. */
typedef struct {
  double time;
  long pha;
  int xi, yi;
  long frame;
  long patnum, patid;
  long pileup;
} WFIEvent;


#endif /* WFIEVENT_H */

