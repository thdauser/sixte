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
  
  /** Valid flag. If the event is a valid event (i.e., not piled-up),
      the valid flag is set to 1. If the event is invalid, the flag is
      set to -1. */
  int f_valid;

} WFIEvent;


#endif /* WFIEVENT_H */

