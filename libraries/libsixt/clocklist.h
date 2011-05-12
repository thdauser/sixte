#ifndef CLOCKLIST_H 
#define CLOCKLIST_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


typedef enum {
  CL_NONE        = 0,
  CL_WAIT        = 1,
  CL_LINESHIFT   = 2,
  CL_READOUTLINE = 3,
  CL_CLEARLINE   = 4
} CLType;


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** ClockList for the operation of time-triggered detectors. */
typedef struct {
  
  /** Number of elements in the ClockList. */
  int nelements;

  /** Type of each ClockList element. */
  CLType* type;

  /** Array of ClockList elements. */
  void** list;

  /** Current element. Index starts at 0. */
  int element;

  /** Current time [s]. */
  double time;

  /** Current frame number. */
  long frame;

} ClockList;


/** ClockList element for wait operation. */
typedef struct {
  /** Time to wait [s]. */
  double time;
} CLWait;


/** ClockList element for line shift. */
typedef struct {
  int dummy;
} CLLineShift;


/** ClockList element for readout line. */
typedef struct {
  /** Detector line to be read out. In CCDs usually this is the 0th
      line. */
  int lineindex;
  
  /** Y-index assigned to the read-out line. Due to the line shifts in
      CCD detectors this value in general is different from
      lineidx. */
  int readoutindex;
} CLReadoutLine;


/** ClockList element for clearing a detector line. */
typedef struct {
  /** Detector line to be cleared. In CCDs usually this is the 0th
      line. */
  int lineindex;
} CLClearLine;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new empty ClockList data
    structure. */
ClockList* newClockList(int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the ClockList data structure to NULL. */
void destroyClockList(ClockList** const list);



/** Append a new element at the end of the ClockList. */
void append2ClockList(ClockList* const list, const CLType type, 
		      void* const element, int* const status);

/** Get the next ClockList element. If the specified list is empty, an
    error is returned via the error status. If the detector currently
    is in a wait status, a CL_NONE type is returned. The CL_WAIT type
    is return, if a wait period has been finished. */
void getClockListElement(ClockList* const list, const double time,
			 CLType* const type, void** const element,
			 int* const status);



/** Constructor for CLWait. */
CLWait* newCLWait(const double time, int* const status);
/** Destructor for CLWait. */
void destroyCLWait(CLWait** const clwait);


/** Constructor for CLLineShift. */
CLLineShift* newCLLineShift(int* const status);
/** Destructor for CLLineShift. */
void destroyCLLineShift(CLLineShift** const cllineshift);


/** Constructor for CLReadoutLine. */
CLReadoutLine* newCLReadoutLine(const int lineindex, const int readoutindex, 
				int* const status);
/** Destructor for CLReadoutLine. */
void destroyCLReadoutLine(CLReadoutLine** const clreadoutline);


/** Constructor for CLClearLine. */
CLClearLine* newCLClearLine(const int lineindex, int* const status);
/** Destructor for CLClearLine. */
void destroyCLClearLine(CLClearLine** const clclearline);


#endif /* CLOCKLIST_H */
