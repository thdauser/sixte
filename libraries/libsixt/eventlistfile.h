#ifndef EVENTLISTFILE_H 
#define EVENTLISTFILE_H 1

#include "sixt.h"
#include "event.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event file for the GenDet generic detector model. */
typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;

  /** Total number of rows in the file. */
  long nrows;

  /** Current row in the file. */
  long row;

  /** Column numbers. */
  int ctime, cpha, ccharge, crawx, crawy, cframe;

} EventListFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty EventListFile data
    structure. */
EventListFile* newEventListFile(int* const status);

/** Destructor. */
void freeEventListFile(EventListFile** const file, int* const status);

/** Create and open a new EventListFile. The new file is generated
    according to the specified template. */
EventListFile* openNewEventListFile(const char* const filename,
				    const char* const template,
				    int* const status);

/** Open an existing EventListFile. */
EventListFile* openEventListFile(const char* const filename,
				 const int mode, int* const status);

/** Append a new event to the event file. */
void addEvent2File(EventListFile* const file, Event* const event, 
		   int* const status);

/** Read the Event at the specified row from the event file. The
    numbering for the rows starts at 1 for the first line. */
void getEventFromFile(const EventListFile* const file,
		      const int row, Event* const event,
		      int* const status);

/** Update the Event at the specified row in the event file. The
    numbering for the rows starts at 1 for the first line. */
void updateEventInFile(const EventListFile* const file,
		       const int row, Event* const event,
		       int* const status);


#endif /* GENEVENTFILE_H */
