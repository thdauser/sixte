#ifndef EVENTFILE_H 
#define EVENTFILE_H 1

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

  /** Column numbers. */
  int ctime, cframe, cpi, csignal, crawx, crawy, cph_id, csrc_id;

} EventFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty EventFile data
    structure. */
EventFile* newEventFile(int* const status);

/** Destructor. */
void freeEventFile(EventFile** const file, int* const status);

/** Create and open a new EventFile. The new file is generated
    according to the specified template. */
EventFile* openNewEventFile(const char* const filename,
			    char* const telescop,
			    char* const instrume,
			    char* const filter,
			    char* const ancrfile,
			    char* const respfile,
			    const double mjdref,
			    const double timezero,
			    const double tstart,
			    const double tstop,
			    const int nxdim,
			    const int nydim,
			    const char clobber,
			    int* const status);

/** Open an existing EventFile. */
EventFile* openEventFile(const char* const filename,
			 const int mode, int* const status);

/** Append a new event at the end of the event file. */
void addEvent2File(EventFile* const file, Event* const event, 
		   int* const status);

/** Read the Event at the specified row from the event file. The
    numbering for the rows starts at 1 for the first line. */
void getEventFromFile(const EventFile* const file,
		      const int row, Event* const event,
		      int* const status);

/** Update the Event at the specified row in the event file. The
    numbering for the rows starts at 1 for the first line. */
void updateEventInFile(const EventFile* const file,
		       const int row, Event* const event,
		       int* const status);


#endif /* EVENTLISTFILE_H */
