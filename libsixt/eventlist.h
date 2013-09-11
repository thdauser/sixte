#ifndef EVENTLIST_H
#define EVENTLIST_H (1)

// OBSOLETE / DEPRECATED file!!

#include "sixt.h"


/** Structure that contains all information, which is necessary to access
 * a generic event list FITS file. */
typedef struct {
  fitsfile *fptr; /**< Pointer to the FITS file containing the event list. */

  /** Current row in the table (starting at 0). 
   * If row is equal to 0, that means that so far no line has been read from the
   * event list file. The next line to read is the first line. After the reading
   * process row will have the value 1. */
  long row; 

  long nrows; /**< Total number of rows in the table. */
} EventList;


////////////////////////////////////////////////////


/** Opens an existing FITS file with a binary table event list.  Apart
    from opening the FITS file the function also determines the number
    of rows in the FITS table and initializes the EventList data
    structure.  The access_mode parameter can be either READONLY or
    READWRITE. */
int openEventList(EventList*, char* filename, int access_mode);

/** Close an open event list FITS file. */
int closeEventList(EventList*);

/** Checks whether the end of the event list is reached.  If the
    internal pointer of the EventList data structure points to the
    last line in the file, i.e., that is the formerly read line, or
    has an even higher value, the function return value is 1,
    otherwise it is 0. */
int EventListEOF(EventList*);

/** Checks whether the specified row is valid with respect to the
    event file.  The return value is '1' if the row exists in the FITS
    file, otherwise it is '0'. */
int EventListRowIsValid(EventList*, long row);


#endif /* EVENTLIST_H */
