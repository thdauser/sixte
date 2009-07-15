#ifndef EVENTLISTFILE_H
#define EVENTLISTFILE_H (1)


#include "sixt.h"
#include "eventlist.types.h"
#include "detectors.types.h"


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
} EventlistFile;


////////////////////////////////////////////////////


/*
// This function creates a new event list table in the specified FITS file.
// It also inserts  header information.
// The function returns '0', if it is run successfully.
// Otherwise the return value is '1'.
struct Eventlist_File* create_Eventlist_File(char* filename, Detector*,
					     double tstart, double tend, int *status);
*/


/** Opens an existing FITS file with a binary table event list.
 * The access_mode parameter can be either READONLY or READWRITE.
 */
struct Eventlist_File* open_EventlistFile(char* filename, int access_mode, int* status);


/** Opens an existing FITS file with a binary table event list.
 * Apart from opening the FITS file the function also determines the number of rows in 
 * the FITS table and initializes the EventlistFile data structure. 
 * The access_mode parameter can be either READONLY or READWRITE.
 */
int openEventlistFile(EventlistFile*, char* filename, int access_mode);

/** Close an open event list FITS file. */
int closeEventlistFile(EventlistFile*);


// This routine inserts one new line in the event list FITS table and writes
// the specified event data.
// The required parameters are:
// * a pointer to the fitsfile,
// * the row, after which the new line is inserted (starting at 0),
// * the event data like time, PHA value, grade and detector coordinates,
// * and the fits status variable for error handling.
void add_eventlist_row(struct Eventlist_File*, struct Event, int *status);


/** This function reads a row of data from the event list FITS file.
 * The function does NOT increment the row counter of the Eventlist_File object. */
int get_eventlist_row(struct Eventlist_File, struct Event*, int *status);


#endif /* EVENTLISTFILE_H */
