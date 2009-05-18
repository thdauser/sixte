#ifndef EVENTLIST_H
#define EVENTLIST_H (1)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#include "fitsio.h"
#include "headas.h"
#include "headas_error.h"

#include "eventlist.types.h"
#include "detectors.types.h"
#include "sixt.h"



// This function creates a new event list table in the specified FITS file.
// It also inserts  header information.
// The function returns '0', if it is run successfully.
// Otherwise the return value is '1'.
int create_eventlist_file(struct Eventlist_File*, Detector*,
			  double tstart, double tend, int *status);


// Open an existing FITS file and try to get the first extension that contains
// a binary table.          !! OBSOLETE !!
int open_eventlist_file(struct Eventlist_File*, int* status);

/** Opens an existing FITS file with a binary table event list.
 * The access_mode parameter can be either READONLY or READWRITE.
 */
struct Eventlist_File* open_EventlistFile(char* filename, int access_mode, int* status);


// This routine inserts one new line in the event list FITS table and writes
// the specified event data.
// The required parameters are:
// * a pointer to the fitsfile,
// * the row, after which the new line is inserted (starting at 0),
// * the event data like time, PHA value, grade and detector coordinates,
// * and the fits status variable for error handling.
void add_eventlist_row(struct Eventlist_File*, struct Event, int *status);


// This function reads a row of data from the event list FITS file.
int get_eventlist_row(struct Eventlist_File, struct Event*, int *status);



#endif /* EVENTLIST_H */
