#ifndef EVENTLISTFILE_H
#define EVENTLISTFILE_H (1)


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
} EventlistFile;


////////////////////////////////////////////////////


/** Opens an existing FITS file with a binary table event list.
 * Apart from opening the FITS file the function also determines the number of rows in 
 * the FITS table and initializes the EventlistFile data structure. 
 * The access_mode parameter can be either READONLY or READWRITE.
 */
int openEventlistFile(EventlistFile*, char* filename, int access_mode);

/** Close an open event list FITS file. */
int closeEventlistFile(EventlistFile*);


/** Checks whether the end of the event list is reached. 
 * If the internal pointer of the EventlistFile data structure points to the last line
 * in the file, i.e. this is the formerly read line, or has an even higher value, the 
 * function return value is 1, otherwise it is 0. */
int EventlistFileEOF(EventlistFile*);


#endif /* EVENTLISTFILE_H */
