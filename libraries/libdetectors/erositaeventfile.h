#ifndef EROSITAEVENTFILE_H
#define EROSITAEVENTFILE_H 1

#include "sixt.h"
#include "erositaevent.h"
#include "eventfile.h"


typedef struct {
  EventFile generic; /**< Generic EventFile data structure. */

  /* Column numbers of the individual eROSITA-specific event list entries.
   * The numbers start at 1. The number 0 means, that there 
   * is no corresponding column in the table. */
  int ctime, cpha, cenergy, crawx, crawy, cframe, cccdnr;
  int cra, cdec, cskyx, cskyy;

} eROSITAEventFile;


/////////////////////////////////////////////////////////////////////


/** Opens an existing FITS file with a binary table event list.
 * Apart from opening the FITS file the function also determines the number of rows in 
 * the FITS table and initializes the eROSITAEventFile data structure. 
 * The access_mode parameter can be either READONLY or READWRITE.
 */
int openeROSITAEventFile(eROSITAEventFile*, char* filename, int access_mode);

/** Create and open a new event list FITS file for the eROSITA detector from 
 * a given FITS template.
 * If the file already exists, the old file is deleted and replaced by an empty one.
 * Apart from opening the FITS file the function also initializes the eROSITAEventFile 
 * data structure. 
 * The access_mode parameter is always READWRITE.
 */
int openNeweROSITAEventFile(eROSITAEventFile*, char* filename, char* template);

/** Close an open eROSITA event list FITS file. */
int closeeROSITAEventFile(eROSITAEventFile*);

/** Append a new eROSITA event to the to event list. 
 * The inserted event has a pixel numbering starting at 0, whereas the numbering
 * in the event file RAWX and RAWY have to start at 1. 
 * So the routine adds a 1 to the raw pixel coordinates. */
int addeROSITAEvent2File(eROSITAEventFile*, eROSITAEvent*);

/** Read the next eROSITAEvent from the eROSITAEventFile.
 * This routine increases the internal counter of the eROSITAEventFile data structure.
 * In the event file the numbering of RAWX and RAXY starts at 1, whereas in the returned
 * eROSITAEvent data structure the numbering of xi and yi starts at 0.
 * The return value is the error status. */
int eROSITAEventFile_getNextRow(eROSITAEventFile*, eROSITAEvent*);


#endif /* EROSITAEVENTFILE */
