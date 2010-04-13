#ifndef WFIEVENTFILE_H
#define WFIEVENTFILE_H 1

#include "sixt.h"
#include "wfievent.h"
#include "eventfile.h"


typedef struct {
  EventFile generic; /**< Generic EventFile data structure. */

  /* Column numbers of the individual WFI-specific event list entries.
     The numbers start at 1. The number 0 means, that there is no
     corresponding column in the table. */
  int ctime, cpha, crawx, crawy, cframe;
  int cpatnum, cpatid, cpileup, cf_valid;

  /** Number of detector columns and rows. */
  int columns, rows;

} WFIEventFile;


/////////////////////////////////////////////////////////////////////


/** Opens an existing FITS file with a binary table event list. Apart
    from opening the FITS file the function also determines the number
    of rows in the FITS table and initializes the WFIEventFile data
    structure. The access_mode parameter can be either READONLY or
    READWRITE.
 */
int openWFIEventFile(WFIEventFile*, char* filename, int access_mode);

/** Create and open a new event list FITS file for the WFI detector
    from a given FITS template. If the file already exists, the old
    file is deleted and replaced by an empty one. Apart from opening
    the FITS file the function also initializes the WFIEventFile data
    structure. The access_mode parameter is always READWRITE.
 */
int openNewWFIEventFile(WFIEventFile*, char* filename, char* template);

/** Close an open WFI event list FITS file. */
int closeWFIEventFile(WFIEventFile*);

/** Append a new WFI event to the to event list. */
int addWFIEvent2File(WFIEventFile*, WFIEvent*);

/** Read the next WFIEvent from the WFIEventFile. This routine
    increases the internal counter of the WFIEventFile data structure.
    In the event file and in the returned WFIEvent data structure the
    numbering of RAWX and RAXY starts at 0. The return value of the
    function is the error status. */
int WFIEventFile_getNextRow(WFIEventFile*, WFIEvent*);


#endif /* WFIEVENTFILE */
