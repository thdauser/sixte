#ifndef XMSEVENTFILE_H
#define XMSEVENTFILE_H 1

#include "sixt.h"
#include "xmsevent.h"
#include "eventfile.h"


typedef struct {

  /** Generic EventFile data structure. */
  EventFile generic; 

  /* Column numbers of the individual XMS-specific event list entries.
      The numbers start at 1. The number 0 means, that there is no
      corresponding column in the table. */
  int ctime, cpha, crawx, crawy, cgrade, carray;

} XMSEventFile;


/////////////////////////////////////////////////////////////////////


/** Opens an existing FITS file with a binary table event list. Apart
    from opening the FITS file the function also determines the number
    of rows in the FITS table and initializes the XMSEventFile data
    structure. The access_mode parameter can be either READONLY or
    READWRITE.
 */
int openXMSEventFile(XMSEventFile*, char* filename, int access_mode);

/** Create and open a new FITS event file for the XMS detector from a
    given FITS template. If the file already exists, the old file is
    deleted and replaced by an empty one. Apart from opening the FITS
    file the function also initializes the XMSEventFile data structure
    by calling openXMSEventFile(). The access_mode parameter is
    always READWRITE.
 */
int openNewXMSEventFile(XMSEventFile*, char* filename, char* template);

/** Close an open XMS event list FITS file. */
int closeXMSEventFile(XMSEventFile*);

/** Append a new XMS event to the to event list. */
int addXMSEvent2File(XMSEventFile*, XMSEvent*);

/** Read the next XMSEvent from the XMSEventFile. This routine
    increases the internal counter of the XMSEventFile data structure.
    In the event file and in the returned XMSEvent data structure the
    numbering of RAWX and RAXY starts at 0. The return value of the
    function is the error status. */
int XMSEventFile_getNextRow(XMSEventFile*, XMSEvent*);

/** Read a specific row from the XMSEventfile. This routine does NOT
    increase the internal counter of the XMSEventFile data structure.
    In the event file and in the returned XMSEvent data structure the
    numbering of RAWX and RAXY starts at 0. The return value of the
    function is the error status. */
int XMSEventFile_getRow(XMSEventFile*, XMSEvent*, long row);

/** Write the specified column in the XMSEventFile. The data to be
    written is given by the XMSEvent. The row already must
    exist. Otherwise an error is returned. */
int XMSEventFile_writeRow(XMSEventFile*, XMSEvent*, long row);


#endif /* XMSEVENTFILE */
