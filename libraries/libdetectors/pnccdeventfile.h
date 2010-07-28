#ifndef PNCCDEVENTFILE_H
#define PNCCDEVENTFILE_H 1

#include "sixt.h"
#include "pnccdevent.h"
#include "eventfile.h"

////////////////////////////////////////////////////////////
// Type declarations
///////////////////////////////////////////////////////////

typedef struct {
	EventFile generic; /**< Generic EventFile data structure */

	/* Column numbers of the individual pnCCD-specific event list 
		 entries. The numbers start at 1. The number 0 means, that there
		 is no corresponding column in the table. */
	int ctime, cpha, cenergy, crawx, crawy, cframe, cccdnr;
	int cpat_typ, cpat_inf;
	int cra, cdec, cskyx, cskyy;

} pnCCDEventFile;


////////////////////////////////////////////////////////////
// Function declarations
///////////////////////////////////////////////////////////

/** Opens an existing FITS file with a binary table event list. Apart
	from opening the FITS file the function also determines the number 
	of rows in the FITS table and initialize the pnCCDEventFile 
	data structure. The access_mode parameter can be either READONLY
	or READWRITE. */
int openpnCCDEventFile(pnCCDEventFile* pef, char* filename, int access_mode);

/** Create and open a new event list FITS file for the pnCCD 
	detector from a given FITS template. If the file already exists,
	the old file is deleted and replaced by an empty one. Apart from
	opening the FITS file the function also initializes the 
	pnCCDEventFile data structure. The access_mode is always READWRITE. */
int openNewpnCCDEventFile(pnCCDEventFile* pef, char* filename, char* template);

/** Close an open pnCCD event list FITS file */
int closepnCCDEventFile(pnCCDEventFile* pef);

/** Append a new pnCCD event to the event list. The inserted 
	event has a pixel numbering starting at 0, whereas the numbering
	in the event file RAWX and RAWY have to start at 1. So the 
	routines adds a 1 to the raw pixel coordinates. */
int addpnCCDEvent2File(pnCCDEventFile* pef, pnCCDEvent* event);

/** Read the next pnCCDEvent from the pnCCDEventFile. This 
	routine increases the internal counter of the pnCCDEventFile
	data structure. In the event file the numbering of RAWX and RAWY 
	starts at 1, whereas in the returned pnCCDEvent data structure
	the numbering of xi and yi starts at 0. The return value is the
	error status. */
int pnCCDEventFile_getNextRow(pnCCDEventFile* pef, pnCCDEvent* event);

/** Read a specific row from the pnCCDEventFile. This routine does 
	NOT increase the internal counter of the pnCCDEventFile data 
	structure. In the event file the numbering of RAWX and RAWY
	starts at 1, whereas in the returned pnCCDEvent data structure
	the numbering of xi and yi starts at 0. The return value of the
	funtion is the error status. */
int pnCCDEventFile_getRow(pnCCDEventFile* pef, pnCCDEvent* event, long row);

#endif /* PNCCDEVENTFILE_H */

