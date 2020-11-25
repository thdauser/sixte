/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2007-2014 Christian Schmid, FAU
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

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

  /** Starting address and size of FITS file (if created in core memory) */
  void* memptr;
  size_t memsize;

  /** Total number of rows in the file. */
  long nrows;

  /** Column numbers. */
  int ctime, cframe, cpha, cpi, csignal, crawx, crawy, cra, cdec,
    cph_id, csrc_id, cnpixels, ctype, cpileup, csignals, cphas;

} EventFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty EventFile data
    structure. */
EventFile* newEventFile(int* const status);

/** Destructor. */
void freeEventFile(EventFile** const file, int* const status);

/** Initializes EventFile column numbers */
void getColNumsFromEventFile(EventFile* const file,	int* const status);

/** Create and open a new EventFile. */
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

/** Add a column to open EventFile. EventFile must be opened in
    READWRITE mode. ttype should be known in EventFile typedef,
    otherwise added column will be ignored by the other
    EventFile routines */
void addCol2EventFile(EventFile* const file,
		int* const colnum, char* const ttype,
		char* const tform, char* const tunit,
		int* const status);

/** Append a new event to the event file. */
void addEvent2File(EventFile* const file,
		   Event* const event,
		   int* const status);

/** Read the Event at the specified row from the file. The
    numbering for the rows starts at 1 for the first line. */
void getEventFromFile(const EventFile* const file,
		      const int row, Event* const event,
		      int* const status);

/** Update the Event at the specified row in the file. The
    numbering for the rows starts at 1 for the first line. */
void updateEventInFile(const EventFile* const file,
		       const int row, Event* const event,
		       int* const status);

/** Fill the destination EventFile with data from the source
    EventFile. The specified thesholds are applied to the transferred
    events. */
void copyEventFile(const EventFile* const src,
		   EventFile* const det,
		   const float threshold_lo_keV,
		   const float threshold_up_keV,
		   int* const status);

/** Creates a copy of infile where fptr is stored in core memory */
EventFile* copyEventFileMemory(EventFile* infile, int* const status);


#endif /* EVENTFILE_H */
