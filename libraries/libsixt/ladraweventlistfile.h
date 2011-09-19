#ifndef LADRAWEVENTLISTFILE_H 
#define LADRAWEVENTLISTFILE_H 1

#include "sixt.h"
#include "ladrawevent.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event file for the GenDet generic detector model. */
typedef struct {
  /** Pointer to the FITS file. */
  fitsfile* fptr;

  /** Total number of rows in the file. */
  long nrows;

  /** Current row in the file. */
  long row;

  /** Column numbers. */
  int ctime, csignal, cpanel, cmodule, celement, canode, cph_id, csrc_id;

} LADRawEventListFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty LADRawEventListFile
    data structure. */
LADRawEventListFile* newLADRawEventListFile(int* const status);

/** Destructor. */
void freeLADRawEventListFile(LADRawEventListFile** const file, 
			     int* const status);

/** Create and open a new LADRawEventListFile. The new file is generated
    according to the specified template. */
LADRawEventListFile* openNewLADRawEventListFile(const char* const filename,
						int* const status);

/** Open an existing LADRawEventListFile. */
LADRawEventListFile* openLADRawEventListFile(const char* const filename,
					     const int mode, 
					     int* const status);

/** Append a new event to the event file. */
void addLADRawEvent2File(LADRawEventListFile* const file, 
			 LADRawEvent* const event, 
			 int* const status);

/** Read the LADRawEvent at the specified row from the event file. The
    numbering for the rows starts at 1 for the first line. */
void getLADRawEventFromFile(const LADRawEventListFile* const file,
			    const int row, 
			    LADRawEvent* const event,
			    int* const status);

/** Update the LADRawEvent at the specified row in the event file. The
    numbering for the rows starts at 1 for the first line. */
void updateLADRawEventInFile(const LADRawEventListFile* const file,
			     const int row, LADRawEvent* const event,
			     int* const status);


#endif /* LADRAWEVENTLISTFILE_H */
