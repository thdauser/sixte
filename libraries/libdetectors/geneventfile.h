#ifndef GENEVENTFILE_H 
#define GENEVENTFILE_H 1

#include "sixt.h"
#include "genevent.h"


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
  int ctime, cpha, ccharge, crawx, crawy, cpileup;
  
} GenEventFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty GenEventFile data
    structure. */
GenEventFile* newGenEventFile(int* const status);

/** Destructor. */
void destroyGenEventFile(GenEventFile** file, int* const status);

/** Create and open a new GenEventFile. The new file is generated
    according to the specified template. */
GenEventFile* openNewGenEventFile(const char* const filename,
				  const char* const template,
				  int* const status);

/** Append a new event to the event file. */
void addGenEvent2File(GenEventFile* const file, GenEvent* const event, 
		      int* const status);


#endif /* GENEVENTFILE_H */
