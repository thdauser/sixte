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
  int ctime, cpha, ccharge, crawx, crawy, cframe, cpileup;

} GenEventFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty GenEventFile data
    structure. */
GenEventFile* newGenEventFile(int* const status);

/** Destructor. */
void destroyGenEventFile(GenEventFile** const file, int* const status);

/** Create and open a new GenEventFile. The new file is generated
    according to the specified template. */
GenEventFile* openNewGenEventFile(const char* const filename,
				  const char* const template,
				  int* const status);

/** Open an existing GenEventFile. */
GenEventFile* openGenEventFile(const char* const filename,
			       const int mode, int* const status);

/** Append a new event to the event file. */
void addGenEvent2File(GenEventFile* const file, GenEvent* const event, 
		      int* const status);

/** Read the GenEvent at the specified row from the event file. The
    numbering for the rows starts at 1 for the first line. */
void getGenEventFromFile(const GenEventFile* const file,
			 const int row, GenEvent* const event,
			 int* const status);

/** Update the GenEvent at the specified row in the event file. The
    numbering for the rows starts at 1 for the first line. */
void updateGenEventFromFile(const GenEventFile* const file,
			    const int row, GenEvent* const event,
			    int* const status);


#endif /* GENEVENTFILE_H */
