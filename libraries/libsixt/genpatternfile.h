#ifndef GENPATTERNFILE_H 
#define GENPATTERNFILE_H 1

#include "sixt.h"
#include "genpattern.h"
#include "eventlistfile.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Event file for the GenDet generic detector model. */
typedef struct {
  
  EventListFile* eventlistfile;

  /** Column numbers. */
  int cpat_type, cpileup, cphas;

} GenPatternFile;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty GenPatternFile data
    structure. */
GenPatternFile* newGenPatternFile(int* const status);

/** Destructor. */
void destroyGenPatternFile(GenPatternFile** const file, int* const status);

/** Create and open a new GenPatternFile. */
GenPatternFile* openNewGenPatternFile(const char* const filename,
				      int* const status);

/** Open an existing EventListFile. */
GenPatternFile* openGenPatternFile(const char* const filename,
				   const int mode, int* const status);

/** Append a new pattern to the event file. */
void addGenPattern2File(GenPatternFile* const file, 
			GenPattern* const pattern, 
			int* const status);


#endif /* GENPATTERNFILE_H */
