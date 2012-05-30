#ifndef GTI_H
#define GTI_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


typedef struct {
  /** Number of entries in the GTI storage. */
  unsigned long nentries;

  /** Array with the start times of the GTIs. */
  double* start;
  /** Array with the stop times of the GTIs. */
  double* stop;

} GTI;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns a pointer to an empty GTI data structure. */
GTI* newGTI(int* const status);

/** Destructor. */
void freeGTI(GTI** const file);

/** Load an existing GTI file. */
GTI* loadGTI(const char* const filename, int* const status);

/** Store the GTI collection in a FITS file. */
void saveGTI(GTI* const gti,
	     const char* const filename,
	     const char clobber,
	     int* const status);

/** Append a new GTI to the GTI collection. */
void appendGTI(GTI* const gti, 
	       const double start, 
	       const double stop, 
	       int* const status);


#endif /* GTI_H */
