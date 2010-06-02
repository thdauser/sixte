#ifndef SOURCELIST_H
#define SOURCELIST_H 1

#include "sixt.h"
#include "vector.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Entry of a list of X-ray sources in 3-dimensional space. */
struct SourceListEntry {
  /* Location of the source in 3-dimensional space. */
  Vector location;

  /* Line number of the X-ray source in the FITS file. The line
     numbering starts at 0 for the first line. */
  long line;
};


/** List of X-ray sources in 3-dimensional space. */
typedef struct SourceListEntry SourceList;


////////////////////////////////////////////////////////////////////////
//   Function declarations
////////////////////////////////////////////////////////////////////////


/** Sort the SourceList with the specified number of entries with
    respect to the requested coordinate axis using a quick sort
    algorithm. */
void quicksortSourceList(SourceList* list, long left, long right, int axis);

void freeSourceList(SourceList* list);


#endif /* SOURCELIST_H */
