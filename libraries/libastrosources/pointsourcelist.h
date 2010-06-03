#ifndef POINTSOURCELIST_H
#define POINTSOURCELIST_H 1

#include "sixt.h"
#include "pointsources.h"
#include "pointsourcefile.h"


#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Contains all PointSources within a certain region (usually the
    FoV). */
typedef struct {
  PointSource* sources;
  long nsources;
} PointSourceList;


/** Linked Point Source List. */
struct structLinkedPointSourceListEntry {
  PointSource* source;
  struct structLinkedPointSourceListEntry* next;
};

typedef struct structLinkedPointSourceListEntry LinkedPointSourceListEntry;
typedef LinkedPointSourceListEntry* LinkedPointSourceList;


////////////////////////////////////////////////////////////////////////
//   Function declarations
////////////////////////////////////////////////////////////////////////


void freePointSourceList(PointSourceList* psl);

/*
PointSourceList* selectFoVPointSourceList(kdNode* tree, 
					  PointSourceFile* psf,
					  Vector* telescope_direction,
					  const double max_align,
					  int* status);

void clearPointSourceList(PointSourceList* psl);
*/

#endif /* POINTSOURCELIST_H */
