#ifndef POINTSOURCELIST_H
#define POINTSOURCELIST_H 1

#include "sixt.h"
#include "sourcelist.h"
#include "pointsources.h"
#include "pointsourcefile.h"
#include "kdtree.h"


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

  /** Number of PointSource objects contained in the
      PointSourceList. */
  long nsources;

  /** Array containing the individual PointSource objects. Length of
      the array is nsources. */
  PointSource* sources;

} PointSourceList;


////////////////////////////////////////////////////////////////////////
//   Function declarations
////////////////////////////////////////////////////////////////////////


PointSourceList* selectFoVPointSourceList(kdNode* tree, 
					  PointSourceFile* psf,
					  Vector* telescope_direction,
					  const double max_align,
					  int* status);

void clearPointSourceList(PointSourceList* psl);

#endif /* POINTSOURCELIST_H */
