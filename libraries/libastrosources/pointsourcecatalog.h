#ifndef POINTSOURCECATALOG_H
#define POINTSOURCECATALOG_H (1)

#include "sixt.h"
#include "pointsources.h"
#include "pointsourcefile.h"
#include "vector.h"


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** Contains all PointSources from one and the same FITS file. This
    data structure is generated from the sources sorted out of the
    PointSourceFile. */
typedef struct {

  /** Number of PointSource objects contained in the PointSourceCatalog. */
  long nsources;

  /** Array containing the individual PointSource objects. Length of
      the array is nsources. */
  PointSource* sources;

} PointSourceCatalog;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Selects point sources along the path of the telescope axis over
    the sky and returns a PointSourceCatalog containing the individual
    PointSource objects. The selection is done based on a normal
    vector perpendicular to the scan plane and the cos(sin?)-value of
    the maximum angle between the source direction and the scan
    plane. To each selected source a spectrum and a light curve are
    assigned. */
PointSourceCatalog* getPointSourceCatalog(PointSourceFile* psf, 
					  Vector normal_vector, 
					  const double max_align, 
					  int* status);

/** Destructor. */
void free_PointSourceCatalog(PointSourceCatalog* psc);


#endif /* POINTSOURCECATALOG_H */

