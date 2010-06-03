#ifndef POINTSOURCECATALOG_H
#define POINTSOURCECATALOG_H (1)

#include "sixt.h"
#include "pointsources.h"
#include "pointsourcefile.h"
#include "kdtree.h"
#include "pointsourcelist.h"
#include "vector.h"
#include "photon.h"

//#define POINTSOURCE_KDTREE 1


////////////////////////////////////////////////////////////////////////
// Type Declarations.
////////////////////////////////////////////////////////////////////////


/** ...  */
typedef struct {

  PointSourceFile* file;

#ifdef POINTSOURCE_KDTREE
  kdTree kdtree;
#else 
  PointSourceList* psl;
#endif

} PointSourceCatalog;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


PointSourceCatalog openPointSourceCatalog(char* filename, int hdu, int* status);

void preselectPointSources(PointSourceCatalog* psc, Vector normal, 
			   double max_align, int* status);

void getFoVPointSources(PointSourceCatalog* psc,
			Vector* ref, double min_align, 
			double time, double dt, 
			struct PhotonOrderedListEntry** list_first,
			struct RMF* rmf,
			int* status);

void clearPointSourceCatalog(PointSourceCatalog* psc);


// OBSOLETE
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

// OBSOLETE
/** Destructor. */
void free_PointSourceCatalog(PointSourceCatalog* psc);


#endif /* POINTSOURCECATALOG_H */

