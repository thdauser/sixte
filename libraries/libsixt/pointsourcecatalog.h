#ifndef POINTSOURCECATALOG_H
#define POINTSOURCECATALOG_H (1)

#include "sixt.h"
#include "pointsources.h"
#include "pointsourcefile.h"
#include "pointsourcelist.h"
#include "kdtree.h"
#include "vector.h"
#include "photon.h"


/** Use kdTree algorithm to sort the point soucres. */
#define POINTSOURCE_KDTREE 1

/** Maximum number of sources in the preselected source catalog. */
#define MAX_N_POINTSOURCES 15000000 // 15 million


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Catalog containing X-ray point sources.  */
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


/** Open a PointSourceCatalog from a FITS file with a specific HDU. */
PointSourceCatalog openPointSourceCatalog(char* filename, int hdu, int* status);

/** Preselect PointSources along the path of the telescope axis among
    all exisiting point sources in order to avoid numerical effort by
    scanning all sources at each step. */
void preselectPointSources(PointSourceCatalog* psc, Vector normal, 
			   double max_align, int* status);

/** Determine the PointSources within / close to the FoV and generate
    photons for them. */
void generateFoVPointSourcePhotons(PointSourceCatalog* psc,
				   Vector* ref, double min_align, 
				   double time, double dt, 
				   struct PhotonOrderedListEntry** list_first,
				   const struct ARF* const arf,
				   int* status);

/** Clear the PointSourceCatalog. */
void clearPointSourceCatalog(PointSourceCatalog* psc);


#endif /* POINTSOURCECATALOG_H */

