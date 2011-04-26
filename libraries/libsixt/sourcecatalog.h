#ifndef SOURCECATALOG_H
#define SOURCECATALOG_H 1

#include "sixt.h"
#include "gendet.h"
#include "kdtreeelement.h"
#include "linkedpholist.h"
#include "simput.h"
#include "source.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Catalog of different X-ray sources. */
typedef struct {
  /** KDTree containing the Source objects. */
  KDTreeElement* tree;

  /** SIMPUT source catalog containing all relevant data. */
  SimputSourceCatalog* simput;

} SourceCatalog;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
SourceCatalog* newSourceCatalog(int* const status);

/** Destructor. */
void freeSourceCatalog(SourceCatalog** cat);

/** Load a SIMPUT source catalog from a FITS file. */
SourceCatalog* loadSourceCatalog(const char* const filename,
				     const GenDet* const det,
				     int* const status);

/** Create photons for all sources in the catalog for the specified
    time interval. Only sources within the FoV (given in [rad])
    defined by the telescope pointing direction are taken into
    account. */
LinkedPhoListElement* genFoVXRayPhotons(SourceCatalog* const cat, 
					const Vector* const pointing, 
					const float fov,
					const double t0, const double t1,
					int* const status);


#endif /* SOURCECATALOG_H */
