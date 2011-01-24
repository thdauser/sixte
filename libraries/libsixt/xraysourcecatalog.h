#ifndef XRAYSOURCECATALOG_H
#define XRAYSOURCECATALOG_H 1

#include "sixt.h"
#include "gendet.h"
#include "kdtreeelement.h"
#include "linkedpholist.h"
#include "xraysource.h"
#include "xraysourcespectrum.h"


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Catalog of different X-ray sources. */
typedef struct {
  /** KDTree containing the XRaySource objects. */
  KDTreeElement* tree;

  /** Library containing all spectra. */
  XRaySourceSpectrum** spectra;
  long nspectra;

} XRaySourceCatalog;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. */
XRaySourceCatalog* newXRaySourceCatalog(int* const status);

/** Destructor. */
void freeXRaySourceCatalog(XRaySourceCatalog** cat);

/** Load a SIMPUT source catalog from a FITS file. */
XRaySourceCatalog* loadSourceCatalog(const char* const filename,
				     const GenDet* const det,
				     int* const status);

/** Create photons for all sources in the catalog for the specified
    time interval. Only sources within the FoV (given in [rad])
    defined by the telescope pointing direction are taken into
    account. */
LinkedPhoListElement* genFoVXRayPhotons(XRaySourceCatalog* const cat, 
					const Vector* const pointing, 
					const float fov,
					const double t0, const double t1,
					int* const status);


#endif /* XRAYSOURCECATALOG_H */
