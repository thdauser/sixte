#ifndef XRAYSOURCECATALOG_H
#define XRAYSOURCECATALOG_H 1

#include "sixt.h"
#include "gendet.h"
#include "kdtreeelement.h"
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


#endif /* XRAYSOURCECATALOG_H */
