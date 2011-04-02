#ifndef PHGEN_H
#define PHGEN_H 1

#include "sixt.h"
#include "attitudecatalog.h"
#include "gendet.h"
#include "photon.h"
#include "photonlistfile.h"
#include "sourcecatalog.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phgen(const GenDet* const det,
	   AttitudeCatalog* const ac,
	   SourceCatalog* const srccat,
	   PhotonListFile* const plf,
	   const double t0, 
	   const double exposure,
	   int* const status);


#endif /* PHGEN_H */
