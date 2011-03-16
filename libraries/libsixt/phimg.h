#ifndef PHIMG_H
#define PHIMG_H 1

#include "sixt.h"
#include "attitudecatalog.h"
#include "check_fov.h"
#include "gendet.h"
#include "photon.h"
#include "photonlistfile.h"
#include "impact.h"
#include "impactlistfile.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phimg(const GenDet* const det,
	   AttitudeCatalog* const ac,
	   PhotonListFile* const plf,
	   ImpactListFile* const ilf,
	   const double t0,
	   const double exposure,
	   int* const status);


#endif /* PHIMG_H */
