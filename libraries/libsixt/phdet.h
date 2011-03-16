#ifndef PHDET_H
#define PHDET_H 1

#include "sixt.h"
#include "gendet.h"
#include "impact.h"
#include "impactlistfile.h"
#include "eventlistfile.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phdetGenDet(GenDet* const det,
		 ImpactListFile* const ilf,
		 EventListFile* const elf,
		 const double t0,
		 const double exposure,
		 int* const status);


#endif /* PHDET_H */
