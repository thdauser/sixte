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
		 Impact* const impact,
		 const double tend,
		 int* const status);


#endif /* PHDET_H */
