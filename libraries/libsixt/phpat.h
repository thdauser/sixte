#ifndef PHPAT_H
#define PHPAT_H 1

#include "sixt.h"
#include "event.h"
#include "eventlistfile.h"
#include "pattern.h"
#include "patternfile.h"
#include "gendet.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phpat(GenDet* const det,
	   EventListFile* const elf,
	   PatternFile* const plf,
	   int* const status);


#endif /* PHPAT_H */
