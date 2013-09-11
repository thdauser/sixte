#ifndef PHPAT_H
#define PHPAT_H 1

#include "sixt.h"
#include "event.h"
#include "eventfile.h"
#include "pattern.h"
#include "patternfile.h"
#include "gendet.h"


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


void phpat(GenDet* const det,
	   EventFile* const elf,
	   PatternFile* const plf,
	   const char skip_invalids,
	   int* const status);


#endif /* PHPAT_H */
