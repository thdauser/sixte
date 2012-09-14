#ifndef ARF_H
#define ARF_H 1

#include "sixt.h"

#ifndef HEASP_H
#define HEASP_H 1
#include "heasp.h"
#endif


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Returns an empty ARF data structure with
    NULL-initialized pointers. */
struct ARF* getARF(int* const status);

/** Load an ARF from a response file. */
struct ARF* loadARF(char* filename, int* const status);

/** Destructor for the ARF data structure. Warning: As there is no
    internal destructor for the ARF data structure in the HEASP
    library, the memory allocated by the function ReadARF is released
    by a self-implemented function, which is not guaranteed to work
    properly. */
void freeARF(struct ARF* const arf);


#endif /* ARF_H */
