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


/** Load an RMF/RSP matrix and the corresponding EBOUNDS from a
    response file. If the compiler flag '-DNORMALIZE_RMF' is set, the
    RSP is renormalized to an RMF in such a way that the sum of each
    matrix row/column(?) is 1. */
struct ARF* loadARF(const char* const filename, int* const status);

/** Destructor for the ARF data structure. Warning: As there is no
    internal destructor for the ARF data structure in the HEASP
    library, the memory allocated by the function ReadARFMatrix is not
    realeased. */
void freeARF(struct ARF* const arf);


#endif /* ARF_H */
