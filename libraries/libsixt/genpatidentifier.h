#ifndef GENPATIDENTIFIER_H 
#define GENPATIDENTIFIER_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////


/** Bad pixel map. */
typedef struct {

  /** Array for the mapping of event patterns to grades. */
  int grade[256];
  // Pattern coding:
  //   32  64  128
  //    8   0   16
  //    1   2    4

  /** Grade for invalid events. */
  int invalid;

  /** Flag whether events touching the border of the detector are
      declared as invalid. */
  int borderinvalid;
  
  /** Flag whether event patterns larger than a 3x3 matrix are
      declared as invalid. */
  int largeinvalid;

} GenPatIdentifier;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new empty GenPatIdentifier
    data structure. The first parameter gives the default event grade
    for all events that are not contained in the map and are therefore
    regarded as invalid. The second parameter is a flag, whether
    events touching the border of the detector are declared as
    invalid, since charge information might have been lost. The third
    parameter is a flag, whether patterns larger than a 3x3 matrix are
    declared as invalid. */
GenPatIdentifier* newGenPatIdentifier(const int invalid,
				      const int borderinvalid,
				      const int largeinvalid,
				      int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenPatIdentifier data structure to NULL. */
void destroyGenPatIdentifier(GenPatIdentifier** const ident);

/** Add a new grade to the GenPatIdentifier list. */
void addGenPatGrade(GenPatIdentifier* const ident,
		    const int code, const int grade);

/** Determine the event grade of a pattern with the particular
    code. The 3rd and 4th parameter indicate, whether the pattern
    touches the border of the detector pixel array and whether the
    pattern is larger than a 3x3 matrix. */
int getGenPatGrade(GenPatIdentifier* const ident,
		   const int code, 
		   const int border,
		   const int large);


#endif /* GENPATIDENTIFIER_H */
