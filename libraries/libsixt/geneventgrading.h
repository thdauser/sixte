#ifndef GENEVENTGRADING_H 
#define GENEVENTGRADING_H 1

#include "sixt.h"


/////////////////////////////////////////////////////////////////
// Constants.
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// Type Declarations.
/////////////////////////////////////////////////////////////////

// Grading scheme:
// <grade p11="0" p12="0" p13="0" p21="0" p23="0" p31="0" p32="0" p33="0" grade="0"/>
// with the following values for the p_nm's:
// 0: contains no charge
// 1: contains charge above threshold
// 2: contains charge above threshold, which is higher than 
//    the charge in pixels with level 1
// and the following matrix indices (nm):
//   31  32  33
//   21  22  23
//   11  12  13

/** Individual event grade. */
typedef struct {
  /** Relative charge distribution in the 3x3 matrix. */
  int p[9];
  /** Assigned event grade (number). */
  int grade;
} GenEventGrade;


/** Information about event grading system. */
typedef struct {

  /** Array of pointers to event grades. */
  GenEventGrade** grades;

  /** Number of grades in the array. */
  int ngrades;

  /** Grade for invalid events. */
  int invalid;

  /** Flag whether events touching the border of the detector are
      declared as invalid. */
  int borderinvalid;
  
  /** Flag whether event patterns larger than a 3x3 matrix are
      declared as invalid. */
  int largeinvalid;

} GenEventGrading;


/////////////////////////////////////////////////////////////////
// Function Declarations.
/////////////////////////////////////////////////////////////////


/** Constructor. Allocates memory for a new empty GenEventGrading
    data structure. The first parameter gives the default event grade
    for all events that are not contained in the map and are therefore
    regarded as invalid. The second parameter is a flag, whether
    events touching the border of the detector are declared as
    invalid, since charge information might have been lost. The third
    parameter is a flag, whether patterns larger than a 3x3 matrix are
    declared as invalid. */
GenEventGrading* newGenEventGrading(const int invalid,
				    const int borderinvalid,
				    const int largeinvalid,
				    int* const status);

/** Destructor. Releases all allocated memory and resets the pointer
    to the GenEventGrading data structure to NULL. */
void destroyGenEventGrading(GenEventGrading** const grading);

/** Constructor. */
GenEventGrade* newGenEventGrade(const int p11, const int p12, const int p13, 
				const int p21, const int p23,
				const int p31, const int p32, const int p33,
				const int grade, int* const status);

/** Destructor for GenEventGrade. */
void destroyGenEventGrade(GenEventGrade** const grade);

/** Add a new grade to the GenEventGrading list. */
void addGenEventGrade(GenEventGrading* const grading,
		      GenEventGrade* const grade,
		      int* const status);

/** Determine the event grade of a pattern with the given charge
    distribution (3x3 matrix starting from the bottom left via the
    bottom line, second line, and upper line to the top right). The
    3rd and 4th parameter indicate, whether the pattern touches the
    border of the detector pixel array and whether the total pattern
    is larger than the specified 3x3 matrix. */
int getGenEventGrade(GenEventGrading* const grading,
		     const float* const charges,
		     const int border,
		     const int large);
		     

#endif /* GENEVENTGRADING_H */
