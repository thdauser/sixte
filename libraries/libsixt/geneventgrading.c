#include "geneventgrading.h"


GenEventGrading* newGenEventGrading(const int invalid,
				    const int borderinvalid,
				    const int largeinvalid,
				    int* const status)
{
  headas_chat(5, "initialize empty GenEventGrading object ...\n");

  GenEventGrading* grading = (GenEventGrading*)malloc(sizeof(GenEventGrading));
  if (NULL==grading) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenEventGrading failed!\n", *status);
    return(grading);
  }

  // Initialize all pointers with NULL.
  grading->grades=NULL;

  // Set the initial values.
  grading->ngrades       = 0;
  grading->invalid       = invalid;
  grading->borderinvalid = borderinvalid;
  grading->largeinvalid  = largeinvalid;

  return(grading);
}



void destroyGenEventGrading(GenEventGrading** const grading) 
{
  if (NULL!=(*grading)) {
    if (NULL!=(*grading)->grades) {
      int ii;
      for (ii=0; ii<(*grading)->ngrades; ii++) {
	destroyGenEventGrade(&(*grading)->grades[ii]);
      }
    }

    free(*grading);
    *grading=NULL;
  }
}


GenEventGrade* newGenEventGrade(const int p11, const int p12, const int p13, 
				const int p21, const int p23,
				const int p31, const int p32, const int p33,
				const int grade, int* const status)
{
  GenEventGrade* ggrade = (GenEventGrade*)malloc(sizeof(GenEventGrade));
  if (NULL==ggrade) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenEventGrade failed!\n", *status);
    return(ggrade);
  }

  // Initialize all pointers with NULL.

  // Set the initial values.
  ggrade->p11 = p11;
  ggrade->p12 = p12;
  ggrade->p13 = p13;
  ggrade->p21 = p21;
  ggrade->p23 = p23;
  ggrade->p31 = p31;
  ggrade->p32 = p32;
  ggrade->p33 = p33;
  ggrade->grade = grade;

  return(ggrade);
}


void destroyGenEventGrade(GenEventGrade** const grade)
{
  if (NULL!=*grade) {
    free(*grade);
    *grade=NULL;
  }
}


void addGenEventGrade(GenEventGrading* const grading,
		      GenEventGrade* const grade,
		      int* const status)
{
  // Resize the allocated amount of memory.
  grading->ngrades++;
  grading->grades = (GenEventGrade**)realloc(grading->grades,
					     (grading->ngrades)*sizeof(GenEventGrade));
  if (NULL==grading->grades) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenEventGrade failed!\n", *status);
    return;
  }

  // Add the pointer to the grading structure to the array.
  grading->grades[grading->ngrades-1] = grade;
}


int getGenEventGrade(GenEventGrading* const grading,
		     const float* const charges,
		     const int border,
		     const int large)
{
  // Check if this a large event and if large events should be
  // handled as invalid.
  if ((large) && (grading->largeinvalid)) return(grading->invalid);

  // Check if this a border event and if border events should be
  // handled as invalid.
  if ((border) && (grading->borderinvalid)) return(grading->invalid);

  // Search for a matching grade pattern.
  int ii;
  for (ii=0; ii<grading->ngrades; ii++) {
    GenEventGrade* grade = grading->grades[ii];
    // Check if the charge distribution matches the current 
    // grade pattern.
    
    if (((charges[0]>0.) == (grade->p11>0)) &&
	((charges[1]>0.) == (grade->p12>0)) &&
	((charges[2]>0.) == (grade->p13>0)) &&
	((charges[3]>0.) == (grade->p21>0)) &&
	((charges[5]>0.) == (grade->p23>0)) &&
	((charges[6]>0.) == (grade->p31>0)) &&
	((charges[7]>0.) == (grade->p32>0)) &&
	((charges[8]>0.) == (grade->p33>0))) {
      return(grade->grade);
    }
  }
  
  // None of the listed pattern matches.
  return(grading->invalid);
}

