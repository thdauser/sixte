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
      free((*grading)->grades);
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
  ggrade->p[0] = p11;
  ggrade->p[1] = p12;
  ggrade->p[2] = p13;
  ggrade->p[3] = p21;
  ggrade->p[4] = 1;
  ggrade->p[5] = p23;
  ggrade->p[6] = p31;
  ggrade->p[7] = p32;
  ggrade->p[8] = p33;
  // Make sure that the central pixel always contains the local maximum.
  int ii;
  for (ii=0; ii<9; ii++) {
    if (4==ii) continue;
    if (ggrade->p[ii] >= ggrade->p[4]) {
      ggrade->p[4] = ggrade->p[ii]+1;
    }
  }
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
    int matching = 1;
    int jj;
    for (jj=0; jj<9; jj++) {
      if ((charges[jj]>0.) != (grade->p[jj]>0)) {
	matching = 0;
	break;
      }
    }
    if (!matching) continue;

    // Check the relative charge distribution among the split partners.
    for (jj=0; jj<9; jj++) {
      int kk;
      for (kk=0; kk<9; kk++) {
	if (jj==kk) continue;
	if (grade->p[jj]>grade->p[kk]) {
	  if (charges[jj]<=charges[kk]) {
	    matching = 0;
	    break;
	  }
	}
      }
      if (!matching) break;
    }
    if (!matching) continue;

    // The charge pattern matches this grade.
    return(grade->grade);
  }
  
  // None of the listed pattern matches.
  return(grading->invalid);
}

