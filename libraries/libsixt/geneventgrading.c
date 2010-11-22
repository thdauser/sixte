#include "geneventgrading.h"


GenEventGrading* newGenEventGrading(const int invalid,
				    const int borderinvalid,
				    const int largeinvalid,
				    int* const status)
{
  headas_chat(5, "initialize empty GenEventGrading object ...\n");

  GenEventGrading* ident = (GenEventGrading*)malloc(sizeof(GenEventGrading));
  if (NULL==ident) {
    *status=EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenEventGrading failed!\n", *status);
    return(ident);
  }

  // Initialize all pointers with NULL.

  // Set the initial values.
  ident->invalid       = invalid;
  ident->borderinvalid = borderinvalid;
  ident->largeinvalid  = largeinvalid;
  int ii;
  for (ii=0; ii<256; ii++) {
    ident->grade[ii] = invalid;
  }

  return(ident);
}



void destroyGenEventGrading(GenEventGrading** const ident) 
{
  if (NULL!=(*ident)) {
    free(*ident);
    *ident=NULL;
  }
}



void addGenPatGrade(GenEventGrading* const ident, 
		    const int code, const int grade)
{
  if (NULL==ident) return;

  if ((code>=0) && (code<256)) {
    ident->grade[code] = grade;
  }
}


int getGenPatGrade(GenEventGrading* const ident,
		   const int code, 
		   const int border,
		   const int large) 
{
  if (NULL==ident) return(0);

  // Check if the pattern codes is within the possible range 0..255.
  if ((code<0) || (code>255)) {
    return(ident->invalid);
  }
  
  // Check if the event is a border event, and whether these
  // are declared as invalid.
  if ((0!=border) && (0!=ident->borderinvalid)) {
    return(ident->invalid);
  }

  // Check if the pattern is larger than 3x3, and if it should
  // therefore be declared as invalid.
  if ((0!=large) && (0!=ident->largeinvalid)) {
    return(ident->invalid);
  }
  
  // Return the corresponding event grade according to the map.
  return(ident->grade[code]);
}

