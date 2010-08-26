#include "gendet.h"


/** Parse the GenDet definition from an XML file. */
static void parseGenDetXML(const char* const filename, int* const status);


GenDet* newGenDet(const char* const filename, int* const status) {
  GenDet* det=NULL;

  // Allocate memory.
  det=(GenDet*)malloc(sizeof(GenDet));
  if (NULL==det) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }

  // Initialize all pointers with NULL.
  det->line=NULL;
  det->rmf =NULL;

  // Read in the XML definition of the detector.
  parseGenDetXML(filename, status);
  
  // TODO RM
  det->xwidth=384;
  det->ywidth=384;

  // Allocate memory for the pixels.
  det->line=(GenDetLine**)malloc(det->ywidth*sizeof(GenDetLine*));
  if (NULL==det->line) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDet failed!\n", *status);
    return(det);
  }
  int i;
  for (i=0; i<det->ywidth; i++) {
    det->line[i] = newGenDetLine(det->xwidth, status);
    if (EXIT_SUCCESS!=*status) return(det);
  }

  return(det);
}



void destroyGenDet(GenDet** det)
{
  if (NULL!=*det) {
    if (NULL!=(*det)->line) {
      int i;
      for (i=0; i<(*det)->ywidth; i++) {
	destroyGenDetLine(&(*det)->line[i]);
      }
      free((*det)->line);
    }

    free(*det);
    *det=NULL;
  }
}



static void parseGenDetXML(const char* const filename, int* const status)
{
  // Open the specified file.
  FILE* xmlfile = fopen(filename, "r");
  if (NULL==xmlfile) {
    *status = EXIT_FAILURE;
    char msg[MAXMSG];
    sprintf(msg, "Error: Failed opening GenDet definition XML "
	    "file '%s' for read access!\n", filename);
    HD_ERROR_THROW(msg, *status);
    return;
  }


  // Close the file handler.
  fclose(xmlfile);
}

