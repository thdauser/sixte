#include "gendetline.h"


GenDetLine* newGenDetLine(int xwidth, int* status)
{
  GenDetLine* line=NULL;

  // Check if the requested line length is positive.
  if (line <= 0) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Invalid line length for GenDetLine!\n", *status);
    return(line);
  }

  line = (GenDetLine*)malloc(sizeof(GenDetLine));
  if (NULL==line) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDetLine failed!\n", *status);
    return(line);
  }
  
  // Allocate memory for the pixels.
  line->xwidth=0;
  line->anycharge=0;
  line->charge = (float*)malloc(xwidth*sizeof(float));
  if (NULL==line->charge) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDetLine failed!\n", *status);
    return(line);
  }
  line->xwidth=xwidth;

  // Clear the pixels. Therefore we have to set the anycharge flag to 1,
  // otherwise the line will not be cleared, because the function assumes
  // that all pixel charges are already set to 0.
  line->anycharge=1;
  clearGenDetLine(line);

  return(line);
}


void destroyGenDetLine(GenDetLine** line)
{
  if (NULL!=*line) {
    if (NULL!=(*line)->charge) {
      free((*line)->charge);
    }
    free(*line);
    *line=NULL;
  }
}


void clearGenDetLine(GenDetLine* line)
{
  // Check if the line contains any charge. If not the clearing
  // is not necessary.
  if (1==line->anycharge) {
    int i;
    for(i=0; i<line->xwidth; i++) {
      line->charge[i] = 0.;
    }
    line->anycharge=0;
  }
}
