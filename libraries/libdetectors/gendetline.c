#include "gendetline.h"


GenDetLine* newGenDetLine(const int xwidth, int* const status)
{
  GenDetLine* line=NULL;

  // Check if the requested line length is positive.
  if (xwidth <= 0) {
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


void clearGenDetLine(GenDetLine* const line)
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


void addGenDetLine(GenDetLine* const line0, GenDetLine* const line1)
{
  // Add the charges.
  int i;
  for(i=0; i<line1->xwidth; i++) {
    line0->charge[i] += line1->charge[i];
  }
  // Set the anycharge flag.
  if (0==line0->anycharge) {
    line0->anycharge = line1->anycharge;
  }
}


void switchGenDetLines(GenDetLine** const line0, GenDetLine** const line1)
{
  GenDetLine* buffer = *line0;
  *line0 = *line1;
  *line1 = buffer;
}


int readoutGenDetLine(GenDetLine* const line, float* charge, int* x)
{
  if (0==line->anycharge) {
    return(0);
  } else {
    int i;
    for (i=0; i<line->xwidth; i++) {
      if (line->charge[i]>0.) {
	// Return the pixel charge.
	*x = i;
	*charge = line->charge[i];
	// Delete the charge in the pixel array.
	line->charge[i] = 0.;
	return(1);
      }
    }
    return(0);
  }
}


inline void addGenDetCharge2Pixel(GenDetLine* const line, 
				  const int column, 
				  float energy)
{
  line->charge[column] += energy;
  line->anycharge       = 1;
}



