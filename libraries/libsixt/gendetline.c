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

  // Set all pointers to NULL.
  line->charge=NULL;
  line->pileup=NULL;

  line->xwidth=0;
  line->anycharge=0;
  // Allocate memory for the pixels.
  line->charge = (float*)malloc(xwidth*sizeof(float));
  if (NULL==line->charge) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Memory allocation for GenDetLine failed!\n", *status);
    return(line);
  }
  // Allocate memory for the pile-up flag for each pixel.
  line->pileup = (GenPileupFlag*)malloc(xwidth*sizeof(GenPileupFlag));
  if (NULL==line->pileup) {
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


void destroyGenDetLine(GenDetLine** const line)
{
  if (NULL!=*line) {
    if (NULL!=(*line)->charge) {
      free((*line)->charge);
    }
    if (NULL!=(*line)->pileup) {
      free((*line)->pileup);
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
    int ii;
#pragma omp for
    for(ii=0; ii<line->xwidth; ii++) {
      line->charge[ii] = 0.;
      line->pileup[ii] = GP_NONE;
    }
    line->anycharge=0;
  }
}


void addGenDetLine(GenDetLine* const line0, const GenDetLine* const line1)
{
  // Add the charges.
  int ii;
#pragma omp for
  for(ii=0; ii<line1->xwidth; ii++) {
    line0->charge[ii] += line1->charge[ii];

    // Update the pile-up flags.
    if (line1->pileup[ii]!=GP_NONE) {
      line0->pileup[ii]=GP_PILEUP;
    }
  }
  // Set the anycharge flag.
  if (0==line0->anycharge) {
    line0->anycharge = line1->anycharge;
  }
}


int readoutGenDetLine(GenDetLine* const line, GenEvent* const event)
{
  if (0==line->anycharge) {
    return(0);
  } else {
    int i;
    for (i=0; i<line->xwidth; i++) {
      if (line->charge[i]>0.) {
	// Return the pixel charge.
	event->rawx = i;
	event->charge = line->charge[i];
	if (GP_PILEUP==line->pileup[i]) {
	  event->pileup = 1;
	} else {
	  event->pileup = 0;
	}
	// Delete the charge in the pixel array.
	line->charge[i] = 0.;
	line->pileup[i] = GP_NONE;
	return(1);
      }
    }
    return(0);
  }
}


void addGenDetCharge2Pixel(GenDetLine* const line, 
			   const int column, 
			   float energy)
{
  line->charge[column] += energy;
  line->anycharge       = 1;
}



