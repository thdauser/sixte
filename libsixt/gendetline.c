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
    HD_ERROR_THROW("Error: Memory allocation for GenDetLine failed!\n", 
		   *status);
    return(line);
  }

  // Set all pointers to NULL.
  line->charge=NULL;
  line->ph_id =NULL;
  line->src_id=NULL;

  line->xwidth=0;
  line->anycharge=0;

  // Allocate memory for the pixel charges.
  line->charge = (float*)malloc(xwidth*sizeof(float));
  CHECK_NULL_RET(line->charge, *status,
		 "memory allocation for GenDetLine failed!\n", line);
  // Allocate memory for the photon and source IDs.
  line->ph_id  = (long**)malloc(xwidth*sizeof(long*));
  CHECK_NULL_RET(line->ph_id, *status,
		 "memory allocation for GenDetLine failed!\n", line);
  line->src_id = (long**)malloc(xwidth*sizeof(long*));
  CHECK_NULL_RET(line->src_id, *status,
		 "memory allocation for GenDetLine failed!\n", line);
  int ii;
  for(ii=0; ii<xwidth; ii++) {
    line->ph_id[ii]  = (long*)malloc(NEVENTPHOTONS*sizeof(long));
    CHECK_NULL_RET(line->ph_id[ii], *status,
		   "memory allocation for GenDetLine failed!\n", line);
    line->src_id[ii] = (long*)malloc(NEVENTPHOTONS*sizeof(long));
    CHECK_NULL_RET(line->src_id, *status,
		   "memory allocation for GenDetLine failed!\n", line);
    // Initialize.
    int jj;
    for (jj=0; jj<NEVENTPHOTONS; jj++) {
      line->ph_id[ii][jj]  = 0;
      line->src_id[ii][jj] = 0;
    }
  }

  line->xwidth=xwidth;

  // Clear the pixels. Therefore we have to set the anycharge flag to 1,
  // otherwise the line will not be cleared, because the function assumes
  // that all pixel charges are already set to 0.
  line->anycharge=1;
  for(ii=0; ii<line->xwidth; ii++) {
    line->charge[ii]=1.;
  }
  clearGenDetLine(line);

  return(line);
}


void destroyGenDetLine(GenDetLine** const line)
{
  if (NULL!=*line) {
    if (NULL!=(*line)->charge) {
      free((*line)->charge);
    }

    if (NULL!=(*line)->ph_id) {
      int ii;
      for (ii=0; ii<(*line)->xwidth; ii++) {
	if (NULL!=(*line)->ph_id[ii]) {
	  free((*line)->ph_id[ii]);
	}
      }
      free((*line)->ph_id);
    }

    if (NULL!=(*line)->src_id) {
      int ii;
      for (ii=0; ii<(*line)->xwidth; ii++) {
	if (NULL!=(*line)->src_id[ii]) {
	  free((*line)->src_id[ii]);
	}
      }
      free((*line)->src_id);
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
    for(ii=0; ii<line->xwidth; ii++) {
      if (line->charge[ii]>0.) {
	line->charge[ii]=0.;
	int jj;
	for (jj=0; jj<NEVENTPHOTONS; jj++) {
	  line->ph_id[ii][jj] = 0;
	  line->src_id[ii][jj]= 0;
	}
      }
    }
    line->anycharge=0;
  }
}


void addGenDetLine(GenDetLine* const line0, const GenDetLine* const line1)
{
  // Check if the 1st line, which is added to the 0th, contains any charges.
  if (0==line1->anycharge) return;

  // Set the anycharge flag.
  line0->anycharge = 1;

  // Add the charges.
  int ii;
  for(ii=0; ii<line1->xwidth; ii++) {
    if (line1->charge[ii]>0.) {
      line0->charge[ii] += line1->charge[ii];

      // Copy the photon and source IDs.
      int jj, kk;
      for (jj=0, kk=0; jj<NEVENTPHOTONS; jj++) {
	if (0==line1->ph_id[ii][kk]) break;
	if (0==line0->ph_id[ii][jj]) {
	  line0->ph_id[ii][jj]  = line1->ph_id[ii][kk];
	  line0->src_id[ii][jj] = line1->src_id[ii][kk];
	  kk++;
	}
      }
    }
  }
}


void addGenDetCharge2Pixel(GenDetLine* const line, 
			   const int column, 
			   const float signal,
			   const long ph_id,
			   const long src_id)
{
  // Set PH_ID and SRC_ID.
  if (line->charge[column]<0.001) {
    // If the charge collect in the pixel up to now is below 1eV,
    // overwrite the old PH_ID and SRC_ID by the new value.
    line->ph_id[column][0]  = ph_id;
    line->src_id[column][0] = src_id;

  } else if (signal>0.001) {
    // Only store the PH_ID and SRC_ID of the new contribution
    // if its signal is above 1eV.
    long ii;
    for (ii=0; ii<NEVENTPHOTONS; ii++) {
      if (0==line->ph_id[column][ii]) {
	line->ph_id[column][ii]  = ph_id;
	line->src_id[column][ii] = src_id;
	break;
      }
    }
  }

  // Add the signal.
  line->charge[column] += signal;
  line->anycharge       = 1;
}



