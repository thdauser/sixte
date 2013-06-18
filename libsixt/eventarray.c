#include "eventarray.h"

ReadEvent* newEventArray(int* const status)
{
  ReadEvent* ea=(ReadEvent*)malloc(sizeof(ReadEvent));
  if (NULL==ea) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: Could not allocate memory for Event Array!\n",
		   *status);
    return(ea);
  }

  //Initialization:
  ea->EventArray=NULL;
  ea->rawx=0;
  ea->rawy=0;
  ea->charge=0.;

  ea->naxis1 = 0;
  ea->naxis2 = 0;

  return(ea);
}

ReadEvent* getEventArray(CodedMask* mask, SquarePixels* detector_pixels, int* status)
{
  ReadEvent* ea=NULL;
  int x,y;

   //Get empty event array-object
  ea=newEventArray(status);
  if (EXIT_SUCCESS!=*status) return(ea);

  //sizes are equal to those of the ReconArray -> mask-width but detector-pixelsize
  //mask has to be >= detector

  /*ea->naxis1=mask->naxis1;
    ea->naxis2=mask->naxis2;*/

    ea->naxis1=2*(mask->naxis1*mask->cdelt1/detector_pixels->xpixelwidth);
      ea->naxis2=2*(mask->naxis2*mask->cdelt2/detector_pixels->ypixelwidth);

  //memory-allocation
  ea->EventArray=(double**)malloc(ea->naxis1*sizeof(double*));
  if(NULL!=ea->EventArray){
    for(x=0; x < ea->naxis1; x++){
      ea->EventArray[x]=(double*)malloc(ea->naxis2*sizeof(double));
	if(NULL==ea->EventArray[x]) {
	  *status=EXIT_FAILURE;
	  HD_ERROR_THROW("Error: could not allocate memory to store the "
			 "EventArray!\n", *status);
	 
	  return(ea);
	}
	//Clear the pixels
	for(y=0; y < ea->naxis2; y++){
	  ea->EventArray[x][y]=0.;
	}
    }
    if (EXIT_SUCCESS!=*status) return(ea);
  } else {
      *status=EXIT_FAILURE;
      HD_ERROR_THROW("Error: could not allocate memory to store the "
		     "EventArray!\n", *status);
      return(ea);
  }//end of memory-allocation 
  return(ea);
  }


int readEventList_nextRow(CoMaEventFile* ef, ReadEvent* ea)
{
  int status=EXIT_SUCCESS;
  int anynul=0;

  // Move counter to next line.
  ef->generic.row++;

  // Check if there is still a row available.
  if (ef->generic.row > ef->generic.nrows) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: event list file contains no further entries!\n", 
		   status);
    return(status);
  }

  // Read in the data.
  ea->rawx = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawx, ef->generic.row, 1, 1, 
		    &ea->rawx, &ea->rawx, &anynul, &status))
    return(status);

  ea->rawy = 0;
  if (fits_read_col(ef->generic.fptr, TINT, ef->crawy, ef->generic.row, 1, 1, 
		    &ea->rawy, &ea->rawy, &anynul, &status))
    return(status);

  ea->charge = 0.;
  if (fits_read_col(ef->generic.fptr, TDOUBLE, ef->ccharge, ef->generic.row, 1, 1,
		    &ea->charge, &ea->charge, &anynul, &status)) 
    return(status);

  // Check if an error occurred during the reading process.
  if (0!=anynul) {
    status = EXIT_FAILURE;
    HD_ERROR_THROW("Error: reading from event list failed!\n", status);
    return(status);
  }

}

double* SaveEventArray1d(ReadEvent* ea, int* status)
{
 double* EventArray1d=NULL;

 //Memory-Allocation for 1d-image
 EventArray1d = (double*)malloc(ea->naxis1*ea->naxis2*sizeof(double));
 if (!EventArray1d) {
    *status = EXIT_FAILURE;
    HD_ERROR_THROW("Error allocating memory for 1d-event-array!\n", *status);
    return(EventArray1d);
 }

    //Create the 1D-image from EventArray 
  int x, y;
  for (x=0; x<ea->naxis1; x++) {
    for (y=0; y<ea->naxis2; y++) {
      EventArray1d[(x+ ea->naxis1*y)] = ea->EventArray[x][y];
   }
  }
  return(EventArray1d); 
}

void FreeEventArray1d(double* EventArray1d)
{
  if (EventArray1d!=NULL) free(EventArray1d);
}

void FreeEventArray(ReadEvent* ea)
{
  if (ea!=NULL) {
    if ((ea->naxis1>0)&&(NULL!=ea->EventArray)) {
      int count;
      for(count=0; count< ea->naxis1; count++) {
	if (NULL!=ea->EventArray[count]) {
	  free(ea->EventArray[count]);
	}
      }
      free(ea->EventArray);
    }
    free(ea);
  }
}
