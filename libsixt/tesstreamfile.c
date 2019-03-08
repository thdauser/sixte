/*
   This file is part of SIXTE.

   SIXTE is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.


   Copyright 2014 Philippe Peille, IRAP
   Copyright 2015-2019 Remeis-Sternwarte, Friedrich-Alexander-Universitaet
                       Erlangen-Nuernberg
*/

#include "tesstreamfile.h"

/** Constructor. Returns a pointer to an empty TesStreamFile data
    structure. */
TesStreamFile* newTesStreamFile(int* const status) {

  TesStreamFile* file=(TesStreamFile*)malloc(sizeof(TesStreamFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TesStreamFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr         =NULL;
  file->columnNumber =NULL;
  file->extNumber    =NULL;

  // Initialize values.
  file->nrows	         =0;
  file->row  	         =0;
  file->Npix	         =0;
  file->timeColumnNumber =0;

  return(file);

}

/** Destructor. */
void freeTesStreamFile(TesStreamFile** const file, int* const status){
  if (NULL!=*file) {
    if((*file)->columnNumber!=NULL){
      free((*file)->columnNumber);
      (*file)->columnNumber=NULL;
    }

    if((*file)->extNumber!=NULL){
      free((*file)->extNumber);
      (*file)->extNumber=NULL;
    }
    (*file)->nrows	      =0;
    (*file)->row  	      =0;
    (*file)->Npix	      =0;
    (*file)->timeColumnNumber =0;

    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
      headas_chat(5, "closed TesStream list file (containing %ld rows).\n",
		  (*file)->nrows);
    }
    free(*file);
    *file=NULL;
  }
}

/** Open an existing TesStreamFile. */
TesStreamFile* openTesStreamFile(const char* const filename,int pixlow,int pixhigh,
				 const int mode, int* const status){
  TesStreamFile* file=newTesStreamFile(status);
  CHECK_STATUS_RET(*status, file);

  headas_chat(5, "open tesstream file '%s' ...\n", filename);

  /////////////////////////////////////////////////////////
  //Iterate over the pixel IDs and check whether they are available
  /////////////////////////////////////////////////////////
  int Npix = pixhigh-pixlow+1;                  //Total number of pixels given
  file->Npix = Npix;
  int nhdus=0;                                  // Number of HDUs

  // Open the FITS file table for reading:
  if (fits_open_file(&(file->fptr), filename, mode, status)) return(file);

  // Get number of hdus in th file
  if (fits_get_num_hdus(file->fptr, &nhdus, status)) return(file);

  // Allocate memory for array countaining extension numbers and column numbers
  file->extNumber = (int*)malloc(Npix*sizeof(int));
  file->columnNumber = (int*)malloc(Npix*sizeof(int));

  //Initialize iteration
  int currentExtensionNumber = 2;
  int hdu_type=0;
  //Move to the first hdu
  if (fits_movabs_hdu(file->fptr,currentExtensionNumber, &hdu_type, status)) return(file);
  //Get time column number
  if (fits_get_colnum(file->fptr, CASEINSEN,"TIME", &(file->timeColumnNumber), status)) return(file);
  //Proper iteration
  for (int pixNumber = 0 ; pixNumber<Npix ; pixNumber++){
    do {
      char columnName[9];
      sprintf(columnName, "PXL%05d", pixlow+pixNumber+1);
      fits_get_colnum(file->fptr, CASEINSEN,columnName, &(file->columnNumber[pixNumber]), status);
      if (*status==COL_NOT_FOUND) {
	currentExtensionNumber++;
	if (currentExtensionNumber>nhdus) {
	  *status=EXIT_FAILURE;
	  char msg[MAXMSG];
	  sprintf(msg, "No data available for pixel %d in file '%s'",pixlow+pixNumber+1,
		  filename);
	  SIXT_ERROR(msg);
	  return(file);
	}
	//Move to next hdu
	if (fits_movabs_hdu(file->fptr,currentExtensionNumber, &hdu_type, status)) return(file);
      }
    } while(*status==COL_NOT_FOUND);
    file->extNumber[pixNumber] = currentExtensionNumber;
  }

  CHECK_STATUS_RET(*status, file);
  return(file);

}

/** Generate a TESDataStream from tesstreamfile. */
TESDataStream* generateTESDataStreamFromFile(TESDataStream* stream,TesStreamFile* file,double tstart_stream,
					     double tstart,double tstop, double sampleFreq,int* const status){

  // Set internal row counter to the beginning correct row in the file.
  file->row = (long)((tstart-tstart_stream)*sampleFreq);
  // Set the number of rows to read
  file->nrows = (long)((tstop-tstart)*sampleFreq);

  //Allocate TesDataStream
  allocateTESDataStream(stream,file->nrows,file->Npix,status);
  CHECK_STATUS_RET(*status,stream);
  ///////////////////////////
  //Populate TESDataStream
  ///////////////////////////

  //Populate adc_value
  int currentExtensionNumber = file->extNumber[0];
  int anynul=0;
  int hdu_type=0;
  int pixNumber;
  int tstep;
  if (fits_movabs_hdu(file->fptr,currentExtensionNumber, &hdu_type, status)) return(stream);
  for (pixNumber=0 ; pixNumber<(file->Npix) ; pixNumber++){
    //Move to correct extension number
    if (currentExtensionNumber!=file->extNumber[pixNumber]){
      if (fits_movabs_hdu(file->fptr,currentExtensionNumber, &hdu_type, status)) return(stream);
    }
    anynul=0;
    for (tstep=0;tstep<file->nrows;tstep++){
      fits_read_col(file->fptr, TUSHORT, file->columnNumber[pixNumber], file->row+tstep+1,1,1,
		    NULL, &(stream->adc_value[tstep][pixNumber]), &anynul, status);
      fflush(stdout);
    }
  }

  //Populate time
  anynul=0;
  fits_read_col(file->fptr, TDOUBLE, file->timeColumnNumber, file->row+1, 1, file->nrows,
		NULL, (stream->time), &anynul, status);

  CHECK_STATUS_RET(*status, stream);
  return(stream);
}
