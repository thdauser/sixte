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
*/

#include "testriggerfile.h"

/** Constructor. Returns a pointer to an empty TesTriggerFile data
    structure. */
TesTriggerFile* newTesTriggerFile(unsigned long triggerSize,int* const status) {
  TesTriggerFile* file=(TesTriggerFile*)malloc(sizeof(TesTriggerFile));
  if (NULL==file) {
    *status=EXIT_FAILURE;
    SIXT_ERROR("memory allocation for TesTriggerFile failed");
    return(file);
  }

  // Initialize pointers with NULL.
  file->fptr    =NULL;

  // Initialize values.
  file->nrows	     =0;
  file->row  	     =1;
  file->timeCol      =1;
  file->trigCol      =2;
  file->pixIDCol     =3;
  file->ph_idCol     =4;
  file->trigger_size =triggerSize;
  file->delta_t=0;

  return(file);
}

/** Destructor. */
void freeTesTriggerFile(TesTriggerFile** const file, int* const status){
  if (NULL!=*file) {
    (*file)->nrows	      =0;
    (*file)->row  	      =0;
    (*file)->timeCol          =0;
    (*file)->trigCol          =0;
    (*file)->ph_idCol         =0;
    (*file)->pixIDCol         =0;
    
    if (NULL!=(*file)->fptr) {
      fits_close_file((*file)->fptr, status);
      headas_chat(5, "closed TesStream list file (containing %ld rows).\n", 
		  (*file)->nrows);
    }
    free(*file);
    *file=NULL;
  }
}

/** Create and open a new TesTriggerFile. */
TesTriggerFile* opennewTesTriggerFile(const char* const filename,
				  SixtStdKeywords * keywords,
				  char* const xmlfile,
				  char* const impactlist,
				  unsigned long triggerSize,
				  int preBufferSize,
				  double sampleFreq,
				  const char clobber,
				  int* const status){

  TesTriggerFile* file = newTesTriggerFile(triggerSize,status);
  CHECK_STATUS_RET(*status, file);

  int exists;
  char buffer[MAXFILENAME];
  sprintf(buffer,"%s",filename);
  fits_file_exists(buffer, &exists, status);
  CHECK_STATUS_RET(*status,file);
  if (0!=exists) {
    if (0!=clobber) {
      // Delete the file.
      remove(buffer);
    } else {
      // Throw an error.
      char msg[MAXMSG];
      sprintf(msg, "file '%s' already exists", buffer);
      SIXT_ERROR(msg);
      *status=EXIT_FAILURE;
      CHECK_STATUS_RET(*status,file);
    }
  }
  fits_create_file(&file->fptr,buffer, status);
  CHECK_STATUS_RET(*status,file);
  int logic=(int)'T';
  int bitpix=8;
  int naxis=0;
  fits_update_key(file->fptr, TLOGICAL, "SIMPLE", &(logic), NULL, status);
  fits_update_key(file->fptr, TINT, "BITPIX", &(bitpix), NULL, status);
  fits_update_key(file->fptr, TINT, "NAXIS", &(naxis), NULL, status);
  sixt_add_fits_stdkeywords(file->fptr, 1,keywords, status);
  fits_update_key(file->fptr, TULONG, "TRIGGSZ", &(triggerSize), "Number of samples in triggers", status);
  fits_update_key(file->fptr, TINT, "PREBUFF", &(preBufferSize), "Number of samples before start of pulse", status);
  double deltat = 1./sampleFreq;
  fits_update_key(file->fptr, TDOUBLE, "DELTAT", &(deltat), "Time resolution of data stream", status);
  CHECK_STATUS_RET(*status,file);
  
  //Write XML into header
  char comment[MAXMSG];
  sprintf(comment, "XMLFILE: %s", xmlfile);
  fits_write_comment(file->fptr, comment, status);
  CHECK_STATUS_RET(*status,file);
  
  //Write pixel impact file into header
  sprintf(comment, "PIXFILE: %s", impactlist);
  fits_write_comment(file->fptr, comment, status);
  CHECK_STATUS_RET(*status,file);

   // Create table
  int tlen=9;
  
  char *ttype[4];
  char *tform[4];
  char *tunit[4];
  
  int ii;
  
  for(ii=0; ii<4; ii++){
    ttype[ii]=(char*)malloc(tlen*sizeof(char));
    if(ttype[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for ttype failed");
      CHECK_STATUS_RET(*status,file);
    }
    tform[ii]=(char*)malloc(tlen*sizeof(char));
    if(tform[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for tform failed");
      CHECK_STATUS_RET(*status,file);
    }
    tunit[ii]=(char*)malloc(tlen*sizeof(char));
    if(tunit[ii]==NULL){
      *status=EXIT_FAILURE;
      SIXT_ERROR("memory allocation for tunit failed");
      CHECK_STATUS_RET(*status,file);
    }
  }
  CHECK_STATUS_RET(*status,file);
  
  //Create first column TIME
  sprintf(ttype[0], "TIME");
  sprintf(tform[0], "1D");
  sprintf(tunit[0], "s");
  
  //Create second column ADC
  sprintf(ttype[1], "ADC");
  sprintf(tform[1], "%ldU",triggerSize);
  sprintf(tunit[1], "ADC");

  //Create third column PIXID
  sprintf(ttype[2], "PIXID");
  sprintf(tform[2], "1J");
  sprintf(tunit[2], "");
  CHECK_STATUS_RET(*status,file);

  //Create fourth column PH_ID
  sprintf(ttype[3], "PH_ID");
  sprintf(tform[3], "1PJ(%d)",MAXIMPACTNUMBER);
  sprintf(tunit[3], "");
  CHECK_STATUS_RET(*status,file);

  char extName[9];
  sprintf(extName,"RECORDS");
  fits_create_tbl(file->fptr, BINARY_TBL, 0, 4,
		  ttype, tform, tunit,extName, status);
  //Add keywords to other extension
  fits_update_key(file->fptr, TULONG, "TRIGGSZ", &(triggerSize), "Number of samples in triggers", status);
  fits_update_key(file->fptr, TINT, "PREBUFF", &(preBufferSize), "Number of samples before start of pulse", status);
  fits_update_key(file->fptr, TDOUBLE, "DELTAT", &(deltat), "Time resolution of data stream", status);
  if(keywords->extname!=NULL){
	  free(keywords->extname);
  }
  keywords->extname=strdup(extName);
  sixt_add_fits_stdkeywords(file->fptr, 2,keywords, status);
  CHECK_STATUS_RET(*status,file);

  int firstpix=0,lastpix=0,numberpix=0;
  float monoen=-1.;
  long nes_tot=0,net_tot=0;
  fits_update_key(file->fptr, TINT, "FIRSTPIX", &firstpix, "First pixel in record file", status);
  fits_update_key(file->fptr, TINT, "LASTPIX", &lastpix, "Last pixel in record file", status);
  fits_update_key(file->fptr, TINT, "NPIX", &numberpix, "Number of pixels in record file", status);
  fits_update_key(file->fptr, TFLOAT, "MONOEN", &monoen, "Monochromatic energy of photons [keV]", status);
  fits_update_key(file->fptr, TLONG, "NESTOT", &nes_tot, "Total number of events simulated", status);
  fits_update_key(file->fptr, TLONG, "NETTOT", &net_tot, "Total number of events actually triggered", status);
  CHECK_STATUS_RET(*status,file);

  //Free memory
  for(ii=0; ii<4; ii++){
    free(ttype[ii]);
    ttype[ii]=NULL;
    free(tform[ii]);
    tform[ii]=NULL;
    free(tunit[ii]);
    tunit[ii]=NULL;
  }

  return(file);

}

/** Open an existing TesTriggerFile. */
TesTriggerFile* openexistingTesTriggerFile(const char* const filename,SixtStdKeywords* keywords,int* const status){
	TesTriggerFile* file = newTesTriggerFile(0,status);
	CHECK_STATUS_RET(*status, file);

	//Open record file in READONLY mode
	fits_open_file(&(file->fptr), filename, READONLY, status);
	CHECK_STATUS_RET(*status, file);

	//Move to the primary hdu
	int hdu_type;
	fits_movabs_hdu(file->fptr,1, &hdu_type, status);
	CHECK_STATUS_RET(*status, file);

	//Read keywords
	sixt_read_fits_stdkeywords(file->fptr,keywords,status);
	CHECK_STATUS_RET(*status, file);

	//Move to the binary table
	fits_movabs_hdu(file->fptr,2, &hdu_type, status);
	CHECK_STATUS_RET(*status, file);

	//Get number of rows
	char comment[MAXMSG];
	fits_read_key(file->fptr, TINT, "NAXIS2", &(file->nrows), comment, status);
	CHECK_STATUS_RET(*status, file);

	//Get trigger_size
	fits_read_key(file->fptr, TULONG, "TRIGGSZ", &(file->trigger_size), comment, status);
	CHECK_STATUS_RET(*status, file);

	//Get delta_t
	fits_read_key(file->fptr, TDOUBLE, "DELTAT", &(file->delta_t), comment, status);
	CHECK_STATUS_RET(*status, file);


	//Associate column numbers
	fits_get_colnum(file->fptr, CASEINSEN,"TIME", &(file->timeCol), status);
	CHECK_STATUS_RET(*status, file);
	fits_get_colnum(file->fptr, CASEINSEN,"ADC", &(file->trigCol), status);
	CHECK_STATUS_RET(*status, file);
	fits_get_colnum(file->fptr, CASEINSEN,"PIXID", &(file->pixIDCol), status);
	CHECK_STATUS_RET(*status, file);
	fits_get_colnum(file->fptr, CASEINSEN,"PH_ID", &(file->ph_idCol), status);
	CHECK_STATUS_RET(*status, file);

	return(file);

}

/** Populates a TesRecord structure with the next record */
int getNextRecord(TesTriggerFile* const file,TesRecord* record,int* const status){
	int anynul=0;
	if (NULL==file || NULL==file->fptr) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("No opened trigger file to read from");
		CHECK_STATUS_RET(*status,0);
	}


	if (file->row<=file->nrows) {
		fits_read_col(file->fptr, TUSHORT, file->trigCol,
					  file->row,1,file->trigger_size,0,record->adc_array, &anynul,status);
		CHECK_STATUS_RET(*status,0);

//		fits_read_col(file->fptr, TLONG, file->ph_idCol,
//					  file->row,1,MAXIMPACTNUMBER,0,record->phid_array, &anynul,status);
//		CHECK_STATUS_RET(*status,0);

		fits_read_col(file->fptr, TLONG, file->pixIDCol,
					  file->row,1,1,0,&(record->pixid), &anynul,status);
		CHECK_STATUS_RET(*status,0);

		fits_read_col(file->fptr, TDOUBLE, file->timeCol,
					  file->row,1,1,0,&(record->time), &anynul,status);
		CHECK_STATUS_RET(*status,0);

		for (unsigned long i=0 ; i < file->trigger_size ; i++) {
			record->adc_double[i]= (double) (record->adc_array[i]);
		}

		file->row++;
		return(1);
	} else {
		return(0);
	}


}

/** Writes a record to a file */
void writeRecord(TesTriggerFile* outputFile,TesRecord* record,int* const status){
	fits_write_col(outputFile->fptr, TDOUBLE, outputFile->timeCol,
			outputFile->row, 1, 1, &(record->time), status);
	fits_write_col(outputFile->fptr, TUSHORT, outputFile->trigCol,
			outputFile->row, 1, record->trigger_size,record->adc_array, status);
	fits_write_col(outputFile->fptr, TLONG, outputFile->pixIDCol,
			outputFile->row, 1, 1, &(record->pixid), status);
	fits_write_col(outputFile->fptr, TLONG, outputFile->ph_idCol,
			outputFile->row, 1,record->phid_list->index,record->phid_list->phid_array, status);
	outputFile->nrows++;
	outputFile->row++;
}

