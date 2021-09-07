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


   Copyright 2015 Philippe Peille, IRAP
*/

#include "teseventlist.h"

/** TesEventList constructor. Returns a pointer to an empty TesEventList data
    structure. */
TesEventList* newTesEventList(int* const status){
	TesEventList* event_list= malloc(sizeof*event_list);
	if (NULL==event_list) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for TesEventList failed");
		return(event_list);
	}
	// Initialize pointers with NULL
	event_list->event_indexes=NULL;
	event_list->pulse_heights=NULL;
	event_list->energies=NULL;
	event_list->grades1=NULL;
	event_list->grades2=NULL;
	event_list->ph_ids=NULL;

	// Initialize values
	event_list->size=0;
	event_list->size_energy=0;
	event_list->index=0;
     
	return(event_list);
}

TesEventList* newTesEventListSIRENA(int* const status){
	TesEventList* event_list= malloc(sizeof*event_list);
	if (NULL==event_list) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for TesEventList failed");
		return(event_list);
	}
	// Initialize pointers with NULL
	event_list->event_indexes=NULL;
	event_list->pulse_heights=NULL;
	event_list->avgs_4samplesDerivative=NULL; 
	event_list->Es_lowres=NULL; 
	event_list->grading=NULL; 
	event_list->phis=NULL; 
	event_list->lagsShifts=NULL; 
	event_list->bsln=NULL; 
	event_list->rmsbsln=NULL; 
	event_list->energies=NULL;
	event_list->grades1=NULL;
	event_list->grades2=NULL;
	event_list->ph_ids=NULL;
    event_list->ph_ids2=NULL;
    event_list->ph_ids3=NULL;
    event_list->pix_ids=NULL;
	event_list->risetimes=NULL;
	event_list->falltimes=NULL;

	// Initialize values
	event_list->size=0;
	event_list->size_energy=0;
	event_list->index=0;

	return(event_list);
}

/** TesEventList Destructor. */
void freeTesEventList(TesEventList* event_list){
	if (NULL!=event_list){
		free(event_list->event_indexes);
		free(event_list->pulse_heights);
		free(event_list->energies);
		free(event_list->grades1);
		free(event_list->grades2);
		free(event_list->ph_ids);
		free(event_list);
		event_list=NULL;
	}
}

void freeTesEventListSIRENA(TesEventList* event_list){
	if (NULL!=event_list){
		free(event_list->event_indexes);
		free(event_list->pulse_heights);
		free(event_list->energies);
        free(event_list->avgs_4samplesDerivative); 
        free(event_list->Es_lowres); 
		free(event_list->grading); 
		free(event_list->phis); 
		free(event_list->lagsShifts); 
		free(event_list->bsln); 
		free(event_list->rmsbsln); 
		free(event_list->grades1);
		free(event_list->grades2);
		free(event_list->ph_ids);
        free(event_list->ph_ids2);
        free(event_list->ph_ids3);
        free(event_list->pix_ids);
		free(event_list->risetimes);
		free(event_list->falltimes);
		free(event_list);
		event_list=NULL;
	}
}

/** Allocates memory for a TesEventList structure for the triggering stage:
 *  only event_index, pulse_height */
void allocateTesEventListTrigger(TesEventList* event_list,int size,int* const status){
	event_list->event_indexes = malloc(size*sizeof*(event_list->event_indexes));
	if (NULL == event_list->event_indexes){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for event_index array in TesEventList failed");
		CHECK_STATUS_VOID(*status);
	}

	event_list->pulse_heights = malloc(size*sizeof*(event_list->pulse_heights));
	if (NULL == event_list->pulse_heights){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for pulse_height array in TesEventList failed");
		CHECK_STATUS_VOID(*status);
	}

	event_list->grades1 = malloc(size*sizeof*(event_list->grades1));
	if (NULL == event_list->grades1){
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for grade1 array in TesEventList failed");
		CHECK_STATUS_VOID(*status);
	}

	event_list->size = size;
}

/** Allocates memory for the energy, grade1, and grade2 arrays according to
 *  current size */
void allocateWholeTesEventList(TesEventList* event_list,unsigned char allocate_ph,int* const status){
	//Only allocate if necessary
	if (event_list->size_energy < event_list->size) {
		//Free previous lists if they were already used
		if (NULL != event_list->energies) {
			free(event_list->energies);
			//free(event_list->grading);
			//free(event_list->phis);
			free(event_list->grades2);
			if(NULL!= event_list->ph_ids){
				free(event_list->ph_ids);
			}
		}
		event_list->energies = malloc(event_list->size*sizeof*(event_list->energies));
		if (NULL == event_list->energies){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		/*event_list->grading = malloc(event_list->size*sizeof*(event_list->grading)); //SIRENA
		if (NULL == event_list->grading){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}*/

		/*event_list->phis = malloc(event_list->size*sizeof*(event_list->phis)); //SIRENA
		if (NULL == event_list->phis){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}*/

		event_list->grades2 = malloc(event_list->size*sizeof*(event_list->grades2));
		if (NULL == event_list->grades2){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for grade2 array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		if (allocate_ph){
			event_list->ph_ids = malloc(event_list->size*sizeof*(event_list->ph_ids));
			if (NULL == event_list->ph_ids){
				*status=EXIT_FAILURE;
				SIXT_ERROR("memory allocation for ph_ids array in TesEventList failed");
				CHECK_STATUS_VOID(*status);
			}
		}

		/*event_list->pix_ids = malloc(event_list->size*sizeof*(event_list->pix_ids));
		if (NULL == event_list->pix_ids){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for pix_ids array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}*/

		event_list->size_energy = event_list->size;
	}
}

/** Appends the index and pulse_height lists to the list */
void addEventToList(TesEventList* event_list,int index,double pulse_height,int grade1,int* const status) {
	//If the list is not big enough, increase size
	if (event_list->index >= event_list->size) {
		//Increase size value
		event_list->size = event_list->size * 2;

		//Allocate new arrays
		int * new_event_index_array = realloc(event_list->event_indexes,(event_list->size)*sizeof(*(event_list->event_indexes)));
		if (!new_event_index_array){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation during size update of TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}
		double * new_pulse_height_array = realloc(event_list->pulse_heights,(event_list->size)*sizeof(*(event_list->pulse_heights)));
		if (!new_pulse_height_array){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation during size update of TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}
		int * new_grades1 = realloc(event_list->grades1,(event_list->size)*sizeof(*(event_list->grades1)));
		if (!new_grades1){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation during size update of TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		//Assign Struct arrays to the new ones
		event_list->event_indexes = new_event_index_array;
		event_list->pulse_heights = new_pulse_height_array;
		event_list->grades1 = new_grades1;
	}
	//Add data to lists
	event_list->event_indexes[event_list->index] = index;
	event_list->pulse_heights[event_list->index] = pulse_height;
	event_list->grades1[event_list->index] = grade1;


	event_list->index++;
}


/** TesEventFile constructor. Returns a pointer to an empty TesEventFile data
    structure. */
TesEventFile* newTesEventFile(int* const status){
	TesEventFile* file=malloc(sizeof(*file));
	if (NULL==file) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for TesEventFile failed");
		return(file);
	}

	// Initialize pointers with NULL.
	file->fptr    =NULL;

	// Initialize values.
	file->row  	    =1;
	file->nrows     =0;
	file->timeCol   =1;
	file->energyCol =2;
	
	file->grade1Col =3;
	file->grade2Col =4;
	file->pixIDCol  =5;
	file->phIDCol   =6;
	file->raCol     =7;
	file->decCol    =8;
	file->detxCol   =9;
	file->detyCol   =10;
	file->gradingCol=11;
	file->srcIDCol  =12;
	file->nxtCol    =13;
	file->extCol    =14;
    
	return(file);
}

TesEventFile* newTesEventFileSIRENA(int* const status){
	TesEventFile* file=malloc(sizeof(*file));
	if (NULL==file) {
		*status=EXIT_FAILURE;
		SIXT_ERROR("memory allocation for TesEventFile failed");
		return(file);
	}

	// Initialize pointers with NULL.
	file->fptr    =NULL;

	// Initialize values.
	file->row  	    =1;
	file->nrows     =0;
	file->timeCol   =1;
	file->energyCol =2;
	file->avg_4samplesDerivativeCol =3; 
	file->E_lowresCol =4; 
	file->grade1Col =5;
	file->grade2Col =6;
	file->phiCol =7; 
	file->lagsShiftCol =8; 
	file->bslnCol =9; 
	file->rmsbslnCol =10; 
	file->pixIDCol  =11;
	file->phIDCol   =12;
	file->riseCol   =13; 
	file->fallCol   =14; 
	file->raCol     =15;
	file->decCol    =16;
	file->detxCol   =17;
	file->detyCol   =18;
	file->gradingCol=19;
	file->srcIDCol  =20;
	file->nxtCol    =21;
	file->extCol    =22;
    
	return(file);
}

/** TesEventFile Destructor. */
void freeTesEventFile(TesEventFile* file, int* const status){
	if (NULL!=file) {
		if (NULL!=file->fptr) {
			fits_close_file(file->fptr, status);
			CHECK_STATUS_VOID(*status);
			headas_chat(5, "closed TesEventFile list file\n");
		}
		free(file);
		file=NULL;
	}
}

/** Create and open a new TesEventFile. */
TesEventFile* opennewTesEventFile(const char* const filename,
				  SixtStdKeywords* keywords,
				  const char clobber,
				  int* const status){
	TesEventFile* file = newTesEventFile(status);
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
	sixt_add_fits_stdkeywords(file->fptr,1,keywords,status);
	CHECK_STATUS_RET(*status,file);

    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    const char * chardate = asctime (timeinfo);
    char keyvalstr[1000];
    strcpy(keyvalstr,chardate);
    fits_write_key(file->fptr,TSTRING,"CREADATE",keyvalstr,NULL,status);
    CHECK_STATUS_RET(*status,file);

	// Create table

	//first column TIME
	char   *ttype[]={"TIME","SIGNAL","GRADE1","GRADE2","PIXID","PH_ID","RA","DEC","DETX","DETY","GRADING","SRC_ID","N_XT","E_XT"};
    char *tform[]={"1D", "1D", "1J", "1J", "1J", "1J", "1D", "1D", "1E", "1E", "1I", "1J", "1I", "1D"};
	char *tunit[]={"s", "keV", "", "", "", "", "deg", "deg", "m", "m", "", "", "", "keV"};

	fits_create_tbl(file->fptr, BINARY_TBL, 0, 14, ttype, tform, tunit,"EVENTS", status);
	sixt_add_fits_stdkeywords(file->fptr,2,keywords,status);
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

	return(file);

}

TesEventFile* opennewTesEventFileSIRENA(const char* const filename,
				  SixtStdKeywords* keywords,
				  const char* const sirenaVersion,
				  const char clobber,
				  int* const status){
	TesEventFile* file = newTesEventFileSIRENA(status);
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
	sixt_add_fits_stdkeywords(file->fptr,1,keywords,status);
	CHECK_STATUS_RET(*status,file);

    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    const char * chardate = asctime (timeinfo);
    char keyvalstr[1000];
    strcpy(keyvalstr,chardate);
    fits_write_key(file->fptr,TSTRING,"CREADATE",keyvalstr,NULL,status);
    CHECK_STATUS_RET(*status,file);

	fits_write_key(file->fptr,TSTRING,"SIRENAV",sirenaVersion,NULL,status);
	CHECK_STATUS_RET(*status,file);
    
	// Create table

	//first column TIME
	char   *ttype[]={"TIME","SIGNAL","AVG4SD","ELOWRES","GRADE1","GRADE2","PHI","LAGS","BSLN","RMSBSLN","PIXID","PH_ID","RISETIME","FALLTIME","RA","DEC","DETX","DETY","GRADING","SRC_ID","N_XT","E_XT"}; 
	char *tform[]={"1D", "1D", "1D", "1D", "1J", "1J", "1D", "1J", "1D", "1D", "1J", "3J", "1D", "1D", "1D", "1D", "1E", "1E", "1I", "1J", "1I", "1D"};
	char *tunit[]={"s", "keV", "", "keV", "", "", "", "", "", "", "", "", "s", "s", "deg", "deg", "m", "m", "", "", "", "keV"};

	fits_create_tbl(file->fptr, BINARY_TBL, 0, 22, ttype, tform, tunit,"EVENTS", status);
	sixt_add_fits_stdkeywords(file->fptr,2,keywords,status);
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

	return(file);

}

static int getColnumIfExists(fitsfile *fptr, int casesen, char* template, int* const status) {
  int colnum = 0;
  fits_get_colnum(fptr, casesen, template, &colnum, status);
  if (*status == COL_NOT_FOUND) {
    *status = 0;
    return -1;
  }

  return colnum;
}

/** Opens a TES event file with the given mode */
TesEventFile* openTesEventFile(const char* const filename,const int mode, int* const status){
	TesEventFile * file = newTesEventFile(status);
	headas_chat(3, "open TES event file '%s' ...\n", filename);
	fits_open_table(&file->fptr, filename, mode, status);
	CHECK_STATUS_RET(*status, NULL);

	// Determine the row numbers.
	fits_get_num_rows(file->fptr, &(file->nrows), status);

	// Determine the column numbers.
	file->timeCol = getColnumIfExists(file->fptr, CASEINSEN, "TIME", status);
	file->energyCol = getColnumIfExists(file->fptr, CASEINSEN, "SIGNAL", status);
	file->grade1Col = getColnumIfExists(file->fptr, CASEINSEN, "GRADE1", status);
	file->grade2Col = getColnumIfExists(file->fptr, CASEINSEN, "GRADE2", status);
	file->pixIDCol = getColnumIfExists(file->fptr, CASEINSEN, "PIXID", status);
	file->phIDCol = getColnumIfExists(file->fptr, CASEINSEN, "PH_ID", status);
	file->raCol = getColnumIfExists(file->fptr, CASEINSEN, "RA", status);
	file->decCol = getColnumIfExists(file->fptr, CASEINSEN, "DEC", status);
	file->detxCol = getColnumIfExists(file->fptr, CASEINSEN, "DETX", status);
	file->detyCol = getColnumIfExists(file->fptr, CASEINSEN, "DETY", status);
	file->gradingCol = getColnumIfExists(file->fptr, CASEINSEN, "GRADING", status);
	file->srcIDCol = getColnumIfExists(file->fptr, CASEINSEN, "SRC_ID", status);
	file->nxtCol = getColnumIfExists(file->fptr, CASEINSEN, "N_XT", status);
	file->extCol = getColnumIfExists(file->fptr, CASEINSEN, "E_XT", status);

	CHECK_STATUS_RET(*status, NULL);

	return(file);
}

/** Adds the data contained in the event list to the given file */
void saveEventListToFile(TesEventFile* file,TesEventList * event_list,
		double start_time,double delta_t,long pixID,int* const status){
	//Save time, PIXID and dummy grading column
	double time;
	int dummy_grading = 0;
	for (int i = 0 ; i<event_list->index ; i++){
		time = start_time + delta_t*event_list->event_indexes[i];
		fits_write_col(file->fptr, TDOUBLE, file->timeCol,
					   file->row, 1, 1, &time, status);
		//fits_write_col(file->fptr, TLONG, file->pixIDCol,
		//			   file->row, 1, 1, &pixID, status);
		fits_write_col(file->fptr, TINT, file->pixIDCol,
					file->row, 1, 1, &event_list->pix_ids[i], status);
		fits_write_col(file->fptr, TINT, file->gradingCol,
				file->row, 1, 1, &dummy_grading, status);
		CHECK_STATUS_VOID(*status);
		file->row++;
	}
	file->row = file->row - event_list->index;
	CHECK_STATUS_VOID(*status);

	//Save energy column
	fits_write_col(file->fptr, TDOUBLE, file->energyCol,
					file->row, 1, event_list->index, event_list->energies, status);
	CHECK_STATUS_VOID(*status);

 	/*//Save grading (GRADING) column   
 	fits_write_col(file->fptr, TINT, file->gradingCol,
					file->row, 1, event_list->index, event_list->grading, status);
	CHECK_STATUS_VOID(*status);*/

	//Save grade1 column
	fits_write_col(file->fptr, TINT, file->grade1Col,
					file->row, 1, event_list->index, event_list->grades1, status);
	CHECK_STATUS_VOID(*status);

	//Save grade2 column
	fits_write_col(file->fptr, TINT, file->grade2Col,
					file->row, 1, event_list->index, event_list->grades2, status);
	CHECK_STATUS_VOID(*status);

	//If PH_ID was computed, save it
	if(NULL!=event_list->ph_ids){
		fits_write_col(file->fptr, TLONG, file->phIDCol,
						file->row, 1, event_list->index, event_list->ph_ids, status);
		CHECK_STATUS_VOID(*status);
	}

	file->row = file->row + event_list->index;
	file->nrows+= event_list->index;
}

/** Adds the data contained in the event list to the given file */
void saveEventListToFileSIRENA(TesEventFile* file,TesEventList * event_list,
		double start_time,double delta_t,long pixID,int* const status){
	//Save time, PIXID and dummy grading column
	double time;
	int dummy_grading = 0;
	for (int i = 0 ; i<event_list->index ; i++){
		time = start_time + delta_t*event_list->event_indexes[i];
		fits_write_col(file->fptr, TDOUBLE, file->timeCol,
					   file->row, 1, 1, &time, status);
		//fits_write_col(file->fptr, TLONG, file->pixIDCol,
		//			   file->row, 1, 1, &pixID, status);
		fits_write_col(file->fptr, TINT, file->pixIDCol,
					file->row, 1, 1, &event_list->pix_ids[i], status);
		fits_write_col(file->fptr, TINT, file->gradingCol,
				file->row, 1, 1, &dummy_grading, status);
		CHECK_STATUS_VOID(*status);
		file->row++;
	}
	file->row = file->row - event_list->index;
	CHECK_STATUS_VOID(*status);

	//Save energy column
	fits_write_col(file->fptr, TDOUBLE, file->energyCol,
					file->row, 1, event_list->index, event_list->energies, status);
	CHECK_STATUS_VOID(*status);

	//Save avgs_4samplesDerivative (AVG4SD) column   
 	fits_write_col(file->fptr, TDOUBLE, file->avg_4samplesDerivativeCol,
					file->row, 1, event_list->index, event_list->avgs_4samplesDerivative, status);
	CHECK_STATUS_VOID(*status);

        //Save Es_lowres (ELOWRES) column   
 	fits_write_col(file->fptr, TDOUBLE, file->E_lowresCol,
					file->row, 1, event_list->index, event_list->Es_lowres, status);
	CHECK_STATUS_VOID(*status);

	//Save phis (PHI) column   
 	fits_write_col(file->fptr, TDOUBLE, file->phiCol,
					file->row, 1, event_list->index, event_list->phis, status);
	CHECK_STATUS_VOID(*status);

	//Save lagsShifts (LAGS) column   
 	fits_write_col(file->fptr, TINT, file->lagsShiftCol,
					file->row, 1, event_list->index, event_list->lagsShifts, status);
	CHECK_STATUS_VOID(*status);

	//Save bsln (BSLN) column   
 	fits_write_col(file->fptr, TDOUBLE, file->bslnCol,
					file->row, 1, event_list->index, event_list->bsln, status);
	CHECK_STATUS_VOID(*status);

	//Save rmsbsln (RMSBSLN) column  
 	fits_write_col(file->fptr, TDOUBLE, file->rmsbslnCol,
					file->row, 1, event_list->index, event_list->rmsbsln, status);
	CHECK_STATUS_VOID(*status);

	//Save grading (GRADING) column  
 	fits_write_col(file->fptr, TINT, file->gradingCol,
					file->row, 1, event_list->index, event_list->grading, status);
	CHECK_STATUS_VOID(*status);

	//Save grade1 column
	fits_write_col(file->fptr, TINT, file->grade1Col,
					file->row, 1, event_list->index, event_list->grades1, status);
	CHECK_STATUS_VOID(*status);

	//Save grade2 column
	fits_write_col(file->fptr, TINT, file->grade2Col,
					file->row, 1, event_list->index, event_list->grades2, status);
	CHECK_STATUS_VOID(*status);

	//If PH_ID was computed, save it
	/*if(NULL!=event_list->ph_ids){
		fits_write_col(file->fptr, TLONG, file->phIDCol,
						file->row, 1, event_list->index, event_list->ph_ids, status);
		CHECK_STATUS_VOID(*status);
	}*/
    //If PH_ID was computed, save it
	if((NULL!=event_list->ph_ids) && (event_list->ph_ids == 0)){
		fits_write_col(file->fptr, TLONG, file->phIDCol,
						file->row, 1, event_list->index, event_list->ph_ids, status);
		CHECK_STATUS_VOID(*status);
	}
	else if ((NULL!=event_list->ph_ids) && (event_list->ph_ids != 0)){
        int dimPH_ID = 3; // Length of the PH_ID column
        int *buffer = (int *) malloc(event_list->index*dimPH_ID*sizeof(int));
        for (int i=0;i<event_list->index;i++)
        {
            buffer[0+i*dimPH_ID] = event_list->ph_ids[i];
            buffer[1+i*dimPH_ID] = event_list->ph_ids2[i];
            buffer[2+i*dimPH_ID] = event_list->ph_ids3[i];
        }
        
        fits_write_col(file->fptr, TINT, file->phIDCol,
                            file->row, 1,event_list->index*3, buffer, status);  // Be careful TINT or TLONG!!!
        free(buffer);
    }
	
	file->row = file->row + event_list->index;
	file->nrows+= event_list->index;
}

/** Updates the RA, DEC and DETX/Y columns with the given coordinates */
void updateRaDecDetXY(TesEventFile* file,double ra, double dec, float detx,float dety,int* const status){
	double dbuffer=ra*180./M_PI;
	fits_write_col(file->fptr, TDOUBLE, file->raCol,
						file->row, 1, 1,&dbuffer, status);
	CHECK_STATUS_VOID(*status);
	dbuffer=dec*180./M_PI;
	fits_write_col(file->fptr, TDOUBLE, file->decCol,
						file->row, 1, 1,&dbuffer, status);
	CHECK_STATUS_VOID(*status);
	fits_write_col(file->fptr, TFLOAT, file->detxCol,
						file->row, 1, 1,&detx, status);
	CHECK_STATUS_VOID(*status);
	fits_write_col(file->fptr, TFLOAT, file->detyCol,
						file->row, 1, 1,&dety, status);
	CHECK_STATUS_VOID(*status);
}

/** Add event as reconstructed with the RMF method */
void addRMFImpact(TesEventFile* file,PixImpact * impact,int grade1,int grade2,int grading,int n_xt,double e_xt,int* const status){
	//Save time column
	fits_write_col(file->fptr, TDOUBLE, file->timeCol,
			file->row, 1, 1, &(impact->time), status);
	CHECK_STATUS_VOID(*status);

	//Save energy column
	double energy = (double)impact->energy;
	fits_write_col(file->fptr, TDOUBLE, file->energyCol,
			file->row, 1, 1, &energy, status);
	CHECK_STATUS_VOID(*status);

	//Save grade1 column
	if (file->grade1Col!=-1) {
	  fits_write_col(file->fptr, TLONG, file->grade1Col,
			 file->row, 1, 1, &grade1, status);
	  CHECK_STATUS_VOID(*status);
	}

	//Save grade2 column
	if (file->grade2Col!=-1) {
	  fits_write_col(file->fptr, TLONG, file->grade2Col,
			 file->row, 1, 1, &grade2, status);
	  CHECK_STATUS_VOID(*status);
	}

	//Save grading column
	if (file->gradingCol!=-1) {
	  fits_write_col(file->fptr, TINT, file->gradingCol,
			 file->row, 1, 1, &grading, status);
	  CHECK_STATUS_VOID(*status);
	}

	//Save PIXID column
	long pixID = impact->pixID+1;
	fits_write_col(file->fptr, TLONG, file->pixIDCol,
			file->row, 1, 1, &pixID, status);
	CHECK_STATUS_VOID(*status);

	//Save PH_ID column
	fits_write_col(file->fptr, TLONG, file->phIDCol,
			file->row, 1, 1, &(impact->ph_id), status);
	CHECK_STATUS_VOID(*status);

	//Save SRC_ID column
	if (file->srcIDCol!=-1) {
	  fits_write_col(file->fptr, TINT, file->srcIDCol,
			 file->row, 1, 1, &(impact->src_id), status);
	  CHECK_STATUS_VOID(*status);
	}

	//Save N_XT and E_XT columns
	fits_write_col(file->fptr, TINT, file->nxtCol,
			file->row, 1, 1, &n_xt, status);
	fits_write_col(file->fptr, TDOUBLE, file->extCol,
			file->row, 1, 1, &e_xt, status);
	CHECK_STATUS_VOID(*status);


	file->row++;
	file->nrows++;
}

/** Adds an event whose signal as not been evaluated yet (necessity in order to keep causality in event file) */
void addEmptyEvent(TesEventFile* file,PixImpact* impact, int* const status){
	//Save time column
	fits_write_col(file->fptr, TDOUBLE, file->timeCol,
			file->row, 1, 1, &(impact->time), status);
	CHECK_STATUS_VOID(*status);

	//Save PIXID column
	long pixID = impact->pixID+1;
	fits_write_col(file->fptr, TLONG, file->pixIDCol,
			file->row, 1, 1, &pixID, status);
	CHECK_STATUS_VOID(*status);

	//Save PH_ID column
	fits_write_col(file->fptr, TLONG, file->phIDCol,
			file->row, 1, 1, &(impact->ph_id), status);
	CHECK_STATUS_VOID(*status);

	//Save SRC_ID column
	fits_write_col(file->fptr, TINT, file->srcIDCol,
			file->row, 1, 1, &(impact->src_id), status);
	CHECK_STATUS_VOID(*status);
	file->row++;
	file->nrows++;
}

/** Update signal and grading columns of an event */
void updateSignal(TesEventFile* file,long row,double energy,long grade1,long grade2,int grading,int n_xt,double e_xt,int* const status){
	//Save energy column
	fits_write_col(file->fptr, TDOUBLE, file->energyCol,
			row, 1, 1, &energy, status);
	CHECK_STATUS_VOID(*status);

	//Save grade1 column
	fits_write_col(file->fptr, TLONG, file->grade1Col,
			row, 1, 1, &grade1, status);
	CHECK_STATUS_VOID(*status);

	//Save grade2 column
	fits_write_col(file->fptr, TLONG, file->grade2Col,
			row, 1, 1, &grade2, status);
	CHECK_STATUS_VOID(*status);

	//Save grading column
	fits_write_col(file->fptr, TINT, file->gradingCol,
			row, 1, 1, &grading, status);
	CHECK_STATUS_VOID(*status);

	//Save N_XT and E_XT columns
	fits_write_col(file->fptr, TINT, file->nxtCol,
			row, 1, 1, &n_xt, status);
	CHECK_STATUS_VOID(*status);

	fits_write_col(file->fptr, TDOUBLE, file->extCol,
			row, 1, 1, &e_xt, status);
	CHECK_STATUS_VOID(*status);
}
