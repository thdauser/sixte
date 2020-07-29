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
	event_list->avgs_4samplesDerivative=NULL; //SIRENA
	event_list->Es_lowres=NULL; //SIRENA
	event_list->grading=NULL; //SIRENA
	event_list->phis=NULL; //SIRENA
	event_list->lagsShifts=NULL; //SIRENA
	event_list->bsln=NULL; //SIRENA
	event_list->rmsbsln=NULL; //SIRENA
	event_list->energies=NULL;
	event_list->grades1=NULL;
	event_list->grades2=NULL;
	event_list->ph_ids=NULL;
        event_list->pix_ids=NULL;//SIRENA
	event_list->risetimes=NULL;//SIRENA
	event_list->falltimes=NULL;//SIRENA

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
                free(event_list->avgs_4samplesDerivative); //SIRENA
                free(event_list->Es_lowres); //SIRENA
		free(event_list->grading); //SIRENA
		free(event_list->phis); //SIRENA
		free(event_list->lagsShifts); //SIRENA
		free(event_list->bsln); //SIRENA
		free(event_list->rmsbsln); //SIRENA
		free(event_list->grades1);
		free(event_list->grades2);
		free(event_list->ph_ids);
                free(event_list->pix_ids);//SIRENA
		free(event_list->risetimes);//SIRENA
		free(event_list->falltimes);//SIRENA
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
			free(event_list->avgs_4samplesDerivative);
			free(event_list->grading);
			free(event_list->phis);
			free(event_list->lagsShifts);
			free(event_list->bsln);
			free(event_list->rmsbsln);
			free(event_list->grades2);
			if(NULL!= event_list->ph_ids){
				free(event_list->ph_ids);
			}
			free(event_list->pix_ids);
			free(event_list->risetimes);
			free(event_list->falltimes);
		}
		event_list->energies = malloc(event_list->size*sizeof*(event_list->energies));
		if (NULL == event_list->energies){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->avgs_4samplesDerivative = malloc(event_list->size*sizeof*(event_list->avgs_4samplesDerivative)); //SIRENA
		if (NULL == event_list->avgs_4samplesDerivative){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}
		
		event_list->Es_lowres = malloc(event_list->size*sizeof*(event_list->Es_lowres)); //SIRENA
		if (NULL == event_list->Es_lowres){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->grading = malloc(event_list->size*sizeof*(event_list->grading)); //SIRENA
		if (NULL == event_list->grading){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->phis = malloc(event_list->size*sizeof*(event_list->phis)); //SIRENA
		if (NULL == event_list->phis){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->lagsShifts = malloc(event_list->size*sizeof*(event_list->lagsShifts)); //SIRENA
		if (NULL == event_list->lagsShifts){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->bsln = malloc(event_list->size*sizeof*(event_list->bsln)); //SIRENA
		if (NULL == event_list->bsln){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->rmsbsln = malloc(event_list->size*sizeof*(event_list->rmsbsln)); //SIRENA
		if (NULL == event_list->rmsbsln){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

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

		event_list->pix_ids = malloc(event_list->size*sizeof*(event_list->pix_ids));
		if (NULL == event_list->pix_ids){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for pix_ids array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->risetimes = malloc(event_list->size*sizeof*(event_list->risetimes)); //SIRENA
		if (NULL == event_list->risetimes){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

		event_list->falltimes = malloc(event_list->size*sizeof*(event_list->falltimes)); //SIRENA
		if (NULL == event_list->falltimes){
			*status=EXIT_FAILURE;
			SIXT_ERROR("memory allocation for energy array in TesEventList failed");
			CHECK_STATUS_VOID(*status);
		}

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
	file->avg_4samplesDerivativeCol =3; //SIRENA	
	file->E_lowresCol =4; //SIRENA	
	file->grade1Col =5;
	file->grade2Col =6;
	file->phiCol =7; //SIRENA
	file->lagsShiftCol =8; //SIRENA
	file->bslnCol =9; //SIRENA
	file->rmsbslnCol =10; //SIRENA
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

	//Write XML into header
	//char comment[MAXMSG];
	//sprintf(comment, "XMLFILE: %s", xmlfile);
	//fits_write_comment(file->fptr, comment, status);
	//CHECK_STATUS_RET(*status,file);

	// Create table

	//first column TIME
	char   *ttype[]={"TIME","SIGNAL","AVG4SD","ELOWRES","GRADE1","GRADE2","PHI","LAGS","BSLN","RMSBSLN","PIXID","PH_ID","RISETIME","FALLTIME","RA","DEC","DETX","DETY","GRADING","SRC_ID","N_XT","E_XT"}; //SIRENA
	char *tform[]={"1D",  "1D",    "1D",    "1D",    "1J",    "1J", "1D", "1J", "1D","1D" , "1J",   "1J", "1D","1D",  "1D","1D","1E","1E", "1I","1J","1I","1D"};
	char *tunit[]={"s",   "keV",   "",   "keV",      "",      "",      "",      "",     "",     "",       "",     "","s","s",     "deg","deg","m","m","","","","keV"};

	fits_create_tbl(file->fptr, BINARY_TBL, 0, 22,		// SIRENA (22 instead of 14)
			ttype, tform, tunit,"EVENTS", status);
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

	fits_write_key(file->fptr,TSTRING,"SIRENAV",sirenaVersion,NULL,status);
	CHECK_STATUS_RET(*status,file);

	//Write XML into header
	//char comment[MAXMSG];
	//sprintf(comment, "XMLFILE: %s", xmlfile);
	//fits_write_comment(file->fptr, comment, status);
	//CHECK_STATUS_RET(*status,file);

	// Create table

	//first column TIME
	char   *ttype[]={"TIME","SIGNAL","AVG4SD","ELOWRES","GRADE1","GRADE2","PHI","LAGS","BSLN","RMSBSLN","PIXID","PH_ID","RISETIME","FALLTIME","RA","DEC","DETX","DETY","GRADING","SRC_ID","N_XT","E_XT"}; //SIRENA
	char *tform[]={"1D",  "1D",    "1D",    "1D",    "1J",    "1J", "1D", "1J" , "1D", "1D","1J",   "1J", "1D","1D",  "1D","1D","1E","1E", "1I","1J","1I","1D"};
	char *tunit[]={"s",   "keV",   "",   "keV",      "",      "",      "",      "",      "",        "",     "",     "","s","s", "deg","deg","m","m","","","","keV"};

	fits_create_tbl(file->fptr, BINARY_TBL, 0, 22,		// SIRENA (22 instead of 14)
			ttype, tform, tunit,"EVENTS", status);
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

/** Opens a TES event file with the given mode */
TesEventFile* openTesEventFile(const char* const filename,const int mode, int* const status){
	TesEventFile * file = newTesEventFile(status);
	headas_chat(3, "open TES event file '%s' ...\n", filename);
	fits_open_table(&file->fptr, filename, mode, status);
	CHECK_STATUS_RET(*status, NULL);

	// Determine the row numbers.
	fits_get_num_rows(file->fptr, &(file->nrows), status);

	// Determine the column numbers.
	fits_get_colnum(file->fptr, CASEINSEN, "TIME", &file->timeCol, status);
	fits_get_colnum(file->fptr, CASEINSEN, "SIGNAL", &file->energyCol, status);
        fits_get_colnum(file->fptr, CASEINSEN, "AVG4SD", &file->avg_4samplesDerivativeCol, status);  //SIRENA
        fits_get_colnum(file->fptr, CASEINSEN, "ELOWRES", &file->E_lowresCol, status);  //SIRENA
	fits_get_colnum(file->fptr, CASEINSEN, "GRADE1", &file->grade1Col, status);
	if (*status==COL_NOT_FOUND) {
	  file->grade1Col=-1;
	  *status=0;
	}

	fits_get_colnum(file->fptr, CASEINSEN, "GRADE2", &file->grade2Col, status);
	if (*status==COL_NOT_FOUND) {
	  file->grade2Col=-1;
	  *status=0;
	}

	fits_get_colnum(file->fptr, CASEINSEN, "PHI", &file->phiCol, status);  //SIRENA
	fits_get_colnum(file->fptr, CASEINSEN, "LAGS", &file->lagsShiftCol, status);  //SIRENA
	fits_get_colnum(file->fptr, CASEINSEN, "BSLN", &file->bslnCol, status);  //SIRENA
	fits_get_colnum(file->fptr, CASEINSEN, "RMSBSLN", &file->rmsbslnCol, status);  //SIRENA

	fits_get_colnum(file->fptr, CASEINSEN, "PIXID", &file->pixIDCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->pixIDCol=-1;
	  *status=0;
	}
	fits_get_colnum(file->fptr, CASEINSEN, "PH_ID", &file->phIDCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->phIDCol=-1;
	  *status=0;
	}
	fits_get_colnum(file->fptr, CASEINSEN, "RISETIME", &file->riseCol, status);  //SIRENA
	fits_get_colnum(file->fptr, CASEINSEN, "FALLTIME", &file->fallCol, status);  //SIRENA
	fits_get_colnum(file->fptr, CASEINSEN, "RA", &file->raCol, status);
	fits_get_colnum(file->fptr, CASEINSEN, "DEC", &file->decCol, status);

	fits_get_colnum(file->fptr, CASEINSEN, "DETX", &file->detxCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->detxCol=-1;
	  *status=0;
	}
	fits_get_colnum(file->fptr, CASEINSEN, "DETY", &file->detyCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->detyCol=-1;
	  *status=0;
	}

	fits_get_colnum(file->fptr, CASEINSEN, "GRADING", &file->gradingCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->gradingCol=-1;
	  *status=0;
	}
	fits_get_colnum(file->fptr, CASEINSEN, "SRC_ID", &file->srcIDCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->srcIDCol=-1;
	  *status=0;
	}
	fits_get_colnum(file->fptr, CASEINSEN, "N_XT", &file->nxtCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->nxtCol=-1;
	  *status=0;
	}
	fits_get_colnum(file->fptr, CASEINSEN, "E_XT", &file->extCol, status);
	if (*status==COL_NOT_FOUND) {
	  file->extCol=-1;
	  *status=0;
	}
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

	//Save avgs_4samplesDerivative (AVG4SD) column   //SIRENA
 	fits_write_col(file->fptr, TDOUBLE, file->avg_4samplesDerivativeCol,
					file->row, 1, event_list->index, event_list->avgs_4samplesDerivative, status);
	CHECK_STATUS_VOID(*status);
        
        //Save Es_lowres (ELOWRES) column   //SIRENA
 	fits_write_col(file->fptr, TDOUBLE, file->E_lowresCol,
					file->row, 1, event_list->index, event_list->Es_lowres, status);
	CHECK_STATUS_VOID(*status);

	//Save phis (PHI) column   //SIRENA
 	fits_write_col(file->fptr, TDOUBLE, file->phiCol,
					file->row, 1, event_list->index, event_list->phis, status);
	CHECK_STATUS_VOID(*status);

	//Save lagsShifts (LAGS) column   //SIRENA
 	fits_write_col(file->fptr, TINT, file->lagsShiftCol,
					file->row, 1, event_list->index, event_list->lagsShifts, status);
	CHECK_STATUS_VOID(*status);

	//Save bsln (BSLN) column   //SIRENA
 	fits_write_col(file->fptr, TDOUBLE, file->bslnCol,
					file->row, 1, event_list->index, event_list->bsln, status);
	CHECK_STATUS_VOID(*status);

	//Save rmsbsln (RMSBSLN) column   //SIRENA
 	fits_write_col(file->fptr, TDOUBLE, file->rmsbslnCol,
					file->row, 1, event_list->index, event_list->rmsbsln, status);
	CHECK_STATUS_VOID(*status);

	//Save grading (GRADING) column   //SIRENA
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
	if(NULL!=event_list->ph_ids){
		fits_write_col(file->fptr, TLONG, file->phIDCol,
						file->row, 1, event_list->index, event_list->ph_ids, status);
		CHECK_STATUS_VOID(*status);
	}
	//printf("%s %d %s","event_list->ph_ids:", event_list->ph_ids[0],"\n");
	/*for (int i = 0 ; i<event_list->index ; i++){
		printf("%s %d %s","event_list->pix_ids:", event_list->pix_ids[i],"\n");
		printf("%s %e %s","event_list->energies:", event_list->energies[i],"\n");
	}*/

	//printf("%s %d %s","event_list->index:", event_list->index,"\n");
	//printf("%s %d %s","file->row:", file->row,"\n");
	//Save pix_ids column
	//fits_write_col(file->fptr, TINT, file->pixIDCol,
	//				file->row, 1, event_list->index, event_list->pix_ids, status);
	//CHECK_STATUS_VOID(*status);

	//Save riseTime (RISETIME) column   //SIRENA
 	fits_write_col(file->fptr, TDOUBLE, file->riseCol,
					file->row, 1, event_list->index, event_list->risetimes, status);
	CHECK_STATUS_VOID(*status);

	//Save fallTime (FALLTIME) column   //SIRENA
 	fits_write_col(file->fptr, TDOUBLE, file->fallCol,
					file->row, 1, event_list->index, event_list->falltimes, status);
	CHECK_STATUS_VOID(*status);

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

        /*//Save avg_4samplesDerivative column  //SIRENA
	double energy = (double)impact->energy;
	fits_write_col(file->fptr, TDOUBLE, file->energyCol,
			file->row, 1, 1, &energy, status);
=======
        //Save avg_4samplesDerivative column  //SIRENA
	/*double avg_4samplesDerivative = (double)impact->avg_4samplesDerivative;
	fits_write_col(file->fptr, TDOUBLE, file->avg_4samplesDerivativeCol,
			file->row, 1, 1, &avg_4samplesDerivative, status);
>>>>>>> bf8ce4341c3d8e276e191e58f674d692a0809835
	CHECK_STATUS_VOID(*status);*/

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
//void updateSignal(TesEventFile* file,long row,double energy,double avg_4samplesDerivative,long grade1,long grade2,int grading,int n_xt,double e_xt,int* const status){ //SIRENA
void updateSignal(TesEventFile* file,long row,double energy,long grade1,long grade2,int grading,int n_xt,double e_xt,int* const status){ 
	//Save energy column
	fits_write_col(file->fptr, TDOUBLE, file->energyCol,
			row, 1, 1, &energy, status);
	CHECK_STATUS_VOID(*status);

	/*//Save avg_4samplesDerivative column   //SIRENA
	fits_write_col(file->fptr, TDOUBLE, file->avg4samplesDerivativeCol,
=======
	//Save avg_4samplesDerivative column   //SIRENA
	/*fits_write_col(file->fptr, TDOUBLE, file->avg4samplesDerivativeCol,
>>>>>>> bf8ce4341c3d8e276e191e58f674d692a0809835
			row, 1, 1, &avg_4samplesDerivative, status);
	CHECK_STATUS_VOID(*status);*/

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
