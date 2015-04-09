/************************************************************************************************

   This file is part of SIXTE/SIRENA software.

   SIXTE/SIRENA is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   SIXTE/SIRENA is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

   Copyright 2014:  This file has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

************************************************************************************************/
#include "integraSIRENA.h"
#include "genutils.h"
#include "tasksSIRENA.h"

/** Initializes the structure ReconstructInitSIRENA with the variables required for SIRENA reconstruction */

extern "C" void initializeReconstructionSIRENA(ReconstructInitSIRENA* reconstruct_init, char* const record_file,
		char* const library_file, double tauFall,int pulse_length, double scaleFactor, double samplesUp,
		double nSgms, int crtLib, int mode, double LrsT, double LbT,char* const noise_file, char* filter_domain, char* filter_method, int calibLQ,  double b_cF, double c_cF, double monoenergy, int interm, char* detectFile, char* filterFile, char* record_file2, double monoenergy2, char clobber, int* const status){

	// Load LibraryCollection structure if library file exists
	int exists=0;
	if(fits_file_exists(library_file, &exists, status)){
		EP_PRINT_ERROR("Error checking if library file exists",*status);
		return;
	}
	if(exists){
		reconstruct_init->library_collection = getLibraryCollection(library_file, status);
		if(*status){
		      EP_PRINT_ERROR((char*)"Error in getLibraryCollection",EPFAIL); 
		      *status=EPFAIL;return;
		}
	}else if(!exists && crtLib==0){
		EP_PRINT_ERROR((char*)"Error accessing library file: it does not exists ",EPFAIL); 
		*status=EPFAIL;return;
	}
	//Load NoiseSpec structure
	reconstruct_init->noise_spectrum = NULL;
	if(crtLib==0){
		exists=0;
		if(fits_file_exists(noise_file, &exists, status)){
			EP_PRINT_ERROR("Error checking if noise file exists",*status);
			return;
		}
		reconstruct_init->noise_spectrum = getNoiseSpec(noise_file, status);
		if(*status){
			EP_PRINT_ERROR((char*)"Error in getNoiseSpec",EPFAIL); 
			*status=EPFAIL;return;
		}
	}

	//Get pulse height  CHECK THAT!!!!!!!!!!!!!!!!!!!!!
	strcpy(reconstruct_init->record_file,record_file);
	strcpy(reconstruct_init->library_file,library_file);
	reconstruct_init->pulse_length 	= pulse_length;
	reconstruct_init->tauFall       = tauFall;
	reconstruct_init->scaleFactor  	= scaleFactor;
	reconstruct_init->samplesUp    	= samplesUp;
	reconstruct_init->nSgms        	= nSgms;
	reconstruct_init->crtLib       	= crtLib;
	reconstruct_init->mode		= mode;
	reconstruct_init->LrsT		= LrsT;
	reconstruct_init->LbT		= LbT;
	reconstruct_init->monoenergy 	= monoenergy;
	strcpy(reconstruct_init->FilterDomain,filter_domain);
	strcpy(reconstruct_init->FilterMethod,filter_method);
	reconstruct_init->calibLQ       = calibLQ;
	reconstruct_init->b_cF          = b_cF;
	reconstruct_init->c_cF          = c_cF;
	reconstruct_init->intermediate  = interm;
	strcpy(reconstruct_init->detectFile,detectFile);
	strcpy(reconstruct_init->filterFile,filterFile);
	strcpy(reconstruct_init->record_file2,record_file2);
	reconstruct_init->monoenergy2 	= monoenergy2;
	if(0!=clobber){
	    reconstruct_init->clobber = 1;
	}else{
	    reconstruct_init->clobber = 0;
	} 
}
extern "C" void reconstructRecordSIRENA(TesRecord* record, TesEventList* event_list, ReconstructInitSIRENA* reconstruct_init,  int lastRecord, int nRecord, PulsesCollection **pulsesAll, OptimalFilterSIRENA **optimalFilter, int* const status){
		// Inititalize structure PulsesCollection
		PulsesCollection* pulsesInRecord = new PulsesCollection;
		pulsesInRecord->ndetpulses = 0;
		pulsesInRecord->pulses_detected = NULL;
		
		PulsesCollection* pulsesAllAux = new PulsesCollection;
		pulsesAllAux->ndetpulses = 0;
		pulsesAllAux->pulses_detected = NULL;

		//
		// Detect pulses in record
		//
		// Check consistency of some input parameters
		if(reconstruct_init->pulse_length > record->trigger_size){
		    EP_PRINT_ERROR("Warning: pulse length is larger than record size. Pulse length set to maximum value (record size)",EPFAIL);
		}
		runDetect(record, lastRecord, *pulsesAll, &reconstruct_init, &pulsesInRecord);
		if (reconstruct_init->crtLib == 0)
		{
			// filter pulses and calculates energy
			runFilter(record, nRecord, lastRecord, &reconstruct_init, *pulsesAll, &pulsesInRecord, optimalFilter);

			if (reconstruct_init->mode == 1)
			{
				// calibrated energy of pulses
				runEnergy(&reconstruct_init, &pulsesInRecord);
			}
		}
		if (nRecord == 1)
		{
			(*pulsesAll)->ndetpulses = pulsesInRecord->ndetpulses;
			if((*pulsesAll)->pulses_detected != NULL) delete [] (*pulsesAll)->pulses_detected;
			(*pulsesAll)->pulses_detected = new PulseDetected[pulsesInRecord->ndetpulses];
			for (int i=0;i<(*pulsesAll)->ndetpulses;i++)
			{
				(*pulsesAll)->pulses_detected[i] = pulsesInRecord->pulses_detected[i];
			}
		}
		else
		{
			pulsesAllAux->ndetpulses = (*pulsesAll)->ndetpulses;
			pulsesAllAux->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];

			for (int i=0;i<(*pulsesAll)->ndetpulses;i++)
			{
				pulsesAllAux->pulses_detected[i] = (*pulsesAll)->pulses_detected[i];
			}
			(*pulsesAll)->ndetpulses = (*pulsesAll)->ndetpulses + pulsesInRecord->ndetpulses;

			if((*pulsesAll)->pulses_detected != NULL) delete [] (*pulsesAll)->pulses_detected;
			(*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];
			// save pulses detected in previous records
			for (int i=0;i<pulsesAllAux->ndetpulses;i++)
			{
				(*pulsesAll)->pulses_detected[i] = pulsesAllAux->pulses_detected[i];
			}
			// save pulses detected in current record
			for (int i=0;i<pulsesInRecord->ndetpulses;i++)
			{
				(*pulsesAll)->pulses_detected[i+pulsesAllAux->ndetpulses] = pulsesInRecord->pulses_detected[i];
			}
		}
		// Fill TesEventList structure
		event_list->index = pulsesInRecord->ndetpulses;
		event_list->energies = new double[event_list->index];
		event_list->grades1  = new int[event_list->index];
		event_list->grades2  = new int[event_list->index];
		event_list->pulse_heights  = new double[event_list->index];
		event_list->ph_ids   = new long[event_list->index];
		
		for (int ip=0; ip<pulsesInRecord->ndetpulses; ip++){
		    event_list->event_indexes[ip] = 
			  (int)((pulsesInRecord->pulses_detected[ip].Tstart - record->time)/record->delta_t);

		    if ((reconstruct_init->crtLib == 0) && (reconstruct_init->mode == 1))
		    {
		    	event_list->energies[ip] = pulsesInRecord->pulses_detected[ip].energy;
		    }
		    else if ((reconstruct_init->crtLib == 0) && (reconstruct_init->mode == 0))
		    {
		    	event_list->energies[ip] = 999.;
		    }

		    event_list->grades1[ip]  = pulsesInRecord->pulses_detected[ip].grade1;
		    event_list->grades2[ip]  = pulsesInRecord->pulses_detected[ip].grade2;
		    event_list->pulse_heights[ip]  = pulsesInRecord->pulses_detected[ip].pulse_height;
		    event_list->ph_ids[ip]   = 0;
		    
		}
		return;
}

/** Constructor. Returns a pointer to an empty ReconstructInitSIRENA data structure. */
extern "C" ReconstructInitSIRENA* newReconstructInitSIRENA(int* const status){
	
	ReconstructInitSIRENA* reconstruct_init = new ReconstructInitSIRENA;
	
	// Initialize pointers with NULL for SIRENA
	reconstruct_init->library_collection =NULL;
	reconstruct_init->noise_spectrum     =NULL;
	
	// Initialize values for SIRENA
	strcpy(reconstruct_init->library_file,"");
	reconstruct_init->pulse_length=0;	
	reconstruct_init->tauFall=0.;
	reconstruct_init->scaleFactor=0.;
	reconstruct_init->samplesUp=0.;
	reconstruct_init->nSgms=0;
	reconstruct_init->crtLib=0;
	reconstruct_init->mode=0;
	reconstruct_init->monoenergy = 0;
	reconstruct_init->LrsT = 0;
	reconstruct_init->LbT = 0;
	strcpy(reconstruct_init->FilterDomain,"");
	strcpy(reconstruct_init->FilterMethod,"");
	reconstruct_init->calibLQ=0;
	reconstruct_init->b_cF=0.;
	reconstruct_init->c_cF=0.;
	reconstruct_init->clobber=0;
	
	return(reconstruct_init);
  }
  
/** Destructor. */
extern "C" void freeReconstructInitSIRENA(ReconstructInitSIRENA* reconstruct_init){
	delete(reconstruct_init);
}

/** Constructor. Returns a pointer to an empty PulsesCollection data structure. */
extern "C" PulsesCollection* newPulsesCollection(int* const status){

	PulsesCollection* PulsesColl = new PulsesCollection;

	// Initialize pointers with NULL for SIRENA
	PulsesColl->pulses_detected =NULL;

	// Initialize values for SIRENA
	PulsesColl->ndetpulses=0;

	return(PulsesColl);
  }

/** Destructor. */
extern "C" void freePulsesCollection(PulsesCollection* PulsesColl){
	delete(PulsesColl);
}

/** Constructor. Returns a pointer to an empty OptimalFilterSIRENA data structure. */
extern "C" OptimalFilterSIRENA* newOptimalFilterSIRENA(int* const status){

	OptimalFilterSIRENA* OFilterColl = new OptimalFilterSIRENA;

	// Initialize pointers with NULL for SIRENA
	OFilterColl->ofilter =NULL;

	// Initialize values for SIRENA
	OFilterColl->ofilter_duration=0;
	OFilterColl->nrmfctr=0.0;

	return(OFilterColl);
  }

/** Destructor. */
extern "C" void freeOptimalFilterSIRENA(OptimalFilterSIRENA* OFilterColl){
	delete(OFilterColl);
}

//extern "C" void writeCalibKeywords(fitsfile* fptr, double const b_cF, double const c_cF, int* const status){
extern "C" void writeCalibKeywords(fitsfile* fptr, double b_cF, double c_cF, int* const status){
    if(fits_write_key(fptr, TDOUBLE, "B_CF", &b_cF, "Calculated Linear Calibration Factor",status))
	EP_PRINT_ERROR("Cannot write B_CF calibration factor keyword to file", EPFAIL);
    if(fits_write_key(fptr, TDOUBLE, "C_CF", &c_cF, "Calculated Quadratic Calibration Factor",status))
	EP_PRINT_ERROR("Cannot write C_CF calibration factor keyword to file", EPFAIL);
    return;
}

/** Create and retrieve a LibraryCollection from a file. */
LibraryCollection* getLibraryCollection(const char* const filename,int* const status){

	//Create LibraryCollection structure
	LibraryCollection* library_collection = new LibraryCollection;
	library_collection->ntemplates=0;
	library_collection->pulse_templates=NULL;
	library_collection->pulse_templates_B0=NULL;
	library_collection->pulse_templates_filder=NULL;
	library_collection->matched_filters=NULL;
	library_collection->matched_filters_B0=NULL;
	
	//Open FITS file in READONLY mode
	fitsfile* fptr = NULL;
	if (fits_open_file(&fptr, filename, READONLY, status)){ 					
		EP_PRINT_ERROR("Error opening library file",*status);
		return(library_collection);}

	//Move to the first hdu
	int extver = 1;
	char HDUname[12];
	strcpy(HDUname,"LIBRARY");
	if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status)){
		EP_PRINT_ERROR("Error moving to HDU LIBRARY in library file",*status);
		return(library_collection);}

	//Get number of templates (rows)
	long ntemplates;

	if (fits_get_num_rows(fptr,&ntemplates, status)){
	      EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
	      return(library_collection);}
	      
	// Allocate library structure
	library_collection->ntemplates = (int)ntemplates;
	library_collection->energies           		  = new double[ntemplates];
	library_collection->pulse_heights      		  = new double[ntemplates];
	library_collection->pulse_templates           = new PulseTemplate[ntemplates];
	library_collection->pulse_templates_filder    = new PulseTemplate[ntemplates];
	library_collection->pulse_templates_B0        = new PulseTemplate[ntemplates];
	library_collection->matched_filters           = new MatchedFilter[ntemplates];
	library_collection->matched_filters_B0        = new MatchedFilter[ntemplates];

	
	//Get energy (ENERGY), template (PULSE) and matched filter ("MF") column numbers
	char column_name[12];
	int template_colnum = 0;
	int mfilter_colnum = 0;
	int template_colnum_B0 = 0;
	int mfilter_colnum_B0 = 0;
	int energy_colnum = 0;
	int pheight_colnum = 0;
	strcpy(column_name,"PULSE");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &template_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for PULSE in library file",*status);
		*status=EPFAIL; return(library_collection);}
		
	strcpy(column_name,"PULSEB0");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &template_colnum_B0, status)){
		EP_PRINT_ERROR("Cannot get column number for PULSEB0 in library file",*status);
		*status=EPFAIL; return(library_collection);}
		
 	strcpy(column_name,"MF");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &mfilter_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for MF in library file",*status);
		*status=EPFAIL; return(library_collection);}
		
	strcpy(column_name,"MFB0");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &mfilter_colnum_B0, status)){
		EP_PRINT_ERROR("Cannot get column number for MFB0 in library file",*status);
		*status=EPFAIL; return(library_collection);}
		
	strcpy(column_name,"ENERGY");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &energy_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for ENERGY in library file",*status);
		*status=EPFAIL; return(library_collection);}
		
	strcpy(column_name,"PHEIGHT");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &pheight_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for PHEIGHT in library file",*status);
		*status=EPFAIL; return(library_collection);}

	//Get template duration
	int naxis = 0;
	long naxes = 0;

	if (fits_read_tdim(fptr, template_colnum, 1, &naxis, &naxes, status)){
	    	EP_PRINT_ERROR("Cannot read dim of column PULSE",*status);
		*status=EPFAIL; return(library_collection);
	}
	int template_duration = naxes;
	
	//Get mfilter duration
	if (fits_read_tdim(fptr, mfilter_colnum, 1, &naxis, &naxes, status)){
	    	EP_PRINT_ERROR("Cannot read dim of column MF",*status);
		*status=EPFAIL; return(library_collection);
	}
	int mfilter_duration = naxes;

	//Check that TEMPLATE & MFILTER have the same duration
	if(template_duration != mfilter_duration){
		EP_PRINT_ERROR("Error: template and matched filter have different durations",EPFAIL);
		*status=EPFAIL; return(library_collection);
	}

	//Iterate over the templates and populate the LibraryCollection structure
	int anynul=0;

	//Read ENERGY column for PulseTemplates and MatchedFilters
	if(fits_read_col(fptr, TDOUBLE, energy_colnum, 1,1, ntemplates, NULL, library_collection->energies, &anynul,status)){
		EP_PRINT_ERROR("Cannot read column ENERGY",*status);
		*status=EPFAIL; return(library_collection);
	}
	//Read PHEIGHT column for PulseTemplates and MatchedFilters
	if(fits_read_col(fptr, TDOUBLE, pheight_colnum, 1,1, ntemplates, NULL, library_collection->pulse_heights, &anynul,status)){
		EP_PRINT_ERROR("Cannot read column PHEIGHT",*status);
		*status=EPFAIL; return(library_collection);
	}
	//Read Templates and matched filter columns and populate LibraryCollection structure
	//library_collection->pulse_templates->ptemplate    = new double[template_duration];
        //library_collection->pulse_templates_B0->ptemplate = new double[template_duration];
	//library_collection->matched_filters->mfilter      = new double[mfilter_duration];
	//library_collection->matched_filters_B0->mfilter   = new double[mfilter_duration];

	for (int it = 0 ; it < ntemplates ; it++){
		library_collection->pulse_templates[it].ptemplate    		= new double[template_duration];
		library_collection->pulse_templates_filder[it].ptemplate    = new double[template_duration];
		library_collection->pulse_templates_B0[it].ptemplate 		= new double[template_duration];
		library_collection->matched_filters[it].mfilter      		= new double[mfilter_duration];
		library_collection->matched_filters_B0[it].mfilter   		= new double[mfilter_duration];

		library_collection->pulse_templates[it].template_duration    		= template_duration;
		library_collection->pulse_templates_filder[it].template_duration    = -1;
		library_collection->pulse_templates_B0[it].template_duration 		= template_duration;
		library_collection->matched_filters[it].mfilter_duration     		= mfilter_duration;
		library_collection->matched_filters_B0[it].mfilter_duration  		= mfilter_duration;


		library_collection->pulse_templates[it].energy    		= library_collection->energies[it];
		library_collection->pulse_templates_filder[it].energy	= library_collection->energies[it];
		library_collection->pulse_templates_B0[it].energy 		= library_collection->energies[it];
		library_collection->matched_filters[it].energy    		= library_collection->energies[it];
		library_collection->matched_filters_B0[it].energy 		= library_collection->energies[it];

		//Actually read template columns and save them into the LibraryCollection structure
		if(fits_read_col(fptr, TDOUBLE, template_colnum, it+1,1, template_duration, 
			NULL, library_collection->pulse_templates[it].ptemplate, &anynul,status)){
		  	EP_PRINT_ERROR("Error: cannot read template from PULSE column",*status);
			*status=EPFAIL; return(library_collection);
		}
		if(fits_read_col(fptr, TDOUBLE, template_colnum_B0, it+1,1, template_duration, 
			NULL, library_collection->pulse_templates_B0[it].ptemplate, &anynul,status)){
			EP_PRINT_ERROR("Error: cannot read template from PULSEB0 column",*status);
			*status=EPFAIL; return(library_collection);
		}

		if(fits_read_col(fptr, TDOUBLE, mfilter_colnum, it+1,1, mfilter_duration, 
			NULL, library_collection->matched_filters[it].mfilter, &anynul,status)){
			EP_PRINT_ERROR("Error: cannot read matched filter from MF column",*status);
			*status=EPFAIL; return(library_collection);
		}
		
		if(fits_read_col(fptr, TDOUBLE, mfilter_colnum_B0, it+1,1, mfilter_duration, 
			NULL, library_collection->matched_filters_B0[it].mfilter, &anynul,status)){
			EP_PRINT_ERROR("Error: cannot read matched filter from MF column",EPFAIL);
			*status=EPFAIL; return(library_collection);
		}
	}

	if (fits_close_file(fptr, status)){
			EP_PRINT_ERROR("Error closing library file",*status);
			return(library_collection);}

	return(library_collection);
}

/** Create and retrieve a NoiseSpec from a file. */
NoiseSpec* getNoiseSpec(const char* const filename,int* const status){

	//Create NoiseSpec structure
	NoiseSpec* noise_spectrum = new NoiseSpec;
	noise_spectrum->noise_duration = 0;
	noise_spectrum->noisespec  = NULL;
	noise_spectrum->noisefreqs = NULL;
	
	//Open FITS file in READONLY mode
	fitsfile* fptr = NULL;
	if(fits_open_file(&fptr, filename, READONLY, status)){
	  	EP_PRINT_ERROR("Error opening noise file",*status);
		*status=EPFAIL; return(noise_spectrum);
	}

	//Move to the first hdu
	int extver = 1;
	char HDUname[12];
	strcpy(HDUname,"NOISEALL");
	if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status)){
	  	EP_PRINT_ERROR("Error moving to HDU NOISE in noise file",*status);
		*status=EPFAIL;return(noise_spectrum);
	}
	
	//Get number of rows
	long noise_duration;
	if (fits_get_num_rows(fptr,&noise_duration, status)){
		EP_PRINT_ERROR("Cannot get number of rows in noise file",*status);
		*status=EPFAIL;return(noise_spectrum);
	}
	noise_spectrum->noise_duration = noise_duration/2;
	
	//Allocate NoiseSpec structure
	noise_spectrum->noisespec  = new double[noise_duration];
	noise_spectrum->noisefreqs = new double[noise_duration];
	
	//Get noise spectrum (CSD), and noise frequencies (FREQ) column numbers
	char column_name[12];
	int CSD_colnum = 0;
	int FREQ_colnum = 0;
	int anynul=0;
	strcpy(column_name,"CSD");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &CSD_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for CSD in noise file",*status);
		*status=EPFAIL; return(noise_spectrum);
	}
	strcpy(column_name,"FREQ");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &FREQ_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for FREQ in noise file", *status);
		*status=EPFAIL; return(noise_spectrum);
	}	
	// Read column CSD and save it into the structure
	if(fits_read_col(fptr, TDOUBLE, CSD_colnum,1,1, noise_duration, NULL, noise_spectrum->noisespec, &anynul,status)){
		EP_PRINT_ERROR("Error: cannot read noise spectrum from CSD column",*status);
		*status=EPFAIL; return(noise_spectrum);
	}

	// Read column FREQ and save it into the structure
	if(fits_read_col(fptr, TDOUBLE, FREQ_colnum,1,1, noise_duration, NULL, noise_spectrum->noisefreqs, &anynul, status)){
		EP_PRINT_ERROR("Error: cannot read noise frequencies from FREQ column",*status);
		*status=EPFAIL; return(noise_spectrum);
	}

	return(noise_spectrum);
}

/***** SECTION CX ************************************************************
* runEnergyCalib: This function...
****************************************************************************/
extern "C" void runEnergyCalib(ReconstructInitSIRENA* reconstruct_init, PulsesCollection* pulsesAll_runEnergyCalib, PulsesCollection* pulsesAll2_runEnergyCalib, double *b_cF_runEnergyCalib, double *c_cF_runEnergyCalib)
{
	const char * create= "runEnergyCalib v.1.0.0";	//Set "CREATOR" keyword of output FITS file

	string message="";
	int status = EPOK;

	long nx, ny;
	double E0x, E0y;
	gsl_vector *xi;
	gsl_vector *yi;

	if (loadUCEnergies(reconstruct_init, pulsesAll_runEnergyCalib, &nx, &xi, &E0x))
	{
		message = "Cannot run loadUCEnergies in calibration mode in runEenergyCalib";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	if (loadUCEnergies(reconstruct_init, pulsesAll2_runEnergyCalib, &ny, &yi, &E0y))
	{
		message = "Cannot run loadUCEnergies in calibration mode in runEenergyCalib";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	if (calculus_bc(reconstruct_init->calibLQ,nx,xi,E0x,ny,yi,E0y, b_cF_runEnergyCalib, c_cF_runEnergyCalib))
	{
		message = "Cannot run calculus_bc in calibration mode in runEenergyCalib";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	gsl_vector_free(xi);
	gsl_vector_free(yi);

	return;
}
/*xxxx end of SECTION C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
