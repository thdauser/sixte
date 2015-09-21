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
		char* const library_file, char* const event_file, double tauFall,int pulse_length, double scaleFactor, double samplesUp,
		double nSgms, int crtLib, int lastELibrary, int mode, double LrsT, double LbT, double baseline, char* const noise_file,
		char* pixel_type, char* filter_domain, char* filter_method, char* energy_method, char* oflength_strategy, int oflength,
		int calibLQ,  double b_cF, double c_cF, double monoenergy, int interm, char* const detectFile, char* const filterFile, char* const record_file2,
		double monoenergy2, char clobber, int maxPulsesPerRecord, int tstartPulse1, int tstartPulse2, int tstartPulse3, int* const status){

	// Load LibraryCollection structure if library file exists
	int exists=0;
	if(fits_file_exists(library_file, &exists, status)){
		EP_PRINT_ERROR("Error checking if library file exists",*status);
		return;
	}
	if(exists){
		reconstruct_init->library_collection = getLibraryCollection(library_file, mode, crtLib, energy_method, status);
		if(*status){
		      EP_PRINT_ERROR((char*)"Error in getLibraryCollection",EPFAIL); 
		      *status=EPFAIL;return;
		}
		if (pulse_length > reconstruct_init->library_collection->pulse_templates[0].template_duration)
		{
			EP_PRINT_ERROR("Templates length in the library file must be at least as the pulse length",EPFAIL);
			*status=EPFAIL;return;
		}
	}else if(!exists && crtLib==0){
		EP_PRINT_ERROR((char*)"Error accessing library file: it does not exists ",EPFAIL); 
		*status=EPFAIL;return;
	}

	//Load NoiseSpec structure
	reconstruct_init->noise_spectrum = NULL;
	if (crtLib == 0)
	{
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
	strncpy(reconstruct_init->record_file,record_file,255);
	reconstruct_init->record_file[255]='\0';
	strncpy(reconstruct_init->library_file,library_file,255);
	reconstruct_init->library_file[255]='\0';
	strncpy(reconstruct_init->noise_file,noise_file,255);
	reconstruct_init->noise_file[255]='\0';
	strncpy(reconstruct_init->event_file,event_file,255);
	reconstruct_init->event_file[255]='\0';
	strncpy(reconstruct_init->detectFile,detectFile,255);
	reconstruct_init->detectFile[255]='\0';
	strncpy(reconstruct_init->filterFile,filterFile,255);
	reconstruct_init->filterFile[255]='\0';
	strncpy(reconstruct_init->record_file2,record_file2,255);
	reconstruct_init->record_file2[255]='\0';
	reconstruct_init->threshold 	= 0.0;
	reconstruct_init->pulse_length 	= pulse_length;
	reconstruct_init->tauFall       = tauFall;
	reconstruct_init->scaleFactor  	= scaleFactor;
	reconstruct_init->samplesUp    	= samplesUp;
	reconstruct_init->nSgms        	= nSgms;
	reconstruct_init->crtLib       	= crtLib;
	reconstruct_init->lastELibrary  = lastELibrary;
	reconstruct_init->mode		= mode;
	reconstruct_init->LrsT		= LrsT;
	reconstruct_init->LbT		= LbT;
	reconstruct_init->baseline		= baseline;
	reconstruct_init->monoenergy 	= monoenergy;
	strcpy(reconstruct_init->PixelType,pixel_type);
	strcpy(reconstruct_init->FilterDomain,filter_domain);
	strcpy(reconstruct_init->FilterMethod,filter_method);
	strcpy(reconstruct_init->EnergyMethod,energy_method);
	strcpy(reconstruct_init->OFStrategy,oflength_strategy);
	if ((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0))		{strcpy(reconstruct_init->FilterMethod,"F0");}
	reconstruct_init->OFLength      = oflength;
	reconstruct_init->calibLQ       = calibLQ;
	reconstruct_init->b_cF          = b_cF;
	reconstruct_init->c_cF          = c_cF;
	reconstruct_init->intermediate  = interm;
	reconstruct_init->monoenergy2 	= monoenergy2;
	reconstruct_init->tstartPulse1 	= tstartPulse1;
	reconstruct_init->tstartPulse2 	= tstartPulse2;
	reconstruct_init->tstartPulse3 	= tstartPulse3;
	if(0!=clobber){
	    reconstruct_init->clobber = 1;
	}else{
	    reconstruct_init->clobber = 0;
	} 
	reconstruct_init->maxPulsesPerRecord = maxPulsesPerRecord;
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
		//cout<<"Acaba runDetect"<<endl;

		if(pulsesInRecord->ndetpulses == 0) // No pulses found in record
		{
			delete(pulsesAllAux);
		    delete(pulsesInRecord);
		    return;
		}
		
		if (reconstruct_init->crtLib == 0)
		{
			// filter pulses and calculates energy
			runFilter(record, nRecord, lastRecord, &reconstruct_init, *pulsesAll, &pulsesInRecord, optimalFilter);
			//cout<<"Acaba runFilter"<<endl;

			if (reconstruct_init->mode == 1)
			{
				// calibrated energy of pulses
				runEnergy(&reconstruct_init, &pulsesInRecord);
			}
			//cout<<"Acaba runEnergy"<<endl;
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
			if(event_list->energies != NULL) delete [] event_list->energies;
			if(event_list->grades1 != NULL) delete [] event_list->grades1;
			if(event_list->grades2 != NULL) delete [] event_list->grades2;
			if(event_list->pulse_heights != NULL) delete [] event_list->pulse_heights;
			if(event_list->ph_ids != NULL) delete [] event_list->ph_ids;
		
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
		
		//cout<<"pulsesAll: "<<(*pulsesAll)->ndetpulses<<endl;
		//cout<<"pulsesInRecord: "<<pulsesInRecord->ndetpulses<<endl;
		//cout<<pulsesInRecord->ndetpulses<<endl;

		// Free & Fill TesEventList structure
		
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
		delete pulsesAllAux;
		delete pulsesInRecord;

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
	strcpy(reconstruct_init->noise_file,"");
	strcpy(reconstruct_init->event_file,"");
	reconstruct_init->threshold=0.;
	reconstruct_init->pulse_length=0;	
	reconstruct_init->tauFall=0.;
	reconstruct_init->scaleFactor=0.;
	reconstruct_init->samplesUp=0.;
	reconstruct_init->nSgms=0;
	reconstruct_init->crtLib=0;
	reconstruct_init->lastELibrary=0;
	reconstruct_init->mode=0;
	reconstruct_init->monoenergy = 0;
	reconstruct_init->LrsT = 0;
	reconstruct_init->LbT = 0;
	strcpy(reconstruct_init->PixelType,"");
	strcpy(reconstruct_init->FilterDomain,"");
	strcpy(reconstruct_init->FilterMethod,"");
	strcpy(reconstruct_init->EnergyMethod,"");
	strcpy(reconstruct_init->OFStrategy,"");
	reconstruct_init->OFLength = 0;
	reconstruct_init->calibLQ=0;
	reconstruct_init->b_cF=0.;
	reconstruct_init->c_cF=0.;
	reconstruct_init->clobber=0;
	reconstruct_init->tstartPulse1=0;
	reconstruct_init->tstartPulse2=0;
	reconstruct_init->tstartPulse3=0;
	
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
LibraryCollection* getLibraryCollection(const char* const filename,int mode, int crtLib,char *energy_method,int* const status){

	//Create LibraryCollection structure
	LibraryCollection* library_collection = new LibraryCollection;
	library_collection->ntemplates=0;
	library_collection->pulse_templates=NULL;
	library_collection->pulse_templates_B0=NULL;
	library_collection->pulse_templates_filder=NULL;
	library_collection->matched_filters=NULL;
	library_collection->matched_filters_B0=NULL;
	library_collection->V=NULL;
	library_collection->W=NULL;
	library_collection->T=NULL;
	library_collection->t=NULL;
	library_collection->X=NULL;
	library_collection->Y=NULL;
	library_collection->Z=NULL;
	library_collection->r=NULL;
	library_collection->PAB=NULL;
	library_collection->DAB=NULL;
	
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
	library_collection->ntemplates                = (int)ntemplates;
	library_collection->energies           		  = gsl_vector_alloc(ntemplates);
	library_collection->pulse_heights      		  = gsl_vector_alloc(ntemplates);
	library_collection->maxDERs      		      = gsl_vector_alloc(ntemplates);
	library_collection->pulse_templates           = new PulseTemplate[ntemplates];
	library_collection->pulse_templates_filder    = new PulseTemplate[ntemplates];
	library_collection->pulse_templates_B0        = new PulseTemplate[ntemplates];
	library_collection->matched_filters           = new MatchedFilter[ntemplates];
	library_collection->matched_filters_B0        = new MatchedFilter[ntemplates];
	
	//Get energy (ENERGY), template (PULSE) and matched filter ("MF") column numbers
	char column_name[12];
	int template_colnum = 0;
	int mfilter_colnum = 0;

	strcpy(column_name,"PULSE");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &template_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for PULSE in library file",*status);
		*status=EPFAIL; return(library_collection);}
		
	/*strcpy(column_name,"PULSEB0");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &template_colnum_B0, status)){
		EP_PRINT_ERROR("Cannot get column number for PULSEB0 in library file",*status);
		*status=EPFAIL; return(library_collection);}
	*/	
 	strcpy(column_name,"MF");
	if(fits_get_colnum(fptr, CASEINSEN,column_name, &mfilter_colnum, status)){
		EP_PRINT_ERROR("Cannot get column number for MF in library file",*status);
		*status=EPFAIL; return(library_collection);}

	//Get template duration
	int naxis = 0;
	long naxes = 0;

	if (fits_read_tdim(fptr, template_colnum, 1, &naxis, &naxes, status)){
	    	EP_PRINT_ERROR("Cannot read dim of column PULSE",*status);
		*status=EPFAIL; return(library_collection);
	}
	int template_duration = naxes;

	library_collection->V = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
	library_collection->W = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
	if (ntemplates > 1)
	{
		library_collection->T = gsl_matrix_alloc(ntemplates-1,template_duration);
		library_collection->t = gsl_vector_alloc(ntemplates-1);
		library_collection->X = gsl_matrix_alloc(ntemplates-1,template_duration*template_duration);
		library_collection->Y = gsl_matrix_alloc(ntemplates-1,template_duration);
		library_collection->Z = gsl_matrix_alloc(ntemplates-1,template_duration);
		library_collection->r = gsl_vector_alloc(ntemplates-1);
		library_collection->PAB = gsl_matrix_alloc(ntemplates-1,template_duration);
		library_collection->DAB = gsl_matrix_alloc(ntemplates-1,template_duration);
	}
	else
	{
		library_collection->T = gsl_matrix_alloc(1,template_duration);
		library_collection->t = gsl_vector_alloc(1);
		library_collection->X = gsl_matrix_alloc(1,template_duration*template_duration);
		library_collection->Y = gsl_matrix_alloc(1,template_duration);
		library_collection->Z = gsl_matrix_alloc(1,template_duration);
		library_collection->r = gsl_vector_alloc(1);
		library_collection->PAB = gsl_matrix_alloc(1,template_duration);
		library_collection->DAB = gsl_matrix_alloc(1,template_duration);
	}

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
	IOData obj;
	obj.inObject = fptr;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"LIBRARY");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.unit = new char [255];
	strcpy(obj.nameCol,"ENERGY");
	obj.type = TDOUBLE;
	obj.iniRow = 1;
	obj.endRow = ntemplates;
	if (readFitsSimple (obj,&library_collection->energies))
	{
	   EP_PRINT_ERROR("Cannot run readFitsSimple in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(library_collection);
	}

	//Read PHEIGHT column for PulseTemplates and MatchedFilters
	strcpy(obj.nameCol,"PHEIGHT");
	if (readFitsSimple (obj,&library_collection->pulse_heights))
	{
	   EP_PRINT_ERROR("Cannot run readFitsSimple in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(library_collection);
	}

	//Read Templates and matched filter columns and populate LibraryCollection structure
	gsl_matrix *matrixAux_PULSE = gsl_matrix_alloc(ntemplates,template_duration);
	strcpy(obj.nameCol,"PULSE");
	if (readFitsComplex (obj,&matrixAux_PULSE))
	{
	   EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(library_collection);
	}

	gsl_matrix *matrixAux_PULSEB0 = gsl_matrix_alloc(ntemplates,template_duration);
	strcpy(obj.nameCol,"PULSEB0");
	if (readFitsComplex (obj,&matrixAux_PULSEB0))
	{
	   EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(library_collection);
	}

	gsl_matrix *matrixAux_MF = gsl_matrix_alloc(ntemplates,mfilter_duration);
	strcpy(obj.nameCol,"MF");
	if (readFitsComplex (obj,&matrixAux_MF))
	{
	   EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(library_collection);
	}

	gsl_matrix *matrixAux_MFB0 = gsl_matrix_alloc(ntemplates,mfilter_duration);
	strcpy(obj.nameCol,"MFB0");
	if (readFitsComplex (obj,&matrixAux_MFB0))
	{
	   EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(library_collection);
	}

	gsl_matrix *matrixAux_V = NULL;
	gsl_vector *vectorAux_V = NULL;
	gsl_matrix *matrixAux_W = NULL;
	gsl_vector *vectorAux_W = NULL;
	gsl_matrix *matrixAux_T = NULL;
	gsl_vector *vectorAux_T = NULL;
	gsl_vector *vectorAux_t = NULL;
	gsl_matrix *matrixAux_X = NULL;
	gsl_vector *vectorAux_X = NULL;
	gsl_matrix *matrixAux_Y = NULL;
	gsl_vector *vectorAux_Y = NULL;
	gsl_matrix *matrixAux_Z = NULL;
	gsl_vector *vectorAux_Z = NULL;
	gsl_vector *vectorAux_r = NULL;
	gsl_matrix *matrixAux_PAB = NULL;
	gsl_vector *vectorAux_PAB = NULL;
	gsl_matrix *matrixAux_DAB = NULL;
	gsl_vector *vectorAux_DAB = NULL;

	if ((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0))
	{
		matrixAux_V = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
		vectorAux_V = gsl_vector_alloc(template_duration*template_duration);
		strcpy(obj.nameCol,"COVARM");
		if (readFitsComplex (obj,&matrixAux_V))
		{
			EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
			*status=EPFAIL; return(library_collection);
		}

		matrixAux_W = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
		vectorAux_W = gsl_vector_alloc(template_duration*template_duration);
		strcpy(obj.nameCol,"WEIGHTM");
		if (readFitsComplex (obj,&matrixAux_W))
		{
			EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
			*status=EPFAIL; return(library_collection);
		}

		if ((mode == 1) || (crtLib == 0))
		{
			int dim;
			if (ntemplates == 1)
			{
				obj.endRow = 1;
				dim = 1;
			}
			else
			{
				obj.endRow = ntemplates-1;
				dim = ntemplates-1;
			}

			if (strcmp(energy_method,"WEIGHT") == 0)
			{
				matrixAux_T = gsl_matrix_alloc(dim,template_duration);
				vectorAux_T = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"TV");
				if (readFitsComplex (obj,&matrixAux_T))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				vectorAux_t = gsl_vector_alloc(dim);
				strcpy(obj.nameCol,"tE");
				if (readFitsSimple (obj,&vectorAux_t))
				{
					EP_PRINT_ERROR("Cannot run readFitsSimplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				matrixAux_X = gsl_matrix_alloc(dim,template_duration*template_duration);
				vectorAux_X = gsl_vector_alloc(template_duration*template_duration);
				strcpy(obj.nameCol,"XM");
				if (readFitsComplex (obj,&matrixAux_X))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				matrixAux_Y = gsl_matrix_alloc(dim,template_duration);
				vectorAux_Y = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"YV");
				if (readFitsComplex (obj,&matrixAux_Y))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				matrixAux_Z = gsl_matrix_alloc(dim,template_duration);
				vectorAux_Z = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"ZV");
				if (readFitsComplex (obj,&matrixAux_Z))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				vectorAux_r = gsl_vector_alloc(dim);
				strcpy(obj.nameCol,"rE");
				if (readFitsSimple (obj,&vectorAux_r))
				{
					EP_PRINT_ERROR("Cannot run readFitsSimplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}
			}
			else if (strcmp(energy_method,"WEIGHTN") == 0)
			{
				matrixAux_PAB = gsl_matrix_alloc(dim,template_duration);
				vectorAux_PAB = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"PAB");
				if (readFitsComplex (obj,&matrixAux_PAB))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				matrixAux_DAB = gsl_matrix_alloc(dim,template_duration);
				vectorAux_DAB = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"DAB");
				if (readFitsComplex (obj,&matrixAux_DAB))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}
			}
		}
	}

	for (int it = 0 ; it < ntemplates ; it++)
	{
		library_collection->pulse_templates[it].ptemplate    		= gsl_vector_alloc(template_duration);
		library_collection->pulse_templates_filder[it].ptemplate    = gsl_vector_alloc(template_duration);
		library_collection->pulse_templates_B0[it].ptemplate 		= gsl_vector_alloc(template_duration);
		library_collection->matched_filters[it].mfilter      		= gsl_vector_alloc(mfilter_duration);
		library_collection->matched_filters_B0[it].mfilter   		= gsl_vector_alloc(mfilter_duration);

		library_collection->pulse_templates[it].template_duration    		= template_duration;
		library_collection->pulse_templates_filder[it].template_duration    = -1;
		library_collection->pulse_templates_B0[it].template_duration 		= template_duration;
		library_collection->matched_filters[it].mfilter_duration     		= mfilter_duration;
		library_collection->matched_filters_B0[it].mfilter_duration  		= mfilter_duration;

		library_collection->pulse_templates[it].energy    		= gsl_vector_get(library_collection->energies,it);
		library_collection->pulse_templates_filder[it].energy	= gsl_vector_get(library_collection->energies,it);
		library_collection->pulse_templates_B0[it].energy 		= gsl_vector_get(library_collection->energies,it);
		library_collection->matched_filters[it].energy    		= gsl_vector_get(library_collection->energies,it);
		library_collection->matched_filters_B0[it].energy 		= gsl_vector_get(library_collection->energies,it);

		gsl_matrix_get_row(library_collection->pulse_templates[it].ptemplate,matrixAux_PULSE,it);
		gsl_matrix_get_row(library_collection->pulse_templates_B0[it].ptemplate,matrixAux_PULSEB0,it);
		gsl_matrix_get_row(library_collection->matched_filters[it].mfilter,matrixAux_MF,it);
		gsl_matrix_get_row(library_collection->matched_filters_B0[it].mfilter,matrixAux_MFB0,it);

		if (strcmp(energy_method,"WEIGHT") == 0)
		{
			gsl_matrix_get_row(vectorAux_V,matrixAux_V,it);
			gsl_matrix_set_row(library_collection->V,it,vectorAux_V);

			gsl_matrix_get_row(vectorAux_W,matrixAux_W,it);
			gsl_matrix_set_row(library_collection->W,it,vectorAux_W);

			if (((mode == 1) || (crtLib == 0)) && (it < ntemplates-1))
			{
				gsl_matrix_get_row(vectorAux_T,matrixAux_T,it);
				gsl_matrix_set_row(library_collection->T,it,vectorAux_T);

				gsl_vector_set(library_collection->t,it,gsl_vector_get(vectorAux_t,it));

				gsl_matrix_get_row(vectorAux_X,matrixAux_X,it);
				gsl_matrix_set_row(library_collection->X,it,vectorAux_X);

				gsl_matrix_get_row(vectorAux_Y,matrixAux_Y,it);
				gsl_matrix_set_row(library_collection->Y,it,vectorAux_Y);

				gsl_matrix_get_row(vectorAux_Z,matrixAux_Z,it);
				gsl_matrix_set_row(library_collection->Z,it,vectorAux_Z);

				gsl_vector_set(library_collection->r,it,gsl_vector_get(vectorAux_r,it));
			}
		}
		else if (strcmp(energy_method,"WEIGHTN") == 0)
		{
			gsl_matrix_get_row(vectorAux_V,matrixAux_V,it);
			gsl_matrix_set_row(library_collection->V,it,vectorAux_V);

			gsl_matrix_get_row(vectorAux_W,matrixAux_W,it);
			gsl_matrix_set_row(library_collection->W,it,vectorAux_W);

			if (((mode == 1) || (crtLib == 0)) && (it < ntemplates-1))
			{
				gsl_matrix_get_row(vectorAux_PAB,matrixAux_PAB,it);
				gsl_matrix_set_row(library_collection->PAB,it,vectorAux_PAB);

				gsl_matrix_get_row(vectorAux_DAB,matrixAux_DAB,it);
				gsl_matrix_set_row(library_collection->DAB,it,vectorAux_DAB);
			}
		}
	}
	gsl_matrix_free(matrixAux_PULSE);
	gsl_matrix_free(matrixAux_PULSEB0);
	gsl_matrix_free(matrixAux_MF);
	gsl_matrix_free(matrixAux_MFB0);
	gsl_matrix_free(matrixAux_V);
	gsl_vector_free(vectorAux_V);
	gsl_matrix_free(matrixAux_W);
	gsl_vector_free(vectorAux_W);
	gsl_matrix_free(matrixAux_T);
	gsl_vector_free(vectorAux_T);
	gsl_vector_free(vectorAux_t);
	gsl_matrix_free(matrixAux_X);
	gsl_vector_free(vectorAux_X);
	gsl_matrix_free(matrixAux_Y);
	gsl_vector_free(vectorAux_Y);
	gsl_matrix_free(matrixAux_Z);
	gsl_vector_free(vectorAux_Z);
	gsl_vector_free(vectorAux_r);
	gsl_matrix_free(matrixAux_PAB);
	gsl_vector_free(vectorAux_PAB);
	gsl_vector_free(vectorAux_DAB);

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
	noise_spectrum->noisespec  = gsl_vector_alloc(noise_duration);
	noise_spectrum->noisefreqs = gsl_vector_alloc(noise_duration);
	
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
	IOData obj;
	obj.inObject = fptr;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"NOISEALL");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.unit = new char [255];
	strcpy(obj.nameCol,"CSD");
	obj.type = TDOUBLE;
	obj.iniRow = 1;
	obj.endRow = noise_duration;
	if (readFitsSimple (obj,&noise_spectrum->noisespec))
	{
	   EP_PRINT_ERROR("Cannot run readFitsSimple in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(noise_spectrum);
	}

	// Read column FREQ and save it into the structure
	strcpy(obj.nameCol,"FREQ");
	if (readFitsSimple (obj,&noise_spectrum->noisefreqs))
	{
	   EP_PRINT_ERROR("Cannot run readFitsSimple in integraSIRENA.cpp",*status);
	   *status=EPFAIL; return(noise_spectrum);
	}

	return(noise_spectrum);
	delete fptr;
}

/***** SECTION CX ************************************************************
* runEnergyCalib: This function...
****************************************************************************/
extern "C" void runEnergyCalib(ReconstructInitSIRENA* reconstruct_init, PulsesCollection* pulsesAll, PulsesCollection* pulsesAll2, double *b_cF, double *c_cF)
{
	const char * create= "runEnergyCalib v.1.0.0";	//Set "CREATOR" keyword of output FITS file

	string message="";
	int status = EPOK;

	long nx, ny;
	double E0x, E0y; // keV
	gsl_vector *xi;  // keV
	gsl_vector *yi;  // keV

	if (loadUCEnergies(reconstruct_init, pulsesAll, &nx, &xi, &E0x))
	{
		message = "Cannot run loadUCEnergies in calibration mode in runEenergyCalib";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	if (loadUCEnergies(reconstruct_init, pulsesAll2, &ny, &yi, &E0y))
	{
		message = "Cannot run loadUCEnergies in calibration mode in runEenergyCalib";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	if (calculus_bc(reconstruct_init->calibLQ,nx,xi,E0x,ny,yi,E0y, b_cF, c_cF))
	{
		message = "Cannot run calculus_bc in calibration mode in runEenergyCalib";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	//cout<<"b_cF: "<<*b_cF<<endl;
	//cout<<"c_cF: "<<*c_cF<<endl;

	gsl_vector_free(xi);
	gsl_vector_free(yi);

	return;
}
/*xxxx end of SECTION C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
