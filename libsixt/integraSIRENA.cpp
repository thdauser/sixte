/***********************************************************************
   This file is part of SIXTE/SIRENA software.

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

   Copyright 2014:  INTEGRASIRENA has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

/***********************************************************************
*                      INTEGRASIRENA
*
*  File:       integraSIRENA.cpp
*  Developers: Beatriz Cobo
* 	       cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

/******************************************************************************
DESCRIPTION:

The purpose of this package is the integration os SIRENA in SIXTE.

MAP OF SECTIONS IN THIS FILE:

 - 1. initializeReconstructionSIRENA
 - 2. reconstructRecordSIRENA
 - 3. newReconstructInitSIRENA
 - 4. freeReconstructInitSIRENA
 - 5. newPulsesCollection
 - 6. freePulsesCollection
 - 7. newOptimalFilterSIRENA
 - 8. freeOptimalFilterSIRENA
 - 9. getLibraryCollection
 - 10. getNoiseSpec

*******************************************************************************/ 

#include "integraSIRENA.h"
#include "genutils.h"
#include "tasksSIRENA.h"

/***** SECTION 1 ************************************************************
* initializeReconstructionSIRENA: This function initializes the structure ReconstructInitSIRENA with the variables required 
*                                 for SIRENA reconstruction. The values are taken from the input parameters.
* 
* - Load LibraryCollection structure if library file exists
* - Load NoiseSpec structure
* - Fill in reconstruct_init
* 
* Parameters:
* - reconstruct_init: Member of ReconstructInitSIRENA structure to initialize the reconstruction parameters (pointer and values)
* - record_file: Filename of input data file with records
* - fptr: FITS object with pointer to data file
* - library_file: File name of calibration library
* - event_file: File name of output events (with reconstructed energy)
* - pulse_length: Pulse length
* - scaleFactor: Detection scale factor for initial filtering
* - samplesUp: Number of samples for threshold trespassing
* - nSgms: Number of standard deviations in the kappa-clipping process for threshold estimation
* - mode: Calibration run (0) or energy reconstruction run (1)
* - LrsT: Running sum length for the RS raw energy estimation (seconds)
* - LbT: Baseline averaging length for the RS raw energy estimation (seconds)
* - noise_file: Noise file
* - filter_domain: Filtering Domain: Time(T) or Frequency(F)
* - filter_method: Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline)
* - energy_method: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED or PCA
* - lagsornot: Lags (1) or no lags (0)
* - ofiter: Iterate (1) or not iterate (0)
* - oflib: Work or not with a library with optimal filters (1/0)
* - ofinterp: Optimal Filter by using the Matched Filter or the DAB as matched filter (MF/DAB)
*             It has been fixed in 'tesreconstruction' as 'DAB'
* - oflength_strategy: Optimal Filter length Strategy: FREE, BASE2, BYGRADE or FIXED
* - oflength: Optimal Filter length (taken into account if :option:`OFStrategy`=FIXED)
* - monoenergy: Monochromatic energy of input file in eV (only for library creation)
* - interm: Write or not intermediate files (1/0)
* - detectFile: Intermediate detections file (if intermediate=1)
* - filterFile: Intermediate filters file (if intermediate=1)
* - clobber: Overwrite or not output files if exist (1/0)
* - maxPulsesPerRecord: Default size of the event list
* - SaturationValue: Saturation level of the ADC curves
* - tstartPulse1: Tstart (samples) of the first pulse (different from 0 if the tstartPulsei input parameters are going to be used)
* - tstartPulse2: Tstart (samples) of the second pulse
* - tstartPulse3: Tstart (samples) of the third pulse (if 0 => PAIRS, if not 0 => TRIOS)
* - energyPCA1: First energy (only for PCA) 
* - energyPCA2: Second energy (only for PCA)
* - XMLFile: File name of the XML input file with instrument definition
* - status: Input/output status
******************************************************************************/
extern "C" void initializeReconstructionSIRENA(ReconstructInitSIRENA* reconstruct_init, char* const record_file, fitsfile *fptr,
		char* const library_file, char* const event_file, int pulse_length, double scaleFactor, double samplesUp,
		double nSgms, int mode, double LrsT, double LbT, char* const noise_file, char* filter_domain, char* filter_method, 
		char* energy_method, int lagsornot, int ofiter, char oflib, char *ofinterp,
		char* oflength_strategy, int oflength,
		double monoenergy, int interm, char* const detectFile, char* const filterFile,
		char clobber, int maxPulsesPerRecord, double SaturationValue,
		int tstartPulse1, int tstartPulse2, int tstartPulse3, double energyPCA1, double energyPCA2, char * const XMLFile, int* const status)
{  
	gsl_set_error_handler_off();
	/*string message = "";
	char valERROR[256];*/
		
	/*gsl_vector *prueba;
	if ((prueba = gsl_vector_alloc(0)) == 0)
	{
		sprintf(valERROR,"%d",__LINE__-2);
		string str(valERROR);
	        message = "Allocating with <= 0 size in line " + str + " (" + __FILE__ + ")";
		EP_EXIT_ERROR(message,EPFAIL);
	}*/
	/*gsl_vector *prueba;
	prueba = gsl_vector_alloc(2);
	if (gsl_vector_set(prueba,-1,1) != 0)	// No puede hacerse pq gsl_vector_set devuelve void, no un entero
	{
	  cout<<"Error"<<endl;
	}*/
	/*gsl_vector *prueba = gsl_vector_alloc(2);
	gsl_vector_set_all(prueba,1);
	gsl_vector *copia = gsl_vector_alloc(2);
	cout<<gsl_vector_memcpy(copia,prueba)<<endl;*/

	// Load LibraryCollection structure if library file exists
	int exists=0;
	if (fits_file_exists(library_file, &exists, status))
	{
		EP_EXIT_ERROR("Error checking if library file exists",*status);
	}
	if (exists)
	{
		reconstruct_init->library_collection = getLibraryCollection(library_file, mode, filter_domain, pulse_length, energy_method, filter_method, oflib, &ofinterp,status);
		if (*status)
		{
			EP_EXIT_ERROR((char*)"Error in getLibraryCollection",EPFAIL); 
		}
	
		if (pulse_length > reconstruct_init->library_collection->pulse_templates[0].template_duration)
		{
			EP_EXIT_ERROR("Templates length in the library file must be at least as the pulse length",EPFAIL);
		}
	}
	else if (!exists && mode==1)
	{
		EP_EXIT_ERROR((char*)"Error accessing library file: it does not exists ",EPFAIL); 
	}
	
	// Load NoiseSpec structure
	reconstruct_init->noise_spectrum = NULL;
	if ((mode == 0) || 
		(((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0) 
		|| (strcmp(energy_method,"I2RFITTED") == 0)) && (mode == 1)) 
		|| ((mode == 1) && (strcmp(energy_method,"WEIGHT") == 0))
		|| ((mode == 1) && (strcmp(energy_method,"WEIGHTN") == 0))) 
	  
	{
		exists=0;
		if(fits_file_exists(noise_file, &exists, status))
		{
			EP_EXIT_ERROR("Error checking if noise file exists",*status);
		}
		if (exists != 1)
		{
			EP_EXIT_ERROR("The necessary noise file does not exist",EPFAIL);
		}
		reconstruct_init->noise_spectrum = getNoiseSpec(noise_file, mode, energy_method, filter_method, status);
		if (*status)
		{
			EP_EXIT_ERROR((char*)"Error in getNoiseSpec",EPFAIL);
		}
	}

	// Fill in reconstruct_init
	strncpy(reconstruct_init->record_file,record_file,255);
	reconstruct_init->record_file[255]='\0';
	reconstruct_init->record_file_fptr = fptr;
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
	reconstruct_init->threshold 	= 0.0;
	reconstruct_init->pulse_length 	= pulse_length;
	reconstruct_init->scaleFactor  	= scaleFactor;
	reconstruct_init->samplesUp    	= samplesUp;
	reconstruct_init->nSgms        	= nSgms;
	reconstruct_init->mode		= mode;
	reconstruct_init->LrsT		= LrsT;
	reconstruct_init->LbT		= LbT;
	reconstruct_init->monoenergy 	= monoenergy;
	strcpy(reconstruct_init->FilterDomain,filter_domain);
	strcpy(reconstruct_init->FilterMethod,filter_method);
	strcpy(reconstruct_init->EnergyMethod,energy_method);
	reconstruct_init->LagsOrNot = lagsornot;
	reconstruct_init->OFIter = ofiter;
	if (0 != oflib)	reconstruct_init->OFLib = 1;
	else		reconstruct_init->OFLib = 0;
	strcpy(reconstruct_init->OFInterp,ofinterp);
	strcpy(reconstruct_init->OFStrategy,oflength_strategy);
	//if ((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0)) strcpy(reconstruct_init->FilterMethod,"F0");
	reconstruct_init->OFLength      = oflength;
	reconstruct_init->intermediate  = interm;
	reconstruct_init->SaturationValue  = SaturationValue;
	if (tstartPulse1 != 0)
	{
		if (strcmp(energy_method,"I2RALL") == 0)	reconstruct_init->tstartPulse1 = tstartPulse1-1-1;	// Because of the derivative
		else 						reconstruct_init->tstartPulse1 = tstartPulse1-1;	// To be consistent in the GSL indexes which start from 0
	}
	else			reconstruct_init->tstartPulse1 = tstartPulse1;
	if (tstartPulse2 != 0)
	{
		if (strcmp(energy_method,"I2RALL") == 0)	reconstruct_init->tstartPulse2 = tstartPulse2-1-1;	// Because of the derivative
		else 						reconstruct_init->tstartPulse2 = tstartPulse2-1;	// To be consistent in the GSL indexes which start from 0
	}
	else			reconstruct_init->tstartPulse2 = tstartPulse2;
	if (tstartPulse3 != 0)
	{
		if (strcmp(energy_method,"I2RALL") == 0)	reconstruct_init->tstartPulse3 = tstartPulse3-1-1;	// Because of the derivative
		else 						reconstruct_init->tstartPulse3 = tstartPulse3-1;	// To be consistent in the GSL indexes which start from 0
	}
	else			reconstruct_init->tstartPulse3 = tstartPulse3;
	reconstruct_init->energyPCA1 	= energyPCA1;
	reconstruct_init->energyPCA2 	= energyPCA2;
	strncpy(reconstruct_init->XMLFile,XMLFile,255);
	reconstruct_init->XMLFile[255]='\0'; 
	if (0 != clobber)	reconstruct_init->clobber = 1;
	else		reconstruct_init->clobber = 0;
	reconstruct_init->maxPulsesPerRecord = maxPulsesPerRecord;
}
/*xxxx end of SECTION 1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 2 ************************************************************
* reconstructRecordSIRENA: This function is the main wrapper function to detect, grade and calculate energy of pulses in input records.
*
* - Inititalize structure PulsesCollection
* - Check consistency of some input parameters
* - Detect pulses in input record (runDetect()). Save information of detected pulses
*      - If PCA, pulses energies are already written in the 'pulsesAll' structures
* - If in RECONSTRUCTION (:option:`mode` = 1) and not PCA:
* 	- Filter record and calculate energy of pulses (runEnergy())
* - Populate output event list with pulses energies, arrival time and grading
*
* Parameters:
* - record: Instance of TesRecord structure that contains the input record
* - event_list: Instance of TesEventList structure that contains the information of the reconstructed pulses
* - reconstruct_init: Member of ReconstructInitSIRENA structure to initialize the reconstruction parameters (pointer and values)
* - lastRecord: If record being analyzed is the last one, lastRecord = 1. Otherwise it is equal to 0
* - nRecord: Input record number
* - pulsesAll: Member of PulsesCollection structure to successively store all the pulses used to create the library. 
*              Re-populated after each processed record.
* - optimalFilter: Optimal filters used in reconstruction
* - status:Input/output status
******************************************************************************/
extern "C" void reconstructRecordSIRENA(TesRecord* record, TesEventList* event_list, ReconstructInitSIRENA* reconstruct_init,  int lastRecord, int nRecord, PulsesCollection **pulsesAll, OptimalFilterSIRENA **optimalFilter, int* const status)
{
	// Inititalize structure PulsesCollection
	PulsesCollection* pulsesInRecord = new PulsesCollection;
	pulsesInRecord->ndetpulses = 0;
	pulsesInRecord->pulses_detected = NULL;
	
	PulsesCollection* pulsesAllAux = new PulsesCollection;
	pulsesAllAux->ndetpulses = 0;
	pulsesAllAux->pulses_detected = NULL;

	// Check consistency of some input parameters
	if (record->trigger_size <= 0)
	{
		EP_EXIT_ERROR("Record size is <= 0",EPFAIL);
	}
	if(reconstruct_init->pulse_length > record->trigger_size)
	{
		EP_EXIT_ERROR("Warning: pulse length is larger than record size. Pulse length set to maximum value (record size)",EPFAIL);
	}

	// Detect pulses in record
	runDetect(record, nRecord, lastRecord, *pulsesAll, &reconstruct_init, &pulsesInRecord);

	if(pulsesInRecord->ndetpulses == 0) // No pulses found in record
	{
		freePulsesCollection(pulsesAllAux);
		freePulsesCollection(pulsesInRecord);
	    
		return;
	}
		
	if ((reconstruct_init->mode == 1) && (strcmp(reconstruct_init->EnergyMethod,"PCA") != 0))
	{
		// Filter pulses and calculates energy
		runEnergy(record, &reconstruct_init, &pulsesInRecord, optimalFilter);
	}

	if (nRecord == 1)
	{
		(*pulsesAll)->ndetpulses = pulsesInRecord->ndetpulses;
		if ((*pulsesAll)->pulses_detected != NULL) delete [] (*pulsesAll)->pulses_detected;
		(*pulsesAll)->pulses_detected = new PulseDetected[pulsesInRecord->ndetpulses];
		for (int i=0;i<(*pulsesAll)->ndetpulses;i++)
		{
			(*pulsesAll)->pulses_detected[i] = pulsesInRecord->pulses_detected[i];
		}
	}
	else
	{
		if (event_list->energies != NULL) 	delete [] event_list->energies;
		if (event_list->grades1 != NULL) 	delete [] event_list->grades1;
		if (event_list->grades2 != NULL) 	delete [] event_list->grades2;
		if (event_list->pulse_heights != NULL) 	delete [] event_list->pulse_heights;
		if (event_list->ph_ids != NULL) 	delete [] event_list->ph_ids;
	
		pulsesAllAux->ndetpulses = (*pulsesAll)->ndetpulses;
		pulsesAllAux->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];

		for (int i=0;i<(*pulsesAll)->ndetpulses;i++)
		{
			pulsesAllAux->pulses_detected[i] = (*pulsesAll)->pulses_detected[i];
		}
		(*pulsesAll)->ndetpulses = (*pulsesAll)->ndetpulses + pulsesInRecord->ndetpulses;

		if ((*pulsesAll)->pulses_detected != NULL) delete [] (*pulsesAll)->pulses_detected; 
		(*pulsesAll)->pulses_detected = new PulseDetected[(*pulsesAll)->ndetpulses];
		// Save pulses detected in previous records
		for (int i=0;i<pulsesAllAux->ndetpulses;i++)
		{
			(*pulsesAll)->pulses_detected[i] = pulsesAllAux->pulses_detected[i];
		}
		// Save pulses detected in current record
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
	
	if (strcmp(reconstruct_init->EnergyMethod,"PCA") != 0)
	{
	   	for (int ip=0; ip<pulsesInRecord->ndetpulses; ip++)
		{	  			
			// '+1' in order to undo the '-1' in initializeReconstructionSIRENA
			event_list->event_indexes[ip] = 
			      (int)((pulsesInRecord->pulses_detected[ip].Tstart - record->time)/record->delta_t + 0.5) + 1;	// '+0.5' to nearest integer (neither 'floor' nor 'ceil') 

			if (reconstruct_init->mode == 1)
			{
				event_list->energies[ip] = pulsesInRecord->pulses_detected[ip].energy;
			}
			else if (reconstruct_init->mode == 0)
			{
				event_list->energies[ip] = 999.;
			}

			event_list->grades1[ip]  = pulsesInRecord->pulses_detected[ip].grade1;
			event_list->grades2[ip]  = pulsesInRecord->pulses_detected[ip].grade2;
			event_list->pulse_heights[ip]  = pulsesInRecord->pulses_detected[ip].pulse_height;
			event_list->ph_ids[ip]   = 0;
		}
	}
	else
	{
		if (lastRecord == 1)
		{
		        // Free & Fill TesEventList structure
			event_list->index = (*pulsesAll)->ndetpulses;
			event_list->event_indexes = new int[event_list->index];
			event_list->energies = new double[event_list->index];
			event_list->grades1  = new int[event_list->index];
			event_list->grades2  = new int[event_list->index];
			event_list->pulse_heights  = new double[event_list->index];
			event_list->ph_ids   = new long[event_list->index];
		
			for (int ip=0; ip<(*pulsesAll)->ndetpulses; ip++)
			{	
				// '+1' in order to undo the '-1' in initializeReconstructionSIRENA
				event_list->event_indexes[ip] = 
				      (int)(((*pulsesAll)->pulses_detected[ip].Tstart - record->time)/record->delta_t + 0.5) + 1;	// '+0.5' to nearest integer (neither 'floor' nor 'ceil')

				if (reconstruct_init->mode == 1)
				{
					event_list->energies[ip] = (*pulsesAll)->pulses_detected[ip].energy;
				}
				else if (reconstruct_init->mode == 0)
				{
					event_list->energies[ip] = 999.;
				}

				event_list->grades1[ip]  = (*pulsesAll)->pulses_detected[ip].grade1;
				event_list->grades2[ip]  = (*pulsesAll)->pulses_detected[ip].grade2;
				event_list->pulse_heights[ip]  = (*pulsesAll)->pulses_detected[ip].pulse_height;
				event_list->ph_ids[ip]   = 0;		    
			}
		}
	}
	
	delete(pulsesAllAux->pulses_detected);
	freePulsesCollection(pulsesAllAux);
	delete(pulsesInRecord->pulses_detected);
	freePulsesCollection(pulsesInRecord);

	return;
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* ReconstructInitSIRENA: Constructor. It returns a pointer to an empty ReconstructInitSIRENA data structure.
* 
* - Initialize pointers with NULL for SIRENA
* - Initialize values for SIRENA
* 
* Parameters:
* - status: Input/output status
******************************************************************************/
extern "C" ReconstructInitSIRENA* newReconstructInitSIRENA(int* const status)
{	
	ReconstructInitSIRENA* reconstruct_init = new ReconstructInitSIRENA;
	
	// Initialize pointers with NULL for SIRENA
	reconstruct_init->library_collection =NULL;
	reconstruct_init->noise_spectrum     =NULL;
	reconstruct_init->grading = NULL;
	
	// Initialize values for SIRENA
	strcpy(reconstruct_init->record_file,"");
	reconstruct_init->record_file_fptr = 0;
	strcpy(reconstruct_init->library_file,"");
	strcpy(reconstruct_init->noise_file,"");
	strcpy(reconstruct_init->event_file,"");
	reconstruct_init->threshold=0.;
	reconstruct_init->pulse_length=0;	
	reconstruct_init->scaleFactor=0.;
	reconstruct_init->samplesUp=0.;
	reconstruct_init->nSgms=0;
	reconstruct_init->mode=0;
	reconstruct_init->monoenergy = 0;
	reconstruct_init->LrsT = 0;
	reconstruct_init->LbT = 0;
	strcpy(reconstruct_init->FilterDomain,"");
	strcpy(reconstruct_init->FilterMethod,"");
	strcpy(reconstruct_init->EnergyMethod,"");
	reconstruct_init->LagsOrNot = 0;
	reconstruct_init->OFIter = 0;
	reconstruct_init->OFLib=0;
	strcpy(reconstruct_init->OFInterp,"");
	strcpy(reconstruct_init->OFStrategy,"");
	reconstruct_init->OFLength = 0;
	reconstruct_init->clobber=0;
	reconstruct_init->SaturationValue=0;
	reconstruct_init->tstartPulse1=0;
	reconstruct_init->tstartPulse2=0;
	reconstruct_init->tstartPulse3=0;
	reconstruct_init->energyPCA1 = 0;
	reconstruct_init->energyPCA2 = 0;
	strcpy(reconstruct_init->XMLFile,"");
	
	return(reconstruct_init);
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* freeReconstructInitSIRENA: Destructor of ReconstructInit structure
* 
* Parameters:
* - reconstruct_init: Member of ReconstructInitSIRENA structure to initialize the reconstruction parameters (pointer and values)
******************************************************************************/
extern "C" void freeReconstructInitSIRENA(ReconstructInitSIRENA* reconstruct_init)
{
	delete(reconstruct_init);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* newPulsesCollection: Constructor. It returns a pointer to an empty PulsesCollection data structure.
* 
* Parameters:
* - status: Input/output status
******************************************************************************/
extern "C" PulsesCollection* newPulsesCollection(int* const status)
{
	PulsesCollection* PulsesColl = new PulsesCollection;

	// Initialize pointers with NULL for SIRENA
	PulsesColl->pulses_detected =NULL;

	// Initialize values for SIRENA
	PulsesColl->ndetpulses=0;

	return(PulsesColl);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
 

/***** SECTION 6 ************************************************************
* freePulsesCollection: Destructor of PulsesCollection structure.
* 
* Parameters:
* - PulsesColl: Instance of PulsesCollection structure
******************************************************************************/
extern "C" void freePulsesCollection(PulsesCollection* PulsesColl)
{
	delete(PulsesColl);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* newOptimalFilterSIRENA: Constructor. It returns a pointer to an empty OptimalFilterSIRENA data structure.
* 
* Parameters:
* - status: Input/output status
******************************************************************************/
extern "C" OptimalFilterSIRENA* newOptimalFilterSIRENA(int* const status)
{
	OptimalFilterSIRENA* OFilterColl = new OptimalFilterSIRENA;

	// Initialize pointers with NULL for SIRENA
	OFilterColl->ofilter =NULL;

	// Initialize values for SIRENA
	OFilterColl->ofilter_duration=0;

	return(OFilterColl);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* freeOptimalFilterSIRENA: Destructor of OptimalFilterSIRENA structure
* 
* Parameters:
* - OFilterColl: Instance of OptimalFilterSIRENA structure
******************************************************************************/
extern "C" void freeOptimalFilterSIRENA(OptimalFilterSIRENA* OFilterColl)
{
	delete(OFilterColl);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* getLibraryCollection: This funtion creates and retrieves a LibraryCollection from a file.
* 
* - Create LibraryCollection structure
* - Open FITS file in READONLY mode (move to the first HDU) and get number of templates (rows)
* - Allocate library structure
* - Get PULSE and MF column numbers (depending on the different options)
* - Get template duration
* - Allocate library structure (cont.)
* - Get matched filter duration
* - Read different columns and populate the *LibraryCollection* structure
* - Added new code to handle the new HDUs FIXFILTF, FIXFILTT and PRECALWN
* - Free allocated GSL vectors and matrices
* 
* Parameters:
* - filename: File with library information
* - mode: Calibration run (0) or energy reconstruction run (1)
* - filter_domain: Time domain ('T') or Frequency domain ('F')
* - pulse_length: Pulse length
* - energy_method: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED or PCA
* - filter_method: Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline)
* - oflib: Work or not with a library with optimal filters (1/0)
* - ofinterp: Optimal Filter by using the Matched Filter or the DAB as matched filter (MF/DAB) 
* 	      It has been fixed in 'tesreconstruction' as 'DAB' (but it would be possible to work with 'MF')
* - status: Input/output status
******************************************************************************/
LibraryCollection* getLibraryCollection(const char* const filename, int mode, char* filter_domain, int pulse_length, char *energy_method, char *filter_method, char oflib, char **ofinterp, int* const status)
{  
	// Create LibraryCollection structure
	LibraryCollection* library_collection = new LibraryCollection;
	library_collection->ntemplates=0;
	library_collection->pulse_templates=NULL;
	library_collection->pulse_templates_B0=NULL;
	library_collection->pulse_templates_filder=NULL;
	library_collection->matched_filters=NULL;
	library_collection->matched_filters_B0=NULL;
	library_collection->optimal_filters=NULL;
	library_collection->V=NULL;
	library_collection->W=NULL;
	library_collection->WAB=NULL;
	library_collection->T=NULL;
	library_collection->t=NULL;
	library_collection->X=NULL;
	library_collection->Y=NULL;
	library_collection->Z=NULL;
	library_collection->r=NULL;
	library_collection->PAB=NULL;
	library_collection->DAB=NULL;
	library_collection->optimal_filtersab=NULL;
	library_collection->optimal_filtersFREQ=NULL;
	library_collection->optimal_filtersTIME=NULL;
	library_collection->optimal_filtersabTIME=NULL;
	library_collection->optimal_filtersabFREQ=NULL;
	
	// Open FITS file in READONLY mode
	fitsfile* fptr = NULL;
	if (fits_open_file(&fptr, filename, READONLY, status))
	{
		EP_PRINT_ERROR("Error opening library file",*status);
		return(library_collection);
	}

	// Move to the first HDU
	int extver = 1;
	char HDUname[12];
	strcpy(HDUname,"LIBRARY");
	if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
	{
		EP_PRINT_ERROR("Error moving to HDU LIBRARY in library file",*status);
		return(library_collection);
	}

	// Get number of templates (rows)
	long ntemplates;
	if (fits_get_num_rows(fptr,&ntemplates, status))
	{
	      EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
	      return(library_collection);
	}
	if (ntemplates == 0)	
	{
		EP_PRINT_ERROR("The library has no rows",EPFAIL); 
		return(library_collection);
	}
	   
	if ((mode == 1) && 
		(((strcmp(energy_method,"OPTFILT") == 0)|| (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0) 
		|| (strcmp(energy_method,"I2RFITTED") == 0))))
	{
		if (ntemplates == 1)
		{	
			if (strcmp(*ofinterp,"DAB") == 0)	strcpy(*ofinterp,"MF");
							
			EP_PRINT_ERROR("The library only has one row, so no interpolation is going to be done (no matter 'OFInterp')",-999); // Only a warning
		}
	}

	// Allocate library structure
	// It is not necessary to check the allocation because 'ntemplates' has been checked previously
	library_collection->ntemplates                = (int)ntemplates;
	library_collection->energies           	      = gsl_vector_alloc(ntemplates);
	library_collection->pulse_heights      	      = gsl_vector_alloc(ntemplates);
	library_collection->maxDERs      	      = gsl_vector_alloc(ntemplates);
	library_collection->samp1DERs      	      = gsl_vector_alloc(ntemplates);
	library_collection->pulse_templates           = new PulseTemplate[ntemplates];
	library_collection->pulse_templates_filder    = new PulseTemplate[ntemplates];
	library_collection->pulse_templates_B0        = new PulseTemplate[ntemplates];
	library_collection->matched_filters           = new MatchedFilter[ntemplates];
	library_collection->matched_filters_B0        = new MatchedFilter[ntemplates];
	library_collection->optimal_filters           = new OptimalFilterSIRENA[ntemplates];
	library_collection->optimal_filtersab         = new OptimalFilterSIRENA[ntemplates];
	library_collection->optimal_filtersFREQ       = new OptimalFilterSIRENA[ntemplates];
	library_collection->optimal_filtersTIME       = new OptimalFilterSIRENA[ntemplates];
	library_collection->optimal_filtersabFREQ     = new OptimalFilterSIRENA[ntemplates];
	library_collection->optimal_filtersabTIME     = new OptimalFilterSIRENA[ntemplates];
	
	// Get PULSE and MF column numbers (depending the different options)
	char column_name[12];
	int template_colnum = 0;
	int mfilter_colnum = 0;
	
	strcpy(column_name,"PULSE");
	if (fits_get_colnum(fptr, CASEINSEN,column_name, &template_colnum, status))
	{
		EP_PRINT_ERROR("Cannot get column number for PULSE in library file",*status);
		*status=EPFAIL; return(library_collection);
	}
		
	if ((mode == 0) || 
		(((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0)
		|| (strcmp(energy_method,"I2RFITTED") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"MF") == 0))) 
	{
		strcpy(column_name,"MF");
		if (fits_get_colnum(fptr, CASEINSEN,column_name, &mfilter_colnum, status))
		{
			EP_PRINT_ERROR("Cannot get column number for MF in library file",*status);
			*status=EPFAIL; return(library_collection);
		  
		}
	}

	// Get template duration
	int naxis = 0;
	long naxes = 0;
	if (fits_read_tdim(fptr, template_colnum, 1, &naxis, &naxes, status))
	{
	    	EP_PRINT_ERROR("Cannot read dim of column PULSE",*status);
		*status=EPFAIL; return(library_collection);
	}
	int template_duration = naxes;
	if (template_duration == 0)	
	{
		EP_PRINT_ERROR("PULSE column vectors length is 0",EPFAIL); 
		return(library_collection);
	}

	// Allocate library structure (cont.)
	// It is not necessary to check the allocation because 'ntemplates' has been checked previously
	library_collection->V = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
	library_collection->W = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
	if (ntemplates > 1)
	{
		library_collection->WAB = gsl_matrix_alloc(ntemplates-1,template_duration*template_duration);
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
		library_collection->WAB = gsl_matrix_alloc(1,template_duration*template_duration);
		library_collection->T = gsl_matrix_alloc(1,template_duration);
		library_collection->t = gsl_vector_alloc(1);
		library_collection->X = gsl_matrix_alloc(1,template_duration*template_duration);
		library_collection->Y = gsl_matrix_alloc(1,template_duration);
		library_collection->Z = gsl_matrix_alloc(1,template_duration);
		library_collection->r = gsl_vector_alloc(1);
		library_collection->PAB = gsl_matrix_alloc(1,template_duration);
		library_collection->DAB = gsl_matrix_alloc(1,template_duration);
	}

	// Get mfilter duration
	int mfilter_duration;
	if ((mode == 0) || 
		(((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0)
		|| (strcmp(energy_method,"I2RFITTED") == 0))&& (oflib == 0) && (strcmp(*ofinterp,"MF") == 0)))
	{
		if (fits_read_tdim(fptr, mfilter_colnum, 1, &naxis, &naxes, status))
		{
			EP_PRINT_ERROR("Cannot read dim of column MF",*status);
			*status=EPFAIL; return(library_collection);
		}
		mfilter_duration = naxes;
		if (mfilter_duration == 0) 	
		{
			EP_PRINT_ERROR("MF column vectors length is 0",EPFAIL); 
			return(library_collection);
		}
		
		//Check that TEMPLATE & MFILTER have the same duration
		if(template_duration != mfilter_duration)
		{
			EP_PRINT_ERROR("Error: template and matched filter have different durations",EPFAIL);
			*status=EPFAIL; return(library_collection);
		}
	}

	// Read different columns and populate the *LibraryCollection* structure
	int anynul=0;

	// Read ENERGY column 
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

	// Read PHEIGHT column
	strcpy(obj.nameCol,"PHEIGHT");
	if (readFitsSimple (obj,&library_collection->pulse_heights))
	{
		EP_PRINT_ERROR("Cannot run readFitsSimple in integraSIRENA.cpp",*status);
		*status=EPFAIL; return(library_collection);
	}

	// It is not necessary to check the allocation because 'ntemplates' and 'template_duration' have been checked previously
	gsl_matrix *matrixAux_PULSE = gsl_matrix_alloc(ntemplates,template_duration);
	strcpy(obj.nameCol,"PULSE");
	if (readFitsComplex (obj,&matrixAux_PULSE))
	{
		EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
		*status=EPFAIL; return(library_collection);
	}

	// It is not necessary to check the allocation because 'ntemplates' and 'template_duration' have been checked previously
	gsl_matrix *matrixAux_PULSEB0 = gsl_matrix_alloc(ntemplates,template_duration);
	strcpy(obj.nameCol,"PULSEB0");
	if (readFitsComplex (obj,&matrixAux_PULSEB0))
	{
		EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
		*status=EPFAIL; return(library_collection);
	}

	gsl_matrix *matrixAux_MF = NULL;
	gsl_matrix *matrixAux_MFB0 = NULL;
	if ((mode == 0) || 
		(((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0) 
		|| (strcmp(energy_method,"I2RFITTED") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"MF") == 0)))
	{
		if ((mode == 0) || (strcmp(filter_method,"F0") == 0))
		{
		        // It is not necessary to check the allocation because 'ntemplates' and 'mfilter_duration' have been checked previously
			matrixAux_MF = gsl_matrix_alloc(ntemplates,mfilter_duration);
			strcpy(obj.nameCol,"MF");
			if (readFitsComplex (obj,&matrixAux_MF))
			{
				EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
				*status=EPFAIL; return(library_collection);
			}
		}
		if ((mode == 0) || (strcmp(filter_method,"B0") == 0))
		{
			// It is not necessary to check the allocation because 'ntemplates' and 'mfilter_duration' have been checked previously
			matrixAux_MFB0 = gsl_matrix_alloc(ntemplates,mfilter_duration);
			strcpy(obj.nameCol,"MFB0");
			if (readFitsComplex (obj,&matrixAux_MFB0))
			{
				EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
				*status=EPFAIL; return(library_collection);
			}  
		}  
	}

	int ofilter_duration;
	
	gsl_matrix *matrixAux_WAB = NULL;
	gsl_vector *vectorAux_WAB = NULL;
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

	int dim;
	if (((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0)) || (mode == 0))
	{
		// It is not necessary to check the allocation because 'ntemplates' and 'template_duration' have been checked previously
		matrixAux_V = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
		vectorAux_V = gsl_vector_alloc(template_duration*template_duration);
		strcpy(obj.nameCol,"COVARM");
		if (readFitsComplex (obj,&matrixAux_V))
		{
			EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
			*status=EPFAIL; return(library_collection);
		}

		// It is not necessary to check the allocation because 'ntemplates' and 'ofilter_durationab' have been checked previously
		matrixAux_W = gsl_matrix_alloc(ntemplates,template_duration*template_duration);
		vectorAux_W = gsl_vector_alloc(template_duration*template_duration);
		strcpy(obj.nameCol,"WEIGHTM");
		if (readFitsComplex (obj,&matrixAux_W))
		{
			EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
			*status=EPFAIL; return(library_collection);
		}

		if ((mode == 1) || ((mode == 0) && (ntemplates >1)))
		{
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

			if (((mode == 1) && (strcmp(energy_method,"WEIGHT") == 0)) || ((mode == 0) && (ntemplates >1)))
			{
				// It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
				matrixAux_T = gsl_matrix_alloc(dim,template_duration);
				vectorAux_T = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"TV");
				if (readFitsComplex (obj,&matrixAux_T))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				// It is not necessary to check the allocation because dim > 0
				vectorAux_t = gsl_vector_alloc(dim);
				strcpy(obj.nameCol,"tE");
				if (readFitsSimple (obj,&vectorAux_t))
				{
					EP_PRINT_ERROR("Cannot run readFitsSimplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				// It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
				matrixAux_X = gsl_matrix_alloc(dim,template_duration*template_duration);
				vectorAux_X = gsl_vector_alloc(template_duration*template_duration);
				strcpy(obj.nameCol,"XM");
				if (readFitsComplex (obj,&matrixAux_X))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				// It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
				matrixAux_Y = gsl_matrix_alloc(dim,template_duration);
				vectorAux_Y = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"YV");
				if (readFitsComplex (obj,&matrixAux_Y))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				// It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
				matrixAux_Z = gsl_matrix_alloc(dim,template_duration);
				vectorAux_Z = gsl_vector_alloc(template_duration);
				strcpy(obj.nameCol,"ZV");
				if (readFitsComplex (obj,&matrixAux_Z))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}

				// It is not necessary to check the allocation because dim > 0
				vectorAux_r = gsl_vector_alloc(dim);
				strcpy(obj.nameCol,"rE");
				if (readFitsSimple (obj,&vectorAux_r))
				{
					EP_PRINT_ERROR("Cannot run readFitsSimplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}
			}
			
			if (((mode == 1) && (strcmp(energy_method,"WEIGHTN") == 0)) || ((mode == 0) && (ntemplates >1)))
			{
				// It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
				matrixAux_WAB = gsl_matrix_alloc(dim,template_duration*template_duration);
				vectorAux_WAB = gsl_vector_alloc(template_duration*template_duration);
				strcpy(obj.nameCol,"WAB");
				if (readFitsComplex (obj,&matrixAux_WAB))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}
			}
		}
	}
		
	if ((mode == 1) && ((strcmp(energy_method,"WEIGHTN") == 0) || (((strcmp(energy_method,"OPTFILT") == 0)|| (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0)
		 || (strcmp(energy_method,"I2RNOL") == 0) || (strcmp(energy_method,"I2RFITTED") == 0)) && (strcmp(*ofinterp,"DAB") == 0))) || ((mode == 0) && (ntemplates >1)))
	{
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
		// It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
		matrixAux_PAB = gsl_matrix_alloc(dim,template_duration);
		vectorAux_PAB = gsl_vector_alloc(template_duration);
		strcpy(obj.nameCol,"PAB");
		if (readFitsComplex (obj,&matrixAux_PAB))
		{
			EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
			*status=EPFAIL; return(library_collection);
		}

		// It is not necessary to check the allocation because dim > 0 and 'template_duration' has been checked previously
		matrixAux_DAB = gsl_matrix_alloc(dim,template_duration);
		vectorAux_DAB = gsl_vector_alloc(template_duration);
		strcpy(obj.nameCol,"DAB");
		if (readFitsComplex (obj,&matrixAux_DAB))
		{
			EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
			*status=EPFAIL; return(library_collection);
		}
	}

	for (int it = 0 ; it < ntemplates ; it++)
	{
		// It is not necessary to check the allocation because 'template_duration' has been checked previously
		library_collection->pulse_templates[it].ptemplate    		= gsl_vector_alloc(template_duration);
		library_collection->pulse_templates_filder[it].ptemplate    	= gsl_vector_alloc(template_duration);
		library_collection->pulse_templates_B0[it].ptemplate 		= gsl_vector_alloc(template_duration);
		library_collection->matched_filters[it].mfilter      		= gsl_vector_alloc(template_duration);
		library_collection->matched_filters_B0[it].mfilter   		= gsl_vector_alloc(template_duration);

		library_collection->pulse_templates[it].template_duration		= template_duration;
		library_collection->pulse_templates_filder[it].template_duration    	= -1;
		library_collection->pulse_templates_B0[it].template_duration 		= template_duration;
		library_collection->matched_filters[it].mfilter_duration     		= template_duration;
		library_collection->matched_filters_B0[it].mfilter_duration  		= template_duration;

		library_collection->pulse_templates[it].energy    	= gsl_vector_get(library_collection->energies,it);
		library_collection->pulse_templates_filder[it].energy	= gsl_vector_get(library_collection->energies,it);
		library_collection->pulse_templates_B0[it].energy 	= gsl_vector_get(library_collection->energies,it);
		library_collection->matched_filters[it].energy    	= gsl_vector_get(library_collection->energies,it);
		library_collection->matched_filters_B0[it].energy 	= gsl_vector_get(library_collection->energies,it);

		gsl_matrix_get_row(library_collection->pulse_templates[it].ptemplate,matrixAux_PULSE,it);
		gsl_matrix_get_row(library_collection->pulse_templates_B0[it].ptemplate,matrixAux_PULSEB0,it);
		if ((mode == 0) || 
			(((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0)
			|| (strcmp(energy_method,"I2RFITTED") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"MF") == 0)))
		{
			if ((mode == 0) || (strcmp(filter_method,"F0") == 0)) 
			{
				gsl_matrix_get_row(library_collection->matched_filters[it].mfilter,matrixAux_MF,it);
			}
			if ((mode == 0) || (strcmp(filter_method,"B0") == 0))
			{
				gsl_matrix_get_row(library_collection->matched_filters_B0[it].mfilter,matrixAux_MFB0,it);
			}
		}
		
		if (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0)
		        || (strcmp(energy_method,"I2RFITTED") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"DAB") == 0))
		{
			if ((mode == 1)  && (it < ntemplates-1))
			{
				gsl_matrix_get_row(vectorAux_DAB,matrixAux_DAB,it);
				gsl_vector_memcpy(library_collection->matched_filters[it].mfilter,vectorAux_DAB);
			}
		}
		
		if (((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0)
			//|| (strcmp(energy_method,"I2RFITTED") == 0)) && (oflib == 0) && (strcmp(*ofinterp,"DAB") == 0) && (mode == 1))
			|| (strcmp(energy_method,"I2RFITTED") == 0)) && (strcmp(*ofinterp,"DAB") == 0) && (mode == 1))
		{
			if ((mode == 1)  && (it < ntemplates-1))
			{
				gsl_matrix_get_row(vectorAux_PAB,matrixAux_PAB,it);
				gsl_matrix_set_row(library_collection->PAB,it,vectorAux_PAB);

				gsl_matrix_get_row(vectorAux_DAB,matrixAux_DAB,it);
				gsl_matrix_set_row(library_collection->DAB,it,vectorAux_DAB);
			}
		}
		
		if ((mode == 0) || ((strcmp(energy_method,"WEIGHT") == 0) || (strcmp(energy_method,"WEIGHTN") == 0)))
		{
			gsl_matrix_get_row(vectorAux_V,matrixAux_V,it);
			gsl_matrix_set_row(library_collection->V,it,vectorAux_V);

			gsl_matrix_get_row(vectorAux_W,matrixAux_W,it);
			gsl_matrix_set_row(library_collection->W,it,vectorAux_W);
		}
			
		if ((mode == 1) && (it < ntemplates-1) && (strcmp(energy_method,"WEIGHT") == 0) ||
			((mode == 0) && (ntemplates > 1) && (it < ntemplates-1)))
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
		
		if ((mode == 1)  && (it < ntemplates-1) &&(strcmp(energy_method,"WEIGHTN") == 0) ||
			((mode == 0) && (ntemplates > 1) && (it < ntemplates-1)))
		{
			gsl_matrix_get_row(vectorAux_WAB,matrixAux_WAB,it);
			gsl_matrix_set_row(library_collection->WAB,it,vectorAux_WAB);
			
			gsl_matrix_get_row(vectorAux_PAB,matrixAux_PAB,it);
			gsl_matrix_set_row(library_collection->PAB,it,vectorAux_PAB);
			
			gsl_matrix_get_row(vectorAux_DAB,matrixAux_DAB,it);
			gsl_matrix_set_row(library_collection->DAB,it,vectorAux_DAB);
		}
	}

	// Free allocated GSL vectors and matrices
	gsl_matrix_free(matrixAux_PULSE);
	gsl_matrix_free(matrixAux_PULSEB0);
	if (matrixAux_MF != NULL) gsl_matrix_free(matrixAux_MF);
	if (matrixAux_MFB0 != NULL) gsl_matrix_free(matrixAux_MFB0);
	if (matrixAux_V != NULL) gsl_matrix_free(matrixAux_V);
	if (vectorAux_V != NULL) gsl_vector_free(vectorAux_V);
	if (matrixAux_W != NULL) gsl_matrix_free(matrixAux_W);
	if (vectorAux_W != NULL) gsl_vector_free(vectorAux_W);
	if (matrixAux_WAB != NULL) gsl_matrix_free(matrixAux_WAB);
	if (vectorAux_WAB != NULL) gsl_vector_free(vectorAux_WAB);
	if (matrixAux_T != NULL) gsl_matrix_free(matrixAux_T);
	if (vectorAux_T != NULL) gsl_vector_free(vectorAux_T);
	if (vectorAux_t != NULL) gsl_vector_free(vectorAux_t);
	if (matrixAux_X != NULL) gsl_matrix_free(matrixAux_X);
	if (vectorAux_X != NULL) gsl_vector_free(vectorAux_X);
	if (matrixAux_Y != NULL) gsl_matrix_free(matrixAux_Y);
	if (vectorAux_Y != NULL) gsl_vector_free(vectorAux_Y);
	if (matrixAux_Z != NULL) gsl_matrix_free(matrixAux_Z);
	if (vectorAux_Z != NULL) gsl_vector_free(vectorAux_Z);
	if (vectorAux_r != NULL) gsl_vector_free(vectorAux_r);
	if (matrixAux_PAB != NULL) gsl_matrix_free(matrixAux_PAB);
	if (vectorAux_PAB != NULL) gsl_vector_free(vectorAux_PAB);
	if (matrixAux_DAB != NULL) gsl_matrix_free(matrixAux_DAB);
	if (vectorAux_DAB != NULL) gsl_vector_free(vectorAux_DAB);
	
	if (mode == 0) 
	{
		// FIXFILTF HDU
		strcpy(HDUname,"FIXFILTF");
		if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
		{
			EP_PRINT_ERROR("Error moving to HDU FIXFILTF in library file",*status);
			return(library_collection);
		}

		// Get number of fixed optimal filters (columnss)
		int nOFs;
		if (fits_get_num_cols(fptr,&nOFs, status))
		{
			EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
			return(library_collection);
                }
		if (ntemplates == 1)	nOFs = nOFs-1;		// -1 because the ENERGYcolumn
		else 			nOFs = (nOFs-1)/2;	// /2 because the AB column
		if (nOFs == 0)	
		{
			EP_PRINT_ERROR("The library has no fixed optimal filters",EPFAIL); 
			return(library_collection);
		}
		library_collection->nfixedfilters = nOFs;
		
		int lengthALL_F = 0;
		int lengthALL_T = 0;
		for (int i=0;i<nOFs;i++)
		{
			lengthALL_F = lengthALL_F + pow(2,floor(log2(pulse_length))-i)*2;
			lengthALL_T = lengthALL_T + pow(2,floor(log2(pulse_length))-i);
		}
		
		gsl_matrix *matrixALL_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
		gsl_matrix *matrixALL_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
		gsl_matrix *matrixALLab_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
		gsl_matrix *matrixALLab_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
		gsl_matrix *matrixALL_PRCLWNx = gsl_matrix_alloc(ntemplates,lengthALL_F);
		
		char str_length[125];
		
		gsl_matrix *matrixAux_OFFx = NULL;
		gsl_matrix *matrixAuxab_OFFx = NULL;
		int index = 0;
		strcpy(obj.nameTable,"FIXFILTF");
		obj.iniRow = 1;
		obj.endRow = ntemplates;
		for (int i=0;i<nOFs;i++)
		{
			snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i)));
			strcpy(obj.nameCol,(string("OFF")+string(str_length)).c_str());
			matrixAux_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i)*2);
			if (readFitsComplex (obj,&matrixAux_OFFx))
			{
				EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
				*status=EPFAIL; return(library_collection);
			}
			for (int j=0;j<matrixAux_OFFx->size1;j++)
			{
				for (int k=0;k<matrixAux_OFFx->size2;k++)
				{
					gsl_matrix_set(matrixALL_OFFx,j,k+index,gsl_matrix_get(matrixAux_OFFx,j,k));
				}
			}
			
			if (ntemplates > 1)
			{
				strcpy(obj.nameCol,(string("OABF")+string(str_length)).c_str());
				matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i)*2);
				if (readFitsComplex (obj,&matrixAuxab_OFFx))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}
				for (int j=0;j<matrixAuxab_OFFx->size1;j++)
				{
					for (int k=0;k<matrixAuxab_OFFx->size2;k++)
					{
						gsl_matrix_set(matrixALLab_OFFx,j,k+index,gsl_matrix_get(matrixAuxab_OFFx,j,k));
					}
				}
			}
			
			index = index + pow(2,floor(log2(pulse_length))-i)*2;
			
			gsl_matrix_free(matrixAux_OFFx);
			gsl_matrix_free(matrixAuxab_OFFx);
		}
		
		// FIXFILTT HDU
		strcpy(HDUname,"FIXFILTT");
		if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
		{
			EP_PRINT_ERROR("Error moving to HDU FIXFILTT in library file",*status);
			return(library_collection);
		}
		
		gsl_matrix *matrixAux_OFTx = NULL;
		gsl_matrix *matrixAuxab_OFTx = NULL;
		index = 0;
		strcpy(obj.nameTable,"FIXFILTT");
		for (int i=0;i<nOFs;i++)
		{
			snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i)));
			strcpy(obj.nameCol,(string("OFT")+string(str_length)).c_str());
			matrixAux_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i));
			if (readFitsComplex (obj,&matrixAux_OFTx))
			{
				EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
				*status=EPFAIL; return(library_collection);
			}
			for (int j=0;j<matrixAux_OFTx->size1;j++)
			{
				for (int k=0;k<matrixAux_OFTx->size2;k++)
				{
					gsl_matrix_set(matrixALL_OFTx,j,k+index,gsl_matrix_get(matrixAux_OFTx,j,k));
				}
			}
			
			if (ntemplates > 1)
			{
				strcpy(obj.nameCol,(string("OABT")+string(str_length)).c_str());
				matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i));
				if (readFitsComplex (obj,&matrixAuxab_OFTx))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}
				for (int j=0;j<matrixAuxab_OFTx->size1;j++)
				{
					for (int k=0;k<matrixAuxab_OFTx->size2;k++)
					{
						gsl_matrix_set(matrixALLab_OFTx,j,k+index,gsl_matrix_get(matrixAuxab_OFTx,j,k));
					}
				}
			}
			
			index = index + pow(2,floor(log2(pulse_length))-i);
			
			gsl_matrix_free(matrixAux_OFTx);
			gsl_matrix_free(matrixAuxab_OFTx);
		}
		
		if (ntemplates > 1)
		{
			// PRECALWN HDU
			strcpy(HDUname,"PRECALWN");
			if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
			{
				EP_PRINT_ERROR("Error moving to HDU PRECALWN in library file",*status);
				return(library_collection);
			}
		
			gsl_matrix *matrixAux_PRCLWNx = NULL;
			index = 0;
			strcpy(obj.nameTable,"PRECALWN");
			for (int i=0;i<nOFs;i++)
			{
				snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(pulse_length))-i)));
				strcpy(obj.nameCol,(string("PRCL")+string(str_length)).c_str());
				matrixAux_PRCLWNx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(pulse_length))-i)*2);
				if (readFitsComplex (obj,&matrixAux_PRCLWNx))
				{
					EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
					*status=EPFAIL; return(library_collection);
				}
				for (int j=0;j<matrixAux_PRCLWNx->size1;j++)
				{
					for (int k=0;k<matrixAux_PRCLWNx->size2;k++)
					{
						gsl_matrix_set(matrixALL_PRCLWNx,j,k+index,gsl_matrix_get(matrixAux_PRCLWNx,j,k));
					}
				}
				
				index = index + pow(2,floor(log2(pulse_length))-i)*2;
				
				gsl_matrix_free(matrixAux_PRCLWNx);
			}
		}
			
		gsl_vector *vectorAux_PRCLWNx = gsl_vector_alloc(lengthALL_F);
		library_collection->PRECALWN = gsl_matrix_alloc(ntemplates,lengthALL_F);
		for (int it=0;it<ntemplates;it++)
		{
			library_collection->optimal_filtersFREQ[it].energy		= gsl_vector_get(library_collection->energies,it);
			library_collection->optimal_filtersFREQ[it].ofilter_duration	= lengthALL_F;
			library_collection->optimal_filtersFREQ[it].ofilter    	= gsl_vector_alloc(lengthALL_F);
			
			gsl_matrix_get_row(library_collection->optimal_filtersFREQ[it].ofilter,matrixALL_OFFx,it);
			
			library_collection->optimal_filtersTIME[it].energy		= gsl_vector_get(library_collection->energies,it);
			library_collection->optimal_filtersTIME[it].ofilter_duration 	= lengthALL_T;
			library_collection->optimal_filtersTIME[it].ofilter    	= gsl_vector_alloc(lengthALL_T);
			
			gsl_matrix_get_row(library_collection->optimal_filtersTIME[it].ofilter,matrixALL_OFTx,it);
			
			if (it < ntemplates-1)
			{
				library_collection->optimal_filtersabFREQ[it].energy		= gsl_vector_get(library_collection->energies,it);
				library_collection->optimal_filtersabFREQ[it].ofilter_duration  = lengthALL_F;
				library_collection->optimal_filtersabFREQ[it].ofilter    	= gsl_vector_alloc(lengthALL_F);
				
				gsl_matrix_get_row(library_collection->optimal_filtersabFREQ[it].ofilter,matrixALLab_OFFx,it);
				
				library_collection->optimal_filtersabTIME[it].energy		= gsl_vector_get(library_collection->energies,it);
				library_collection->optimal_filtersabTIME[it].ofilter_duration	= lengthALL_F;
				library_collection->optimal_filtersabTIME[it].ofilter    	= gsl_vector_alloc(lengthALL_T);
				
				gsl_matrix_get_row(library_collection->optimal_filtersabTIME[it].ofilter,matrixALLab_OFTx,it);
				
				//library_collection->PRECALWN = gsl_matrix_alloc(ntemplates,lengthALL_F);
				
				gsl_matrix_get_row(vectorAux_PRCLWNx,matrixALL_PRCLWNx,it);
				gsl_matrix_set_row(library_collection->PRECALWN,it,vectorAux_PRCLWNx);
			}
			
		}
		gsl_matrix_free(matrixALL_OFFx);
		gsl_matrix_free(matrixALL_OFTx);
		gsl_matrix_free(matrixALLab_OFFx);
		gsl_matrix_free(matrixALLab_OFTx);
		gsl_matrix_free(matrixALL_PRCLWNx);
		gsl_vector_free(vectorAux_PRCLWNx);
	}
	
	if ((mode == 1) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0) || (strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0)
	     || (strcmp(energy_method,"I2RFITTED") == 0)) && (oflib == 1))
	{
		char str_length[125];
		obj.iniRow = 1;
		obj.endRow = ntemplates;
		int index = 0;
			
		if (strcmp(filter_domain,"F") == 0)
		{
			// FIXFILTF HDU
			strcpy(HDUname,"FIXFILTF");
			if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
			{
				EP_PRINT_ERROR("Error moving to HDU FIXFILTF in library file",*status);
				return(library_collection);
			}

			// Get number of fixed optimal filters (columnss)
			int nOFs;
			if (fits_get_num_cols(fptr,&nOFs, status))
			{
				EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
				return(library_collection);
			}
			
			if (ntemplates == 1)
			{
				//nOFs = nOFs-1;		// -1 because the ENERGYcolumn
				int nOFs_aux;
				nOFs_aux = nOFs_aux-1;		// -1 because the ENERGYcolumn
				if (nOFs_aux == floor(log2(template_duration)))		nOFs = nOFs-1;		// -1 because the ENERGYcolumn
				else 							nOFs = (nOFs-1)/2;	// /2 because the AB column
			}
			else 								nOFs = (nOFs-1)/2;	// /2 because the AB column
			if (nOFs == 0)	
			{
				EP_PRINT_ERROR("The library has no fixed optimal filters",EPFAIL); 
				return(library_collection);
			}
			library_collection->nfixedfilters = nOFs;
		
			int lengthALL_F = 0;
			for (int i=0;i<nOFs;i++)
			{
				lengthALL_F = lengthALL_F + pow(2,floor(log2(template_duration))-i)*2;
			}
						
			strcpy(obj.nameTable,"FIXFILTF");
			
			if (strcmp(*ofinterp,"MF") == 0)
			{
				gsl_matrix *matrixALL_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
				gsl_matrix *matrixAux_OFFx = NULL;
				
				for (int i=0;i<nOFs;i++)
				{
					snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
					strcpy(obj.nameCol,(string("OFF")+string(str_length)).c_str());
					matrixAux_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i)*2);
					if (readFitsComplex (obj,&matrixAux_OFFx))
					{
						EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
						*status=EPFAIL; return(library_collection);
					}
					for (int j=0;j<matrixAux_OFFx->size1;j++)
					{
						for (int k=0;k<matrixAux_OFFx->size2;k++)
						{
							gsl_matrix_set(matrixALL_OFFx,j,k+index,gsl_matrix_get(matrixAux_OFFx,j,k));
						}
					}
					
					index = index + pow(2,floor(log2(template_duration))-i)*2;
					
					gsl_matrix_free(matrixAux_OFFx);
				}
				
				for (int it=0;it<ntemplates;it++)
				{
					library_collection->optimal_filters[it].energy			= gsl_vector_get(library_collection->energies,it);
					library_collection->optimal_filters[it].ofilter_duration 	= lengthALL_F;
					library_collection->optimal_filters[it].ofilter    		= gsl_vector_alloc(lengthALL_F);
					
					gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALL_OFFx,it);
				}
				
				gsl_matrix_free(matrixALL_OFFx);
			}
			else if (strcmp(*ofinterp,"DAB") == 0)
			{
				gsl_matrix *matrixALLab_OFFx = gsl_matrix_alloc(ntemplates,lengthALL_F);
				gsl_matrix *matrixAuxab_OFFx = NULL;
				
				for (int i=0;i<nOFs;i++)
				{
					snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
					strcpy(obj.nameCol,(string("OABF")+string(str_length)).c_str());
					matrixAuxab_OFFx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i)*2);
					if (readFitsComplex (obj,&matrixAuxab_OFFx))
					{
						EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
						*status=EPFAIL; return(library_collection);
					}
					for (int j=0;j<matrixAuxab_OFFx->size1;j++)
					{
						for (int k=0;k<matrixAuxab_OFFx->size2;k++)
						{
							gsl_matrix_set(matrixALLab_OFFx,j,k+index,gsl_matrix_get(matrixAuxab_OFFx,j,k));
						}
					}
					
					index = index + pow(2,floor(log2(template_duration))-i)*2;
					
					gsl_matrix_free(matrixAuxab_OFFx);
				}
				
				for (int it=0;it<ntemplates;it++)
				{
					library_collection->optimal_filters[it].energy			= gsl_vector_get(library_collection->energies,it);
					library_collection->optimal_filters[it].ofilter_duration 	= lengthALL_F;
					library_collection->optimal_filters[it].ofilter    		= gsl_vector_alloc(lengthALL_F);
					
					gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALLab_OFFx,it);
				}
				
				gsl_matrix_free(matrixALLab_OFFx);
			}
		}
		else if (strcmp(filter_domain,"T") == 0)
		{
			// FIXFILTT HDU
			strcpy(HDUname,"FIXFILTT");
			if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
			{
				EP_PRINT_ERROR("Error moving to HDU FIXFILTT in library file",*status);
				return(library_collection);
			}

			// Get number of fixed optimal filters (columnss)
			int nOFs;
			if (fits_get_num_cols(fptr,&nOFs, status))
			{
				EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
				return(library_collection);
			}
			if (ntemplates == 1)	
			{
				  //nOFs = nOFs-1;		// -1 because the ENERGYcolumn
				  int nOFs_aux;
				  nOFs_aux = nOFs_aux-1;		// -1 because the ENERGYcolumn
				  if (nOFs_aux == floor(log2(template_duration)))		nOFs = nOFs-1;		// -1 because the ENERGYcolumn
				  else 							nOFs = (nOFs-1)/2;	// /2 because the AB column
			}
			else 			nOFs = (nOFs-1)/2;	// /2 because the AB column
			if (nOFs == 0)	
			{
				EP_PRINT_ERROR("The library has no fixed optimal filters",EPFAIL); 
				return(library_collection);
			}
			library_collection->nfixedfilters = nOFs;
		  
			int lengthALL_T = 0;
			for (int i=0;i<nOFs;i++)
			{
				lengthALL_T = lengthALL_T + pow(2,floor(log2(template_duration))-i);
			}
		
			strcpy(obj.nameTable,"FIXFILTT");
			
			if (strcmp(*ofinterp,"MF") == 0)
			{
				gsl_matrix *matrixALL_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
				gsl_matrix *matrixAux_OFTx = NULL;
				
				for (int i=0;i<nOFs;i++)
				{
					snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
					strcpy(obj.nameCol,(string("OFT")+string(str_length)).c_str());
					matrixAux_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i));
					if (readFitsComplex (obj,&matrixAux_OFTx))
					{
						EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
						*status=EPFAIL; return(library_collection);
					}
					for (int j=0;j<matrixAux_OFTx->size1;j++)
					{
						for (int k=0;k<matrixAux_OFTx->size2;k++)
						{
							gsl_matrix_set(matrixALL_OFTx,j,k+index,gsl_matrix_get(matrixAux_OFTx,j,k));
						}
					}
					
					index = index + pow(2,floor(log2(template_duration))-i);
					
					gsl_matrix_free(matrixAux_OFTx);
				}
				
				for (int it=0;it<ntemplates;it++)
				{
					library_collection->optimal_filters[it].energy			= gsl_vector_get(library_collection->energies,it);
					library_collection->optimal_filters[it].ofilter_duration 	= lengthALL_T;
					library_collection->optimal_filters[it].ofilter    		= gsl_vector_alloc(lengthALL_T);
					
					gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALL_OFTx,it);
				}
				
				gsl_matrix_free(matrixALL_OFTx);
			}
			else if (strcmp(*ofinterp,"DAB") == 0)
			{
				gsl_matrix *matrixALLab_OFTx = gsl_matrix_alloc(ntemplates,lengthALL_T);
				gsl_matrix *matrixAuxab_OFTx = NULL;
				
				for (int i=0;i<nOFs;i++)
				{
					snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
					strcpy(obj.nameCol,(string("OABT")+string(str_length)).c_str());
					matrixAuxab_OFTx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i));
					if (readFitsComplex (obj,&matrixAuxab_OFTx))
					{
						EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
						*status=EPFAIL; return(library_collection);
					}
					for (int j=0;j<matrixAuxab_OFTx->size1;j++)
					{
						for (int k=0;k<matrixAuxab_OFTx->size2;k++)
						{
							gsl_matrix_set(matrixALLab_OFTx,j,k+index,gsl_matrix_get(matrixAuxab_OFTx,j,k));
						}
					}
										
					index = index + pow(2,floor(log2(template_duration))-i);
					
					gsl_matrix_free(matrixAuxab_OFTx);
				}	
				
				for (int it=0;it<ntemplates;it++)
				{
					library_collection->optimal_filters[it].energy			= gsl_vector_get(library_collection->energies,it);
					library_collection->optimal_filters[it].ofilter_duration 	= lengthALL_T;
					library_collection->optimal_filters[it].ofilter    		= gsl_vector_alloc(lengthALL_T);
					
					gsl_matrix_get_row(library_collection->optimal_filters[it].ofilter,matrixALLab_OFTx,it);
				}
				
				gsl_matrix_free(matrixALLab_OFTx);
			}	
		}  
	}
	
	if ((mode == 1) && (strcmp(energy_method,"WEIGHTN") == 0) && (ntemplates > 1))
	{
		char str_length[125];
		obj.iniRow = 1;
		obj.endRow = ntemplates;
		int index = 0;
		
		// PRECALWN HDU
		strcpy(HDUname,"PRECALWN");
		if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
		{
			EP_PRINT_ERROR("Error moving to HDU PRECALWN in library file",*status);
			return(library_collection);
		}

		// Get number of fixed optimal filters (columnss)
		int nOFs;
		if (fits_get_num_cols(fptr,&nOFs, status))
		{
			EP_PRINT_ERROR("Cannot get number of rows in library file",*status);
			return(library_collection);
		}
		nOFs = nOFs-1;	// -1 because the ENERGYcolumn
		if (nOFs == 0)	
		{
			EP_PRINT_ERROR("The library has no fixed optimal filters",EPFAIL); 
			return(library_collection);
		}
		library_collection->nfixedfilters = nOFs;
	
		int lengthALL_F = 0;
		for (int i=0;i<nOFs;i++)
		{
			lengthALL_F = lengthALL_F + pow(2,floor(log2(template_duration))-i)*2;
		}
					
		strcpy(obj.nameTable,"PRECALWN");
		
		gsl_matrix *matrixALL_PRCLWNx = gsl_matrix_alloc(ntemplates,lengthALL_F);
		gsl_matrix *matrixAux_PRCLWNx = NULL;
		for (int i=0;i<nOFs;i++)
		{
			snprintf(str_length,125,"%d",(int) (pow(2,floor(log2(template_duration))-i)));
			strcpy(obj.nameCol,(string("PRCL")+string(str_length)).c_str());
			matrixAux_PRCLWNx = gsl_matrix_alloc(ntemplates,pow(2,floor(log2(template_duration))-i)*2);
			if (readFitsComplex (obj,&matrixAux_PRCLWNx))
			{
				EP_PRINT_ERROR("Cannot run readFitsComplex in integraSIRENA.cpp",*status);
				*status=EPFAIL; return(library_collection);
			}
		
			for (int j=0;j<matrixAux_PRCLWNx->size1;j++)
			{
				for (int k=0;k<matrixAux_PRCLWNx->size2;k++)
				{
					gsl_matrix_set(matrixALL_PRCLWNx,j,k+index,gsl_matrix_get(matrixAux_PRCLWNx,j,k));
				}
			}
			
			index = index + pow(2,floor(log2(template_duration))-i)*2;
			
			gsl_matrix_free(matrixAux_PRCLWNx);
		}
		
		library_collection->PRECALWN = gsl_matrix_alloc(ntemplates, lengthALL_F);
		gsl_matrix_memcpy(library_collection->PRECALWN,matrixALL_PRCLWNx);
		gsl_matrix_free(matrixALL_PRCLWNx);
	}
			
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	if (fits_close_file(fptr, status))
	{
		EP_PRINT_ERROR("Error closing library file",*status);
		return(library_collection);
	}
				
	return(library_collection);
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ************************************************************
* getNoiseSpec: This function creates and retrieves a NoiseSpec from a file.
* 
* - Create *NoiseSpec* structure
* - Open FITS file, move to the NOISE and NOISEALL HDUs and get interesting keywords
* - Allocate *NoiseSpec* structure
* - Get noise spectrum (CSD), and noise frequencies (FREQ) column numbers
* - Read column CSD and save it into the structure
* - Read column FREQ and save it into the structure
* - Return noise spectrum
    
* Parameters:
* - filename: File name with noise
* - mode: Calibration run (0) or energy reconstruction run (1)
* - energy_method: Energy calculation Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RALL, I2RNOL, I2RFITTED or PCA
* - filter_method: Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline)
* - status: Input/output status
******************************************************************************/
NoiseSpec* getNoiseSpec(const char* const filename, int mode, char *energy_method, char *filter_method, int* const status)
{
	// Create NoiseSpec structure
	NoiseSpec* noise_spectrum = new NoiseSpec;
	noise_spectrum->noise_duration = 0;
	noise_spectrum->noisespec  = NULL;
	noise_spectrum->noisefreqs = NULL;
	
	// Open FITS file in READONLY mode
	fitsfile* fptr = NULL;
	if(fits_open_file(&fptr, filename, READONLY, status))
	{
	  	EP_PRINT_ERROR("Error opening noise file",*status);
		*status=EPFAIL; return(noise_spectrum);
	}
	
	// Move to the NOISE HDU
	int extver = 1;
	char HDUname[12];
	char keyname[10];
	strcpy(HDUname,"NOISE");
	if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
	{
	  	EP_PRINT_ERROR("Error moving to HDU NOISE in noise file",*status);
		*status=EPFAIL;return(noise_spectrum);
	}
	strcpy(keyname,"NOISESTD");
	if (fits_read_key(fptr,TDOUBLE,keyname, &noise_spectrum->noiseStd,NULL,status))
	{
		EP_PRINT_ERROR("Cannot read keyword NOISESTD",*status);
		*status=EPFAIL;return(noise_spectrum);
	}
	
	if ((mode == 0) ||
                ((mode == 1) && (strcmp(filter_method,"B0") == 0) && ((strcmp(energy_method,"OPTFILT") == 0) || (strcmp(energy_method,"I2R") == 0)))
                || ((mode == 1) && ((strcmp(energy_method,"I2RALL") == 0) || (strcmp(energy_method,"I2RNOL") == 0) || (strcmp(energy_method,"I2RFITTED") == 0)))
		|| ((mode == 1) && (strcmp(energy_method,"WEIGHT") == 0)))
	{
		strcpy(keyname,"BASELINE");
		
		if (fits_read_key(fptr,TDOUBLE,keyname, &noise_spectrum->baseline,NULL,status))
		{
			EP_PRINT_ERROR("Cannot read keyword BASELINE",*status);
			*status=EPFAIL;return(noise_spectrum);
		}
	}
	
	if (strcmp(energy_method,"WEIGHT") != 0)
	{
		// Move to the NOISEALL hdu
		strcpy(HDUname,"NOISEALL");
		if (fits_movnam_hdu(fptr, ANY_HDU,HDUname, extver, status))
		{
			EP_PRINT_ERROR("Error moving to HDU NOISEALL in noise file",*status);
			*status=EPFAIL;return(noise_spectrum);
		}
		
		// Get number of rows
		long noise_duration;
		if (fits_get_num_rows(fptr,&noise_duration, status))
		{
			EP_PRINT_ERROR("Cannot get number of rows in noise file",*status);
			*status=EPFAIL;return(noise_spectrum);
		}
		if (noise_duration == 0)	
		{	
			EP_PRINT_ERROR("Number of rows in noise file is 0",EPFAIL); 
			return(noise_spectrum);
		  
		}
		noise_spectrum->noise_duration = noise_duration/2;
		
		// Allocate NoiseSpec structure
		// It is not necessary to check the allocation because 'noise_duration' has been checked previously
		noise_spectrum->noisespec  = gsl_vector_alloc(noise_duration);
		noise_spectrum->noisefreqs = gsl_vector_alloc(noise_duration);
		
		// Get noise spectrum (CSD), and noise frequencies (FREQ) column numbers
		char column_name[12];
		int CSD_colnum = 0;
		int FREQ_colnum = 0;
		int anynul=0;
		strcpy(column_name,"CSD");
		if (fits_get_colnum(fptr, CASEINSEN,column_name, &CSD_colnum, status))
		{
			EP_PRINT_ERROR("Cannot get column number for CSD in noise file",*status);
			*status=EPFAIL; return(noise_spectrum);
		}
		strcpy(column_name,"FREQ");
		if (fits_get_colnum(fptr, CASEINSEN,column_name, &FREQ_colnum, status))
		{
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
		
		delete [] obj.nameTable;
		delete [] obj.nameCol;
		delete [] obj.unit;
	}

	if (fits_close_file(fptr, status))
	{
		EP_PRINT_ERROR("Error closing noise file",*status);
		return(noise_spectrum);
	  
	}
	
	// Return noise spectrum
	return(noise_spectrum);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/