#include "tasksSIRENA.h"
//#include <iomanip>
#include "assert.h"
#include <gsl/gsl_eigen.h>

/***** SECTION A ************************************************************
* runDetect: This function is responsible for the detection in SIRENA, record by record.
*            There are two run operation modes, CALIBRATION ('mode'=0) and the usual PRODUCTION ('mode'=1).
*            In CALIBRATION mode the purpose is the library creation.
*
* If first record and PRODUCTION => 'filderLibrary'
* If last record and CALIBRATION => 'calculateTemplate' and 'writeLibrary'
*                                   
* If 'reconstruct_init->intermediate'=1 => 'writeTestInfo' and 'writePulses'
* 
* If CALIBRATION => Find pulses by using 'findPulsesCAL'
* If PRODUCTION => Find pulses by 'InitialTriggering' and 'FindSecondaries'
*
* - Create library if it is necessary, CALIBRATION and 'lastRecord'=1 ('createLibrary')
* - Create intermediate output FITS file if 'reconstruct_init->intermediate'=1 ('createDetectFile')
* - Filter and derive the 'models' of the library (only for the first record in PRODUCTION mode) ('filderLibrary')
* - Store the input record in 'invector' ('loadRecord')
* - Convert I into R if 'EnergyMethod' = I2R or I2RBIS ('convertI2R')
* - Process each record ('proceRecord')
* 	- Low-pass filter and derive
* 	- Find pulses
* 	- Load the found pulses data in the input/output 'foundPulses' structure
* 	- Write test info in intermediate output FITS file if 'reconstruct_init->intermediate'=1 ('writeTestInfo')
* 	- Write pulses info in intermediate output FITS file if 'reconstruct_init->intermediate'=1 ('writePulses')
* - If last record in CALIBRATION mode:
* 	- 'calculateTemplate'
* 	- 'writeLibrary'
* - Close intermediate output FITS file if it is necessary
*
* Parameters:
* - record: Structure that contains the input record
* - lastRecord: Last record => 1, not last record => 0
* - pulsesAll: Structure containing the info of the found pulses in the previous records
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
* - pulsesInRecord: Structure containing the info of the found pulses in the input record
******************************************************************************/
void runDetect(TesRecord* record, int lastRecord, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
{
 	const char * create= "runDetect v.17.0.0";	// In order to set "CREATOR" keyword

	/*gsl_matrix *matrix = gsl_matrix_alloc(2,2);
	gsl_matrix_set(matrix,0,0,0.616555556);
	gsl_matrix_set(matrix,0,1,0.615444444);
	gsl_matrix_set(matrix,1,0,0.615444444);
	gsl_matrix_set(matrix,1,1,0.716555556);
	gsl_matrix *eigenvectors = gsl_matrix_alloc(2,2);
	gsl_vector *eigenvalues = gsl_vector_alloc(2);
	eigenVV(matrix,eigenvectors,eigenvalues);*/
	
	string message="";
	int status=EPOK;

	// Declare variables
	fitsfile *inLibObject = NULL;	// Object which contains information of the library FITS file
	bool appendToLibrary = false;	// Pulse templates library FITS file new (appendToLibrary=false) or not (appendToLibrary=true)

	fitsfile *dtcObject = NULL;	// Object which contains information of the intermediate FITS file ('dtc' comes from 'detectFile')
	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	int eventsz = record->trigger_size;
	double tstartRecord;
	gsl_vector *invector = gsl_vector_alloc(eventsz);	// Record
	
	if (((*reconstruct_init)->tstartPulse1 != 0) && (((*reconstruct_init)->tstartPulse1 > record->trigger_size) || ((*reconstruct_init)->tstartPulse2 > record->trigger_size) || ((*reconstruct_init)->tstartPulse3 > record->trigger_size)))
	{
		message = "tstartPulsx can not be greater than the record size";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	//assert(((*reconstruct_init)->tstartPulse1 <= record->trigger_size) && ((*reconstruct_init)->tstartPulse2 <= record->trigger_size) && ((*reconstruct_init)->tstartPulse3 <= record->trigger_size));
	
	// Create library if it is necessary
	if (((*reconstruct_init)->mode == 0) && (lastRecord == 1))
	{
		if (createLibrary(*reconstruct_init, &appendToLibrary, &inLibObject, create))
		{
			message = "Cannot run routine createLibrary to create pulses library";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Create intermediate output FITS file if 'reconstruct_init->intermediate'=1
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (createDetectFile(*reconstruct_init,1/record->delta_t,create, &dtcObject))
		{
			message = "Cannot create file " +  string((*reconstruct_init)->detectFile);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Filter and derive the 'models' of the library (only for the first record)
	if ((*reconstruct_init)->mode == 1)
	{
		if (filderLibrary(reconstruct_init,1/record->delta_t))
		{
			message = "Cannot run routine filderLibrary to filter & derive library if the 1st record";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Store the input record in 'invector'
	if (loadRecord(record, &tstartRecord, &invector))
	{
		message = "Cannot run routine loadRecord";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	eventsz = invector->size;	// Just in case the last record has been filled in with 0's => Re-allocate invector

	// Convert I into R if 'EnergyMethod' = I2R or I2RBIS
	gsl_vector *invectorWithoutConvert2R = gsl_vector_alloc(invector->size);
	gsl_vector_memcpy(invectorWithoutConvert2R,invector);
	if ((strcmp((*reconstruct_init)->EnergyMethod,"I2R") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"I2RBIS") == 0))
	{
		if (convertI2R(reconstruct_init, &record, &invector))
		{
			message = "Cannot run routine convertI2R";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Process each record
	if (procRecord(reconstruct_init, tstartRecord, 1/record->delta_t, dtcObject, invector, invectorWithoutConvert2R, *pulsesInRecord))
	{
		message = "Cannot run routine procRecord for record processing";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_vector_free(invectorWithoutConvert2R);
	
	if (((*reconstruct_init)->intermediate == 1) && (lastRecord == 1))
	{		
		// Write output keywords (their values have been previously checked)
		char keyname[10];
		
		char extname[10];
		strncpy(extname,"PULSES",9);
		extname[9]='\0';
		//strcpy(extname,"PULSES");
		if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU " + string(extname) +" in " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		long totalpulses;
		if (fits_get_num_rows(dtcObject,&totalpulses, &status))
		{
			message = "Cannot get number of rows in " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}

		int ttpls1 = (int) totalpulses;
		strcpy(keyname,"EVENTCNT");
		if (ttpls1 < 0)
		{
			message = "Legal values for EVENTCNT (PULSES) are integer numbers greater than or equal to 0";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		//if(fits_write_key(dtcObject,TINT,keyname,&ttpls1,NULL,&status))
		if(fits_write_key(dtcObject,TINT,keyname,&ttpls1,NULL,&status))
		{
			message = "Cannot write keyword " + string(keyname) +" in " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	if ((lastRecord == 1) && (*reconstruct_init)->mode == 0 && pulsesAll->ndetpulses>0)	// CALIBRATION mode => Calculate the pulse template by averaging some found pulses
	{
		gsl_vector *pulsetemplate = gsl_vector_alloc((*reconstruct_init)->pulse_length);
		double pulseheighttemplate = 0;
		gsl_matrix *weight = gsl_matrix_alloc((*reconstruct_init)->pulse_length,(*reconstruct_init)->pulse_length);
		gsl_matrix *covariance = gsl_matrix_alloc((*reconstruct_init)->pulse_length,(*reconstruct_init)->pulse_length);
		gsl_matrix_set_zero(weight);
		gsl_matrix_set_zero(covariance);

		if (calculateTemplate (*reconstruct_init, pulsesAll, *pulsesInRecord, 1/record->delta_t, &pulsetemplate, &pulseheighttemplate, &covariance, &weight))
		{
			message = "Cannot run routine calculateTemplate in CALIBRATION mode";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		if (writeLibrary(*reconstruct_init, 1/record->delta_t, pulseheighttemplate, pulsetemplate, covariance, weight, appendToLibrary, &inLibObject))
		{
			message = "Cannot run routine writeLibrary in CALIBRATION mode";
			EP_EXIT_ERROR(message,EPFAIL);
		}

		gsl_vector_free(pulsetemplate);
		gsl_matrix_free(weight);
		gsl_matrix_free(covariance);
	}

	// Close intermediate output FITS file if it is necessary
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (fits_close_file(dtcObject,&status))
		{
			message = "Cannot close file " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	gsl_vector_free(invector);

	return;
}
/*xxxx end of SECTION A xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A1 ************************************************************
* createLibrary: This function creates the pulse templates library FITS file, if it does not exist yet, or open it.
*
* - If it exists => Open it and 'appendToLibrary' = true
* - If it does not exist => Create it and 'appendToLibrary' = false
* 	- Write 'EVENTCNT'=1 in the LIBRARY extension
* 	- Write info keywords about the running parameters in the Primary extension
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses parameters to call the library ('library_file') and to write input parameters
*                     info in the Primary extension of the library file
* - appendToLibrary: It is going to be used by the function 'writeLibrary'
* - inLibObject: Object which contains information of the library FITS file (it is going to be used also by 'writeLibrary') 
* - create: To write the 'CREATOR' keyword
****************************************************************************/
int createLibrary(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, fitsfile **inLibObject, const char * create)
{
	int status = EPOK;
	string message = "";

	char inLibName[256];
	strncpy(inLibName, reconstruct_init->library_file,255);
	inLibName[255]='\0';
	
	char extname[10];
	strncpy(extname,"LIBRARY",9);
	extname[9]='\0';

	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];
	char keyname[10];

	long eventcntLib;

	// If the library exists => Open it
	if (fileExists(string(inLibName)))
	{
		*appendToLibrary = true;

		if (fits_open_file(inLibObject,inLibName,READWRITE,&status))
		{
		    message = "Cannot open library file " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
		    message = "Cannot move to HDU " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}
	else	// If the library does not exist yet => Create it
	{
		*appendToLibrary = false;
		if (fits_create_file(inLibObject, inLibName, &status))
		{
		    message = "Cannot create library file " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_open_file(inLibObject,inLibName,READWRITE,&status))
		{
		    message = "Cannot open library file " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Create LIBRARY extension
		if (fits_create_tbl(*inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
		{
		    message = "Cannot create table " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		eventcntLib = 1;
		strncpy(keyname,"EVENTCNT",9);
		keyname[9]='\0';
		assert(eventcntLib > 0);
		if (fits_write_key(*inLibObject,TLONG,keyname,&eventcntLib,NULL,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Primary HDU
		strncpy(extname,"Primary",9);
		extname[9]='\0';
		if (fits_movabs_hdu(*inLibObject, 1, NULL, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strncpy(keyname,"CREATOR",9);
		keyname[9]='\0';
		string creator (string("File CREATED by") + ' ' + (string) create);
		strncpy(keyvalstr,creator.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		//strcpy(keyname,"PROC0");
		strncpy(keyname,"PROC0",9);
		keyname[9]='\0';
		const char * charproc= "PROC0 Starting parameter list";
		strncpy(keyvalstr,charproc,999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		string strproc (string("RecordFile = ") + reconstruct_init->record_file);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		strproc=string("TesEventFile = ") + reconstruct_init->event_file;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		strproc=string("LibraryFile = ") + reconstruct_init->library_file;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		strproc=string("NoiseFile = ") + reconstruct_init->noise_file;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_mode[125];		snprintf(str_mode,125,"%d",reconstruct_init->mode);
		strproc = string("mode = ") + string(str_mode);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		strproc=string("PixelType = ") + reconstruct_init->PixelType;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		strproc=string("FilterDomain = ") + reconstruct_init->FilterDomain;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		strproc=string("FilterMethod = ") + reconstruct_init->FilterMethod;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strproc=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_LagsOrNot[125];	snprintf(str_LagsOrNot,125,"%d",reconstruct_init->LagsOrNot);
		strproc=string("LagsOrNot = ") + string(str_LagsOrNot);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_OFIter[125];		snprintf(str_OFIter,125,"%d",reconstruct_init->OFIter);
		strproc=string("OFIter = ") + string(str_OFIter);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_OFLib[125];      	snprintf(str_OFLib,125,"%d",reconstruct_init->OFLib);
		strproc=string("OFLib = ") + string(str_OFLib);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strproc=string("OFInterp = ") + reconstruct_init->OFInterp;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strproc=string("OFStrategy = ") + reconstruct_init->OFStrategy;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_OFLength[125];		snprintf(str_OFLength,125,"%d",reconstruct_init->OFLength);
		strproc=string("OFLength = ") + string(str_OFLength);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_maxPulsesPerRecord[125];	snprintf(str_maxPulsesPerRecord,125,"%d",reconstruct_init->maxPulsesPerRecord);
		strproc=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_pulse_length[125];	snprintf(str_pulse_length,125,"%d",reconstruct_init->pulse_length);
		strproc=string("PulseLength = ") + string(str_pulse_length);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tauFall[125];		snprintf(str_tauFall,125,"%e",reconstruct_init->tauFall);
		strproc=string("tauFall = ") + string(str_tauFall);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_scaleFactor[125];	snprintf(str_scaleFactor,125,"%f",reconstruct_init->scaleFactor);
		strproc=string("scaleFactor = ") + string(str_scaleFactor);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_samplesUp[125];	snprintf(str_samplesUp,125,"%f",reconstruct_init->samplesUp);
		strproc=string("samplesUp = ") + string(str_samplesUp);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_nSgms[125];	    	snprintf(str_nSgms,125,"%f",reconstruct_init->nSgms);
		strproc=string("nSgms = ") + string(str_nSgms);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LrsT[125];		snprintf(str_LrsT,125,"%e",reconstruct_init->LrsT);
		strproc=string("LrsT = ") + string(str_LrsT);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LbT[125];		snprintf(str_LbT,125,"%e",reconstruct_init->LbT);
		strproc=string("LbT = ") + string(str_LbT);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_monoenergy[125];	snprintf(str_monoenergy,125,"%f",reconstruct_init->monoenergy);
		strproc=string("monoenergy = ") + string(str_monoenergy);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_intermediate[125];     snprintf(str_intermediate,125,"%d",reconstruct_init->intermediate);
		strproc=string("intermediate = ") + string(str_intermediate);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strproc=string("detectFile = ") + reconstruct_init->detectFile;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strproc=string("filterFile = ") + reconstruct_init->filterFile;
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_clobber[125];      	snprintf(str_clobber,125,"%d",reconstruct_init->clobber);
		strproc=string("clobber = ") + string(str_clobber);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_SaturationValue[125];	snprintf(str_SaturationValue,125,"%e",reconstruct_init->SaturationValue);
		strproc=string("SaturationValue = ") + string(str_SaturationValue);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse1[125];	snprintf(str_tstartPulse1,125,"%d",reconstruct_init->tstartPulse1);
		strproc=string("tstartPulse1 = ") + string(str_tstartPulse1);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse2[125];	snprintf(str_tstartPulse2,125,"%d",reconstruct_init->tstartPulse2);
		strproc=string("tstartPulse2 = ") + string(str_tstartPulse2);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse3[125];	snprintf(str_tstartPulse3,125,"%d",reconstruct_init->tstartPulse3);
		strproc=string("tstartPulse3 = ") + string(str_tstartPulse3);
		strncpy(keyvalstr,strproc.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		charproc= "PROC0 Ending parameter list";
		strncpy(keyvalstr,charproc,999);
		keyvalstr[999]='\0';
		fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		if (status != 0)
		{
			message = "Cannot write keyword in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return (EPOK);
}
/*xxxx end of SECTION A1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A2 ************************************************************
* createDetectFile function: This function creates an intermediate FITS file with some useful info (during the developing) if the 'intermediate' input parameter is 1.
*                            The intermediate FITS file will contain 2 extensions, PULSES and TESTINFO.
*                            The PULSES extension contains some found pulses info: TSTART, I0 (the pulse itself), TEND, TAURISE, TAUFALL and QUALITY.
*                            The TESTINFO extension contains FILDER (the low-pass filtered and derived records) and THRESHOLD.
*
* - If it exists => Check 'clobber' in order to overwrite or not
* - If it does not exist => Create it
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses parameters to call the intermediate file ('detectFile'), to write input parameters
*                     info in the Primary extension of the intermediate file and also 'clobber'
* - samprate: Sampling rate
* - create: To write the 'CREATOR' keyword
* - dtcObject: Object which contains information of the intermediate FITS file (it is going to be used also by 'writeTestInfo' and 'writePulses') 
***************************************************************************/
int createDetectFile(ReconstructInitSIRENA* reconstruct_init, double samprate, const char * create, fitsfile **dtcObject)
{
	int status = EPOK;
	string message = "";

	char dtcName[256];
	strncpy(dtcName,reconstruct_init->detectFile,255);
	dtcName[255]='\0'; // enforce zero ending string in case of buffer overflows

	// Create intermediate output FITS file: If it does not exist yet
	// If dtcName does not finish as '.fits' and the file dtcName+'.fits' already exists =>
	// => Data are appended to dtcName file => Must not be allowed
	// '.fits' => 5 characters
	if (strlen(dtcName)<6)
	{
		// dtcName has 5 or less characters => Does not contain '.fits' =>Append '.fits' to dtcName
		strcat(dtcName,".fits");
	}
	else
	{
		// Check if dtcName has '.fits' and if not, append it
		if (strncmp(strndup(dtcName+strlen(dtcName)-5, 5),".fits",5) != 0)
		{
			// dtcName does not finish as '.fits' => Append '.fits' to dtcName
			char dtcNameaux[256];
			snprintf(dtcNameaux,256,dtcName);
			strcat(dtcNameaux,".fits");
			strncpy(dtcName,dtcNameaux,255);
			dtcName[255]='\0';
		}
	}

	// Create intermediate file (if file already exists => check clobber)
	if (fileExists(string(dtcName)) && (reconstruct_init->clobber == 1))
	{
		if (remove(dtcName))
		{
			message = "Output intermediate file already exists & cannot be deleted ("+string(strerror(errno))+")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	    }
	}
	else if (fileExists(string(dtcName)) && (reconstruct_init->clobber == 0))
	{
		message = "Output intermediate file already exists: must not be overwritten";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	if (!fileExists(string(dtcName)))
	{
		if(fits_create_file(dtcObject, dtcName, &status))
		{
			message = "Cannot create output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		message = "Create Detect Fits File: " + string(dtcName);
	}

	// Create extension PULSES
	// To work with tables (extensions)
	char *tt[1];
	char *tf[1];
	char *tu[1];
	
	// To write keywords
	char keyname[10];
	char keyvalstr[1000];
	char extname[10];

	if (fits_open_file(dtcObject,dtcName,READWRITE,&status))
	{
		message = "Cannot open output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (reconstruct_init->clobber == 1)
	{
		strncpy(extname,"PULSES",9);
		extname[9]='\0';
		
		// PULSES HDU
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_movnam_hdu(*dtcObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strncpy(keyname,"MODE",9);
		keyname[9]='\0';
		fits_write_key(*dtcObject,TINT,keyname,&(reconstruct_init->mode),NULL,&status);
		
		strncpy(keyname,"EVENTSZ",9);
		keyname[9]='\0';
		assert(reconstruct_init->pulse_length > 0);
		fits_write_key(*dtcObject,TINT,keyname,&(reconstruct_init->pulse_length),NULL,&status);
		if (reconstruct_init->mode == 0)
		{
			strcpy(keyname,"ENERGY");
			fits_write_key(*dtcObject,TDOUBLE,keyname,&(reconstruct_init-> monoenergy),NULL,&status);
		}
		strncpy(keyname,"SAMPRATE",9);
		keyname[9]='\0';
		assert(samprate > 0);
		fits_write_key(*dtcObject,TDOUBLE,keyname,&samprate,NULL,&status);
		if (status != 0)
		{
			message = "Cannot write keyword some keyword in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// TESTINFO HDU
		strncpy(extname,"TESTINFO",9);
		extname[9]='\0';
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Primary HDU
		strncpy(extname,"Primary",9);
		extname[9]='\0';
		if (fits_movabs_hdu(*dtcObject, 1, NULL, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strncpy(keyname,"CREATOR",9);
		keyname[9]='\0';
		string creator (string("File CREATED by") + ' ' + (string) create);
		strncpy(keyvalstr,creator.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strncpy(keyname,"HISTORY",9);
		keyname[9]='\0';
		const char * charhistory= "HISTORY Starting parameter list";
		strncpy(keyvalstr,charhistory,999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		string strhistory (string("RecordFile = ") + reconstruct_init->record_file);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key_longstr(*dtcObject,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("TesEventFile = ") + reconstruct_init->event_file;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key_longstr(*dtcObject,keyname,keyvalstr,NULL,&status);

		strhistory=string("LibraryFile = ") + reconstruct_init->library_file;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key_longstr(*dtcObject,keyname,keyvalstr,NULL,&status);

		strhistory=string("NoiseFile = ") + reconstruct_init->noise_file;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key_longstr(*dtcObject,keyname,keyvalstr,NULL,&status);

		char str_mode[125];		snprintf(str_mode,125,"%d",reconstruct_init->mode);
		strhistory = string("mode = ") + string(str_mode);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("PixelType = ") + reconstruct_init->PixelType;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("FilterDomain = ") + reconstruct_init->FilterDomain;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("FilterMethod = ") + reconstruct_init->FilterMethod;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LagsOrNot[125];      	snprintf(str_LagsOrNot,125,"%d",reconstruct_init->LagsOrNot);
		strhistory=string("LagsOrNot = ") + string(str_LagsOrNot);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_OFIter[125];     	snprintf(str_OFIter,125,"%d",reconstruct_init->OFIter);
		strhistory=string("OFIter = ") + string(str_OFIter);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_OFLib[125];      	snprintf(str_OFLib,125,"%d",reconstruct_init->OFLib);
		strhistory=string("OFLib = ") + string(str_OFLib);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
				
		strhistory=string("OFInterp = ") + reconstruct_init->OFInterp;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
				
		strhistory=string("OFStrategy = ") + reconstruct_init->OFStrategy;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
				
		char str_OFLength[125];		snprintf(str_OFLength,125,"%d",reconstruct_init->OFLength);
		strhistory=string("OFLength = ") + string(str_OFLength);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);

		char str_maxPulsesPerRecord[125];	snprintf(str_maxPulsesPerRecord,125,"%d",reconstruct_init->maxPulsesPerRecord);
		strhistory=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_pulse_length[125];	snprintf(str_pulse_length,125,"%d",reconstruct_init->pulse_length);
		strhistory=string("PulseLength = ") + string(str_pulse_length);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tauFall[125];		snprintf(str_tauFall,125,"%e",reconstruct_init->tauFall);
		strhistory=string("tauFall = ") + string(str_tauFall);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_scaleFactor[125];	snprintf(str_scaleFactor,125,"%f",reconstruct_init->scaleFactor);
		strhistory=string("scaleFactor = ") + string(str_scaleFactor);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_samplesUp[125];	snprintf(str_samplesUp,125,"%f",reconstruct_init->samplesUp);
		strhistory=string("samplesUp = ") + string(str_samplesUp);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_nSgms[125];	    	snprintf(str_nSgms,125,"%f",reconstruct_init->nSgms);
		strhistory=string("nSgms = ") + string(str_nSgms);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LrsT[125];		snprintf(str_LrsT,125,"%e",reconstruct_init->LrsT);
		strhistory=string("LrsT = ") + string(str_LrsT);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_LbT[125];		snprintf(str_LbT,125,"%e",reconstruct_init->LbT);
		strhistory=string("LbT = ") + string(str_LbT);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_monoenergy[125];	snprintf(str_monoenergy,125,"%f",reconstruct_init->monoenergy);
		strhistory=string("monoenergy = ") + string(str_monoenergy);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_intermediate[125];      snprintf(str_intermediate,125,"%d",reconstruct_init->intermediate);
		strhistory=string("intermediate = ") + string(str_intermediate);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("detectFile = ") + reconstruct_init->detectFile;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key_longstr(*dtcObject,keyname,keyvalstr,NULL,&status);
		
		strhistory=string("filterFile = ") + reconstruct_init->filterFile;
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key_longstr(*dtcObject,keyname,keyvalstr,NULL,&status);
		
		char str_clobber[125];      	snprintf(str_clobber,125,"%d",reconstruct_init->clobber);
		strhistory=string("clobber = ") + string(str_clobber);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_SaturationValue[125];	snprintf(str_SaturationValue,125,"%e",reconstruct_init->SaturationValue);
		strhistory=string("SaturationValue = ") + string(str_SaturationValue);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse1[125];	snprintf(str_tstartPulse1,125,"%d",reconstruct_init->tstartPulse1);
		strhistory=string("tstartPulse1 = ") + string(str_tstartPulse1);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse2[125];	snprintf(str_tstartPulse2,125,"%d",reconstruct_init->tstartPulse2);
		strhistory=string("tstartPulse2 = ") + string(str_tstartPulse2);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_tstartPulse3[125];	snprintf(str_tstartPulse3,125,"%d",reconstruct_init->tstartPulse3);
		strhistory=string("tstartPulse3 = ") + string(str_tstartPulse3);
		strncpy(keyvalstr,strhistory.c_str(),999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		charhistory= "HISTORY Ending parameter list";
		strncpy(keyvalstr,charhistory,999);
		keyvalstr[999]='\0';
		fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		if (status != 0)
		{
			message = "Cannot write some keyword in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return EPOK;
}
/*xxxx end of SECTION A2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A3 ************************************************************
* filderLibrary: This function calculates the low-pass filtered and derivative of the models ('pulse_templates') of the library (only necessary
*                if first record), and it stores the 'pulse_templates_filder' and the 'maxDERs' in the 'reconstruct_init' structure.
*                The maximum of the low-pass filtered and derived pulse has to be compared to the 'maxDERs' to select the appropriate model.
*
* - Check if it is the first record
* - Low-pass filtered and derive the models ('pulse_templates') of the library
* - Store the low-pass filtered derivatives in 'pulse_templates_filder'
* - Calculate the maximum of the low-pass filtered and derived models ('maxDERs')
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses parameters to filter ('scaleFactor' and 'tauFall') and others to handle the pulse templates
* 		      and their derivatives ('ntemplates', 'pulse_templates', 'pulse_templates_filder' and 'maxDERs')
* - samprate: Sampling rate
******************************************************************************/
int filderLibrary(ReconstructInitSIRENA** reconstruct_init, double samprate)
{
	int status = EPOK;
	string message = "";

	if ((*reconstruct_init)->library_collection->pulse_templates_filder[0].template_duration == -1)		// First record
	{
		double scaleFactor = (*reconstruct_init)->scaleFactor;
		double tauFALL = (*reconstruct_init)->tauFall;

		// Check boxLength
		double cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
		int boxLength = (int) ((1/cutFreq) * samprate);
		if (boxLength <= 1)
		{
			message = "lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
			EP_PRINT_ERROR(message,-999);
		}

		// 'pulse_templates' are filtered and derived
		gsl_vector *model = gsl_vector_alloc((*reconstruct_init)->library_collection->pulse_templates[0].template_duration);
		for (int i=0; i<(*reconstruct_init)->library_collection->ntemplates; i++)
		{
			gsl_vector_memcpy(model, (*reconstruct_init)->library_collection->pulse_templates[i].ptemplate);

			// PULSE TEMPLATE: Low-pass filtering
			status = lpf_boxcar(&model,model->size,tauFALL*scaleFactor,samprate);
			if (status == 1)
			{
				message = "Cannot run routine lpf_boxcar for low-pass filtering";
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}
			if (status == 3)
			{
				status = EPOK;
			}
			if (status == 4)
			{
				message = "lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low";
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}

			// PULSE TEMPLATE: Derivative after filtering
			if (derivative (&model,model->size))
			{
				message = "Cannot run routine derivative in filderLibrary";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			// Store the low-pass filtered derivatives in 'pulse_templates_filder'
			gsl_vector_memcpy((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,model);

			(*reconstruct_init)->library_collection->pulse_templates_filder[i].template_duration = (*reconstruct_init)->library_collection->pulse_templates[i].template_duration;

			// Calculate the maximum of the low-pass filtered and derived models
			gsl_vector_set((*reconstruct_init)->library_collection->maxDERs,i,gsl_vector_max(model));
		}

		gsl_vector_free(model);
	}

	return(EPOK);
}
/*xxxx end of SECTION A3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A4 ************************************************************
* loadRecord: This function loads the struture 'record' in the 'adc_double' vector.
*
* It checks if the record has been filled in with 0's => It only loads the first values (which are different from 0).
* 
* Parameters:
* - record: Structure that contains the input record
* - time_record: Starting time of the record (output)
* - adc_double: The record to process is stored in 'adc_double'(input/output)
******************************************************************************/
int loadRecord(TesRecord* record, double *time_record, gsl_vector **adc_double)
{
	*time_record = record->time;
	for (int i=0;i<record->trigger_size;i++)
	{
		gsl_vector_set(*adc_double,i,record->adc_double[i]);
	}

	// Just in case the last record has been filled in with 0's => Re-allocate 'invector'
	//if (gsl_vector_ispos(*adc_double) != 1)
	if ((gsl_vector_get(*adc_double,record->trigger_size-1) == 0.0) && (gsl_vector_get(*adc_double,record->trigger_size-2) == 0.0))
	// Changed because your baseline could be near to zero, even having negative values of the record
	{
		// Know the new dimension of the last record (elements different from 0)
		long eventszLastRecord;
		eventszLastRecord = gsl_vector_min_index(*adc_double);

		// Auxiliary vector
		gsl_vector *vector_aux = gsl_vector_alloc(eventszLastRecord);
		gsl_vector_view temp;
		temp = gsl_vector_subvector(*adc_double,0,eventszLastRecord);
		gsl_vector_memcpy(vector_aux,&temp.vector);

		// Free invector, allocate with the new dimension and free the auxiliary vector
		gsl_vector_free(*adc_double);
		*adc_double = gsl_vector_alloc(eventszLastRecord);
		gsl_vector_memcpy(*adc_double,vector_aux);
		gsl_vector_free(vector_aux);
	}

	return (EPOK);
}
/*xxxx end of SECTION A4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A5 ************************************************************
* procRecord function:  This function processes the input record.
*
* - Declare and initialize variables
* - Allocate GSL vectors
* - Low-pass filtering and derivative
* - Find the events (pulses) in the record
* - Calculate the tend of the found pulses and check if the pulse is saturated
* - Obtain the approximate rise and fall times of each pulse (to be done)
* - Load the found pulses data in the input/output 'foundPulses' structure
* - Write test info (if 'reconstruct_init->intermediate'=1)
* - Write pulses info in intermediate output FITS file (if 'reconstruct_init->intermediate'=1)
* - Free allocate of GSL vectors
*
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses parameters to filter ('scaleFactor', 'tauFall'),
*                     to find pulses ('pulse_length', 'samplesUp', 'nSgms', 'LrsT', 'LbT', 'maxPulsesPerRecord')
*                     and to wrtite info ('detectFile') if 'reconstruct_init->intermediate'=1
* - tstartRecord: Starting time of the record (in order to calculate absolute times)
* - samprate: Sampling rate (in order to low-pass filter)
* - dtcObject: Object which contains information of the intermediate FITS file (to be written if 'intermediate'=1)
* - record: Input record
* - foundPulses: Input/output structure where the found pulses info is stored
****************************************************************************/
int procRecord(ReconstructInitSIRENA** reconstruct_init, double tstartRecord, double samprate, fitsfile *dtcObject, gsl_vector *record, gsl_vector *recordWithoutCovert2R,PulsesCollection *foundPulses)
{
	int status = EPOK;
	string message = "";

	// Declare and initialize variables
	int numPulses = 0;
	double threshold = 0.0;

	double stopCriteriaMKC = 1.0;	// Used in medianKappaClipping
	                               	// Given in %
	double kappaMKC = 3.0;		// Used in medianKappaClipping
	double levelPrvPulse = 100.0;  	// Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse

	gsl_vector_view temp;

	double scaleFactor = (*reconstruct_init)->scaleFactor;
	double tauFALL = (*reconstruct_init)->tauFall;
	int sizePulse_b = (*reconstruct_init)->pulse_length;
	double samplesUp = (*reconstruct_init)->samplesUp;
	double nSgms = (*reconstruct_init)->nSgms;
	double Lrs = (int) ((*reconstruct_init)->LrsT*samprate);	// Running sum length (in the RS filter case): 'LrsT' in samples
	double Lb = (int) ((*reconstruct_init)->LbT*samprate); 		// Baseline averaging length (in the RS filter case): 'LbT' in samples

	// Allocate GSL vectors
	gsl_vector *recordNOTFILTERED = gsl_vector_alloc(record->size); // Record without having been filtered
	gsl_vector *recordDERIVATIVE = gsl_vector_alloc(record->size);  // Derivative of 'invectorFILTERED'

	// To look for pulses
	gsl_vector *tstartgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *tendgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *qualitygsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *pulseHeightsgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *maxDERgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector_set_zero(qualitygsl);
	gsl_vector_set_zero(pulseHeightsgsl);		// In order to choose the proper pulse model to calculate
	                                                // the adjusted derivative and to fill in the ESTENRGY column
	                                                // in the output FITS file
	gsl_vector_set_zero(maxDERgsl);

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	// Check boxLength
	double cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
	int boxLength = (int) ((1/cutFreq) * samprate);
	if (boxLength <= 1)
	{
		message = "lpf_boxcar: tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
		//?? Too many maessages
	}

	// Low-pass filtering
	gsl_vector_memcpy(recordNOTFILTERED,record);
	status = lpf_boxcar(&record,record->size,scaleFactor*tauFALL,samprate);
	if (status == EPFAIL)
	{
		message = "Cannot run routine lpf_boxcar for low pass filtering";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (status == 3)
	{
		status = EPOK;
	}
	if (status == 4)
	{
		message = "lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low";
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	// Derivative after filtering
	if (derivative (&record, record->size))
	{
		message = "Cannot run routine derivative for derivative after filtering";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_vector_memcpy(recordDERIVATIVE,record);

	gsl_vector *recordDERIVATIVEOriginal = gsl_vector_alloc(recordDERIVATIVE->size);
	gsl_vector_memcpy(recordDERIVATIVEOriginal,recordDERIVATIVE);
	/*cout<<"Not Filtered * Derived"<<endl;
	for (int i=20998;i<20998+50;i++)
	{
	      cout<<gsl_vector_get(recordNOTFILTERED,i)<<" "<<gsl_vector_get(recordDERIVATIVE,i)<<endl;
	}*/

	// Find the events (pulses) in the record
	if ((*reconstruct_init)->mode == 1)	// In PRODUCTION mode
	{
		int tstartFirstEvent = 0;
		bool triggerCondition;
		int flagTruncated;
		int tstartProvided;

		if (InitialTriggering (recordDERIVATIVE, samplesUp, nSgms,
			tauFALL, scaleFactor, samprate, stopCriteriaMKC, kappaMKC,
			&triggerCondition, &tstartFirstEvent, &flagTruncated,
			&threshold, (*reconstruct_init)->tstartPulse1))
		{
			message = "Cannot run routine InitialTriggering";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		if (FindSecondaries ((*reconstruct_init)->maxPulsesPerRecord,
			recordDERIVATIVE, threshold,
			samplesUp,(*reconstruct_init),
			tstartFirstEvent,
			&numPulses,&tstartgsl,&qualitygsl, &maxDERgsl))
		{
			message = "Cannot run routine FindSecondaries";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	else if ((*reconstruct_init)->mode == 0)	// In CALIBRATION mode
	{
		if (findPulsesCAL (recordNOTFILTERED, recordDERIVATIVE, &tstartgsl, &qualitygsl, &pulseHeightsgsl, &maxDERgsl,
			&numPulses, &threshold,
			tauFALL, scaleFactor, samprate,
			samplesUp, nSgms,
			Lb, Lrs,
			(*reconstruct_init),
			stopCriteriaMKC,
			kappaMKC))
		{
			message = "Cannot run routine findPulsesCAL";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	(*reconstruct_init)->threshold = threshold;

	// Write test info
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (writeTestInfo((*reconstruct_init), recordDERIVATIVEOriginal, threshold, dtcObject))
		{
			message = "Cannot run routine writeTestInfo";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	gsl_vector_free(recordDERIVATIVEOriginal);

	// Calculate the tend of the found pulses and check if the pulse is saturated
	for (int i=0;i<numPulses;i++)
	{
		gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i)+sizePulse_b);	//tend_i = tstart_i + (ntaus*tauFALL*samprate)

		if (gsl_vector_get(tendgsl,i) >= recordDERIVATIVE->size)		// Truncated pulses at the end of the record
		{
			gsl_vector_set(tendgsl,i,(recordDERIVATIVE->size)-1);
			gsl_vector_set (qualitygsl,i,2);
		}

		if ((numPulses != 1) && (i != numPulses-1)) 				// More than one pulse in the record and not the last one
		{
			if (gsl_vector_get(tendgsl,i) > gsl_vector_get(tstartgsl,i+1))
			{
				gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i+1));
			}
		}
		
		temp = gsl_vector_subvector(recordWithoutCovert2R,gsl_vector_get(tstartgsl,i),gsl_vector_get(tendgsl,i)-gsl_vector_get(tstartgsl,i));
		if (gsl_vector_max(&temp.vector) > (*reconstruct_init)->SaturationValue)	gsl_vector_set(qualitygsl,i,gsl_vector_get(qualitygsl,i)+10);
	}
	
	// Obtain the approximate rise and fall times of each pulse
	gsl_vector *tauRisegsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *tauFallgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector_set_zero(tauRisegsl);
	gsl_vector_set_zero(tauFallgsl);
	/*if (obtainTau (recordNOTFILTERED, tstartgsl, tendgsl, *numPulses, &tauRisegsl, &tauFallgsl))
	{
	    message = "Cannot run routine obtainTau to calculate pulses slopes";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}*/

	// Load the found pulses data in the input/output 'foundPulses' structure
	foundPulses->ndetpulses = numPulses;
	foundPulses->pulses_detected = new PulseDetected[numPulses];
	for (int i=0;i<numPulses;i++)
	{
		foundPulses->pulses_detected[i].pulse_duration = gsl_vector_get(tendgsl,i)-gsl_vector_get(tstartgsl,i);

		// 'grade1' will be known after running 'runEnergy' (but initialize for library creation!)
		foundPulses->pulses_detected[i].grade1 = 0;
		if (i!= 0)
		{
			foundPulses->pulses_detected[i].grade2 = (int)(gsl_vector_get(tstartgsl,i)-gsl_vector_get(tendgsl,i-1));
			foundPulses->pulses_detected[i].grade2_1 = (int)(gsl_vector_get(tstartgsl,i)-gsl_vector_get(tstartgsl,i-1));
		}
		else
		{
			foundPulses->pulses_detected[i].grade2 = (*reconstruct_init)->pulse_length;
			foundPulses->pulses_detected[i].grade2_1 = (*reconstruct_init)->pulse_length;
		}

		foundPulses->pulses_detected[i].pulse_adc = gsl_vector_alloc(foundPulses->pulses_detected[i].pulse_duration);
		temp = gsl_vector_subvector(recordNOTFILTERED,gsl_vector_get(tstartgsl,i),foundPulses->pulses_detected[i].pulse_duration);
		//foundPulses->pulses_detected[i].pulse_adc = gsl_vector_alloc(foundPulses->pulses_detected[i].pulse_duration+safetyMargin);
		//temp = gsl_vector_subvector(recordNOTFILTERED,gsl_vector_get(tstartgsl,i)-safetyMargin,foundPulses->pulses_detected[i].pulse_duration+safetyMargin);

		gsl_vector_memcpy(foundPulses->pulses_detected[i].pulse_adc,&temp.vector);

		foundPulses->pulses_detected[i].Tstart = gsl_vector_get(tstartgsl,i)/samprate+tstartRecord;
		foundPulses->pulses_detected[i].Tend = gsl_vector_get(tendgsl,i)/samprate+tstartRecord;
		foundPulses->pulses_detected[i].riseTime = gsl_vector_get(tauRisegsl,i);
		foundPulses->pulses_detected[i].fallTime = gsl_vector_get(tauFallgsl,i);
		foundPulses->pulses_detected[i].pulse_height = gsl_vector_get(pulseHeightsgsl,i);
		foundPulses->pulses_detected[i].maxDER = gsl_vector_get(maxDERgsl,i);
		// 'energy' will be known after running 'runEnergy'
		foundPulses->pulses_detected[i].quality = gsl_vector_get(qualitygsl,i);
		//cout<<"Pulse "<<i<<" tstart="<<gsl_vector_get(tstartgsl,i)<<", maxDER= "<<foundPulses->pulses_detected[i].maxDER<<", pulse_duration= "<<foundPulses->pulses_detected[i].pulse_duration<<",quality= "<<foundPulses->pulses_detected[i].quality<<endl;
		//cout<<gsl_vector_get(tstartgsl,i)<<endl;
	}

	// Write pulses info in intermediate output FITS file
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (writePulses (reconstruct_init, samprate, tstartRecord, recordNOTFILTERED, numPulses, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, dtcObject))
		{
			message = "Cannot run routine writePulses to write pulses in intermediate output file";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	// Free allocate of GSL vectors
	gsl_vector_free(recordNOTFILTERED);
	gsl_vector_free(recordDERIVATIVE);

	gsl_vector_free(tstartgsl);
	gsl_vector_free(tendgsl);
	gsl_vector_free(qualitygsl);
	gsl_vector_free(pulseHeightsgsl);
	gsl_vector_free(maxDERgsl);

	gsl_vector_free(tauRisegsl);
	gsl_vector_free(tauFallgsl);

	return EPOK;
}
/*xxxx end of SECTION A5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A6 ************************************************************
* writePulses function: This function writes the data of the found pulses in the record in the intermediate FITS file.
*                       It writes the pulses info in the PULSES extension.
*                       The pulses info is: TSTART, I0 (the pulse itself), TEND, TAURISE, TAUFALL and QUALITY.
*
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'detectFile', 'pulse_length' and 'clobber'
* - samprate: Sampling rate (to convert samples to seconds)
* - initialtime: Starting time of the record (in order to calculate absolute times)
* - invectorNOTFIL: Original record (neither low-pass filtered nor derived)
* - numPulsesRecord: Number of pulses found in the record
* - tstart: Tstarts of the found pulses
* - tend: Tends of the found pulses
* - quality: Quality of the found pulses
* - taurise: Taurise of the forund pulses (to be done)
* - taufall: Taufall of the found pulses (to be done)
* - dtcObject: Object which contains information of the intermediate FITS file 
******************************************************************************/
int writePulses(ReconstructInitSIRENA** reconstruct_init, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, fitsfile *dtcObject)
{
	int status = EPOK;
	string message = "";

	// Declare variables
	int t0;		// First value of index of pulse
	gsl_matrix *vgslout2;

	gsl_vector_view temp;

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	// If intermediate=1 => First record, createDetectFile
	//	                Change clobber to 2
	//                   => Not first record, append info to output intermediate file (because clobber is 2)
	long totalpulses;	// It is necessary to know the row where the info is going to be written
	if ((*reconstruct_init)->clobber == 1)
	{
		totalpulses = 0;
		if ((*reconstruct_init)->mode == 0)	(*reconstruct_init)->clobber = 2;
	}
	else if ((*reconstruct_init)->clobber == 2)
	{
		if (fits_get_num_rows(dtcObject,&totalpulses, &status))
		{
			message = "Cannot get number of rows in " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
	}

	if (numPulsesRecord!=0)
	{
		vgslout2 = gsl_matrix_alloc(numPulsesRecord,(*reconstruct_init)->pulse_length);

		// Converting samples to seconds
		for (int i=0; i<numPulsesRecord; i++)
		{
			t0 = gsl_vector_get (tstart,i);
			gsl_vector_set(tstart,i,initialtime + (gsl_vector_get (tstart,	i) * (1/samprate)));
			gsl_vector_set(tend,i,initialtime + (gsl_vector_get (tend,	i) * (1/samprate)));

			if (invectorNOTFIL->size - t0 > (*reconstruct_init)->pulse_length)	//The invectorNOTFIL has more sampless than sizePulse
			{
				temp = gsl_vector_subvector(invectorNOTFIL,t0,(*reconstruct_init)->pulse_length);
				gsl_matrix_set_row(vgslout2, i, &temp.vector);
			}
			else 									// The invectorNOTFIL has less bins than sizePulse (truncated)
			{
				for (int j=0; j<(invectorNOTFIL->size) - t0; j++)
				{
					if (t0 == -1) t0 = 0;
					gsl_matrix_set (vgslout2,i,j,gsl_vector_get(invectorNOTFIL,j+t0));
				}

				for (int j=(invectorNOTFIL->size)-t0; j< (*reconstruct_init)->pulse_length; j++) {gsl_matrix_set (vgslout2,i,j,0.0);}
			}
		}

		IOData obj;

		// PULSES HDU
		// Creating TSTART Column
		obj.inObject = dtcObject;
		obj.nameTable = new char [255];
		strcpy(obj.nameTable,"PULSES");
		obj.iniRow = totalpulses+1;
		obj.endRow = totalpulses+numPulsesRecord-1+1;
		obj.iniCol = 0;
		obj.nameCol = new char [255];
		strcpy(obj.nameCol,"TSTART");
		obj.type = TDOUBLE;
		obj.unit = new char [255];
		strcpy(obj.unit,"seconds");

		temp = gsl_vector_subvector(tstart,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TSTART";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating I0 Column
		strcpy(obj.nameCol,"I0");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"ADC");
		if (writeFitsComplex(obj, vgslout2))
		{
			message = "Cannot run routine writeFitsComplex for column IO";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating TEND Column
		strcpy(obj.nameCol,"TEND");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"seconds");
		temp = gsl_vector_subvector(tend,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TEND";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating TAURISE Column
		strcpy(obj.nameCol,"TAURISE");
		temp = gsl_vector_subvector(taurise,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TAURISE";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating TAUFALL Column
		strcpy(obj.nameCol,"TAUFALL");
		temp = gsl_vector_subvector(taufall,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column TAUFALL";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		// Creating QUALITY Column
		strcpy(obj.nameCol,"QUALITY");
		obj.type = TSHORT;
		strcpy(obj.unit," ");
		temp = gsl_vector_subvector(quality,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column QUALITY";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		totalpulses = totalpulses + numPulsesRecord;

		// Free allocate GSL vectors
		gsl_matrix_free (vgslout2);

		delete [] obj.nameTable;
		delete [] obj.nameCol;
		delete [] obj.unit;
	}

	return (EPOK);
}
/*xxxx end of SECTION A6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A7 ************************************************************
* writeTestInfo function: This function writes the TESTINFO extension in the intermediate FITS file.
*                         The written columns are FILDER (low-pass filtered and derived record) and THRESHOLD.
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'detectFile'
* - recordDERIVATIVE: Input record low-pass filtered and derived
* - threshold: Threshold used to find pulses
* - dtcObject: Object which contains information of the intermediate FITS file
******************************************************************************/
int writeTestInfo(ReconstructInitSIRENA* reconstruct_init, gsl_vector *recordDERIVATIVE, double threshold, fitsfile *dtcObject)
{
	int status = EPOK;
	string message = "";

	long totalrecords;

	char dtcName[256];
	strncpy(dtcName,reconstruct_init->detectFile,255);
	dtcName[255]='\0'; // enforce zero ending string in case of buffer overflows

	// To work with tables (extensions)
	char extname[20];

	strcpy(extname,"TESTINFO");
	if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, 0, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(dtcObject,&totalrecords, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName) + " (TESTINFO HDU)";
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	IOData obj;

	// Creating FILDER Column
	obj.inObject = dtcObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"TESTINFO");
	obj.iniRow = totalrecords+1;
	obj.endRow = totalrecords+1;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit," ");

	gsl_matrix *matrixToWrite = gsl_matrix_alloc(1,recordDERIVATIVE->size);
	strcpy(obj.nameCol,"FILDER");
	gsl_matrix_set_row(matrixToWrite,0,recordDERIVATIVE);
	if (writeFitsComplex(obj, matrixToWrite))
	{
		message = "Cannot run routine writeFitsComplex for recordDERIVATIVE (TESTINFO HDU)";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating THRESHOLD Column
	strcpy(obj.nameCol,"Threshold");
	gsl_vector *vectorToWrite = gsl_vector_alloc(1);
	gsl_vector_set(vectorToWrite,0,threshold);
	if (writeFitsSimple(obj, vectorToWrite))
	{
		message = "Cannot run routine writeFitsSimple for threshold (TESTINFO HDU)";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	gsl_matrix_free(matrixToWrite);
	gsl_vector_free(vectorToWrite);

	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	return (EPOK);
}
/*xxxx end of SECTION A7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION  A8 ************************************************************
* calculateTemplate function: This function calculates the template (PULSE column in the library) of non piled-up pulses.
*                             Just in case in the detection process some piled-up pulses have not been distinguished as different pulses => Build a pulseheights histogram.
*                             Use the pulseheights histogram (built by using the PHEIGGHT column of the library), 'Tstart' and 'quality' to select the non piled-up pulses.
*
* - Declare and initialize variables
* - Before building the histogram, select the pulseheihts of the pulses which are enough separated from others and whose 'quality' is 0
* - Create the pulseheights histogram
* - Calculate the pulseaverage only taking into account the valid pulses
* 	- Check if the pulse is piled-up or not
* 	- Non piled-up pulses => Align and average them (currently 'align' is not used because the starting time of the pulses are provided as input parameters)
* //- Just in case due to the noise influence in the alignment, the first sample of the pulseaverage is not around the baseline but higher => Correction
* - Calculate covariance and weight matrices
* - Free allocate of GSL vectors
*
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'pulse_length' and 'EnergyMethod'
* - pulsesAll: Found pulses in the previous records
* - pulsesInRecord: Found pulses in the current record
* - samprate: Sampling rate
* - pulseaverage: Pulseaverage (template) of the non piled-up pulses
* - pulseaverageHeight: Height of the pulseaverage
* - covariance: Covariance matrix
* - weight: Weight matrix
******************************************************************************/
int calculateTemplate (ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight, gsl_matrix **covariance, gsl_matrix **weight)
{
	string message = "";

	// Declare and initialize variables
	int totalPulses = pulsesAll->ndetpulses + pulsesInRecord->ndetpulses;
	gsl_vector *tstart = gsl_vector_alloc(totalPulses);
	gsl_vector *pulseheight = gsl_vector_alloc(totalPulses);
	gsl_vector *quality = gsl_vector_alloc(totalPulses);
	for (int i=0;i<pulsesAll->ndetpulses;i++)
	{
		gsl_vector_set(tstart,i,pulsesAll->pulses_detected[i].Tstart);
		gsl_vector_set(pulseheight,i,pulsesAll->pulses_detected[i].pulse_height);
		gsl_vector_set(quality,i,pulsesAll->pulses_detected[i].quality);
	}
	for (int i=0;i<pulsesInRecord->ndetpulses;i++)
	{
		gsl_vector_set(tstart,i+pulsesAll->ndetpulses,pulsesInRecord->pulses_detected[i].Tstart);
		gsl_vector_set(pulseheight,i+pulsesAll->ndetpulses,pulsesInRecord->pulses_detected[i].pulse_height);
		gsl_vector_set(quality,i+pulsesAll->ndetpulses,pulsesInRecord->pulses_detected[i].quality);
	}

	gsl_vector *nonpileup = gsl_vector_alloc(totalPulses);	// Piled-up pulse => Not taken into account to calculate the template
	long nonpileupPulses = totalPulses;			// A priori, all the found pulses are considered as non piled-up
	gsl_vector_set_all(nonpileup,1);

	int nBins;						// Square-root choice (used by Excel and many others)
	gsl_vector *xhisto;					// X-axis of the pulseheights histogram
	gsl_vector *yhisto;					// Y-axis of the pulseheights histogram
	int index_maximumpulseheight;				// Index where the maximum of the pulseheights histogram is
	double maximumpulseheight;				// Maximum of the pulseheights histogram

	bool firstnonpileupPulse = true;
	gsl_vector *pulse = gsl_vector_alloc(reconstruct_init->pulse_length);
	//gsl_vector *pulse = gsl_vector_alloc(reconstruct_init->pulse_length+safetyMargin);
	//gsl_vector *pulseSHORT = gsl_vector_alloc(reconstruct_init->pulse_length);

	double tstartnext;

	gsl_vector_view temp;

	gsl_vector_scale(tstart,samprate); 			//tstarts not in seconds but in samples
	
	//gsl_vector *pulseaverageSM = gsl_vector_alloc((*pulseaverage)->size+safetyMargin);

	// Before building the histogram, select the pulseheihts of the pulses which are enough separated from others and whose 'quality' is 0
	gsl_vector *pulseheightAUX = gsl_vector_alloc(totalPulses);
	int cnt = 0;
	for (int i=0;i<totalPulses;i++)
	{
		if (i == totalPulses-1)		tstartnext = gsl_vector_get(tstart,i)+2*reconstruct_init->pulse_length;
		else				tstartnext = gsl_vector_get(tstart,i+1);

		//if ((tstartnext-gsl_vector_get(tstart,i) > reconstruct_init->pulse_length) && (gsl_vector_get(quality,i) == 0))
		if ((tstartnext-gsl_vector_get(tstart,i) > reconstruct_init->pulse_length) && (gsl_vector_get(quality,i) != 1) && (gsl_vector_get(quality,i) != 11))
		{
			 gsl_vector_set(pulseheightAUX,cnt,gsl_vector_get(pulseheight,i));
			 cnt = cnt +1;
		}
	}
	temp = gsl_vector_subvector(pulseheightAUX,0,cnt);
	gsl_vector *pulseheightAUX2 = gsl_vector_alloc(cnt);
	gsl_vector_memcpy(pulseheightAUX2,&temp.vector);
	gsl_vector_free(pulseheightAUX);

	// Create the pulseheights histogram
	nBins = floor(sqrt(cnt));
	xhisto = gsl_vector_alloc(nBins);	// X-axis of the pulseheights histogram
	yhisto = gsl_vector_alloc(nBins);	// Y-axis of the pulseheights histogram
	if (createHisto(pulseheightAUX2, nBins, &xhisto, &yhisto))
	{
	    message = "Cannot run createHisto routine";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(pulseheightAUX2);

	index_maximumpulseheight = gsl_vector_max_index(yhisto);
	maximumpulseheight = gsl_vector_get(xhisto,index_maximumpulseheight);

	// Calculate the pulseaverage only taking into account the valid pulses
	gsl_vector_set_all(*pulseaverage,0.0);
	//gsl_vector_set_all(pulseaverageSM,0.0);
	for (int i=0;i<totalPulses;i++)
	{
		if (i == totalPulses-1)		tstartnext = gsl_vector_get(tstart,i)+2*reconstruct_init->pulse_length;
		else 				tstartnext = gsl_vector_get(tstart,i+1);

		// Check if the pulse is piled-up or not
		if ((gsl_vector_get(pulseheight,i) < maximumpulseheight-0.1*maximumpulseheight) || (gsl_vector_get(pulseheight,i) > maximumpulseheight+0.1*maximumpulseheight)
		//	|| (tstartnext-gsl_vector_get(tstart,i) <= reconstruct_init->pulse_length) || (gsl_vector_get(quality,i) >= 1))
			|| (tstartnext-gsl_vector_get(tstart,i) <= reconstruct_init->pulse_length) || ((gsl_vector_get(quality,i) != 0) && (gsl_vector_get(quality,i) != 10)))
		{
 			gsl_vector_set(nonpileup,i,0);
			nonpileupPulses --;
		}
		else
		{
			if (i < pulsesAll->ndetpulses)
			{
				gsl_vector_memcpy(pulse,pulsesAll->pulses_detected[i].pulse_adc);
			}
			else
			{
				gsl_vector_memcpy(pulse,pulsesInRecord->pulses_detected[i-pulsesAll->ndetpulses].pulse_adc);
			}
			//temp = gsl_vector_subvector(pulse,50,reconstruct_init->pulse_length);
			//gsl_vector_memcpy(pulseSHORT,&temp.vector);

			// Non piled-up pulses => Align and average them
			if (firstnonpileupPulse == true)
			{
				gsl_vector_memcpy(*pulseaverage,pulse);
				//gsl_vector_memcpy(pulseaverageSM,pulse);
				//gsl_vector_memcpy(*pulseaverage,pulseSHORT);
				*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
			}
			else
			{
				/*if (align(samprate, pulseaverage,&pulse))
				//if (align(samprate, &pulseaverageSM,&pulse))
				{
					message = "Cannot run align for pulse " + boost::lexical_cast<std::string>(i) + " when 1st pulse is piled-up";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}*/
				/*if (align(samprate, pulseaverage,&pulseSHORT))
				{
					message = "Cannot run align for pulse " + boost::lexical_cast<std::string>(i) + " when 1st pulse is piled-up";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}*/
				gsl_vector_add(*pulseaverage,pulse);
				//gsl_vector_add(pulseaverageSM,pulse);
				*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
			}
			if (firstnonpileupPulse == true)	firstnonpileupPulse = false;
		}
	}

	gsl_vector_scale(*pulseaverage,1.0/(nonpileupPulses));
	//gsl_vector_scale(pulseaverageSM,1.0/(nonpileupPulses));
	
	// In order to delete the samples of the safetyMargin or less (the first one must be buried in the baseline)
	/*double noiseStd = reconstruct_init->noise_spectrum->noiseStd;
	double noiseBsln = reconstruct_init->baseline;
	int numSamplesToBeDeleted = 0;
	int index=gsl_vector_max_index(pulseaverageSM);
	numSamplesToBeDeleted = index;
	//temp = gsl_vector_subvector(pulseaverageSM,numSamplesToBeDeleted,reconstruct_init->pulse_length);
	temp = gsl_vector_subvector(pulseaverageSM,0,reconstruct_init->pulse_length);
	gsl_vector_memcpy(*pulseaverage,&temp.vector);*/

	// Calculate covariance and weight matrix
	bool saturatedPulses = false;
	if (pulsesInRecord->pulses_detected[0].quality >= 10)		saturatedPulses = true;
	if (weightMatrix(reconstruct_init, saturatedPulses, pulsesAll, pulsesInRecord, nonpileupPulses, nonpileup, *pulseaverage, covariance, weight))
	{
		message = "Cannot run weightMatrix routine";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	
	//cout<<"Antes de eigenVV"<<endl;
	//gsl_vector *eigenvalues = gsl_vector_alloc((*covariance)->size1);
	//gsl_matrix *eigenvectors = gsl_matrix_alloc((*covariance)->size1,(*covariance)->size2);
	/*gsl_vector *eigenvalues;
	gsl_matrix *eigenvectors;
	if (eigenVV(*covariance,&eigenvectors,&eigenvalues))
	{
		message = "Cannot run eigenVV routine";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(eigenvalues);
	gsl_matrix_free(eigenvectors);*/
	//cout<<"Despues de eigenVV"<<endl;
	
	//Just in case due to the noise influence in the alignment, the first sample of the pulseaverage is not around the baseline but higher
	/*double meanLast200points, sgLast200points;
	temp = gsl_vector_subvector(*pulseaverage,reconstruct_init->pulse_length-200-1,200);
	if (findMeanSigma (&temp.vector, &meanLast200points, &sgLast200points))
	{
		message = "Cannot run findMeanSigma routine for kappa-sigma iteration";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	if (fabs(gsl_vector_get(*pulseaverage,0))>fabs(meanLast200points)+3*sgLast200points)
	{
		gsl_vector_set(*pulseaverage,0,meanLast200points);
	}*/

	*pulseaverageHeight = *pulseaverageHeight/nonpileupPulses;

	// Free allocate of GSL vectors
	gsl_vector_free(tstart);
	gsl_vector_free(pulseheight);
	gsl_vector_free(quality);
	gsl_vector_free(nonpileup);
	gsl_vector_free(xhisto);
	gsl_vector_free(yhisto);
	gsl_vector_free(pulse);
	//gsl_vector_free(pulseaverageSM);

	return (EPOK);
}
/*xxxx end of SECTION A8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A9 ************************************************************
* createHisto function: This function builds the histogram of the input vector.
*                       Histogram x-axis values are the different input vector values (pulseheights).
*                       Histogram y-axis values are the the number of cases per unit of the variable on the horizontal axis.
*
* - Declare variables
* - It is only going to work with the positive elements of the input vector -> 'invectoraux2'
* - Check if all the values of 'invector' are the same => Histogram of only one bin
* - Obtain 'invector_max 'and 'invector_min'
* - Obtain 'binSize'
* - Create histogram axis
* - Free allocate of GSL vectors
*
* Parameters:
* - invector: Input vector
* - nbins: Number of bins to build the histogram
* - xhistogsl: Output histogram x-axis
* - yhistogsl: Output histogram y-axis
******************************************************************************/
int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhistogsl, gsl_vector **yhistogsl)
{
	string message = "";

	// Declare variables
	int size = invector->size;
	double invectormax= 0;      				// Maximum of 'invector'
	double invectormin=1e10;  				// Minimum of 'invector'
	double binSize;						// Size in samples of each bin
	int ind = 0;                				// Index of the bin which contains each 'invector' element
	gsl_vector *invectoraux = gsl_vector_alloc(size);	// Auxiliary variable
	gsl_vector *invectoraux2;				// Auxiliary variable
	gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)

	// It is only going to work with the positive elements of the input vector -> 'invectoraux2'
	for (int i=0;i<size;i++)
	{
		if (gsl_vector_get(invector,i) > 0)
		{
			gsl_vector_set(invectoraux,ind,gsl_vector_get(invector,i));
			ind = ind+1;
		}
	}
	temp = gsl_vector_subvector(invectoraux,0,ind);
	invectoraux2 = gsl_vector_alloc(ind);
	gsl_vector_memcpy(invectoraux2, &temp.vector);
	size = invectoraux2->size;
	gsl_vector_free(invectoraux);

	// To check if all the values of 'invector' are the same
	// For example, high energies and no noise => Saturated pulses
	invectoraux = gsl_vector_alloc(size);
	gsl_vector *invectoraux1 = gsl_vector_alloc(size);
	gsl_vector_memcpy(invectoraux,invectoraux2);
	gsl_vector_set_all(invectoraux1,gsl_vector_get(invectoraux2,0));
	gsl_vector_scale(invectoraux1,-1.0);
	gsl_vector_add(invectoraux,invectoraux1);
	gsl_vector_free(invectoraux1);

	if (gsl_vector_isnull(invectoraux) == 1) // All the values are the same
	{
		gsl_vector_free(*xhistogsl);
		gsl_vector_free(*yhistogsl);
		*xhistogsl = gsl_vector_alloc(1);
		*yhistogsl = gsl_vector_alloc(1);
		gsl_vector_set(*xhistogsl,0,gsl_vector_get(invectoraux2,0));
		gsl_vector_set(*yhistogsl,0,1); //It doesn't matter the value because the maximum is going to look for
	}
	else
	{
		// Obtain 'invector_max'
		for (int i=0; i<size; i++)
		{
			if (invectormax < gsl_vector_get (invectoraux2,i))	invectormax = gsl_vector_get (invectoraux2,i);
		}
		// Obtain 'invector_min'
		for (int i=0; i<size; i++)
		{
			if (invectormin > gsl_vector_get (invectoraux2,i))	invectormin = gsl_vector_get (invectoraux2,i);
		}

		if (invectormax == invectormin)
		{
			message = "Invectormax == invectormin";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		if ((invectormax != 0) && (invectormin != 1e10))
		{
			// Obtain binSize
			binSize = (invectormax - invectormin) / nbins;        // Size of each bin of the histogram

			// Create histogram x-axis
			for (int i=0; i<nbins; i++)	{gsl_vector_set (*xhistogsl,i,binSize*i+invectormin+binSize/2);}

			// Create histogram y-axis
			gsl_vector_set_zero (*yhistogsl);
			for (int i=0; i<size; i++)
			{
				ind = (int) ((gsl_vector_get(invectoraux2,i) - invectormin) / binSize);
				if (ind == nbins) ind--;
				gsl_vector_set (*yhistogsl,ind, gsl_vector_get(*yhistogsl,ind) + 1);
			}
		}
	}

	// Free allocate of GSL vectors
	gsl_vector_free(invectoraux);
	gsl_vector_free(invectoraux2);

	return EPOK;
}
/*xxxx end of SECTION A9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A10 ************************************************************
* align function: This function aligns 'vector1' with 'vector2' (by delaying or moving forward 'vector2').
*                 It is supposed that 'vector1' and 'vector2' are shifted replicas of the same function.
*
* From the discrete function x[n] (n=0,...,N-1 => Length = N) and according the time shifting property of the Fourier transform:
*
*  x[n]   <------> X[f]
*  x[n-m] <------> X[f]·exp(-j2·pi·m/N)
*
*  Shift = m => Phase due to the shift = 2·pi·m/N => m = Phase due to the shift·N/(2·pi)
*
* - Declare variables
* - FFT of 'vector1'
* - FFT of 'vector2'
* - Phases of the FFT_vector1 and FFT_vector2, *size/(2*pi)
* - Shift between the input vectors
* - 'shiftdouble' into 'shiftint' (because we are working with samples)
* - Move forward or delay 'vector2' depending on positive or negative shift
*
* Parameters:
* - samprate: Sampling rate
* - vector1: Input vector
* - vector2: Input vector which is delayed or moved forward to be aligned with 'vector1'
******************************************************************************/
int align(double samprate, gsl_vector **vector1, gsl_vector **vector2)
{
	const double pi = 4.0 * atan(1.0);
	string message = "";

	// Declare variables
	int size = (*vector1)->size;

	double SelectedTimeDuration = size/samprate;
	gsl_vector_complex *vector1fft = gsl_vector_complex_alloc(size);
	gsl_vector_complex *vector2fft = gsl_vector_complex_alloc(size);
	double vector1fft_ph;
	double vector2fft_ph;

	double shiftdouble;
	int shiftint;

	gsl_vector *vector2shifted = gsl_vector_alloc(size);

	// FFT of 'vector1'
	if (FFT(*vector1,vector1fft,SelectedTimeDuration))
	{
		message = "Cannot run FFT routine for vector1";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// FFT of 'vector2'
	if (FFT(*vector2,vector2fft,SelectedTimeDuration))
	{
		message = "Cannot run FFT routine for vector2";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Phases of the FFT_vector1 and FFT_vector2, *size/(2*pi)
	vector1fft_ph= gsl_complex_arg(gsl_vector_complex_get(vector1fft,1))*size/(2*pi);
	vector2fft_ph= gsl_complex_arg(gsl_vector_complex_get(vector2fft,1))*size/(2*pi);
	gsl_vector_complex_free(vector1fft);
	gsl_vector_complex_free(vector2fft);

	// Shift between the input vectors
	shiftdouble = vector1fft_ph-vector2fft_ph;

	// 'shiftdouble' into 'shiftint' (because we are working with samples)
	if ((shiftdouble > -1) && (shiftdouble < 1)) 	shiftint = 0;
	else if (shiftdouble > 1)			shiftint = floor(shiftdouble);
	else if (shiftdouble < -1)			shiftint = ceil(shiftdouble);

	// Move forward or delay vector2 depending on positive or negative shift
	if (shiftint > 0)
	{
		if (shift_m(*vector2,vector2shifted,shiftint))
		{
			message = "Cannot run shift_m routine for vector2";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		gsl_vector_memcpy(*vector2,vector2shifted);
	}
	else if (shiftint < 0)
	{
		if (shiftm(*vector2,vector2shifted,fabs(shiftint)))
		{
			message = "Cannot run shiftm routine for vector2";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		gsl_vector_memcpy(*vector2,vector2shifted);
	}

	gsl_vector_free(vector2shifted);

	return (EPOK);
}
/*xxxx end of SECTION A10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A11 ************************************************************
* shiftm function: This function returns (in 'vectorout') the 'vectorin' delayed m samples.
*
* Parameters:
* - vectorin: Input vector
* - vectorout: Output vector which is the 'vectorin' delayed m samples
******************************************************************************/
int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	int size = vectorin->size;

	for (int i=0;i<size-m;i++)
	{
		gsl_vector_set(vectorout,i+m,gsl_vector_get(vectorin,i));
	}
	for (int i=0;i<m;i++)
	{
		gsl_vector_set(vectorout,i,gsl_vector_get(vectorin,0));
	}

	return (EPOK);
}
/*xxxx end of SECTION A11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A12 ************************************************************
* shift_m function: This function returns (in 'vectorout') the 'vectorin' moved forward m samples.
*
* Parameters:
* - vectorin: Input vector
* - vectorout: Output vector which is the 'vectorin' moved forward m samples
******************************************************************************/
int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	int size = vectorin->size;

	for (int i=m;i<size;i++)
	{
		gsl_vector_set(vectorout,i-m,gsl_vector_get(vectorin,i));
	}
	for (int i=size-m;i<size;i++)
	{
		gsl_vector_set(vectorout,i,gsl_vector_get(vectorin,size-1));
	}

	return (EPOK);
}
/*xxxx end of SECTION A12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A13 ************************************************************
* weightMatrix function: This function calculates the weight matrix by using the non piled-up pulses found in all the records, stored
*                        in 'pulsesAll' (previous records) and 'pulsesInRecord' (current record). The weight matrix of each energy
*                        (and other intermediate values) will be stored in the library by the function 'fillInLibraryData'.
*
* Si^p: Value of the i-sample of the pulse p
* Mi^p: Value of the i-sample of the model p (model='pulseaverage')  Mi = <Si> = (1/N)sum(p=1,N){Si^p}
* N: number of non piled-up pulses
* Di = Si - Mi
* <DiDj> = E[(Si-Mi)(Sj-Mj)] = (1/N)sum(p=1,N){(Si^p-Mi^p)(Sj^p-Mj^p)}
* 		  |<D1D1> <D1D2>...<D1Dn>|
*  Vij = <DiDj> = |<D2D1> <D2D2>...<D2Dn>|	where n is the pulse length     V => (nxn)
*                 |...                   |
*                 |<DnD1> <DnD2>...<DnDn>|
*  W = 1/V
*
* - Calculate the elements of the diagonal of the covariance matrix
* - Calculate the elements out of the diagonal of the covariance matrix
* - If saturated pulses => Covariance matrix is a singular matrix => Non invertible 
*   In order to allow the covariance matrix to be inverted => Replacing 0's (0's are due to the saturated values, equal in the pulse and in the model)
* 	- Elements of the diagonal: Generating a random double f1 between a range fMin,fMax (-NoiseStd,NoiseStd) to replace 0's with f1*f1 
*       - Elements out of the diagonal: Generating two random doubles f1 and f2 between a range fMin,fMax (-NoiseStd,NoiseStd) to replace 0's with f1*f2 
* - Calculate the weight matrix
*
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'pulse_length'
* - saturatedPulses: If 'true', all the pulses (CALIBRATION => all the pulses have the same energy) are saturated
* - pulsesAll: Found pulses in the previous records
* - pulsesInRecord: Found pulses in the current record
* - nonpileupPulses: Number of non piled-up pulses
* - nonpileup: Vector containing info about all the pulses informing if they are piled-up or not
* - pulseaverage: Average=Model=Template of all the found pulses
* - covariance: Covariance matrix
* - weight: Weight matrix
******************************************************************************/
int weightMatrix (ReconstructInitSIRENA *reconstruct_init, bool saturatedPulses, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, long nonpileupPulses, gsl_vector *nonpileup, gsl_vector *pulseaverage, gsl_matrix **covariance, gsl_matrix **weight)
{
	double elementValue = 0.0;
	gsl_permutation *perm = gsl_permutation_alloc(reconstruct_init->pulse_length);
	int s=0;

	// Elements of the diagonal of the covariance matrix
	for (int i=0;i<reconstruct_init->pulse_length;i++)
	{
		for (int p=0;p<pulsesAll->ndetpulses;p++)
		{
			if (gsl_vector_get(nonpileup,p) == 1)
			{
				elementValue = elementValue +
					pow(gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i),2.0);
			}
		}
		for (int p=0;p<pulsesInRecord->ndetpulses;p++)
		{
			if (gsl_vector_get(nonpileup,pulsesAll->ndetpulses+p) == 1)
			{
				elementValue = elementValue +
					pow(gsl_vector_get(pulsesInRecord->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i),2.0);
			}
		}
		elementValue = elementValue/nonpileupPulses;
		
		//cout<<i<<" "<<elementValue<<endl;

		gsl_matrix_set(*covariance,i,i,elementValue);

		elementValue = 0.0;
	}

	// Other elements
	for (int i=0;i<reconstruct_init->pulse_length;i++)
	{
		for (int j=i+1;j<reconstruct_init->pulse_length;j++)
		{
			for (int p=0;p<pulsesAll->ndetpulses;p++)
			{
				if (gsl_vector_get(nonpileup,p) == 1)
				{
					elementValue = elementValue +
						(gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i))*
						(gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,j)-gsl_vector_get(pulseaverage,j));	
				}
			}
			for (int p=0;p<pulsesInRecord->ndetpulses;p++)
			{
				if (gsl_vector_get(nonpileup,pulsesAll->ndetpulses+p) == 1)
				{
					elementValue = elementValue +
						(gsl_vector_get(pulsesInRecord->pulses_detected[p].pulse_adc,i)-gsl_vector_get(pulseaverage,i))*
						(gsl_vector_get(pulsesInRecord->pulses_detected[p].pulse_adc,j)-gsl_vector_get(pulseaverage,j));
				}
			}
			elementValue = elementValue/nonpileupPulses;

			gsl_matrix_set(*covariance,i,j,elementValue);
			gsl_matrix_set(*covariance,j,i,elementValue);

			elementValue = 0.0;
		}
	}
	
	// If saturated pulses => Covariance matrix is a singular matrix => Non invertible 
	// In order to allow the covariance matrix to be inverted => Replacing 0's (0's are due to the saturated values, equal in the pulse and in the model)
	//    - Elements of the diagonal: Generating a random double f1 between a range fMin,fMax (-NoiseStd,NoiseStd) to replace 0's with f1*f1 
	//    - Elements out of the diagonal: Generating two random doubles f1 and f2 between a range fMin,fMax (-NoiseStd,NoiseStd) to replace 0's with f1*f2
	int cnt = 0;
	if (saturatedPulses == true)
	{
		//double noiseStd = 8.5195838175;
		double noiseStd = reconstruct_init->noise_spectrum->noiseStd;
		double f1 = 0;
		double f2 = 0;
		double fMin = (-1.0)*noiseStd;
		double fMax = noiseStd;
			
		for (int i=0;i<reconstruct_init->pulse_length;i++)
		{
			//for (int j=0;j<reconstruct_init->pulse_length;j++)
			for (int j=i;j<reconstruct_init->pulse_length;j++)	// In order to have a symmetric matrix (j=i)
			{
				if (fabs(gsl_matrix_get(*covariance,i,j)) < 1e-24)	// Meaning 'equal to 0'
				{	
					cnt++;
					srand(time(NULL)+i+j+cnt);	// In order to change each run time the seed
					f1 = (double)rand()/RAND_MAX; 
					f1 = fMin + f1 * (fMax - fMin);
					if (i == j)
					{
						srand(time(NULL)+i+j+cnt+1);	// In order to get two different numbers (f1 and f2)
						f2 = (double)rand()/RAND_MAX; 
						f2 = fMin + f2 * (fMax - fMin);
						gsl_matrix_set(*covariance,i,j,f1*f2);
					}
					else	
					{
						gsl_matrix_set(*covariance,i,j,f1*f1);	
						gsl_matrix_set(*covariance,j,i,f1*f1);	// In order to have a symmetric matrix
					}
					//cout<<"i+j+cnt= "<<i+j+cnt<<" f1="<<f1<<" f2="<<f2<<endl;
				}
				//if (i == j) cout<<i<<" "<<i+j+cnt<<" "<<gsl_matrix_get(*covariance,i,j)<<" f1="<<f1<<" f2="<<f2<<" f1*f2="<<f1*f2<<endl;
				//if (i == j) cout<<i<<" "<<gsl_matrix_get(*covariance,i,j)<<endl;
			}
		}
	}
	
	// Calculate the weight matrix
	gsl_matrix *covarianceaux = gsl_matrix_alloc((*covariance)->size1,(*covariance)->size2);
	gsl_matrix_memcpy(covarianceaux,*covariance);
	gsl_linalg_LU_decomp(covarianceaux, perm, &s);
	gsl_linalg_LU_invert(covarianceaux, perm, *weight);
	gsl_matrix_free(covarianceaux);
	/*cout<<gsl_matrix_get(*covariance,0,0)<<" "<<gsl_matrix_get(*covariance,0,1)<<" "<<gsl_matrix_get(*covariance,0,2)<<" "<<gsl_matrix_get(*covariance,0,3)<<"..."<<gsl_matrix_get(*covariance,0,1023)<<endl;
	cout<<gsl_matrix_get(*covariance,1,0)<<" "<<gsl_matrix_get(*covariance,1,1)<<" "<<gsl_matrix_get(*covariance,1,2)<<" "<<gsl_matrix_get(*covariance,1,3)<<"..."<<gsl_matrix_get(*covariance,1,1023)<<endl;
	cout<<gsl_matrix_get(*weight,0,0)<<" "<<gsl_matrix_get(*weight,0,1)<<" "<<gsl_matrix_get(*weight,0,2)<<" "<<gsl_matrix_get(*weight,0,3)<<"..."<<gsl_matrix_get(*weight,0,1023)<<endl;
	cout<<gsl_matrix_get(*weight,1,0)<<" "<<gsl_matrix_get(*weight,1,1)<<" "<<gsl_matrix_get(*weight,1,2)<<" "<<gsl_matrix_get(*weight,1,3)<<"..."<<gsl_matrix_get(*weight,1,1023)<<endl;*/

	gsl_permutation_free(perm);

	return (EPOK);
}
/*xxxx end of SECTION A13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION Axx ************************************************************
* eigenVV function: This function provides the principal eigenvectors and eigenvalues of the input matrix
*                   The eigenvalues and eigenvectors are ordered and only the principal components are provided
* - 
* 
* Parameters:
* - 
******************************************************************************/
int eigenVV (gsl_matrix *matrixin, gsl_matrix **eigenvectors, gsl_vector **eigenvalues)
{
	int status = EPOK;
	
	gsl_vector *eigenvaluesAll = gsl_vector_alloc(matrixin->size1);
	gsl_matrix *eigenvectorsAll = gsl_matrix_alloc(matrixin->size1,matrixin->size2);
	
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (matrixin->size1);
  
	/*cout<<eigenvalues->size<<endl;
	cout<<eigenvectors->size1<<" "<<eigenvectors->size2<<endl;
	cout<<matrixin->size1<<" "<<matrixin->size2<<endl;*/
	gsl_eigen_symmv (matrixin, eigenvaluesAll, eigenvectorsAll, w);

	gsl_eigen_symmv_free (w);

	gsl_eigen_symmv_sort (eigenvaluesAll, eigenvectorsAll, GSL_EIGEN_SORT_ABS_ASC);
	
	// Choose eigenvectors whose abs(eigenvalues) are grater than 1
	int indexToStartToTakeAccount;
	for (int i = 0; i < eigenvaluesAll->size; i++)
	{
		if (abs(gsl_vector_get(eigenvaluesAll,i)) > 1.0)
		{
			indexToStartToTakeAccount = i;
			break;
		}
	}
	gsl_vector_view temp;
	*eigenvalues = gsl_vector_alloc(eigenvaluesAll->size-indexToStartToTakeAccount);
	*eigenvectors = gsl_matrix_alloc(eigenvectorsAll->size1,(*eigenvalues)->size);
	temp = gsl_vector_subvector(eigenvaluesAll,indexToStartToTakeAccount,eigenvaluesAll->size-indexToStartToTakeAccount);
	
	gsl_vector_memcpy(*eigenvalues,&temp.vector);
	

	/*gsl_vector *column = gsl_vector_alloc((*eigenvectors)->size1);
	for (int i = 0; i < matrixin->size1; i++)
	{
		cout<<"eigenvvalue("<<i<<"): "<<gsl_vector_get(*eigenvalues,i)<<endl;
		//gsl_matrix_get_col(column,eigenvectors,i);
		//for (int j = 0; j < eigenvectors->size1; j++)
		//{
		//	cout<<"eigenvector("<<i<<"): "<<gsl_vector_get(column,j)<<endl;
		//}
	}
	gsl_vector_free(column);*/
	
	gsl_vector_free(eigenvaluesAll);
	gsl_matrix_free(eigenvectorsAll);
	
	return (EPOK);  
}
/*xxxx end of SECTION Axx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A14 ************************************************************
* writeLibrary function: This function writes the library (sorting if it is necesary and calculating some intermediate parameters)
*
* - Adding a new row to the library
* - Write the first row of the library
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
* - samprate: Sampling rate
* - estenergy: Pulse height of the template whose energy is going to be added to the library
* - pulsetemplate: Pulse template whose energy is going to be added to the library
* - covariance: Covariance matrix of the energy which is going to be added to the library
* - weight: Weight matrix of the energy which is going to be added to the library
* - appendToLibrary: 'true' if adding a new row to the library
*                    'false' if it is the first row
* - inLibObject: Object which contains information of the library FITS file 
******************************************************************************/
int writeLibrary(ReconstructInitSIRENA *reconstruct_init, double samprate, double estenergy, gsl_vector *pulsetemplate, gsl_matrix *covariance, gsl_matrix *weight, bool appendToLibrary, fitsfile **inLibObject)
{
	int status = EPOK;
	string message = "";

	char inLibName[256];
	strncpy(inLibName, reconstruct_init->library_file,255);
	inLibName[255]='\0';
	
	char keyname[10];
	char extname[20];
	//IOData obj;
	
	int runF0orB0val;
	if (strcmp(reconstruct_init->FilterMethod,"F0") == 0)		// Deleting the frequency-zero bin
	{
		runF0orB0val = 0;
	}
	else if (strcmp(reconstruct_init->FilterMethod,"B0") == 0)	// Working without baseline
	{
		runF0orB0val = 1;
	}

	if (appendToLibrary == true)
	{
		long eventcntLib;
		if (fits_get_num_rows(*inLibObject,&eventcntLib, &status))
		{
			message = "Cannot get number of rows in " + string(inLibName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		assert(eventcntLib > 0);
		long eventcntLib1 = eventcntLib + 1;
		strcpy(keyname,"EVENTCNT");
		
		if (fits_update_key(*inLibObject,TLONG,keyname, &eventcntLib1,NULL,&status))
		{
		    message = "Cannot update keyword " + string(keyname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		
		if (readAddSortParams(reconstruct_init,inLibObject,samprate,eventcntLib,estenergy,pulsetemplate,covariance,weight))
		{
			message = "Cannot run routine readAddSortParams in writeLibrary";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		// Primary HDU
                strcpy(extname,"Primary");
                if (fits_movabs_hdu(*inLibObject, 1, NULL, &status))
                {
                        message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                }

                char keyvalstr[1000];

                char str_procnumber[125];               snprintf(str_procnumber,125,"%d",eventcntLib);
                string strprocname (string("PROC") + string(str_procnumber));
                strcpy(keyname,strprocname.c_str());
                string strprocval (string("PROC") + string(str_procnumber) + string(" Starting parameter list"));
                strcpy(keyvalstr,strprocval.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

                string strproc (string("RecordFile = ") + reconstruct_init->record_file);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);

                strproc=string("TesEventFile = ") + reconstruct_init->event_file;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("LibraryFile = ") + reconstruct_init->library_file;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("NoiseFile = ") + reconstruct_init->noise_file;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_mode[125];                     snprintf(str_mode,125,"%d",reconstruct_init->mode);
                strproc = string("mode = ") + string(str_mode);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("PixelType = ") + reconstruct_init->PixelType;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("FilterDomain = ") + reconstruct_init->FilterDomain;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("FilterMethod = ") + reconstruct_init->FilterMethod;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_LagsOrNot[125];                snprintf(str_LagsOrNot,125,"%d",reconstruct_init->LagsOrNot);
                strproc=string("LagsOrNot = ") + string(str_LagsOrNot);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_OFIter[125];                   snprintf(str_OFIter,125,"%d",reconstruct_init->OFIter);
                strproc=string("OFIter = ") + string(str_OFIter);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_OFLib[125];                    snprintf(str_OFLib,125,"%d",reconstruct_init->OFLib);
                strproc=string("OFLib = ") + string(str_OFLib);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
                               
                strproc=string("OFInterp = ") + reconstruct_init->OFInterp;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("OFStrategy = ") + reconstruct_init->OFStrategy;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
                               
                char str_OFLength[125];                 snprintf(str_OFLength,125,"%d",reconstruct_init->OFLength);
                strproc=string("OFLength = ") + string(str_OFLength);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_maxPulsesPerRecord[125];       snprintf(str_maxPulsesPerRecord,125,"%d",reconstruct_init->maxPulsesPerRecord);
                strproc=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_pulse_length[125];             snprintf(str_pulse_length,125,"%d",reconstruct_init->pulse_length);
                strproc=string("PulseLength = ") + string(str_pulse_length);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_tauFall[125];                  snprintf(str_tauFall,125,"%e",reconstruct_init->tauFall);
                strproc=string("tauFall = ") + string(str_tauFall);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_scaleFactor[125];              snprintf(str_scaleFactor,125,"%f",reconstruct_init->scaleFactor);
                strproc=string("scaleFactor = ") + string(str_scaleFactor);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_samplesUp[125];                snprintf(str_samplesUp,125,"%f",reconstruct_init->samplesUp);
                strproc=string("samplesUp = ") + string(str_samplesUp);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_nSgms[125];                    snprintf(str_nSgms,125,"%f",reconstruct_init->nSgms);
                strproc=string("nSgms = ") + string(str_nSgms);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_LrsT[125];                     snprintf(str_LrsT,125,"%e",reconstruct_init->LrsT);
                strproc=string("LrsT = ") + string(str_LrsT);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_LbT[125];                      snprintf(str_LbT,125,"%e",reconstruct_init->LbT);
                strproc=string("LbT = ") + string(str_LbT);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_monoenergy[125];               snprintf(str_monoenergy,125,"%f",reconstruct_init->monoenergy);
                strproc=string("monoenergy = ") + string(str_monoenergy);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_intermediate[125];             snprintf(str_intermediate,125,"%d",reconstruct_init->intermediate);
                strproc=string("intermediate = ") + string(str_intermediate);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("detectFile = ") + reconstruct_init->detectFile;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strproc=string("filterFile = ") + reconstruct_init->filterFile;
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_clobber[125];                  snprintf(str_clobber,125,"%d",reconstruct_init->clobber);
                strproc=string("clobber = ") + string(str_clobber);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
		
		char str_SaturationValue[125];          snprintf(str_SaturationValue,125,"%e",reconstruct_init->SaturationValue);
                strproc=string("SaturationValue = ") + string(str_SaturationValue);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_tstartPulse1[125];             snprintf(str_tstartPulse1,125,"%d",reconstruct_init->tstartPulse1);
                strproc=string("tstartPulse1 = ") + string(str_tstartPulse1);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_tstartPulse2[125];             snprintf(str_tstartPulse2,125,"%d",reconstruct_init->tstartPulse2);
                strproc=string("tstartPulse2 = ") + string(str_tstartPulse2);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                char str_tstartPulse3[125];             snprintf(str_tstartPulse3,125,"%d",reconstruct_init->tstartPulse3);
                strproc=string("tstartPulse3 = ") + string(str_tstartPulse3);
                strcpy(keyvalstr,strproc.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                strcpy(keyname,strprocname.c_str());
                strprocval = string("PROC") + string(str_procnumber) + string(" Ending parameter list");
                strcpy(keyvalstr,strprocval.c_str());
                fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,NULL,&status);
               
                if (status != 0)
                {
                        message = "Cannot write some keyword in library file " + string(inLibName);
                        EP_PRINT_ERROR(message,status); return(EPFAIL);
                }
	}
	else
	{
		gsl_vector *energyoutgsl = gsl_vector_alloc(1);
		gsl_vector *estenergyoutgsl = gsl_vector_alloc(1);
		gsl_matrix *pulsetemplates_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
		gsl_matrix *pulsetemplatesb0_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
		gsl_matrix *matchedfilters_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
		gsl_matrix *matchedfiltersb0_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	
		strcpy(extname,"LIBRARY");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, 0, &status))
		{
			message = "Cannot move to HDU  " + string(extname) + " in library";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}

		if (fits_write_key(*inLibObject,TINT,"EVENTSZ",&reconstruct_init->pulse_length,NULL,&status))
		{
			message = "Cannot write keyword EVENTSZ in library";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}

		gsl_vector_set (energyoutgsl,0,reconstruct_init->monoenergy);
		gsl_vector_set (estenergyoutgsl,0,estenergy);

		gsl_matrix_set_row(pulsetemplates_matrix,0,pulsetemplate);
		gsl_matrix_set_row(matchedfilters_matrix,0,pulsetemplate);

		gsl_vector *baselinegsl = gsl_vector_alloc(pulsetemplate->size);
		gsl_vector_set_all(baselinegsl,-1.0*reconstruct_init->noise_spectrum->baseline);
		gsl_vector_add(pulsetemplate,baselinegsl);
		gsl_vector_free(baselinegsl);
		gsl_matrix_set_row(pulsetemplatesb0_matrix,0,pulsetemplate);
		gsl_matrix_set_row(matchedfiltersb0_matrix,0,pulsetemplate);
		
		gsl_matrix_scale(matchedfilters_matrix,1./reconstruct_init->monoenergy);
		gsl_matrix_scale(matchedfiltersb0_matrix,1./reconstruct_init->monoenergy);
		
		if (addFirstRow(reconstruct_init,inLibObject,samprate,runF0orB0val,energyoutgsl,estenergyoutgsl,pulsetemplates_matrix,pulsetemplatesb0_matrix,matchedfilters_matrix, matchedfiltersb0_matrix,covariance,weight))
		{
			message = "Cannot run addFirstRow in writeLibrary";
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		
		// Free allocate of GSL vectors and matrices
		gsl_vector_free(energyoutgsl);
		gsl_vector_free(estenergyoutgsl);
		gsl_matrix_free (pulsetemplates_matrix);
		gsl_matrix_free (pulsetemplatesb0_matrix);
		gsl_matrix_free (matchedfilters_matrix);
		gsl_matrix_free (matchedfiltersb0_matrix);
	}

	if (fits_close_file(*inLibObject,&status))
	{
	    message = "Cannot close file " + string(inLibName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	
	return (EPOK);
}
/*xxxx end of SECTION A14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A15 ************************************************************
* addFirstRow function: This function writes the first row of the library (without intermediate parameters because it would be necessary to have at least two rows=energies in the library)
* 
* Parameters: 
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'mode' and 'noise_spectrum' in order to run 'calculus_optimalFilter'
* - inLibObject: Object which contains information of the library FITS file
* - samprate: Sampling rate
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - E:	First energy to be included in the library
* - PHEIGHT: Pulse height associated to the first energy to be included in the library
* - PULSE: Pulse template associated to the first energy to be included in the library
* - PULSEB0: Pulse template without baseline associated to the first energy to be included in the library
* - MF: Matched filter associated to the first energy to be included in the library
* - MFB0: Matched filter without baseline associated to the first energy to be included in the library
* - COVAR: Covariance matrix associated to the first energy to be included in the library
* - WEIGHT: Weight matrix associated to the first energy to be included in the library
******************************************************************************/
int addFirstRow(ReconstructInitSIRENA *reconstruct_init, fitsfile **inLibObject, double samprate, int runF0orB0val, gsl_vector *E, gsl_vector *PHEIGHT, gsl_matrix *PULSE, gsl_matrix *PULSEB0, gsl_matrix *MF, gsl_matrix *MFB0, gsl_matrix *COVAR, gsl_matrix *WEIGHT)
{
	int status = EPOK;
	string message = "";
	
	gsl_vector *optimalfilter = NULL;
	gsl_vector *optimalfilter_f = NULL;
	gsl_vector *optimalfilter_FFT = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	double normalizationFactor;
	
	IOData obj;
	
	obj.inObject = *inLibObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"LIBRARY");
	obj.iniRow = 1;
	obj.endRow = 1;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	
	// Creating ENERGY Column
	strcpy(obj.nameCol,"ENERGY");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"eV");
	if (writeFitsSimple(obj, E))
	{
		message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating PHEIGHT Column
	strcpy(obj.nameCol,"PHEIGHT");
	strcpy(obj.unit,"ADC");
	if (writeFitsSimple(obj, PHEIGHT))
	{
		message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating PULSE Column
	strcpy(obj.nameCol,"PULSE");
	strcpy(obj.unit,"ADC");
	if (writeFitsComplex(obj, PULSE))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating PULSEB0 Column
	strcpy(obj.nameCol,"PULSEB0");
	strcpy(obj.unit,"ADC");
	if (writeFitsComplex(obj, PULSEB0))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating MF Column
	strcpy(obj.nameCol,"MF");
	strcpy(obj.unit," ");
	//gsl_matrix_scale(MF,1.0/reconstruct_init->monoenergy);
	if (writeFitsComplex(obj, MF))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating MFB0 Column
	strcpy(obj.nameCol,"MFB0");
	strcpy(obj.unit," ");
	//gsl_matrix_scale(MFB0,1.0/reconstruct_init->monoenergy);
	if (writeFitsComplex(obj, MFB0))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	
	strcpy(obj.nameCol,"OF");
	strcpy(obj.unit," ");
	gsl_vector *matchedfilters_row = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix_get_row(matchedfilters_row,MF,0);	//Matched filter
	// Calculate the optimal filter
	if (calculus_optimalFilter (0, 0, reconstruct_init->mode, matchedfilters_row, matchedfilters_row->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex, &normalizationFactor))
	{
		message = "Cannot run routine calculus_optimalFilter in writeLibrary";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_matrix *optimalfilters_matrix = gsl_matrix_alloc(1,optimalfilter->size);
	gsl_matrix_set_row(optimalfilters_matrix,0,optimalfilter);
	if (writeFitsComplex(obj, optimalfilters_matrix))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_matrix_free(optimalfilters_matrix);

	// Creating NFCTR Column
	strcpy(obj.nameCol,"NFCTR");
	strcpy(obj.unit," ");
	gsl_vector *nrmfctrgsl = gsl_vector_alloc(1);
	gsl_vector_set(nrmfctrgsl,0,normalizationFactor);
	if (writeFitsSimple(obj, nrmfctrgsl))
	{
		message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(nrmfctrgsl);
	
	// Creating COVARM Column
	strcpy(obj.nameCol,"COVARM");
	strcpy(obj.unit," ");
	if (writeFitsComplex(obj, COVAR))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating WEIGHTM Column
	strcpy(obj.nameCol,"WEIGHTM");
	strcpy(obj.unit," ");
	if (writeFitsComplex(obj, WEIGHT))
	{
		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;
	
	
	return (EPOK);
}
/*xxxx end of SECTION A15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A16 ************************************************************
* readAddSortParams function: This function read the library data, add new data and sort the data according an ascending energy order
* 
* - Declare variables
* - Load values already in the library
* - Add new values
* - Re-sort
* - Add intermeadiate values
* - Recalculate intermediate values of some new pairs
* - Write values in the library
* - Free allocate of GSL vectors
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'FilterMethod', 'pulse_length', 'library_collection', 'monoenergy', 'mode' and 'noise_spectrum'
* - inLibObject: Object which contains information of the library FITS file
* - samprate: Sampling rate
* - eventcntLib: Number of templates in the library
* - estenergy: Pulse height of the template whose energy is going to be added to the library
* - pulsetemplate: Pulse template whose energy is going to be added to the library
* - covariance: Covariance matrix of the energy which is going to be added to the library
* - weight: Weight matrix of the energy which is going to be added to the library
******************************************************************************/
int readAddSortParams(ReconstructInitSIRENA *reconstruct_init,fitsfile **inLibObject,double samprate,int eventcntLib, double estenergy, gsl_vector *pulsetemplate,gsl_matrix *covariance, gsl_matrix *weight)
{
	int status = EPOK;
	string message = "";
	
	IOData obj;
	
	int runF0orB0val;
	if (strcmp(reconstruct_init->FilterMethod,"F0") == 0)		// Deleting the frequency-zero bin
	{
		runF0orB0val = 0;
	}
	else if (strcmp(reconstruct_init->FilterMethod,"B0") == 0)	// Working without baseline
	{
		runF0orB0val = 1;
	}
	
	// Declare variables
	gsl_vector *energycolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *estenergycolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_matrix *modelsaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *modelsb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix *optimalfiltersaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filters->ofilter_duration);
	gsl_vector *nrmfctrcolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_matrix *weightaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *covarianceaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	
	gsl_matrix *TVaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(TVaux);
	gsl_vector *tEaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(tEaux);
	gsl_matrix *XMaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix_set_zero(XMaux);
	gsl_matrix *YVaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(YVaux);
	gsl_matrix *ZVaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(ZVaux);
	gsl_vector *rEaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(rEaux);
	gsl_matrix *Pabaux = gsl_matrix_alloc(eventcntLib+1, reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Pabaux);
	gsl_matrix *Dabaux = gsl_matrix_alloc(eventcntLib+1, reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Dabaux);
	gsl_matrix *optimalfiltersabaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filters->ofilter_duration);
	gsl_vector *nrmfctrabcolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(nrmfctrabcolumn);
	
	gsl_vector *vectoraux = gsl_vector_alloc(1);
	gsl_vector *vectoraux1 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *vectoraux2 = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	
	// Load values already in the library
	for (int i=0;i<eventcntLib;i++)
	{
		gsl_vector_set(energycolumn,i,gsl_vector_get(reconstruct_init->library_collection->energies,i));
		gsl_vector_set(estenergycolumn,i,gsl_vector_get(reconstruct_init->library_collection->pulse_heights,i));
			
		gsl_matrix_set_row(modelsaux,i,reconstruct_init->library_collection->pulse_templates[i].ptemplate);
		gsl_matrix_set_row(modelsb0aux,i,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate);
		gsl_matrix_set_row(matchedfiltersaux,i,reconstruct_init->library_collection->matched_filters[i].mfilter);
		gsl_matrix_set_row(matchedfiltersb0aux,i,reconstruct_init->library_collection->matched_filters_B0[i].mfilter);
		    
		gsl_matrix_set_row(optimalfiltersaux,i,reconstruct_init->library_collection->optimal_filters[i].ofilter);
		gsl_vector_set(nrmfctrcolumn,i,gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,i));
		      
		gsl_matrix_get_row(vectoraux2,reconstruct_init->library_collection->W,i);
		gsl_matrix_set_row(weightaux,i,vectoraux2);

		gsl_matrix_get_row(vectoraux2,reconstruct_init->library_collection->V,i);
		gsl_matrix_set_row(covarianceaux,i,vectoraux2);
		
		if ((eventcntLib > 1) && (i < eventcntLib-1))
		{
			gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->T,i);
			gsl_matrix_set_row(TVaux,i,vectoraux1);
			gsl_vector_set(tEaux,i,gsl_vector_get(reconstruct_init->library_collection->t,i));  			
			gsl_matrix_get_row(vectoraux2,reconstruct_init->library_collection->X,i);
			gsl_matrix_set_row(XMaux,i,vectoraux2);
			gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->Y,i);
			gsl_matrix_set_row(YVaux,i,vectoraux1);
			gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->Z,i);
			gsl_matrix_set_row(ZVaux,i,vectoraux1);
			gsl_vector_set(rEaux,i,gsl_vector_get(reconstruct_init->library_collection->r,i));  
			gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->PAB,i);
			gsl_matrix_set_row(Pabaux,i,vectoraux1);
			gsl_matrix_get_row(vectoraux1,reconstruct_init->library_collection->DAB,i);
			gsl_matrix_set_row(Dabaux,i,vectoraux1);
			gsl_matrix_set_row(optimalfiltersabaux,i,reconstruct_init->library_collection->optimal_filtersab[i].ofilter);
			gsl_vector_set(nrmfctrabcolumn,i,gsl_vector_get(reconstruct_init->library_collection->nrmfctrsab,i));
		}
	}
	
	// Add new values
	gsl_vector_set(energycolumn,eventcntLib,reconstruct_init->monoenergy);	
	gsl_vector_set(estenergycolumn,eventcntLib,estenergy);
	gsl_matrix_set_row(modelsaux,eventcntLib,pulsetemplate);
	
	gsl_vector_memcpy(vectoraux1,pulsetemplate);
	gsl_vector_scale(vectoraux1,1/reconstruct_init->monoenergy);
	gsl_matrix_set_row(matchedfiltersaux,eventcntLib,vectoraux1);
	
	gsl_vector *optimalfilter = NULL;
	gsl_vector *optimalfilter_f = NULL;
	gsl_vector *optimalfilter_FFT = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	double normalizationFactor;
	if (calculus_optimalFilter (0, 0, reconstruct_init->mode, vectoraux1, vectoraux1->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex, &normalizationFactor))
	{
		message = "Cannot run routine calculus_optimalFilter in writeLibrary";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_matrix_set_row(optimalfiltersaux,eventcntLib,optimalfilter);
	gsl_vector_set(nrmfctrcolumn,eventcntLib,normalizationFactor);
	gsl_vector_free(optimalfilter_f);
	gsl_vector_free(optimalfilter_FFT);
	gsl_vector_complex_free(optimalfilter_FFT_complex);

	gsl_vector *baselinegsl = gsl_vector_alloc(pulsetemplate->size);
	gsl_vector_set_all(baselinegsl,-1.0*reconstruct_init->noise_spectrum->baseline);
	gsl_vector_add(pulsetemplate,baselinegsl);
	gsl_vector_free(baselinegsl);
	gsl_matrix_set_row(modelsb0aux,eventcntLib,pulsetemplate);

	gsl_vector_memcpy(vectoraux1,pulsetemplate);
	gsl_vector_scale(vectoraux1,1/reconstruct_init->monoenergy);
	gsl_matrix_set_row(matchedfiltersb0aux,eventcntLib,vectoraux1);

	matrix2vector(weight,&vectoraux2);
	gsl_matrix_set_row(weightaux,eventcntLib,vectoraux2);
	matrix2vector(covariance,&vectoraux2);
	gsl_matrix_set_row(covarianceaux,eventcntLib,vectoraux2);
	
	// Re-sort
	gsl_vector *energycolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *estenergycolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *modelsrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *modelsaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_vector *modelsrowb0 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *modelsb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_vector *matchedfiltersrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_vector *matchedfiltersrowb0 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_vector *optimalfiltersrow = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filters->ofilter_duration);
	gsl_matrix *optimalfiltersaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filters->ofilter_duration);
	gsl_vector *nrmfctrcolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_matrix *weightaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *covarianceaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_vector *TVrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *TVaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(TVaux1);
	gsl_vector *tEcolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *tEcolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(tEcolumnaux);
	gsl_vector *XMrow = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *XMaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_vector *YVrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *YVaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(YVaux1);
	gsl_vector *ZVrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *ZVaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(ZVaux1);
	gsl_vector *rEcolumn = gsl_vector_alloc(eventcntLib+1);
	gsl_vector *rEcolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(rEcolumnaux);
	gsl_vector *Pabrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *Pabaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Pabaux1);
	gsl_vector *Dabrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *Dabaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
	gsl_matrix_set_zero(Dabaux1);
	gsl_vector *optimalfiltersabrow = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filters->ofilter_duration);
	gsl_matrix *optimalfiltersabaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->library_collection->optimal_filters->ofilter_duration);
	gsl_vector *nrmfctrabcolumnaux = gsl_vector_alloc(eventcntLib+1);
	gsl_vector_set_zero(nrmfctrabcolumnaux);
	
	gsl_permutation *perm = gsl_permutation_alloc(eventcntLib+1);
	// 'gsl_sort_vector_index' indirectly sorts the elements of the vector v into ascending order, storing the resulting
	// permutation in p. The elements of p give the index of the vector element which would have been stored in that position
	// if the vector had been sorted in place. The first element of p gives the index of the least element in v, and the last
	// element of p gives the index of the greatest element in v. The vector v is not changed.
	// Example: tstartaux=(5200 6000 200 3000) tauxsorted=(200 3000 5200 6000) perm=(2 3 0 1)
	gsl_sort_vector_index(perm,energycolumn);
	for (int i=0;i<eventcntLib+1;i++)
	{
		gsl_vector_set(energycolumnaux,i,gsl_vector_get(energycolumn,gsl_permutation_get(perm,i)));
		gsl_vector_set(estenergycolumnaux,i,gsl_vector_get(estenergycolumn,gsl_permutation_get(perm,i)));
		gsl_matrix_get_row(modelsrow,modelsaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(modelsaux1,i,modelsrow);
		gsl_matrix_get_row(modelsrowb0,modelsb0aux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(modelsb0aux1,i,modelsrowb0);
		gsl_matrix_get_row(matchedfiltersrow,matchedfiltersaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(matchedfiltersaux1,i,matchedfiltersrow);
		gsl_matrix_get_row(matchedfiltersrowb0,matchedfiltersb0aux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(matchedfiltersb0aux1,i,matchedfiltersrowb0);
		
		gsl_matrix_get_row(optimalfiltersrow,optimalfiltersaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(optimalfiltersaux1,i,optimalfiltersrow);
		gsl_vector_set(nrmfctrcolumnaux,i,gsl_vector_get(nrmfctrcolumn,gsl_permutation_get(perm,i)));
		
		gsl_matrix_get_row(vectoraux2,weightaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(weightaux1,i,vectoraux2);
		gsl_matrix_get_row(vectoraux2,covarianceaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(covarianceaux1,i,vectoraux2);
		
		gsl_matrix_get_row(TVrow,TVaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(TVaux1,i,TVrow);
				
		gsl_vector_set(tEcolumnaux,i,gsl_vector_get(tEaux,gsl_permutation_get(perm,i)));
		
		gsl_matrix_get_row(XMrow,XMaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(XMaux1,i,XMrow);
		
		gsl_matrix_get_row(YVrow,YVaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(YVaux1,i,YVrow);
		
		gsl_matrix_get_row(ZVrow,ZVaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(ZVaux1,i,ZVrow);
		
		gsl_vector_set(rEcolumnaux,i,gsl_vector_get(rEaux,gsl_permutation_get(perm,i)));
		
		gsl_matrix_get_row(Pabrow,Pabaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(Pabaux1,i,Pabrow);
		
		gsl_matrix_get_row(Dabrow,Dabaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(Dabaux1,i,Dabrow);
		
		gsl_matrix_get_row(optimalfiltersabrow,optimalfiltersabaux,gsl_permutation_get(perm,i));
		gsl_matrix_set_row(optimalfiltersabaux1,i,optimalfiltersabrow);
		gsl_vector_set(nrmfctrabcolumnaux,i,gsl_vector_get(nrmfctrabcolumn,gsl_permutation_get(perm,i)));
	}
	gsl_vector_memcpy(energycolumn,energycolumnaux);
	gsl_vector_memcpy(estenergycolumn,estenergycolumnaux);
	gsl_matrix_memcpy(modelsaux,modelsaux1);

	gsl_matrix_memcpy(modelsb0aux,modelsb0aux1);
	gsl_matrix_memcpy(matchedfiltersaux,matchedfiltersaux1);
	gsl_matrix_memcpy(matchedfiltersb0aux,matchedfiltersb0aux1);
	
	gsl_matrix_memcpy(optimalfiltersaux,optimalfiltersaux1);
	gsl_vector_memcpy(nrmfctrcolumn,nrmfctrcolumnaux);
	
	gsl_matrix_memcpy(weightaux,weightaux1);
	gsl_matrix_memcpy(covarianceaux,covarianceaux1);
	
	gsl_matrix_memcpy(TVaux,TVaux1);
	gsl_vector_memcpy(tEcolumn,tEcolumnaux);
	gsl_matrix_memcpy(XMaux,XMaux1);
	gsl_matrix_memcpy(YVaux,YVaux1);
	gsl_matrix_memcpy(ZVaux,ZVaux1);
	gsl_vector_memcpy(rEcolumn,rEcolumnaux);
	gsl_matrix_memcpy(Pabaux,Pabaux1);
	gsl_matrix_memcpy(Dabaux,Dabaux1);
	gsl_matrix_memcpy(optimalfiltersabaux,optimalfiltersabaux1);
	gsl_vector_memcpy(nrmfctrabcolumn,nrmfctrabcolumnaux);
			
	gsl_permutation_free(perm);
	gsl_vector_free(energycolumnaux);
	gsl_vector_free(estenergycolumnaux);
	gsl_vector_free(modelsrow);
	gsl_matrix_free(modelsaux1);
	gsl_vector_free(modelsrowb0);
	gsl_matrix_free(modelsb0aux1);
	gsl_vector_free(matchedfiltersrow);
	gsl_matrix_free(matchedfiltersaux1);
	gsl_vector_free(matchedfiltersrowb0);
	gsl_matrix_free(matchedfiltersb0aux1);
	gsl_matrix_free(optimalfiltersaux1);
	gsl_vector_free(nrmfctrcolumnaux);
	gsl_matrix_free(weightaux1);
	gsl_matrix_free(covarianceaux1);
	gsl_vector_free(TVrow);
	gsl_matrix_free(TVaux1);
	gsl_vector_free(tEaux);
	gsl_vector_free(tEcolumnaux);
	gsl_vector_free(XMrow);
	gsl_matrix_free(XMaux1);
	gsl_vector_free(YVrow);
	gsl_matrix_free(YVaux1);
	gsl_vector_free(ZVrow);
	gsl_matrix_free(ZVaux1);
	gsl_vector_free(rEaux);
	gsl_vector_free(rEcolumnaux);
	gsl_vector_free(Pabrow);
	gsl_matrix_free(Pabaux1);
	gsl_vector_free(Dabrow);
	gsl_matrix_free(Dabaux1);
	gsl_matrix_free(optimalfiltersabaux1);
	gsl_vector_free(nrmfctrabcolumnaux);
			
	// Add intermeadiate values
	int option;
	int indexa, indexb;
	if (reconstruct_init->monoenergy > gsl_vector_get(energycolumn,eventcntLib-1))	
	{	
		// Intermediate values of the last pair of energies are going to be calculated
		option = 2;
		indexa = eventcntLib-1;
		indexb = eventcntLib;
	}
	else if  (reconstruct_init->monoenergy < gsl_vector_get(energycolumn,0))
	{
		// Intermediate values of the first pair of energies are going to be calculated
		option = 0;
		indexa = 0;
		indexb = 1;
	}
	else
	{
		for (int i=0;i<eventcntLib;i++)
		{
			if ((reconstruct_init->monoenergy > gsl_vector_get(energycolumn,i)) && (reconstruct_init->monoenergy < gsl_vector_get(energycolumn,i+1)))	
			{
				// Intermediate values of the (i,i+1) pair of energies are going to be calculated
				option = 1;
				indexa = i;
				indexb = i+1;
			      
				break;
			}
		}
	}
	
	// Recalculate intermediate values of some new pairs
	if (calculateIntParams(reconstruct_init,indexa,indexb,samprate,runF0orB0val,modelsaux,weightaux,energycolumn, 
		&TVaux,&tEcolumn,&XMaux,&YVaux,&ZVaux,&rEcolumn, &Pabaux, &Dabaux, &optimalfiltersabaux, &nrmfctrabcolumn))
	{
		message = "Cannot run calculateIntParams routine for option 2 in readAddSortParams";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	
	if (option == 1)
	{
		if (calculateIntParams(reconstruct_init,indexa+1,indexb+1,samprate,runF0orB0val,modelsaux,weightaux,energycolumn, 
		&TVaux,&tEcolumn,&XMaux,&YVaux,&ZVaux,&rEcolumn, &Pabaux, &Dabaux, &optimalfiltersabaux, &nrmfctrabcolumn))
		{
			message = "Cannot run calculateIntParams routine for option 2 in readAddSortParams";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
		
	// Write values
	gsl_matrix *matrixaux = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	gsl_matrix *matrixaux2 = gsl_matrix_alloc(1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	
	obj.inObject = *inLibObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"LIBRARY");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.unit = new char [255];
	for (int i=0;i<eventcntLib+1;i++)
	{
	  	obj.iniRow = i+1;
		obj.endRow = i+1;
		strcpy(obj.nameCol,"ENERGY");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"eV");
		gsl_vector_set (vectoraux,0,gsl_vector_get(energycolumn,i));
		if (writeFitsSimple(obj, vectoraux))
		{
		    message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"PHEIGHT");
		strcpy(obj.unit,"ADC");
		gsl_vector_set (vectoraux,0,gsl_vector_get(estenergycolumn,i));
		if (writeFitsSimple(obj, vectoraux))
		{
		    message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"PULSE");
		strcpy(obj.unit,"ADC");
		gsl_matrix_get_row(vectoraux1,modelsaux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
		    message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"PULSEB0");
		strcpy(obj.unit,"ADC");
		gsl_matrix_get_row(vectoraux1,modelsb0aux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
		    message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"MF");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(vectoraux1,matchedfiltersaux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
		    message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"MFB0");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(vectoraux1,matchedfiltersb0aux,i);
		gsl_matrix_set_row(matrixaux,0,vectoraux1);
		if (writeFitsComplex(obj, matrixaux))
		{
		    message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"OF");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(optimalfiltersrow,optimalfiltersaux,i);
		gsl_matrix *optimalfilters_matrix = gsl_matrix_alloc(1,optimalfilter->size);
		gsl_matrix_set_row(optimalfilters_matrix,0,optimalfiltersrow);
		if (writeFitsComplex(obj, optimalfilters_matrix))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"NFCTR");
		strcpy(obj.unit,"ADC");
		gsl_vector_set (vectoraux,0,gsl_vector_get(nrmfctrcolumn,i));
		if (writeFitsSimple(obj, vectoraux))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"COVARM");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(vectoraux2,covarianceaux,i);
		gsl_matrix_set_row(matrixaux2,0,vectoraux2);
		if (writeFitsComplex(obj, matrixaux2))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"WEIGHTM");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(vectoraux2,weightaux,i);
		gsl_matrix_set_row(matrixaux2,0,vectoraux2);
		if (writeFitsComplex(obj, matrixaux2))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		
		if (i < eventcntLib)
		{
			strcpy(obj.nameCol,"TV");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux1,TVaux,i);
			gsl_matrix_set_row(matrixaux,0,vectoraux1);
			if (writeFitsComplex(obj, matrixaux))
			{
			    message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
			    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			strcpy(obj.nameCol,"tE");
			strcpy(obj.unit," ");
			gsl_vector_set (vectoraux,0,gsl_vector_get(tEcolumn,i));
			if (writeFitsSimple(obj, vectoraux))
			{
				message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			strcpy(obj.nameCol,"XM");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux2,XMaux,i);
			gsl_matrix_set_row(matrixaux2,0,vectoraux2);
			if (writeFitsComplex(obj, matrixaux2))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			strcpy(obj.nameCol,"YV");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux1,YVaux,i);
			gsl_matrix_set_row(matrixaux,0,vectoraux1);
			if (writeFitsComplex(obj, matrixaux))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			strcpy(obj.nameCol,"ZV");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux1,ZVaux,i);
			gsl_matrix_set_row(matrixaux,0,vectoraux1);
			if (writeFitsComplex(obj, matrixaux))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			strcpy(obj.nameCol,"rE");
			strcpy(obj.unit," ");
			gsl_vector_set (vectoraux,0,gsl_vector_get(rEcolumn,i));
			if (writeFitsSimple(obj, vectoraux))
			{
				message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			strcpy(obj.nameCol,"PAB");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux1,Pabaux,i);
			gsl_matrix_set_row(matrixaux,0,vectoraux1);
			if (writeFitsComplex(obj, matrixaux))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			
			strcpy(obj.nameCol,"DAB");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(vectoraux1,Dabaux,i);
			gsl_matrix_set_row(matrixaux,0,vectoraux1);
			if (writeFitsComplex(obj, matrixaux))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}

			strcpy(obj.nameCol,"OFAB");
			strcpy(obj.unit," ");
			gsl_matrix_get_row(optimalfiltersabrow,optimalfiltersabaux,i);
			gsl_matrix_set_row(optimalfilters_matrix,0,optimalfiltersabrow);
			if (writeFitsComplex(obj, optimalfilters_matrix))
			{
				message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
			gsl_matrix_free(optimalfilters_matrix);
					
			strcpy(obj.nameCol,"NFCTRAB");
			strcpy(obj.unit," ");
			gsl_vector_set (vectoraux,0,gsl_vector_get(nrmfctrabcolumn,i));
			if (writeFitsSimple(obj, vectoraux))
			{
				message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
	}
	
	// Free allocate of GSL vectors
	gsl_vector_free(optimalfilter);
		
	gsl_vector_free(vectoraux);
	gsl_vector_free(vectoraux1);
	gsl_vector_free(vectoraux2);
	gsl_matrix_free(matrixaux);
	gsl_matrix_free(matrixaux2);
		
	gsl_vector_free(energycolumn);
	gsl_vector_free(estenergycolumn);
	gsl_matrix_free(modelsaux);
	gsl_matrix_free(modelsb0aux);
	gsl_matrix_free(matchedfiltersaux);
	gsl_matrix_free(matchedfiltersb0aux);
	gsl_matrix_free(optimalfiltersaux);
	gsl_vector_free(optimalfiltersrow);
	gsl_vector_free(nrmfctrcolumn);
	gsl_matrix_free(weightaux);
	gsl_matrix_free(covarianceaux);
	gsl_matrix_free(TVaux);
	gsl_vector_free(tEcolumn);
	gsl_matrix_free(XMaux);
	gsl_matrix_free(YVaux);
	gsl_matrix_free(ZVaux);
	gsl_vector_free(rEcolumn);
	gsl_matrix_free(Pabaux);
	gsl_matrix_free(Dabaux);
	gsl_matrix_free(optimalfiltersabaux);
	gsl_vector_free(optimalfiltersabrow);
	gsl_vector_free(nrmfctrabcolumn);
		
	return (EPOK);
}
/*xxxx end of SECTION A16 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A17 ************************************************************
* calculateIntParams function: This function calculates some intermediate scalars, vectors and matrices (TV, tE, XM, YV, ZV, rE, PAB, DAB, OFAB and NRMFCTRAB).
*
* - Declare variables and allocate GSL vectors and matrices
* - Calculate intermediate scalars, vectors and matrices 
* - Free allocate of GSL vectors and matrices
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'pulse_length'
* - indexa: Low index of the library to calculate the intermediate params
* - indexb: High index of the library to calculate the intermediate params
* - samprate: Sampling rate
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - modelsaux, weightaux, energycolumn: Input data
* - TVaux, tEcolumn, XMaux, YVaux, ZVaux, rEcolumn, Pabaux, Dabaux, optimalfiltersabaux, nrmfctrabcolumn: Input/output intermediate parameters
******************************************************************************/
int calculateIntParams(ReconstructInitSIRENA *reconstruct_init, int indexa, int indexb, double samprate, int runF0orB0val, 
	gsl_matrix *modelsaux,gsl_matrix *weightaux,gsl_vector *energycolumn, 
	gsl_matrix **TVaux,gsl_vector **tEcolumn,gsl_matrix **XMaux,gsl_matrix **YVaux,gsl_matrix **ZVaux,gsl_vector **rEcolumn, gsl_matrix **Pabaux, gsl_matrix **Dabaux, gsl_matrix **optimalfiltersabaux, gsl_vector **nrmfctrabcolumn)
{
	int status = EPOK;
        string message = "";
	
	// Declare variables and allocate GSL vectors and matrices
	gsl_vector *T = gsl_vector_alloc(reconstruct_init->pulse_length);
	double t;
	gsl_vector *Y = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *Z = gsl_vector_alloc(reconstruct_init->pulse_length);
	double r;
	gsl_vector *Pab = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *Dab = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *Walpha = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
	gsl_vector *Xvector = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *X = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
	double Ea,Eb;
	
	gsl_vector *vectoraux = gsl_vector_alloc(1);
	gsl_vector *vectoraux1 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *vectoraux2 = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	
	gsl_vector *optimalfilter = NULL;
	gsl_vector *optimalfilter_f = NULL;
	gsl_vector *optimalfilter_FFT = NULL;
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	double normalizationFactor;
		
	// Calculate
	gsl_matrix_get_row(vectoraux1,modelsaux,indexb);
	gsl_vector_memcpy(T,vectoraux1);
	gsl_matrix_get_row(vectoraux1,modelsaux,indexa);
	gsl_vector_sub(T,vectoraux1);
	gsl_matrix_set_row(*TVaux,indexa,T);
	
	gsl_matrix_get_row(vectoraux2,weightaux,indexa);
	vector2matrix(vectoraux2,&Walpha);
	gsl_blas_dgemv(CblasNoTrans, 1.0, Walpha, T, 0.0, vectoraux1);
	gsl_blas_ddot(T,vectoraux1,&t);
	gsl_vector_set(*tEcolumn,indexa,t);
	
	gsl_matrix_get_row(vectoraux2,weightaux,indexb);
	gsl_vector_memcpy(Xvector,vectoraux2);
	gsl_matrix_get_row(vectoraux2,weightaux,indexa);
	gsl_vector_sub(Xvector,vectoraux2);
	gsl_vector_scale(Xvector,1/t);
	gsl_matrix_set_row(*XMaux,indexa,Xvector);
	
	gsl_blas_dgemv(CblasNoTrans, 1.0, Walpha, T, 0.0, Y);
	gsl_vector_scale(Y,1/t);
	gsl_matrix_set_row(*YVaux,indexa,Y);
	
	vector2matrix(Xvector,&X);
	gsl_blas_dgemv(CblasNoTrans, 1.0, X, T, 0.0, Z);
	gsl_matrix_set_row(*ZVaux,indexa,Z);
	
	gsl_blas_ddot(Z,T,&r);
	r=1/r;
	gsl_vector_set(*rEcolumn,indexa,r);
	
	Ea = gsl_vector_get(energycolumn,indexa);
	Eb = gsl_vector_get(energycolumn,indexb);
	gsl_vector_memcpy(Pab,T);
	gsl_vector_scale(Pab,-Ea/(Eb-Ea));
	gsl_matrix_get_row(vectoraux1,modelsaux,indexa);
	gsl_vector_add(Pab,vectoraux1);
	gsl_matrix_set_row(*Pabaux,indexa,Pab);
	
	gsl_vector_memcpy(Dab,T);
	gsl_vector_scale(Dab,1/(Eb-Ea));
	gsl_matrix_set_row(*Dabaux,indexa,Dab);
	
	if (calculus_optimalFilter (0, 0, 0, Dab, Dab->size, samprate, runF0orB0val, reconstruct_init->noise_spectrum->noisefreqs, reconstruct_init->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex, &normalizationFactor))
	{
		message = "Cannot run routine calculus_optimalFilter in writeLibrary";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_matrix_set_row(*optimalfiltersabaux,indexa,optimalfilter);
	gsl_vector_set(*nrmfctrabcolumn,indexa,normalizationFactor);
	
	// Free allocate of GSL
	gsl_vector_free(T);
	gsl_vector_free(Y);
	gsl_vector_free(Z);
	gsl_vector_free(Pab);
	gsl_vector_free(Dab);
	gsl_matrix_free(Walpha);
	gsl_vector_free(Xvector);
	gsl_matrix_free(X);
	
	gsl_vector_free(vectoraux);
	gsl_vector_free(vectoraux1);
	gsl_vector_free(vectoraux2);
	
	gsl_vector_free(optimalfilter);
	gsl_vector_free(optimalfilter_f);
	gsl_vector_free(optimalfilter_FFT);
	gsl_vector_complex_free(optimalfilter_FFT_complex);
		
	return (EPOK);
}
/*xxxx end of SECTION A17 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A18 ************************************************************
* matrix2vector function: This function converts an input n x n square matrix into an output n^2 vector.
*                         It puts the first row of the matrix (n elements) in the first n elements of the vector (from 0 to n-1),
*                         the second row of the matrix in the (n,2n-1)elements of the vector and so on.
*
* Parameters:
* - matrixin: Input square matrix whose dimensions are nxn
* - vectorout: Output vector whose length is n^2
******************************************************************************/
int matrix2vector (gsl_matrix *matrixin, gsl_vector **vectorout)
{
	// matrixin is a square matrix
	int dim = matrixin->size1;
	gsl_vector *vectoraux = gsl_vector_alloc(dim);

	for (int i=0;i<dim;i++)
	{
		gsl_matrix_get_row(vectoraux,matrixin,i);
		for (int j=0;j<dim;j++)
		{
			gsl_vector_set(*vectorout,i*dim+j,gsl_vector_get(vectoraux,j));
		}
	}

	gsl_vector_free(vectoraux);

	return (EPOK);
}
/*xxxx end of SECTION A18 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A19 ************************************************************
* vector2matrix function: This function converts an input n^2 vector into an output n x n square matrix.
*                         It puts the first n elements of the vector in the first row of the matrix,
*                         the second group of n elements (from n to 2n-1) of the vector in the second row and so on.
*
* Parameters:
* - vectorin: Input vector whose length is n^2
* - matrixout: Output matrix whose dimensions are nxn
******************************************************************************/
int vector2matrix (gsl_vector *vectorin, gsl_matrix **matrixout)
{
	// matrixin is a square matrix
	double dim = sqrt(vectorin->size);

	for (int i=0;i<dim;i++)
	{
		for (int j=0;j<dim;j++)
		{
			gsl_matrix_set(*matrixout,i,j,gsl_vector_get(vectorin,i*dim+j));
		}
	}

	return (EPOK);
}
/*xxxx end of SECTION A19 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A20 ************************************************************
* convertI2R: This funcion converts the current space into a resistance space. 
*             The 'invector' filled with currents is filled in with resistances at the output.
* 
* - Read R0, I0_START=Ibias, ADUCNV and IMIN from the 'RecordFile' by using the pointer '(*reconstruct_init)->record_file_fptr'
* - Conversion according to 'EnergyMethod'=I2R
* - Conversion according to 'EnergyMethod'=I2RBIS
* 	- Read also RL, LFILTER and TTR from the 'RecordFile'
* - From this point forward, I2R and I2RBIS are completely equivalent to OPTFILT
* - Because in 'runEnergy' the record (TesRecord) is used => The I2R or I2RBIS transformed record has to be used
* 
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
* - record: Structure that contains the input record
* - invector: Input current vector & output resistance vector
******************************************************************************/
int convertI2R (ReconstructInitSIRENA** reconstruct_init, TesRecord **record, gsl_vector **invector)
{
	int status = EPOK;
	string message="";
  
	double R0;		// Operating point resistance [Ohm]
	double Ibias;		// I0_STARRT = Initial bias current [A]
	double aducnv;		// ADU conversion factor [A/ADU]
	//double baselineNew;
	double baseline;
	double Imin;		// Current corresponding to 0 ADU [A]

	// Read R0, I0_START=Ibias, ADUCNV and IMIN from the 'RecordFile' by using the pointer '(*reconstruct_init)->record_file_fptr'
	// The 'RecordFile' can not be openned because it is already open from tesrecontruction.c
	char recordFileName[256];
	strncpy(recordFileName, (*reconstruct_init)->record_file,255);
	recordFileName[255]='\0';
	char extname[20];
	char keyname[10];
	strcpy(extname,"RECORDS");
	//strcpy(extname,"TESDATASTREAM");
	if (fits_movnam_hdu((*reconstruct_init)->record_file_fptr, ANY_HDU,extname, 0, &status))
	{
		message = "Cannot move to HDU " + string(extname);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	strcpy(keyname,"R0");
	fits_read_key((*reconstruct_init)->record_file_fptr,TDOUBLE,keyname, &R0,NULL,&status);
	strcpy(keyname,"I0_START");
	fits_read_key((*reconstruct_init)->record_file_fptr,TDOUBLE,keyname, &Ibias,NULL,&status);
	strcpy(keyname,"ADUCNV");
	fits_read_key((*reconstruct_init)->record_file_fptr,TDOUBLE,keyname, &aducnv,NULL,&status);
	strcpy(keyname,"IMIN");
	fits_read_key((*reconstruct_init)->record_file_fptr,TDOUBLE,keyname, &Imin,NULL,&status);
	if (status != 0)
	{
		message = "Cannot read keyword in convertI2R";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	if (strcmp((*reconstruct_init)->EnergyMethod,"I2R") == 0)
	{
		double DeltaI = abs(((*reconstruct_init)->noise_spectrum->baseline*aducnv+Imin) - Ibias);
		//baselineNew = R0 -R0*(DeltaI/Ibias)/(1.+(DeltaI/Ibias));

		// DeltaI <- I-Ibias
		// R <- R0 - R0*(abs(DeltaI)/Ibias)/(1+abs(DeltaI)/Ibias)

		gsl_vector *invector_modified = gsl_vector_alloc((*invector)->size);
		gsl_vector_scale(*invector,aducnv);             // invector = I(ADC)*ADUCNV
		gsl_vector_add_constant(*invector,Imin);	// invector = I(Amp) = I(ADC)*ADUCNV+Imin
		gsl_vector_add_constant(*invector,-1*Ibias); 	// invector = DeltaI = abs(I(Amp)-Ibias)
		for (int i=0;i<(*invector)->size;i++)
		{
			if (gsl_vector_get(*invector,i)<0) gsl_vector_set(*invector,i,(-1.*gsl_vector_get(*invector,i)));
		}
		gsl_vector_scale(*invector,1./Ibias); 			// invector = DeltaI/Ibias = (I(Amp)-Ibias)/Ibias
		gsl_vector_memcpy(invector_modified,*invector);  	// invector_modified = invector = DeltaI/Ibias
		gsl_vector_add_constant(invector_modified,+1.0);	// invector_modified = 1 + DeltaI/Ibias
		gsl_vector_div(*invector,invector_modified);     	// invector = invector/invector_modified = (DeltaI/Ibias)/(1+DeltaI/Ibias)
		gsl_vector_scale(*invector,-1.*R0);			// invector = -R0*(DeltaI/Ibias)/(1+DeltaI/Ibias)
		gsl_vector_add_constant(*invector,R0); 			// invector = R0 - R0*(DeltaI/Ibias)/(1+DeltaI/Ibias)
		gsl_vector_free(invector_modified);
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RBIS") == 0)
	{
		double RL;		// Shunt/load resistor value [oHM]
		double TTR;		// Transformer Turns Ratio
		double LFILTER;		// Circuit filter inductance [H]
		double V0;
		double L;
		gsl_vector *I = gsl_vector_alloc((*invector)->size);
		gsl_vector *dI = gsl_vector_alloc((*invector)->size);
		gsl_vector *invector_modified = gsl_vector_alloc((*invector)->size);

		strcpy(keyname,"RL");
		fits_read_key((*reconstruct_init)->record_file_fptr,TDOUBLE,keyname, &RL,NULL,&status);
		strcpy(keyname,"TTR");
		fits_read_key((*reconstruct_init)->record_file_fptr,TDOUBLE,keyname, &TTR,NULL,&status);
		strcpy(keyname,"LFILTER");
		fits_read_key((*reconstruct_init)->record_file_fptr,TDOUBLE,keyname, &LFILTER,NULL,&status);
		if (status != 0)
		{
			message = "Cannot read keyword in convertI2R and I2RBIS";
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		V0 = Ibias*(R0+RL);			// V0 = I0(R0+RL)
		L = 2*LFILTER/(pow(TTR,2));		// L = 2LFILTER/(TTR)\B2

		// I
		gsl_vector_memcpy(invector_modified,*invector);
		gsl_vector_scale(*invector,aducnv);             // invector = I(ADC)*ADUCNV
		gsl_vector_add_constant(*invector,Imin);	// invector = I(ADC)*ADUCNV+Imin
		gsl_vector_scale(*invector,-1.0); 		// invector = Ibias-(I(ADC)*ADUCNV+Imin)
		gsl_vector_add_constant(*invector,Ibias);
		gsl_vector_memcpy(I,*invector);

		// dI/dt
		gsl_vector_memcpy(dI,I);
		if (derivative (&dI, dI->size))
		{
			message = "Cannot run routine derivative for convertI2R and I2RBIS";
			EP_EXIT_ERROR(message,EPFAIL);
		}

		// R = (V0-LdI/dt)/I
		gsl_vector_memcpy(*invector,dI);
		gsl_vector_scale(*invector,-L);
		gsl_vector_add_constant(*invector,V0);
		gsl_vector_div(*invector,I);
				
		// Change to positive polarity
		//gsl_vector_scale(*invector,-1.0);
		//gsl_vector_add_constant(*invector,2*(*reconstruct_init)->noise_spectrum->baseline);
		
		gsl_vector_free(I);
		gsl_vector_free(dI);
		gsl_vector_free(invector_modified);
	}
	
	strcpy((*reconstruct_init)->EnergyMethod,"OPTFILT"); // From this point forward, I2R and I2RBIS are completely equivalent to OPTFILT

	for (int i=0;i<(*invector)->size;i++)	// Because in 'runEnergy' the record (TesRecord) is used => The I2R or I2RBIS transformed record has to be used
	{
		(*record)->adc_double[i] = gsl_vector_get(*invector,i);
	}
  
	return(EPOK);
}
/*xxxx end of SECTION A20 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B ************************************************************
* runEnergy: This function is responsible for calculating the pulse energy applying different methods ('EnergyMethod'):
* 				- Optimal filtering ('OPTFILT')
* 				- Covariance matrix ('WEIGHT')
* 				- Covariance matrix simplified ('WEIGHTN')
* 				- Conversion from current to resistance space ('I2R' or 'I2RBIS')
* 	     It only runs in PRODUCTION mode
*
* - Declare variables
* - Store the record in 'invector' ('loadRecord')
* - Check Quality: If there are no valid pulses in the record (not unsaturated pulses) => The task finishes
* - For each not saturated pulse:
* 	- Establish the pulse grade (HighRes=0, MidRes=1, LowRes=-1) and the optimal filter length
* 	- Pulse: Load the proper piece of the record in 'pulse'
* 	- If 'OPTFILT' (or 'I2R') or 'I2RBIS':
* 	  If the 'OPIter'=1, in the first iteration ('numiteration'=0) the values of 'maxDER' and 'maxDERs' are used in 'find_matchedfilter' or 'dind_optimalfilter'
* 	  getting the values of the 'energies' which straddle the 'maxDER' ('Ealpha' and 'Ebeta'). In there are more iterations (if it is necessary because the
* 	  calculated 'energy' is out of ['Ealpha','Ebeta'], the values of 'energy' and 'energies' are used in 'find_matchedfilter' or 'find_optimalfilter').
* 	  If 'energy' is in ['Elapha','Ebeta'] the iterative process stops.
* 	  	- If 'OFLib'=0:
*	 		- Filter (find the matched filter and load it in 'filter') ('find_matchedfilter')
* 		    - Calculate the optimal filter
* 		- If 'OFLib'=1:
*	 		- Filter (find the optimal filter and load it in 'optimalfilter') ('find_optimalfilter')
*	 		- If 'FilterDomain'=F => FFT('optimalfilter')
*		- Calculate the energy of each pulse
*	- else if 'WEIGHT' or 'WEIGHTN':
*		- Get the indexes of the two energies which straddle the pulse
*		- Calculate the energy of each pulse
*	- Subtract the pulse model from the record
*	- Write info of the pulse in the output intemediate file if 'intermediate'=1
* - Not valid pulse => Its info has to be also stored in the intermediate file (if 'intermediate'=1) and in the structure 'pulsesInRecord'
* - Free allocate of GSL vectors
*
* Parameters:
* - record:
* - reconstruct_init:
* - pulsesInRecord:
* - optimalFilter:
******************************************************************************/
void runEnergy(TesRecord* record,ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter)
{
	const char * create= "runEnergy v.15.0.0";	// In order to set "CREATOR" keyword

	// Declare variables
	string message="";
	int status = EPOK;

	fitsfile *dtcObject = NULL;	    // Object which contains information of the intermediate output FITS file
	if ((*reconstruct_init)->intermediate == 1)
	{
		char dtcName[256];
		strncpy(dtcName,(*reconstruct_init)->detectFile,255);
		dtcName[255]='\0';
	}

	int TorF;
	if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)			// Time
	{
		TorF=0;
	}
	else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)	// Frequency
	{
		TorF=1;
	}

	int runF0orB0val;
	if (strcmp((*reconstruct_init)->FilterMethod,"F0") == 0)		// Deleting the frequency-zero bin
	{
		runF0orB0val = 0;
	}
	else if (strcmp((*reconstruct_init)->FilterMethod,"B0") == 0)	// Working without baseline
	{
		runF0orB0val = 1;
	}

	// I2R or I2RBIS methods convert I into R at the beginnig and after that 'I2R' or 'I2RBIS' are 'OPTFILT'
	int runEMethod;
	if (strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0)
	{
		runEMethod = 0;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHT") == 0)
	{
		runEMethod = 1;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0)
	{
		runEMethod = 2;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"I2RBIS") == 0)
	{
		runEMethod = 3;
	}

	int OFlength_strategy;
	if (strcmp((*reconstruct_init)->OFStrategy,"FREE") == 0)
	{
		OFlength_strategy = 0;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"BASE2") == 0)
	{
		OFlength_strategy = 1;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"BYGRADE") == 0)
	{
		OFlength_strategy = 2;
	}
	else if (strcmp((*reconstruct_init)->OFStrategy,"FIXED") == 0)
	{
		OFlength_strategy = 3;
	}

	long energyInLibrary_row;

	double normalizationFactor;
	double energy;

	double Ealpha, Ebeta;
	gsl_vector *optimalfilter = NULL;	// Resized optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector *optimalfilter_f = NULL;		// Resized optimal filter f's when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector *optimalfilter_FFT = NULL;	// Resized optimal filter magnitudes when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;
	
	gsl_vector *Pab = NULL;

	gsl_vector *model;

	gsl_vector *recordAux;

	gsl_vector *pulse = NULL;

	gsl_vector *filtergsl = NULL;			// Matched filter values (time domain)

	int iter;
	gsl_vector_view(temp);

	double tstartSamplesRecord;
	double tstartRecord;
	double tstartRecordSamples = floor(record->time/record->delta_t+0.5);	// Close integer

	int indexEalpha = 0;
	int indexEbeta = 0;

	long resize_mf;
	int pulseGrade; // HighRes=0, MidRes=1, LowRes=-1

	// Store the record in 'invector'
	recordAux = gsl_vector_alloc(record->trigger_size);
	if (loadRecord(record, &tstartRecord, &recordAux))
	{
		message = "Cannot run routine loadRecord";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	//int eventsz = recordAux->size;	// Just in case the last record has been filled in with 0's => Re-allocate invector

	// Subtract the baseline if 'runF0orB0val'= 1 ('FilterMethod'=B0)
	if (runF0orB0val == 1)
	{
		gsl_vector *baselinegsl = gsl_vector_alloc(recordAux->size);
		//gsl_vector_set_all(baselinegsl,-1.0*(*reconstruct_init)->baseline);
		gsl_vector_set_all(baselinegsl,-1.0*(*reconstruct_init)->noise_spectrum->baseline);
		gsl_vector_add(recordAux,baselinegsl);
		gsl_vector_free(baselinegsl);
	}

	// Check Quality: If there are no valid pulses in the record => The task finishes
	iter = 0;
	for (int i=0; i<(*pulsesInRecord)->ndetpulses;i++)
	{
		if ((*pulsesInRecord)->pulses_detected[i].quality >= 10)  //saturated
		{
			iter++;
		}
	}
	if (iter == (*pulsesInRecord)->ndetpulses)
	{
		message = "There are no unsaturated pulses in one record";
		EP_PRINT_ERROR(message,-999);
	}

	model =gsl_vector_alloc((*reconstruct_init)->pulse_length);
	gsl_vector *filtergsl_aux;

	for (int i=0; i<(*pulsesInRecord)->ndetpulses ;i++)
	{
		// Pulses are going to be validated by checking its quality
		//if ((*pulsesInRecord)->pulses_detected[i].quality < 10)	// Not saturated
		//{
			//resize_mf = (*pulsesInRecord)->pulses_detected[i].pulse_duration; // Resize the matched filter by using the tstarts
			//resize_mf = pow(2,floor(log2(resize_mf)));

			// Establish the pulse grade (HighRes=0, MidRes=1, LowRes=-1) and the optimal filter length
			if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
			else								(*pulsesInRecord)->pulses_detected[i].grade1 = (*pulsesInRecord)->pulses_detected[i].pulse_duration;

			if (pulseGrading(*reconstruct_init,(*pulsesInRecord)->pulses_detected[i].grade1,(*pulsesInRecord)->pulses_detected[i].grade2_1,OFlength_strategy,&pulseGrade,&resize_mf))
			{
				message = "Cannot run routine pulseGrading";
				EP_EXIT_ERROR(message,EPFAIL);
			}
			//cout<<"resize_mf="<<resize_mf<<endl;

			// Pulse: Load the proper piece of the record in 'pulse'
			tstartSamplesRecord = floor((*pulsesInRecord)->pulses_detected[i].Tstart/record->delta_t+0.5)-tstartRecordSamples;
			pulse = gsl_vector_alloc(resize_mf);
			temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf);
			gsl_vector_memcpy(pulse,&temp.vector);

			// EnergyMethod != WEIGHT/WEIGHTN => Filter and optimal filter
			if (strcmp((*reconstruct_init)->EnergyMethod,"OPTFILT") == 0)
			{
				bool iterate = true;
				gsl_matrix *Estraddle = gsl_matrix_alloc(2,(*reconstruct_init)->library_collection->ntemplates);
				int numiteration = -1;

				do
				{
					numiteration++;

					// Filter (find the matched filter and load it in 'filter')
					if ((*reconstruct_init)->OFLib == 0)
					{
						filtergsl_aux= gsl_vector_alloc((*reconstruct_init)->library_collection->matched_filters[0].mfilter_duration);
						Pab = gsl_vector_alloc((*reconstruct_init)->library_collection->matched_filters[0].mfilter_duration);
						if (numiteration == 0)
						{
							if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
							{
								if (find_matchedfilter(runF0orB0val, (*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_aux, &Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_matchedfilter for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
							else //(*reconstruct_init)->OFInterp = DAB
							{
								if (find_matchedfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_aux, &Pab, &Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_matchedfilterDAB for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
						}
						else
						{
							if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
							{
								if (find_matchedfilter(runF0orB0val, energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl_aux, &Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_matchedfilter for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
							else //(*reconstruct_init)->OFInterp = DAB
							{
								if (find_matchedfilterDAB(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl_aux, &Pab, &Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_matchedfilterDAB for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
						}

					}
					else if ((*reconstruct_init)->OFLib == 1)
					{
						filtergsl_aux= gsl_vector_alloc((*reconstruct_init)->library_collection->optimal_filters[0].ofilter_duration);
						Pab = gsl_vector_alloc((*reconstruct_init)->library_collection->optimal_filters[0].ofilter_duration);
						if (numiteration == 0)
						{
							if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
							{
								if (find_optimalfilter((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_aux, &normalizationFactor, &Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_optimalfilter for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
							else //(*reconstruct_init)->OFInterp = DAB
							{
								if (find_optimalfilterDAB((*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_aux, &normalizationFactor, &Pab,&Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_optimalfilter for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
						}
						else
						{	
						        if (strcmp((*reconstruct_init)->OFInterp,"MF") == 0)
							{
								if (find_optimalfilter(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl_aux, &normalizationFactor, &Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_optimalfilter for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
							else //(*reconstruct_init)->OFInterp = DAB
							{
								if (find_optimalfilterDAB(energy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl_aux, &normalizationFactor, &Pab, &Ealpha, &Ebeta))
								{
									message = "Cannot run routine find_optimalfilter for filter interpolation";
									EP_EXIT_ERROR(message,EPFAIL);
								}
							}
						}
					}

					gsl_matrix_set(Estraddle,0,numiteration,Ealpha);
					gsl_matrix_set(Estraddle,1,numiteration,Ebeta);

					filtergsl = gsl_vector_alloc(resize_mf);			// Filter values
					temp = gsl_vector_subvector(filtergsl_aux,0,resize_mf);
					gsl_vector_memcpy(filtergsl,&temp.vector);
					gsl_vector_free(filtergsl_aux);

					if ((*reconstruct_init)->OFLib == 0)
					{
						// Calculate the optimal filter
						if (calculus_optimalFilter (TorF, (*reconstruct_init)->intermediate, (*reconstruct_init)->mode, filtergsl, filtergsl->size, 1.0/record->delta_t, runF0orB0val, (*reconstruct_init)->noise_spectrum->noisefreqs, (*reconstruct_init)->noise_spectrum->noisespec, &optimalfilter, &optimalfilter_f, &optimalfilter_FFT, &optimalfilter_FFT_complex, &normalizationFactor))
						{
							message = "Cannot run routine calculus_optimalFilter";
							EP_EXIT_ERROR(message,EPFAIL);
						}
					}
					else if ((*reconstruct_init)->OFLib == 1)
					{
						optimalfilter = gsl_vector_alloc(filtergsl->size);
						gsl_vector_memcpy(optimalfilter,filtergsl);
						if (TorF == 1)
						{
							optimalfilter_FFT_complex = gsl_vector_complex_alloc(filtergsl->size);
							if (FFT(optimalfilter,optimalfilter_FFT_complex,filtergsl->size/(1/record->delta_t)))
							{
								message = "Cannot run FFT routine to calculate optimal filter trnsformation";
								EP_EXIT_ERROR(message,EPFAIL);
							}
							gsl_vector_complex_set(optimalfilter_FFT_complex,0,gsl_complex_rect(0.0,0.0));
							for (int k=0;k<filtergsl->size;k++)
							{
								gsl_vector_complex_set(optimalfilter_FFT_complex,k,gsl_complex_conjugate(gsl_vector_complex_get(optimalfilter_FFT_complex,k)));
							}
						}
					}

					// Calculate the energy of each pulse
					if (calculateEnergy(pulse,pulseGrade,optimalfilter,optimalfilter_FFT_complex,runEMethod,indexEalpha,indexEbeta,(*reconstruct_init),TorF,normalizationFactor,1/record->delta_t,Pab,&energy))
					{
						message = "Cannot run calculateEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_EXIT_ERROR(message,EPFAIL);
					}

					if (((Ealpha == 0) && (energy >= Ebeta)) || ((Ebeta == 0) && (energy <= Ealpha)) || ((Ealpha != Ebeta) && ((energy < Ealpha) || (energy > Ebeta))))
					{
						if (numiteration != 0)
						{
							for (int j=0; j<numiteration; j++)
							{
								if ((Ealpha == gsl_matrix_get(Estraddle,0,j)) && (Ebeta == gsl_matrix_get(Estraddle,1,j)))
								{
									iterate = false;
									break;
								}
							}
						}
						if ((*reconstruct_init)-> OFIter == 1)	iterate = true;
						else iterate = false;
					}
					else iterate = false;

				} while (iterate);

				gsl_matrix_free(Estraddle);
			}
			else
			{
				// Get the indexes of the two energies which straddle the pulse
				if (find_Esboundary((*pulsesInRecord)->pulses_detected[i].maxDER,(*reconstruct_init)->library_collection->maxDERs,&indexEalpha,&indexEbeta))
				{
					message = "Cannot run routine find_Esboundary for filter interpolation";
					EP_EXIT_ERROR(message,EPFAIL);
				}

				if ((indexEalpha == indexEbeta) && ((*reconstruct_init)->library_collection->ntemplates-1 == indexEalpha))
				{
					message = "Not eneough info in the library. At least is necessary to add a new last row with energy higher than the energy of the pulses in the input FITS file";
					EP_EXIT_ERROR(message,EPFAIL);
				}

				// Calculate the energy of each pulse
				if (calculateEnergy(pulse,pulseGrade,optimalfilter,optimalfilter_FFT_complex,runEMethod,indexEalpha,indexEbeta,(*reconstruct_init),TorF,normalizationFactor,1/record->delta_t,Pab,&energy))
				{
					message = "Cannot run calculateEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
					EP_EXIT_ERROR(message,EPFAIL);
				}
			}

			// Subtract the pulse model from the record
			if (find_model_energies(energy, (*reconstruct_init), &model))
			{
				message = "Cannot run find_model_energies routine for pulse i=" + boost::lexical_cast<std::string>(i);
				EP_EXIT_ERROR(message,EPFAIL);
			}

			for (int j=tstartSamplesRecord;j<(tstartSamplesRecord+(model->size));j++)
			{
				gsl_vector_set(recordAux,j,gsl_vector_get(recordAux,j)-gsl_vector_get(model,j-tstartSamplesRecord));
			}

			// Write info of the pulse in the output intemediate file if 'intermediate'=1
			if ((*reconstruct_init)->intermediate == 1)
			{
				if (writeFilterHDU(reconstruct_init, i, normalizationFactor, energy, optimalfilter, optimalfilter_f, optimalfilter_FFT, &dtcObject, create))
				{
					message = "Cannot run writeFilterHDU routine for pulse i=" + boost::lexical_cast<std::string>(i);
					EP_EXIT_ERROR(message,EPFAIL);
				}
			}

			(*pulsesInRecord)->pulses_detected[i].energy = energy/1e3;	// In SIXTE, SIGNAL is in keV
		/*}
		else	// Not valid pulse => Its info has to be also stored in the intermediate file (if 'intermediate'=1) and in the structure 'pulsesInRecord'
		{
			optimalfilter = gsl_vector_alloc((*reconstruct_init)->pulse_length);
			optimalfilter_f = gsl_vector_alloc((*reconstruct_init)->pulse_length);
			optimalfilter_FFT = gsl_vector_alloc((*reconstruct_init)->pulse_length);
			gsl_vector_set_zero(optimalfilter);
			gsl_vector_set_zero(optimalfilter_f);
			gsl_vector_set_zero(optimalfilter_FFT);

			if ((*reconstruct_init)->intermediate == 1)
			{
				if (writeFilterHDU(reconstruct_init, i, 0.0, 0.0, optimalfilter, optimalfilter_f, optimalfilter_FFT, &dtcObject, create))
				{
					message = "Cannot run writeFilterHDU routine for pulse i=" + boost::lexical_cast<std::string>(i);
					EP_EXIT_ERROR(message,EPFAIL);
				}
			}

			if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
			else								(*pulsesInRecord)->pulses_detected[i].grade1 = 0;
			(*pulsesInRecord)->pulses_detected[i].energy = 0.0;
		} // End if*/

		// Free allocate of GSL vectors
		gsl_vector_free(optimalfilter);
		gsl_vector_free(optimalfilter_f);
		gsl_vector_free(optimalfilter_FFT);
		if ((*pulsesInRecord)->pulses_detected[i].quality < 10)
		{
			gsl_vector_free(pulse);
			gsl_vector_free(filtergsl);
			gsl_vector_complex_free(optimalfilter_FFT_complex);
			gsl_vector_free(Pab);
		}
	} // End for

	return;
}
/*xxxx end of SECTION B xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B1 ************************************************************
* calculus_optimalFilter: This function calculates the optimal filter for a pulse whose matched filter is provided as input parameter.
*                         An optimal filter is just a matched filter that has been adjusted based on the noise spectrum of the system.
*
* It is assumed that the are modelled as a form H·s(t), where s(t) is constant and H is simply an amplitude scalar for the photon energy.
* That is, it is assumed all pulses are scaled versions of a template, s(t).
* This assumption is incorrect, and the variation of pulse shape with energy is a large source of nonlinearity in the calculated energy.
*
* It is also assumed that the noise is stationary.
*
*                                                                                                          (P(f)-E·S(f))²
* For a given pulse p(t), a least-squares minimization is performed to find the best value of H: X² = SUM ---------------- where N²(f) is the power
* spectrum of the noise.                                                                                         N²(f)
*
*                                                       P(f)·S(f)
* Minimizing yields the optimal value for H: H = k·SUM ----------- where k is a normalization factor to make the estimate in units of energy
*                                                         N²(f)
*
* Transforming back to the time domain: H = k·SUM p(t)*op(t) where op(t) is the inverse Fourier Transform of S(f)/N²(f)
* If the power spectrum of the noise were white (constant), the optimal filter would be identical to the pulse template (average pulse shape).
*
*                    S(f)                1                               |S(f)|²
* OptimalFilter(f)= -------  			--- = NormalizationFactor = SUM ---------
*                    N²(f)               k                               |N(f)|²
*
* (*) Boyce et all, "Design and performance of the Astro-E/XRS signal processing system"
*
* - FFT calculus of the matched filter (filter template)
* 	- Declare variables
* 	- Complex FFT values for positive and negative f's
* 	- FFT calculus
* 	- Generation of the f's (positive and negative)
* 	- Magnitude and argument for positive and negative f's
* - N(f)
* - To divide MatchedFilter(f)/N^2(f) => MatchedFilter(f) and N(f) must have the same number of points
* 	- 'if (mf_size' < freqgsl->size)' => Decimate noise samples
* 		- 'if ((freqgsl->size)%mf_size == 0)'
* 		- 'else' => It is necessary to work only with the positive frequencies in order to not handle the f=0
* 		            N(f) interpolation ('interpolatePOS')
* 	- 'else if (mf_size > freqgsl->size)' => Error: Noise spectrum must have more samples than pulse spectrum
*	- 'else if (mf_size == freqgsl->size)' => It is not necessary to do anything
* - OptimalFilter = MatchedFilter'(f)/N^2(f)
* - Calculus of the normalization factor
* - Inverse FFT (to get the expression of the optimal filter in time domain)
*	- Complex OptimalFilter(f) => Taking into account magnitude (MatchedFilter(f)/N^2(f)) and phase (given by MatchedFilter(f))
* - Free of GSL vectors
*
* Parameters:
* - TorF: 'FilterDomain' = T => 'TorF' = 0
*         'FilterDomain' = F => 'TorF' = 1
* - intermediate: 'intermediate' = 0 => Not write an intermediate file
*                 'intermediate' = 1 => Write an intermediate file
* - mode: 'mode' = 0 => CALIBRATION
*         'mode' = 1 => PRODUCTION
* - matchedfiltergsl: Matched filter associated to the pulse (in general, from the interpolation between two matched filters of the library)
* - mf_size: Matched filter size
* - samprate: Sampling rate
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - freqgsl: Frequency axis of the current noise spectral density (input)
* - csdgsl: Current noise spectral density (input)
* - optimal_filtergsl: Optimal filter in time domain (output)
* - of_f: Frequency axis of the optimal filter spectrum (output)
* - of_FFT: Optimal filter spectrum (absolute values) (output)
* - of_FFT_complex: Optimal filter spectrum (complex values) (output)
* - normalizationFactor: Normalization factor (output)
****************************************/
int calculus_optimalFilter(int TorF, int intermediate, int mode, gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val, gsl_vector *freqgsl, gsl_vector *csdgsl, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, gsl_vector_complex **of_FFT_complex, double *normalizationFactor)
{
	// FFT calculus of the matched filter (filter template)
	// Declare variables
	double SelectedTimeDuration = mf_size/samprate;
	string message = "";

	// Complex FFT values for positive and negative f's
	gsl_vector_complex *mfFFTcomp = gsl_vector_complex_alloc(mf_size);
	gsl_vector_complex *mfFFTcomp_conj = gsl_vector_complex_alloc(mf_size);
	gsl_vector *mf_arg = gsl_vector_alloc(mf_size);				// Argument for positive and negative f's
	gsl_vector *mf_f = gsl_vector_alloc(mf_size);				// Sorted f's according [-fmax,...,0,...,fmax]
	gsl_vector *mf_FFT = gsl_vector_alloc(mf_size);				// Sorted magnitude according [-fmax,...,0,...,fmax]

	// FFT calculus
	if (FFT(matchedfiltergsl,mfFFTcomp,SelectedTimeDuration))
	{
		message = "Cannot run FFT routine to calculate filter template";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (runF0orB0val == 0)		gsl_vector_complex_set(mfFFTcomp,0,gsl_complex_rect(0.0,0.0));

	// Generation of the f's (positive and negative)
	// Here is a table which shows the layout of the array data, and the correspondence between the time-domain data z,
	// and the frequency-domain data x

	// index    z               x = FFT(z)

	// 0        z(t = 0)        x(f = 0)
	// 1        z(t = 1)        x(f = 1/(n Delta))
	// 2        z(t = 2)        x(f = 2/(n Delta))
	// .        ........        ..................
	// n/2      z(t = n/2)      x(f = +1/(2 Delta),
	//                                -1/(2 Delta))
	// .        ........        ..................
	// n-3      z(t = n-3)      x(f = -3/(n Delta))
	// n-2      z(t = n-2)      x(f = -2/(n Delta))
	// n-1      z(t = n-1)      x(f = -1/(n Delta))

	if (mf_size == 1)
	{
		gsl_vector_set(mf_f,0,0);
	}
	else if (mf_size%2 == 0)	//Even
	{
		for (int i=0; i<=(mf_size)/2; i++)
		{
			gsl_vector_set(mf_f,i,i/SelectedTimeDuration);
		}
		for (int i=1; i<(mf_size/2); i++)
		{
			gsl_vector_set(mf_f,i+mf_size/2,(-1.0)*(i+mf_size/2-i*2)/SelectedTimeDuration);
		}
	}
	else	//Odd
	{
		for (int i=0; i<=(mf_size)/2; i++)
		{
			gsl_vector_set(mf_f,i,i/SelectedTimeDuration);
		}
		gsl_vector_set(mf_f,mf_size/2+1,(-1.0)*(mf_size/2)/SelectedTimeDuration);
		for (int i=2; i<=(mf_size/2); i++)
		{
			gsl_vector_set(mf_f,i+mf_size/2,(-1.0)*(1+mf_size/2-i)/SelectedTimeDuration);
		}
	}

	// Magnitude and argument for positive and negative f's
	gsl_vector_complex_absRIFCA(mf_FFT,mfFFTcomp);		// Magnitude
	for (int i=0;i<mf_size;i++)
	{
		gsl_vector_complex_set(mfFFTcomp_conj,i,gsl_complex_conjugate(gsl_vector_complex_get(mfFFTcomp,i)));
	}

	gsl_vector_complex_argIFCA(mf_arg,mfFFTcomp);	    // Argument

	// Free allocate of GSL vectors
	gsl_vector_complex_free(mfFFTcomp);

	// N(f)
	gsl_vector *n_f;
	gsl_vector *n_FFT;

	// To divide MatchedFilter(f)/N^2(f) => MatchedFilter(f) and N(f) must have the same number of points
	gsl_vector_view temp;
	if (mf_size < freqgsl->size)			// Decimate noise samples
	{
		if ((freqgsl->size)%mf_size == 0)
		{
			int timesNoverMF = freqgsl->size/mf_size;
			n_f = gsl_vector_alloc(mf_size);
			n_FFT = gsl_vector_alloc(mf_size);
			for (int i=0;i<n_f->size;i++)
			{
				gsl_vector_set(n_f,i,gsl_vector_get(freqgsl,i*timesNoverMF));
				gsl_vector_set(n_FFT,i,gsl_vector_get(csdgsl,i*timesNoverMF));
			}
		}
		else
		{
			// It is necessary to work only with the positive frequencies in order to not handle the f=0
			int noisePOS_size = floor(freqgsl->size/2);
			gsl_vector *freqgsl_POS = gsl_vector_alloc(noisePOS_size);
			temp = gsl_vector_subvector(freqgsl,1,noisePOS_size);
			gsl_vector_memcpy(freqgsl_POS,&temp.vector);
			gsl_vector *csdgsl_POS = gsl_vector_alloc(noisePOS_size);
			temp = gsl_vector_subvector(csdgsl,1,noisePOS_size);
			gsl_vector_memcpy(csdgsl_POS,&temp.vector);

			// N(f) interpolation
			int mfPOS_size = floor(mf_size/2);
			gsl_vector *n_f_interp = gsl_vector_alloc(mfPOS_size);
			gsl_vector *n_FFT_interp = gsl_vector_alloc(mfPOS_size);
		   	if (interpolatePOS (freqgsl_POS, csdgsl_POS, mfPOS_size, gsl_vector_get(mf_f,1)-gsl_vector_get(mf_f,0), &n_f_interp, &n_FFT_interp))
		   	{
				message = "Cannot run routine interpolate for matched filter interpolation";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

		   	n_f = gsl_vector_alloc(mf_size);
		   	n_FFT = gsl_vector_alloc(mf_size);
		   	gsl_vector_set(n_f,0,0);
		   	gsl_vector_set(n_FFT,0,gsl_vector_get(csdgsl,0));
			for (int i=0;i<n_f_interp->size;i++)
		   	{
		   		gsl_vector_set(n_f,i+1,gsl_vector_get(n_f_interp,i));
		   		gsl_vector_set(n_FFT,i+1,gsl_vector_get(n_FFT_interp,i));
		   	}
			if (mf_size%2 == 0)
			{
				for (int i=0;i<n_f_interp->size-1;i++)
				{
					gsl_vector_set(n_f,n_f->size-i-1,-1.0*gsl_vector_get(n_f_interp,i));
					gsl_vector_set(n_FFT,n_f->size-i-1,gsl_vector_get(n_FFT_interp,i));
				}
			}
			else
			{
				for (int i=0;i<n_f_interp->size;i++)
				{
					gsl_vector_set(n_f,n_f->size-i-1,-1.0*gsl_vector_get(n_f_interp,i));
					gsl_vector_set(n_FFT,n_f->size-i-1,gsl_vector_get(n_FFT_interp,i));
				}
			}

			gsl_vector_free(n_f_interp);
		   	gsl_vector_free(n_FFT_interp);
		   	gsl_vector_free(freqgsl_POS);
		   	gsl_vector_free(csdgsl_POS);
		}
	}
	else if (mf_size > freqgsl->size)		// Error
	{
		message = "Noise must have more samples than pulse";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	else if (mf_size == freqgsl->size)		// It is not necessary to do anything
	{
		n_f = gsl_vector_alloc(freqgsl->size);
		n_FFT = gsl_vector_alloc(freqgsl->size);
		gsl_vector_memcpy(n_f,freqgsl);
		gsl_vector_memcpy(n_FFT,csdgsl);
	}

	// OptimalFilter = MatchedFilter'(f)/N^2(f)
	*of_f = gsl_vector_alloc(mf_f->size);
	*of_FFT = gsl_vector_alloc(mf_f->size);
	gsl_vector_memcpy(*of_f,mf_f);

	gsl_vector *mf_FFT_2 = gsl_vector_alloc(mf_f->size);
	gsl_vector *n_FFT_2 = gsl_vector_alloc(mf_f->size);
	gsl_vector_memcpy(mf_FFT_2,mf_FFT);
	gsl_vector_mul(mf_FFT_2,mf_FFT_2);
	gsl_vector_memcpy(n_FFT_2,n_FFT);
	gsl_vector_mul(n_FFT_2,n_FFT_2);

	gsl_vector_memcpy(*of_FFT,mf_FFT);
	*of_FFT_complex = gsl_vector_complex_alloc(mf_size);
	for (int i=0;i<mf_size;i++)
	{
		gsl_vector_complex_set(*of_FFT_complex,i,gsl_complex_div_real(gsl_vector_complex_get(mfFFTcomp_conj,i),gsl_vector_get(n_FFT_2,i)));
	}

	gsl_vector_complex_absRIFCA(*of_FFT,*of_FFT_complex);

	// Calculus of the normalization factor
	*normalizationFactor = 0;
	for (int i=0; i<mf_f->size; i++)
	{
		*normalizationFactor = *normalizationFactor + gsl_vector_get(mf_FFT_2,i)/gsl_vector_get(n_FFT_2,i);
	}

	gsl_vector_free(mf_FFT_2);
	gsl_vector_free(n_FFT_2);

	if ((TorF == 0) || (intermediate == 1) || (mode == 0))
	{
		// Inverse FFT (to get the expression of the optimal filter in time domain)
		// Complex OptimalFilter(f) => Taking into account magnitude (MatchedFilter(f)/N^2(f)) and phase (given by MatchedFilter(f))
		gsl_vector_complex *of_FFTcomp = gsl_vector_complex_alloc(mf_size);
		*optimal_filtergsl = gsl_vector_alloc(mf_size);
		for (int i=0;i<mf_size;i++)
		{
			gsl_vector_complex_set(of_FFTcomp,i,gsl_complex_polar(gsl_complex_abs(gsl_vector_complex_get(*of_FFT_complex,i)),gsl_vector_get(mf_arg,i)));
		}
		if (FFTinverse(of_FFTcomp,*optimal_filtergsl,SelectedTimeDuration))
		{
			message = "Cannot run routine FFTinverse to get optimal filter in time domain";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		gsl_vector_complex_free(of_FFTcomp);
	}
	else
	{
		*optimal_filtergsl = gsl_vector_alloc(mf_size);
		gsl_vector_set_zero(*optimal_filtergsl);
	}

	// Free of GSL vectors
	gsl_vector_complex_free(mfFFTcomp_conj);
	gsl_vector_free(mf_f);
	gsl_vector_free(mf_FFT);
	gsl_vector_free(mf_arg);
	gsl_vector_free(n_f);
	gsl_vector_free(n_FFT);

	return(EPOK);
}
/*xxxx end of SECTION B1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B2 ************************************************************
* interpolatePOS: This function interpolates an input vector ('x_in', 'y_in'), creating an output vector ('x_out', 'y_out') with the size and
*                 frequency step given.
*                 POS is due to the fact that the input spectrum only has positive frequencies (in order to not handle the f=0).
*
* - Declare and initialize variables
* - Method applied to interpolate
* - Generate the interpolated output vector
* - Free memory
*
* Parameters:
* - x_in, y_in: Input vector which is going to be interpolated
* - size: Size of the interpolated output vector
* - step: Frequency step of the interpolated output vector
* - x_out, y_out: Interpolated output vector
****************************************/
int interpolatePOS (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out)
{
	// Declare variables
	gsl_vector_set_zero(*x_out);
	gsl_vector_set_zero(*y_out);
	int N = x_in->size;
	double x[N], y[N];
	for (int i=0; i<N; i++)
	{
		x[i] = gsl_vector_get(x_in,i);
		y[i] = gsl_vector_get(y_in,i);
	}
	double xi, yi;

	// Method applied to interpolate
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	const gsl_interp_type *t = gsl_interp_cspline;

	gsl_spline *spline = gsl_spline_alloc (t, N);

	gsl_spline_init (spline, x, y, N);

	// Generate the interpolated output vector
	xi = step;
	for (int i=0; i<size; i++)
	{
		if ((i == size-1) && (xi <= x[N-1]+x[N-1]/1000)) xi = x[N-1];
		if ((xi >= x[0]) && (xi <= x[N-1]))
		{
		    yi = gsl_spline_eval (spline, xi, acc);

		    gsl_vector_set(*x_out,i,xi);
		    gsl_vector_set(*y_out,i,yi);
		}

	    xi = xi+step;
	}

	// Free memory
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

	return EPOK;
}
/*xxxx end of SECTION B2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B3 ************************************************************
* find_matchedfilter: This function chooses the proper matched filter of the pulse templates library ('MF' or 'MFB0')
*                     by comparing the 'maxDER' to the 'maxDERs'.
*
* It finds the two 'maxDERs' closer in the pulse models library and interpolates ('interpolate_model')
*   - If 'maxDER' is lower than the lower 'maxDERs' in the pulse models library => The model with
*     a lower 'maxDERs' in the pulse models library is chosen
*   - If 'maxDER' is higher than the higher 'maxDERs' in the pulse models library => The model with
*     a higher 'maxDERs' in the pulse models library is chosen
*
* Parameters:
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - maxDER: 'maxDER' of the filtered pulse whose matched filter is looking for
* - maxDERs: 'maxDERs' to compare with
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses the info in the library ('matched_filters')
* - matchedfilterFound: Found matched filter
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_matchedfilter(int runF0orB0val, double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	long nummodels = maxDERs->size;

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		if (runF0orB0val == 0)	gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[0].mfilter);
		else    		gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters_B0[0].mfilter);
		
		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
		if (runF0orB0val == 0)	gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[nummodels-1].mfilter);
		
		*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
		*Ebeta = 0.0;
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				if (runF0orB0val == 0)	gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[i].mfilter);
				else			gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters_B0[i].mfilter);

				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);
				//cout<<"Ealpha: "<<*Ealpha<<endl;
				//cout<<"Ebeta: "<<*Ebeta<<endl;

				// Interpolate between the two corresponding rows in "models"
				gsl_vector *matchedfilterAux = gsl_vector_alloc(reconstruct_init->library_collection->matched_filters[0].mfilter_duration);
				gsl_vector_set_zero(matchedfilterAux);
				if (runF0orB0val == 0)
				{
					if (interpolate_model(&matchedfilterAux,maxDER,reconstruct_init->library_collection->matched_filters[i].mfilter,gsl_vector_get(maxDERs,i),
						reconstruct_init->library_collection->matched_filters[i+1].mfilter,gsl_vector_get(maxDERs,i+1)))
					{
						message = "Cannot run interpolate_model routine for model interpolation";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}
				}
				else
				{
					if (interpolate_model(&matchedfilterAux,maxDER,reconstruct_init->library_collection->matched_filters_B0[i].mfilter,gsl_vector_get(maxDERs,i),
						reconstruct_init->library_collection->matched_filters_B0[i+1].mfilter,gsl_vector_get(maxDERs,i+1)))
					{
						message = "Cannot run interpolate_model routine for model interpolation";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}
				}
				gsl_vector_memcpy(*matchedfilterFound,matchedfilterAux);

				gsl_vector_free(matchedfilterAux);
				
				break;
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION B3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B3 ************************************************************
* find_matchedfilterDAB: This function chooses the proper matched filter of the pulse templates library ('DAB')
*                        by comparing the 'maxDER' to the 'maxDERs'.
*
* It finds the two 'maxDERs' closer in the pulse models library and interpolates ('interpolate_model')
*   - If 'maxDER' is lower than the lower 'maxDERs' in the pulse models library => The model with
*     a lower 'maxDERs' in the pulse models library is chosen
*   - If 'maxDER' is higher than the higher 'maxDERs' in the pulse models library => The model with
*     a higher 'maxDERs' in the pulse models library is chosen
*
* Parameters:
* - runF0orB0val: 'FilterMethod' = F0 => 'runF0orB0val' = 1
*                 'FilterMethod' = B0 => 'runF0orB0val' = 0
* - maxDER: 'maxDER' of the filtered pulse whose matched filter is looking for
* - maxDERs: 'maxDERs' to compare with
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses the info in the library ('matched_filters')
* - matchedfilterFound: Found matched filter
* - PabFound: PAB column from the library (only provided if 'OFInterp'=DAB )
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_matchedfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	long nummodels = maxDERs->size;

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[0].mfilter);
	
		gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,0);
				
		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
		//cout<<"maxDER: "<<maxDER<<endl;
		//cout<<"gsl_vector_get(maxDERs,nummodels-1): "<<gsl_vector_get(maxDERs,nummodels-1)<<endl;
		gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[nummodels-1].mfilter);
		
		message = "In spite of supposedly having a higher energy than the last row of the library, the PAB data of the penultimate row of the library is going to be used";
		EP_PRINT_ERROR(message,-999);
		
		gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,nummodels-2);	//!!!!!! -2 !!!!!!!

		*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-2);
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[i].mfilter);
				
				gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,i);

				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);

				gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[i].mfilter);
				
				gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,i);

				break;
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION B3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B4 ************************************************************
* find_optimalfilter: This function chooses the proper optimal filter of the pulse templates library (OF)
*                     by comparing the 'maxDER' to the 'maxDERs'.
*
* It finds the two 'maxDERs' closer in the pulse models library and interpolates ('interpolate_model')
*   - If 'maxDER' is lower than the lower 'maxDERs' in the pulse models library => The model with
*     a lower 'maxDERs' in the pulse models library is chosen
*   - If 'maxDER' is higher than the higher 'maxDERs' in the pulse models library => The model with
*     a higher 'maxDERs' in the pulse models library is chosen
*
* Parameters:
* - maxDER: 'maxDER' of the filtered pulse whose optimal filter is looking for
* - maxDERs: 'maxDERs' to compare with
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses the info in the library ('optimal_filters')
* - optimalfilterFound: Found optimal filter
* - nrmfctrFound: Found normalization factor
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_optimalfilter(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, double *nrmfctrFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	long nummodels = maxDERs->size;

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		gsl_vector_memcpy(*optimalfilterFound,reconstruct_init->library_collection->optimal_filters[0].ofilter);
		*nrmfctrFound = gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,0);

		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
		gsl_vector_memcpy(*optimalfilterFound,reconstruct_init->library_collection->optimal_filters[nummodels-1].ofilter);
		*nrmfctrFound = gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,nummodels-1);

		*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
		*Ebeta = 0.0;
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				gsl_vector_memcpy(*optimalfilterFound,reconstruct_init->library_collection->optimal_filters[i].ofilter);
				*nrmfctrFound = gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,i);

				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);

				double factor1 = (gsl_vector_get(maxDERs,i+1)-maxDER)/(gsl_vector_get(maxDERs,i+1)-gsl_vector_get(maxDERs,i));
				double factor2 = (maxDER-gsl_vector_get(maxDERs,i))/(gsl_vector_get(maxDERs,i+1)-gsl_vector_get(maxDERs,i));
				*nrmfctrFound = factor1*gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,i)+factor2*gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,i+1);

				// Interpolate between the two corresponding rows in "models"
				gsl_vector *optimalfilterAux = gsl_vector_alloc(reconstruct_init->library_collection->optimal_filters[0].ofilter_duration);
				gsl_vector_set_zero(optimalfilterAux);
				if (interpolate_model(&optimalfilterAux,maxDER,reconstruct_init->library_collection->optimal_filters[i].ofilter,gsl_vector_get(maxDERs,i),
					reconstruct_init->library_collection->optimal_filters[i+1].ofilter,gsl_vector_get(maxDERs,i+1)))
				{
					message = "Cannot run interpolate_model routine for model interpolation";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}

				gsl_vector_memcpy(*optimalfilterFound,optimalfilterAux);

				gsl_vector_free(optimalfilterAux);
							
				break;
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION B4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B4 ************************************************************
* find_optimalfilterDAB: This function chooses the proper optimal filter of the pulse templates library ('OFAB')
*                        by comparing the 'maxDER' to the 'maxDERs'.
*
* It finds the two 'maxDERs' closer in the pulse models library and interpolates ('interpolate_model')
*   - If 'maxDER' is lower than the lower 'maxDERs' in the pulse models library => The model with
*     a lower 'maxDERs' in the pulse models library is chosen
*   - If 'maxDER' is higher than the higher 'maxDERs' in the pulse models library => The model with
*     a higher 'maxDERs' in the pulse models library is chosen
*
* Parameters:
* - maxDER: 'maxDER' of the filtered pulse whose optimal filter is looking for
* - maxDERs: 'maxDERs' to compare with
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses the info in the library ('optimal_filters')
* - optimalfilterFound: Found optimal filter
* - nrmfctrFound: Found normalization factor
* - PabFound: PAB column from the library (only provided if 'OFInterp'=DAB )
* - Ealpha: Energy (in eV) which straddle the 'maxDER' in the lower limit
* - Ebeta: Energy (in eV) which straddle the 'maxDER' in the higher limit
****************************************/
int find_optimalfilterDAB(double maxDER, gsl_vector *maxDERs, ReconstructInitSIRENA *reconstruct_init, gsl_vector **optimalfilterFound, double *nrmfctrFound, gsl_vector **PabFound, double *Ealpha, double *Ebeta)
{
	string message = "";
	long nummodels = maxDERs->size;

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		gsl_vector_memcpy(*optimalfilterFound,reconstruct_init->library_collection->optimal_filters[0].ofilter);
		*nrmfctrFound = gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,0);

		gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,0);
		
		*Ealpha = 0.0;
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,0);
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
		//cout<<"maxDER: "<<maxDER<<endl;
		//cout<<"gsl_vector_get(maxDERs,nummodels-1): "<<gsl_vector_get(maxDERs,nummodels-1)<<endl;
		gsl_vector_memcpy(*optimalfilterFound,reconstruct_init->library_collection->optimal_filters[nummodels-1].ofilter);
		
		message = "In spite of supposedly having a higher energy than the last row of the library, the PAB data of the penultimate row of the library is going to be used";
		EP_PRINT_ERROR(message,-999);
		
		gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,nummodels-2);	//!!!!!! -2 !!!!!!!

		*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-2);
		*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,nummodels-1);
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				gsl_vector_memcpy(*optimalfilterFound,reconstruct_init->library_collection->optimal_filters[i].ofilter);
				*nrmfctrFound = gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,i);
				
				gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,i);
				
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i);

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*Ealpha = gsl_vector_get(reconstruct_init->library_collection->energies,i);
				*Ebeta = gsl_vector_get(reconstruct_init->library_collection->energies,i+1);

				gsl_vector_memcpy(*optimalfilterFound,reconstruct_init->library_collection->optimal_filters[i].ofilter);
				*nrmfctrFound = gsl_vector_get(reconstruct_init->library_collection->nrmfctrs,i);

				gsl_matrix_get_row(*PabFound,reconstruct_init->library_collection->PAB,i);
								
				break;
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION B4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B5 ************************************************************
* find_Esboundary: This function provides the indexes of the two energies which straddle the pulse, by
*                  comparing 'maxDER' of the pulse to the 'maxDERs' of the temples in the library.
*
* It finds the two 'maxDERs' closer in the pulse models library and interpolates ('interpolate_model')
* 	- If 'maxDER' is lower than the lower 'maxDERs' in the pulse models library => indexEalpha = indexEbeta = 0
*   - If 'maxDER' is higher than the higher 'maxDERs' in the pulse models library => indexEalpha = indexEbeta = Number of templates -1
*
* Parameters:
* - maxDER: Maximum of the derivative of the filtered pulse which is being analyzed
* - maxDERs: Maximum of the derivative of the filtered templates in the library
* - indexEalpha: Index of the energy lower than the energy of the pulse which is being analyzed
* - indexEbeta: Index of the energy higher than the energy of the pulse which is being analyzed
****************************************/
int find_Esboundary(double maxDER, gsl_vector *maxDERs, int *indexEalpha, int *indexEbeta)
{
	string message = "";
	int nummodels = maxDERs->size;

	if (maxDER < gsl_vector_get(maxDERs,0))
	{
		*indexEalpha = 0;
		*indexEbeta = 0;
	}
	else if (maxDER > gsl_vector_get(maxDERs,nummodels-1))
	{
		*indexEalpha = nummodels - 1;
		*indexEbeta = nummodels - 1;
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (maxDER == gsl_vector_get(maxDERs,i))
			{
				*indexEalpha = i;
				*indexEbeta = i;

				break;
			}
			else if ((maxDER > gsl_vector_get(maxDERs,i)) && (maxDER < gsl_vector_get(maxDERs,i+1)))
			{
				*indexEalpha = i;
				*indexEbeta = i+1;
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION B5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B6 ************************************************************
* pulseGrading: This function provides the pulse grade (HighRes=0, MidRes=1, LowRes=-1) and the optimal filter length
*               by taking into account the 'PixelType' and the 'OFlength_strategy' (FREE, BASE2, BYGRADE or FIXED).
*
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses the 'PixelType' and the 'OFLength'
* - grade1: Pulse duration
* - grade2_1: tstart of the pulse minus the tstart of the previous pulse
* - OFlength_strategy: 'OFStrategy' (input)
* - pulseGrade: Pulse grade (output)
* - OFlength: Optimal filter length (='OFLength' only if 'OFStrategy'=FIXED) (output)
****************************************/
int pulseGrading (ReconstructInitSIRENA *reconstruct_init, int grade1, int grade2_1, int OFlength_strategy, int *pulseGrade, long *OFlength)
{
	// SPA
	//			 Next pulse  Previous pulse
	// HighRes    >= 512(H1)   >= 128(H2)
	// MidRes     >= 128(M1)   >= 128(M2)
	// LowRes                  >= 128(L2)
	// Rejected                < 128

	int L2;
	int H1;
	int M1;
	if (strcmp(reconstruct_init->PixelType,"SPA") == 0)
	{
		L2 = 128;
		H1 = 512;
		M1 = 128;
	}
	else if (strcmp(reconstruct_init->PixelType,"LPA1") == 0)
	{
		L2 = 400;
		H1 = 1024;
		M1 = 256;
	}
	else if (strcmp(reconstruct_init->PixelType,"LPA2") == 0)
	{
		L2 = 800;
		H1 = 16384;
		M1 = 512;
	}
	else if (strcmp(reconstruct_init->PixelType,"LPA3") == 0)
	{
		L2 = 1400;
		H1 = 16384;
		M1 = 1024;
	}

	// pulseGrade
	// HighRes=0, MidRes=1, LowRes=-1

	if (grade2_1 < L2)
	{
		*pulseGrade = -1;
		if (OFlength_strategy == 0) 		*OFlength = grade1;
		else if (OFlength_strategy == 1) 	*OFlength = pow(2,floor(log2(grade1)));
		//else if (OFlength_strategy == 2) 	*OFlength = L2;
		else if (OFlength_strategy == 2) 	*OFlength = grade1;
		else if (OFlength_strategy == 3) 	*OFlength = reconstruct_init->OFLength;
	}
	else
	{
		if (grade1 >= H1)
		{
			*pulseGrade = 0;
			if (OFlength_strategy == 0)		*OFlength = grade1;
			else if (OFlength_strategy == 1) 	*OFlength = pow(2,floor(log2(grade1)));
			else if (OFlength_strategy == 2) 	*OFlength = H1;
			else if (OFlength_strategy == 3) 	*OFlength = reconstruct_init->OFLength;
		}
		else if (grade1 >= M1)
		{
			*pulseGrade = 1;
			if (OFlength_strategy == 0) 		*OFlength = grade1;
			else if (OFlength_strategy == 1) 	*OFlength = pow(2,floor(log2(grade1)));
			else if (OFlength_strategy == 2)	*OFlength = H1/4;
			else if (OFlength_strategy == 3) 	*OFlength = reconstruct_init->OFLength;
		}
		else
		{
			*pulseGrade = -1;
			if (OFlength_strategy == 0)		*OFlength = grade1;
			else if (OFlength_strategy == 1) 	*OFlength = pow(2,floor(log2(grade1)));
			//else if (OFlength_strategy == 2) 	*OFlength = L2;
			else if (OFlength_strategy == 2) 	*OFlength = grade1;
			else if (OFlength_strategy == 3) 	*OFlength = reconstruct_init->OFLength;
		}
	}

	return EPOK;
}
/*xxxx end of SECTION B6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B7 ************************************************************
* calculateEnergy function: This function calculates the energy of a pulse ('vector') depending on
*                           the 'EnergyMethod' and the 'FilterDomain' basically (and 'OFInterp' if 'OPTFILT', 'I2R' or 'I2RBIS').
*
* OPTFILT (= I2R) or I2RBIS: Optimal filter = Wiener filter
*
*   Once the filter template has been created ('filter' or 'filterFFT'), pulse height analysis is performed by aligning the template
*   with a pulse and multiplying each point in the template by the corresponding point in the pulse. The sum of these products is the energy.
*
* 	Time domain: E = (1/nrmfctr)SUM p(t)·op(t)
*
* 	Frequency domain: E = (1/nrmfctr)SUM P(f)·OF(f)
*
* 	In practice, the alignment of the pulse relative to the trigger is not completely accurate, so a number of n lags could be used in
* 	order to find the peak value of the energy. The n peak values are fitted to a parabola to find the most accurate energy.
*
* 	(*) IXO Onboard processing Trade-Off
* 
* 	If 'OFInterp'=DAB, E = (1/nrmfctr)SUM {(p(t)-PAB(t))·op(t)} or E = (1/nrmfctr)SUM {(P(f)-PAB(f))·OF(f)}
*
* WEIGHT:
*
* 	It is an algorithm for the least-squares measure of pulse size in the presence of non-stationary system noise through the signal event, which
* 	is optimal in both the linear and nonlinear regime. In the case of a linear detector with stationary noise, the algorithm reduces to the Wiener filter.
*
* 	On the transition the TES response is approximately linear, but large changes in resistance can result in a significant change in noise level.
* 	So,	modelling the response is not enough. To get the optimum performace, we must also model the system noise. We describe the noise by its covariance
* 	matrix <didj>, which for an ensemble of time series of length n, is a nxn matrix.
*
* 	Consider a set of measurements of the current in the dectector element, Si. An average of these is used as a model template, Mi=<Si>=(1/N)SUM(p=1,N){Sip}
* 	(Mi is the i-sample of the pulseaverage). And the deviationsfrom this mean, Di=Si-Mi, are used to construct a covariance matrix, Vij=<DiDj> (weight
* 	matrix W=1/V).
*
*	Given a number of models, M, with their associated weight matrices W, only a crude estimation of signal size is sufficient to determine which two
*	calibration points alpha and beta straddle the unkonwn signal U (the pulse whose energy we want to calculate).
*
*	With a linear interpolation of the signal and weight matrix the best energy estimate is:
*
*   E = Ealpha + (Ebeta-Ealpha)(r/3){(2DZ-1)+sqrt[(2DZ-1)²+3(2DY-DXD)/r]}
*
*   where D=U-Salpha. The terms T=Sbeta-Salpha, t=TWalphaT, X=(Wbeta-Walpha)/t, Y=WalphaT/t, Z=XT and r=1/(ZT) can be precalculated with the calibration
*   data alone.
*
* 	(*) Fixsen et all, "Pulse estimation in nonllinear detectors with nonstationary noise"
*
* WEIGHTN:
*
*	The starting idea is the same as in WEIGHT, i.e. minimizing (S-M)W(S-M) but without interpolating W. To do that, a first order development of the pulse
*	shape between 2 calibration points:
*
*	P(E) = P(Ea) + (E-Ea)/(Eb-Ea) . (P(Eb)-P(Ea)) => P(E) - P(Ea) + Ea/(Eb-Ea) . (P(Eb)-P(Ea)) = E . (P(Eb) - P(Ea))/(Eb - Ea)
*
*   => P(E) - Pab = E . Dab => Datos = E · Modelo => y = E·x (condition equation)
*
*   where Pab = P(Ea) - Ea/(Eb-Ea) * (P(Eb)-P(Ea))  and Dab = (P(Eb) - P(Ea))/(Eb - Ea) can be pre-calculated. That way, you see that if you subtract Pab
*   to your data, you end up again in an optimal filter like situation with your data modeled by something that is proportional to a template
*   (that here would be more accurately called a "differential" template).
*   Anyway, this is now a linear chi^2 problem which has an exact solution:
*
*   y = E·X => y0 = E·x0 => SUM{xiyi} = E·SUM{xi²} => Matricial expression => chi^2 taking into account the errors
*              y1 = E·x1                              X'Y = E·X'X             chi^2 = SUM{(data-model)²}/sigma²
*              ...                                    E = (1/(X'X))X'Y        y/sigma = E·x/sigma
*              yN = E·xN                                                      E = (1/(X'WX))X'WY
*
*   E = (Dab'.W.Dab)^(-1) . Dab'.W.(P-Pab) with ' being the transposition operation.
*
* Parameters:
* - vector: Pulse whose energy is looking for
* - pulseGrade: Grade of the input pulse (in order to apply a more complicated method or not to calculate the energy)
* - filter: Optimal filter in time domain
* - filterFFT: Optimal filter in frequency domain
* - runEMethod: 'EnergyMethod' = OPTFILT => 'runEMethod' = 0
* 				'EnergyMethod' = I2R => 'runEMethod' = 0
* 				'EnergyMethod' = I2RBIS => 'runEMethod' = 0
* 				'EnergyMethod' = WEIGHT => 'runEMethod' = 1
* 				'EnergyMethod' = WEIGHTN => 'runEMethod' = 2
* - indexEalpha: Index of the energy lower than the energy of the pulse which is being analyzed
* - indexEbeta: Index of the energy higher than the energy of the pulse which is being analyzed
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
* - domain: 'FilterDomain' = T => 'domain' = 0
*           'FilterDomain' = F => 'domain' = 1
* - nrmfctr: Normalization factor to make the estimate of H in units of energy
* - samprate: Sampling rate
* - calculatedEnergy: Calculated energy in eV
****************************************************************************/
int calculateEnergy (gsl_vector *vector, int pulseGrade, gsl_vector *filter, gsl_vector_complex *filterFFT,int runEMethod, int indexEalpha, int indexEbeta, ReconstructInitSIRENA *reconstruct_init, int domain, double nrmfctr, double samprate, gsl_vector *Pab, double *calculatedEnergy)
{
	string message = "";

	double LagsOrNot = reconstruct_init->LagsOrNot;

	if ((runEMethod == 0) || (runEMethod == 3))
	//       OPTFILT           I2R => OPTFILT, I2RBIS is only slighhtly different from OPTFILT
	{
		if (strcmp(reconstruct_init->OFInterp,"DAB") == 0)	// DAB
		{
			gsl_vector_scale(Pab,-1.0);
			gsl_vector_add(vector,Pab);
		}
			
		int numlags;
		gsl_vector *lags_vector;
		if ((LagsOrNot == 0) || (pulseGrade < 0))	//NOLAGS or LowRes
		{
			numlags = 1;
			lags_vector = gsl_vector_alloc(numlags);
			gsl_vector_set(lags_vector,0,0);
		}
		else
		{
			numlags = 5; // Must be odd
			lags_vector = gsl_vector_alloc(numlags);
			for (int i=0;i<numlags;i++)
			{
				gsl_vector_set(lags_vector,i,-numlags/2+i);
			}
		}
		//cout<<"numlags: "<<numlags<<endl;

		int index_vector;
		gsl_vector *calculatedEnergy_vector = gsl_vector_alloc(numlags);
		gsl_vector_set_zero(calculatedEnergy_vector);
		double a,b,c;
		double xmax;

		double SelectedTimeDuration = vector->size/samprate;
		*calculatedEnergy = 0.0;

		if (domain == 0)	// Time domain filtering
		{
			if (vector->size != filter->size) *calculatedEnergy = 0.0;
			else
			{
				for (int j=0;j<numlags;j++)
				{
					for (int i=gsl_vector_get(lags_vector,j); i<vector->size+gsl_vector_get(lags_vector,j); i++)
					{
						if (i < 0)						index_vector = 0;
						else if (i > vector->size-1)	index_vector = vector->size-1;
						else 							index_vector = i;

						gsl_vector_set(calculatedEnergy_vector,j,gsl_vector_get(calculatedEnergy_vector,j)+gsl_vector_get(vector,index_vector)*gsl_vector_get(filter,i-gsl_vector_get(lags_vector,j)));
						//cout<<i<<" "<<gsl_vector_get(vector,index_vector)<<" "<<gsl_vector_get(filter,i-gsl_vector_get(lags_vector,j))<<" "<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
					}
					gsl_vector_set(calculatedEnergy_vector,j,(gsl_vector_get(calculatedEnergy_vector,j)/nrmfctr)*2*SelectedTimeDuration);
					//cout<<"lag="<<gsl_vector_get(lags_vector,j)<<", E="<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
				}

				if ((pulseGrade < 0) || (LagsOrNot == 0))		*calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
				else
				{
					if (polyFit (lags_vector, calculatedEnergy_vector, &a, &b, &c))
					{
						message = "Cannot run routine polyFit";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}
					//cout<<"a="<<a<<", b="<<b<<" , c="<<c<<endl;
					xmax = -b/(2*a);
					*calculatedEnergy = a*pow(xmax,2.0) + b*xmax +c;
					//cout<<"xmax="<<xmax<<endl;
				}
				//cout<<"calculatedEnergy="<<*calculatedEnergy<<endl;
			}
		}
		else if (domain == 1)	// Frequency domain filtering (multiply vectorFFT and filterFFT)
		{
			if (vector->size != filter->size) *calculatedEnergy = 0.0;
			else
			{
				// Declare variables
				double argExp;
				gsl_vector_complex *calculatedEnergy_vectorcomplex = gsl_vector_complex_alloc(numlags);
				gsl_vector_complex_set_zero(calculatedEnergy_vectorcomplex);

				gsl_vector_complex *vectorFFT = gsl_vector_complex_alloc(vector->size);

				if (FFT(vector,vectorFFT,SelectedTimeDuration))
				{
					message = "Cannot run routine FFT";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}

				*calculatedEnergy = 0.0;

				for (int j=0;j<numlags;j++)
				{
					for (int i=0; i<vector->size; i++)
					{
						argExp = -2*pi*gsl_vector_get(lags_vector,j)*i/vector->size;
						gsl_vector_complex_set(calculatedEnergy_vectorcomplex,j,gsl_complex_add(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j),gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_complex_mul(gsl_vector_complex_get(filterFFT,i),gsl_complex_rect(cos(argExp),sin(argExp))))));
						//cout<<i<<" "<<argExp<<" "<<gsl_complex_abs(gsl_vector_complex_get(vectorFFT,index_vector))<<" "<<gsl_complex_abs(gsl_vector_complex_get(filterFFT,index_vector))<<" "<<GSL_REAL(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))<<","<<GSL_IMAG(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))<<endl;
					}
					gsl_vector_set(calculatedEnergy_vector,j,gsl_complex_abs(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))/nrmfctr);
					//cout<<"lag="<<gsl_vector_get(lags_vector,j)<<", E="<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
				}

				if ((pulseGrade < 0) || (LagsOrNot == 0))		*calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
				else
				{
					if (polyFit (lags_vector, calculatedEnergy_vector, &a, &b, &c))
					{
						message = "Cannot run routine polyFit";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}
					//cout<<"a="<<a<<", b="<<b<<" , c="<<c<<endl;
					xmax = -b/(2*a);
					*calculatedEnergy = a*pow(xmax,2.0) + b*xmax +c;
					//cout<<"xmax="<<xmax<<endl;
				}
				//cout<<"calculatedEnergy="<<*calculatedEnergy<<endl;
				
				gsl_vector_complex_free(calculatedEnergy_vectorcomplex);
				gsl_vector_complex_free(vectorFFT);
			}
		}

		gsl_vector_free(lags_vector);
		gsl_vector_free(calculatedEnergy_vector);
	}
	else if (runEMethod == 1) //WEIGHT
	{
		int pulselength = reconstruct_init->pulse_length;
		gsl_vector *D = gsl_vector_alloc(vector->size);
		gsl_vector *Salpha = gsl_vector_alloc(vector->size);
		gsl_matrix *X = gsl_matrix_alloc(vector->size,vector->size);
		gsl_vector *Y = gsl_vector_alloc(vector->size);
		gsl_vector *Z = gsl_vector_alloc(vector->size);
		gsl_vector_view tempv;
		gsl_matrix_view tempm;
		double scalar_aux1;
		double scalar_aux2;
		double scalar_aux3;
		gsl_vector *vector_aux = gsl_vector_alloc(vector->size);
		gsl_vector *vector_aux1 = gsl_vector_alloc(reconstruct_init->pulse_length);
		gsl_vector *vector_auxlong = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
		gsl_matrix *matrix_auxlong = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);

		//cout<<"indexEalpha: "<<indexEalpha<<endl;
		//cout<<"indexEbeta: "<<indexEbeta<<endl;
		tempv = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates[indexEalpha].ptemplate,0,vector->size);
		gsl_vector_memcpy(Salpha,&tempv.vector);	// Salpha
		gsl_vector_memcpy(D,vector);
		gsl_vector_sub(D,Salpha);		// D = U - Salpha

		gsl_matrix_get_row(vector_aux1,reconstruct_init->library_collection->Y,indexEalpha);
		tempv = gsl_vector_subvector(vector_aux1,0,vector->size);
		gsl_vector_memcpy(Y,&tempv.vector);														// Y
		gsl_blas_ddot(D,Y,&scalar_aux1);	// DY
		scalar_aux1 = 2*scalar_aux1;		// 2DY
		//cout<<"scalar_aux1_0: "<<scalar_aux1<<endl;

		gsl_matrix_get_row(vector_auxlong,reconstruct_init->library_collection->X,indexEalpha);
		vector2matrix(vector_auxlong,&matrix_auxlong);
		tempm = gsl_matrix_submatrix(matrix_auxlong,0,0,vector->size,vector->size);
		gsl_matrix_memcpy(X,&tempm.matrix);														// X
		gsl_blas_dgemv(CblasNoTrans, 1.0, X, D, 0.0, vector_aux);	// XD
		gsl_blas_ddot(D,vector_aux,&scalar_aux2);	// DXD
		//cout<<"scalar_aux2_0: "<<scalar_aux2<<endl;

		scalar_aux1 = scalar_aux1-scalar_aux2;		// 2DY-DXD
		scalar_aux1 = 3*scalar_aux1/gsl_vector_get(reconstruct_init->library_collection->r,indexEalpha);	// 3(2DY-DXD)/r
		//cout<<"scalar_aux1_1: "<<scalar_aux1<<endl;

		gsl_matrix_get_row(vector_aux1,reconstruct_init->library_collection->Z,indexEalpha);
		tempv = gsl_vector_subvector(vector_aux1,0,vector->size);
		gsl_vector_memcpy(Z,&tempv.vector);		// Z
		gsl_blas_ddot(D,Z,&scalar_aux3);		// DZ
		scalar_aux2 = pow(2*scalar_aux3-1,2.0);		// (2DZ-1)²
		//cout<<"scalar_aux2_1: "<<scalar_aux2<<endl;

		scalar_aux1 = sqrt(scalar_aux1+scalar_aux2);	// sqrt[(2DZ-1)² + 3(2DY-DXD)/r]

		scalar_aux3 = 2*scalar_aux3-1;			// (2DZ-1)

		// (Ebeta-Ealpha)*(r/3)*{(2DZ-1) + sqrt[(2DZ-1)² + 3(2DY-DXD)/r]}
		//cout<<"scalar_aux1_2: "<<scalar_aux1<<endl;
		//cout<<"scalar_aux3: "<<scalar_aux3<<endl;
		//cout<<"gsl_vector_get(reconstruct_init->library_collection->r,indexEalpha): "<<gsl_vector_get(reconstruct_init->library_collection->r,indexEalpha)<<endl;
		//cout<<"gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha): "<<gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha)<<" gsl_vector_get(reconstruct_init->library_collection->energies,indexEbeta):"<<gsl_vector_get(reconstruct_init->library_collection->energies,indexEbeta)<<endl;
		scalar_aux1 = (scalar_aux3 + scalar_aux1)*(gsl_vector_get(reconstruct_init->library_collection->r,indexEalpha)/3)*(gsl_vector_get(reconstruct_init->library_collection->energies,indexEbeta)-gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha));
		//cout<<"scalar_aux1: "<<scalar_aux1<<endl;

		// Ealpha + (Ebeta-Ealpha)*(r/3)*{(2DZ-1) + sqrt[(2DZ-1)² + 3(2DY-DXD)/r]}
		*calculatedEnergy = gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha) + scalar_aux1;

		gsl_vector_free(D);
		gsl_vector_free(Salpha);
		gsl_matrix_free(X);
		gsl_vector_free(Y);
		gsl_vector_free(Z);
		gsl_vector_free(vector_aux);
		gsl_vector_free(vector_aux1);
		gsl_vector_free(vector_auxlong);
		gsl_matrix_free(matrix_auxlong);
	}
	else if (runEMethod == 2) //WEIGHTN
	{
		gsl_vector *Pab = gsl_vector_alloc(reconstruct_init->pulse_length);
		gsl_vector *P_Pab = gsl_vector_alloc(reconstruct_init->pulse_length);
		gsl_vector *Dab = gsl_vector_alloc(reconstruct_init->pulse_length);
		gsl_matrix *Wmalpha = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
		gsl_vector *Wvalpha = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
		gsl_matrix *Wmbeta = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
		gsl_vector *Wvbeta = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
		gsl_matrix *Wm = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
		gsl_vector *WDab = gsl_vector_alloc(reconstruct_init->pulse_length);
		double Dab_transWDab;
		gsl_vector *Dab_mod = gsl_vector_alloc(reconstruct_init->pulse_length);	// MOD = modified
		gsl_vector *WP_Pab = gsl_vector_alloc(reconstruct_init->pulse_length);

		gsl_matrix_get_row(Dab,reconstruct_init->library_collection->DAB,indexEalpha);	// Dab

		gsl_vector_memcpy(P_Pab,vector);
		gsl_matrix_get_row(Pab,reconstruct_init->library_collection->PAB,indexEalpha);	// Pab
		gsl_vector_sub(P_Pab,Pab);							// P-Pab

		gsl_matrix_get_row(Wvalpha,reconstruct_init->library_collection->W,indexEalpha);// Walpha vector
		vector2matrix(Wvalpha,&Wmalpha);						// Walpha matrix
		gsl_matrix_get_row(Wvbeta,reconstruct_init->library_collection->W,indexEbeta);	// Wbeta vector
		vector2matrix(Wvbeta,&Wmbeta);							// Wbeta matrix
		gsl_matrix_memcpy(Wm,Wmalpha);
		gsl_matrix_add(Wm,Wmbeta);
		gsl_matrix_scale(Wm,0.5);							// Wmatrix = (Walpha + Wbeta)/2
		gsl_blas_dgemv(CblasNoTrans, 1.0, Wm, Dab, 0.0, WDab);				// WDab => W·Dab
		gsl_blas_ddot(Dab,WDab,&Dab_transWDab);						// Dab_transWDab = Dab*·W·Dab

		gsl_vector_memcpy(Dab_mod,Dab);
		gsl_vector_scale(Dab_mod,1/Dab_transWDab);					// Dab·(Dab*·W·Dab)^(-1)

		gsl_blas_dgemv(CblasNoTrans, 1.0, Wm, P_Pab, 0.0, WP_Pab);			// (Dab*·W·Dab)^(-1)·Dab*·W·(P-Pab)
		gsl_blas_ddot(Dab_mod,WP_Pab,calculatedEnergy);

		gsl_vector_free(Pab);
		gsl_vector_free(P_Pab);
		gsl_vector_free(Dab);
		gsl_matrix_free(Wmalpha);
		gsl_vector_free(Wvalpha);
		gsl_matrix_free(Wmbeta);
		gsl_vector_free(Wvbeta);
		gsl_matrix_free(Wm);
		gsl_vector_free(WDab);
		gsl_vector_free(Dab_mod);
		gsl_vector_free(WP_Pab);
	}
	//cout<<"*calculatedEnergy: "<<*calculatedEnergy<<endl;

	return EPOK;
}
/*xxxx end of SECTION B7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B8 ************************************************************
* writeFilterHDU: This function runs in PRODUCTION mode and it writes the optimal filter info (in the FILTER HDU) for each pulse
*                 if 'intermediate'=1.
*
* - Declare variables
* - Open intermediate FITS file
* - Create the FILTER HDU is if it is the first pulse
* - Write data
* 	- OPTIMALF column
* 	- OFLENGTH column
* 	- NFCTR column
* 	- FREQ column
* 	- OPTIMALFF column
* 	- ENERGY column
* - Modify the Primary HDU
* - Close intermediate output FITS file if it is necessary
* - Free memory
*
* Parameters:
* - reconstruct_init: Structure containing all the pointers and values to run the reconstruction stage
*                     This function uses 'detectFile', 'clobber' and 'pulse_length'
* - pulse_index: Index of the pulse whose info is going to be written (to know if it is the first pulse)
* - normalizationFactor: Normalization factor
* - energy: Calculate energy
* - optimalfilter: Optimal filter in time domain
* - optimalfilter_f: Frequency axis of the optimal filter spectrum
* - optimalfilter_FFT: Optimal filter spectrum
* - dtcObject: Intermediate file name
* - create: To write the 'CREATOR' keyword
******************************************************************************/
int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init, int pulse_index, double normalizationFactor, double energy, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject, const char *create)
{
	// Declare variables
	string message = "";
	int status = EPOK;

	long totalpulses = 0;

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];

	// Open intermediate FITS file
	if (fits_open_file(dtcObject,dtcName,READWRITE,&status))
	{
		message = "Cannot open output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	// Create the FILTER HDU is if it is the first pulse
	if ( ((*reconstruct_init)->clobber == 1) && (pulse_index == 0))
	{
		strcpy(extname,"FILTER");
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
	}

	strcpy(extname,"FILTER");
	if (fits_movnam_hdu(*dtcObject, ANY_HDU,extname, 0, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(*dtcObject,&totalpulses, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	gsl_matrix *optimalfilter_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
	gsl_matrix *optimalfilter_f_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);
	gsl_matrix *optimalfilter_FFT_matrix = gsl_matrix_alloc(1,(*reconstruct_init)->pulse_length);

	gsl_matrix_set_zero(optimalfilter_matrix);
	gsl_matrix_set_zero(optimalfilter_f_matrix);
	gsl_matrix_set_zero(optimalfilter_FFT_matrix);

	for (int i=0;i<optimalfilter->size;i++)
	{
		gsl_matrix_set(optimalfilter_matrix,0,i,gsl_vector_get(optimalfilter,i));
		gsl_matrix_set(optimalfilter_f_matrix,0,i,gsl_vector_get(optimalfilter_f,i));
		gsl_matrix_set(optimalfilter_FFT_matrix,0,i,gsl_vector_get(optimalfilter_FFT,i));
	}

	// Write data
	IOData obj;
	obj.inObject = *dtcObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"FILTER");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];

	// OPTIMALF column
	obj.iniRow = totalpulses+1;
	obj.endRow = totalpulses+1;
	strcpy(obj.nameCol,"OPTIMALF");
	strcpy(obj.unit," ");
	if (writeFitsComplex (obj,optimalfilter_matrix))
	{
		message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_matrix_free(optimalfilter_matrix);

	// OFLENGTH column
	gsl_vector *oflength = gsl_vector_alloc(1);
	gsl_vector_set(oflength,0,optimalfilter->size);
	strcpy(obj.nameCol,"OFLENGTH");
	if (writeFitsSimple (obj,oflength))
	{
		message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) +
			    " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(oflength);

	// NFCTR column
	gsl_vector *nrmfctr = gsl_vector_alloc(1);
	gsl_vector_set(nrmfctr,0,normalizationFactor);
	strcpy(obj.nameCol,"NFCTR");
	if (writeFitsSimple (obj,nrmfctr))
	{
		message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(nrmfctr);

	// FREQ column
	strcpy(obj.nameCol,"FREQ");
	strcpy(obj.unit,"Hz");
	if (writeFitsComplex (obj,optimalfilter_f_matrix))
	{
		message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_matrix_free(optimalfilter_f_matrix);

	// OPTIMALFF column
	strcpy(obj.nameCol,"OPTIMALFF");
	strcpy(obj.unit," ");
	if (writeFitsComplex (obj,optimalfilter_FFT_matrix))
	{
		message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_matrix_free(optimalfilter_FFT_matrix);

	strcpy(obj.nameTable,"PULSES");

	// ENERGY column
	gsl_vector *energygsl = gsl_vector_alloc(1);
	gsl_vector_set(energygsl,0,energy);
	gsl_vector_scale(energygsl,1.0/1e3);
	strcpy(obj.nameCol,"ENERGY");
	strcpy(obj.unit,"keV");
	if (writeFitsSimple (obj,energygsl))
	{
		message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in PULSES";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(energygsl);

	// Modify the Primary HDU
	if (((*reconstruct_init)->clobber == 1) && (pulse_index == 0))
	{
		strcpy(extname,"Primary");
		if (fits_movabs_hdu(*dtcObject, 1, NULL, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output intermediate file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		string mod1 (string("File MODIFIED by") + ' ' +	(string) create);

		strcpy(keyname,"MOD0");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,NULL,&status))
		{
			message = "Cannot write key " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		if (fits_update_key_longstr(*dtcObject,keyname,keyvalstr,NULL,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	// Close intermediate output FITS file if it is necessary
	if (fits_close_file(*dtcObject,&status))
	{
		message = "Cannot close file " + string(dtcName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	// Free memory
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	return(EPOK);
}
/*xxxx end of SECTION B8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
