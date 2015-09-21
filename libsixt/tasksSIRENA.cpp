#include "tasksSIRENA.h"

/***** SECTION A ************************************************************
* runDetect: This function ... Run detection routine in SIRENA, record by record
*
******************************************************************************/
void runDetect(TesRecord* record, int lastRecord, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
{
	const char * create= "runDetect v.17.0.0";	//Set "CREATOR" keyword of output FITS file

	string message="";
	int status=EPOK;

	// Declare variables
	fitsfile *inLibObject = NULL;	// Object which contains information of the library FITS file
	bool appendToLibrary = false;	// Pulse templates library FITS file new (appendToLibrary=false) or not (appendToLibrary=true)

	fitsfile *dtcObject = NULL;	    // Object which contains information of the output FITS file
	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	int eventsz = record->trigger_size;
	double tstartRecord;
	gsl_vector *invector = gsl_vector_alloc(eventsz);	// Record

	// Handle the library data or the library file
	if (handleLibraryDetect(*reconstruct_init,&appendToLibrary, lastRecord, &inLibObject,create))
	{
		message = "Cannot run routine handleLibraryDetect to handle library";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Create output FITS file: Detect file (*_dtc.fits) if reconstruct_init->intermediate=1
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (createDetectFile(*reconstruct_init,1/record->delta_t,create, &dtcObject))
		{
			message = "Cannot create file " +  string((*reconstruct_init)->detectFile);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Filter and derive the 'models'
	if ((*reconstruct_init)->crtLib == 0)
	{
		if (filderLibrary(reconstruct_init,1/record->delta_t))
		{
			message = "Cannot run routine filderLibrary to filter & derive library if the 1st record";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Store the record in 'invector'
	if (loadRecord(record, &tstartRecord, &invector))
	{
		message = "Cannot run routine loadRecord";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	eventsz = invector->size;	// Just in case the last record has been filled in with 0's => Re-allocate invector

	// Process each record
	if (procRecord(reconstruct_init, tstartRecord, 1/record->delta_t, dtcObject, invector, *pulsesInRecord))
	{
	    message = "Cannot run routine procRecord for record processing";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	int extver=0;
	char extname[20];
	
	if (((*reconstruct_init)->intermediate == 1) && (lastRecord == 1))
	{
		// Write output keywords (their values have been previously checked)
		char keyname[10];
		char *comment=NULL;
		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, extver, &status))
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
		if(fits_write_key(dtcObject,TINT,keyname,&ttpls1,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) +" in " + string(dtcName);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	if ((lastRecord == 1) && (*reconstruct_init)->crtLib == 1 && pulsesAll->ndetpulses>0)	// CREATIONLIB run mode => Calculate the pulse template by averaging some found pulses
	{
		gsl_vector *pulsetemplate = gsl_vector_alloc((*reconstruct_init)->pulse_length);
		double pulseheighttemplate = 0;
		gsl_matrix *weight = gsl_matrix_alloc((*reconstruct_init)->pulse_length,(*reconstruct_init)->pulse_length);
		gsl_matrix *covariance = gsl_matrix_alloc((*reconstruct_init)->pulse_length,(*reconstruct_init)->pulse_length);
		gsl_matrix_set_zero(weight);
		gsl_matrix_set_zero(covariance);

		char extname[20];
		char keyname[10];
		char *comment=NULL;

		if (calculateTemplate (*reconstruct_init, pulsesAll, *pulsesInRecord, 1/record->delta_t, &pulsetemplate, &pulseheighttemplate, &covariance, &weight))
		{
		    message = "Cannot run routine calculateTemplate in creationlib run mode";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		//cout<<"COVAR2 "<<gsl_matrix_get(covariance,0,1)<<" "<<gsl_matrix_get(covariance,0,2)<<" "<<gsl_matrix_get(covariance,0,3)<<endl;
		//cout<<"COVAR2 "<<gsl_matrix_get(covariance,1,0)<<" "<<gsl_matrix_get(covariance,2,0)<<" "<<gsl_matrix_get(covariance,3,0)<<endl;

		if (writeLibrary(*reconstruct_init, pulseheighttemplate, &pulsetemplate, covariance, weight, appendToLibrary, &inLibObject))
		{
		    message = "Cannot run routine writeLibrary in creationlib run mode";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		if (((*reconstruct_init)->lastELibrary == 1) && ((strcmp((*reconstruct_init)->EnergyMethod,"WEIGHT") == 0) || (strcmp((*reconstruct_init)->EnergyMethod,"WN") == 0)))
		{
			if (fillInLibraryData(*reconstruct_init))
			{
				message = "Cannot run routine Library in crationlib run mode";
				EP_EXIT_ERROR(message,EPFAIL);
			}
		}

		gsl_vector_free(pulsetemplate);
		gsl_matrix_free(weight);
		gsl_matrix_free(covariance);
	}

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
* handleLibraryDetect: This function ...
*
******************************************************************************/
int handleLibraryDetect(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, int lastRecord, fitsfile **inLibObject, const char* create)
{
	// Declare variables
	string message = "";
	int crtLib = reconstruct_init->crtLib;

	if (crtLib == 1)
	{
		if (lastRecord == 1)
		{
			if (createLibraryDetect(reconstruct_init, appendToLibrary, inLibObject, create))
			{
				message = "Cannot run routine createLibraryDetect to create pulses library";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION A1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A4 ************************************************************
* createLibraryDetect: This function...
****************************************************************************/
int createLibraryDetect(ReconstructInitSIRENA* reconstruct_init, bool *appendToLibrary, fitsfile **inLibObject, const char * create)
{
	int status = EPOK;
	int extver=0;
	string message = "";

	char extname[20];
	char keyname[10];
	char *comment=NULL;

	char inLibName[256];
	strncpy(inLibName, reconstruct_init->library_file,255);
	inLibName[255]='\0';

	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];
	//int evtcnt=0, ttpls1=0, modeval=0,eventcntLib1=0, lib_id=0;

	long eventcntLib;

	// Create pulse templates library FITS file: If it does not exist yet
	// Create pulse templates library file or open it
	if (fileExists(string(inLibName)))
	{
		*appendToLibrary = true;

		if (fits_open_file(inLibObject,inLibName,READWRITE,&status))
		{
		    message = "Cannot open library file " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(extname,"LIBRARY");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}
	else
	{
		*appendToLibrary = false;
		status = EPOK;
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
		strcpy(extname,"LIBRARY");
		if (fits_create_tbl(*inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
		{
		    message = "Cannot create table " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(extname,"LIBRARY");
		if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		eventcntLib = 1;
		strcpy(keyname,"EVENTCNT");
		if (eventcntLib <= 0)
		{
		    message = "Legal values for EVENTCNT (LIBRARY) are integer numbers greater than 0";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_write_key(*inLibObject,TLONG,keyname,&eventcntLib,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		/*// WEIGHT image extension
		strcpy(extname,"LIBRARY");
		long naxes[3]={reconstruct_init->pulse_length,reconstruct_init->pulse_length,1}
		if (fits_create_img(*inLibObject,TDOUBLE,3,naxes,&status))
		{
		    message = "Cannot create table " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}*/

		// Primary HDU
		strcpy(extname,"Primary");
		int *hdutype;
		if (fits_movabs_hdu(*inLibObject, 1, hdutype, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"CREATOR");
		string creator (string("File CREATED by") + ' ' + (string) create);
		strcpy(keyvalstr,creator.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"PROC0");
		const char * charproc= "PROC0 Starting parameter list";
		strcpy(keyvalstr,charproc);
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		string strproc (string("RecordFile = ") + reconstruct_init->record_file);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("TesEventFile = ") + reconstruct_init->event_file;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("LibraryFile = ") + reconstruct_init->library_file;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("NoiseFile = ") + reconstruct_init->noise_file;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_mode[125];			sprintf(str_mode,"%d",reconstruct_init->mode);
		strproc = string("mode = ") + string(str_mode);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_crtLib[125];		sprintf(str_crtLib,"%d",reconstruct_init->crtLib);
		strproc=string("crtLib = ") + string(str_crtLib);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_lastELibrary[125];		sprintf(str_lastELibrary,"%d",reconstruct_init->lastELibrary);
		strproc=string("lastELibrary = ") + string(str_lastELibrary);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("PixelType = ") + reconstruct_init->PixelType;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("FilterDomain = ") + reconstruct_init->FilterDomain;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("FilterMethod = ") + reconstruct_init->FilterMethod;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_calibLQ[125];      sprintf(str_calibLQ,"%d",reconstruct_init->calibLQ);
		strproc=string("calibLQ = ") + string(str_calibLQ);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_b_cF[125];	sprintf(str_b_cF,"%f",reconstruct_init->b_cF);
		strproc=string("b_cF = ") + string(str_b_cF);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_c_cF[125];	sprintf(str_c_cF,"%f",reconstruct_init->c_cF);
		strproc=string("c_cF = ") + string(str_c_cF);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_maxPulsesPerRecord[125];	sprintf(str_maxPulsesPerRecord,"%d",reconstruct_init->maxPulsesPerRecord);
		strproc=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_pulse_length[125];	sprintf(str_pulse_length,"%d",reconstruct_init->pulse_length);
		strproc=string("PulseLength = ") + string(str_pulse_length);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tauFall[125];		sprintf(str_tauFall,"%e",reconstruct_init->tauFall);
		strproc=string("tauFall = ") + string(str_tauFall);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",reconstruct_init->scaleFactor);
		strproc=string("scaleFactor = ") + string(str_scaleFactor);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_samplesUp[125];	sprintf(str_samplesUp,"%f",reconstruct_init->samplesUp);
		strproc=string("samplesUp = ") + string(str_samplesUp);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_nSgms[125];	    sprintf(str_nSgms,"%f",reconstruct_init->nSgms);
		strproc=string("nSgms = ") + string(str_nSgms);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_baseline[125];	sprintf(str_baseline,"%f",reconstruct_init->baseline);
		strproc=string("baseline = ") + string(str_baseline);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_LrsT[125];			sprintf(str_LrsT,"%e",reconstruct_init->LrsT);
		strproc=string("LrsT = ") + string(str_LrsT);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_LbT[125];			sprintf(str_LbT,"%e",reconstruct_init->LbT);
		strproc=string("LbT = ") + string(str_LbT);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_monoenergy[125];	sprintf(str_monoenergy,"%f",reconstruct_init->monoenergy);
		strproc=string("monoenergy = ") + string(str_monoenergy);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_intermediate[125];      sprintf(str_intermediate,"%d",reconstruct_init->intermediate);
		strproc=string("intermediate = ") + string(str_intermediate);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("detectFile = ") + reconstruct_init->detectFile;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("RecordFileCalib2 = ") + reconstruct_init->record_file2;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_monoenergy2[125];	sprintf(str_monoenergy2,"%f",reconstruct_init->monoenergy2);
		strproc=string("monoenergy2 = ") + string(str_monoenergy2);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strproc=string("filterFile = ") + reconstruct_init->filterFile;
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_clobber[125];      sprintf(str_clobber,"%d",reconstruct_init->clobber);
		strproc=string("clobber = ") + string(str_clobber);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tstartPulse1[125];	sprintf(str_tstartPulse1,"%d",reconstruct_init->tstartPulse1);
		strproc=string("tstartPulse1 = ") + string(str_tstartPulse1);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tstartPulse2[125];	sprintf(str_tstartPulse2,"%d",reconstruct_init->tstartPulse2);
		strproc=string("tstartPulse2 = ") + string(str_tstartPulse2);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tstartPulse3[125];	sprintf(str_tstartPulse3,"%d",reconstruct_init->tstartPulse3);
		strproc=string("tstartPulse3 = ") + string(str_tstartPulse3);
		strcpy(keyvalstr,strproc.c_str());
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		charproc= "PROC0 Ending parameter list";
		strcpy(keyvalstr,charproc);
		if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return (EPOK);
}
/*xxxx end of SECTION A4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A5 ************************************************************
* createDetectFile function: This function ...
***************************************************************************/
int createDetectFile(ReconstructInitSIRENA* reconstruct_init, double samprate, const char * create, fitsfile **dtcObject)
{
	int status = EPOK;
	string message = "";

	char dtcName[256];
	strncpy(dtcName,reconstruct_init->detectFile,255);
	dtcName[255]='\0'; // enforce zero ending string in case of buffer overflows

	// Create output FITS file: If it does not exist yet
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
			char dtcNameaux[255];
			sprintf(dtcNameaux,dtcName);
			strcat(dtcNameaux,".fits");
			strcpy(dtcName,dtcNameaux);
		}
	}

	// Create _dtc file (if file already exists => check clobber)
	if (fileExists(string(dtcName)) && (reconstruct_init->clobber == 1))
	{
		if (remove(dtcName))
		{
			message = "Output detect file already exists & cannot be deleted ("+string(strerror(errno))+")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	    }
	}
	else if (fileExists(string(dtcName)) && (reconstruct_init->clobber == 0))
	{
		message = "Output detect file already exists: must not be overwritten";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	if (!fileExists(string(dtcName)))
	{
		if(fits_create_file(dtcObject, dtcName, &status))
		{
			message = "Cannot create output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		message = "Create Detect Fits File: " + string(dtcName);
	}

	// Create extension PULSES
		// To work with tables (extensions)
	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	int extver = 0;
		// To write keywords
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	if (fits_open_file(dtcObject,dtcName,READWRITE,&status))
	{
	    message = "Cannot open output detect file " + string(dtcName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (reconstruct_init->clobber == 1)
	{
		// PULSES HDU
		strcpy(extname,"PULSES");
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(*dtcObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"MODE");
		if ((reconstruct_init->mode != 0) && (reconstruct_init->mode != 1))
		{
			message = "Legal values for MODE (PULSES) are 0 or 1";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_write_key(*dtcObject,TINT,keyname,&(reconstruct_init->mode),comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(keyname,"EVENTSZ");
		if (reconstruct_init->pulse_length <= 0)
		{
			message = "Legal values for EVENTSZ (PULSES) are integer numbers greater than 0";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_write_key(*dtcObject,TINT,keyname,&(reconstruct_init->pulse_length),comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		if (reconstruct_init->crtLib == 1)
		{
			strcpy(keyname,"ENERGY");
			if (reconstruct_init-> monoenergy < 0)
			{
				message = "Legal values for ENERGY (PULSES) are non negative real numbers";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			if (fits_write_key(*dtcObject,TDOUBLE,keyname,&(reconstruct_init-> monoenergy),comment,&status))
			{
				message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}
		}
		strcpy(keyname,"SAMPRATE");
		if (samprate <= 0)
		{
			message = "Legal values for SAMPRATE (PULSES) are non negative real numbers";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_write_key(*dtcObject,TDOUBLE,keyname,&samprate,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// TEST HDU
		strcpy(extname,"TESTINFO");
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Primary HDU
		strcpy(extname,"Primary");
		int *hdutype;
		if (fits_movabs_hdu(*dtcObject, 1, hdutype, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"CREATOR");
		string creator (string("File CREATED by") + ' ' + (string) create);
		strcpy(keyvalstr,creator.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"HISTORY");
		const char * charhistory= "HISTORY Starting parameter list";
		strcpy(keyvalstr,charhistory);
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		string strhistory (string("RecordFile = ") + reconstruct_init->record_file);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("TesEventFile = ") + reconstruct_init->event_file;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("LibraryFile = ") + reconstruct_init->library_file;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("NoiseFile = ") + reconstruct_init->noise_file;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_mode[125];			sprintf(str_mode,"%d",reconstruct_init->mode);
		strhistory = string("mode = ") + string(str_mode);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_crtLib[125];		sprintf(str_crtLib,"%d",reconstruct_init->crtLib);
		strhistory=string("crtLib = ") + string(str_crtLib);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_lastELibrary[125];		sprintf(str_lastELibrary,"%d",reconstruct_init->lastELibrary);
		strhistory=string("lastELibrary = ") + string(str_lastELibrary);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("PixelType = ") + reconstruct_init->PixelType;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("FilterDomain = ") + reconstruct_init->FilterDomain;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("FilterMethod = ") + reconstruct_init->FilterMethod;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_calibLQ[125];      sprintf(str_calibLQ,"%d",reconstruct_init->calibLQ);
		strhistory=string("calibLQ = ") + string(str_calibLQ);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_b_cF[125];	sprintf(str_b_cF,"%f",reconstruct_init->b_cF);
		strhistory=string("b_cF = ") + string(str_b_cF);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_c_cF[125];	sprintf(str_c_cF,"%f",reconstruct_init->c_cF);
		strhistory=string("c_cF = ") + string(str_c_cF);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_maxPulsesPerRecord[125];	sprintf(str_maxPulsesPerRecord,"%d",reconstruct_init->maxPulsesPerRecord);
		strhistory=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_pulse_length[125];	sprintf(str_pulse_length,"%d",reconstruct_init->pulse_length);
		strhistory=string("PulseLength = ") + string(str_pulse_length);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tauFall[125];		sprintf(str_tauFall,"%e",reconstruct_init->tauFall);
		strhistory=string("tauFall = ") + string(str_tauFall);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",reconstruct_init->scaleFactor);
		strhistory=string("scaleFactor = ") + string(str_scaleFactor);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_samplesUp[125];	sprintf(str_samplesUp,"%f",reconstruct_init->samplesUp);
		strhistory=string("samplesUp = ") + string(str_samplesUp);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_nSgms[125];	    sprintf(str_nSgms,"%f",reconstruct_init->nSgms);
		strhistory=string("nSgms = ") + string(str_nSgms);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_baseline[125];	sprintf(str_baseline,"%f",reconstruct_init->baseline);
		strhistory=string("baseline = ") + string(str_baseline);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_LrsT[125];			sprintf(str_LrsT,"%e",reconstruct_init->LrsT);
		strhistory=string("LrsT = ") + string(str_LrsT);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_LbT[125];			sprintf(str_LbT,"%e",reconstruct_init->LbT);
		strhistory=string("LbT = ") + string(str_LbT);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_monoenergy[125];	sprintf(str_monoenergy,"%f",reconstruct_init->monoenergy);
		strhistory=string("monoenergy = ") + string(str_monoenergy);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_intermediate[125];      sprintf(str_intermediate,"%d",reconstruct_init->intermediate);
		strhistory=string("intermediate = ") + string(str_intermediate);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("detectFile = ") + reconstruct_init->detectFile;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("RecordFileCalib2 = ") + reconstruct_init->record_file2;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_monoenergy2[125];	sprintf(str_monoenergy2,"%f",reconstruct_init->monoenergy2);
		strhistory=string("monoenergy2 = ") + string(str_monoenergy2);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strhistory=string("filterFile = ") + reconstruct_init->filterFile;
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_clobber[125];      sprintf(str_clobber,"%d",reconstruct_init->clobber);
		strhistory=string("clobber = ") + string(str_clobber);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tstartPulse1[125];	sprintf(str_tstartPulse1,"%d",reconstruct_init->tstartPulse1);
		strhistory=string("tstartPulse1 = ") + string(str_tstartPulse1);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tstartPulse2[125];	sprintf(str_tstartPulse2,"%d",reconstruct_init->tstartPulse2);
		strhistory=string("tstartPulse2 = ") + string(str_tstartPulse2);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		char str_tstartPulse3[125];	sprintf(str_tstartPulse3,"%d",reconstruct_init->tstartPulse3);
		strhistory=string("tstartPulse3 = ") + string(str_tstartPulse3);
		strcpy(keyvalstr,strhistory.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in library file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		charhistory= "HISTORY Ending parameter list";
		strcpy(keyvalstr,charhistory);
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

	}

	return EPOK;
}
/*xxxx end of SECTION A5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A6 ************************************************************
* filderLibrary: This function ...
*
******************************************************************************/
int filderLibrary(ReconstructInitSIRENA** reconstruct_init, double samprate)
{
	int status = EPOK;

	// NOTCREATIONLIB mode and first record
	if (((*reconstruct_init)->crtLib == 0) && ((*reconstruct_init)->library_collection->pulse_templates_filder[0].template_duration == -1))
	{
		string message = "";

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

		// 'models' is filtered and derived
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
				message = "Cannot run routine derMTHSimple to calculate derivative";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			gsl_vector_memcpy((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,model);

			(*reconstruct_init)->library_collection->pulse_templates_filder[i].template_duration = (*reconstruct_init)->library_collection->pulse_templates[i].template_duration;

			gsl_vector_set((*reconstruct_init)->library_collection->maxDERs,i,gsl_vector_max(model));
			//cout<<"maxDERs["<<i<<"]="<<gsl_vector_get((*reconstruct_init)->library_collection->maxDERs,i)<<endl;
			//cout<<"maxDERs["<<i<<"]="<<gsl_vector_get((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,0)<<" "<<gsl_vector_get((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,1)<<" "<<gsl_vector_get((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,2)<<" "<<endl;
			/*for (int k=0;k<(*reconstruct_init)->pulse_length;k++)
			{
				cout<<k<<" "<<gsl_vector_get((*reconstruct_init)->library_collection->pulse_templates_filder[i].ptemplate,k)<<endl;
			}*/
		}
		/*cout<<"W: "<<endl;
		cout<<gsl_matrix_get((*reconstruct_init)->library_collection->W,0,0)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->W,0,1)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->W,0,2)<<" "<<endl;
		cout<<gsl_matrix_get((*reconstruct_init)->library_collection->W,1,0)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->W,1,1)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->W,1,2)<<" "<<endl;
		cout<<"T: "<<endl;
		cout<<gsl_matrix_get((*reconstruct_init)->library_collection->T,0,0)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->T,0,1)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->T,0,2)<<" "<<endl;
		cout<<"t: "<<endl;
		cout<<gsl_vector_get((*reconstruct_init)->library_collection->t,0)<<endl;
		cout<<"X: "<<endl;
		cout<<gsl_matrix_get((*reconstruct_init)->library_collection->X,0,0)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->X,0,1)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->X,0,2)<<" "<<endl;
		cout<<"Y: "<<endl;
		cout<<gsl_matrix_get((*reconstruct_init)->library_collection->Y,0,0)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->Y,0,1)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->Y,0,2)<<" "<<endl;
		cout<<"Z: "<<endl;
		cout<<gsl_matrix_get((*reconstruct_init)->library_collection->Z,0,0)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->Z,0,1)<<" "<<gsl_matrix_get((*reconstruct_init)->library_collection->Z,0,2)<<" "<<endl;
		cout<<"r: "<<endl;
		cout<<gsl_vector_get((*reconstruct_init)->library_collection->r,0)<<endl;*/

		gsl_vector_free(model);
	}

	return(EPOK);
}
/*xxxx end of SECTION A6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A7 ************************************************************
* loadRecord: This function loads the struture 'record' in the 'adc_double' vector.
*             The time length of the record is stored in 'time_record'.
*
* It checks if the record has been filled in with 0's => It only loads the first values (which are different from 0).
******************************************************************************/
int loadRecord(TesRecord* record, double *time_record, gsl_vector **adc_double)
{
	*time_record = record->time;
	for (int i=0;i<record->trigger_size;i++)
	{
		gsl_vector_set(*adc_double,i,record->adc_double[i]);
	}

	// Just in case the last record has been filled in with 0's => Re-allocate 'invector'
	if (gsl_vector_ispos(*adc_double) != 1)
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
/*xxxx end of SECTION A7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A8 ************************************************************
* procRecord function:  This function ...
****************************************************************************/
int procRecord(ReconstructInitSIRENA** reconstruct_init, double tstartRecord, double samprate, fitsfile *dtcObject, gsl_vector *record, PulsesCollection *foundPulses)
{
	int status = EPOK;
	string message = "";

	// Initialize variables
	int numPulses = 0;
	double threshold = 0.0;

	double asquid = 1.0;
	double plspolar = 1.0;
	double ivcal=1.0;
	int chngplrt = 0;// 1 => Polarity changed (pulses multiplied by -1)
    				// 0 => Polarity not changed

	double stopCriteriaMKC = 1.0;  			// Used in medianKappaClipping
	                               	   	   	// Given in %
	double kappaMKC = 3.0;					// Used in medianKappaClipping
	double levelPrvPulse = 100.0;  		    // Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse

	gsl_vector_view temp;

	double scaleFactor = (*reconstruct_init)->scaleFactor;
	double tauFALL = (*reconstruct_init)->tauFall;
	int sizePulse_b = (*reconstruct_init)->pulse_length;
	double samplesUp = (*reconstruct_init)->samplesUp;
	double nSgms = (*reconstruct_init)->nSgms;
	double Lrs = (int) ((*reconstruct_init)->LrsT*samprate);			// Running sum length (in the RS filter case): LrsT in samples
	double Lb = (int) ((*reconstruct_init)->LbT*samprate); 				// Baseline averaging length (in the RS filter case): LbT in samples

	// Allocate GSL vectors
	gsl_vector *recordNOTFILTERED = gsl_vector_alloc(record->size); 	// Record without having been filtered
	gsl_vector *recordDERIVATIVE = gsl_vector_alloc(record->size);  	// Derivative of invectorFILTERED

	// To look for pulses
	gsl_vector *tstartgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *tendgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *qualitygsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *pulseHeightsgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector *maxDERgsl = gsl_vector_alloc((*reconstruct_init)->maxPulsesPerRecord);
	gsl_vector_set_zero(qualitygsl);
	gsl_vector_set_zero(pulseHeightsgsl);						// In order to choose the proper pulse model to calculate
	                                                            // the adjusted derivative and to fill in the ESTENRGY column
	                                                            // in the output FITS file
	gsl_vector_set_zero(maxDERgsl);

	gsl_vector_scale(record,ivcal);			// IVCAL to change arbitrary units of voltage to non-arbitrary
			                                // units of current (Amps)

	// Assign positive polarity to the pulses
	if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))
	{
		gsl_vector_scale(record,-1);
		chngplrt = 1;
	}

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	if ((*reconstruct_init)->intermediate == 1)
	{
		if (fits_open_file(&dtcObject,dtcName,1,&status))
		{
			message = "Cannot open file " +  string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		char extname[20];
		char keyname[10];
		int extver=0;
		char *comment=NULL;
		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(keyname,"CHNGPLRT");
		if ((chngplrt != 0) && (chngplrt != 1))
		{
			message = "Legal values for CHNGPLRT (PULSES) are 0 or 1";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_update_key(dtcObject,TINT,keyname,&chngplrt,comment,&status))
		{
			message = "Cannot update keyword " + string(keyname) + " in output detect file ";
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

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
	//if ((*reconstruct_init)->intermediate == 1)
	//{
		gsl_vector *recordDERIVATIVEOriginal = gsl_vector_alloc(recordDERIVATIVE->size);
		gsl_vector_memcpy(recordDERIVATIVEOriginal,recordDERIVATIVE);
	//}
	//cout<<"recordDerivative"<<endl;
	//for (int j=950;j<2010;j++) cout<<j<<" "<<gsl_vector_get(recordNOTFILTERED,j)<<" "<<gsl_vector_get(recordDERIVATIVE,j)<<endl;

	// Find pulses of the record
	if ((*reconstruct_init)->crtLib == 0)
	{
		/*if (findPulsesPROD (recordDERIVATIVE, &tstartgsl, &qualitygsl, &maxDERgsl,
						&numPulses, &threshold,
						tauFALL, scaleFactor, samprate,
						samplesUp, nSgms,
						(*reconstruct_init),
						stopCriteriaMKC,
						kappaMKC))
				{
					message = "Cannot run routine findPulsesPROD";
					EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}*/
		int tstartFirstEvent = 0;
		bool triggerCondition;
		//int tstart = 0;
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

		if (FindSecondaries	((*reconstruct_init)->maxPulsesPerRecord,
				recordDERIVATIVE, threshold,
				samplesUp,(*reconstruct_init),
				tstartFirstEvent,
				&numPulses,&tstartgsl,&qualitygsl, &maxDERgsl))
		{
			message = "Cannot run routine FindSecondaries";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	else if ((*reconstruct_init)->crtLib == 1)
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
	//cout<<"threshold: "<<threshold<<endl;

	// To write test info
	if ((*reconstruct_init)->intermediate == 1)
	{
		if (writeTestInfo((*reconstruct_init), recordNOTFILTERED, recordDERIVATIVEOriginal, threshold, dtcObject))
		{
			message = "Cannot run routine writeTestInfo";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	gsl_vector_free(recordDERIVATIVEOriginal);

	for (int i=0;i<numPulses;i++)
	{
		gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i)+sizePulse_b); 	//tend_i = tstart_i + (ntaus*tauFALL*samprate)

		if (gsl_vector_get(tendgsl,i) >= recordDERIVATIVE->size)	// Truncated pulses at the end of the record
		{
			gsl_vector_set(tendgsl,i,(recordDERIVATIVE->size)-1);
			gsl_vector_set (qualitygsl,i,2);
		}

		if ((numPulses !=1) && (i != numPulses-1)) 		// More than one pulse in the record and not the last one
		{
			if (gsl_vector_get(tendgsl,i) > gsl_vector_get(tstartgsl,i+1))
			{
				gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i+1));
			}
		}
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

		// 'grade1' will be known after running runFilter (but initialize for library creation!)
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

		gsl_vector_memcpy(foundPulses->pulses_detected[i].pulse_adc,&temp.vector);

		foundPulses->pulses_detected[i].Tstart = gsl_vector_get(tstartgsl,i)/samprate+tstartRecord;
		foundPulses->pulses_detected[i].Tend = gsl_vector_get(tendgsl,i)/samprate+tstartRecord;
		foundPulses->pulses_detected[i].riseTime = gsl_vector_get(tauRisegsl,i);
		foundPulses->pulses_detected[i].fallTime = gsl_vector_get(tauFallgsl,i);
		foundPulses->pulses_detected[i].pulse_height = gsl_vector_get(pulseHeightsgsl,i);
		foundPulses->pulses_detected[i].maxDER = gsl_vector_get(maxDERgsl,i);
		// 'ucenergy' and 'energy' will be known after running runFilter and runEnergy respectively
		foundPulses->pulses_detected[i].quality = gsl_vector_get(qualitygsl,i);
		//cout<<"Pulse "<<i<<" tstart="<<gsl_vector_get(tstartgsl,i)<<", maxDER= "<<foundPulses->pulses_detected[i].maxDER<<", pulse_duration= "<<foundPulses->pulses_detected[i].pulse_duration<<",quality= "<<foundPulses->pulses_detected[i].quality<<endl;
		//cout<<gsl_vector_get(tstartgsl,i)<<endl;
	}

	// Write pulses info in output FITS file
	if ((*reconstruct_init)->intermediate == 1)
	{
		//if (writePulses (reconstruct_init, samprate, tstartRecord, recordNOTFILTERED, numPulses, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, pulseHeightsgsl, dtcObject))
		if (writePulses (reconstruct_init, samprate, tstartRecord, recordNOTFILTERED, numPulses, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, dtcObject))
		{
			message = "Cannot run routine writePulses to write pulses in output file";
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
/*xxxx end of SECTION A8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A9 ************************************************************
* writePulses function:
*
******************************************************************************/
int writePulses(ReconstructInitSIRENA** reconstruct_init, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, fitsfile *dtcObject)
{
	int status = EPOK;
	string message = "";

	// Declare variables
	int t0;   					// First value of index of pulse
	gsl_matrix *vgslout2;

	gsl_vector_view temp;

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

    // If intermediate=1 => First record, createDetectFile
    //	                    Change clobber to 2
	//                   => Not first record, append info to output dtc file (because clobber is 2)
	long totalpulses;	// It is necessary to know the row to write info
	if ((*reconstruct_init)->clobber == 1)
	{
		totalpulses = 0;
		if ((*reconstruct_init)->crtLib == 1)	(*reconstruct_init)->clobber = 2;
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

		// Converting bins to time
		for (int i=0; i<numPulsesRecord; i++)
		{
			t0 = gsl_vector_get (tstart,i);
			gsl_vector_set(tstart,i,initialtime + (gsl_vector_get (tstart,	i) * (1/samprate)));
			gsl_vector_set(tend,i,initialtime + (gsl_vector_get (tend,	i) * (1/samprate)));

			if (invectorNOTFIL->size - t0 > (*reconstruct_init)->pulse_length)	//The invectorNOTFIL has more bins than sizePulse
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
		strcpy(obj.unit,"a.u.");
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
/*xxxx end of SECTION A9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A10 ************************************************************
* writeTestInfo function:
*
******************************************************************************/
int writeTestInfo(ReconstructInitSIRENA* reconstruct_init, gsl_vector *recordNOTFILTERED, gsl_vector *recordDERIVATIVE, double threshold, fitsfile *dtcObject)
{
	int status = EPOK;
	string message = "";

	long totalrecords;

	char dtcName[256];
	strncpy(dtcName,reconstruct_init->detectFile,255);
	dtcName[255]='\0'; // enforce zero ending string in case of buffer overflows

	// To work with tables (extensions)
	char extname[20];
	int extver = 0;

	strcpy(extname,"TESTINFO");
	if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(dtcObject,&totalrecords, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName) + " (TESTINFO HDU)";
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	IOData obj;

	// Creating ORIGINAL Column
	obj.inObject = dtcObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"TESTINFO");
	obj.iniRow = totalrecords+1;
	obj.endRow = totalrecords+1;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	//strcpy(obj.nameCol,"ORIGINAL");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit," ");

	gsl_matrix *matrixToWrite = gsl_matrix_alloc(1,recordDERIVATIVE->size);
	/*gsl_matrix_set_row(matrixToWrite,0,recordNOTFILTERED);
	if (writeFitsComplex(obj, matrixToWrite))
	{
		message = "Cannot run routine writeFitsComplex for recordNOTFILTERED (TESTINFO HDU)";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}*/

	// Creating FILDER Column
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
/*xxxx end of SECTION A10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION  A11 ************************************************************
* calculateTemplate function:
******************************************************************************/
int calculateTemplate (ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight, gsl_matrix **covariance, gsl_matrix **weight)
{
	string message = "";

	// Declare (and initialize) variables
	int totalPulses = pulsesAll->ndetpulses + pulsesInRecord->ndetpulses;
	gsl_vector *tstart = gsl_vector_alloc(totalPulses);			// Tstart column from the output dtc FITS file
	gsl_vector *pulseheight = gsl_vector_alloc(totalPulses);	// EstEnrgy column from the output dtc FITS file
	gsl_vector *quality = gsl_vector_alloc(totalPulses);		// Quality column from the output dtc FITS file
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

	gsl_vector *nonpileup = gsl_vector_alloc(totalPulses);	// Piled up pulse => Not taking into account to calculate the template
	long nonpileupPulses = totalPulses;						// A priori, all the found pulses are considered as non piled up
	gsl_vector_set_all(nonpileup,1);

	int nBins;										// Square-root choice (used by Excel and many others)
	gsl_vector *xhisto;								// X-axis of the pulseheights histogram
	gsl_vector *yhisto;								// Y-axis of the pulseheights histogram
	int index_maximumpulseheight;					// Index where the maximum of the pulseheights histogram is
	double maximumpulseheight;						// Maximum of the pulseheights histogram

	bool firstnonpileupPulse = true;
	gsl_vector *pulse = gsl_vector_alloc(reconstruct_init->pulse_length);

	double tstartnext;

	gsl_vector_view temp;

	// Only if one of the found pulses is validated as non piled up pulse by using EstEnrgy
	// (pulseheights histogram) and Tstart, and Quality => The pulse, I0, is going
	// to be read from the dtc output FITS file (in order to not handle a long matrix of
	// pulses, found pulses x sizePulse_b)

	gsl_vector_scale(tstart,samprate); 	//tstarts not in sec but in samples

	gsl_vector *pulseheightAUX = gsl_vector_alloc(totalPulses);
	int cnt = 0;
	for (int i=0;i<totalPulses;i++)
	{
		if (i == totalPulses-1)		tstartnext = gsl_vector_get(tstart,i)+2*reconstruct_init->pulse_length;
		else						tstartnext = gsl_vector_get(tstart,i+1);

		 if ((tstartnext-gsl_vector_get(tstart,i) > reconstruct_init->pulse_length) && (gsl_vector_get(quality,i) == 0))
		 {
			 gsl_vector_set(pulseheightAUX,cnt,gsl_vector_get(pulseheight,i));
			 cnt = cnt +1;
		 }
	}
	temp = gsl_vector_subvector(pulseheightAUX,0,cnt);
	gsl_vector *pulseheightAUX2 = gsl_vector_alloc(cnt);
	gsl_vector_memcpy(pulseheightAUX2,&temp.vector);
	gsl_vector_free(pulseheightAUX);

	// Create pulseheights histogram
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

	for (int i=0;i<totalPulses;i++)
	{
		if (i == totalPulses-1)		tstartnext = gsl_vector_get(tstart,i)+2*reconstruct_init->pulse_length;
		else 						tstartnext = gsl_vector_get(tstart,i+1);

		// Check if the pulse is piled up or not
		if ((gsl_vector_get(pulseheight,i) < maximumpulseheight-0.1*maximumpulseheight) || (gsl_vector_get(pulseheight,i) > maximumpulseheight+0.1*maximumpulseheight)
			|| (tstartnext-gsl_vector_get(tstart,i) <= reconstruct_init->pulse_length) || (gsl_vector_get(quality,i) >= 1))
		{
 			gsl_vector_set(nonpileup,i,0);
			nonpileupPulses --;
			//cout<<"Pile-up: "<<i<<endl;
		}
		else
		{
			//cout<<"NO pile-up: "<<i<<endl;
			if (i < pulsesAll->ndetpulses)
			{
				gsl_vector_memcpy(pulse,pulsesAll->pulses_detected[i].pulse_adc);
			}
			else
			{
				gsl_vector_memcpy(pulse,pulsesInRecord->pulses_detected[i-pulsesAll->ndetpulses].pulse_adc);
			}

			// Non piled up pulses => Read, align and average them
			if (firstnonpileupPulse == true)
			{
				gsl_vector_memcpy(*pulseaverage,pulse);
				*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
			}
			else
			{
				if (align(samprate, pulseaverage,&pulse))
				{
				    message = "Cannot run align for pulse " + boost::lexical_cast<std::string>(i) + " when 1st pulse is piled-up";
				    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_vector_add(*pulseaverage,pulse);
				*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
			}
			if (firstnonpileupPulse == true)	firstnonpileupPulse = false;
		}
	}

	gsl_vector_scale(*pulseaverage,1.0/(nonpileupPulses));

	//cout<<"Pulsos usados para promediar "<<nonpileupPulses<<" de un total de "<<totalPulses<<" detectados"<<endl;

	if ((strcmp(reconstruct_init->EnergyMethod,"WEIGHT") == 0) || (strcmp(reconstruct_init->EnergyMethod,"WN") == 0))
	{
		if (weightMatrix(reconstruct_init, pulsesAll, pulsesInRecord, nonpileupPulses, nonpileup, *pulseaverage, covariance, weight))
		{
			message = "Cannot run weightMatrix routine";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	//cout<<"COVAR1 "<<gsl_matrix_get(*covariance,0,1)<<" "<<gsl_matrix_get(*covariance,0,2)<<" "<<gsl_matrix_get(*covariance,0,3)<<endl;
	//cout<<"COVAR1 "<<gsl_matrix_get(*covariance,1,0)<<" "<<gsl_matrix_get(*covariance,2,0)<<" "<<gsl_matrix_get(*covariance,3,0)<<endl;

	//Just in case due to the noise influence in the alignment, the first sample of the pulseaverage is not around the baseline but higher
	double meanLast200points, sgLast200points;
	temp = gsl_vector_subvector(*pulseaverage,reconstruct_init->pulse_length-200-1,200);
	if (findMeanSigma (&temp.vector, &meanLast200points, &sgLast200points))
	{
		message = "Cannot run findMeanSigma routine for kappa-sigma iteration";
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	//cout<<"pulseaverage(0): "<<gsl_vector_get(*pulseaverage,0)<<endl;
	if (fabs(gsl_vector_get(*pulseaverage,0))>fabs(meanLast200points)+3*sgLast200points)
	{
		gsl_vector_set(*pulseaverage,0,meanLast200points);
	}
	//cout<<"meanLast200points: "<<meanLast200points<<endl;
	//cout<<"pulseaverage(1): "<<gsl_vector_get(*pulseaverage,1)<<endl;

	*pulseaverageHeight = *pulseaverageHeight/nonpileupPulses;

	// Free allocate of GSL vectors
	gsl_vector_free(tstart);
	gsl_vector_free(pulseheight);
	gsl_vector_free(quality);
	gsl_vector_free(nonpileup);
	gsl_vector_free(xhisto);
	gsl_vector_free(yhisto);
	gsl_vector_free(pulse);

	return (EPOK);
}
/*xxxx end of SECTION A11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A12 ************************************************************
* createHisto function: This function builds ...
*
* - Declare variables
* - It is only going to work with the positive elements of the input vector -> invectoraux2
* - Check if all the values of the invector are the same => Histogrm of only one bin
* - Obtain invector_max and invector_min
* - Obtain binSize
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
    double invectormax= 0;      						// Maximum of invector
    double invectormin=1e10;  							// Minimum of invector
    double binSize;										// Size in samples of each bin
    int ind = 0;                						// Index of the bin which contains each invector element
    gsl_vector *invectoraux = gsl_vector_alloc(size);	// Auxiliary variable
    gsl_vector *invectoraux2;							// Auxiliary variable
    gsl_vector_view temp;								// In order to handle with gsl_vector_view (subvectors)

    // It is only going to work with the positive elements of the input vector -> invectoraux2
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

    // To check if all the values of the invector are the same
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
    	// Obtain invector_max
    	for (int i=0; i<size; i++)
    	{
    		if (invectormax < gsl_vector_get (invectoraux2,i))	invectormax = gsl_vector_get (invectoraux2,i);
    	}
    	// Obtain invector_min
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
/*xxxx end of SECTION A12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A13 ************************************************************
* align function:
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

	// FFT of vector1
	if (FFT(*vector1,vector1fft,SelectedTimeDuration))
	{
	    message = "Cannot run FFT routine for vector1";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	// FFT of vector 2
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

	// shiftdouble into shiftint
	if ((shiftdouble > -1) && (shiftdouble < 1)) shiftint = 0;
	else if (shiftdouble > 1)	shiftint = floor(shiftdouble);
	else if (shiftdouble < -1)	shiftint = ceil(shiftdouble);

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
/*xxxx end of SECTION A13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A14 ************************************************************
* shiftm function: This function returns (in vectorout) the vectorin delayed m samples.
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
/*xxxx end of SECTION A14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A15 ************************************************************
* shift_m function: This function returns (in vectorout) the vectorin moved forward m samples.
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
/*xxxx end of SECTION A15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A16 ************************************************************
* weightMatrix function:
*
******************************************************************************/
int weightMatrix (ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, long nonpileupPulses, gsl_vector *nonpileup, gsl_vector *pulseaverage, gsl_matrix **covariance, gsl_matrix **weight)
{
	double elementValue = 0.0;
	gsl_permutation *perm = gsl_permutation_alloc(reconstruct_init->pulse_length);
	int s=0;

	// Elements of the diagonal of the matrix
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
						(gsl_vector_get(pulsesAll->pulses_detected[p].pulse_adc,j)-gsl_vector_get(pulseaverage,j));	//KeV
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

	gsl_matrix *covarianceaux = gsl_matrix_alloc((*covariance)->size1,(*covariance)->size2);
	gsl_matrix_memcpy(covarianceaux,*covariance);
	gsl_linalg_LU_decomp(covarianceaux, perm, &s);
	gsl_linalg_LU_invert(covarianceaux, perm, *weight);
	gsl_matrix_free(covarianceaux);

	gsl_permutation_free(perm);

	return (EPOK);
}
/*xxxx end of SECTION A16 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A17 ************************************************************
* writeLibrary function:
*
******************************************************************************/
int writeLibrary(ReconstructInitSIRENA *reconstruct_init, double estenergy, gsl_vector **pulsetemplate, gsl_matrix *covariance, gsl_matrix *weight, bool appendToLibrary, fitsfile **inLibObject)
{
	int status = EPOK;
	int extver=0;
	string message = "";

	char inLibName[256];
	strncpy(inLibName, reconstruct_init->library_file,255);
	inLibName[255]='\0';

	// Declare variables
	gsl_vector *energyoutgsl = gsl_vector_alloc(1);
	gsl_vector *estenergyoutgsl = gsl_vector_alloc(1);
	gsl_matrix *pulsetemplates_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	gsl_matrix *pulsetemplatesb0_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	gsl_matrix *matchedfilters_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	gsl_matrix *matchedfiltersb0_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	gsl_matrix *weight_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *covariance_matrix = gsl_matrix_alloc(1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);

	char keyname[10];
	char *comment=NULL;
	char extname[20];
	IOData obj;

    if (appendToLibrary == true)
    {
    	long eventcntLib;
    	strcpy(keyname,"EVENTCNT");
    	if (fits_read_key(*inLibObject,TLONG,keyname, &eventcntLib,comment,&status))
    	{
    	    message = "Cannot read keyword " + string(keyname);
    	    EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}
    	if (eventcntLib <= 0)
    	{
    		message = "Legal values for read EVENTCNT (LIBRARY) are integer numbers greater than 0";
    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
    	}
    	long eventcntLib1 = eventcntLib + 1;
    	strcpy(keyname,"EVENTCNT");
	
    	if (fits_update_key(*inLibObject,TLONG,keyname, &eventcntLib1,comment,&status))
    	{
    	    message = "Cannot update keyword " + string(keyname);
    	    EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	gsl_vector *energycolumn = gsl_vector_alloc(eventcntLib+1);
    	gsl_vector *estenergycolumn = gsl_vector_alloc(eventcntLib+1);

    	gsl_vector *modelsrow = gsl_vector_alloc(reconstruct_init->pulse_length);
    	gsl_vector *modelsrowb0 = gsl_vector_alloc(reconstruct_init->pulse_length);
    	gsl_matrix *modelsaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
    	gsl_matrix *modelsb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
    	gsl_matrix *modelsaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
    	gsl_matrix *modelsb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);

    	gsl_vector *matchedfiltersrow = gsl_vector_alloc(reconstruct_init->pulse_length);
    	gsl_vector *matchedfiltersrowb0 = gsl_vector_alloc(reconstruct_init->pulse_length);
    	gsl_matrix *matchedfiltersaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
    	gsl_matrix *matchedfiltersb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
    	gsl_matrix *matchedfiltersaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);
    	gsl_matrix *matchedfiltersb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length);

    	gsl_vector *weightrow = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
    	gsl_matrix *weightaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
    	gsl_matrix *weightaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);

    	gsl_vector *covariancerow = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
    	gsl_matrix *covarianceaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
    	gsl_matrix *covarianceaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);

    	for (int i=0;i<eventcntLib;i++)
    	{
    		gsl_vector_set(energycolumn,i,gsl_vector_get(reconstruct_init->library_collection->energies,i));
            gsl_vector_set(estenergycolumn,i,gsl_vector_get(reconstruct_init->library_collection->pulse_heights,i));
		
    		gsl_matrix_set_row(modelsaux,i,reconstruct_init->library_collection->pulse_templates[i].ptemplate);
    		gsl_matrix_set_row(modelsb0aux,i,reconstruct_init->library_collection->pulse_templates_B0[i].ptemplate);
    		gsl_matrix_set_row(matchedfiltersaux,i,reconstruct_init->library_collection->matched_filters[i].mfilter);
    		gsl_matrix_set_row(matchedfiltersb0aux,i,reconstruct_init->library_collection->matched_filters_B0[i].mfilter);

    		if ((strcmp(reconstruct_init->EnergyMethod,"WEIGHT") == 0) || (strcmp(reconstruct_init->EnergyMethod,"WN") == 0))
    		{
    			gsl_matrix_get_row(weightrow,reconstruct_init->library_collection->W,i);
    			gsl_matrix_set_row(weightaux,i,weightrow);

    			gsl_matrix_get_row(covariancerow,reconstruct_init->library_collection->V,i);
    			gsl_matrix_set_row(covarianceaux,i,covariancerow);
    		}
    	}

    	gsl_vector_set(energycolumn,eventcntLib,reconstruct_init->monoenergy);
    	gsl_vector_set(estenergycolumn,eventcntLib,estenergy);
    	gsl_matrix_set_row(modelsaux,eventcntLib,*pulsetemplate);

    	gsl_vector *row_aux= gsl_vector_alloc((*pulsetemplate)->size);
    	gsl_vector_memcpy(row_aux,*pulsetemplate);
    	gsl_vector_scale(row_aux,1/reconstruct_init->monoenergy);
    	gsl_matrix_set_row(matchedfiltersaux,eventcntLib,row_aux);

    	gsl_vector *baselinegsl = gsl_vector_alloc((*pulsetemplate)->size);
    	gsl_vector_set_all(baselinegsl,-1.0*reconstruct_init->baseline);
    	gsl_vector_add(*pulsetemplate,baselinegsl);
    	gsl_vector_free(baselinegsl);
    	gsl_matrix_set_row(pulsetemplatesb0_matrix,0,*pulsetemplate);
    	gsl_matrix_set_row(modelsb0aux,eventcntLib,*pulsetemplate);

    	gsl_vector_memcpy(row_aux,*pulsetemplate);
    	gsl_vector_scale(row_aux,1/reconstruct_init->monoenergy);
    	gsl_matrix_set_row(matchedfiltersb0aux,eventcntLib,row_aux);
    	gsl_vector_free(row_aux);

    	if ((strcmp(reconstruct_init->EnergyMethod,"WEIGHT") == 0) || (strcmp(reconstruct_init->EnergyMethod,"WN") == 0))
    	{
    		matrix2vector(weight,&weightrow);
    		gsl_matrix_set_row(weightaux,eventcntLib,weightrow);

    		matrix2vector(covariance,&covariancerow);
    		gsl_matrix_set_row(covarianceaux,eventcntLib,covariancerow);
    	}

    	gsl_permutation *perm = gsl_permutation_alloc(eventcntLib+1);
    	// 'gsl_sort_vector_index' indirectly sorts the elements of the vector v into ascending order, storing the resulting
    	// permutation in p. The elements of p give the index of the vector element which would have been stored in that position
    	// if the vector had been sorted in place. The first element of p gives the index of the least element in v, and the last
    	// element of p gives the index of the greatest element in v. The vector v is not changed.
    	// Example: tstartaux=(5200 6000 200 3000) tauxsorted=(200 3000 5200 6000) perm=(2 3 0 1)
    	gsl_sort_vector_index(perm,energycolumn);
    	gsl_vector *energycolumnaux = gsl_vector_alloc(eventcntLib+1);
    	gsl_vector *estenergycolumnaux = gsl_vector_alloc(eventcntLib+1);
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

    		if ((strcmp(reconstruct_init->EnergyMethod,"WEIGHT") == 0) || (strcmp(reconstruct_init->EnergyMethod,"WN") == 0))
    		{
    			gsl_matrix_get_row(weightrow,weightaux,gsl_permutation_get(perm,i));
    			gsl_matrix_set_row(weightaux1,i,weightrow);

    			gsl_matrix_get_row(covariancerow,covarianceaux,gsl_permutation_get(perm,i));
    			gsl_matrix_set_row(covarianceaux1,i,covariancerow);
    		}
    	}

    	gsl_vector_memcpy(energycolumn,energycolumnaux);
    	gsl_vector_memcpy(estenergycolumn,estenergycolumnaux);
    	gsl_matrix_memcpy(modelsaux,modelsaux1);

    	gsl_matrix_memcpy(modelsb0aux,modelsb0aux1);
    	gsl_matrix_memcpy(matchedfiltersaux,matchedfiltersaux1);
    	gsl_matrix_memcpy(matchedfiltersb0aux,matchedfiltersb0aux1);

    	if ((strcmp(reconstruct_init->EnergyMethod,"WEIGHT") == 0) || (strcmp(reconstruct_init->EnergyMethod,"WN") == 0))
    	{
    		gsl_matrix_memcpy(weightaux,weightaux1);

    		gsl_matrix_memcpy(covarianceaux,covarianceaux1);
    	}

    	gsl_permutation_free(perm);
    	gsl_vector_free(energycolumnaux);
    	gsl_vector_free(estenergycolumnaux);
    	gsl_matrix_free(modelsaux1);
    	gsl_matrix_free(modelsb0aux1);
    	gsl_matrix_free(matchedfiltersaux1);
    	gsl_matrix_free(matchedfiltersb0aux1);
    	gsl_matrix_free(weightaux1);
    	gsl_matrix_free(covarianceaux1);

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
    	    gsl_vector_set (energyoutgsl,0,gsl_vector_get(energycolumn,i));
    	    if (writeFitsSimple(obj, energyoutgsl))
    	    {
    	    	message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
    	    	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	    }

    	    strcpy(obj.nameCol,"PHEIGHT");
    	    strcpy(obj.unit,"a.u.");
    	    gsl_vector_set (estenergyoutgsl,0,gsl_vector_get(estenergycolumn,i));
    	    if (writeFitsSimple(obj, estenergyoutgsl))
    	    {
    	    	message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
    	    	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	    }

    	    strcpy(obj.nameCol,"PULSE");
    	    strcpy(obj.unit,"a.u.");
    	    gsl_matrix_get_row(modelsrow,modelsaux,i);
    	    gsl_matrix_set_row(pulsetemplates_matrix,0,modelsrow);
    	    if (writeFitsComplex(obj, pulsetemplates_matrix))
    	    {
    	    	message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    	    	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	    }

    	    strcpy(obj.nameCol,"PULSEB0");
    	    strcpy(obj.unit,"a.u.");
    	    gsl_matrix_get_row(modelsrowb0,modelsb0aux,i);
    	    gsl_matrix_set_row(pulsetemplatesb0_matrix,0,modelsrowb0);
    	    if (writeFitsComplex(obj, pulsetemplatesb0_matrix))
    	    {
    	      	message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    	      	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	    }

    	    strcpy(obj.nameCol,"MF");
       	    strcpy(obj.unit," ");
       	    gsl_matrix_get_row(matchedfiltersrow,matchedfiltersaux,i);
       	    //gsl_vector_scale(matchedfiltersrow,1./reconstruct_init->monoenergy);
       	    gsl_matrix_set_row(matchedfilters_matrix,0,matchedfiltersrow);
       	    if (writeFitsComplex(obj, matchedfilters_matrix))
       	    {
      	    	message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
       	    	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
       	    }

       	    strcpy(obj.nameCol,"MFB0");
       	    strcpy(obj.unit," ");
       	    gsl_matrix_get_row(matchedfiltersrowb0,matchedfiltersb0aux,i);
       	    gsl_matrix_set_row(matchedfiltersb0_matrix,0,matchedfiltersrowb0);
       	    if (writeFitsComplex(obj, matchedfiltersb0_matrix))
       	    {
       	      	message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
       	      	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
       	    }

       	    if ((strcmp(reconstruct_init->EnergyMethod,"WEIGHT") == 0) || (strcmp(reconstruct_init->EnergyMethod,"WN") == 0))
       	    {
       	    	strcpy(obj.nameCol,"COVARM");
       	    	strcpy(obj.unit," ");
       	    	gsl_matrix_get_row(covariancerow,covarianceaux,i);
       	    	gsl_matrix_set_row(covariance_matrix,0,covariancerow);
       	    	if (writeFitsComplex(obj, covariance_matrix))
       	    	{
       	    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
       	    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
       	    	}

       	    	strcpy(obj.nameCol,"WEIGHTM");
       	    	strcpy(obj.unit," ");
       	    	gsl_matrix_get_row(weightrow,weightaux,i);
       	    	gsl_matrix_set_row(weight_matrix,0,weightrow);
       	    	if (writeFitsComplex(obj, weight_matrix))
       	    	{
       	    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
       	    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
       	    	}
       	    }
    	}

    	gsl_vector_free(modelsrow);
    	gsl_vector_free(modelsrowb0);
    	gsl_vector_free(matchedfiltersrow);
    	gsl_vector_free(matchedfiltersrowb0);

    	gsl_vector_free(weightrow);
    	gsl_vector_free(covariancerow);

    	gsl_vector_free (energycolumn);
    	gsl_vector_free (estenergycolumn);
    	gsl_matrix_free(modelsaux);
    	gsl_matrix_free(modelsb0aux);
    	gsl_matrix_free(matchedfiltersaux);
    	gsl_matrix_free(matchedfiltersb0aux);

    	gsl_matrix_free(weightaux);
    	gsl_matrix_free(covarianceaux);

    	delete [] obj.nameTable;
    	delete [] obj.nameCol;
    	delete [] obj.unit;

    	// Primary HDU
    	strcpy(extname,"Primary");
    	int *hdutype;
    	if (fits_movabs_hdu(*inLibObject, 1, hdutype, &status))
    	{
    		message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char keyvalstr[1000];

    	char str_procnumber[125];			sprintf(str_procnumber,"%d",eventcntLib);
    	string strprocname (string("PROC") + string(str_procnumber));
    	strcpy(keyname,strprocname.c_str());
    	string strprocval (string("PROC") + string(str_procnumber) + string(" Starting parameter list"));
    	strcpy(keyvalstr,strprocval.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	string strproc (string("RecordFile = ") + reconstruct_init->record_file);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("TesEventFile = ") + reconstruct_init->event_file;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("LibraryFile = ") + reconstruct_init->library_file;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("NoiseFile = ") + reconstruct_init->noise_file;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_mode[125];			sprintf(str_mode,"%d",reconstruct_init->mode);
    	strproc = string("mode = ") + string(str_mode);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_crtLib[125];		sprintf(str_crtLib,"%d",reconstruct_init->crtLib);
    	strproc=string("crtLib = ") + string(str_crtLib);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_lastELibrary[125];		sprintf(str_lastELibrary,"%d",reconstruct_init->lastELibrary);
    	strproc=string("lastELibrary = ") + string(str_lastELibrary);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("PixelType = ") + reconstruct_init->PixelType;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("FilterDomain = ") + reconstruct_init->FilterDomain;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("FilterMethod = ") + reconstruct_init->FilterMethod;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_calibLQ[125];      sprintf(str_calibLQ,"%d",reconstruct_init->calibLQ);
    	strproc=string("calibLQ = ") + string(str_calibLQ);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_b_cF[125];	sprintf(str_b_cF,"%f",reconstruct_init->b_cF);
    	strproc=string("b_cF = ") + string(str_b_cF);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_c_cF[125];	sprintf(str_c_cF,"%f",reconstruct_init->c_cF);
    	strproc=string("c_cF = ") + string(str_c_cF);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_maxPulsesPerRecord[125];	sprintf(str_maxPulsesPerRecord,"%d",reconstruct_init->maxPulsesPerRecord);
    	strproc=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_pulse_length[125];	sprintf(str_pulse_length,"%d",reconstruct_init->pulse_length);
    	strproc=string("PulseLength = ") + string(str_pulse_length);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_tauFall[125];		sprintf(str_tauFall,"%e",reconstruct_init->tauFall);
    	strproc=string("tauFall = ") + string(str_tauFall);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",reconstruct_init->scaleFactor);
    	strproc=string("scaleFactor = ") + string(str_scaleFactor);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_samplesUp[125];	sprintf(str_samplesUp,"%f",reconstruct_init->samplesUp);
    	strproc=string("samplesUp = ") + string(str_samplesUp);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_nSgms[125];	    sprintf(str_nSgms,"%f",reconstruct_init->nSgms);
    	strproc=string("nSgms = ") + string(str_nSgms);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_baseline[125];	sprintf(str_baseline,"%f",reconstruct_init->baseline);
    	strproc=string("baseline = ") + string(str_baseline);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_LrsT[125];			sprintf(str_LrsT,"%e",reconstruct_init->LrsT);
    	strproc=string("LrsT = ") + string(str_LrsT);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_LbT[125];			sprintf(str_LbT,"%e",reconstruct_init->LbT);
    	strproc=string("LbT = ") + string(str_LbT);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_monoenergy[125];	sprintf(str_monoenergy,"%f",reconstruct_init->monoenergy);
    	strproc=string("monoenergy = ") + string(str_monoenergy);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_intermediate[125];      sprintf(str_intermediate,"%d",reconstruct_init->intermediate);
    	strproc=string("intermediate = ") + string(str_intermediate);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("detectFile = ") + reconstruct_init->detectFile;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("RecordFileCalib2 = ") + reconstruct_init->record_file2;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_monoenergy2[125];	sprintf(str_monoenergy2,"%f",reconstruct_init->monoenergy2);
    	strproc=string("monoenergy2 = ") + string(str_monoenergy2);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strproc=string("filterFile = ") + reconstruct_init->filterFile;
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_clobber[125];      sprintf(str_clobber,"%d",reconstruct_init->clobber);
    	strproc=string("clobber = ") + string(str_clobber);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_tstartPulse1[125];	sprintf(str_tstartPulse1,"%d",reconstruct_init->tstartPulse1);
    	strproc=string("tstartPulse1 = ") + string(str_tstartPulse1);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_tstartPulse2[125];	sprintf(str_tstartPulse2,"%d",reconstruct_init->tstartPulse2);
    	strproc=string("tstartPulse2 = ") + string(str_tstartPulse2);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	char str_tstartPulse3[125];	sprintf(str_tstartPulse3,"%d",reconstruct_init->tstartPulse3);
    	strproc=string("tstartPulse3 = ") + string(str_tstartPulse3);
    	strcpy(keyvalstr,strproc.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	strcpy(keyname,strprocname.c_str());
    	strprocval = string("PROC") + string(str_procnumber) + string(" Ending parameter list");
    	strcpy(keyvalstr,strprocval.c_str());
    	if (fits_write_key(*inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
    	{
    		message = "Cannot write keyword " + string(keyname) + " in library file " + string(inLibName);
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}
    }
    else
    {
    	strcpy(extname,"LIBRARY");
    	if (fits_movnam_hdu(*inLibObject, ANY_HDU,extname, extver, &status))
    	{
    		message = "Cannot move to HDU  " + string(extname) + " in library";
    		EP_PRINT_ERROR(message,status);return(EPFAIL);
    	}

    	if (reconstruct_init->pulse_length <= 0)
    	{
    		message = "Legal values for EVENTSZ (PULSES) are integer numbers greater than 0";
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}
    	if (fits_write_key(*inLibObject,TINT,"EVENTSZ",&reconstruct_init->pulse_length,comment,&status))
    	{
    		message = "Cannot write keyword EVENTSZ in library";
    		EP_PRINT_ERROR(message,status);return(EPFAIL);
    	}

    	gsl_vector_set (energyoutgsl,0,reconstruct_init->monoenergy);
    	gsl_vector_set (estenergyoutgsl,0,estenergy);

    	gsl_matrix_set_row(pulsetemplates_matrix,0,*pulsetemplate);

    	gsl_vector *baselinegsl = gsl_vector_alloc((*pulsetemplate)->size);
    	gsl_vector_set_all(baselinegsl,-1.0*reconstruct_init->baseline);
    	gsl_vector_add(*pulsetemplate,baselinegsl);
    	gsl_vector_free(baselinegsl);
    	gsl_matrix_set_row(pulsetemplatesb0_matrix,0,*pulsetemplate);

    	// Creating ENERGY Column
    	obj.inObject = *inLibObject;
    	obj.nameTable = new char [255];
    	strcpy(obj.nameTable,"LIBRARY");
    	obj.iniRow = 1;
    	obj.endRow = 1;
    	obj.iniCol = 0;
    	obj.nameCol = new char [255];
    	strcpy(obj.nameCol,"ENERGY");
    	obj.type = TDOUBLE;
    	obj.unit = new char [255];
    	strcpy(obj.unit,"eV");
    	if (writeFitsSimple(obj, energyoutgsl))
    	{
    		message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating PHEIGHT Column
    	strcpy(obj.nameCol,"PHEIGHT");
    	strcpy(obj.unit,"a.u.");
    	if (writeFitsSimple(obj, estenergyoutgsl))
    	{
    		message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating PULSE Column
    	strcpy(obj.nameCol,"PULSE");
    	strcpy(obj.unit,"a.u.");
    	if (writeFitsComplex(obj, pulsetemplates_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating PULSEB0 Column
    	strcpy(obj.nameCol,"PULSEB0");
    	strcpy(obj.unit,"a.u.");
    	if (writeFitsComplex(obj, pulsetemplatesb0_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating MF Column
    	strcpy(obj.nameCol,"MF");
    	strcpy(obj.unit," ");
    	gsl_matrix_scale(pulsetemplates_matrix,1.0/reconstruct_init->monoenergy);
    	if (writeFitsComplex(obj, pulsetemplates_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating MFB0 Column
    	strcpy(obj.nameCol,"MFB0");
    	strcpy(obj.unit," ");
    	gsl_matrix_scale(pulsetemplatesb0_matrix,1.0/reconstruct_init->monoenergy);
    	if (writeFitsComplex(obj, pulsetemplatesb0_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	if ((strcmp(reconstruct_init->EnergyMethod,"WEIGHT") == 0) || (strcmp(reconstruct_init->EnergyMethod,"WN") == 0))
    	{
    		// Creating COVARM Column
    		strcpy(obj.nameCol,"COVARM");
    		strcpy(obj.unit," ");
    		if (writeFitsComplex(obj, covariance))
    		{
    			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    		}

    		// Creating WEIGHTM Column
    		strcpy(obj.nameCol,"WEIGHTM");
    		strcpy(obj.unit," ");
    		if (writeFitsComplex(obj, weight))
    		{
    			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    		}
    	}

    	delete [] obj.nameTable;
    	delete [] obj.nameCol;
    	delete [] obj.unit;
    }

    if (fits_close_file(*inLibObject,&status))
    {
    	message = "Cannot close file " + string(inLibName);
    	EP_PRINT_ERROR(message,status);return(EPFAIL);
    }

    // Free allocate of GSL vectors
    gsl_vector_free(energyoutgsl);
    gsl_vector_free(estenergyoutgsl);

    gsl_matrix_free (pulsetemplates_matrix);
    gsl_matrix_free (pulsetemplatesb0_matrix);
    gsl_matrix_free (matchedfilters_matrix);
    gsl_matrix_free (matchedfiltersb0_matrix);

    return (EPOK);
}
/*xxxx end of SECTION A17 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A18 ************************************************************
* matrix2vector function: This function transforms an input n x n square matrix into an output n^2 vector.
*                         It puts the first row of the matrix (n elements) in the first n elements of the vector (from 0 to n-1),
*                         the second row of the matrix in the (n,2n-1)elements of the vector and so on.
******************************************************************************/
int matrix2vector (gsl_matrix *matrixin, gsl_vector **vectorout)
{
	// matrixin is a square matrix
	int dim = matrixin->size1;
	gsl_vector *vectoraux = gsl_vector_alloc(dim);

	//*vectorout = gsl_vector_alloc(dim*dim);

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
* vector2matrix function: This function transforms an input n^2 vector into an output n x n square matrix.
*                         It puts the first n elements of the vector in the first row of the matrix,
*                         the second group of n elements (from n to 2n-1) of the vector in the second row and so on.
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
* fillInLibraryData: This function ...
*
******************************************************************************/
int fillInLibraryData (ReconstructInitSIRENA* reconstruct_init)
{
	string message="";

	char extname[20];
	//char keyname[10];
	//char *comment=NULL;
	int extver=0;
	int status = EPOK;

	fitsfile *inLibObject = NULL;

	char inLibName[256];
	strncpy(inLibName, reconstruct_init->library_file,255);
	inLibName[255]='\0';

	long totalrows;

	IOData obj;

	gsl_vector *vectorAux_E;
	gsl_matrix *matrixAux_PULSE;
	gsl_matrix *matrixAux_W;
	gsl_vector *matrixrow = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *matrixrow2 = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_matrix *Walpha = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
	//gsl_matrix *Wbeta = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
	gsl_vector *T = gsl_vector_alloc(reconstruct_init->pulse_length);
	double t;
	gsl_matrix *X = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
	gsl_vector *Xvector = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	gsl_vector *Y = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *Z = gsl_vector_alloc(reconstruct_init->pulse_length);
	double r;
	double Ea, Eb;
	gsl_vector *Pab = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *Dab = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_vector *vector_aux = gsl_vector_alloc(1);
	gsl_vector *vector_aux2 = gsl_vector_alloc(reconstruct_init->pulse_length);
	gsl_matrix *matrix_aux = gsl_matrix_alloc(1,reconstruct_init->pulse_length);
	gsl_matrix *matrix_aux2 = gsl_matrix_alloc(1,reconstruct_init->pulse_length*reconstruct_init->pulse_length);

	// Open the library FITS file
	if (fits_open_file(&inLibObject,inLibName,READWRITE,&status))
	{
	    message = "Cannot open library file " + string(inLibName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(extname,"LIBRARY");
	if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(inLibObject,&totalrows, &status))
	{
		message = "Cannot get number of rows in " + string(inLibName);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	obj.inObject = inLibObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"LIBRARY");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.unit = new char [255];
	obj.type = TDOUBLE;

	obj.iniRow = 1;
	obj.endRow = totalrows;
	strcpy(obj.nameCol,"ENERGY");
	vectorAux_E = gsl_vector_alloc(totalrows);
	if (readFitsSimple (obj,&vectorAux_E))
	{
		message = "Cannot run readFitsSimple routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	strcpy(obj.nameCol,"PULSE");
	matrixAux_PULSE = gsl_matrix_alloc(totalrows,reconstruct_init->pulse_length);
	if (readFitsComplex (obj,&matrixAux_PULSE))
	{
		message = "Cannot run readFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	strcpy(obj.nameCol,"WEIGHTM");
	matrixAux_W = gsl_matrix_alloc(totalrows,reconstruct_init->pulse_length*reconstruct_init->pulse_length);
	if (readFitsComplex (obj,&matrixAux_W))
	{
		message = "Cannot run readFitsComplex routine for column " + string(obj.nameCol);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	for (int i=0;i<totalrows-1;i++)
	{
		obj.iniRow = i+1;
		obj.endRow = i+1;

		strcpy(obj.nameCol,"TV");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(matrixrow,matrixAux_PULSE,i+1);
		gsl_vector_memcpy(T,matrixrow);
		gsl_matrix_get_row(matrixrow,matrixAux_PULSE,i);
		gsl_vector_sub(T,matrixrow);
		gsl_matrix_set_row(matrix_aux,0,T);
		if (writeFitsComplex(obj, matrix_aux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"tE");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(matrixrow2,matrixAux_W,i);
		vector2matrix(matrixrow2,&Walpha);
		gsl_blas_dgemv(CblasNoTrans, 1.0, Walpha, T, 0.0, vector_aux2);
		gsl_blas_ddot(T,vector_aux2,&t);
		gsl_vector_set(vector_aux,0,t);
		if (writeFitsSimple(obj, vector_aux))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"XM");
		strcpy(obj.unit," ");
		gsl_matrix_get_row(matrixrow2,matrixAux_W,i+1);
		gsl_vector_memcpy(Xvector,matrixrow2);
		gsl_matrix_get_row(matrixrow2,matrixAux_W,i);
		gsl_vector_sub(Xvector,matrixrow2);
		gsl_vector_scale(Xvector,1/t);
		gsl_matrix_set_row(matrix_aux2,0,Xvector);
		if (writeFitsComplex(obj, matrix_aux2))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"YV");
		strcpy(obj.unit," ");
		gsl_blas_dgemv(CblasNoTrans, 1.0, Walpha, T, 0.0, Y);
		gsl_vector_scale(Y,1/t);
		gsl_matrix_set_row(matrix_aux,0,Y);
		if (writeFitsComplex(obj, matrix_aux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"ZV");
		strcpy(obj.unit," ");
		vector2matrix(Xvector,&X);
		gsl_blas_dgemv(CblasNoTrans, 1.0, X, T, 0.0, Z);
		gsl_matrix_set_row(matrix_aux,0,Z);
		if (writeFitsComplex(obj, matrix_aux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"rE");
		strcpy(obj.unit," ");
		gsl_blas_ddot(Z,T,&r);
		r=1/r;
		gsl_vector_set(vector_aux,0,r);
		if (writeFitsSimple(obj, vector_aux))
		{
			message = "Cannot run writeFitsSimple routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"PAB");
		strcpy(obj.unit," ");
		Ea = gsl_vector_get(vectorAux_E,i);
		Eb = gsl_vector_get(vectorAux_E,i+1);
		gsl_vector_memcpy(Pab,T);
		gsl_vector_scale(Pab,-Ea/(Eb-Ea));
		gsl_matrix_get_row(matrixrow,matrixAux_PULSE,i);
		gsl_vector_add(Pab,matrixrow);
		gsl_matrix_set_row(matrix_aux,0,Pab);
		if (writeFitsComplex(obj, matrix_aux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		strcpy(obj.nameCol,"DAB");
		strcpy(obj.unit," ");
		gsl_vector_memcpy(Dab,T);
		gsl_vector_scale(Dab,1/(Eb-Ea));
		gsl_matrix_set_row(matrix_aux,0,Dab);
		if (writeFitsComplex(obj, matrix_aux))
		{
			message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	// Close the library FITS file
	if (fits_close_file(inLibObject,&status))
	{
		message = "Cannot close file " + string(inLibName);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	gsl_vector_free(vectorAux_E);
	gsl_matrix_free(matrixAux_PULSE);
	gsl_matrix_free(matrixAux_W);
	gsl_vector_free(matrixrow);
	gsl_vector_free(matrixrow2);
	gsl_matrix_free(Walpha);
	gsl_vector_free(T);
	gsl_matrix_free(X);
	gsl_vector_free(Xvector);
	gsl_vector_free(Y);
	gsl_vector_free(Z);
	gsl_vector_free(Pab);
	gsl_vector_free(Dab);
	gsl_vector_free(vector_aux);
	gsl_vector_free(vector_aux2);
	gsl_matrix_free(matrix_aux);
	gsl_matrix_free(matrix_aux2);

	return(EPOK);
}
/*xxxx end of SECTION A20 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B ************************************************************
* runFilter: This function ...
*
******************************************************************************/
void runFilter(TesRecord* record, int nRecord, int lastRecord, ReconstructInitSIRENA** reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter)
{
	const char * create= "runFilter v.15.0.0";	//Set "CREATOR" keyword of output FITS file

	string message="";
	int status = EPOK;

	fitsfile *dtcObject = NULL;	    // Object which contains information of the output FITS file
	if ((*reconstruct_init)->intermediate == 1)
	{
		char dtcName[256];
		strncpy(dtcName,(*reconstruct_init)->detectFile,255);
		dtcName[255]='\0';
	}

	int TorF;
	if (strcmp((*reconstruct_init)->FilterDomain,"T") == 0)
	{
	    TorF=0;
	}
	else if (strcmp((*reconstruct_init)->FilterDomain,"F") == 0)
	{
		TorF=1;
	}
	else
	{
	    message = "Parameter reconstruct_init->FilterDomain out of range: [T/F]";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	int runF0orB0val;
	if (strcmp((*reconstruct_init)->FilterMethod,"F0") == 0)
	{
		runF0orB0val = 0;
	}
	else if (strcmp((*reconstruct_init)->FilterMethod,"B0") == 0)
	{
		runF0orB0val = 1;
	}
	else
	{
	    message = "Parameter reconstruct_init->FilterMethod out of range: [F0/B0]";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	int runEMethod;
	if (strcmp((*reconstruct_init)->EnergyMethod,"NOLAGS") == 0)
	{
		runEMethod = 0;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"LAGS") == 0)
	{
		runEMethod = 1;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHT") == 0)
	{
		runEMethod = 2;
	}
	else if (strcmp((*reconstruct_init)->EnergyMethod,"WEIGHTN") == 0)
	{
		runEMethod = 3;
	}
	else
	{
	    message = "Parameter reconstruct_init->EnergyMethod out of range: [NOLAGS/LAGS/WEIGHT/WEIGHTN]";
	    EP_EXIT_ERROR(message,EPFAIL);
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
	else
	{
	    message = "Parameter reconstruct_init->OFStrategy out of range: [FREE/BASE2/BYGRADE/FIXED]";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	long energyInLibrary_row;

	double normalizationFactor;
	double uncE;

	// Declare variables
	gsl_vector *optimalfilter_SHORT = NULL;	// Resized optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector *optimalfilter_f_SHORT = NULL;		// Resized optimal filter f's when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector *optimalfilter_FFT_SHORT = NULL;	// Resized optimal filter magnitudes when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector_complex *optimalfilter_FFT_complex = NULL;

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

	// Store the record in 'invector'
	recordAux = gsl_vector_alloc(record->trigger_size);
	if (loadRecord(record, &tstartRecord, &recordAux))
	{
		message = "Cannot run routine loadRecord";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	//int eventsz = recordAux->size;	// Just in case the last record has been filled in with 0's => Re-allocate invector

	if (runF0orB0val == 1)
	{
	    //cout<<"Elimina la baseline"<<endl;
		gsl_vector *baselinegsl = gsl_vector_alloc(recordAux->size);
		gsl_vector_set_all(baselinegsl,-1.0*(*reconstruct_init)->baseline);
		gsl_vector_add(recordAux,baselinegsl);
		gsl_vector_free(baselinegsl);
	}

	if ((*reconstruct_init)->mode == 0)
	{
		filtergsl = gsl_vector_alloc((*reconstruct_init)->pulse_length);			// Filter values

		if ((strcmp((*reconstruct_init)-> EnergyMethod,"WEIGHT") != 0) && (strcmp((*reconstruct_init)-> EnergyMethod,"WEIGHTN") != 0))
		{
		if (nRecord == 1)
		{
			// Matched filter
			if (find_matchedfilter(runF0orB0val, (*reconstruct_init)->monoenergy, (*reconstruct_init)->library_collection->energies, (*reconstruct_init), &filtergsl))
			{
				message = "Cannot run routine find_matchedfilter for filter interpolation";
				EP_EXIT_ERROR(message,EPFAIL);
			}

			// Optimal filter
			if (calculus_optimalFilter (TorF, (*reconstruct_init)->intermediate, (*reconstruct_init)->mode, filtergsl, filtergsl->size, 1/record->delta_t, runF0orB0val, (*reconstruct_init)->noise_spectrum->noisefreqs, (*reconstruct_init)->noise_spectrum->noisespec, &optimalfilter_SHORT, &optimalfilter_f_SHORT, &optimalfilter_FFT_SHORT, &optimalfilter_FFT_complex, &normalizationFactor))
			{
				message = "Cannot run routine calculus_optimalFilter";
				EP_EXIT_ERROR(message,EPFAIL);
			}

			(*optimalFilter)->ofilter_duration = optimalfilter_SHORT->size;
			(*optimalFilter)->nrmfctr = normalizationFactor;
			(*optimalFilter)->ofilter = gsl_vector_alloc((*optimalFilter)->ofilter_duration);
			gsl_vector_memcpy((*optimalFilter)->ofilter,optimalfilter_SHORT);

			// Write output _flt FITS file
			if ((*reconstruct_init)->intermediate == 1)
			{
				if (writeFilter(*reconstruct_init, normalizationFactor, optimalfilter_SHORT, optimalfilter_f_SHORT, optimalfilter_FFT_SHORT, &dtcObject, create))
				{
					message = "Cannot run routine writeFilter";
					EP_EXIT_ERROR(message,EPFAIL);
				}
			}
		}
		else
		{
			optimalfilter_SHORT = gsl_vector_alloc((*optimalFilter)->ofilter_duration);
			optimalfilter_FFT_complex = gsl_vector_complex_alloc((*optimalFilter)->ofilter_duration);
			gsl_vector_complex_set_zero(optimalfilter_FFT_complex);
			normalizationFactor = (*optimalFilter)->nrmfctr;
			gsl_vector_memcpy(optimalfilter_SHORT,(*optimalFilter)->ofilter);
		}
		}

		// Check Quality: If there are no valid pulses in the DETECT FITS file => The task finishes
		iter = 0;
		for (int i=0; i<(*pulsesInRecord)->ndetpulses;i++)
		{
			if ((*pulsesInRecord)->pulses_detected[i].quality != 0)
			{
				iter++;
			}
		}
		if (iter == (*pulsesInRecord)->ndetpulses)
		{
			message = "There are no valid pulses (quality == 0) in one record";
			EP_PRINT_ERROR(message,-999);
		}

		for (int i=0; i<(*pulsesInRecord)->ndetpulses ;i++)
		{
			// Pulses are going to be validated by checking its quality
			if ((*pulsesInRecord)->pulses_detected[i].quality == 0)	// Neither truncated or saturated
			{
				// Pulse
				tstartSamplesRecord = floor((*pulsesInRecord)->pulses_detected[i].Tstart/record->delta_t+0.5)-tstartRecordSamples;
				pulse = gsl_vector_alloc((*pulsesInRecord)->pulses_detected[i].pulse_duration);
				temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,(*pulsesInRecord)->pulses_detected[i].pulse_duration);
				gsl_vector_memcpy(pulse,&temp.vector);

				// Calculate the uncalibrated energy of each pulse
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OJO q he puesto pulseGrade=0
				if (calculateUCEnergy(pulse,0,optimalfilter_SHORT,optimalfilter_FFT_complex,runEMethod,indexEalpha,indexEbeta,(*reconstruct_init),TorF,(*reconstruct_init)->mode,normalizationFactor,1/record->delta_t,&uncE))
				{
					message = "Cannot run calculateUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
					EP_EXIT_ERROR(message,EPFAIL);
				}
				gsl_vector_free(pulse);
				//cout<<"Pulso "<<i<<" "<<uncE<<endl;

				if ((*reconstruct_init)->intermediate == 1)
				{
					if (writeUCEnergy(reconstruct_init, pulsesAll, *pulsesInRecord, i, uncE, &dtcObject, create))
					{
						message = "Cannot run writeUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_EXIT_ERROR(message,EPFAIL);
					}
					//reconstruct_init->clobber = 2;
				}

				if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
				else								 (*pulsesInRecord)->pulses_detected[i].grade1 = optimalfilter_SHORT->size;
				//(*pulsesInRecord)->pulses_detected[i].grade1 = optimalfilter_SHORT->size;
				(*pulsesInRecord)->pulses_detected[i].ucenergy = uncE/1e3;	// In SIXTE, SIGNAL is in keV
				(*pulsesInRecord)->pulses_detected[i].energy = uncE/1e3;
			}
			else
			{
				if ((*reconstruct_init)->intermediate == 1)
				{
					if (writeUCEnergy(reconstruct_init, pulsesAll, *pulsesInRecord, i, uncE, &dtcObject, create))
					{
						message = "Cannot run writeUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_EXIT_ERROR(message,EPFAIL);
					}
					//reconstruct_init->clobber = 2;
				}

				if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
				else														(*pulsesInRecord)->pulses_detected[i].grade1 = 0;
				//(*pulsesInRecord)->pulses_detected[i].grade1 = 0;
				(*pulsesInRecord)->pulses_detected[i].ucenergy = 0.0;
				(*pulsesInRecord)->pulses_detected[i].energy = 0.0;
			}
		}

		gsl_vector_free(optimalfilter_SHORT);
		if (optimalfilter_f_SHORT != NULL) gsl_vector_free(optimalfilter_f_SHORT);
		if (optimalfilter_FFT_SHORT != NULL) gsl_vector_free(optimalfilter_FFT_SHORT);
		gsl_vector_complex_free(optimalfilter_FFT_complex);
		gsl_vector_free(filtergsl);
	}
	else if ((*reconstruct_init)->mode == 1)
	{
		long resize_mf;
		int pulseGrade; // HighRes=0, MidRes=1, LowRes=-1

		// Check Quality: If there are no valid pulses in the DETECT FITS file => The task finishes
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
			if ((*pulsesInRecord)->pulses_detected[i].quality < 10)	// Neither truncated or saturated
			{
				//resize_mf = (*pulsesInRecord)->pulses_detected[i].pulse_duration; // Resize the matched filter by using the tstarts
				//resize_mf = pow(2,floor(log2(resize_mf)));

				if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
				//else												 		(*pulsesInRecord)->pulses_detected[i].grade1 = resize_mf;
				else												 		(*pulsesInRecord)->pulses_detected[i].grade1 = (*pulsesInRecord)->pulses_detected[i].pulse_duration;

				if (pulseGrading(*reconstruct_init,(*pulsesInRecord)->pulses_detected[i].grade1,(*pulsesInRecord)->pulses_detected[i].grade2_1,OFlength_strategy,&pulseGrade,&resize_mf))
				{
					message = "Cannot run routine pulseGrading";
					EP_EXIT_ERROR(message,EPFAIL);
				}
				//cout<<"resize_mf="<<resize_mf<<endl;

				// Pulse
				tstartSamplesRecord = floor((*pulsesInRecord)->pulses_detected[i].Tstart/record->delta_t+0.5)-tstartRecordSamples;
				pulse = gsl_vector_alloc(resize_mf);
				temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf);
				gsl_vector_memcpy(pulse,&temp.vector);

				// EnergyMethod != WEIGHT/WEIGHTN => Filter and optimal filter
				if ((strcmp((*reconstruct_init)-> EnergyMethod,"WEIGHT") != 0) && (strcmp((*reconstruct_init)-> EnergyMethod,"WEIGHTN") != 0))
				{
					// Filter
					filtergsl_aux= gsl_vector_alloc((*reconstruct_init)->library_collection->matched_filters[0].mfilter_duration);
					if (find_matchedfilter(runF0orB0val, (*pulsesInRecord)->pulses_detected[i].maxDER, (*reconstruct_init)->library_collection->maxDERs, (*reconstruct_init), &filtergsl_aux))
					{
						message = "Cannot run routine find_matchedfilter for filter interpolation";
						EP_EXIT_ERROR(message,EPFAIL);
					}

					filtergsl = gsl_vector_alloc(resize_mf);			// Filter values
					temp = gsl_vector_subvector(filtergsl_aux,0,resize_mf);
					gsl_vector_memcpy(filtergsl,&temp.vector);
					gsl_vector_free(filtergsl_aux);

					// Calculate the optimal filter
					if (calculus_optimalFilter (TorF, (*reconstruct_init)->intermediate, (*reconstruct_init)->mode, filtergsl, filtergsl->size, 1/record->delta_t, runF0orB0val, (*reconstruct_init)->noise_spectrum->noisefreqs, (*reconstruct_init)->noise_spectrum->noisespec, &optimalfilter_SHORT, &optimalfilter_f_SHORT, &optimalfilter_FFT_SHORT, &optimalfilter_FFT_complex, &normalizationFactor))
					{
						message = "Cannot run routine calculus_optimalFilter";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				else
				{
					if (find_Esboundary((*pulsesInRecord)->pulses_detected[i].maxDER,(*reconstruct_init)->library_collection->maxDERs,(*reconstruct_init),&indexEalpha,&indexEbeta))
					{
						message = "Cannot run routine find_Esboundary for filter interpolation";
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}

				// Calculate the uncalibrated energy of each pulse
				if (calculateUCEnergy(pulse,pulseGrade,optimalfilter_SHORT,optimalfilter_FFT_complex,runEMethod,indexEalpha,indexEbeta,(*reconstruct_init),TorF,(*reconstruct_init)->mode,normalizationFactor,1/record->delta_t,&uncE))
				{
					message = "Cannot run calculateUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
					EP_EXIT_ERROR(message,EPFAIL);
				}

				// Subtract pulse model from the record
				if (find_model_energies(uncE, (*reconstruct_init), &model))
				{
				    message = "Cannot run find_model_energies routine for pulse i=" + boost::lexical_cast<std::string>(i);
				    EP_EXIT_ERROR(message,EPFAIL);
				}

				for (int j=tstartSamplesRecord;j<(tstartSamplesRecord+(model->size));j++)
				{
					gsl_vector_set(recordAux,j,gsl_vector_get(recordAux,j)-gsl_vector_get(model,j-tstartSamplesRecord));
				}

				if ((*reconstruct_init)->intermediate == 1)
				{
					if (writeFilterHDU(reconstruct_init, i, normalizationFactor, uncE, optimalfilter_SHORT, optimalfilter_f_SHORT, optimalfilter_FFT_SHORT, &dtcObject, create))
					{
						message = "Cannot run writeFilterHDU routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}

				(*pulsesInRecord)->pulses_detected[i].ucenergy = uncE/1e3; // In SIXTE, SIGNAL is in keV
				(*pulsesInRecord)->pulses_detected[i].energy = uncE/1e3;
			}
			else
			{
				optimalfilter_SHORT = gsl_vector_alloc((*reconstruct_init)->pulse_length);
				optimalfilter_f_SHORT = gsl_vector_alloc((*reconstruct_init)->pulse_length);
				optimalfilter_FFT_SHORT = gsl_vector_alloc((*reconstruct_init)->pulse_length);
				gsl_vector_set_zero(optimalfilter_SHORT);
				gsl_vector_set_zero(optimalfilter_f_SHORT);
				gsl_vector_set_zero(optimalfilter_FFT_SHORT);

				if ((*reconstruct_init)->intermediate == 1)
				{
					if (writeFilterHDU(reconstruct_init, i, 0.0, 0.0, optimalfilter_SHORT, optimalfilter_f_SHORT, optimalfilter_FFT_SHORT, &dtcObject, create))
					{
						message = "Cannot run writeFilterHDU routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_EXIT_ERROR(message,EPFAIL);
					}
				}

				if ((*pulsesInRecord)->pulses_detected[i].quality == 1)		(*pulsesInRecord)->pulses_detected[i].grade1 = -1;
				else														(*pulsesInRecord)->pulses_detected[i].grade1 = 0;
				(*pulsesInRecord)->pulses_detected[i].ucenergy = 0.0;
				(*pulsesInRecord)->pulses_detected[i].energy = 0.0;
			} // End if

			gsl_vector_free(optimalfilter_SHORT);
			gsl_vector_free(optimalfilter_f_SHORT);
			gsl_vector_free(optimalfilter_FFT_SHORT);
			if((*pulsesInRecord)->pulses_detected[i].quality < 10)
			{
				gsl_vector_free(pulse);
				gsl_vector_free(filtergsl);
				gsl_vector_complex_free(optimalfilter_FFT_complex);
			}
		} // End for

		gsl_vector_free(recordAux);
		gsl_vector_free(model);
	} // End if

	return;
}
/*xxxx end of SECTION B xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* find_energy: This function
*******************************************************************************/
int find_energy(double energyKeyword, gsl_vector *energygsl, long *rowFound)
{
	long nummodels = energygsl->size;
	*rowFound = -1;

	for (int i=0;i<nummodels;i++)
	{
		if (energyKeyword == gsl_vector_get(energygsl,i))
		{
			*rowFound = i;

			break;
		}
	}

    return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* calculus_optimalFilter:
****************************************/
int calculus_optimalFilter(int TorF, int intermediate, int mode, gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val, gsl_vector *freqgsl, gsl_vector *csdgsl, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, gsl_vector_complex **of_FFT_complex, double *normalizationFactor)
{
	// FFT calculus of the filter template (MATCHED FILTER->matchedfiltergsl)
	// Declare variables
	double SelectedTimeDuration = mf_size/samprate;
	string message = "";

	//Complex FFT values for positive and negative f's
	gsl_vector_complex *mfFFTcomp = gsl_vector_complex_alloc(mf_size);
	gsl_vector_complex *mfFFTcomp_conj = gsl_vector_complex_alloc(mf_size);
	gsl_vector *mf_arg = gsl_vector_alloc(mf_size);	// Argument for positive and negative f's
	gsl_vector *mf_f = gsl_vector_alloc(mf_size);				// Ordered f's according [-fmax,...,0,...,fmax]
	gsl_vector *mf_FFT = gsl_vector_alloc(mf_size);				// Ordered magnitude according [-fmax,...,0,...,fmax]

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
			long n0start = 0;
			long n0end = 0;
		   	if (interpolatePOS (freqgsl_POS, csdgsl_POS, mfPOS_size, gsl_vector_get(mf_f,1)-gsl_vector_get(mf_f,0), &n_f_interp, &n_FFT_interp, &n0start, &n0end))
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
	//gsl_vector_complex *of_FFT_complex = gsl_vector_complex_alloc(mf_size);
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

	//gsl_vector_complex_free(of_FFT_complex);
	gsl_vector_complex_free(mfFFTcomp_conj);
	gsl_vector_free(mf_f);
	gsl_vector_free(mf_FFT);
	gsl_vector_free(mf_arg);
	gsl_vector_free(n_f);
	gsl_vector_free(n_FFT);

    return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* interpolatePOS: This function interpolates an input vector, creating an output vector with the size and
*              frequency step given
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
* - numzerosstart: Number of zeros at the beginning of 'xout'
* - numzerosend: Number of zeros at the end of 'xout'
****************************************/
int interpolatePOS (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out, long *numzerosstart, long *numzerosend)
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

	*numzerosstart = 0;
	*numzerosend = 0;
	for (int i=0;i<size;i++)
	{
	    if ((i == 0) && (gsl_vector_get(*x_out,i) == 0))	    *numzerosstart = *numzerosstart+1;
	    else if ((i == 0) && (gsl_vector_get(*x_out,i) != 0)) 	break;
	    else if ((i != 0) && (gsl_vector_get(*x_out,i) == 0))	*numzerosstart = *numzerosstart+1;
	    else if ((i != 0) && (gsl_vector_get(*x_out,i) != 0))	break;
	}
	for (int i=(size-1);i>=0;i--)
	{
	    if ((i == size-1) && (gsl_vector_get(*x_out,i) == 0))	{*numzerosend = *numzerosend+1;}
	    else if ((i == size-1) && (gsl_vector_get(*x_out,i) != 0))	{break;}
	    else if ((i != size-1) && (gsl_vector_get(*x_out,i) == 0))	{*numzerosend = *numzerosend+1;}
	    else if ((i != size-1) && (gsl_vector_get(*x_out,i) != 0))	{break;}
	}

	return EPOK;
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* find_matchedfilter:
****************************************/
int find_matchedfilter(int runF0orB0val, double parameterToUse, gsl_vector *valuesToCompare, ReconstructInitSIRENA *reconstruct_init, gsl_vector **matchedfilterFound)
{
	string message = "";
	long nummodels = valuesToCompare->size;

	if (parameterToUse < gsl_vector_get(valuesToCompare,0))
	{
		if (runF0orB0val == 0)	gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[0].mfilter);
		else					gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters_B0[0].mfilter);

		//gsl_vector_scale(*matchedfilterFound,parameterToUse/gsl_vector_get(valuesToCompare,0));
	}
	else if (parameterToUse > gsl_vector_get(valuesToCompare,nummodels-1))
	{
		if (runF0orB0val == 0)	gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[nummodels-1].mfilter);
		else					gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters_B0[nummodels-1].mfilter);

		//gsl_vector_scale(*matchedfilterFound,parameterToUse/gsl_vector_get(valuesToCompare,nummodels-1));
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (parameterToUse == gsl_vector_get(valuesToCompare,i))
			{
				if (runF0orB0val == 0)	gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters[i].mfilter);
				else					gsl_vector_memcpy(*matchedfilterFound,reconstruct_init->library_collection->matched_filters_B0[i].mfilter);

				//gsl_vector_scale(*matchedfilterFound,parameterToUse/gsl_vector_get(valuesToCompare,i));

				break;
			}
			else if ((parameterToUse > gsl_vector_get(valuesToCompare,i)) && (parameterToUse < gsl_vector_get(valuesToCompare,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector *matchedfilterAux = gsl_vector_alloc(reconstruct_init->library_collection->matched_filters[0].mfilter_duration);
				gsl_vector_set_zero(matchedfilterAux);
				if (runF0orB0val == 0)
				{
					if (interpolate_model(&matchedfilterAux,parameterToUse,reconstruct_init->library_collection->matched_filters[i].mfilter,gsl_vector_get(valuesToCompare,i),
							reconstruct_init->library_collection->matched_filters[i+1].mfilter,gsl_vector_get(valuesToCompare,i+1)))
					{
						message = "Cannot run interpolate_model routine for model interpolation";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}
				}
				else
				{
					if (interpolate_model(&matchedfilterAux,parameterToUse,reconstruct_init->library_collection->matched_filters_B0[i].mfilter,gsl_vector_get(valuesToCompare,i),
							reconstruct_init->library_collection->matched_filters_B0[i+1].mfilter,gsl_vector_get(valuesToCompare,i+1)))
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
/*xxxx end of SECTION Bx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* find_Esboundary:
****************************************/
int find_Esboundary(double parameterToUse, gsl_vector *valuesToCompare, ReconstructInitSIRENA *reconstruct_init, int *indexEalpha, int *indexEbeta)
{
	string message = "";
	int nummodels = valuesToCompare->size;

	/*cout<<"parameterToUse: "<<parameterToUse<<endl;
	cout<<"valuesToCompare0: "<<gsl_vector_get(valuesToCompare,0)<<endl;
	cout<<"valuesToCompare1: "<<gsl_vector_get(valuesToCompare,1)<<endl;
	cout<<"valuesToCompare2: "<<gsl_vector_get(valuesToCompare,2)<<endl;*/

	if (parameterToUse < gsl_vector_get(valuesToCompare,0))
	{
		*indexEalpha = 0;
		*indexEbeta = 0;
	}
	else if (parameterToUse > gsl_vector_get(valuesToCompare,nummodels-1))
	{
		*indexEalpha = nummodels - 1;
		*indexEbeta = nummodels - 1;
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (parameterToUse == gsl_vector_get(valuesToCompare,i))
			{
				*indexEalpha = i;
				*indexEbeta = i;

				break;
			}
			else if ((parameterToUse > gsl_vector_get(valuesToCompare,i)) && (parameterToUse < gsl_vector_get(valuesToCompare,i+1)))
			{
				*indexEalpha = i;
				*indexEbeta = i+1;
			}
		}
	}

	return(EPOK);
}

/*xxxx end of SECTION Bx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* pulseGrading:
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

	//cout<<"grade1: "<<grade1<<endl;
	//cout<<"grade2_1: "<<grade2_1<<endl;
	if (grade2_1 < L2)
	{
		*pulseGrade = -1;
		if (OFlength_strategy == 0) 		*OFlength = grade1;
		else if (OFlength_strategy == 1) 	*OFlength = pow(2,floor(log2(grade1)));
		//else if (OFlength_strategy == 2) 	*OFlength = L2;
		else if (OFlength_strategy == 2) 	*OFlength = grade1;
	}
	else
	{
		if (grade1 >= H1)
		{
			*pulseGrade = 0;
			if (OFlength_strategy == 0)			*OFlength = grade1;
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
			if (OFlength_strategy == 0)			*OFlength = grade1;
			else if (OFlength_strategy == 1) 	*OFlength = pow(2,floor(log2(grade1)));
			//else if (OFlength_strategy == 2) 	*OFlength = L2;
			else if (OFlength_strategy == 2) 	*OFlength = grade1;
			else if (OFlength_strategy == 3) 	*OFlength = reconstruct_init->OFLength;
		}
	}

	//cout<<"pulseGrade: "<<*pulseGrade<<endl;
	//cout<<"OFlength: "<<*OFlength<<endl;

	return EPOK;
}
/*xxxx end of SECTION Bx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* calculateUCEnergy function:
****************************************************************************/
int calculateUCEnergy (gsl_vector *vector, int pulseGrade, gsl_vector *filter, gsl_vector_complex *filterFFT,int runEMethod, int indexEalpha, int indexEbeta, ReconstructInitSIRENA *reconstruct_init, int domain, int mode, double nrmfctr, double samprate, double *calculatedEnergy)
{
	string message = "";

	if ((runEMethod == 0) || (runEMethod == 1))
	{
		int numlags;
		gsl_vector *lags_vector;
		if ((runEMethod == 0) || (pulseGrade < 0))	//NOLAGS or LowRes
		//if (pulseGrade < 0) // LowRes
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

				if ((pulseGrade < 0) || (runEMethod == 0))		*calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
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

				if (mode == 0)
				{
					gsl_vector_complex *filterFFT_aux = gsl_vector_complex_alloc(filter->size);
					if (FFT(filter,filterFFT_aux,SelectedTimeDuration))
					{
						message = "Cannot run routine FFT when domain=1)";
						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
					}

					for (int j=0;j<numlags;j++)
					{
						for (int i=0; i<vector->size; i++)
						{
							argExp = -2*pi*gsl_vector_get(lags_vector,j)*i/vector->size;
							gsl_vector_complex_set(calculatedEnergy_vectorcomplex,j,gsl_complex_add(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j),gsl_complex_mul(gsl_vector_complex_get(vectorFFT,i),gsl_complex_mul(gsl_vector_complex_get(filterFFT_aux,i),gsl_complex_rect(cos(argExp),sin(argExp))))));
							//cout<<i<<" "<<argExp<<" "<<gsl_complex_abs(gsl_vector_complex_get(vectorFFT,index_vector))<<" "<<gsl_complex_abs(gsl_vector_complex_get(filterFFT_aux,index_vector))<<" "<<GSL_REAL(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))<<","<<GSL_IMAG(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))<<endl;
						}
						gsl_vector_set(calculatedEnergy_vector,j,gsl_complex_abs(gsl_vector_complex_get(calculatedEnergy_vectorcomplex,j))/nrmfctr);
						//cout<<"lag="<<gsl_vector_get(lags_vector,j)<<", E="<<gsl_vector_get(calculatedEnergy_vector,j)<<endl;
					}

					if ((pulseGrade < 0) || (runEMethod == 0))		*calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
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

					gsl_vector_complex_free(filterFFT_aux);
				}
				else if (mode == 1)
				{
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

					if ((pulseGrade < 0) || (runEMethod == 0))		*calculatedEnergy = gsl_vector_get(calculatedEnergy_vector,0);
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

				gsl_vector_complex_free(calculatedEnergy_vectorcomplex);
				gsl_vector_complex_free(vectorFFT);
			}
		}

		gsl_vector_free(lags_vector);
		gsl_vector_free(calculatedEnergy_vector);
	}
	else if (runEMethod == 2) //WEIGHT
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

		tempv = gsl_vector_subvector(reconstruct_init->library_collection->pulse_templates[indexEalpha].ptemplate,0,vector->size);
		gsl_vector_memcpy(Salpha,&tempv.vector);
		gsl_vector_memcpy(D,vector);
		gsl_vector_sub(D,Salpha);

		gsl_matrix_get_row(vector_aux1,reconstruct_init->library_collection->Y,indexEalpha);
		tempv = gsl_vector_subvector(vector_aux1,0,vector->size);
		gsl_vector_memcpy(Y,&tempv.vector);
		gsl_blas_ddot(D,Y,&scalar_aux1);
		scalar_aux1 = 2*scalar_aux1;

		gsl_matrix_get_row(vector_auxlong,reconstruct_init->library_collection->X,indexEalpha);
		vector2matrix(vector_auxlong,&matrix_auxlong);
		tempm = gsl_matrix_submatrix(matrix_auxlong,0,0,vector->size,vector->size);
		gsl_matrix_memcpy(X,&tempm.matrix);
		gsl_blas_dgemv(CblasNoTrans, 1.0, X, D, 0.0, vector_aux);
		gsl_blas_ddot(D,vector_aux,&scalar_aux2);

		scalar_aux1 = scalar_aux1-scalar_aux2;
		scalar_aux1 = 3*scalar_aux1/gsl_vector_get(reconstruct_init->library_collection->r,indexEalpha);

		gsl_matrix_get_row(vector_aux1,reconstruct_init->library_collection->Z,indexEalpha);
		tempv = gsl_vector_subvector(vector_aux1,0,vector->size);
		gsl_vector_memcpy(Z,&tempv.vector);
		gsl_blas_ddot(D,Z,&scalar_aux3);
		scalar_aux2 = pow(2*scalar_aux3-1,2.0);

		scalar_aux1 = sqrt(scalar_aux1+scalar_aux2);

		scalar_aux3 = 2*scalar_aux3-1;

		scalar_aux1 = (scalar_aux3 + scalar_aux1)*(gsl_vector_get(reconstruct_init->library_collection->r,indexEalpha)/3)*(gsl_vector_get(reconstruct_init->library_collection->energies,indexEbeta)-gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha));

		*calculatedEnergy = gsl_vector_get(reconstruct_init->library_collection->energies,indexEalpha) + scalar_aux1;
		cout<<"*calculatedEnergy: "<<*calculatedEnergy<<endl;

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
	else if (runEMethod == 3) //WEIGHTN
	{
		cout<<"Paso1"<<endl;
		gsl_vector *P_Pab = gsl_vector_alloc(reconstruct_init->pulse_length);
		gsl_vector *Dab = gsl_vector_alloc(reconstruct_init->pulse_length);
		gsl_matrix *Wm = gsl_matrix_alloc(reconstruct_init->pulse_length,reconstruct_init->pulse_length);
		gsl_vector *Wv = gsl_vector_alloc(reconstruct_init->pulse_length*reconstruct_init->pulse_length);
		gsl_vector *WDab = gsl_vector_alloc(reconstruct_init->pulse_length);
		double DabWDab;
		gsl_vector *Dab_mod = gsl_vector_alloc(reconstruct_init->pulse_length);
		gsl_vector *WP_Pab = gsl_vector_alloc(reconstruct_init->pulse_length);

		gsl_matrix_get_row(Dab,reconstruct_init->library_collection->DAB,indexEalpha);
		cout<<"Paso3"<<endl;

		gsl_matrix_get_row(P_Pab,reconstruct_init->library_collection->PAB,indexEalpha);
		cout<<"Paso4"<<endl;

		gsl_matrix_get_row(Wv,reconstruct_init->library_collection->W,indexEalpha);
		vector2matrix(Wv,&Wm);
		gsl_blas_dgemv(CblasNoTrans, 1.0, Wm, Dab, 0.0, WDab);
		gsl_blas_ddot(Dab,WDab,&DabWDab);
		cout<<"Paso5"<<endl;

		gsl_vector_memcpy(Dab_mod,Dab);
		gsl_vector_scale(Dab_mod,1/DabWDab);
		cout<<"Paso6"<<endl;

		gsl_blas_dgemv(CblasNoTrans, 1.0, Wm, P_Pab, 0.0, WDab);
		gsl_blas_ddot(Dab_mod,WP_Pab,calculatedEnergy);
		cout<<"Paso7"<<endl;
		cout<<"*calculatedEnergy: "<<*calculatedEnergy<<endl;

		gsl_vector_free(P_Pab);
		gsl_vector_free(Dab);
		gsl_matrix_free(Wm);
		gsl_vector_free(Wv);
		gsl_vector_free(WDab);
		gsl_vector_free(Dab_mod);
		gsl_vector_free(WP_Pab);
		cout<<"Paso8"<<endl;
	}

    return EPOK;
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* writeFilter: This function ...
*
******************************************************************************/
int writeFilter(ReconstructInitSIRENA *reconstruct_init, double normalizationFactor, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject, const char *create)
{
	string message = "";
	int status = EPOK;

	fitsfile *fltObject;
	char fltName[256];
	strncpy(fltName,reconstruct_init->filterFile,255);
	fltName[255]='\0';

	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	// Create _dtc file (if file already exists => check clobber)
	if (fileExists(string(fltName)) && (reconstruct_init->clobber == 1))
	{
		if (remove(fltName))
		{
			message = "Output filter file already exists & cannot be deleted ("+string(strerror(errno))+")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	    }
	}
	else if (fileExists(string(fltName)) && (reconstruct_init->clobber == 0))
	{
		message = "Output filter file already exists: must not be overwritten";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	if(!fileExists(string(fltName)))
	{
	    if(fits_create_file(&fltObject, fltName, &status))
	    {
	    	message = "Cannot create file " + string(fltName);
	    	EP_PRINT_ERROR(message,status); return(EPFAIL);
	    }
	}

	if (fits_open_file(&fltObject,fltName,READWRITE,&status))
	{
		message = "Cannot open output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// FILTER HDU
	strcpy(extname,"FILTER");
	if (fits_create_tbl(fltObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
	{
		message = "Cannot create table " + string(extname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	strcpy(extname,"FILTER");
	if (fits_movnam_hdu(fltObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(keyname,"NRMFCTR");
	if (fits_write_key(fltObject,TDOUBLE,keyname,&normalizationFactor,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Primary HDU
	strcpy(extname,"Primary");
	int *hdutype;
	if (fits_movabs_hdu(fltObject, 1, hdutype, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(keyname,"CREATOR");
	string creator (string("File CREATED by") + ' ' + (string) create);
	strcpy(keyvalstr,creator.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(keyname,"HISTORY");
	const char * charhistory= "HISTORY Starting parameter list";
	strcpy(keyvalstr,charhistory);
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	string strhistory (string("RecordFile = ") + reconstruct_init->record_file);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("TesEventFile = ") + reconstruct_init->event_file;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("LibraryFile = ") + reconstruct_init->library_file;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("NoiseFile = ") + reconstruct_init->noise_file;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_mode[125];			sprintf(str_mode,"%d",reconstruct_init->mode);
	strhistory = string("mode = ") + string(str_mode);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_crtLib[125];		sprintf(str_crtLib,"%d",reconstruct_init->crtLib);
	strhistory=string("crtLib = ") + string(str_crtLib);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_lastELibrary[125];		sprintf(str_lastELibrary,"%d",reconstruct_init->lastELibrary);
	strhistory=string("lastELibrary = ") + string(str_lastELibrary);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("FilterDomain = ") + reconstruct_init->FilterDomain;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("FilterMethod = ") + reconstruct_init->FilterMethod;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("EnergyMethod = ") + reconstruct_init->EnergyMethod;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_calibLQ[125];      sprintf(str_calibLQ,"%d",reconstruct_init->calibLQ);
	strhistory=string("calibLQ = ") + string(str_calibLQ);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_b_cF[125];	sprintf(str_b_cF,"%f",reconstruct_init->b_cF);
	strhistory=string("b_cF = ") + string(str_b_cF);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_c_cF[125];	sprintf(str_c_cF,"%f",reconstruct_init->c_cF);
	strhistory=string("c_cF = ") + string(str_c_cF);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_maxPulsesPerRecord[125];	sprintf(str_maxPulsesPerRecord,"%d",reconstruct_init->maxPulsesPerRecord);
	strhistory=string("maxPulsesPerRecord = ") + string(str_maxPulsesPerRecord);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_pulse_length[125];	sprintf(str_pulse_length,"%d",reconstruct_init->pulse_length);
	strhistory=string("PulseLength = ") + string(str_pulse_length);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_tauFall[125];		sprintf(str_tauFall,"%e",reconstruct_init->tauFall);
	strhistory=string("tauFall = ") + string(str_tauFall);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",reconstruct_init->scaleFactor);
	strhistory=string("scaleFactor = ") + string(str_scaleFactor);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_samplesUp[125];	sprintf(str_samplesUp,"%f",reconstruct_init->samplesUp);
	strhistory=string("samplesUp = ") + string(str_samplesUp);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_nSgms[125];	    sprintf(str_nSgms,"%f",reconstruct_init->nSgms);
	strhistory=string("nSgms = ") + string(str_nSgms);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_baseline[125];	sprintf(str_baseline,"%f",reconstruct_init->baseline);
	strhistory=string("baseline = ") + string(str_baseline);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_LrsT[125];			sprintf(str_LrsT,"%e",reconstruct_init->LrsT);
	strhistory=string("LrsT = ") + string(str_LrsT);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_LbT[125];			sprintf(str_LbT,"%e",reconstruct_init->LbT);
	strhistory=string("LbT = ") + string(str_LbT);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_monoenergy[125];	sprintf(str_monoenergy,"%f",reconstruct_init->monoenergy);
	strhistory=string("monoenergy = ") + string(str_monoenergy);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_intermediate[125];      sprintf(str_intermediate,"%d",reconstruct_init->intermediate);
	strhistory=string("intermediate = ") + string(str_intermediate);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("detectFile = ") + reconstruct_init->detectFile;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("RecordFileCalib2 = ") + reconstruct_init->record_file2;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_monoenergy2[125];	sprintf(str_monoenergy2,"%f",reconstruct_init->monoenergy2);
	strhistory=string("monoenergy2 = ") + string(str_monoenergy2);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strhistory=string("filterFile = ") + reconstruct_init->filterFile;
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	char str_clobber[125];      sprintf(str_clobber,"%d",reconstruct_init->clobber);
	strhistory=string("clobber = ") + string(str_clobber);
	strcpy(keyvalstr,strhistory.c_str());
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	charhistory= "HISTORY Ending parameter list";
	strcpy(keyvalstr,charhistory);
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	IOData obj;
	obj.inObject = fltObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"FILTER");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];

	// OPTIMALF column
	obj.iniRow = 1;
	obj.endRow = optimalfilter->size;
	strcpy(obj.nameCol,"OPTIMALF");
	strcpy(obj.unit," ");
	if (writeFitsSimple (obj,optimalfilter))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// FREQ column
	strcpy(obj.nameCol,"Freq");
	strcpy(obj.unit,"Hz");
	if (writeFitsSimple (obj,optimalfilter_f))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// OPTIMALFF column
	strcpy(obj.nameCol,"OptimalFF");
	strcpy(obj.unit," ");
	if (writeFitsSimple (obj,optimalfilter_FFT))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Write output keywords
	strcpy(keyname,"EVENTCNT");
	long keyvalint = optimalfilter->size;
	if (fits_write_key(fltObject,TLONG,keyname,&keyvalint,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(keyname,"EVENTSZ");
	keyvalint = 1;
	if (fits_write_key(fltObject,TLONG,keyname,&keyvalint,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(keyname,"NRMFCTR");
	if (fits_write_key(fltObject,TDOUBLE,keyname,&normalizationFactor,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(fltName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	if (fits_close_file(fltObject,&status))
	{
	    message = "Cannot close file " + string(fltName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* writeUCEnergy: This function ...
*
******************************************************************************/
int writeUCEnergy(ReconstructInitSIRENA **reconstruct_init, PulsesCollection *pulsesAll, PulsesCollection *pulsesInRecord, int pulse_index, double uncE, fitsfile **dtcObject, const char *create)
{
	string message = "";
	int status = EPOK;

	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	char extname[20];
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	if (fits_open_file(dtcObject,dtcName,READWRITE,&status))
	{
		message = "Cannot open output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	IOData obj;
	obj.inObject = *dtcObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"PULSES");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"keV");

	// UNCE column
	obj.iniRow = pulsesAll->ndetpulses + pulse_index + 1;
	obj.endRow = pulsesAll->ndetpulses + pulse_index + 1;
	gsl_vector *uncEgsl = gsl_vector_alloc(1);
	gsl_vector_set(uncEgsl,0,uncE);
	gsl_vector_scale(uncEgsl,1.0/1e3);
	strcpy(obj.nameCol,"UNCE");
	if (writeFitsSimple (obj,uncEgsl))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in PULSES";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(uncEgsl);

	if (((*reconstruct_init)->clobber == 1) && (pulse_index == 0))
	{
		strcpy(extname,"Primary");
		int *hdutype;
		if (fits_movabs_hdu(*dtcObject, 1, hdutype, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		string mod1 (string("File MODIFIED by") + ' ' +	(string) create);

		strcpy(keyname,"MOD0");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		if (fits_update_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if ((*reconstruct_init)->mode == 0) (*reconstruct_init)->clobber = 2;
	}

	if (fits_close_file(*dtcObject,&status))
	{
	    message = "Cannot close file " + string(dtcName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* writeFilterHDU: This function ...
*
******************************************************************************/
int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init, int pulse_index, double normalizationFactor, double uncE, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject, const char *create)
{
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
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	if (fits_open_file(dtcObject,dtcName,READWRITE,&status))
	{
		message = "Cannot open output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	if ( ((*reconstruct_init)->clobber == 1) && (pulse_index == 0))
	{
		strcpy(extname,"FILTER");
		if (fits_create_tbl(*dtcObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
	}

	strcpy(extname,"FILTER");
	if (fits_movnam_hdu(*dtcObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
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

	// NRMFCTR column
	gsl_vector *nrmfctr = gsl_vector_alloc(1);
	gsl_vector_set(nrmfctr,0,normalizationFactor);
	strcpy(obj.nameCol,"NRMFCTR");
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

	// UNCE column
	gsl_vector *uncEgsl = gsl_vector_alloc(1);
	gsl_vector_set(uncEgsl,0,uncE);
	gsl_vector_scale(uncEgsl,1.0/1e3);
	strcpy(obj.nameCol,"UNCE");
	strcpy(obj.unit,"keV");
	if (writeFitsSimple (obj,uncEgsl))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in PULSES";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	gsl_vector_free(uncEgsl);

	if (fits_get_num_rows(*dtcObject,&totalpulses, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (((*reconstruct_init)->clobber == 1) && (pulse_index == 0))
	{
		strcpy(extname,"Primary");
		int *hdutype;
		if (fits_movabs_hdu(*dtcObject, 1, hdutype, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		string mod1 (string("File MODIFIED by") + ' ' +	(string) create);

		strcpy(keyname,"MOD0");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(*dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		if (fits_update_key_longstr(*dtcObject,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	if (fits_close_file(*dtcObject,&status))
	{
	    message = "Cannot close file " + string(dtcName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION C ************************************************************
* runEnergy: This function ...
*
******************************************************************************/
void runEnergy(ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
{
	const char * create= "runEnergy v.13.0.0";	//Set "CREATOR" keyword of output FITS file

	string message="";

	if (calculateEnergy(*reconstruct_init,pulsesInRecord))
	{
		message = "Cannot run calculateEnergy in runEnergy";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	if ((*reconstruct_init)->intermediate == 1)
	{
		if (writeEnergy(reconstruct_init,*pulsesInRecord, create))
		{
			message = "Cannot run writeEnergy in runEnergy";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	return;
}
/*xxxx end of SECTION C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION CX ************************************************************
* loadUCEnergies: This function...
****************************************************************************/
int loadUCEnergies(ReconstructInitSIRENA *reconstruct_init, PulsesCollection *pulsesAll, long *nz, gsl_vector **zi, double *E0z)
{
        *E0z = reconstruct_init->monoenergy/1e3;  // keV

        *nz = 0;

        gsl_vector *ziAUX = gsl_vector_alloc(pulsesAll->ndetpulses);
        gsl_vector_set_zero(ziAUX);

        for (int i=0;i<pulsesAll->ndetpulses;i++)
        {
                if (pulsesAll->pulses_detected[i].quality == 0)
                {
                        gsl_vector_set(ziAUX, *nz, pulsesAll->pulses_detected[i].ucenergy);  // Already in keV
                        *nz = *nz + 1;
                }
        }

        gsl_vector_view temp;
        *zi = gsl_vector_alloc(*nz);
        gsl_vector_subvector(ziAUX,0,*nz);
        gsl_vector_memcpy(*zi,&temp.vector);

        gsl_vector_free(ziAUX);

        return(EPOK);
}
/*xxxx end of SECTION CX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION CX ************************************************************
* calculus_bc function:
****************************************************************************/
int calculus_bc (int calibLQ, long nx, gsl_vector *xi, double E0x, long ny, gsl_vector *yj, double E0y, double *b_cF, double *c_cF)
{
	// Declare variables
	long nx_OK = 0;
	long ny_OK = 0;
	double x_= 0.0;
	double x2_= 0.0;
	double y_= 0.0;
	double y2_ = 0.0;

	//int relEe = 0;  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (calibLQ == 1)
	{
		for (int i=0; i < nx; i++)
		{
			if (gsl_vector_get(xi,i) != 0.0)
			{
				x_ = x_ + gsl_vector_get(xi,i);
				nx_OK = nx_OK +1;
			}
		}
		x_ = x_/nx_OK;

		*b_cF = E0x/x_;
		*c_cF = 0.0;
	}
	else if (calibLQ == 2)
	{
		for (int i=0; i < nx; i++)
		{
			if (gsl_vector_get(xi,i) != 0.0)
			{
				x_ = x_ + gsl_vector_get(xi,i);
				x2_ = x2_ + pow(gsl_vector_get(xi,i),2.0);
				nx_OK = nx_OK +1;
			}
		}
		x_ = x_/nx_OK;
		x2_ = x2_/nx_OK;

		for (int j=0; j < ny; j++)
		{
			if (gsl_vector_get(yj,j) != 0.0)
			{
				y_ = y_ + gsl_vector_get(yj,j);
				y2_ = y2_ + pow(gsl_vector_get(yj,j),2.0);
				ny_OK = ny_OK +1;
			}
		}
		y_ = y_/ny_OK;
		y2_ = y2_/ny_OK;

		*b_cF = (E0x*y2_ - E0y*x2_)/(x_*y2_ - y_*x2_);
		*c_cF = (E0y*x_ - E0x*y_)/(x_*y2_ - y_*x2_);
	}

	return (EPOK);
}
/*xxxx end of SECTION CX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION CX ************************************************************
* calculateEnergy: This function converts uncalibrated energies (e) into energies (E).
*
* It reads the uncalibrated energies of the found pulses from the 'pulses' structure and it write the energies in the same structure.
*
* It uses the info in the 'reconstruct_init' structure: calilbration factors ('b_cF' and 'c_cF') and the 'calibLQ' element to know
* if the calibration of the energies is linear or quadratic:
*
* Linear: E = b·e
* Quadratic: E = b·e + c·e²
****************************************************************************/
int calculateEnergy(ReconstructInitSIRENA *reconstruct_init, PulsesCollection **pulses)
{
	double ucenergy;
	double b = reconstruct_init->b_cF;
	double c = reconstruct_init->c_cF;
	int calibLQ = reconstruct_init->calibLQ;

	for (int i=0;i<(*pulses)->ndetpulses;i++)
	{
		ucenergy = (*pulses)->pulses_detected[i].ucenergy;

		if (calibLQ == 1)	// Linear
		{
			(*pulses)->pulses_detected[i].energy = b*ucenergy;	// In SIXTE, SIGNAL is in keV
		}
		else if (calibLQ == 2)	// Quadratic
		{
			(*pulses)->pulses_detected[i].energy = b*ucenergy + c*pow(ucenergy,2.0);	// In SIXTE, SIGNAL is in keV
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION CX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION CX ************************************************************
* writeEnergy: This function...
****************************************************************************/
int writeEnergy(ReconstructInitSIRENA **reconstruct_init, PulsesCollection *pulsesInRecord, const char *create)
{
	string message = "";
	int status = EPOK;

	long totalpulses = 0;

	fitsfile *dtcObject;
	char dtcName[256];
	strncpy(dtcName,(*reconstruct_init)->detectFile,255);
	dtcName[255]='\0';

	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	if (fits_open_file(&dtcObject,dtcName,READWRITE,&status))
	{
		message = "Cannot open output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(extname,"PULSES");
	if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(dtcObject,&totalpulses, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	IOData obj;
	obj.inObject = dtcObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"PULSES");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	for (int i=0;i<pulsesInRecord->ndetpulses;i++)
	{
		// ENERGY column
		obj.iniRow = totalpulses + i + 1 -pulsesInRecord->ndetpulses;
		obj.endRow = totalpulses + i + 1 -pulsesInRecord->ndetpulses;
		strcpy(obj.nameCol,"ENERGY");
		strcpy(obj.unit,"keV");
		gsl_vector *energy = gsl_vector_alloc(1);
		gsl_vector_set(energy,0,pulsesInRecord->pulses_detected[i].energy);
		strcpy(obj.nameCol,"ENERGY");
		if (writeFitsSimple (obj,energy))
		{
		    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) +
		    		" column in FILTER";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		gsl_vector_free(energy);
	}

	if ((*reconstruct_init)->clobber == 1)
	{
		strcpy(extname,"Primary");
		int *hdutype;
		if (fits_movabs_hdu(dtcObject, 1, hdutype, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		string mod1 (string("File MODIFIED by") + ' ' +	(string) create);

		strcpy(keyname,"MOD1");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(dtcObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(dtcName);
		    EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		if (fits_update_key_longstr(dtcObject,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		(*reconstruct_init)->clobber = 2;
	}

	if (fits_close_file(dtcObject,&status))
	{
	    message = "Cannot close file " + string(dtcName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	return(EPOK);
}
/*xxxx end of SECTION CX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
