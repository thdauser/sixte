#include "tasksSIRENA.h"

/***** SECTION A ************************************************************
* runDetect: This function ... Run detection routine in SIRENA, record by record
*
******************************************************************************/
void runDetect(TesRecord* record, int lastRecord_runDetect, PulsesCollection *pulsesAll, ReconstructInitSIRENA** reconstruct_init, PulsesCollection** pulsesInRecord)
{
	const char * create= "runDetect v.17.0.0";	//Set "CREATOR" keyword of output FITS file

	string message="";
	int status=EPOK;

	// Declare variables
	fitsfile *inLibObject = NULL;	// Object which contains information of the library FITS file
	bool appendToLibrary = false;	// Pulse templates library FITS file new (appendToLibrary=false) or not (appendToLibrary=true)

	fitsfile *dtcObject = NULL;	    // Object which contains information of the output FITS file
	char dtcName[256];
	strcpy(dtcName,(*reconstruct_init)->detectFile);

	gsl_matrix *library=NULL;		// ENERGY - PHEIGHT
	gsl_matrix *models=NULL;		// Matrix to store all the pulse templates (PULSES) of the library

	int eventsz = record->trigger_size;
	double tstartRecord;
	gsl_vector *invector = gsl_vector_alloc(eventsz);	// Record

	// Handle the library data or the library file
	if (handleLibraryDetect(*reconstruct_init,&appendToLibrary, lastRecord_runDetect, &inLibObject,create,&library,&models))
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
		if (filderLibrary(reconstruct_init,1/record->delta_t,&models))
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
	if (procRecord(reconstruct_init, tstartRecord, 1/record->delta_t, dtcObject, invector, library, models, *pulsesInRecord))
	{
	    message = "Cannot run routine procRecord for record processing";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	int extver=0;
	char extname[20];

	
	if (((*reconstruct_init)->intermediate == 1) && (lastRecord_runDetect == 1))
	{
		// Write output keywords (their values have been previously checked)
		char keyname[10];
		char *comment=NULL;
		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(dtcObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) +" in " + string(dtcName);
			EP_EXIT_ERROR(message,status);
		}
		
		long totalpulses;
		if (fits_get_num_rows(dtcObject,&totalpulses, &status))
		{
			message = "Cannot get number of rows in " + string(dtcName);
			EP_EXIT_ERROR(message,status);
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
			EP_EXIT_ERROR(message,status);
		}
	}

	if ((lastRecord_runDetect == 1) && (*reconstruct_init)->crtLib == 1)	// CREATIONLIB run mode => Calculate the pulse template by averaging some found pulses
	{
		gsl_vector *pulsetemplate = gsl_vector_alloc((*reconstruct_init)->pulse_length);
		double pulseheighttemplate = 0;

		char extname[20];
		char keyname[10];
		char *comment=NULL;

		if (calculateTemplate (*reconstruct_init, pulsesAll, *pulsesInRecord, 1/record->delta_t, &pulsetemplate,&pulseheighttemplate))
		{
		    message = "Cannot run routine calculateTemplate in creationlib run mode";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
		
		if (writeLibrary(*reconstruct_init, pulseheighttemplate, &pulsetemplate, appendToLibrary, &inLibObject))
		{
		    message = "Cannot run routine writeLibrary in crationlib run mode";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		gsl_vector_free(pulsetemplate);
	}

	if ((*reconstruct_init)->intermediate == 1)
	{
		if (fits_close_file(dtcObject,&status))
		{
			message = "Cannot close file " + string(dtcName);
			EP_EXIT_ERROR(message,status);
		}
	}
	if(library != NULL) gsl_matrix_free(library);
	if(models != NULL) gsl_matrix_free(models);
	gsl_vector_free(invector);
	return;
}
/*xxxx end of SECTION A xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/***** SECTION A1 ************************************************************
* handleLibraryDetect: This function ...
*
******************************************************************************/
int handleLibraryDetect(ReconstructInitSIRENA* reconstruct_init_handleLibraryDetect, bool *appendToLibrary_handleLibraryDetect, int lastRecord_handleLibrary, fitsfile **inLibObject_handleLibraryDetect, const char* create_handleLibraryDetect, gsl_matrix **library_handleLibraryDetect, gsl_matrix **models_handleLibraryDetect)
{
	// Declare variables
	string message = "";
	int crtLib = reconstruct_init_handleLibraryDetect->crtLib;

	// If CREATIONLIB run mode (crtLib = 1) => Secondary pulses are not searched for => Pulse templates library not used => Library has to be created
	// If NOTCREATIONLIB mode (crtLib = 0) => Read the pulse templates library input FITS file
	if (crtLib == 0)
	{
		// Initialize the vectors/matrix of the library
		if (initLibraryDetect(reconstruct_init_handleLibraryDetect, library_handleLibraryDetect, models_handleLibraryDetect))
		{
			message = "Cannot run routine initLibraryDetect to initialize library";
			EP_EXIT_ERROR(message,EPFAIL);
		}

	    // Store the library data, which are stored in reconstruct_init->library_collection, in gsl vectors/matrix
		if (loadLibraryDetect(reconstruct_init_handleLibraryDetect,library_handleLibraryDetect,models_handleLibraryDetect))
		{
		   	message = "Cannot run routine readLib to read pulses library";
		  	EP_EXIT_ERROR(message,EPFAIL);
		}
	}
	else
	{
		if (lastRecord_handleLibrary == 1)
		{
			if (createLibraryDetect(reconstruct_init_handleLibraryDetect, appendToLibrary_handleLibraryDetect, inLibObject_handleLibraryDetect, create_handleLibraryDetect, library_handleLibraryDetect, models_handleLibraryDetect))
			{
				message = "Cannot run routine createLibraryDetect to create pulses library";
				EP_EXIT_ERROR(message,EPFAIL);
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION A1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A2 ************************************************************
* initLibraryDetect: This function ...
*
******************************************************************************/
int initLibraryDetect(ReconstructInitSIRENA* reconstruct_init_initLibraryDetect, gsl_matrix **library_initLibraryDetect, gsl_matrix **models_initLibraryDetect)
{
	int nummodels = reconstruct_init_initLibraryDetect->library_collection->ntemplates;
	int eventsz_initLibraryDetect = reconstruct_init_initLibraryDetect->library_collection->pulse_templates->template_duration;

	*library_initLibraryDetect = gsl_matrix_alloc(nummodels,2);					// reconstruct_init->LibraryCollection->energies + reconstruct_init->LibraryCollection->pulse_heights
	*models_initLibraryDetect = gsl_matrix_alloc(nummodels,eventsz_initLibraryDetect);	// All the pulse templates

	return(EPOK);
}
/*xxxx end of SECTION A2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A3 ************************************************************
* loadLibraryDetect:
*
******************************************************************************/
int loadLibraryDetect(ReconstructInitSIRENA* reconstruct_init_loadLibraryDetect, gsl_matrix **library_loadLibraryDetect, gsl_matrix **models_loadLibraryDetect)
{
	for (int i=0;i<reconstruct_init_loadLibraryDetect->library_collection->ntemplates;i++)
	{
		gsl_matrix_set(*library_loadLibraryDetect,i,0,(reconstruct_init_loadLibraryDetect->library_collection->energies[i]));
		gsl_matrix_set(*library_loadLibraryDetect,i,1,(reconstruct_init_loadLibraryDetect->library_collection->pulse_heights[i]));
		for (int j=0;j<reconstruct_init_loadLibraryDetect->library_collection->pulse_templates->template_duration;j++)
		{
			gsl_matrix_set(*models_loadLibraryDetect,i,j,(reconstruct_init_loadLibraryDetect->library_collection->pulse_templates[i].ptemplate[j]));
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION A3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A4 ************************************************************
* createLibraryDetect: This function...
****************************************************************************/
int createLibraryDetect(ReconstructInitSIRENA* reconstruct_init_createLibraryDetect, bool *appendToLibrary_createLibraryDetect, fitsfile **inLibObject_createLibraryDetect, const char * create_createLibraryDetect, gsl_matrix **library_createLibraryDetect, gsl_matrix **models_createLibraryDetect)
{
	int status = EPOK;
	int extver=0;
	string message = "";

	char extname[20];
	char keyname[10];
	char *comment=NULL;

	//fitsfile *inLibObject;
	char inLibName[200];
	strcpy(inLibName, reconstruct_init_createLibraryDetect->library_file);

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
		*appendToLibrary_createLibraryDetect = true;

		if (fits_open_file(inLibObject_createLibraryDetect,inLibName,READWRITE,&status))
		{
		    message = "Cannot open library file " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(extname,"LIBRARY");
		if (fits_movnam_hdu(*inLibObject_createLibraryDetect, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Initialize the vectors/matrix of the library
		//if (initLibraryDetect(reconstruct_init_createLibraryDetect, librarygsl, modelsgsl, modelsb0gsl))
		if (initLibraryDetect(reconstruct_init_createLibraryDetect, library_createLibraryDetect, models_createLibraryDetect))
		{
			message = "Cannot run routine initLibraryDetect to initialize library";
			EP_EXIT_ERROR(message,EPFAIL);
		}

		// Store library, which is stored in reconstruct_init->library_collection, in gsl vectors/matrix
		//if (loadLibraryDetect(reconstruct_init_createLibraryDetect, &librarygsl, &modelsgsl, &modelsb0gsl))
		if (loadLibraryDetect(reconstruct_init_createLibraryDetect, library_createLibraryDetect, models_createLibraryDetect))
		{
		    message = "Cannot run routine loadLibraryDetect to read pulse templates library";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
	}
	else
	{
		*appendToLibrary_createLibraryDetect = false;
		status = EPOK;
		if (fits_create_file(inLibObject_createLibraryDetect, inLibName, &status))
		{
		    message = "Cannot create library file " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (fits_open_file(inLibObject_createLibraryDetect,inLibName,READWRITE,&status))
		{
		    message = "Cannot open library file " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Create extension LIBRARY
		strcpy(extname,"LIBRARY");
		if (fits_create_tbl(*inLibObject_createLibraryDetect,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
		{
		    message = "Cannot create table " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(extname,"LIBRARY");
		if (fits_movnam_hdu(*inLibObject_createLibraryDetect, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		eventcntLib = 1;
		strcpy(keyname,"EVENTCNT");
		if (eventcntLib <= 0)
		{
		    message = "Legal values for EVENTCNT (LIBRARY) are integer numbers greater than 0";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		if (fits_write_key(*inLibObject_createLibraryDetect,TLONG,keyname,&eventcntLib,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"CREATOR");
		strcpy(keyvalstr,create_createLibraryDetect);
		if (fits_write_key(*inLibObject_createLibraryDetect,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	char str_energy[125];       sprintf(str_energy,"%f",reconstruct_init_createLibraryDetect->monoenergy);

	string process (string("") 	+ ' ' +
	string(inLibName) 	  + ' ' +
	string(str_energy)      + ' ' +
	string("(")				+      (string) create_createLibraryDetect 		  +   	  string(")"));

	strcpy(keyname,"PROC0");
	strcpy(keyvalstr,process.c_str());
	if (fits_write_key_longwarn(*inLibObject_createLibraryDetect,&status))
	{
	    message = "Cannot write keyword long warn in library";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_write_key_longstr(*inLibObject_createLibraryDetect,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in library";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	return (EPOK);
}
/*xxxx end of SECTION A4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A5 ************************************************************
* createDetectFile function: This function ...
***************************************************************************/
int createDetectFile(ReconstructInitSIRENA* reconstruct_init_createDetectFile, double samprate, const char * create_createDetectFile, fitsfile **dtcObject_createDetectFile)
{
	int status = EPOK;
	string message = "";

	char dtcName[256];
	strcpy(dtcName,reconstruct_init_createDetectFile->detectFile);

	// Create output FITS file: If it does not exist yet
	// If dtcName does not finish as '.fits' and the file dtcName+'.fits' already exists =>
	// => Data are appended to dtcName file => Must not be allowed
	// '.fits' => 5 characters
	if (strlen(dtcName)<6)
	{
		// dtcName has 5 or less characters => Does not contain '.fits' =>Append '.fits' to dtcName
		char dtcNameaux[255];
		sprintf(dtcNameaux,dtcName);
		strcat(dtcNameaux,".fits");
		strcpy(dtcName,dtcNameaux);
	}
	else if (strlen(dtcName)>=6)
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
	if (fileExists(string(dtcName)) && (reconstruct_init_createDetectFile->clobber == 1))
	{
		if (remove(dtcName))
		{
			message = "Output detect file already exists & cannot be deleted ("+string(strerror(errno))+")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	    }
	}
	else if (fileExists(string(dtcName)) && (reconstruct_init_createDetectFile->clobber == 0))
	{
		message = "Output detect file already exists: must not be overwritten";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	if (!fileExists(string(dtcName)))
	{
		if(fits_create_file(dtcObject_createDetectFile, dtcName, &status))
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

	if (fits_open_file(dtcObject_createDetectFile,dtcName,READWRITE,&status))
	{
	    message = "Cannot open output detect file " + string(dtcName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (reconstruct_init_createDetectFile->clobber == 1)
	{
		strcpy(extname,"PULSES");
		if (fits_create_tbl(*dtcObject_createDetectFile, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Create keywords
		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(*dtcObject_createDetectFile, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"MODE");
		if ((reconstruct_init_createDetectFile->mode != 0) && (reconstruct_init_createDetectFile->mode != 1))
		{
			message = "Legal values for MODE (PULSES) are 0 or 1";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_write_key(*dtcObject_createDetectFile,TINT,keyname,&(reconstruct_init_createDetectFile->mode),comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(keyname,"EVENTSZ");
		//if (sizePulse_b <= 0)
		if (reconstruct_init_createDetectFile->pulse_length <= 0)
		{
			message = "Legal values for EVENTSZ (PULSES) are integer numbers greater than 0";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		//if (fits_write_key(dtcObject,TINT,keyname,&sizePulse_b,comment,&status))
		if (fits_write_key(*dtcObject_createDetectFile,TINT,keyname,&(reconstruct_init_createDetectFile->pulse_length),comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		if (reconstruct_init_createDetectFile->crtLib == 1)
		{
			strcpy(keyname,"ENERGY");
			if (reconstruct_init_createDetectFile-> monoenergy < 0)
			{
				message = "Legal values for ENERGY (PULSES) are non negative real numbers";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			if (fits_write_key(*dtcObject_createDetectFile,TDOUBLE,keyname,&(reconstruct_init_createDetectFile-> monoenergy),comment,&status))
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
		if (fits_write_key(*dtcObject_createDetectFile,TDOUBLE,keyname,&samprate,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(keyname,"CREATOR");
		strcpy(keyvalstr,create_createDetectFile);
		if (fits_write_key(*dtcObject_createDetectFile,TSTRING,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Set PROCESS keyword
		char str_pulse_length[125];	sprintf(str_pulse_length,"%d",reconstruct_init_createDetectFile->pulse_length);
		char str_tauFall[125];		sprintf(str_tauFall,"%e",reconstruct_init_createDetectFile->tauFall);
		char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",reconstruct_init_createDetectFile->scaleFactor);
		char str_samplesUp[125];	sprintf(str_samplesUp,"%d",reconstruct_init_createDetectFile->samplesUp);
		char str_nSgms[125];	    sprintf(str_nSgms,"%f",reconstruct_init_createDetectFile->nSgms);
		char str_mode[125];			sprintf(str_mode,"%d",reconstruct_init_createDetectFile->mode);
		char str_crtLib[125];		sprintf(str_crtLib,"%d",reconstruct_init_createDetectFile->crtLib);
		char str_LrsT[125];			sprintf(str_LrsT,"%e",reconstruct_init_createDetectFile->LrsT);
		char str_LbT[125];			sprintf(str_LbT,"%e",reconstruct_init_createDetectFile->LbT);
		char str_clobber[125];      sprintf(str_clobber,"%d",reconstruct_init_createDetectFile->clobber);

		string process (string(" runDetect") 		+ ' ' +
		string(reconstruct_init_createDetectFile->record_file)          + ' ' + string(reconstruct_init_createDetectFile->library_file)		+ ' ' + string(dtcName)		+ ' ' +
		string(str_mode) 				+ ' ' + string(str_crtLib) 		+ ' ' +
		string(str_LrsT)   	    		+ ' ' + string(str_LbT)         + ' ' +
		string(str_tauFall)     		+ ' ' + string(str_scaleFactor) + ' ' +
		string(str_samplesUp)   		+ ' ' + string(str_nSgms)       + ' ' +
		string(str_pulse_length) 		+ ' ' +
		string(str_clobber) + ' ' +
		string("(") + (string) create_createDetectFile + string(")"));

		strcpy(keyname,"PROC0");
		strcpy(keyvalstr,process.c_str());
		if (fits_write_key_longwarn(*dtcObject_createDetectFile,&status))
		{
			message = "Cannot write long warn in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		if (fits_write_key_longstr(*dtcObject_createDetectFile,keyname,keyvalstr,comment,&status))
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
int filderLibrary(ReconstructInitSIRENA** reconstruct_init_filderLibrary, double samprate, gsl_matrix **models_filderLibrary)
{
	int status = EPOK;

	// NOTCREATIONLIB mode and first record
	if (((*reconstruct_init_filderLibrary)->crtLib == 0) && ((*reconstruct_init_filderLibrary)->library_collection->pulse_templates_filder->template_duration == -1))
	{
		string message = "";

		double scaleFactor = (*reconstruct_init_filderLibrary)->scaleFactor;
		double tauFALL = (*reconstruct_init_filderLibrary)->tauFall;

		// Check boxLength
		double cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
		int boxLength = (int) ((1/cutFreq) * samprate);
		if (boxLength <= 1)
		{
			message = "lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
	        // Cómo se escribe un warning desde una función?
		}

		// 'models' is filtered and derived
		gsl_vector *model = gsl_vector_alloc((*reconstruct_init_filderLibrary)->library_collection->pulse_templates->template_duration);
		for (int i=0; i<(*reconstruct_init_filderLibrary)->library_collection->ntemplates; i++)
		{
			gsl_matrix_get_row(model,*models_filderLibrary,i);

			// PULSE TEMPLATE: Low-pass filtering
			status = lpf_boxcar(&model,model->size,tauFALL*scaleFactor,samprate);
			if (status == 1)
			{
				message = "Cannot run routine lpf_boxcar for low-pass filtering";
				EP_EXIT_ERROR(message,status); return(EPFAIL);
			}
			if (status == 3)
			{
				status = EPOK;
			}
			if (status == 4)
			{
				message = "lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low";
				EP_EXIT_ERROR(message,status); return(EPFAIL);
			}

			// PULSE TEMPLATE: Derivative after filtering
			if (derivative (&model,model->size))
			{
				message = "Cannot run routine derMTHSimple to calculate derivative";
				EP_EXIT_ERROR(message,EPFAIL);
			}

			for (int j=0; j<(*reconstruct_init_filderLibrary)->library_collection->pulse_templates->template_duration;j++)
			{
				(*reconstruct_init_filderLibrary)->library_collection->pulse_templates_filder[i].ptemplate[j] = gsl_vector_get(model,j);
			}

			(*reconstruct_init_filderLibrary)->library_collection->pulse_templates_filder[i].template_duration = (*reconstruct_init_filderLibrary)->library_collection->pulse_templates->template_duration;

			gsl_matrix_set_row(*models_filderLibrary,i,model);
		}

		gsl_vector_free(model);
	}
	else
	{
		for (int i=0; i<(*reconstruct_init_filderLibrary)->library_collection->ntemplates; i++)
		{
			for (int j=0; j<(*reconstruct_init_filderLibrary)->library_collection->pulse_templates->template_duration;j++)
			{
				gsl_matrix_set(*models_filderLibrary,i,j,(*reconstruct_init_filderLibrary)->library_collection->pulse_templates_filder[i].ptemplate[j]);
			}
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION A6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A7 ************************************************************
* loadRecord:
*
******************************************************************************/
int loadRecord(TesRecord* record_loadRecord, double *time_record, gsl_vector **adc_double)
{
	*time_record = record_loadRecord->time;
	for (int i=0;i<record_loadRecord->trigger_size;i++)
	{
		gsl_vector_set(*adc_double,i,record_loadRecord->adc_double[i]);
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
int procRecord(ReconstructInitSIRENA** reconstruct_init_procRecord, double tstartRecord_procRecord, double samprate, fitsfile *dtcObject_procRecord, gsl_vector *record, gsl_matrix *library_procRecord, gsl_matrix *models_procRecord, PulsesCollection *foundPulses)
{
	int status = EPOK;
	string message = "";

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Provisional => To be deleted in future
	FILE * temporalFile;
	char temporalFileName[256];
	sprintf(temporalFileName,"auxfile");
	strcat(temporalFileName,".txt");
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// Initialize variables
	int numPulses = 0;

	double asquid = 1.0;
	double plspolar = 1.0;
	double ivcal=1.0;
	int chngplrt = 0;// 1 => Polarity changed (pulses multiplied by -1)
    				// 0 => Polarity not changed

	double stopCriteriaMKC = 1.0;  			// Used in medianKappaClipping
	                               	   	   	// Given in %
	double kappaMKC = 3.0;					// Used in medianKappaClipping
	double levelPrvPulse = 100.0;  		    // Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse

	double scaleFactor = (*reconstruct_init_procRecord)->scaleFactor;
	double tauFALL = (*reconstruct_init_procRecord)->tauFall;
	int sizePulse_b = (*reconstruct_init_procRecord)->pulse_length;
	double samplesUp = (*reconstruct_init_procRecord)->samplesUp;
	double nSgms = (*reconstruct_init_procRecord)->nSgms;
	double Lrs = (int) ((*reconstruct_init_procRecord)->LrsT*samprate);			// Running sum length (in the RS filter case): LrsT in samples
	double Lb = (int) ((*reconstruct_init_procRecord)->LbT*samprate); 				// Baseline averaging length (in the RS filter case): LbT in samples

	// Allocate GSL vectors
	gsl_vector *recordNOTFILTERED = gsl_vector_alloc(record->size); 	// Record without having been filtered
	gsl_vector *recordFILTERED = gsl_vector_alloc(record->size);	 	// Filtered (LPF) record
	gsl_vector *recordDERIVATIVE = gsl_vector_alloc(record->size);  	// Derivative of invectorFILTERED

	// To look for pulses
	gsl_vector *tstartgsl = gsl_vector_alloc(record->size);
	gsl_vector *tendgsl = gsl_vector_alloc(record->size);
	gsl_vector *qualitygsl = gsl_vector_alloc(record->size);
	gsl_vector *pulseHeightsgsl = gsl_vector_alloc(record->size);
	gsl_vector_set_zero(qualitygsl);
	gsl_vector_set_zero(pulseHeightsgsl);						// In order to choose the proper pulse model to calculate
	                                                            // the adjusted derivative and to fill in the ESTENRGY column
	                                                            // in the output FITS file

	gsl_vector_scale(record,ivcal);			// IVCAL to change arbitrary units of voltage to non-arbitrary
			                                // units of current (Amps)

	// Assign positive polarity to the pulses
	if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))
	{
		gsl_vector_scale(record,-1);
		chngplrt = 1;
	}

	char dtcName[256];
	strcpy(dtcName,(*reconstruct_init_procRecord)->detectFile);
	if ((*reconstruct_init_procRecord)->intermediate == 1)
	{
		if (fits_open_file(&dtcObject_procRecord,dtcName,1,&status))
		{
			message = "Cannot open file " +  string(dtcName);
			EP_EXIT_ERROR(message,status);
		}
		char extname[20];
		char keyname[10];
		int extver=0;
		char *comment=NULL;
		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(dtcObject_procRecord, ANY_HDU,extname, extver, &status))
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
		if (fits_update_key(dtcObject_procRecord,TINT,keyname,&chngplrt,comment,&status))
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
		//????????????????????
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
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_vector_memcpy(recordFILTERED,record);

	// Derivative after filtering
	if (derivative (&record, record->size))
	{
	    message = "Cannot run routine derivative for derivative after filtering";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gsl_vector_memcpy(recordDERIVATIVE,record);

	// Find pulses of the record
	if ((*reconstruct_init_procRecord)->crtLib == 0)
	{
		if (findPulses (recordNOTFILTERED, recordDERIVATIVE, &tstartgsl, &qualitygsl, &pulseHeightsgsl,
				&numPulses,
				1,
				tauFALL, scaleFactor, sizePulse_b, samprate,
				samplesUp, nSgms,
				Lb, Lrs,
				library_procRecord, models_procRecord,
				stopCriteriaMKC,
				kappaMKC,
				levelPrvPulse,
				temporalFile))
		{
			message = "Cannot run routine findPulses";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	else if ((*reconstruct_init_procRecord)->crtLib == 1)
	{
		if (findPulses (recordNOTFILTERED, recordDERIVATIVE, &tstartgsl, &qualitygsl, &pulseHeightsgsl,
				&numPulses,
				0,
				tauFALL, scaleFactor, sizePulse_b, samprate,
				samplesUp, nSgms,
				Lb, Lrs,
				library_procRecord, models_procRecord,
				stopCriteriaMKC,
				kappaMKC,
				levelPrvPulse,
				temporalFile))
		{
			message = "Cannot run routine findPulses";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	for (int i=0;i<numPulses;i++)
	{
		gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i)+sizePulse_b); 	//tend_i = tstart_i + (ntaus*tauFALL*samprate)

		if ((gsl_vector_get(qualitygsl,i) != 1) && (gsl_vector_get(qualitygsl,i) != 3))
		// If it is already truncated at the beginning, it is not taken into account to classify it again as truncated (at the end)
		{
			if (gsl_vector_get(tendgsl,i) >= recordDERIVATIVE->size)	// Truncated pulses at the end of the record
			{
				gsl_vector_set(tendgsl,i,(recordDERIVATIVE->size)-1);
				gsl_vector_set (qualitygsl,i,1);
			}
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
	gsl_vector *tauRisegsl = gsl_vector_alloc(recordDERIVATIVE->size);
	gsl_vector *tauFallgsl = gsl_vector_alloc(recordDERIVATIVE->size);
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
		}
		else
		{
			foundPulses->pulses_detected[i].grade2 = (*reconstruct_init_procRecord)->pulse_length;
		}
		foundPulses->pulses_detected[i].pulse_adc = new double[sizePulse_b];
		for (int j=0;j<foundPulses->pulses_detected[i].pulse_duration;j++)
		{

			foundPulses->pulses_detected[i].pulse_adc[j]=gsl_vector_get(recordNOTFILTERED,gsl_vector_get(tstartgsl,i)+j);
		}
		foundPulses->pulses_detected[i].Tstart = gsl_vector_get(tstartgsl,i)/samprate+tstartRecord_procRecord;
		foundPulses->pulses_detected[i].Tend = gsl_vector_get(tendgsl,i)/samprate+tstartRecord_procRecord;
		foundPulses->pulses_detected[i].riseTime = gsl_vector_get(tauRisegsl,i);
		foundPulses->pulses_detected[i].fallTime = gsl_vector_get(tauFallgsl,i);
		foundPulses->pulses_detected[i].pulse_height = gsl_vector_get(pulseHeightsgsl,i);
		// 'ucenergy' and 'energy' will be known after running runFilter and runEnergy respectively
		foundPulses->pulses_detected[i].quality = gsl_vector_get(qualitygsl,i);
	}

	// Write pulses info in output FITS file
	if ((*reconstruct_init_procRecord)->intermediate == 1)
	{
		if (writePulses (reconstruct_init_procRecord, samprate, tstartRecord_procRecord, recordNOTFILTERED, numPulses, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, pulseHeightsgsl, dtcObject_procRecord))
		{
			message = "Cannot run routine writePulses to write pulses in output file";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	// Free allocate of GSL vectors
	gsl_vector_free(recordNOTFILTERED);
	gsl_vector_free(recordFILTERED);
	gsl_vector_free(recordDERIVATIVE);

	gsl_vector_free(tstartgsl);
	gsl_vector_free(tendgsl);
	gsl_vector_free(qualitygsl);
	gsl_vector_free(pulseHeightsgsl);

	gsl_vector_free(tauRisegsl);
	gsl_vector_free(tauFallgsl);

	return EPOK;
}
/*xxxx end of SECTION A8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A9 ************************************************************
* writePulses function:
*
******************************************************************************/
int writePulses(ReconstructInitSIRENA** reconstruct_init_writePulses, double samprate, double initialtime, gsl_vector *invectorNOTFIL, int numPulsesRecord, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, gsl_vector *pulseheights, fitsfile *dtcObject_writePulses)
{
	int status = EPOK;
	string message = "";

	// Declare variables
	int t0;   					// First value of index of pulse
	gsl_matrix *vgslout2;

	gsl_vector_view temp;

	char dtcName[256];
	strcpy(dtcName,(*reconstruct_init_writePulses)->detectFile);

    // If intermediate=1 => First record, createDetectFile
    //	                    Change clobber to 2
	//                   => Not first record, append info to output dtc file (because clobber is 2)
	long totalpulses;	// It is necessary to know the row to write info
	if ((*reconstruct_init_writePulses)->clobber == 1)
	{
		totalpulses = 0;
		if ((*reconstruct_init_writePulses)->crtLib == 1 )	(*reconstruct_init_writePulses)->clobber = 2;
	}
	else if ((*reconstruct_init_writePulses)->clobber == 2)
	{
		if (fits_get_num_rows(dtcObject_writePulses,&totalpulses, &status))
		{
			message = "Cannot get number of rows in " + string(dtcName);
			EP_EXIT_ERROR(message,status);
		}
	}

	if (numPulsesRecord!=0)
	{
		vgslout2 = gsl_matrix_alloc(numPulsesRecord,(*reconstruct_init_writePulses)->pulse_length);

		// Converting bins to time
		for (int i=0; i<numPulsesRecord; i++)
		{
			t0 = gsl_vector_get (tstart,i);
			gsl_vector_set(tstart,i,initialtime + (gsl_vector_get (tstart,	i) * (1/samprate)));
			gsl_vector_set(tend,i,initialtime + (gsl_vector_get (tend,	i) * (1/samprate)));

			if (invectorNOTFIL->size - t0 > (*reconstruct_init_writePulses)->pulse_length)	//The invectorNOTFIL has more bins than sizePulse
			{
				temp = gsl_vector_subvector(invectorNOTFIL,t0,(*reconstruct_init_writePulses)->pulse_length);
				gsl_matrix_set_row(vgslout2, i, &temp.vector);
			}
			else 									// The invectorNOTFIL has less bins than sizePulse (truncated)
			{
				for (int j=0; j<(invectorNOTFIL->size) - t0; j++)
				{
					if (t0 == -1) t0 = 0;
					gsl_matrix_set (vgslout2,i,j,gsl_vector_get(invectorNOTFIL,j+t0));
				}

				for (int j=(invectorNOTFIL->size)-t0; j< (*reconstruct_init_writePulses)->pulse_length; j++) {gsl_matrix_set (vgslout2,i,j,0.0);}
			}
		}

		IOData obj;

		// Creating TSTART Column
		obj.inObject = dtcObject_writePulses;
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
		strcpy(obj.unit,"Amps");
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

		// Creating PHEIGHT Column
		strcpy(obj.nameCol,"PHEIGHT");
		strcpy(obj.unit,"a.u.");
		temp = gsl_vector_subvector(pulseheights,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column PHEIGHT";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating QUALITY Column
		strcpy(obj.nameCol,"QUALITY");
		obj.type = TSHORT;
		strcpy(obj.unit,"bits");
		temp = gsl_vector_subvector(quality,0,numPulsesRecord);
		if (writeFitsSimple(obj, &temp.vector))
		{
		    message = "Cannot run routine writeFitsSimple for column QUALITY";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		totalpulses = totalpulses + numPulsesRecord;

		// Free allocate GSL vectors
		gsl_matrix_free (vgslout2);
	}

	return (EPOK);
}
/*xxxx end of SECTION A9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION  A10 ************************************************************
* calculateTemplate function:
******************************************************************************/
int calculateTemplate (ReconstructInitSIRENA *reconstruct_init_calculateTemplate, PulsesCollection *pulsesAll_calculateTemplate, PulsesCollection *pulsesInRecord_calculateTemplate, double samprate, gsl_vector **pulseaverage, double *pulseaverageHeight)
{
	int status = EPOK;
	string message = "";
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Provisional => To be deleted in future
	FILE * temporalFile;
	char temporalFileName[256];
	sprintf(temporalFileName,"auxfile");
	strcat(temporalFileName,".txt");
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// Declare (and initialize) variables
	int totalPulses = pulsesAll_calculateTemplate->ndetpulses + pulsesInRecord_calculateTemplate->ndetpulses;
	gsl_vector *tstart = gsl_vector_alloc(totalPulses);			// Tstart column from the output dtc FITS file
	gsl_vector *pulseheight = gsl_vector_alloc(totalPulses);	// EstEnrgy column from the output dtc FITS file
	gsl_vector *quality = gsl_vector_alloc(totalPulses);		// Quality column from the output dtc FITS file
	for (int i=0;i<pulsesAll_calculateTemplate->ndetpulses;i++)
	{
		gsl_vector_set(tstart,i,pulsesAll_calculateTemplate->pulses_detected[i].Tstart);
		gsl_vector_set(pulseheight,i,pulsesAll_calculateTemplate->pulses_detected[i].pulse_height);
		gsl_vector_set(quality,i,pulsesAll_calculateTemplate->pulses_detected[i].quality);
	}
	for (int i=0;i<pulsesInRecord_calculateTemplate->ndetpulses;i++)
	{
		gsl_vector_set(tstart,i+pulsesAll_calculateTemplate->ndetpulses,pulsesInRecord_calculateTemplate->pulses_detected[i].Tstart);
		gsl_vector_set(pulseheight,i+pulsesAll_calculateTemplate->ndetpulses,pulsesInRecord_calculateTemplate->pulses_detected[i].pulse_height);
		gsl_vector_set(quality,i+pulsesAll_calculateTemplate->ndetpulses,pulsesInRecord_calculateTemplate->pulses_detected[i].quality);
	}

	gsl_vector *nonpileup = gsl_vector_alloc(totalPulses);	// Piled up pulse => Not taking into account to calculate the template
	long nonpileupPulses = totalPulses;						// A priori, all the found pulses are considered as non piled up
	gsl_vector_set_all(nonpileup,1);

	int nBins = floor(sqrt(totalPulses));			// Number of bins of the pulseheights histogram
													// Square-root choice (used by Excel and many others)
	gsl_vector *xhisto = gsl_vector_alloc(nBins);	// X-axis of the pulseheights histogram
	gsl_vector *yhisto = gsl_vector_alloc(nBins);	// Y-axis of the pulseheights histogram
	int index_maximumpulseheight;					// Index where the maximum of the pulseheights histogram is
	double maximumpulseheight;						// Maximum of the pulseheights histogram

	bool firstnonpileupPulse = true;
	gsl_vector *pulse = gsl_vector_alloc(reconstruct_init_calculateTemplate->pulse_length);

	double tstartnext;

	// Only if one of the found pulses is validated as non piled up pulse by using EstEnrgy
	// (pulseheights histogram) and Tstart, and Quality => The pulse, I0, is going
	// to be read from the dtc output FITS file (in order to not handle a long matrix of
	// pulses, found pulses x sizePulse_b)

	gsl_vector_scale(tstart,samprate); 	//tstarts not in sec but in samples

	// Create pulseheights histogram
	if (createHisto(pulseheight, nBins, &xhisto, &yhisto))
	{
	    message = "Cannot run createHisto routine";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	index_maximumpulseheight = gsl_vector_max_index(yhisto);
	maximumpulseheight = gsl_vector_get(xhisto,index_maximumpulseheight);

	for (int i=0;i<totalPulses;i++)
	{
		if (i == totalPulses-1)		tstartnext = gsl_vector_get(tstart,i)+2*reconstruct_init_calculateTemplate->pulse_length;
		else 						tstartnext = gsl_vector_get(tstart,i+1);

		// Check if the pulse is piled up or not
		if ((gsl_vector_get(pulseheight,i) < maximumpulseheight-0.1*maximumpulseheight) || (gsl_vector_get(pulseheight,i) > maximumpulseheight+0.1*maximumpulseheight)
			|| (tstartnext-gsl_vector_get(tstart,i) <= reconstruct_init_calculateTemplate->pulse_length) || (gsl_vector_get(quality,i) >= 1))
		{
 			gsl_vector_set(nonpileup,i,0);
			nonpileupPulses --;
		}
		else
		{
			if (i < pulsesAll_calculateTemplate->ndetpulses)
			{
				for (int j=0;j<reconstruct_init_calculateTemplate->pulse_length;j++)
				{
					gsl_vector_set(pulse,j,(pulsesAll_calculateTemplate->pulses_detected[i].pulse_adc[j]));
				}
			}
			else
			{
				for (int j=0;j<reconstruct_init_calculateTemplate->pulse_length;j++)
				{
					gsl_vector_set(pulse,j,(pulsesInRecord_calculateTemplate->pulses_detected[i-pulsesAll_calculateTemplate->ndetpulses].pulse_adc[j]));
				}
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

	//Just in case due to the noise influence in the alignment, the first sample of the pulseaverage is not around the baseline but higher
	double meanLast200points, sgLast200points;
	gsl_vector_view temp;
	temp = gsl_vector_subvector(*pulseaverage,reconstruct_init_calculateTemplate->pulse_length-200-1,200);
	if (findMeanSigma (&temp.vector, &meanLast200points, &sgLast200points, temporalFile))
	{
		message = "Cannot run findMeanSigma routine for kappa-sigma iteration";
		EP_PRINT_ERROR(message,EPFAIL);
	}
	if (fabs(gsl_vector_get(*pulseaverage,0))>fabs(meanLast200points)+3*sgLast200points)
	{
		gsl_vector * pulseaverageAUX = gsl_vector_alloc(reconstruct_init_calculateTemplate->pulse_length);
		gsl_vector_set(pulseaverageAUX,0,meanLast200points);
		for (int i=1;i<reconstruct_init_calculateTemplate->pulse_length;i++)
		{
			gsl_vector_set(pulseaverageAUX,i,gsl_vector_get(*pulseaverage,i-1));
		}
		gsl_vector_memcpy(*pulseaverage,pulseaverageAUX);
		gsl_vector_free(pulseaverageAUX);
	}

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
/*xxxx end of SECTION A10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A11 ************************************************************
* createHisto function: This function builds ...
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

    // Free allocate of GSL vectors
    gsl_vector_free(invectoraux);
    gsl_vector_free(invectoraux2);

    return EPOK;
}
/*xxxx end of SECTION A11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A12 ************************************************************
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

	return (EPOK);
}
/*xxxx end of SECTION A12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A13 ************************************************************
* shiftm function: This function returns (in vectorout) the vectorin delayed m samples
*
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
/*xxxx end of SECTION A13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A14 ************************************************************
* shift_m function: This function returns (in vectorout) the vectorin moved forward m samples
*
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
/*xxxx end of SECTION A14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION A15 ************************************************************
* writeLibrary function:
*
******************************************************************************/
int writeLibrary(ReconstructInitSIRENA *reconstruct_init_writeLibrary, double estenergy, gsl_vector **pulsetemplate, bool appendToLibrary_writeLibrary, fitsfile **inLibObject_writeLibrary)
{
	int status = EPOK;
	int extver=0;
	string message = "";

	char inLibName[200];
	strcpy(inLibName, reconstruct_init_writeLibrary->library_file);

	// Declare variables
	gsl_vector *energyoutgsl = gsl_vector_alloc(1);
	gsl_vector *estenergyoutgsl = gsl_vector_alloc(1);
	gsl_matrix *pulsetemplates_matrix = gsl_matrix_alloc(1,reconstruct_init_writeLibrary->pulse_length);
	gsl_matrix *pulsetemplatesb0_matrix = gsl_matrix_alloc(1,reconstruct_init_writeLibrary->pulse_length);
	gsl_matrix *matchedfilters_matrix = gsl_matrix_alloc(1,reconstruct_init_writeLibrary->pulse_length);
	gsl_matrix *matchedfiltersb0_matrix = gsl_matrix_alloc(1,reconstruct_init_writeLibrary->pulse_length);

	char keyname[10];
	char *comment=NULL;
	char extname[20];
	IOData obj;

    if (appendToLibrary_writeLibrary == true)
    {
    	long eventcntLib;
    	strcpy(keyname,"EVENTCNT");
    	if (fits_read_key(*inLibObject_writeLibrary,TLONG,keyname, &eventcntLib,comment,&status))
    	{
    	    message = "Cannot read keyword " + string(keyname);
    	    EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}
    	if (eventcntLib <= 0)
    	{
    		message = "Legal values for read EVENTCNT (LIBRARY) are integer numbers greater than 0";
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}
    	int eventcntLib1 = eventcntLib + 1;
    	strcpy(keyname,"EVENTCNT");
    	if (eventcntLib1 <= eventcntLib)
    	{
    		message = "Legal value for written EVENTCNT (LIBRARY) is read EVENTCNT (LIBRARY) plus 1";
    		EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}
    	if (fits_update_key(*inLibObject_writeLibrary,TLONG,keyname, &eventcntLib1,comment,&status))
    	{
    	    message = "Cannot update keyword " + string(keyname);
    	    EP_PRINT_ERROR(message,status); return(EPFAIL);
    	}

    	gsl_vector *energycolumn = gsl_vector_alloc(eventcntLib+1);
    	gsl_vector *estenergycolumn = gsl_vector_alloc(eventcntLib+1);

    	gsl_vector *modelsrow = gsl_vector_alloc(reconstruct_init_writeLibrary->pulse_length);
    	gsl_vector *modelsrowb0 = gsl_vector_alloc(reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *modelsaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *modelsb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *modelsaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *modelsb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);

    	gsl_vector *matchedfiltersrow = gsl_vector_alloc(reconstruct_init_writeLibrary->pulse_length);
    	gsl_vector *matchedfiltersrowb0 = gsl_vector_alloc(reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *matchedfiltersaux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *matchedfiltersb0aux = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *matchedfiltersaux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);
    	gsl_matrix *matchedfiltersb0aux1 = gsl_matrix_alloc(eventcntLib+1,reconstruct_init_writeLibrary->pulse_length);

    	for (int i=0;i<eventcntLib;i++)
    	{
    		gsl_vector_set(energycolumn,i,reconstruct_init_writeLibrary->library_collection->energies[i]);
    		gsl_vector_set(estenergycolumn,i,reconstruct_init_writeLibrary->library_collection->pulse_heights[i]);
    		for (int j=0;j<reconstruct_init_writeLibrary->pulse_length;j++)
    		{
    			gsl_matrix_set(modelsaux,i,j,reconstruct_init_writeLibrary->library_collection->pulse_templates[i].ptemplate[j]);
    			gsl_matrix_set(modelsb0aux,i,j,reconstruct_init_writeLibrary->library_collection->pulse_templates_B0[i].ptemplate[j]);
    			gsl_matrix_set(matchedfiltersaux,i,j,reconstruct_init_writeLibrary->library_collection->matched_filters[i].mfilter[j]);
    		    gsl_matrix_set(matchedfiltersb0aux,i,j,reconstruct_init_writeLibrary->library_collection->matched_filters_B0[i].mfilter[j]);
    		}
    	}

    	gsl_vector_set(energycolumn,eventcntLib,reconstruct_init_writeLibrary->monoenergy);
    	gsl_vector_set(estenergycolumn,eventcntLib,estenergy);
    	gsl_matrix_set_row(modelsaux,eventcntLib,*pulsetemplate);

    	gsl_vector *row_aux= gsl_vector_alloc((*pulsetemplate)->size);
    	gsl_vector_memcpy(row_aux,*pulsetemplate);
    	gsl_vector_scale(row_aux,1/reconstruct_init_writeLibrary->monoenergy);
    	gsl_matrix_set_row(matchedfiltersaux,eventcntLib,row_aux);

    	gsl_vector *baseline = gsl_vector_alloc((*pulsetemplate)->size);
    	gsl_vector_set_all(baseline,-200.0);
    	gsl_vector_add(*pulsetemplate,baseline);
    	gsl_vector_free(baseline);
    	gsl_matrix_set_row(pulsetemplatesb0_matrix,0,*pulsetemplate);
    	gsl_matrix_set_row(modelsb0aux,eventcntLib,*pulsetemplate);

    	gsl_vector_memcpy(row_aux,*pulsetemplate);
    	gsl_vector_scale(row_aux,1/reconstruct_init_writeLibrary->monoenergy);
    	gsl_matrix_set_row(matchedfiltersb0aux,eventcntLib,row_aux);
    	gsl_vector_free(row_aux);

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
    	}
    	gsl_vector_memcpy(energycolumn,energycolumnaux);
    	gsl_vector_memcpy(estenergycolumn,estenergycolumnaux);
    	gsl_matrix_memcpy(modelsaux,modelsaux1);
    	gsl_matrix_memcpy(modelsb0aux,modelsb0aux1);
    	gsl_matrix_memcpy(matchedfiltersaux,matchedfiltersaux1);
    	gsl_matrix_memcpy(matchedfiltersb0aux,matchedfiltersb0aux1);

    	gsl_permutation_free(perm);
    	gsl_vector_free(energycolumnaux);
    	gsl_vector_free(estenergycolumnaux);
    	gsl_matrix_free(modelsaux1);
    	gsl_matrix_free(modelsb0aux1);
    	gsl_matrix_free(matchedfiltersaux1);
    	gsl_matrix_free(matchedfiltersb0aux1);

    	obj.inObject = *inLibObject_writeLibrary;
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
    	    strcpy(obj.unit,"Amps");
    	    gsl_matrix_get_row(modelsrow,modelsaux,i);
    	    gsl_matrix_set_row(pulsetemplates_matrix,0,modelsrow);
    	    if (writeFitsComplex(obj, pulsetemplates_matrix))
    	    {
    	    	message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    	    	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	    }

    	    strcpy(obj.nameCol,"PULSEB0");
    	    strcpy(obj.unit,"Amps");
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
       	    //gsl_vector_scale(matchedfiltersrow,1./reconstruct_init_writeLibrary->monoenergy);
       	    gsl_matrix_set_row(matchedfilters_matrix,0,matchedfiltersrow);
       	    if (writeFitsComplex(obj, matchedfilters_matrix))
       	    {
      	    	message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
       	    	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
       	    }

       	    strcpy(obj.nameCol,"MFB0");
       	    strcpy(obj.unit," ");
       	    gsl_matrix_get_row(matchedfiltersrowb0,matchedfiltersb0aux,i);
       	    //gsl_vector_scale(matchedfiltersrowb0,1./reconstruct_init_writeLibrary->monoenergy);
       	    gsl_matrix_set_row(matchedfiltersb0_matrix,0,matchedfiltersrowb0);
       	    if (writeFitsComplex(obj, matchedfiltersb0_matrix))
       	    {
       	      	message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
       	      	EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
       	    }

    	}
    	gsl_vector_free(modelsrow);
    	gsl_vector_free(modelsrowb0);
    	gsl_vector_free(matchedfiltersrow);
    	gsl_vector_free(matchedfiltersrowb0);

    	gsl_vector_free (energycolumn);
    	gsl_vector_free (estenergycolumn);
    	gsl_matrix_free(modelsaux);
    	gsl_matrix_free(modelsb0aux);
    	gsl_matrix_free(matchedfiltersaux);
    	gsl_matrix_free(matchedfiltersb0aux);
    }
    else
    {
    	strcpy(extname,"LIBRARY");
    	if (fits_movnam_hdu(*inLibObject_writeLibrary, ANY_HDU,extname, extver, &status))
    	{
    		message = "Cannot move to HDU  " + string(extname) + " in library";
    		EP_PRINT_ERROR(message,status);return(EPFAIL);
    	}

    	if (reconstruct_init_writeLibrary->pulse_length <= 0)
    	{
    		message = "Legal values for EVENTSZ (PULSES) are integer numbers greater than 0";
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}
    	if (fits_write_key(*inLibObject_writeLibrary,TINT,"EVENTSZ",&reconstruct_init_writeLibrary->pulse_length,comment,&status))
    	{
    		message = "Cannot write keyword EVENTSZ in library";
    		EP_PRINT_ERROR(message,status);return(EPFAIL);
    	}

    	gsl_vector_set (energyoutgsl,0,reconstruct_init_writeLibrary->monoenergy);
    	gsl_vector_set (estenergyoutgsl,0,estenergy);

    	gsl_matrix_set_row(pulsetemplates_matrix,0,*pulsetemplate);

    	gsl_vector *baseline = gsl_vector_alloc((*pulsetemplate)->size);
    	gsl_vector_set_all(baseline,-200.0);
    	gsl_vector_add(*pulsetemplate,baseline);
    	gsl_vector_free(baseline);
    	gsl_matrix_set_row(pulsetemplatesb0_matrix,0,*pulsetemplate);

    	// Creating ENERGY Column
    	obj.inObject = *inLibObject_writeLibrary;
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
    	strcpy(obj.unit,"Amps");
    	if (writeFitsComplex(obj, pulsetemplates_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating PULSEB0 Column
    	strcpy(obj.nameCol,"PULSEB0");
    	strcpy(obj.unit,"Amps");
    	if (writeFitsComplex(obj, pulsetemplatesb0_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating MF Column
    	strcpy(obj.nameCol,"MF");
    	strcpy(obj.unit," ");
    	gsl_matrix_scale(pulsetemplates_matrix,1.0/reconstruct_init_writeLibrary->monoenergy);
    	if (writeFitsComplex(obj, pulsetemplates_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}

    	// Creating MFB0 Column
    	strcpy(obj.nameCol,"MFB0");
    	strcpy(obj.unit," ");
    	gsl_matrix_scale(pulsetemplatesb0_matrix,1.0/reconstruct_init_writeLibrary->monoenergy);
    	if (writeFitsComplex(obj, pulsetemplatesb0_matrix))
    	{
    		message = "Cannot run writeFitsComplex routine for column " + string(obj.nameCol);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}
    }

    if (fits_close_file(*inLibObject_writeLibrary,&status))
    {
    	message = "Cannot close file " + string(inLibName);
        EP_EXIT_ERROR(message,status);
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
/*xxxx end of SECTION 15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION B ************************************************************
* runFilter: This function ...
*
******************************************************************************/
void runFilter(TesRecord* record, int nRecord_runFilter, int lastRecord_runFilter, ReconstructInitSIRENA** reconstruct_init, PulsesCollection *pulsesAll_runFilter, PulsesCollection** pulsesInRecord, OptimalFilterSIRENA **optimalFilter_runFilter)
{
	const char * create= "runFilter v.15.0.0";	//Set "CREATOR" keyword of output FITS file

	string message="";
	int status = EPOK;

	fitsfile *dtcObject = NULL;	    // Object which contains information of the output FITS file
	if((*reconstruct_init)->intermediate==1){
	    char dtcName[256];
	    strcpy(dtcName,(*reconstruct_init)->detectFile);
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

	long energyInLibrary_row;

	double normalizationFactor;
	double uncE;

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Provisional => To be deleted in future
	FILE * temporalFile;
	char temporalFileName[256];
	sprintf(temporalFileName,"auxfile");
	strcat(temporalFileName,".txt");
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// Declare variables
	gsl_vector *energylibrary;
	gsl_vector *estenergylibrary;
	gsl_matrix *templateslibrary;
	gsl_matrix *templateslibraryb0;
	gsl_matrix *matchedfilters;
	gsl_matrix *matchedfiltersb0;

	gsl_vector *freqgsl;
	gsl_vector *csdgsl;

	gsl_vector *optimalfilter_SHORT;	// Resized optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector *optimalfilter_f_SHORT;		// Resized optimal filter f's when f's are according to [0,...fmax,-fmax,...] (frequency domain)
	gsl_vector *optimalfilter_FFT_SHORT;	// Resized optimal filter magnitudes when f's are according to [0,...fmax,-fmax,...] (frequency domain)

	gsl_vector *model;

	gsl_vector *recordAux;

	gsl_vector *pulse=NULL;

	gsl_vector *filtergsl=NULL;			// Matched filter values (time domain)

	int iter;
	gsl_vector_view(temp);

	double tstartSamplesRecord;
	double tstartRecord;
	double tstartRecordSamples = floor(record->time/record->delta_t+0.5);	// Close integer

	// Initialize the vectors/matrix of the library
	if (initLibraryFilter(*reconstruct_init,&energylibrary, &estenergylibrary, &templateslibrary, &templateslibraryb0,&matchedfilters, &matchedfiltersb0))
	{
		message = "Cannot run routine initLibraryFilter to initialize library";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Store the library data, which are stored in reconstruct_init->library_collection, in gsl vectors/matrix
	if (loadLibraryFilter(*reconstruct_init, &energylibrary, &estenergylibrary, &templateslibrary, &templateslibraryb0,&matchedfilters, &matchedfiltersb0))
	{
	   	message = "Cannot run routine readLib to read pulses library";
	  	EP_EXIT_ERROR(message,EPFAIL);
	}

	// Initialize the vectors of the noise spectrum
	if (initNoisespectrum(*reconstruct_init,&freqgsl, &csdgsl))
	{
		message = "Cannot run routine initNoisespectrum to initialize noise spectrum";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Store the noise spectrum data, which are stored in reconstruct_init->library_collection, in gsl vectors
	if (loadNoisespectrum(*reconstruct_init, &freqgsl, &csdgsl))
	{
	   	message = "Cannot run routine loadNoisespectrum to read noise spectrum";
	  	EP_EXIT_ERROR(message,EPFAIL);
	}

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
		gsl_vector *baseline = gsl_vector_alloc(recordAux->size);
		gsl_vector_set_all(baseline,-200.0);
		gsl_vector_add(recordAux,baseline);
		gsl_vector_free(baseline);
	}

	if ((*reconstruct_init)->mode == 0)
	{
		gsl_vector *filtergsl = gsl_vector_alloc((*reconstruct_init)->pulse_length);			// Filter values

		if (nRecord_runFilter == 1)
		{
			// Matched filter
			if (getMatchedFilter(runF0orB0val, (*reconstruct_init)->monoenergy, energylibrary, matchedfilters, matchedfiltersb0, &filtergsl))
			{
				message = "Cannot run routine getMatchedFilter";
				EP_EXIT_ERROR(message,EPFAIL);
			}

			// Optimal filter
			if (calculus_optimalFilter (filtergsl, filtergsl->size, 1/record->delta_t, runF0orB0val, freqgsl, csdgsl, &optimalfilter_SHORT, &optimalfilter_f_SHORT, &optimalfilter_FFT_SHORT, &normalizationFactor))
			{
				message = "Cannot run routine calculus_optimalFilter";
				EP_EXIT_ERROR(message,EPFAIL);
			}

			(*optimalFilter_runFilter)->ofilter_duration = optimalfilter_SHORT->size;
			(*optimalFilter_runFilter)->nrmfctr = normalizationFactor;
			(*optimalFilter_runFilter)->ofilter = new double[(*optimalFilter_runFilter)->ofilter_duration];
			for (int i=0;i<optimalfilter_SHORT->size;i++)
			{
				(*optimalFilter_runFilter)->ofilter[i] = gsl_vector_get(optimalfilter_SHORT,i);
			}

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
			optimalfilter_SHORT = gsl_vector_alloc((*optimalFilter_runFilter)->ofilter_duration);
			normalizationFactor = (*optimalFilter_runFilter)->nrmfctr;
			for (int i=0;i<(*optimalFilter_runFilter)->ofilter_duration;i++)
			{
				gsl_vector_set(optimalfilter_SHORT,i,(*optimalFilter_runFilter)->ofilter[i]);
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
			message = "There are no valid pulses (quality == 0) in the found pulses";
			EP_EXIT_ERROR(message,EPFAIL);
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
				if (calculateUCEnergy(pulse,optimalfilter_SHORT,TorF,normalizationFactor,1/record->delta_t,&uncE))
				{
					message = "Cannot run calculateUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
					EP_PRINT_ERROR(message,EPFAIL);
				}

				if ((*reconstruct_init)->intermediate == 1)
				{
					if (writeUCEnergy(reconstruct_init, pulsesAll_runFilter, *pulsesInRecord, i, uncE, &dtcObject, create))
					{
						message = "Cannot run writeUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_PRINT_ERROR(message,EPFAIL);
					}
				}

				(*pulsesInRecord)->pulses_detected[i].grade1 = optimalfilter_SHORT->size;
				(*pulsesInRecord)->pulses_detected[i].ucenergy = uncE;
				(*pulsesInRecord)->pulses_detected[i].energy = uncE;
			}
			else
			{
				if ((*reconstruct_init)->intermediate == 1)
				{
					if (writeUCEnergy(reconstruct_init, pulsesAll_runFilter, *pulsesInRecord, i, uncE, &dtcObject, create))
					{
						message = "Cannot run writeUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_PRINT_ERROR(message,EPFAIL);
					}
				}
				(*pulsesInRecord)->pulses_detected[i].grade1 = 0;
				(*pulsesInRecord)->pulses_detected[i].ucenergy = 0.0;
				(*pulsesInRecord)->pulses_detected[i].energy = 0.0;
			}
		}

		gsl_vector_free(optimalfilter_SHORT);
	}
	else if ((*reconstruct_init)->mode == 1)
	{
		long resize_mf;

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
			message = "There are no valid pulses (quality == 0) in the found pulses";
			EP_EXIT_ERROR(message,EPFAIL);
		}

		model =gsl_vector_alloc((*reconstruct_init)->pulse_length);
		gsl_vector *filtergsl_aux = gsl_vector_alloc(matchedfilters->size2);

		for (int i=0; i<(*pulsesInRecord)->ndetpulses ;i++)
		{
			// Pulses are going to be validated by checking its quality
			if ((*pulsesInRecord)->pulses_detected[i].quality == 0)	// Neither truncated or saturated
			{
				// Filter
				if (getMatchedFilter(runF0orB0val, (*pulsesInRecord)->pulses_detected[i].pulse_height, estenergylibrary, matchedfilters, matchedfiltersb0, &filtergsl_aux))
				{
					message = "Cannot run routine getMatchedFilter";
					EP_EXIT_ERROR(message,EPFAIL);
				}

				resize_mf = (*pulsesInRecord)->pulses_detected[i].pulse_duration; // Resize the matched filter by using the tstarts
				resize_mf = pow(2,floor(log2(resize_mf)));

				filtergsl = gsl_vector_alloc(resize_mf);			// Filter values
				temp = gsl_vector_subvector(filtergsl_aux,0,resize_mf);
				gsl_vector_memcpy(filtergsl,&temp.vector);
				gsl_vector_free(filtergsl_aux);

				// Pulse
				tstartSamplesRecord = floor((*pulsesInRecord)->pulses_detected[i].Tstart/record->delta_t+0.5)-tstartRecordSamples;
				pulse = gsl_vector_alloc(resize_mf);
				temp = gsl_vector_subvector(recordAux,tstartSamplesRecord,resize_mf);
				gsl_vector_memcpy(pulse,&temp.vector);

				// Calculate the optimal filter
				if (calculus_optimalFilter (filtergsl, filtergsl->size, 1/record->delta_t, runF0orB0val, freqgsl, csdgsl, &optimalfilter_SHORT, &optimalfilter_f_SHORT, &optimalfilter_FFT_SHORT, &normalizationFactor))
				{
					message = "Cannot run routine calculus_optimalFilter";
				    EP_EXIT_ERROR(message,EPFAIL);
				}

				// Calculate the uncalibrated energy of each pulse
				if (calculateUCEnergy(pulse,optimalfilter_SHORT,TorF,normalizationFactor,1/record->delta_t,&uncE))
				{
					message = "Cannot run calculateUCEnergy routine for pulse i=" + boost::lexical_cast<std::string>(i);
					EP_PRINT_ERROR(message,EPFAIL);
				}

				// Subtract pulse model from the record
				if (find_model(uncE, energylibrary, templateslibraryb0, &model, temporalFile))
				{
				    message = "Cannot run find_model routine for pulse i=" + boost::lexical_cast<std::string>(i);
				    EP_PRINT_ERROR(message,EPFAIL);
				}
				for (int j=tstartSamplesRecord;j<tstartSamplesRecord+(*reconstruct_init)->pulse_length;j++)
				{
					gsl_vector_set(recordAux,j,gsl_vector_get(recordAux,j)-gsl_vector_get(model,j-tstartSamplesRecord));
				}

				if ((*reconstruct_init)->intermediate == 1)
				{
					if (writeFilterHDU(reconstruct_init, i, normalizationFactor, uncE, optimalfilter_SHORT, optimalfilter_f_SHORT, optimalfilter_FFT_SHORT, &dtcObject, create))
					{
						message = "Cannot run writeFilterHDU routine for pulse i=" + boost::lexical_cast<std::string>(i);
						EP_PRINT_ERROR(message,EPFAIL);
					}
				}

				(*pulsesInRecord)->pulses_detected[i].grade1 = resize_mf;
				(*pulsesInRecord)->pulses_detected[i].ucenergy = uncE;
				(*pulsesInRecord)->pulses_detected[i].energy = uncE;
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
						EP_PRINT_ERROR(message,EPFAIL);
					}
				}

				(*pulsesInRecord)->pulses_detected[i].grade1 = 0;
				(*pulsesInRecord)->pulses_detected[i].ucenergy = 0.0;
				(*pulsesInRecord)->pulses_detected[i].energy = 0.0;
			} // End if
			gsl_vector_free(optimalfilter_SHORT);
			gsl_vector_free(optimalfilter_f_SHORT);
			gsl_vector_free(optimalfilter_FFT_SHORT);
			if((*pulsesInRecord)->pulses_detected[i].quality == 0)gsl_vector_free(pulse);
			if((*pulsesInRecord)->pulses_detected[i].quality == 0)gsl_vector_free(filtergsl);
		} // End for
		gsl_vector_free(recordAux);
		gsl_vector_free(model);
	} // End if
	gsl_vector_free(energylibrary);
	gsl_vector_free(estenergylibrary);
	gsl_matrix_free(templateslibrary);
	gsl_matrix_free(templateslibraryb0);
	gsl_matrix_free(matchedfilters);
	gsl_matrix_free(matchedfiltersb0);

	gsl_vector_free(freqgsl);
	gsl_vector_free(csdgsl);

	return;
}
/*xxxx end of SECTION B xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* initLibraryFilter: This function ...
*
******************************************************************************/
int initLibraryFilter(ReconstructInitSIRENA* reconstruct_init_initLibraryFilter, gsl_vector **energylibrary_initLibraryFilter, gsl_vector **estenergylibrary_initLibraryFilter,
		gsl_matrix **templateslibrary_initLibraryFilter, gsl_matrix **templateslibraryb0_initLibraryFilter,
		gsl_matrix **matchedfilters_initLibraryFilter, gsl_matrix **matchedfiltersb0_initLibraryFilter)
{
	int nummodels = reconstruct_init_initLibraryFilter->library_collection->ntemplates;
	int eventsz_initLibraryFilter = reconstruct_init_initLibraryFilter->library_collection->pulse_templates->template_duration;

	*energylibrary_initLibraryFilter = gsl_vector_alloc(nummodels);					// reconstruct_init->LibraryCollection->energies
	*estenergylibrary_initLibraryFilter = gsl_vector_alloc(nummodels);				// reconstruct_init->LibraryCollection->pulse_heights
	*templateslibrary_initLibraryFilter = gsl_matrix_alloc(nummodels,eventsz_initLibraryFilter);	// All the pulse templates
	*templateslibraryb0_initLibraryFilter = gsl_matrix_alloc(nummodels,eventsz_initLibraryFilter);	// All the pulse templates_B0
	*matchedfilters_initLibraryFilter = gsl_matrix_alloc(nummodels,eventsz_initLibraryFilter);		// All the matched_filters
	*matchedfiltersb0_initLibraryFilter = gsl_matrix_alloc(nummodels,eventsz_initLibraryFilter);	// All the matched_filters_B0

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* loadLibraryFilter:
*
******************************************************************************/
int loadLibraryFilter(ReconstructInitSIRENA* reconstruct_init_loadLibraryFilter, gsl_vector **energylibrary_loadLibraryFilter, gsl_vector **estenergylibrary_loadLibraryFilter,
		gsl_matrix **templateslibrary_loadLibraryFilter, gsl_matrix **templateslibraryb0_loadLibraryFilter,
		gsl_matrix **matchedfilters_loadLibraryFilter, gsl_matrix **matchedfiltersb0_loadLibraryFilter)
{
	for (int i=0;i<reconstruct_init_loadLibraryFilter->library_collection->ntemplates;i++)
	{
		gsl_vector_set(*energylibrary_loadLibraryFilter,i,(reconstruct_init_loadLibraryFilter->library_collection->energies[i]));
		gsl_vector_set(*estenergylibrary_loadLibraryFilter,i,(reconstruct_init_loadLibraryFilter->library_collection->pulse_heights[i]));
		for (int j=0;j<reconstruct_init_loadLibraryFilter->library_collection->pulse_templates->template_duration;j++)
		{
			gsl_matrix_set(*templateslibrary_loadLibraryFilter,i,j,(reconstruct_init_loadLibraryFilter->library_collection->pulse_templates[i].ptemplate[j]));
			gsl_matrix_set(*templateslibraryb0_loadLibraryFilter,i,j,(reconstruct_init_loadLibraryFilter->library_collection->pulse_templates_B0[i].ptemplate[j]));
			gsl_matrix_set(*matchedfilters_loadLibraryFilter,i,j,(reconstruct_init_loadLibraryFilter->library_collection->matched_filters[i].mfilter[j]));
			gsl_matrix_set(*matchedfiltersb0_loadLibraryFilter,i,j,(reconstruct_init_loadLibraryFilter->library_collection->matched_filters_B0[i].mfilter[j]));
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* initNoisespectrum: This function ...
*
******************************************************************************/
int initNoisespectrum(ReconstructInitSIRENA* reconstruct_init_initNoisespectrum, gsl_vector **freqgsl_initNoisespectrum, gsl_vector **csdgsl_initNoisespectrum)
{
	int eventcnt_initNoisespectrum = reconstruct_init_initNoisespectrum->noise_spectrum->noise_duration*2;

	*freqgsl_initNoisespectrum = gsl_vector_alloc(eventcnt_initNoisespectrum);
	*csdgsl_initNoisespectrum = gsl_vector_alloc(eventcnt_initNoisespectrum);

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* loadNoisespectrum:
*
******************************************************************************/
int loadNoisespectrum(ReconstructInitSIRENA* reconstruct_init_loadNoisespectrum, gsl_vector **freqgsl_loadNoisespectrum, gsl_vector **csdgsl_loadNoisespectrum)
{
	for (int i=0;i<reconstruct_init_loadNoisespectrum->noise_spectrum->noise_duration*2;i++)
	{
		gsl_vector_set(*freqgsl_loadNoisespectrum,i,(reconstruct_init_loadNoisespectrum->noise_spectrum->noisefreqs[i]));
		gsl_vector_set(*csdgsl_loadNoisespectrum,i,(reconstruct_init_loadNoisespectrum->noise_spectrum->noisespec[i]));
	}

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


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
int calculus_optimalFilter(gsl_vector *matchedfiltergsl, long mf_size, double samprate, int runF0orB0val_calculusOptimalFilter, gsl_vector *freqgsl_calculusOptimalFilter, gsl_vector *csdgsl_calculusOptimalFilter, gsl_vector **optimal_filtergsl, gsl_vector **of_f, gsl_vector **of_FFT, double *normalizationFactor_calculusOptimalFilter)
{
	// FFT calculus of the filter template (MATCHED FILTER->matchedfiltergsl)
	// Declare variables
	double SelectedTimeDuration = mf_size/samprate;
	int status = EPOK;
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
	if (runF0orB0val_calculusOptimalFilter == 0)		gsl_vector_complex_set(mfFFTcomp,0,gsl_complex_rect(0.0,0.0));

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

	if (mf_size%2 == 0)	//Even
	{
		for (int i=0; i< (mf_size)/2; i++)
		{
			gsl_vector_set(mf_f,i,i/SelectedTimeDuration);
		}
		gsl_vector_set(mf_f,mf_size/2,(mf_size/2)/SelectedTimeDuration);
		for (int i=1; i<(mf_size/2); i++)
		{
			gsl_vector_set(mf_f,i+mf_size/2,(-1.0)*(i+mf_size/2-i*2)/SelectedTimeDuration);
		}
	}
	else	//Odd
	{
		for (int i=0; i< (mf_size)/2; i++)
		{
			gsl_vector_set(mf_f,i,i/SelectedTimeDuration);
		}
		gsl_vector_set(mf_f,mf_size/2,(mf_size/2)/SelectedTimeDuration);
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
	if (mf_size < freqgsl_calculusOptimalFilter->size)			// Decimate noise samples
	{
		if ((freqgsl_calculusOptimalFilter->size)%mf_size == 0)
		{
			int timesNoverMF = freqgsl_calculusOptimalFilter->size/mf_size;
			n_f = gsl_vector_alloc(mf_size);
			n_FFT = gsl_vector_alloc(mf_size);
			for (int i=0;i<n_f->size;i++)
			{
				gsl_vector_set(n_f,i,gsl_vector_get(freqgsl_calculusOptimalFilter,i*timesNoverMF));
				gsl_vector_set(n_FFT,i,gsl_vector_get(csdgsl_calculusOptimalFilter,i*timesNoverMF));
			}
		}
		else
		{
			// It is necessary to work only with the positive frequencies in order to not handle the f=0
			int noisePOS_size = floor(freqgsl_calculusOptimalFilter->size/2);
			gsl_vector *freqgsl_POS = gsl_vector_alloc(noisePOS_size);
			temp = gsl_vector_subvector(freqgsl_calculusOptimalFilter,1,noisePOS_size);
			gsl_vector_memcpy(freqgsl_POS,&temp.vector);
			gsl_vector *csdgsl_POS = gsl_vector_alloc(noisePOS_size);
			temp = gsl_vector_subvector(csdgsl_calculusOptimalFilter,1,noisePOS_size);
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
		   	gsl_vector_set(n_FFT,0,gsl_vector_get(csdgsl_calculusOptimalFilter,0));
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
	else if (mf_size > freqgsl_calculusOptimalFilter->size)		// Error
	{
		message = "Noise must have more samples than pulse";
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	else if (mf_size == freqgsl_calculusOptimalFilter->size)		// It is not necessary to do anything
	{
		n_f = gsl_vector_alloc(freqgsl_calculusOptimalFilter->size);
		n_FFT = gsl_vector_alloc(freqgsl_calculusOptimalFilter->size);
		gsl_vector_memcpy(n_f,freqgsl_calculusOptimalFilter);
		gsl_vector_memcpy(n_FFT,csdgsl_calculusOptimalFilter);
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
	gsl_vector_complex *of_FFT_complex = gsl_vector_complex_alloc(mf_size);
	for (int i=0;i<mf_size;i++)
	{
		gsl_vector_complex_set(of_FFT_complex,i,gsl_complex_div_real(gsl_vector_complex_get(mfFFTcomp_conj,i),gsl_vector_get(n_FFT_2,i)));
	}

	gsl_vector_complex_absRIFCA(*of_FFT,of_FFT_complex);

	// Calculus of the normalization factor
	*normalizationFactor_calculusOptimalFilter = 0;
	for (int i=0; i<mf_f->size; i++)
	{
		*normalizationFactor_calculusOptimalFilter = *normalizationFactor_calculusOptimalFilter + gsl_vector_get(mf_FFT_2,i)/gsl_vector_get(n_FFT_2,i);
	}

	gsl_vector_free(mf_FFT_2);
	gsl_vector_free(n_FFT_2);

	// Calculus of the pseudoenergy in frequency domain
	/*if ((mode == 1) && (TorF == 1))
	{
		gsl_vector_complex *pulseFFTcomp = gsl_vector_complex_alloc(mf_f->size);
		gsl_vector *pulse_FFT = gsl_vector_alloc(mf_f->size);				// Ordered magnitude according [-fmax,...,0,...,fmax]
		if (FFT(pulse,pulseFFTcomp,SelectedTimeDuration))
		{
		    message = "Cannot run FFT routine to calculate pulse FFT";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		gsl_vector_complex_absRIFCA(pulse_FFT,pulseFFTcomp);

		uncE = 0.0;
		for (int i=0;i<mf_f->size;i++)
		{
			uncE = uncE + gsl_vector_get(pulse_FFT,i)*gsl_complex_abs(gsl_vector_complex_get(of_FFT_complex,i));
		}

		uncE = uncE/normalizationFactor;
		//uncE=uncE/1.61403237679;
		//uncE=uncE*1000;
		//cout<<"uncE_Freq="<<uncE<<endl;

		gsl_vector_complex_free(pulseFFTcomp);
		gsl_vector_free(pulse_FFT);
	}*/

	// Inverse FFT (to get the expression of the optimal filter in time domain)
	// Complex OptimalFilter(f) => Taking into account magnitude (MatchedFilter(f)/N^2(f)) and phase (given by MatchedFilter(f))
	gsl_vector_complex *of_FFTcomp = gsl_vector_complex_alloc(mf_size);
	*optimal_filtergsl = gsl_vector_alloc(mf_size);
	for (int i=0;i<mf_size;i++)
	{
		gsl_vector_complex_set(of_FFTcomp,i,gsl_complex_polar(gsl_complex_abs(gsl_vector_complex_get(of_FFT_complex,i)),gsl_vector_get(mf_arg,i)));

	}
	if (FFTinverse(of_FFTcomp,*optimal_filtergsl,SelectedTimeDuration))
	{
	    message = "Cannot run routine FFTinverse to get optimal filter in time domain";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Calculus of the pseudoenergy in time domain
	/*if ((mode == 1) && (TorF == 0))
	{
		uncE = 0.0;
		for (int i=0;i<mf_f->size;i++)
		{
			uncE = uncE + gsl_vector_get(pulse,i)*gsl_vector_get(*optimal_filtergsl,i);
			cout<<i<<" "<<gsl_vector_get(pulse,i)<<" "<<gsl_vector_get(*optimal_filtergsl,i)<<" "<<uncE<<endl;
		}
		cout<<"uncE: "<<uncE<<endl;
		uncE = uncE/normalizationFactor;
		uncE=uncE*2*SelectedTimeDuration;
	}*/

	gsl_vector_complex_free(of_FFTcomp);
	gsl_vector_complex_free(of_FFT_complex);
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
	int status = EPOK;
	string message = "";

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
* getMatchedFilter: This function ...
*
******************************************************************************/
int getMatchedFilter(int runF0orB0val_getMatchedFilter, double pheight, gsl_vector *estenergylibrary_getMatchedFilter, gsl_matrix *matchedfilters_getMatchedFilter, gsl_matrix *matchedfiltersb0_getMatchedFilter, gsl_vector **matchedfilter_getMatchedFilter)
{
	string message = "";

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Provisional => To be deleted in future
	FILE * temporalFile;
	char temporalFileName[256];
	sprintf(temporalFileName,"auxfile");
	strcat(temporalFileName,".txt");
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// Obtain the matched filter by interpolating
	if (runF0orB0val_getMatchedFilter == 0)
	{
		if (find_matchedfilter(pheight, estenergylibrary_getMatchedFilter, matchedfilters_getMatchedFilter, matchedfilter_getMatchedFilter, temporalFile))
		{
			message = "Cannot run routine find_matchedfilter for filter interpolation";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}
	else if (runF0orB0val_getMatchedFilter == 1)
	{
		if (find_matchedfilter(pheight, estenergylibrary_getMatchedFilter, matchedfiltersb0_getMatchedFilter, matchedfilter_getMatchedFilter, temporalFile))
		{
			message = "Cannot run routine find_matchedfilter for filter interpolation";
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* find_matchedfilter:
****************************************/
int find_matchedfilter(double ph, gsl_vector *energiesvalues, gsl_matrix *matchedfilters, gsl_vector **matchedfilterFound, FILE * temporalFile)
{
	int status = EPOK;
	string message = "";
	long nummodels = energiesvalues->size;

	if (ph < gsl_vector_get(energiesvalues,0))
	{
		gsl_matrix_get_row(*matchedfilterFound,matchedfilters,0);
	}
	else if (ph > gsl_vector_get(energiesvalues,nummodels-1))
	{
		gsl_matrix_get_row(*matchedfilterFound,matchedfilters,nummodels-1);
	}
	else
	{
		for (int i=0;i<nummodels;i++)
		{
			if (ph == gsl_vector_get(energiesvalues,i))
			{
				gsl_matrix_get_row(*matchedfilterFound,matchedfilters,i);

				break;
			}
			else if ((ph > gsl_vector_get(energiesvalues,i)) && (ph < gsl_vector_get(energiesvalues,i+1)))
			{
				// Interpolate between the two corresponding rows in "models"
				gsl_vector *matchedfilterAux = gsl_vector_alloc(matchedfilters->size2);
				gsl_vector_set_zero(matchedfilterAux);
				gsl_vector *model1 = gsl_vector_alloc(matchedfilters->size2);
				gsl_vector *model2 = gsl_vector_alloc(matchedfilters->size2);
				gsl_matrix_get_row(model1,matchedfilters,i);
				gsl_matrix_get_row(model2,matchedfilters,i+1);

				if (interpolate_model(&matchedfilterAux,ph,model1,gsl_vector_get(energiesvalues,i), model2,gsl_vector_get(energiesvalues,i+1), temporalFile))
				{
				    message = "Cannot run interpolate_model routine for model interpolation";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				*matchedfilterFound = gsl_vector_alloc(matchedfilterAux->size);
				gsl_vector_memcpy(*matchedfilterFound,matchedfilterAux);

				gsl_vector_free(matchedfilterAux);
				gsl_vector_free(model1);
				gsl_vector_free(model2);

				break;
			}
		}
	}

    return(EPOK);
}
/*xxxx end of SECTION Bx xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* calculateUCEnergy function:
****************************************************************************/
int calculateUCEnergy (gsl_vector *vector, gsl_vector *filter, int domain, double nrmfctr, double samprate, double *calculatedEnergy)
{
	string message = "";

	double SelectedTimeDuration = vector->size/samprate;
	*calculatedEnergy = 0.0;

	if (domain == 0)	// Time domain filtering
	{
		if (vector->size != filter->size) *calculatedEnergy = 0.0;
		else
		{
			for (int i=0; i<vector->size; i++)
			{
				*calculatedEnergy = *calculatedEnergy + gsl_vector_get(vector,i)*gsl_vector_get(filter,i);
			}
			*calculatedEnergy = *calculatedEnergy/nrmfctr;
			*calculatedEnergy = *calculatedEnergy*2*SelectedTimeDuration;
		}
	}
	else if (domain == 1)	// Frequency domain filtering (multiply vectorFFT and filterFFT)
    {
		if (vector->size != filter->size) *calculatedEnergy = 0.0;
		else
		{
			// Declare variables
			gsl_vector_complex *vectorFFT = gsl_vector_complex_alloc(vector->size);
			gsl_vector_complex *filterFFT = gsl_vector_complex_alloc(filter->size);

			if (FFT(vector,vectorFFT,SelectedTimeDuration))
			{
				message = "Cannot run routine FFT";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
			if (FFT(filter,filterFFT,SelectedTimeDuration))
			{
				message = "Cannot run routine FFT when domain=1)";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
			*calculatedEnergy = 0.0;
			for (int i=0;i<vector->size;i++)
			{
				*calculatedEnergy = *calculatedEnergy + gsl_complex_abs(gsl_vector_complex_get(vectorFFT,i))*gsl_complex_abs(gsl_vector_complex_get(filterFFT,i));
			}
			*calculatedEnergy = *calculatedEnergy/nrmfctr;

			gsl_vector_complex_free(vectorFFT);
			gsl_vector_complex_free(filterFFT);
		}
    }

    return EPOK;
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* writeFilter: This function ...
*
******************************************************************************/
int writeFilter(ReconstructInitSIRENA *reconstruct_init_writeFilter, double normalizationFactor_writeFilter, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject_writeUCEnergy, const char *create_writeFilter)
{
	string message = "";
	int status = EPOK;

	fitsfile *fltObject;
	char fltName[256];
	strcpy(fltName,reconstruct_init_writeFilter->filterFile);

	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	// Create _dtc file (if file already exists => check clobber)
	if (fileExists(string(fltName)) && (reconstruct_init_writeFilter->clobber == 1))
	{
		if (remove(fltName))
		{
			message = "Output filter file already exists & cannot be deleted ("+string(strerror(errno))+")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	    }
	}
	else if (fileExists(string(fltName)) && (reconstruct_init_writeFilter->clobber == 0))
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

	strcpy(extname,"FILTER");
	if (fits_create_tbl(fltObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
	{
		message = "Cannot create table " + string(extname) + " in output detect file " + string(fltName);
		EP_EXIT_ERROR(message,status);
	}

	strcpy(extname,"FILTER");
	if (fits_movnam_hdu(fltObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,create_writeFilter);
	if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(keyname,"NRMFCTR");
	if (fits_write_key(fltObject,TDOUBLE,keyname,&normalizationFactor_writeFilter,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Set PROCESS keyword
	//char str_TorF[125];			sprintf(str_TorF,"%d",TorF);
	//char str_runF0orB0[125];	sprintf(str_runF0orB0,"%d",runF0orB0val);
	char str_clobber[125];      sprintf(str_clobber,"%d",reconstruct_init_writeFilter->clobber);

	string processoutFLT (string(" runFilter") 	       + ' ' +
	string(reconstruct_init_writeFilter->record_file)  + ' ' + string(reconstruct_init_writeFilter->library_file) + ' ' +
	string(fltName)	                                   + ' ' + string(reconstruct_init_writeFilter->noise_file)   + ' ' +
	string(reconstruct_init_writeFilter->FilterDomain) + ' ' + string(reconstruct_init_writeFilter->FilterMethod) 	    + ' ' +
	string(str_clobber) + ' ' +
	string("(") + (string) create_writeFilter + string(")"));

	strcpy(keyname,"PROC0");
	if (fits_write_key_longwarn(fltObject,&status))
	{
	    message = "Cannot write long warn in output file " + string(fltName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(keyvalstr,processoutFLT.c_str());
	if (fits_write_key_longstr(fltObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
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
	strcpy(obj.unit,"--");
	if (writeFitsSimple (obj,optimalfilter))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// FREQ column
	strcpy(obj.nameCol,"Freq");
	strcpy(obj.unit,"Hz");
	if (writeFitsSimple (obj,optimalfilter_f))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// OPTIMALFF column
	strcpy(obj.nameCol,"OptimalFF");
	strcpy(obj.unit,"--");
	if (writeFitsSimple (obj,optimalfilter_FFT))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Write output keywords
	strcpy(keyname,"EVENTCNT");
	long keyvalint = optimalfilter->size;
	if (fits_write_key(fltObject,TLONG,keyname,&keyvalint,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(fltName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(keyname,"EVENTSZ");
	keyvalint = 1;
	if (fits_write_key(fltObject,TLONG,keyname,&keyvalint,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(fltName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(keyname,"NRMFCTR");
	if (fits_write_key(fltObject,TDOUBLE,keyname,&normalizationFactor_writeFilter,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(fltName);
	    EP_EXIT_ERROR(message,status);
	}

	if (fits_close_file(fltObject,&status))
	{
	    message = "Cannot close file " + string(fltName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/



/***** SECTION BX ************************************************************
* writeUCEnergy: This function ...
*
******************************************************************************/
int writeUCEnergy(ReconstructInitSIRENA **reconstruct_init_writeUCEnergy, PulsesCollection *pulsesAll_writeUCEnergy, PulsesCollection *pulsesInRecord_writeUCEnergy, int pulse_index, double uncE_writeUCEnergy, fitsfile **dtcObject_writeUCEnergy, const char *create_writeUCEnergy)
{
	string message = "";
	int status = EPOK;

	char dtcName[256];
	strcpy(dtcName,(*reconstruct_init_writeUCEnergy)->detectFile);

	char extname[20];
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	if (fits_open_file(dtcObject_writeUCEnergy,dtcName,READWRITE,&status))
	{
		message = "Cannot open output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	IOData obj;
	obj.inObject = *dtcObject_writeUCEnergy;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"PULSES");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit," ");

	// UNCE column
	obj.iniRow = pulsesAll_writeUCEnergy->ndetpulses + pulse_index + 1;
	obj.endRow = pulsesAll_writeUCEnergy->ndetpulses + pulse_index + 1;
	gsl_vector *uncE = gsl_vector_alloc(1);
	gsl_vector_set(uncE,0,uncE_writeUCEnergy);
	strcpy(obj.nameCol,"UNCE");
	if (writeFitsSimple (obj,uncE))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in PULSES";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_vector_free(uncE);

	if ((*reconstruct_init_writeUCEnergy)->clobber == 1)
	{
		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(*dtcObject_writeUCEnergy, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in file " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}

		string mod1 (string("File MODIFIED by") + ' ' +	(string) create_writeUCEnergy);

		strcpy(keyname,"MOD0");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(*dtcObject_writeUCEnergy,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(dtcName);
		    EP_EXIT_ERROR(message,status);
		}
		if (fits_update_key_longstr(*dtcObject_writeUCEnergy,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		//(*reconstruct_init_writeUCEnergy)->clobber = 2;
		//if ((*reconstruct_init_writeUCEnergy)->mode == 0) (*reconstruct_init_writeUCEnergy)->clobber = 2;
	}

	if (fits_close_file(*dtcObject_writeUCEnergy,&status))
	{
	    message = "Cannot close file " + string(dtcName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	return(EPOK);
}
/*xxxx end of SECTION BX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION BX ************************************************************
* writeFilterHDU: This function ...
*
******************************************************************************/
int writeFilterHDU(ReconstructInitSIRENA **reconstruct_init_writeFilterHDU, int pulse_index, double normalizationFactor_writeFilterHDU, double uncE_writeFilterHDU, gsl_vector *optimalfilter, gsl_vector *optimalfilter_f, gsl_vector *optimalfilter_FFT, fitsfile **dtcObject_writeFilterHDU, const char *create_writeFilterHDU)
{
	string message = "";
	int status = EPOK;

	long totalpulses = 0;

	char dtcName[256];
	strcpy(dtcName,(*reconstruct_init_writeFilterHDU)->detectFile);

	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	if (fits_open_file(dtcObject_writeFilterHDU,dtcName,READWRITE,&status))
	{
		message = "Cannot open output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (((*reconstruct_init_writeFilterHDU)->clobber == 1) && (pulse_index == 0))
	{
		strcpy(extname,"FILTER");
		if (fits_create_tbl(*dtcObject_writeFilterHDU, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
			message = "Cannot create table " + string(extname) + " in output detect file " + string(dtcName);
			EP_EXIT_ERROR(message,status);
		}
		totalpulses = 0;

		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(*dtcObject_writeFilterHDU, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}
	else
	{
		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(*dtcObject_writeFilterHDU, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		if (pulse_index == 0)
		{
			if (fits_get_num_rows(*dtcObject_writeFilterHDU,&totalpulses, &status))
			{
				message = "Cannot get number of rows in " + string(dtcName);
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}
		}
		else
		{
			if (fits_get_num_rows(*dtcObject_writeFilterHDU,&totalpulses, &status))
			{
				message = "Cannot get number of rows in " + string(dtcName);
				EP_PRINT_ERROR(message,status); return(EPFAIL);
			}
			totalpulses = totalpulses-1;
		}
	}

	gsl_matrix *optimalfilter_matrix = gsl_matrix_alloc(1,(*reconstruct_init_writeFilterHDU)->pulse_length);
	gsl_matrix *optimalfilter_f_matrix = gsl_matrix_alloc(1,(*reconstruct_init_writeFilterHDU)->pulse_length);
	gsl_matrix *optimalfilter_FFT_matrix = gsl_matrix_alloc(1,(*reconstruct_init_writeFilterHDU)->pulse_length);

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
	obj.inObject = *dtcObject_writeFilterHDU;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"FILTER");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];

	// OPTIMALF column
	obj.iniRow = pulse_index+totalpulses+1;
	obj.endRow = pulse_index+totalpulses+1;
	strcpy(obj.nameCol,"OPTIMALF");
	strcpy(obj.unit,"--");
	if (writeFitsComplex (obj,optimalfilter_matrix))
	{
	    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
	    EP_EXIT_ERROR(message,EPFAIL);
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
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_vector_free(oflength);

	// NRMFCTR column
	gsl_vector *nrmfctr = gsl_vector_alloc(1);
	gsl_vector_set(nrmfctr,0,normalizationFactor_writeFilterHDU);
	strcpy(obj.nameCol,"NRMFCTR");
	if (writeFitsSimple (obj,nrmfctr))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in FILTER";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_vector_free(nrmfctr);

	// FREQ column
	strcpy(obj.nameCol,"FREQ");
	strcpy(obj.unit,"Hz");
	if (writeFitsComplex (obj,optimalfilter_f_matrix))
	{
	    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_matrix_free(optimalfilter_f_matrix);

	// OPTIMALFF column
	strcpy(obj.nameCol,"OPTIMALFF");
	strcpy(obj.unit,"--");
	if (writeFitsComplex (obj,optimalfilter_FFT_matrix))
	{
	    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_matrix_free(optimalfilter_FFT_matrix);

	strcpy(obj.nameTable,"PULSES");

	// UNCE column
	gsl_vector *uncE = gsl_vector_alloc(1);
	gsl_vector_set(uncE,0,uncE_writeFilterHDU);
	strcpy(obj.nameCol,"UNCE");
	if (writeFitsSimple (obj,uncE))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in PULSES";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	gsl_vector_free(uncE);

	if ((*reconstruct_init_writeFilterHDU)->clobber == 1)
	{
		strcpy(extname,"PULSES");
		if (fits_movnam_hdu(*dtcObject_writeFilterHDU, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in file " + string(dtcName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}

		string mod1 (string("File MODIFIED by") + ' ' +	(string) create_writeFilterHDU);

		strcpy(keyname,"MOD0");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(*dtcObject_writeFilterHDU,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(dtcName);
		    EP_EXIT_ERROR(message,status);
		}
		if (fits_update_key_longstr(*dtcObject_writeFilterHDU,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		//(*reconstruct_init_writeFilterHDU)->clobber = 2;
		/*if ((*reconstruct_init_writeFilterHDU)->mode == 0) (*reconstruct_init_writeFilterHDU)->clobber = 2;
		else
		{
			if (pulse_index == 0) (*reconstruct_init_writeFilterHDU)->clobber = 2;
		}*/
	}

	if (fits_get_num_rows(*dtcObject_writeFilterHDU,&totalpulses, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_close_file(*dtcObject_writeFilterHDU,&status))
	{
	    message = "Cannot close file " + string(dtcName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

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
	int status = EPOK;

	if (calculateEnergy(*reconstruct_init,pulsesInRecord))
	{
		message = "Cannot run calculateEnergy in runEnergy";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	if((*reconstruct_init)->intermediate==1){
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
int loadUCEnergies(ReconstructInitSIRENA *reconstruct_init_loadUCEnergies, PulsesCollection *pulsesAll_loadUCEnergies, long *nz_loadUCEnergies, gsl_vector **zi_loadUCEnergies, double *E0z_loadUCEnergies)
{
	*E0z_loadUCEnergies = reconstruct_init_loadUCEnergies->monoenergy;

	*nz_loadUCEnergies = 0;

	gsl_vector *zi_loadUCEnergiesAUX = gsl_vector_alloc(pulsesAll_loadUCEnergies->ndetpulses);
	gsl_vector_set_zero(zi_loadUCEnergiesAUX);

	for (int i=0;i<pulsesAll_loadUCEnergies->ndetpulses;i++)
	{
		if (pulsesAll_loadUCEnergies->pulses_detected[i].quality == 0)
		{
			gsl_vector_set(zi_loadUCEnergiesAUX, *nz_loadUCEnergies, pulsesAll_loadUCEnergies->pulses_detected[i].ucenergy);
			*nz_loadUCEnergies = *nz_loadUCEnergies + 1;
		}
	}

	gsl_vector_view temp;
	*zi_loadUCEnergies = gsl_vector_alloc(*nz_loadUCEnergies);
	gsl_vector_subvector(zi_loadUCEnergiesAUX,0,*nz_loadUCEnergies);
	gsl_vector_memcpy(*zi_loadUCEnergies,&temp.vector);

	gsl_vector_free(zi_loadUCEnergiesAUX);

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

	//cout<<"b_cF = "<<*b_cF<<endl;
	//cout<<"c_cF = "<<*c_cF<<endl;
	return (EPOK);
}
/*xxxx end of SECTION CX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION CX ************************************************************
* calculateEnergy: This function...
****************************************************************************/
int calculateEnergy(ReconstructInitSIRENA *reconstruct_init_calculateEnergy, PulsesCollection **pulses)
{
	double ucenergy;
	double b = reconstruct_init_calculateEnergy->b_cF;
	double c = reconstruct_init_calculateEnergy->c_cF;
	int calibLQ = reconstruct_init_calculateEnergy->calibLQ;

	for (int i=0;i<(*pulses)->ndetpulses;i++)
	{
		ucenergy = (*pulses)->pulses_detected[i].ucenergy;

		if (calibLQ == 1)	// Linear
		{
			(*pulses)->pulses_detected[i].energy = b*ucenergy;
		}
		else if (calibLQ == 2)	// Quadratic
		{
			(*pulses)->pulses_detected[i].energy = b*ucenergy + c*pow(ucenergy,2.0);
		}
	}

	return(EPOK);
}
/*xxxx end of SECTION CX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION CX ************************************************************
* writeEnergy: This function...
****************************************************************************/
int writeEnergy(ReconstructInitSIRENA **reconstruct_init_writeEnergy, PulsesCollection *pulsesInRecord_writeEnergy, const char *create_writeEnergy)
{
	string message = "";
	int status = EPOK;

	long totalpulses = 0;

	fitsfile *dtcObject_writeEnergy;
	char dtcName[256];
	strcpy(dtcName,(*reconstruct_init_writeEnergy)->detectFile);

	char *tt[1];
	char *tf[1];
	char *tu[1];
	char extname[20];
	int extver = 0;
	char keyname[10];
	char keyvalstr[1000];
	char *comment=NULL;

	if (fits_open_file(&dtcObject_writeEnergy,dtcName,READWRITE,&status))
	{
		message = "Cannot open output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	strcpy(extname,"PULSES");
	if (fits_movnam_hdu(dtcObject_writeEnergy, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output detect file " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	if (fits_get_num_rows(dtcObject_writeEnergy,&totalpulses, &status))
	{
		message = "Cannot get number of rows in " + string(dtcName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	IOData obj;
	obj.inObject = dtcObject_writeEnergy;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"PULSES");
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	for (int i=0;i<pulsesInRecord_writeEnergy->ndetpulses;i++)
	{

		// ENERGY column
		obj.iniRow = totalpulses + i + 1 -pulsesInRecord_writeEnergy->ndetpulses;
		obj.endRow = totalpulses + i + 1 -pulsesInRecord_writeEnergy->ndetpulses;
		strcpy(obj.nameCol,"ENERGY");
		strcpy(obj.unit,"eV");
		gsl_vector *energy = gsl_vector_alloc(1);
		gsl_vector_set(energy,0,pulsesInRecord_writeEnergy->pulses_detected[i].energy);
		strcpy(obj.nameCol,"ENERGY");
		if (writeFitsSimple (obj,energy))
		{
		    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) +
		    		" column in FILTER";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
		gsl_vector_free(energy);
	}

	if ((*reconstruct_init_writeEnergy)->clobber == 1)
	{
		string mod1 (string("File MODIFIED by") + ' ' +	(string) create_writeEnergy);

		strcpy(keyname,"MOD1");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(dtcObject_writeEnergy,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(dtcName);
		    EP_EXIT_ERROR(message,status);
		}
		if (fits_update_key_longstr(dtcObject_writeEnergy,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(dtcName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		(*reconstruct_init_writeEnergy)->clobber = 2;
	}

	if (fits_close_file(dtcObject_writeEnergy,&status))
	{
	    message = "Cannot close file " + string(dtcName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	return(EPOK);
}
/*xxxx end of SECTION CX xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
