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

   Copyright 2014:  Trigger has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

/***********************************************************************
*                      PULSEGRADE
*
*  File:      pulsegrade.cpp
*  Version:   13.0.0
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene Gonzélez Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
*	version 1.0.0 	10/10/06	First version
*  	version 1.0.1  	08/03/07	Change the name of Group Extension of input Fits file
* 							   	for "EUR-TRG"
* 							   	Change the name of Group Extension of output Fits file
* 							   	for "EUR-PSH" and "EUR-PULSESHAPE"
*  	version 1.0.2	23/04/07 	Delete environment value "FITS_DIR", it's incompatible with the browse of param_gui.
* 							   	"FIT_DIR" imposes a input directory which contains the input fits file (That is incompatible)
* 							   	with the browse).
*  								Change input Name of the fits file. Now, the Name include the .fits extension. Now, the file
* 							   	Name is compatible with the browse of param_gui.
* 								Change in the Call to "readFitsComplex", "readFitsSimple", "writeFitsComplex" and
* 								"writeFitsSimple". Now is used the struct "OIData". See documentation of Utils module. (inoutUtils).
*  	version 1.0.3	24/04/07 	Include new input parameter in the module. Adding in .par file "output File". Now,
* 								the user can choose the name of the output .fits file.
*								-jjs-
* 								Change maximum number of data channels in an input file (4)
*  								Parameterise the names for the input ADC channels, and
*          					 	the names for the result data column names.
*  								If inputfile is NOT modulated, do not demodulate, but do all else, eg. scaling the data
* 								ADC channel names are now parameterised
* 	version	1.1.0	23/05/07 	Modification of the module using the new module inoutFits included in Utils. Now, Reads and Writes
* 								of Fits use block of rows.
* 								Including coments using RIL logs.
* 	version 1.1.1 	29/05/07	Free allocate of memory
*
*  	version 2.0.0   07/04/08	Pulse Shape reads from trg fits file.
* 								Reading column "Quality" of input fits file
* 								Changing the Group Extension "ADC_RAW" for "EUR-PSH" of the output fits file.
* 								Filtering of pulse using the column "Quality"
* 								Adding the output column "quality"
* 	version 2.1.0	26/06/08	Adding new input keyword "NBQUAL"
* 								Adding findMultipulse function
* 								Adding transformDatas function
* 								Adding fitNorris function
* 								Adding input/output keyword "ql"
* 								Adding Column "Quality" in PSH FITS file.
* 								Adding output Keyword in PSH FITS "NP-FIT"
* 								Create a HISTORY in the trigger file.
* 								Deleting input/output keyword "ql";
* 								Adding input parameter "ql"
* 								Adding baseline and fit_baseline columns in output file
* 								Adding getPulseParams function
*  								Adding input parameter "wRise" and "wFall"
* 								Adding new columns "tstart" and "tend" from output FITS file.
* 								D//+eleting column "chisq " from output FITS file
*
*  version 3.0.0	10/09/08	Changing documentacion of file
* 								Adding variable called "version". It contains the number version of this module.
* 								Adding input parametrer "nameLog"
* 								Adding input parametrer "verbosity" and processing it using RIL libraries.
* 								Adding new function "readInputKeywords"
* 								Deleting in/output keyword "channelcount"
* 								Including error codes.
*  version 3.1.0 	29/09/08	Included "CREATE" and "PROCESS" keyword in the output FITS file
* 								Include error processing of "INPUT_PARAMETER_NON_DEFINE"
*  version 3.2.0	24/10/08	If the module fail. It will be warning in the output.
* 								Included parameter time to measure.
* version 4.0.0		14/01/09 	Changed the input column name of Current from "V" to "I0"
* 								Changed the unit of wRise and wFall input parameter from bins to time
* 								Included new input keywords
* version 4.0.1		19/01/09	Changed fixed value of wFall to input parameter wFall.
* version 4.1.0		22/01/09	Included new error process:
* 									-63804 PSH_NON_FOUND_PULSES: Non found pulses. Check the input paramenters of trigger task.
* 								Modificated error code of input parameters
* version 4.2.0		02/02/09	Included new Warning: Non Valid Norris Fit and new codes of quality columns.
* 								Changed name of output keyword "CREATE" to "CREATOR"
* version 4.2.1		04/02/09	Resolved bug of norris fit (nan values).
* 								Resolved bug of baseline. the TRG pulses have baseline = 0;
* version 5.0.0		20/02/09 	Included, modificated and removed input keywords
* version 5.1.0		24/02/09	Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
* 								Solved bug 4: stops on zero value for the keyword TIMEZERO, error -62436.
* 								Included FTYPE keywords in output FITS file.
* version 5.2.0		12/03/09	Pulses with quality values non zero is not processed.
* 								Included new function called findSaturatedPulses.
* version 6.0.0						Deleted inDataIteratorBl function.
* 								Read baseline and sigma columns of EUR-TRG extension.
* 								Create baseline and sigma columns in the output FITS file.
* 								Changed the method to obtain trise and tfall
*                               Used the keywords libraries
* version 6.1.0		15/04/09    Delete "PILInit(argc, argv);" and "PILClose(PIL_OK);"
*                               The input parameter fitabstol is now used during the fitting instead a hardpoint
* 								Corrected the calculus of Trise and Tfall in getPulseParams
*                               Documentation updated
* version 7.0.0     08/06/09    Norris fit columns not written in Quick Look mode
* version 8.0.0                 Changed order in delete's
* 								Assumed Jaap's suggestions
* 								findMultipulses: Warning instead of error.
*                               getPulseParams runs if gsl_vector_get(qualoutgsl,0) == 0 (not qual == 0)
* version 8.0.1		28/09/09    Columns units corrected
* version x.x.x		xx/xx/09	If error file exists, its contests will be deleted
* 					30/11/09    The PROCESS keyword written in the modified TRG file keeps the PROCESS keyword read from the
*                               TRG input FITS file and the PROCESS generated in PULSESHAPE and the PULSESHAPE version are added
*                               The HISTORY keyword is written in the _trg FITS file, not in the _psh FITS file
*                               63804 => If there is no pulses in the _trg input FITS file
* version 9.0.0		10/09/10	Deleted Norris function fitting
*								Deleted output columns obtained from the fit: MaxTimeFit, MaxCurrentFit, TriseFit, NriseFit, TfallFit, NfallFit, ChisqFit, BaselineFit.
*								Deleted other output columns: TIME, Baseline, Sigma
*								Now, no differences between Quick and Exhaustive methods of analysis
*								Saturated pulses are now found. New function added findSaturatedPulses
*								Deleted input parameter fitabstol. Added input parameter CurrentLevel
*								Not necessary to cut most part of the pulse (tstart-Tend)
*								Deleted extension BASELINE in input FITS file
*								Input column difTstrt added. Input columns TIME, EndPulse, Baseline, Sigma deleted
*								Deleted functions fit_norris and norris
* version 9.1.0		28/09/10	"fit_linear2" function has been modified and moved to Utils library with name "fit_linear"
* version 9.2.0		30/09/10	The value of the keyword "EVENTSZ" is changed
* 					08/11/10	Adding keyword "FLT_ID" (FLT_ID=0)
* 								Documentation updated
* version 10.0.0	10/12/10	Resolved bug 49 --> Multiple pulses with the maximum at the beginning or at the end are now found
* 					16/12/10	Deleted input parameters "wRise" & "wfall"
* 								'Log' argument instead 'OK' in some "writeLog" functions
* 								Deleted output columns: "Trise", "Nrise", "Tfall" & "Nfall"
* 								Deleted functions: "getPulseParams" & "transformDatas"
* 								Not necessary to read "tauRise", "tauFall" and "difTstrt" columns from xay_t.fits
* 								Documentation updated
*                   18/02/11    findMultipulses not used (inDataIteratorTrg and pulseShape modified)
*                               In TRG there is not an I0_Filtered column => Updated the PULSESHAPE code
*                   28/02/11	PROCESS also includes PULSESHAPE version (CREATOR)
* version 11.0.0    22/03/11    "columnName" input parameter not used
*                   24/03/11    "status" allocation and "writeLog" order changed in some cases
*                   25/03/11    Quick Look mode no longer used
*                   05/04/11    No pulses in _trg input FITS file => Empty EUR-PSH in _psh output FITS file
*                   22/12/11    Event grading (new function "eventGrade" in "pulseShape") => New column GRADE in EUR-PSH
*                   13/02/12    New column GRADE in EUR-TRG
*                   06/02/13    New output keyword "MODE"
*                               Saturated pulses are no longer searched for (already done in TRIGGER):
*									'CurrentLevel' variable deleted
*									'findSaturatedPulses' function deleted
*                               Pulses are no longer classified as multipulses (redundant with the grading)
*								All pulses are going to be graded (including truncated or saturated)
*					            Deleted error handling (63801) related to non 'goodPulses'
*					15/05/13    New output keywords "B_CF" and "C_CF"
*					15/10/13    Deleted keyword "FLT_ID"
*					27/11/13    ANNALS keyword substitutes HISTORY 'pseudokeyword'
*					23/01/14    'if (verbosity<0 || verbosity>3)' instead of 'if (verbosity <= 0)' in initModule
*		            05/07/14  	Removed PIL/RIL/Common dependency
*  version 11.0.1   07/07/14    Modified some comments & parameters description
*  version 11.1.0   07/07/14    Adapted for new parameter in interativePars function for task name
*  version 12.0.0   10/07/14    Features of the pulses are no longer got => Output MaxTime and MaxCurrent columns are no longer written
*                   11/07/14    Solved bug preventing the task from reading the full command line
*  version 12.0.1   11/07/14    Comments of the 'initModule' function modified
*                   15/07/14    'initModule' modified in order to accept negative int or double (if parameters are read from command line)
*                   31/07/14    Comments referring to PIL and RIL deleted
*                   10/09/14    Differences when reading or writing some keywords
*                   24/09/14    inTaus and outTaus input parameters have been changed from 'int' to 'double'
*                   15/10/14    TailBfr column is read from _trg and if there is a tail before the pulse, its grade is fixed at 32
*  version 13.0.0     Dec/14    DAL -> CFITSIO migration (MC)
*								'pulsegrade' instead of 'pulseshape'
*								Deleted some unnecessary input/output keywords
***********************************************************************/

/******************************************************************************
DESCRIPTION:

The goal of this task is grading the pulses (also called events in the literature).

The user must supply the following input parameters:

- trgFile: Name of the input _trg FITS file which contains the found pulses
- pshFile: Name of the output FITS file
- tauFALL: Fall time of the pulses in seconds
- inTaus: Times the fall time of the pulses to establish the inner interval when grading the events
- outTaus: Times the fall time of the pulses to establish the outer interval when grading the events
- namelog: Output log file name
- verbosity: Verbosity level of the output log file

The results are written into an _psh output FITS file (RunGroupIDRunIDSubrunID_psh.fits) and also
the _trg input FITS file is modified.

MAP OF SECTIONS IN THIS FILE:

 - 1. INCLUDE's
 - 2. MAIN
 - 3. initModule
 - 4. createPulseShapeFile
 - 5. inDataIteratorTrg
 - 6. pulseShape
 - 7. eventGrade

*******************************************************************************/

/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include "pulsegrade.h"

/***** SECTION 2 ************************************
*  MAIN function: This function is the main function of the PULSEGRADE task
*
* - Read input parameters
* - Open _trg input FITS file
* - Read input keywords and check their values
* - If there are no pulses in _trg input FITS file (empty EUR-TRG extension)
* 	- Create output FITS files (_psh.fits file) (call createPulseShapeFile)
*	- Create ANNALS in the TRIGGER file
*	- Close input and output FITS files
*	- Finalize the task
* - There are pulses in _trg input FITS file
* 	- Get structure of input FITS file columns
*   - Create output FITS files (_psh.fits file) (call createPulseShapeFile)
*   - Create structure to run Iteration: inDataIteratorTrg
*		- Read columns (Tstart, I0, Tend, Quality, TailBfr)
*		- Called iteration function: inDataIteratorTrg (which it is going to call pulseShape)
* 		- Create an ANNALS in the TRIGGER file
* 		- Write output keywords
*		- Close input and output FITS files
* 		- Finalize the task
****************************************************/
int main (int argc, char **argv)
{
	create = "pulsegrade v.13.0.0";			//Set "CREATOR" keyword of output FITS file
	time_t t_start = time(NULL);
	string message;
	int status=EPOK, extver=0, one=1;
	
	sprintf(temporalFileName,"PULSEGRADEauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if (temporalFile == NULL)
	{
	  message = "Cannot open auxiliary file PULSEGRADEauxfile.txt";
	  EP_EXIT_ERROR(message,EPFAIL);
	}

	// Read input parameters
	if (initModule(argc, argv))
	{
	  message = "Error in initModule";
	  EP_EXIT_ERROR(message,EPFAIL);
	}
	
	writeLog(fileRef,"Log", verbosity,"Into Pulsegrade task");

	// Open _trg input FITS file
	if (fits_open_file(&trgObject, trgName,READWRITE,&status))
	{
	    message = "Cannot open file " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(extname,"EUR-TRG");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(trgName) + " (before keyword reading)";
	    EP_EXIT_ERROR(message,status);
	}

	// Read	input keywords and check their values
	strcpy(keyname,"EVENTCNT");
	if (fits_read_key(trgObject,TLONG,keyname, &eventcnt,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	if (eventcnt < 0)
	{
	    message = "Legal values for EVENTCNT (EUR-TRG) are non negative integer numbers";
	    writeLog(fileRef, "Error", verbosity, message);
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	strcpy(keyname,"EVENTSZ");
	if (fits_read_key(trgObject,TLONG,keyname, &eventsz,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	if (eventsz < 0)
	{
	    message = "Legal values for EVENTSZ (EUR-TRG) are non negative integer numbers";
	    writeLog(fileRef, "Error", verbosity, message);
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	strcpy(keyname,"PROCESS");
	if (fits_read_key_longstr(trgObject,keyname,&processin,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	strcpy(keyname,"ENERGY");
	if(fits_read_key(trgObject,TDOUBLE,keyname, &energy,comment,&status)){
	    message = "Cannot read key " + string(keyname) + "in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	
	if ((eventcnt == 0) && (eventsz == 0))	// There are no pulses in the _trg input FITS file
	{
		// Create output FITS file
		if (createPulseShapeFile())		// Pulsegrade file (*_psh.fits)
		{
		    message = "Error creating pulsegrade file ";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		strcpy(keyname,"EVENTCNT");
		if (fits_update_key(trgObject,TLONG,keyname, 0,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(trgName);
		    EP_EXIT_ERROR(message,status);	  
		}
		strcpy(keyname,"EVENTSZ");
		if (fits_update_key(trgObject,TLONG,keyname, 0,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(trgName);
		    EP_EXIT_ERROR(message,status);	  
		}
		string annals (string("File MODIFIED by") + ' ' +	(string) create);
		strcpy(keyname,"ANNALS");
		strcpy(keyvalstr,annals.c_str());
		if (fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(trgName);
		    EP_EXIT_ERROR(message,status);	  
		}

		if (fits_close_file(trgObject,&status) || fits_close_file(pshObject,&status))
		{
	  	    message = "Cannot close file " + string(trgName);
		    EP_EXIT_ERROR(message,status);	  
		}

		time_t t_end = time(NULL);
		message = "There are no pulses in the " + string(trgName) + " input FITS file ";
		writeLog(fileRef,"Warning",verbosity,message);
		message = "The FITS file " + string(trgName) + " will have an empty extension EUR-PSH";
		writeLog(fileRef,"Warning",verbosity,message);
		writeLog(fileRef,"OK", verbosity,"Pulsegrade Module OK");
		sprintf(straux,"%f",(double) (t_end - t_start));
		message = "Time: " + string(straux);
		writeLog(fileRef,"Log", verbosity,message);

		if (fclose(fileRef))
		{
	  	    message = "Cannot close log file ";
		    EP_EXIT_ERROR(message,EPFAIL);	  
		}	
	}

	// There are pulses in the _trg input FITS file
	// Get structure of input FITS file columns
	strcpy(extname,"EUR-TRG");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in " + string(trgName) + " (before column access)";
	    EP_EXIT_ERROR(message,status);	  
	}
	strcpy(straux,"TSTART");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
  	    message = "Cannot access column " + string(straux) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	strcpy(straux,"I0");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
        message = "Cannot access column " + string(straux) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	strcpy(straux,"TEND");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
        message = "Cannot access column " + string(straux) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	strcpy(straux,"QUALITY");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot access column " + string(straux) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	strcpy(straux,"TAILBFR");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot access column " + string(straux) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}

	message = "Open TriggerFits: " + string(trgName);
	writeLog(fileRef,"Log", verbosity,message);

	// Create output FITS File
	if (createPulseShapeFile())		// Pulsegrade file (*_psh.fits)
	{
	    message = "Error creating Pulsegrade file ";
	    EP_EXIT_ERROR(message,EPFAIL);	  
	}
	if (fits_open_file(&pshObject,pshName,READWRITE,&status))
	{
	    message = "Cannot open file " + string(pshName);
	    EP_EXIT_ERROR(message,status);	  
	}

	extern int inDataIteratorTrg(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols,void *user_strct );

	// Create structure to run Iteration
	iteratorCol cols [5];			// Structure of Iteration
	int n_cols = 5; 				// Number of columns: TSTART + I0 + TEND + QUALITY + TAILBFR
	long rows_per_loop = 0; 		// 0: Use default: Optimum number of rows
	long offset=0; 					// 0: Process all the rows

	strcpy(extname,"EUR-TRG");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in " + string(trgName) + " (before iteration process)";
	    EP_EXIT_ERROR(message,status);
	}
		// Read Tstart
	strcpy(straux,"Tstart");
	status = fits_iter_set_by_name(&cols[0], trgObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate column " + string(straux) + " in file " + string(trgName);
	    EP_EXIT_ERROR(message,status);	
	}
		// Read I
	strcpy(straux,"I0");
	status = fits_iter_set_by_name(&cols[1], trgObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate column " + string(straux) + " in file " + string(trgName);
	    EP_EXIT_ERROR(message,status);	
	}
		// Read Other Columns
	strcpy(straux,"Tend");
	status = fits_iter_set_by_name(&cols[2], trgObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate column " + string(straux) + " in file " + string(trgName);
	    EP_EXIT_ERROR(message,status);	
	}
	strcpy(straux,"Quality");
	status = fits_iter_set_by_name(&cols[3], trgObject, straux, TSHORT, InputCol);
	if (status)
	{
	    message = "Cannot iterate column " + string(straux) + " in file " + string(trgName);
	    EP_EXIT_ERROR(message,status);	
	}
	strcpy(straux,"TailBfr");
	status = fits_iter_set_by_name(&cols[4], trgObject, straux, TSHORT, InputCol);
	if (status)
	{
	    message = "Cannot iterate column " + string(straux) + " in file " + string(trgName);
	    EP_EXIT_ERROR(message,status);	
	}

	// Called iteration function
	pulseF = gsl_vector_alloc(eventsz);
	pulseG = gsl_vector_alloc(eventsz);
	if (fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIteratorTrg,0L,&status))
	{
	    message = "Cannot iterate data through inDataIteratorTrg ";
	    EP_EXIT_ERROR(message,status);	
	}
	gsl_vector_free(pulseF);
	gsl_vector_free(pulseG);

	// Create ANNALS in trigger file
	string annals (string("File MODIFIED by") + ' ' +	(string) create);
	strcpy(keyname,"ANNALS");
	strcpy(keyvalstr,annals.c_str());
	if (fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(trgName);
	    EP_EXIT_ERROR(message,status);	
	}

	// Write output keywords
	strcpy(keyname,"EVENTSZ");
	if (fits_write_key(pshObject,TINT,keyname,&one,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in file " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}

	// Close input and output FITS files
	if (fits_close_file(trgObject,&status) || fits_close_file(pshObject,&status))
	{
	    message = "Cannot close file(s) " + string(trgName) + "/" + string(pshName);
	    EP_EXIT_ERROR(message,status);
	}
	if (fclose(temporalFile))
	{
	    message = "Cannot close temp file " + string(temporalFileName);
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Finalize the task
	time_t t_end = time(NULL);
	sprintf(straux,"%f",(double) (t_end - t_start));
	message = "Time:" + string(straux);
	writeLog(fileRef,"Log", verbosity,message);
	writeLog(fileRef,"OK", verbosity,"Pulsegrade Module OK");
	
	if (fclose(fileRef))
	{
	    message = "Cannot close log file ";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	
	return EPOK;
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
* initModule function: This function reads and processes the input parameters
*
* - It defines the input parameters and assigns the values to variables
* - It defines the structure for the command line options
* - It gets the parameter values from the command line or user input
* - It saves the parameter values into meaningful variables
*   (checking also if the parameter value is in the allowed range)
* - It opens the log file of the task
****************************************************************************/
int initModule (int argc, char **argv)
{
	int status = EPOK;
	
	// Define PULSEGRADE input parameters and assign values to variables
	// Parameter definition and assignation of default values
	const int npars = 7, npars1 = 8;
	inparam pshapePars[npars];
	int optidx =0, par=0, fst=0, ipar; 
	string message, task="pulsegrade";

	pshapePars[0].name = "trgFile";
    pshapePars[0].description = " TRG Input file name";
    pshapePars[0].defValStr = "a_trg.fits";
    pshapePars[0].type =  "char";
	pshapePars[0].ValStr = pshapePars[0].defValStr;

	pshapePars[1].name = "pshFile";
    pshapePars[1].description = "Output file name";
    pshapePars[1].defValStr = "a_psh.fits";
    pshapePars[1].type =  "char";
	pshapePars[1].ValStr = pshapePars[1].defValStr;
	
	pshapePars[2].name = "tauFALL"; 
    pshapePars[2].description = "Fall time of the pulses (seconds, [>0])";
    pshapePars[2].defValReal = 3.5E-4;
    pshapePars[2].type = "double";
    pshapePars[2].minValReal = 1.E-10;
    pshapePars[2].maxValReal = 1.E+10;
	pshapePars[2].ValReal = pshapePars[2].defValReal;

	pshapePars[3].name = "inTaus";
    pshapePars[3].description = "Times the fall time of the pulses to establish the inner interval for event grading [>0]";
    pshapePars[3].defValReal = 10.;
    pshapePars[3].type = "double";
    pshapePars[3].minValReal = 1.E-10;
    pshapePars[3].maxValReal = 1.E+10;
	pshapePars[3].ValReal = pshapePars[3].defValReal;

	pshapePars[4].name = "outTaus";
    pshapePars[4].description = "Times the fall time of the pulses to establish the outer interval for event grading [>0, >inTaus]";
    pshapePars[4].defValReal = 40.;
    pshapePars[4].type = "double";
    pshapePars[4].minValReal = 1.E-10;
    pshapePars[4].maxValReal = 1.E+10;	
	pshapePars[4].ValReal = pshapePars[4].defValReal;

	pshapePars[5].name = "nameLog"; 
    pshapePars[5].description = "Output log file Name";
    pshapePars[5].defValStr = "psh_log.txt";
    pshapePars[5].type = "char";
	pshapePars[5].ValStr = pshapePars[5].defValStr;
	
	pshapePars[6].name = "verbosity"; 
    pshapePars[6].description = "Verbosity Level of the output log file ( [0,3] )";
    pshapePars[6].defValInt = 3;
    pshapePars[6].type = "int";
    pshapePars[6].minValInt = 0;
    pshapePars[6].maxValInt = 3;
	pshapePars[6].ValInt = pshapePars[6].defValInt;
	
	// Define structure for command line options
	static struct option long_options[npars1];
	for (optidx = 0; optidx<npars; optidx++)
	{
		long_options[optidx].name= pshapePars[optidx].name.c_str();
		long_options[optidx].has_arg = 2;
		long_options[optidx].flag = 0;
		long_options[optidx].val = 0;
	}
	long_options[npars].name=0;
	long_options[npars].has_arg = 0;
	long_options[npars].flag = 0;
	long_options[npars].val = 0;

	// Get parameter values from command line or user input
	optidx = 0;
	int commandLine = 0; // initialize command line as being empty
	while ((par=getopt_long(argc, argv,"",long_options, &optidx ))!= -1) // while there are still parameters in command line
	{
		// Params have been specified in the command line
		commandLine = 1;
		switch(par)
		{
			case 0:
				// identify parameter read from command line
				for(int i=0;i<npars; i++)
				{
					if (long_options[optidx].name == pshapePars[i].name.c_str())
					{
						if (pshapePars[i].type == "char") //save char value for par
						{
							pshapePars[i].ValStr = optarg;
						}
						else // check if numeric value
						{
							if ((!isdigit(optarg[0]) && (optarg[0] != '-')) ||
									(!isdigit(optarg[0]) && (optarg[0] == '-') && (!isdigit(optarg[1]))))
							{
								message = "Invalid value for input argument " + string(pshapePars[i].name);
								EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
							}
							if (pshapePars[i].type == "int")
							{
								pshapePars[i].ValInt = atoi(optarg);
							}
							else
							{
								pshapePars[i].ValReal= atof(optarg);
							}
						}
						break;
					} // endif
				} // endfor
				break;
			default:
				message = "Invalid parameter name " + string(long_options[optidx].name);
				EP_PRINT_ERROR(message,EPFAIL);
		}//switch
	}//while
	// If command line is empty: ask for params interactively
	if(commandLine == 0)
	{
		if (interactivePars(pshapePars,npars,task))
		{
		  message = "Error reading parameters interactively";
		  EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
	}

	// Save parameter values into meaningful variables
	for(int i=0;i<npars; i++)
	{
		if(pshapePars[i].name == "trgFile")
		{
			strcpy(trgName, pshapePars[i].ValStr.c_str());
		}
		else if(pshapePars[i].name == "pshFile")
		{
			strcpy(pshName, pshapePars[i].ValStr.c_str());
		}
		else if(pshapePars[i].name == "tauFALL")
		{
			tauFALL = pshapePars[i].ValReal;
		}
		else if(pshapePars[i].name == "inTaus")
		{
			inTaus = pshapePars[i].ValReal;
		}
		else if(pshapePars[i].name == "outTaus")
		{
			outTaus = pshapePars[i].ValReal;
		}
		else if(pshapePars[i].name == "nameLog")
		{
			strcpy(nameLog,pshapePars[i].ValStr.c_str());
		}
		else if(pshapePars[i].name == "verbosity")
		{
			verbosity = pshapePars[i].ValInt;
		}

		// Check if parameter value is in allowed range
		if (pshapePars[i].type == "int" &&
	        (pshapePars[i].ValInt < pshapePars[i].minValInt || 
	         pshapePars[i].ValInt > pshapePars[i].maxValInt))
		{
	    
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(pshapePars[i].minValInt) + "," + boost::lexical_cast<std::string>(pshapePars[i].maxValInt) + "]";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		else if (pshapePars[i].type == "double" &&
	        (pshapePars[i].ValReal < pshapePars[i].minValReal || 
	        		pshapePars[i].ValReal > pshapePars[i].maxValReal))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(pshapePars[i].minValReal) + "," + boost::lexical_cast<std::string>(pshapePars[i].maxValReal) + "]";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	// Open the log file of the task
	fileRef = fopen(nameLog,"w+");	// Remove file if it already exists and open a new file to save log messages
	if (fileRef == NULL)
	{
	  message = "Cannot open log file";
	  EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	return (EPOK);
}	
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* createPulseShapeFile function: This function creates the structure of the _psh output FITS file
*
* - Create _psh output FITS file (if it does not exist yet)
* - Create extension EUR-PSH
* - Create keywords
***************************************************************************/
int createPulseShapeFile ()
{
	string message;
	int status = EPOK, extver=0, one=1;
	
	message = "Into createPulseGradeFile";
	writeLog(fileRef,"Log", verbosity,message);

	// Create output FITS file: If it does not exist yet
	// Create Pulsegrade file
	if (fits_open_file(&pshObject,pshName,READWRITE,&status) == EPOK)
	{
  	    message = "Input file (" + string(pshName) + ") already exists, it cannot be overwritten ";
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	status = EPOK;
	if (fits_create_file(&pshObject, pshName, &status))
	{
	    message = "Cannot create file " + string(pshName);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	
	message = "Create Pulse Shape Fits File: " + string(pshName);
	writeLog(fileRef,"Log", verbosity,message);

	// Create extension: EUR-PSH
	char *tt[1];
	char *tf[1];
	char *tu[1];
	strcpy(extname,"EUR-PSH");
	if (fits_open_file(&pshObject,pshName,READWRITE,&status))
	{
		message = "Cannot open input file " + string(pshName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	if (fits_create_tbl(pshObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
	{
		message = "Cannot crate table " + string(extname) + " in file " + string(pshName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	if (fits_movnam_hdu(pshObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in file " + string(pshName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	// Create keywords
	strcpy(keyname,"EVENTCNT");
	if (eventcnt < 0)
	{
		message = "Legal values for EVENTCNT (EUR-PSH) are integer numbers greater than or equal to 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	if (fits_write_key(pshObject,TLONG,keyname,&eventcnt,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(pshName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,create);
	if (fits_write_key(pshObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(pshName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(keyname,"ENERGY");
	if(fits_write_key(pshObject,TDOUBLE,keyname,&energy,comment,&status)){
	      message = "Cannot write keyword " + string(keyname) + "in " + string(pshName);
	      EP_PRINT_ERROR(message,status); return(EPFAIL);  
	}
	
	// Set PROCESS keyword
	char str_verb[125];	sprintf(str_verb,"%d",verbosity);
	char str_tauFALL[125];	sprintf(str_tauFALL,"%f",tauFALL);
	char str_inTaus[125];	sprintf(str_inTaus,"%f",inTaus);
	char str_outTaus[125];	sprintf(str_outTaus,"%f",outTaus);

	string processoutPSH (string("PULSEGRADE") 	+ ' ' +
	string(trgName)				+ ' ' + string(pshName)	    + ' ' +
	string(str_tauFALL)	        + ' ' + string(str_inTaus)	+ ' ' + string(str_outTaus)	+ ' ' +
	string(nameLog)             + ' ' + string(str_verb)	+ ' ' +
	string("(")				    +      (string) create 	+   	string(")"));

	string processoutTRG ((string) processin + ' ' + processoutPSH);

	strcpy(keyname,"PROCESS");
	strcpy(keyvalstr,processoutPSH.c_str());
	if (fits_write_key_longwarn(pshObject,&status))
	{
		message = "Cannot write long warn in output file " + string(pshName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	if (fits_update_key_longstr(pshObject,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(pshName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(extname,"EUR-TRG");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in file " + string(trgName);
		EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(keyname,"PROCESS");
	strcpy(keyvalstr,processoutTRG.c_str());
	
	if (fits_update_key_longstr(trgObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot update keyword " + string(keyname) + " in " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	return EPOK;
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* inDataIteratorTrg: This function takes the optimum number of rows to read the input FITS file
*                 and works iteratively
*
* - Declaration of variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each found pulse
* 	- Call pulseShape
* - Free allocate of GSL vectors
*****************************************************************************/
int inDataIteratorTrg (long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct)
{
	char val[256];
	char val_aux[256];

 	int status=EPOK;
	string message="";

	// Declaration of variables
	double *v, *vin;			// Vector of V column
	double *tstartv,*tstartvin;		// Vector of Tstart column
	double *tend,*tendin;			// Vector of Tend column
	int *qualityP,*qualityPin;		// Vector of QUALITY column
	int *tailbefore, *tailbeforein;
	
	// Info about a pulse (i) and the next one (nexti) and the previous one (previ) is necessary to event grade
	gsl_vector *pulseprevi = gsl_vector_alloc(eventsz);
	gsl_vector *pulsei = gsl_vector_alloc(eventsz);
	gsl_vector *pulsenexti = gsl_vector_alloc(eventsz);
	double tstartprevi,tstarti, tstartnexti;
	double tendprevi, tendi, tendnexti;
	double qualityprevi,qualityi, qualitynexti;
	double tailbeforeprevi,tailbeforei, tailbeforenexti;

	// Allocate input GSL vectors
	gsl_matrix *pulses = gsl_matrix_alloc(nrows,eventsz);				// GSL of input pulse column
	gsl_vector *tstartgsl = gsl_vector_alloc(nrows);
	gsl_vector *tendgsl  = gsl_vector_alloc(nrows);
	gsl_vector *pulse = gsl_vector_alloc(eventsz);						// GSL of each row (event)
	gsl_vector *qualitygsl = gsl_vector_alloc(nrows);					// GSL of input QUALITY column
	gsl_vector *tailbeforegsl = gsl_vector_alloc(nrows);

	// Read iterator
	// NOTE: fits_iter_get_array: in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	tstartvin = (double *) fits_iter_get_array(&cols[0]);
	tstartv = &tstartvin[1];
	if (toGslVector(((void **)&tstartv),&tstartgsl,nrows,0,TDOUBLE))
	{
	    message = "Error converting Tstart column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	vin = (double *) fits_iter_get_array(&cols[1]);
	v = &vin[1];
	if (toGslMatrix(((void **)&v),&pulses,eventsz,nrows,(int)TDOUBLE,0))
	{
	    message = "Error converting I0 column to GSL matrix";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}	  
	tendin = (double *) fits_iter_get_array(&cols[2]);
	tend = &tendin[1];
	if (toGslVector(((void **)&tend),&tendgsl,nrows,0,TDOUBLE))
	{
	    message = "Error converting Tend column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	qualityPin = (int *) fits_iter_get_array(&cols[3]);
	qualityP = &qualityPin[1];
	if (toGslVector(((void **)&qualityP),&qualitygsl,nrows,0,TSHORT))
	{
	    message = "Error converting Quality column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	tailbeforein = (int *) fits_iter_get_array(&cols[4]);
	tailbefore = &tailbeforein[1];
	if (toGslVector(((void **)&tailbefore),&tailbeforegsl,nrows,0,TSHORT))
	{
	    message = "Error converting TailBfr column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	//  Processing each found pulse
	for (int i=0; i< nrows; i++){

		sprintf(val,"-----------> Pulse: ");
		sprintf(val_aux,"%d",ntotalrows);
		strcat(val,val_aux);
		sprintf(val_aux," of ");
		strcat(val,val_aux);
		sprintf(val_aux,"%d",eventcnt);
		strcat(val,val_aux);
		message = string(val) + "<------------------ ";
		writeLog(fileRef,"Log", verbosity,message);
		strcat(val,"\n");
		fputs(val,temporalFile);

		if (nrows == 1)
		{
			if (totalpulses == 1)
			// First 1-inDataIteratorTrg-set
			{
				// Read pulse-prev-i and stored in *F
				gsl_matrix_get_row(pulseF,pulses,i);
				tstartF = gsl_vector_get(tstartgsl,i);
				tendF = gsl_vector_get(tendgsl,i);
				qualityF = gsl_vector_get(qualitygsl,i);
				tailbeforeF = gsl_vector_get(tailbeforegsl,i);
			}
			else if (totalpulses == 2)
			// Second 1-inDataIteratorTrg-set
			{
				// Not read pulse-i data because it is already in *F
				gsl_vector_memcpy(pulsei,pulseF);
				tstarti = tstartF;
				tendi = tendF;
				qualityi = qualityF;
				tailbeforei = tailbeforeF;
				// Read pulse-next-i
				gsl_matrix_get_row(pulsenexti,pulses,i);
				tstartnexti = gsl_vector_get(tstartgsl,i);
				tendnexti = gsl_vector_get(tendgsl,i);
				qualitynexti = gsl_vector_get(qualitygsl,i);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				// First pulse of the file (tendprevi = -1 to pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,-1.0, -1.0,tstartnexti))
				{
				    message = "Error in pulseShape under nrows=1 & totalpulses=2 condition";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				indexToWrite++;
				// Pulse-next-i data stored in *G
				gsl_vector_memcpy(pulseG,pulsenexti);
				tstartG = tstartnexti;
				tendG = tendnexti;
				qualityG = qualitynexti;
				tailbeforeG = tailbeforenexti;
			}
			else if ((totalpulses != 1) && (totalpulses != 2) && (totalpulses != eventcnt))
			// Neither first or second or last 1-inDataIteratorTrg-set
			{
				// Not read pulse-prev-i data because it is already in *F
				gsl_vector_memcpy(pulseprevi,pulseF);
				tstartprevi = tstartF;
				tendprevi = tendF;
				qualityprevi = qualityF;
				tailbeforeprevi = tailbeforeF;
				// Not read pulse-i data because it is already in *G
				gsl_vector_memcpy(pulsei,pulseG);
				tstarti = tstartG;
				tendi = tendG;
				qualityi = qualityG;
				tailbeforei = tailbeforeG;
				// Read pulse-next-i
				gsl_matrix_get_row(pulsenexti,pulses,i);
				tstartnexti = gsl_vector_get(tstartgsl,i);
				tendnexti = gsl_vector_get(tendgsl,i);
				qualitynexti = gsl_vector_get(qualitygsl,i);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,tstartnexti))
				{
				    message = "Error in pulseShape under nrows=1 & totalpulses!=1 ... condition";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				indexToWrite++;
				// Pulse-i data stored in *F
				gsl_vector_memcpy(pulseF,pulsei);
				tstartF = tstarti;
				tendF = tendi;
				qualityF = qualityi;
				tailbeforeF = tailbeforei;
				// Pulse-next-i data stored in *G
				gsl_vector_memcpy(pulseG,pulsenexti);
				tstartG = tstartnexti;
				tendG = tendnexti;
				qualityG = qualitynexti;
				tailbeforeG = tailbeforenexti;
			}
			else if (totalpulses == eventcnt)
			// Last 1-inDataIteratorTrg-set
			{
				// PULSE OF THE PREVIOUS 1-INDATAITERATORTRG-SET
				// Not read pulse-prev-i data because it is already in *F
				gsl_vector_memcpy(pulseprevi,pulseF);
				tstartprevi = tstartF;
				tendprevi = tendF;
				qualityprevi = qualityF;
				tailbeforeprevi = tailbeforeF;
				// Not read pulse-i data because it is already in *G
				gsl_vector_memcpy(pulsei,pulseG);
				tstarti = tstartG;
				tendi = tendG;
				qualityi = qualityG;
				tailbeforei = tailbeforeG;
				// Read pulse-next-i
				gsl_matrix_get_row(pulsenexti,pulses,i);
				tstartnexti = gsl_vector_get(tstartgsl,i);
				tendnexti = gsl_vector_get(tendgsl,i);
				qualitynexti = gsl_vector_get(qualitygsl,i);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,tstartnexti))
				{
				    message = "Error in pulseShape under nrows=1 & totalpulses=eventcnt condition";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
				indexToWrite++;
				//PULSE OF THE ACTUAL 1-INDATAITERATORTRG-SET
				// Pulse-i data transferred to pulse-prev-i
				gsl_vector_memcpy(pulseprevi,pulsei);
				tstartprevi = tstarti;
				tendprevi = tendi;
				qualityprevi = qualityi;
				tailbeforeprevi = tailbeforei;
				// Pulse-next-i data transferred to pulse-i
				gsl_vector_memcpy(pulsei,pulsenexti);
				tstarti = tstartnexti;
				tendi = tendnexti;
				qualityi = qualitynexti;
				tailbeforei = tailbeforenexti;
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				// Last pulse of the file (tstartnexti = -1 to pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,-1.0))
				{
				    message = "Error in pulseShape under nrows =1 & last pulse condition";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
			}
		}
		else if (nrows == 2)
		{
			if ((i == 0) && (totalpulses == 1))
			// First row of the first 2-inDataIteratorTrg-set
			{
				// Read pulse-i
				gsl_matrix_get_row(pulsei,pulses,i);
				tstarti = gsl_vector_get(tstartgsl,i);
				tendi = gsl_vector_get(tendgsl,i);
				qualityi = gsl_vector_get(qualitygsl,i);
				tailbeforei = gsl_vector_get(tailbeforegsl,i);
			}
			else if ((i == 1) && (totalpulses == 2))
			// Second row of the first 2-indataIteratorTrg-set
			{
				// Read pulse-next-i
				gsl_matrix_get_row(pulsenexti,pulses,i);
				tstartnexti = gsl_vector_get(tstartgsl,i);
				tendnexti = gsl_vector_get(tendgsl,i);
				qualitynexti = gsl_vector_get(qualitygsl,i);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				// First pulse of the file (tendprevi = -1 to pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,-1.0,-1.0,tstartnexti))
				{
				    message = "Error in pulseShape under nrows=2 & totalpulses=2 condition";
				    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
				}
				indexToWrite++;
				// Pulse-i stored in *F
				gsl_vector_memcpy(pulseF,pulsei);
				tstartF = tstarti;
				tendF = tendi;
				qualityF = qualityi;
				tailbeforeF = tailbeforei;
				// Pulse-next-i stored in *G
				gsl_vector_memcpy(pulseG,pulsenexti);
				tstartG = tstartnexti;
				tendG = tendnexti;
				qualityG = qualitynexti;
				tailbeforeG = tailbeforenexti;
			}
			else if ((i == 0) && (totalpulses != 1))
			// First row of a 2-inDataIteratorTrg-set different from the first one
			{
				// Not read pulse-prev-i because it is already in *F
				gsl_vector_memcpy(pulseprevi,pulseF);
				tstartprevi = tstartF;
				tendprevi = tendF;
				qualityprevi = qualityF;
				tailbeforeprevi = tailbeforeF;
				// Not read pulse-i because it is already in *G
				gsl_vector_memcpy(pulsei,pulseG);
				tstarti = tstartG;
				tendi = tendG;
				qualityi = qualityG;
				tailbeforei = tailbeforeG;
				// Read pulse-next-i
				gsl_matrix_get_row(pulsenexti,pulses,i);
				tstartnexti = gsl_vector_get(tstartgsl,i);
				tendnexti = gsl_vector_get(tendgsl,i);
				qualitynexti = gsl_vector_get(qualitygsl,i);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,tstartnexti))
				{
				    message = "Error in pulseShape under nrows=2 & totalpulses!=1 condition";
				    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
				}
				indexToWrite++;
			}
			else if ((i == 1) && (totalpulses != 2))
			// Second row of a 2-inDataIteratorTrg-set different from the first one
			{
				// Pulse-i data transferred to pulse-prev-i
				gsl_vector_memcpy(pulseprevi,pulsei);
				tstartprevi = tstarti;
				tendprevi = tendi;
				qualityprevi = qualityi;
				tailbeforeprevi = tailbeforei;
				// Pulse-next-i data transferred to pulse-i
				gsl_vector_memcpy(pulsei,pulsenexti);
				tstarti = tstartnexti;
				tendi = tendnexti;
				qualityi = qualitynexti;
				tailbeforei = tailbeforenexti;
				// Read pulse-next-i
				gsl_matrix_get_row(pulsenexti,pulses,i);
				tstartnexti = gsl_vector_get(tstartgsl,i);
				tendnexti = gsl_vector_get(tendgsl,i);
				qualitynexti = gsl_vector_get(qualitygsl,i);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,tstartnexti))
				{
				    message = "Error in pulseShape under nrows=2 & totalpulses!=2 condition";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL); 
				}
				indexToWrite++;
				if (totalpulses == eventcnt)
				// The 2-inDataIteratorTrg-set is the last one
				{
					// Pulse-i data transferred to pulse-prev-i
					gsl_vector_memcpy(pulseprevi,pulsei);
					tstartprevi = tstarti;
					tendprevi = tendi;
					qualityprevi = qualityi;
					tailbeforeprevi = tailbeforei;
					// Pulse-next-i data transferred to pulse-i
					gsl_vector_memcpy(pulsei,pulsenexti);
					tstarti = tstartnexti;
					tendi = tendnexti;
					qualityi = qualitynexti;
					tailbeforei = tailbeforenexti;
					// Analyze pulse-i (and write pulse-i inside pulseShape)
					// Last pulse of the file (tstartnexti = -1 to pulseShape)
					if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,-1.0))
					{
					    message = "Error in pulseShape under nrows=2 & totalpulses=eventcnt condition";
					    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
					}
				}
				else
				// The 2-inDataIteratorTrg-set is not the last one
				{
					// Pulse-i data stored in *F
					gsl_vector_memcpy(pulseF,pulsei);
					tstartF = tstarti;
					tendF = tendi;
					qualityF = qualityi;
					tailbeforeF = tailbeforei;
					// Pulse-next-i data stored in *G
					gsl_vector_memcpy(pulseG,pulsenexti);
					tstartG = tstartnexti;
					tendG = tendnexti;
					qualityG = qualitynexti;
					tailbeforeG = tailbeforenexti;
				}
			}
		}
		else
		// nrows different from 1 or 2
		{
			if ((i == 0) && (totalpulses == 1))
			// First row of the first nrows-inDataIteratorTrg-set
			{
				// Read pulse-i data
				gsl_matrix_get_row(pulsei,pulses,i);
				tstarti = gsl_vector_get(tstartgsl,i);
				tendi = gsl_vector_get(tendgsl,i);
				qualityi = gsl_vector_get(qualitygsl,i);
				tailbeforei = gsl_vector_get(tailbeforegsl,i);
				// Read pulse-next-i data
				gsl_matrix_get_row(pulsenexti,pulses,i+1);
				tstartnexti = gsl_vector_get(tstartgsl,i+1);
				tendnexti = gsl_vector_get(tendgsl,i+1);
				qualitynexti = gsl_vector_get(qualitygsl,i+1);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i+1);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				// First pulse of the file (tendprevi = -1 to pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,-1.0,-1.0,tstartnexti))
				{
				    message = "Error in pulseShape under nrows>2 & totalpulses=1 condition";
				    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
				}
				indexToWrite++;
			}
			else if ((i == 0) && (totalpulses != 1))
			// First row of a nrows-inDataIteratorTrg-set different from the first one
			{
				// LAST PULSE OF THE PREVIOUS NROWS-INDATAITERATORTRG-SET
				// Not read pulse-prev-i data because it is already in *F
				gsl_vector_memcpy(pulseprevi,pulseF);
				tstartprevi = tstartF;
				tendprevi = tendF;
				qualityprevi = qualityF;
				tailbeforeprevi = tailbeforeF;
				// Not read pulse-i data because it is already in *G
				gsl_vector_memcpy(pulsei,pulseG);
				tstarti = tstartG;
				tendi = tendG;
				qualityi = qualityG;
				tailbeforei = tailbeforeG;
				// Read pulse-next-i data
				gsl_matrix_get_row(pulsenexti,pulses,i);
				tstartnexti = gsl_vector_get(tstartgsl,i);
				tendnexti = gsl_vector_get(tendgsl,i);
				qualitynexti = gsl_vector_get(qualitygsl,i);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,tstartnexti))
				{
				    message = "Error in pulseShape under nrows>2 & Last Pulse condition";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL); 
				}
				indexToWrite++;
					//FIRST PULSE OF THE NROWS-INDATAITERATORTRG-SET
				// Pulse-i data transferred to pulse-prev-i
				gsl_vector_memcpy(pulseprevi,pulsei);
				tstartprevi = tstarti;
				tendprevi = tendi;
				qualityprevi = qualityi;
				tailbeforeprevi = tailbeforei;
				// Pulse-next-i data transferred to pulse-i
				gsl_vector_memcpy(pulsei,pulsenexti);
				tstarti = tstartnexti;
				tendi = tendnexti;
				qualityi = qualitynexti;
				tailbeforei = tailbeforenexti;
				// Read pulse-next-i data
				gsl_matrix_get_row(pulsenexti,pulses,i+1);
				tstartnexti = gsl_vector_get(tstartgsl,i+1);
				tendnexti = gsl_vector_get(tendgsl,i+1);
				qualitynexti = gsl_vector_get(qualitygsl,i+1);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i+1);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,tstartnexti))
				{
				    message = "Error in pulseShape under nrows>2 & First Pulse condition";
				    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
				}
				indexToWrite++;
			}
			else if ((i != 0) && (i < nrows-1))
			// Neither first pulse or last pulse of a nrows-indataIteratorTrg-set
			{
				// Pulse-i data transferred to pulse-prev-i
				gsl_vector_memcpy(pulseprevi,pulsei);
				tstartprevi = tstarti;
				tendprevi = tendi;
				qualityprevi = qualityi;
				tailbeforeprevi = tailbeforei;
				// Pulse-next-i data transferred to pulse-i
				gsl_vector_memcpy(pulsei,pulsenexti);
				tstarti = tstartnexti;
				tendi = tendnexti;
				qualityi = qualitynexti;
				tailbeforei = tailbeforenexti;
				// Read pulse-next-i data
				gsl_matrix_get_row(pulsenexti,pulses,i+1);
				tstartnexti = gsl_vector_get(tstartgsl,i+1);
				tendnexti = gsl_vector_get(tendgsl,i+1);
				qualitynexti = gsl_vector_get(qualitygsl,i+1);
				tailbeforenexti = gsl_vector_get(tailbeforegsl,i+1);
				// Analyze pulse-i (and write pulse-i inside pulseShape)
				if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,tstartnexti))
				{
				    message = "Error in pulseShape under nrows>2 & not first/last pulse condition";
				    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL); 
				}
				indexToWrite++;
			}
			else if ((i != 0) && (i == nrows-1))
			// Last pulse of a nrows-inDataIteratorTrg-set
			{
				// Pulse-i data transferred to pulse-prev-i
				gsl_vector_memcpy(pulseprevi,pulsei);
				tstartprevi = tstarti;
				tendprevi = tendi;
				qualityprevi = qualityi;
				tailbeforeprevi = tailbeforei;
				// Pulse-next-i data transferred to pulse-i
				gsl_vector_memcpy(pulsei,pulsenexti);
				tstarti = tstartnexti;
				tendi = tendnexti;
				qualityi = qualitynexti;
				tailbeforei = tailbeforenexti;
				if (totalpulses == eventcnt)
				// Last pulse of the file (tstartnexti = -1 to pulseShape)
				{
					if (pulseShape (pulsei,tstarti,tendi,qualityi,tailbeforei,tstartprevi,tendprevi,-1.0))
					{
					    message = "Error in pulseShape under nrows>2 & Last Pulse in file condition";
					    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
					}
				}
				else
				// Last pulse of a nrows-inDataIteratorTrg-set
				{
					// Pulse-prev-i data stored in *F
					gsl_vector_memcpy(pulseF,pulseprevi);
					tstartF = tstartprevi;
					tendF = tendprevi;
					qualityF = qualityprevi;
					tailbeforeF = tailbeforeprevi;
					// Pulse-i data stored in *G
					gsl_vector_memcpy(pulseG,pulsei);
					tstartG = tstarti;
					tendG = tendi;
					qualityG = qualityi;
					tailbeforeG = tailbeforei;
				}
			}
		}
		totalpulses ++;
		ntotalrows ++;
	}

	iteration = iteration + nrows;

	// Free allocate of GSL vectors
	gsl_matrix_free (pulses);
	gsl_vector_free (pulse);
	gsl_vector_free (tstartgsl);
	gsl_vector_free (tendgsl);
	gsl_vector_free (qualitygsl);
	gsl_vector_free (tailbeforegsl);

	gsl_vector_free (pulsei);
	gsl_vector_free (pulsenexti);

	return(EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ***************************************
* pulseShape function: This function obtains the features of each pulse and writes data in the output FITS file
*                      and in the _trg.fits file
*
* - Declare variables
* - Allocate GSL vectors
* - Store data in vectors
* - Grading the pulses (called 'events' in this case)
* - Write each pulse into the output FITS file
* 		- Output columns: Tstart, Tend, Quality, TailBfr and Grade
* - Add Grade column to the TRIGGER FITS file
* - Free allocate of GSL vectors
* - Free memory
* *********************************************/
int pulseShape (gsl_vector *pulse, double tstart, double tend, double qual, double tail, double tstartprev, double tendprev, double tstartnext)
{
	int status=EPOK;

	// Declare variables
	int size = pulse->size;
	int index_peak;						// Index of the maximum voltage of the pulse
	string message = "";

	// Allocate GSL vectors
	// It is necessary that the data to write in the output FITS file (by using writeFitsSimple) is stored in a vector
	// Here, the vector has an only element because it is going to write pulse by pulse (row by row)
	gsl_vector *tstartoutgsl = gsl_vector_alloc(1);
	gsl_vector *tendoutgsl = gsl_vector_alloc(1);
	gsl_vector *qualoutgsl = gsl_vector_alloc(1);
	gsl_vector *tailoutgsl = gsl_vector_alloc(1);
	gsl_vector *gradegsl = gsl_vector_alloc(1);

	index_peak = gsl_vector_max_index (pulse);

	// Store data in vectors
	gsl_vector_set(tstartoutgsl,0,tstart);
	gsl_vector_set(tendoutgsl,0,tend);
	gsl_vector_set(qualoutgsl,0,(int) (qual));
	gsl_vector_set(tailoutgsl,0,(int) (tail));
	gsl_vector_set(gradegsl,0,0);

	// Grading the pulses (called 'events' in this case)
	char val[256];
	if (eventGrade(tstart,tstartprev,tstartnext,&gradegsl))
	{
	    message = "Error during Event Grading";
	    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
	}

	if ((qual == 0) && (tail ==1))	gsl_vector_set(gradegsl,0,32);

 	// Creating Tstart Column
 	IOData obj;
 	obj.inObject = pshObject;
 	obj.nameTable = new char [255];
 	strcpy(obj.nameTable,"EUR-PSH");
 	obj.iniRow = indexToWrite;
 	obj.endRow = indexToWrite;
 	obj.iniCol = 0;
 	obj.nameCol = new char [255];
 	strcpy(obj.nameCol,"Tstart");
 	obj.type = TDOUBLE;
 	obj.unit = new char [255];
 	strcpy(obj.unit,"seconds");
	
 	if (writeFitsSimple (obj,tstartoutgsl))
 	{
	    message = "Error creating " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
	}

 	// Creating Tend Column
 	strcpy(obj.nameCol,"Tend");
 	if (writeFitsSimple (obj,tendoutgsl))
 	{
	    message = "Error creating " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
	}

 	// Creating Quality Column
	strcpy(obj.nameCol,"Quality"); obj.type = TSHORT; strcpy(obj.unit," ");
	if (writeFitsSimple (obj,qualoutgsl))
	{
	    message = "Error creating " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
	}

	// Creating TailBfr Column
	strcpy(obj.nameCol,"TailBfr"); obj.type = TSHORT; strcpy(obj.unit," ");
	if (writeFitsSimple (obj,tailoutgsl))
	{
	    message = "Error creating " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
	}

	// Creating Grade Column
	strcpy(obj.nameCol,"Grade"); obj.type = TSHORT; strcpy(obj.unit," ");
	if (writeFitsSimple (obj,gradegsl))
	{
	    message = "Error creating " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
	}

	// Add Grade Column to TRIGGER FITS
 	IOData objTrg;
 	objTrg.inObject = trgObject;
 	objTrg.nameTable = new char [255];
 	strcpy(objTrg.nameTable,"EUR-TRG");
 	objTrg.iniRow = indexToWrite;
 	objTrg.endRow = indexToWrite;
 	objTrg.iniCol = 0;
 	objTrg.nameCol = new char [255];
 	strcpy(objTrg.nameCol,"Grade");
 	objTrg.type = TSHORT;
 	objTrg.unit = new char [255];
 	strcpy(objTrg.unit," ");
 	if (writeFitsSimple(objTrg,gradegsl))
 	{
	    message = "Error creating " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);  return(EPFAIL);
	}

	// Free allocate of GSL vectors
	gsl_vector_free(tstartoutgsl);
	gsl_vector_free(tendoutgsl);
	gsl_vector_free(qualoutgsl);
	gsl_vector_free(tailoutgsl);

 	// Free memory
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;
	delete [] objTrg.nameTable;
	delete [] objTrg.nameCol;
	delete [] objTrg.unit;

	return EPOK;
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* eventGrade: This function judges the event grade as defined in the following table.
*
* 								 	tp<=innerInterval	innerInterval<tp<=outerInterval	  outerInterval<tp
*               tn<=innerInterval          Ls                        Ls                         Lp
* innerInterval<tn<=outerInterval		   Ls	                     Ms                         Mp
* outerInterval<tn                         Ls                        Ms                         Hp
*
* tp = Tstart-Tstartprev (arrival time of the preceding pulse)
* tn = Tstartnext-Tstart (arrival time of the next pulse)
*
* Hp: High primary  (Grade=1)
* Mp: Med primary   (Grade=21)
* Ms: Med secondary (Grade=22)
* Lp: Low primary   (Grade=31)
* Ls: Low secondary (Grade=32)
*
****************************************/
int eventGrade (double Tstart, double Tstartprev, double Tstartnext, gsl_vector **Gradinggsl)
{
	char val[256];
	char val_aux[256];

	int status = EPOK;
	string message = "";

	double innerInterval, outerInterval;

	innerInterval = inTaus*tauFALL;
	outerInterval = outTaus*tauFALL;

	if (Tstartnext == -1)	Tstartnext = Tstart*1e6;	// Last pulse of the file

	if (((Tstart-Tstartprev) > outerInterval) && ((Tstartnext-Tstart) > outerInterval))
	{
		// High primary (Hp)
		gsl_vector_set(*Gradinggsl,0,1);
	}
	else if (((Tstart-Tstartprev) > innerInterval) && ((Tstartnext-Tstart) > innerInterval))
	{
		if ((Tstart-Tstartprev) > outerInterval)
		{
			// Med primary (Mp)
			gsl_vector_set(*Gradinggsl,0,21);
		}
		else
		{
			// Med secondary (Ms)
			gsl_vector_set(*Gradinggsl,0,22);
		}
	}
	else if (((Tstart-Tstartprev) < innerInterval) || ((Tstartnext-Tstart) < innerInterval))
	{
		if ((Tstart-Tstartprev) > outerInterval)
		{
			// Low primary (Lp)
			gsl_vector_set(*Gradinggsl,0,31);
		}
		else
		{
			// Low secondary (Ls)
			gsl_vector_set(*Gradinggsl,0,32);
		}
	}
	else if (((Tstart-Tstartprev) < innerInterval) && ((Tstartnext-Tstart) < innerInterval))
	{
		// Low secondary (Ls)
		gsl_vector_set(*Gradinggsl,0,32);
	}

	return (EPOK);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
