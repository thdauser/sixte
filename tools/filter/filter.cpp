/************************************************************************************************

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

/************************************************************************************************
*               			  FILTER
*
*  File:      filter.cpp
*  Version:   14.0.0
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Prez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
* version 1.0.0		12/03/08    First version
* version 2.0.0		07/04/08	Filtering of pulse using the column "Quality"
* 								Adding the output column "quality"
* version 3.1.0		26/06/08	Adding input/output keyword "ql"
* 								Adding Error Procesing: "The input FITS file hasn't run in PSH module."
* 								New function to filter pulses not valid of TRG file
*								Adding baseline and normalize function
* 								Adding input parameters "dt_before", "cut_filter" and "dt_after"
* 								Deleting input/output keyword "ql";
* 								Adding input parameter "ql"
* 								New input FITS file "Pulse Shape"
* version 3.2.0		30/09/08	Adding new function "readInputKeywords"
* 								Deleting in/output keyword "channelcount"
* 								Including error codes.
* version 4.0.0 	29/09/08	Included "CREATE" and "PROCESS" keyword in the output FITS file
* 								Include error processing of "INPUT_PARAMETER_NON_DEFINE"
* 								Include input parameters: "nameLog" and "verbosity"
* version 4.1.0		24/09/08	If the module fail. It will be warning in the output.
* 								Included parameter time to measure.
* version 4.2.0		19/09/08	Included Warning about the size of the filter
* version 5.0.0		14/01/09 	Included new Warnings of Input Parameters values.
* 								Changed the output column name of Current from "V" to "I0"
* 								Included new input keywords
* version 5.1.0		22/01/09	Changed the verbosity of any logs.
* 								Modified error code of input parameters
* version 5.2.0		02/02/09	Changed name of output keyword "CREATE" to "CREATOR"
* 								Create a HISTORY in pulse shape file.
* version 5.2.1		04/02/09	Resolved bug in output FIT file.(PSH)
* 								Checked free memory allocate.
* version 6.0.0		20/02/09	Included, modified and removed input keywords
* version 6.1.0		24/02/09	Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
* 								Solved bug 4: stops on zero value for the keyword TIMEZERO, error -62436.
* 								Included FTYPE keywords in output FITS file.
* version 6.2.0		12/03/09	Included new column called "filterError" in the output FITS file.
* 								Created obtainFilterError function.
* version 7.0.0					Resolver bug about max_size of the filter variable.
* 								Used the keywords libraries
* version 7.0.1     16/04/09    Delete "PILInit(argc, argv);" and "PILClose(PIL_OK);"
* 								Documentation updated
* version 8.0.0     08/06/09    If Norris fit columns not in PSH file, not read here
* version 8.0.1     18/06/09    Repaired comment "PulseShape has been run in QuickLook Mode, then EnergyResol too"
* 								Changed order in delete's
* version 9.0.0		24/11/09	If error file exists, its contests will be deleted
* 								endPulse is not a hardpoint but a keyword read from _trg input FITS file (EVENTSZ)
* 								createFilter: When cutting the pulse, if it is necessary, dt_beforeCF and dt_afterCF
* 											  (local variables) are changed with respect the dt_before and dt_after
*                                             given by the user as input parameters (global variables). If they have
*                                             to be changed, the change only affects to the particular pulse.
* version 9.0.1		25/11/09	The EVENTSZ keyword read from _trg input FITS file is not written in the output file
*                   30/11/09    The PROCESS keyword written in the modified PSH file keeps the PROCESS keyword read from the
*                               PSH input FITS file and the PROCESS generated in FILTER and the FILTER version are added
*                   02/12/09    If there are no valid pulses in the input FITS file => The task finishes (but no in an uncontrolled way)
* version 9.0.2     04/12/09    Bug 28 solved
* version 10.0.0	24/09/10	Added input parameter "tmin"
* 								Deleted input parameters "dt_after" & "dt_before"
* 								Normalize the filter template ("normalize" function added)
*								Not necessary cut the pulse from "tstart-dt_before" to "tend+dt_after" ("cut_filter" function deleted)
*								Deleted extension BASELINE in input TRIGGER FITS file
*								From TRIGGER FITS file: Necessary column "difTstrt" to consider valid or non-valid pulses and columns "TIME",
*								"EndPulse", "Baseline" and "Sigma" deleted
*								From PULSESHAPE FITS file: columns "TIME", "Baseline", "Sigma", "MaxTimeFit", "MaxCurrentFit", "TriseFit",
*								"NriseFit", "TfallFit", "NfallFit", "ChisqFit" and "BaselineFit" deleted
*								Deleted "writeFitsSimple2" function
*								Documentation updated
* version 10.1.0	30/09/10	The value of the keywords "EVENTCNT" and "EVENTSZ" are changed
* 					08/11/10	FILTER task obtains the filter template using non-filtered pulses instead filtered pulses
* 								Documentation updated
* version 11.0.0	14/12/10	'Log' argument instead 'OK' in some "writeLog" functions
* 					21/12/10	From PULSESHAPE FITS file: columns "Trise", "Nrise", "Tfall" and "Nfall" deleted
* 								Documentation updated
* 					??/??/??	Filter template with equal areas below and above the x-axis
* 					23/02/11    I0 instead I0_NotFiltered
*                   28/02/11    PROCESS also includes FILTER version (CREATOR)
*                   01/03/11    Neither HISTORY, FLT_ID or PROCESS written in _psh, because psh is read but not modified
*                   17/03/11    _psh input FITS file is no longer used: QUALITY will be read from _trg
*                               	inDataIteratorPsh deleted and inDataIteratorTrg modified
*                                   getKh deleted => getKhaux renamed as getKh
*                               "columnName" input parameter not used
* version 12.0.0    25/03/11    Quick Look mode no longer used
*                   04/04/11    No pulses in the _trg input FITS file => Warning and the task finishes
*                   31/01/12    New _noisespec input FITS file
*                               New "inDataIteratorTPS" function
*                               "obtainFilterError" changed to handle local variables
*                               Optimal filter in time and frequency domain as output
*                               	"interpolate", "reorderFFT" and "disorderFFT" new functions
*                   05/03/12    Nor "difTstrt" or "tmin"(input parameter) are not going to be used => Delete
*                               Grade column of _trg is read and used to select the pulses to average in order to obtain the optimal filter
*                   12/02/13    New output keyword "MODE"
*                   31/07/13    Added the "OPTFSIZE" keyword in the extension EUR-FLT of the _trg FITS file
*                   25/03/13    New input parameters to take into account an energy window: Emin, Emax
*                               Necessary changes to divide the thread between CALIBRATION and NORMAL mode
*                               NORMAL mode: Read EstEnrgy column from _trg input FITS file
*                               Added as input file the pulse model library
*                   xx/03/13    Cut the ending part of the sum of all the pulses to avoid an undesirable minimum (an artifact)
*                   10/05/13    New function "find_matchedfilter"
*                   15/05/13    New output keywords "B_CF" and "C_CF"
*                   05/06/13    To avoid undesirable maximums, not use the maximum of the pulse to sum the maximum of the pulse and the maximum of the filter (create_filter)
*                               Instead, use an average of the previous and the following samples to the maximum of the pulse.
*                               'pulse_value = (gsl_vector_get(pulse,index_peak-1) + gsl_vector_get(pulse,index_peak+1))/2;'
*                   06/06/13    Improved the cutting of the ending part of sum of all the pulses to avoid an undesirable minimum('filterlength')
*                   12/06/13    Rewrote the normal operation mode code
*                   27/11/13    ANNALS in EUR-TRG from _trg FITS file is modified (in case of normal mode)
*                   30/01/14    The pulse model library input FITS file is read (readLib) in a different point of the main thread
	                   			Added the library FITS file name to the output keyword PROCESS
*                   31/01/14    EVENTCNT in EUR-FLT (_flt) is written before optimalfiltergsl was free allocated
*                               (and it is not the 'filterlength' but the 'optimalfiltergsl->size' because of n0start and n0end)
*                   07/01/14    New output column MF_L (matched filter length) to know where the artifact at the end of the
*                               matched filter starts (mode = 0)
*                   20/02/14    New output column OF_L (optimal filter length) to know where the artifact at the end of the
*                               optimal filter starts
*                   14/05/14    New output column NRMFCTR (normalization factor)
*                   01/07/14    Deleted the "OPTFSIZE" keyword in the extension EUR-FLT of the _trg FITS file (not used in ERS)
*                               Deleted the 'normalize' function
*                               New 'area0' function
*                               calculus_matchedFilter: - Not 'normalize'
*                                                       - Divide by the energy(eV)
*                                                       - 'area0' used
*                               calculus_optimalFilter: - Calculate the normalization factor
*                                                       - Not normalize the optimal filter
*                               'area0' used after obtaining the matched filter by interpolating (normal mode)
* version 13.0.0    07/07/14	Removed PIL/RIL/Common dependencies. Adapted use of writeLog routine
* version 13.1.0    08/07/14    Adapted for new parameter in interativePars function for task name
* version 13.2.1    11/07/14    Solved bug preventing the task from reading the full command line
* version 13.2.2    11/07/14    Comments of the 'initModule' function modified
*                   15/07/14    'initModule' modified in order to accept negative int or double (if parameters are read from command line)
*                   31/07/14    Comments referring to PIL and RIL deleted
*                   12/09/14    Differences when reading or writing some keywords
*                   16/09/14    In 'calculus_optimalFilter', when  Af{MatchedFilter(f)} < Af{N(f)} => Interpolation of N(f):
*	    						'gsl_vector_memcpy(f2,&temp.vector);' instead 'gsl_vector_memcpy(f1,&temp.vector);'
*	   	                	   	'gsl_vector_memcpy(FFT2,&temp.vector);' instead 'gsl_vector_memcpy(FFT1,&temp.vector);'
*	            21/10/14    In 'calculus_optimalFilter',  when you go backwards and Af{MatchedFilter(f)} < Af{N(f)} => Interpolation of N(f):
*	   	                        Vector allocation
* version 14.0.0    Dec/14    DAL->CFITSIO migration
*   		                    Deleted some input/output keywords
*   		                    Deleted some unnecessary input parameters (Emin, Emax)
*   		                    In calibration mode, it is not necessary to read from TRIGGER file any more because
*   		                    the average of the pulses is already in the pulse templates library (PULSE)
*   		                    (deleted inDataIteratorTrg0)
*                   09/01/15	'short' instead of 'int' when QUALITY and TAILBFR columns are read (in 'inDataIteratorTrg')
*                   			'writeLib' modified (TAURISE, TAUFALL and PRETRIGS columns in ANX-LIB extension are not going to be written)
*                   12/01/15    'PROCESS' input keyword from _trg input FITS file is not read
*                               'PROCESS' output keyword in _flt output FITS file renamed as 'PROC0'
*					            'PROCESS' output keyword in _trg output FITS file divided in three ones ('PROC0' generated by TRIGGER, 'PROC1'
*					            generated by PULSEGRADE and 'PROC2' generated by FILTER)
*					            'ANNALS' input keyword from _trg input FITS file is not read
*					            'ANNALS' output keyword in _trg output FITS file divided in two ones ('MOD0' generated by PULSEGRADE and 'MOD1'
*            					generated by FILTER)
*            		20/01/15    Clobber parameter added. Some renameing of variables/FITS extensions (EUR-TRG -> TRIGGER, EUR-LIB --> LIBRARY, EUR-TEST->TEST)
*					22/01/15    Added error condition if in calibration mode, energy is not found in calibration library
*
****************************************************************************************/

/******************************************************************************
DESCRIPTION:

This document describes the task FILTER .

The goal of the FILTER task is obtain an optimal filtering both in time and in frequency domain.
There are two operation modes of the FILTER task: calibration mode and production mode.

- If the calibration mode is selected, the average of the pulses is already in the pulse templates library input FITS file.
  Therefore, to have the matched filter, the filter template read from the library is vertically shifted (in time domain)
  to set its area to zero (with equal areas below and above the x-axis). In the frequency domain it means to remove the
  bin 0, i.e., the baseline essentially in order to calculate properly the pulse energy in the next tasks.

  The matched filter will be store in the pulse templates library FITS file.

  Moreover, the matched filter must be also normalized with the power spectrum of the noise in order to obtain an
  optimal filtering.  The optimal filtering will be stored in the _flt output FITS file.

- If the production mode is selected, the proper matched filter to each valid pulse is chosen from the pulse templates
  library FITS file. Then, an optimal filter is calculated to each pulse and it will be stored in the proper row
  of a new extension of the _trg input FITS file (FILTER).

  To classify pulses as valid ones, quality and grade are used.

The users must supply the following input parameters:

- trgFile: Name of the _trg.fits file which contains all the found pulses
- noisespecFile: Name of the _noisespec.fits file which contains the current noise spectral density
- inLibFile: Name of the pulse templates library FITS file
- fltFile: Name of the _flt.fits output file to write the optimal filter (only in calibration mode)
- Hp, Mp, Ms, Lp, Ls: Some parameters to choose the pulses according to its grade to be processed to obtain the matched filter
- namelog: Output log file name
- verbosity: Verbosity level of the output log file

 MAP OF SECTIONS IN THIS FILE:

 - 1. INCLUDE's
 - 2. MAIN
 - 3. initModule
 - 4. readLib
 - 5. inDataIteratorLib0
 - 6. inDataIteratorLib1
 - 7. createFilterFile
 - 8. inDataIteratorTPS
 - 9. find_estenergy
 - 10. inDataIteratorTrg
 - 11. calculus_matchedFilter
 - 12. writeLib
 - 13. calculus_optimalFilter
 - 14. reorderFFT
 - 15. interpolate
 - 16. disorderFFT
 - 17. find_matchedfilter

*******************************************************************************/

/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include "filter.h"

/***** SECTION 2 ************************************
* MAIN function: This function is the main function of the FILTER task
*
* - Read input parameters (call initModule)
* - Open input TRIGGER FITS file
* - Open input NOISESPEC FITS file (NOISEALL extension)
* - Read input keywords and check their values
* - If there are no pulses in _trg input FITS file (empty TRIGGER extension) => Provide a warning and finish
* - Get structure of input FITS fileS columns (from _trg and _noisespec)
* - Read the pulse templates library input FITS file (call readLib)
* - Handle the output FITS File (_flt.fits in calibration mode or _trg.fits in production mode) (call createFilterFile)
* - Create structure to run Iteration: inDataIteratorTPS
* - Read columns (FREQ and CSD)
* - Called iteration function: inDataIteratorTPS
* - If CALIBRATION mode
* 	- Read ENERGY keyword from the _trg input FITS file (the input file is monochromatic)
* 	  If ENERGY does not match to any energy of the pulse templates library, a new row in the ANX-LIB extension
*     of the library FITS file is created. The matched filter will be written in this new row, although this
*     one will not have all the info (neither ESTENERGY nor TAURISE nor TAUFALL nor PULSE)
* 	- Call find_energy (to know in which row from the library the calculated matched filter is going to be written)
* 	- Call calculus_matchedFilter
* 	  Averaged of pulses have already been done in PULSE (library), then it is only necessary vertically shift the
* 	  pulse template in order to have equal areas below and above the x-axis (in order to be baseline insensitive
* 	  when calculating the pulses energy)
* 	- Write library FITS file (MacthedF column) (called writeLib)
* 	- Close pulses templates library FITS file (LIBRARY extension)
* 	- Calculate the optimal filter (call calculus_optimalFilter)
* 		- FFT calculus of the matched filter:
* 			- FFT calculus
* 			- Generation of the f's (positive and negative)
* 			- Magnitude and argument for positive and negative f's
* 			- Order the f's and magnitudes according [-fmax,...,0,...,fmax]
* 		- Calculus of the optimal filter:
* 			- Also order the f's and csd's of N(f) according [-fmax,...,0,...,fmax]
* 			- To divide MatchedFilter(f)/N^2(f) => Interpolate the one what has a frequency step bigger (matched filter or N(f)):
* 				- Af{MatchedFilter(f)} > Af{N(f)} => Interpolation of MatchedFilter(f)
* 				- Af{MatchedFilter(f)} < Af{N(f)} => Interpolation of N(f)
* 			- Calculation: MatchedFilter(f)/N^2(f)
* 			- Calculus of the normalization factor
* 			- It is necessary to go backwards in order to have the same number of points than MatchedFilter(f):
* 				- Af{MatchedFilter(f)} > Af{N(f)} => Interpolation of OptimalFilter(f)
*			 	- Af{MatchedFilter(f)} < Af{N(f)} => Already the same number of points than MatchedFilter(f)
* 			- Disorder the f's and magnitudes according [0,...,fmax,-fmax,...] in order to store the OptimalFilter(f) in the output file
* 			- Inverse FFT (to get the expression of the optimal filter in the time domain)
*				- Complex OptimalFilter(f) => Taking into account the magnitude (MatchedFilter(f)/N^2(f)) and the phase
*				  (given by MatchedFilter(f))
*				- OptimalFilter(t) (It is necessary to delete the last argfilterFFT->size/2 elements)
* 	- Write the output FITS file (_flt.fits)
* 	- Write output keywords
* 	- Close output FITS file (_flt.fits)
* 	- Free allocate of GSL vectors
* - elseif PRODUCTION mode
*   - Create structure to run Iteration: inDataIteratorTrg
*	- Read columns: EstEnrgy, Quality and Grade
*	- Call iteration function: inDataIteratorTrg
*   - Check QUALITY: If there are no valid pulses in the input FITS file => The task finishes
*   - Check Hp, Mp, Ms, Lp, Ls and Grade: If there are no valid pulses in the input FITS file => The task finishes
*   - Operate to create the optimal filter:
*   	- Create the extension FILTER in _trg.fits file
*		- For each detected pulse:
*			- If it is a valid pulse according to Quality and Grade:
*				- Obtain the matched filter by interpolating (call find_matchedFilter)
*				- Call area0 (just in case the interpolation modifies the 0 area)
*				- Call calculus_optimalFilter
*				- Write the FILTER extension in _trg FITS file
*			- elseif it is not a valid pulse:
*				- Write 0's in the corresponding row in the FILTER extension in _trg FITS file
* - Close input FITS files
* - Free allocate of memory
* - Finalize the task
****************************************************/
int main (int argc, char **argv)
{
	create = "filter v.14.0.0";			//Set "CREATOR" keyword of output FITS file
	time_t t_start = time(NULL);
	string message = "";
	int status = EPOK, extver=0;
	char *tt[1];
	char *tf[1];
	char *tu[1];
	
	sprintf(temporalFileName,"FILTERauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if (temporalFile == NULL)
	{
	    message = "Cannot open auxiliary file FILTERauxfile.txt";
	    writeLog(fileRef,"Error",verbosity,message);
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Read input parameters
	if (initModule(argc, argv))
	{
		message = "Error in initModule";
		EP_EXIT_ERROR(message,EPFAIL);
	}
	writeLog(fileRef,"Log", verbosity,"Into Filter task");

	//   Open input TRIGGER FITS file
	if (fits_open_file(&trgObject, trgName,READWRITE,&status))
	{
	    message = "Cannot open file " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(extname,"TRIGGER");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(trgName) + " (before keyword reading)";
        EP_EXIT_ERROR(message,status);
	}

	message = "Open TriggerFits: " + string(trgName);
	writeLog(fileRef,"Log", verbosity,message);
	
	// Open input NOISESPEC FITS file (NOISEALL extension)
	// NOISEALL extension contains all the frequencies (positive and negative)
	if (fits_open_file(&noisespecObject, noisespecName,READWRITE,&status))
	{
	    message = "Cannot open file " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(extname,"NOISEALL");
	if (fits_movnam_hdu(noisespecObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(noisespecName) + " (before keyword reading)";
	    EP_EXIT_ERROR(message,status);
	}
	
	message = "Open NoisespecFits: " +  string(noisespecName);
	writeLog(fileRef,"Log", verbosity,message);

	//	Read input keywords and check their values
	strcpy(keyname,"EVENTCNT");
	if (fits_read_key(trgObject,TLONG,keyname, &eventcnt,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	if (eventcnt < 0)
	{
		message = "Legal values for EVENTCNT (TRIGGER) are non negative integer numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	if (eventcnt == 0) // There are no pulses in _trg input FITS file (empty TRIGGER extension)
	{
		time_t t_end = time(NULL);
		message = "There are no pulses in the input FITS file " + string(trgName) +" There is no output FITS file ";
		writeLog(fileRef,"Warning", verbosity,message);
		writeLog(fileRef,"OK", verbosity,"Filter Module OK");

		sprintf(straux,"%f",(double) (t_end - t_start));
		message = "Time: " + string(straux);
		writeLog(fileRef,"Log", verbosity,message);

		if (fclose(fileRef))
		{
		    message = "Cannot close log file";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
	}
	strcpy(keyname,"EVENTSZ");
	if (fits_read_key(trgObject,TLONG,keyname, &eventsz,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	if (eventsz < 0)
	{
		message = "Legal values for EVENTSZ (TRIGGER) are non negative integer numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	strcpy(keyname,"SAMPRATE");
	if (fits_read_key(trgObject,TDOUBLE,keyname, &samprate,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	if (samprate <= 0)
	{
		message = "Legal values for SAMPRATE (TRIGGER) are real numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	/*strcpy(keyname,"ANNALS");
	if (fits_read_key(trgObject,TSTRING,keyname, annalsin,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}*/
	strcpy(keyname,"MODE");
	if (fits_read_key(trgObject,TINT,keyname, &mode,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	if ((mode != 0) && (mode !=1))
	{
		message = "Legal values for MODE (TRIGGER) are 0 or 1";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	strcpy(keyname,"EVENTCNT");
	if (fits_read_key(noisespecObject,TLONG,keyname, &eventcnt_noisespec,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);	  
	}
	if (eventcnt_noisespec <= 0)
	{
		message = "Legal values for EVENTCNT (NOISEALL) are integer numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Get structure of input FITS fileS columns
	//From _trg
	strcpy(extname,"TRIGGER");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) +" in " + string(trgName) + " to get columns structure";
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"I0");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"QUALITY");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"GRADE");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"ESTENRGY");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}

	//From _noisespec
	strcpy(straux,"FREQ");
	if (fits_get_colnum(noisespecObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"CSD");
	if (fits_get_colnum(noisespecObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);
	}
	
	// Read the pulse templates library input FITS file
	if (readLib())
	{
	    message = "Cannot run routine readLib to read pulses library";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Handle the output FITS File (_flt.fits in calibration mode or _trg.fits in production mode)
	if (createFilterFile())
	{
	    message = "Cannot run routine createFilterFile to create output file";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	
	// Iteration NOISEALL extension (NOISESPEC file)
	ntotalrows = 0;
	freqgsl = gsl_vector_alloc(eventcnt_noisespec);
	csdgsl = gsl_vector_alloc(eventcnt_noisespec);

	extern int inDataIteratorTPS(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols,void *user_strct );

	// Create structure to run Iteration: inDatatIteratorTPS
	iteratorCol cols_tps[2];	// Structure of Iteration
	int n_cols_tps = 2; 		// Number of columns: FREQ, CSD
	long rows_per_loop = 0;		// 0: Use default: Optimum number of rows
	long offset = 0;			// 0: Process all the rows

		// Read Columns
	strcpy(extname,"NOISEALL");
	if (fits_movnam_hdu(noisespecObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) +" in " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"FREQ");
	status = fits_iter_set_by_name(&cols_tps[0], noisespecObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"CSD");
	status = fits_iter_set_by_name(&cols_tps[1], noisespecObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);
	}

	 // Called iteration function: inDataIteratorTPS
	if (fits_iterate_data(n_cols_tps, cols_tps, offset, rows_per_loop, inDataIteratorTPS,0L,&status))
	{
	    message = "Cannot iterate data using InDataIteratorTPS";
	    EP_EXIT_ERROR(message,status);
	}

	// Initialize parameters
	qualitygsl	= gsl_vector_alloc(eventcnt);		// GSL vector that contains values of quality column
	gradegsl	= gsl_vector_alloc(eventcnt);		// GSL vector that contains values of grade column
	energygsl	= gsl_vector_alloc(eventcnt);		// GSL vector that contains values of energy column

	IOData obj;
	if (mode == 0)	//Calibration mode
	{
		filtergsl = gsl_vector_alloc(eventsz);			// Filter values
		gsl_vector_set_all(filtergsl,0);
		gsl_vector *row_aux = gsl_vector_alloc(eventsz);

		// Read ENERGY keyword from the _trg input FITS file (the input file is monochromatic)
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(keyname,"ENERGY");
		if (fits_read_key(trgObject,TDOUBLE,keyname, &energy,comment,&status))
		{
		    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
		    EP_EXIT_ERROR(message,status);	  
		}
		// If ENERGY does not match to any energy of the pulse templates library, a new row in the ANX-LIB extension
		// of the library FITS file is created. The matched filter will be written in this new row, although this
		// one will not have all the info (neither ESTENERGY nor TAURISE nor TAUFALL nor PULSE).

		// Call find_energy (to know in which row from the library the calculated matched filter is going to be written)
		if (find_energy (energy, energylibrary, &energyInLibrary_row))
		{
		    message = "Cannot run routine find_energy to find row in library to write match filter";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
		if(energyInLibrary_row==-1)
		{
		    message="Calibration Energy="+boost::lexical_cast<std::string>(energy) + " (eV) must be present in Calibration Library";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
					
		gsl_matrix_get_row(row_aux,templateslibrary,energyInLibrary_row);
		gsl_vector_memcpy(filtergsl,row_aux);

		if (calculus_matchedFilter (&filtergsl, eventsz))
		{
		    message = "Cannot run routine calculus_matchedFilter";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		// Write library FITS file (MacthedF column)
		matchedfilters = gsl_matrix_alloc(1,eventsz);
		gsl_matrix_set_zero(matchedfilters);
		for (int i=0;i<eventsz;i++)
		{
			gsl_matrix_set(matchedfilters,0,i,gsl_vector_get(filtergsl,i));
		}
		if (energyInLibrary_row != -1)
		{
			// Creating MATCHEDF Column
			obj.inObject = inLibObject;
			obj.nameTable = new char [255];
			strcpy(obj.nameTable,"LIBRARY");
			obj.iniCol = 0;
			obj.nameCol = new char [255];
			strcpy(obj.nameCol,"MatchedF");
			obj.type = TDOUBLE;
			obj.unit = new char [255];
			strcpy(obj.unit," ");

			if (matchedf_exist == 1)
			{
				obj.iniRow = energyInLibrary_row+1;
				obj.endRow = energyInLibrary_row+1;

				if (writeFitsComplex(obj, matchedfilters))
				{
				    message = "Cannot run routine writeFitsComplex to write MATCHEDF if already present";
				    EP_EXIT_ERROR(message,EPFAIL);
				}
				
			}
			else if (matchedf_exist == 0)
			{
				gsl_matrix *matrixaux = gsl_matrix_alloc(1,eventsz);
				gsl_matrix_set_zero(matrixaux);
				for (int i=0;i<nummodels;i++)
				{
					obj.iniRow = i+1;
					obj.endRow = i+1;

					if (writeFitsComplex(obj, matrixaux))
					{
					    message = "Cannot run routine writeFitsComplex to write MATCHEDF for model " + 
					    		boost::lexical_cast<std::string>(i);
					    EP_EXIT_ERROR(message,EPFAIL);
					}
				}
				gsl_matrix_free(matrixaux);
				obj.iniRow = energyInLibrary_row+1;
				obj.endRow = energyInLibrary_row+1;

				if (writeFitsComplex(obj, matchedfilters))
				{
				    message = "Cannot run routine writeFitsComplex to write MATCHEDF if already present";
				    EP_EXIT_ERROR(message,EPFAIL);
				}
			}
		}
		else if (energyInLibrary_row == -1)
		{
			if (writeLib(energy, matchedfilters))
			{
			    message = "Cannot run routine writeLib when energyInLibrary_row == -1";
			    EP_EXIT_ERROR(message,EPFAIL);
			}
		}

		// Close pulses templates library FITS file (LIBRARY extension)
		if (fits_close_file(inLibObject,&status))
		{
	  	    message = "Cannot close file " + string(inLibName);
		    EP_EXIT_ERROR(message,status);	  
		}

		// Calculate the optimal filter

		if (calculus_optimalFilter (filtergsl, eventsz, &optimalfiltergsl, &foptimalfiltergslaux, &optimalfilterFFTgslaux))
		{
		    message = "Cannot run routine calculus_optimalFilter";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		// Write the output FITS file (_flt.fits)
		f_optimal_filtergsl = gsl_vector_alloc(foptimalfiltergslaux->size/2);
		optimal_filterFFTgsl = gsl_vector_alloc(foptimalfiltergslaux->size/2);
		for (int i=0;i<foptimalfiltergslaux->size/2;i++)
		{
		    gsl_vector_set(f_optimal_filtergsl,i,gsl_vector_get(foptimalfiltergslaux,i));
		    gsl_vector_set(optimal_filterFFTgsl,i,gsl_vector_get(optimalfilterFFTgslaux,i));
		}
		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(fltObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) + " in " + string(fltName) + " to write columns";
		    EP_EXIT_ERROR(message,status);
		}

		// OptimalF column
		obj.inObject = fltObject;
		obj.nameTable = new char [255];
		strcpy(obj.nameTable,"FILTER");
		obj.iniRow = 1;
		obj.endRow = optimalfiltergsl->size;
		obj.iniCol = 0;
		obj.nameCol = new char [255];
		strcpy(obj.nameCol,"OptimalF");
		obj.type = TDOUBLE;
		obj.unit = new char [255];
		strcpy(obj.unit,"--");
		if (writeFitsSimple (obj,optimalfiltergsl))
		{
		    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		// Freq column
		obj.endRow = foptimalfiltergslaux->size/2;
		strcpy(obj.nameCol,"Freq");
		strcpy(obj.unit,"Hz");
		if (writeFitsSimple (obj,f_optimal_filtergsl))
		{
		    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		// OptimalFF column
		strcpy(obj.nameCol,"OptimalFF");
		strcpy(obj.unit,"--");
		if (writeFitsSimple (obj,optimal_filterFFTgsl))
		{
		    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		// Write output keywords
		strcpy(keyname,"EVENTCNT");
		keyvalint = optimalfiltergsl->size;
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
		if (fits_write_key(fltObject,TDOUBLE,keyname,&normalizationFactor,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(fltName);
		    EP_EXIT_ERROR(message,status);	  
		}
	    
		// Close output FITS files
		if (fits_close_file(fltObject,&status))
		{
		    message = "Cannot close file " + string(fltName);
		    EP_EXIT_ERROR(message,status);	  
		}

		// Free allocate of GSL vectors
		gsl_vector_free(filtergsl);
		gsl_vector_free(optimalfiltergsl);
		gsl_vector_free(f_optimal_filtergsl);
		gsl_vector_free(optimal_filterFFTgsl);
	}
	else if (mode == 1)		// Production mode
	{
		extern int inDataIteratorTrg(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols,void *user_strct );

		// Create structure to run Iteration
		iteratorCol cols_trg1 [3];	// Structure of Iteration
		int n_cols_trg1 = 3; 		// Number of columns: EstEnrgy, Quality and Grade
		int rows_per_loop = 0;		// 0: Use default: Optimum number of rows
		int offset = 0;			    // 0: Process all the rows

			// Read Columns
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(trgName) + " to get columns structure";
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"EstEnrgy");
		status = fits_iter_set_by_name(&cols_trg1[0], trgObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"Quality");
		status = fits_iter_set_by_name(&cols_trg1[1], trgObject, straux, TSHORT, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"Grade");
		status = fits_iter_set_by_name(&cols_trg1[2], trgObject, straux, TSHORT, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}

		// Called iteration function
		ntotalrows = 0;
		if (fits_iterate_data(n_cols_trg1,cols_trg1,offset,rows_per_loop,inDataIteratorTrg,0L,&status))
		{
		    message = "Cannot iterate data through inDataIteratorTrg ";
		    EP_EXIT_ERROR(message,status);	
		}

		// Check Quality: If there are no valid pulses in the TRIGGER FITS file => The task finishes
		iter = 0;
		for (int i=0; i<eventcnt ;i++)
		{
			if (gsl_vector_get(qualitygsl,i) != 0)
			{
				iter++;
			}
		}
		if (iter == eventcnt)
		{
			message = "There are no valid pulses (quality == 0) in the input FITS file " + string(trgName) + ". There is no output FITS file ";
			writeLog(fileRef,"Warning", verbosity,message);
			EP_EXIT_ERROR(message,EPFAIL);
		}

		// Check Hp, Mp, Ms, Lp, Ls and Grade: If there are no valid pulses in the input FITS file => The task finishes
		iter = 0;
		iter1 = 0;
		iter2 = 0;
		iter3 = 0;
		iter4 = 0;
		iter5 = 0;
		for (int i=0; i<eventcnt ;i++)
		{
			if ((Hp == 1) && (gsl_vector_get(gradegsl,i) == 1))	    	iter1++;
			else if ((Mp == 1) && (gsl_vector_get(gradegsl,i) == 21))	iter2++;
			else if ((Ms == 1) && (gsl_vector_get(gradegsl,i) == 22))	iter3++;
			else if ((Lp == 1) && (gsl_vector_get(gradegsl,i) == 31))	iter4++;
			else if ((Ls == 1) && (gsl_vector_get(gradegsl,i) == 32))	iter5++;
		}
		iter = iter1 + iter2 + iter3 + iter4 + iter5;
		if (iter == 0)
		{
			message = "There are no valid pulses (due to GRADE) in the TRIGGER FITS file " + string(trgName) + ". There is no output FITS file ";
			writeLog(fileRef,"Warning", verbosity,message);
			EP_EXIT_ERROR(message,EPFAIL);
		}

		// Operate to create the optimal filter
		// Declare variables
		double gradei;				// In order to validate pulses according its grade
		// In order to write in FILTER (_trg.fits) matrix are necessary
		gsl_matrix *optimalfiltergsl_matrix;
		gsl_matrix *f_optimal_filtergsl_matrix;
		gsl_matrix *optimal_filterFFTgsl_matrix;
		gsl_matrix *optimalfiltergsl_matrix0;		// *_matrix0 variables necessaries to write non valid pulses at the beginning (all set to 0)
		gsl_matrix *f_optimal_filtergsl_matrix0;
		gsl_matrix *optimal_filterFFTgsl_matrix0;
		gsl_vector *validPulses = gsl_vector_alloc(eventcnt);
		gsl_vector_set_all(validPulses,0);

		// Create the extension FILTER in _trg.fits file
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(trgName) + " to add a new table FILTER";
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(extname,"FILTER");
		if (fits_create_tbl(trgObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
		    message = "Cannot crate table " + string(extname) + " in file " + string(trgName);
		    EP_EXIT_ERROR(message,status); 
		}

		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}

		for (int i=0; i<eventcnt ;i++)
		{
			sprintf(val,"-------------> Pulse: %d of %d<------------------ ",i+1, eventcnt);
			writeLog(fileRef,"Log", verbosity,string(val));
			strcat(val,"\n");
			fputs(val,temporalFile);

			filtergsl = gsl_vector_alloc(matchedfilters->size2);			// Filter values
			gsl_vector_set_all(filtergsl,0);

			gradei = gsl_vector_get(gradegsl,i);

			// Pulses are going to be validated by checking its quality and grade
			if (gsl_vector_get(qualitygsl,i) == 0)	// Neither truncated or saturated
			{
				if (((gradei == 1) && (Hp == 1)) || ((gradei == 21) && (Mp == 1)) || ((gradei == 22) && (Ms == 1)) || ((gradei == 31) && (Lp == 1)) || ((gradei == 32) && (Ls == 1)))
				{
					nValidPulses++;
					gsl_vector_set(validPulses,i,1);
					if (alreadyValidPulse == 0)		alreadyValidPulse = 1;

					// Obtain the matched filter by interpolating
					if (find_matchedfilter (gsl_vector_get(energygsl,i), estenergylibrary, matchedfilters, &filtergsl, temporalFile))
					{
					      message = "Cannot run routine find_matchedfilter for filter interpolation";
					      EP_EXIT_ERROR(message,EPFAIL);
					}

					if (area0 (&filtergsl))		// Just in case the interpolation modifies the 0 area
					{
					      message = "Cannot run routine area0";
					      EP_EXIT_ERROR(message,EPFAIL);
					}

					// Calculate the optimal filter
					if (calculus_optimalFilter (filtergsl, filtergsl->size, &optimalfiltergsl, &foptimalfiltergslaux, &optimalfilterFFTgslaux))
					{
					      message = "Cannot run routine calculus_optimalFilter";
					      EP_EXIT_ERROR(message,EPFAIL);
					}

					optimalfiltergsl_matrix = gsl_matrix_alloc(1,matchedfilters->size2);
					gsl_vector *optimalfiltergsl_matrix_row = gsl_vector_alloc(matchedfilters->size2);
					gsl_vector_set_zero(optimalfiltergsl_matrix_row);
					for (int j=0; j<optimalfiltergsl->size; j++)
					{
						gsl_vector_set(optimalfiltergsl_matrix_row,j,gsl_vector_get(optimalfiltergsl,j));
					}
					gsl_matrix_set_row(optimalfiltergsl_matrix,0,optimalfiltergsl_matrix_row);
					gsl_vector_free(optimalfiltergsl_matrix_row);
					sizeoptimalfiltergsl = matchedfilters->size2;

					// Write the FILTER in _trg FITS file
					f_optimal_filtergsl = gsl_vector_alloc(matchedfilters->size2/2);
					optimal_filterFFTgsl = gsl_vector_alloc(matchedfilters->size2/2);
					gsl_vector_set_zero(f_optimal_filtergsl);
					gsl_vector_set_zero(optimal_filterFFTgsl);
					for (int ii=0;ii<filtergsl->size/2;ii++)
					{
						gsl_vector_set(f_optimal_filtergsl,ii,gsl_vector_get(foptimalfiltergslaux,ii));
						gsl_vector_set(optimal_filterFFTgsl,ii,gsl_vector_get(optimalfilterFFTgslaux,ii));
					}
					f_optimal_filtergsl_matrix = gsl_matrix_alloc (1,f_optimal_filtergsl->size);
					optimal_filterFFTgsl_matrix = gsl_matrix_alloc (1,f_optimal_filtergsl->size);
					gsl_matrix_set_row(f_optimal_filtergsl_matrix,0,f_optimal_filtergsl);
					gsl_matrix_set_row(optimal_filterFFTgsl_matrix,0,optimal_filterFFTgsl);
					sizeoptimalfilterFFTgsl = optimal_filterFFTgsl->size;

					obj.inObject = trgObject;
					obj.nameTable = new char [255];
					strcpy(obj.nameTable,"FILTER");
					obj.iniCol = 0;
					obj.nameCol = new char [255];
					obj.type = TDOUBLE;
					obj.unit = new char [255];

					if (nonValidPulsesBeginning != 0)
					{
						optimalfiltergsl_matrix0 = gsl_matrix_alloc(1,sizeoptimalfiltergsl);
						f_optimal_filtergsl_matrix0 = gsl_matrix_alloc(1,sizeoptimalfilterFFTgsl);
						optimal_filterFFTgsl_matrix0 = gsl_matrix_alloc(1,sizeoptimalfilterFFTgsl);
						gsl_vector_set(nrmfctrgsl,0,0);

						gsl_matrix_set_all(optimalfiltergsl_matrix0,0);
						gsl_matrix_set_all(f_optimal_filtergsl_matrix0,0);
						gsl_matrix_set_all(optimal_filterFFTgsl_matrix0,0);

						for (int j=0;j<nonValidPulsesBeginning;j++)
						{
							obj.iniRow = j+1;
							obj.endRow = j+1;

							// OptimalF column
							strcpy(obj.nameCol,"OptimalF");
							strcpy(obj.unit,"--");
							if (writeFitsComplex (obj,optimalfiltergsl_matrix0))
							{
							    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) +
							    		" column in FILTER (non-valid pulse number " + boost::lexical_cast<std::string>(j) + ")";
							    EP_EXIT_ERROR(message,EPFAIL);
							}

							// NRMFCTR column
							strcpy(obj.nameCol,"NRMFCTR");
							if (writeFitsSimple (obj,nrmfctrgsl))
							{
							    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) +
							    		" column in FILTER (non-valid pulse number " + boost::lexical_cast<std::string>(j) + ")";
							    EP_EXIT_ERROR(message,EPFAIL);
							}

							// Freq column
							strcpy(obj.nameCol,"Freq");
							strcpy(obj.unit,"Hz");
							if (writeFitsComplex (obj,f_optimal_filtergsl_matrix0))
							{
							    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) +
							    " column in FILTER (non-valid pulse number " + boost::lexical_cast<std::string>(j) + ")";
							    EP_EXIT_ERROR(message,EPFAIL);
							}

							// OptimalFF column
							strcpy(obj.nameCol,"OptimalFF");
							strcpy(obj.unit,"--");
							if (writeFitsComplex (obj,optimal_filterFFTgsl_matrix0))
							{
							    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) +
							    		" column in FILTER (non-valid pulse number " + boost::lexical_cast<std::string>(j) + ")";
							    EP_EXIT_ERROR(message,EPFAIL);
							}
						}

						nonValidPulsesBeginning = 0;

						gsl_matrix_free(optimalfiltergsl_matrix0);
						gsl_matrix_free(f_optimal_filtergsl_matrix0);
						gsl_matrix_free(optimal_filterFFTgsl_matrix0);
					}

					// OptimalF column
					obj.iniRow = i+1;
					obj.endRow = i+1;
					strcpy(obj.nameCol,"OptimalF");
					strcpy(obj.unit,"--");
					if (writeFitsComplex (obj,optimalfiltergsl_matrix))
					{
					    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					// NRMFCTR column
					gsl_vector_set(nrmfctrgsl,0,normalizationFactor);
					strcpy(obj.nameCol,"NRMFCTR");
					if (writeFitsSimple (obj,nrmfctrgsl))
					{
					    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in FILTER";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					// Freq column
					strcpy(obj.nameCol,"Freq");
					strcpy(obj.unit,"Hz");
					if (writeFitsComplex (obj,f_optimal_filtergsl_matrix))
					{
					    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					// OptimalFF column
					strcpy(obj.nameCol,"OptimalFF");
					strcpy(obj.unit,"--");
					if (writeFitsComplex (obj,optimal_filterFFTgsl_matrix))
					{
					    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					// Free allocate
					gsl_vector_free(optimalfiltergsl);
					gsl_matrix_free(optimalfiltergsl_matrix);
					gsl_vector_free(f_optimal_filtergsl);
					gsl_matrix_free(f_optimal_filtergsl_matrix);
					gsl_vector_free(optimal_filterFFTgsl);
					gsl_matrix_free(optimal_filterFFTgsl_matrix);
				}
			}

			//If a pulse is not a valid one (and in this case is not at the beginning), columns in FILTER in the corresponding row are filled with 0's'
			if (gsl_vector_get(validPulses,i) == 0)
			{
				if (alreadyValidPulse == 1)
				{
					optimalfiltergsl_matrix = gsl_matrix_alloc(1,sizeoptimalfiltergsl);
					f_optimal_filtergsl_matrix = gsl_matrix_alloc(1,sizeoptimalfilterFFTgsl);
					optimal_filterFFTgsl_matrix = gsl_matrix_alloc(1,sizeoptimalfilterFFTgsl);

					gsl_matrix_set_all(optimalfiltergsl_matrix,0);
					gsl_matrix_set_all(f_optimal_filtergsl_matrix,0);
					gsl_matrix_set_all(optimal_filterFFTgsl_matrix,0);

					obj.inObject = trgObject;
					obj.nameTable = new char [255];
					strcpy(obj.nameTable,"FILTER");
					obj.iniCol = 0;
					obj.nameCol = new char [255];
					obj.type = TDOUBLE;
					obj.unit = new char [255];

					// OptimalF column
					obj.iniRow = i+1;
					obj.endRow = i+1;
					strcpy(obj.nameCol,"OptimalF");
					strcpy(obj.unit,"--");
					if (writeFitsComplex (obj,optimalfiltergsl_matrix))
					{
					    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER (valid pulses = 0)";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					// NRMFCTR column
					gsl_vector_set(nrmfctrgsl,0,0);
					strcpy(obj.nameCol,"NRMFCTR");
					if (writeFitsSimple (obj,nrmfctrgsl))
					{
					    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " column in FILTER (valid pulses = 0)";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					// Freq column
					strcpy(obj.nameCol,"Freq");
					strcpy(obj.unit,"Hz");
					if (writeFitsComplex (obj,f_optimal_filtergsl_matrix))
					{
					    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER (valid pulses = 0)";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					// OptimalFF column
					strcpy(obj.nameCol,"OptimalFF");
					strcpy(obj.unit,"--");
					if (writeFitsComplex (obj,optimal_filterFFTgsl_matrix))
					{
					    message = "Cannot run routine writeFitsComplex to write " + string(obj.nameCol) + " column in FILTER (valid pulses = 0)";
					    EP_EXIT_ERROR(message,EPFAIL);
					}

					gsl_matrix_free(optimalfiltergsl_matrix);
					gsl_matrix_free(f_optimal_filtergsl_matrix);
					gsl_matrix_free(optimal_filterFFTgsl_matrix);
				}
				else
				{
					nonValidPulsesBeginning++;
				}
			}

			gsl_vector_free(filtergsl);
		}

		//string annals (string("and by")+ ' ' +(string) create);
		//string annalsall ((string) annalsin + ' ' + annals);
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname) + " in file " + string(trgName);
			EP_PRINT_ERROR(message,status);return(EPFAIL);
		}
		string mod1 (string("File MODIFIED by") + ' ' +	(string) create);

		strcpy(keyname,"MOD1");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(trgName);
		    EP_EXIT_ERROR(message,status);	  
		}
		if (fits_update_key_longstr(trgObject,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(trgName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	// Close input FITS files
	if (fits_close_file(trgObject,&status))
	{
	    message = "Cannot close file " + string(trgName);
	    EP_EXIT_ERROR(message,status);	  
	}
	if (fits_close_file(noisespecObject,&status))
	{
	    message = "Cannot close file " + string(noisespecName);
	    EP_EXIT_ERROR(message,status);	  
	}

	// Free allocate of memory
	gsl_vector_free(nrmfctrgsl);
	//delete [] straux;
	/*for (int i=0;i<khNumber;i++){
		delete [] kh[i].limitMin;
		delete [] kh[i].limitMax;
	}
	delete [] kh;*/
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	// To delete
	if (fclose(temporalFile))
	{
	    message = "Cannot close temporalFile " + string(temporalFileName);
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	
	// Finalize the task
	time_t t_end = time(NULL);
	sprintf(straux,"%f",(double) (t_end - t_start));
	message = "Time:" + string(straux);
	writeLog(fileRef,"Log", verbosity,message);
	
	writeLog(fileRef,"OK", verbosity,"Filter Module OK");

	if (fclose(fileRef))
	{
	    message = "Cannot close log file ";
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	
	return EPOK;
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 *************************************************************
* initModule function: This function reads and processes the input parameters
*
* - It defines the input parameters and assigns the values to variables
* - It defines the structure for the command line options
* - It gets the parameter values from the command line or user input
* - It saves the parameter values into meaningful variables
*   (checking also if the parameter value is in the allowed range)
* - It opens the log file of the task
****************************************************************************/
int initModule(int argc, char **argv)
{
	int status = EPOK;
	
	// Define TRIGGER input parameters and assign values to variables
	// Parameter definition and assignation of default values
	const int npars = 12, npars1 = 13;
	inparam filterPars[npars];
	int optidx =0, par=0, fst=0, ipar; 
	string message, task="filter";
	
	filterPars[0].name = "trgFile";
    filterPars[0].description = "TRG fits filename";
    filterPars[0].defValStr = "a_trg.fits";
    filterPars[0].type =  "char";
	filterPars[0].ValStr = filterPars[0].defValStr;

	filterPars[1].name = "noisespecFile";
    filterPars[1].description = "NOISESPEC fits filename";
    filterPars[1].defValStr = "a_noisespec.fits";
    filterPars[1].type =  "char";
	filterPars[1].ValStr = filterPars[1].defValStr;
	
	filterPars[2].name = "inLibFile";
    filterPars[2].description = "Pulse models library file name";
    filterPars[2].defValStr = "library.fits";
    filterPars[2].type =  "char";
	filterPars[2].ValStr = filterPars[2].defValStr;
	
	filterPars[3].name = "fltFile"; 
    filterPars[3].description = "Output file name (only for calibration mode)";
    filterPars[3].defValStr = "a_flt.fits";
    filterPars[3].type =  "char";
	filterPars[3].ValStr = filterPars[3].defValStr;

	filterPars[4].name = "Hp"; 
    filterPars[4].description = "Include High-Res Primary pulses in optimal filter calculation? (1:TRUE, 0:FALSE)";
    filterPars[4].defValInt = 1;
    filterPars[4].type = "int";
    filterPars[4].minValInt = 0;
    filterPars[4].maxValInt = 1;
	filterPars[4].ValInt = filterPars[4].defValInt;
	
	filterPars[5].name = "Mp"; 
    filterPars[5].description = "Include Mid-Res Primary pulses in optimal filter calculation? (1:TRUE, 0:FALSE)";
    filterPars[5].defValInt = 1;
    filterPars[5].type = "int";
    filterPars[5].minValInt = 0;
    filterPars[5].maxValInt = 1;
	filterPars[5].ValInt = filterPars[5].defValInt;
	
	filterPars[6].name = "Ms"; 
    filterPars[6].description = "Include Mid-Res Secondary pulses in optimal filter calculation? (1:TRUE, 0:FALSE)";
    filterPars[6].defValInt = 1;
    filterPars[6].type = "int";
    filterPars[6].minValInt = 0;
    filterPars[6].maxValInt = 1;    
	filterPars[6].ValInt = filterPars[6].defValInt;
	
	filterPars[7].name = "Lp"; 
    filterPars[7].description = "Include Low-Res Primary pulses in optimal filter calculation? (1:TRUE, 0:FALSE)";
    filterPars[7].defValInt = 1;
    filterPars[7].type = "int";
    filterPars[7].minValInt = 0;
    filterPars[7].maxValInt = 1;    
	filterPars[7].ValInt = filterPars[7].defValInt;
	
	filterPars[8].name = "Ls"; 
    filterPars[8].description = "Include Low-Res Secondary pulses in optimal filter calculation? (1:TRUE, 0:FALSE)";
    filterPars[8].defValInt = 1;
    filterPars[8].type = "int";
    filterPars[8].minValInt = 0;
    filterPars[8].maxValInt = 1;    
	filterPars[8].ValInt = filterPars[8].defValInt;

	filterPars[9].name = "nameLog";
    filterPars[9].description = "Output log file name";
    filterPars[9].defValStr = "flt_log.txt";
    filterPars[9].type = "char";
	filterPars[9].ValStr = filterPars[9].defValStr;
		
	filterPars[10].name = "verbosity";
    filterPars[10].description = "Verbosity level of the output log file (in [0,3])";
    filterPars[10].defValInt = 3;
    filterPars[10].type = "int";
    filterPars[10].minValInt = 0;
    filterPars[10].maxValInt = 3;
	filterPars[10].ValInt = filterPars[10].defValInt;
	
	filterPars[11].name = "clobber";
    filterPars[11].description = "Re-write output files if clobber=yes";
    filterPars[11].defValStr = "no";
    filterPars[11].type = "char";
	filterPars[11].ValStr = filterPars[11].defValStr;
	
	// Define structure for command line options
	static struct option long_options[npars1];
	for (optidx = 0; optidx<npars; optidx++)
	{
		long_options[optidx].name= filterPars[optidx].name.c_str();
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
			for (int i=0;i<npars; i++)
			{
				if (long_options[optidx].name == filterPars[i].name.c_str())
				{
				    if (filterPars[i].type == "char") //save char value for par
				    {
				    	filterPars[i].ValStr = optarg;
				    }
				    else // check if numeric value
				    {
				    	if ((!isdigit(optarg[0]) && (optarg[0] != '-')) ||
				    			(!isdigit(optarg[0]) && (optarg[0] == '-') && (!isdigit(optarg[1]))))
				    	{
				    		message = "Invalid value for input argument " + string(filterPars[i].name);
				    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				    	}
				    	if (filterPars[i].type == "int")
				    	{
				    		filterPars[i].ValInt = atoi(optarg);
				    	}
				    	else
				    	{
				    		filterPars[i].ValReal= atof(optarg);
				    	}
				    }
				    break;
				} // endif
			} // endfor
			break;
		    default:
		    	message = "Invalid parameter name ";
		    	EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}//switch
	}//while

	// If command line is empty: ask for params interactively
	if (commandLine == 0)
	{
		if (interactivePars(filterPars,npars,task))
		{
		    message = "Error reading parameters interactively";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL); 
		}
	}

	for(int i=0;i<npars; i++)
	{
		if(filterPars[i].name == "trgFile")
		{
			strcpy(trgName, filterPars[i].ValStr.c_str());
		}
		else if(filterPars[i].name == "noisespecFile")
		{
			strcpy(noisespecName, filterPars[i].ValStr.c_str());
		}
		else if(filterPars[i].name == "inLibFile")
		{
			strcpy(inLibName, filterPars[i].ValStr.c_str());
		}
		else if(filterPars[i].name == "fltFile")
		{
			strcpy(fltName, filterPars[i].ValStr.c_str());
		}
		// Parameters to consider a pulse as a valid one
		else if(filterPars[i].name == "Hp")
		{
			Hp = filterPars[i].ValInt;
		}
		else if(filterPars[i].name == "Mp")
		{
			Mp = filterPars[i].ValInt;
		}
		else if(filterPars[i].name == "Ms")
		{
			Ms = filterPars[i].ValInt;
		}
		else if(filterPars[i].name == "Lp")
		{
			Lp = filterPars[i].ValInt;
		}
		else if(filterPars[i].name == "Ls")
		{
			Ls = filterPars[i].ValInt;

		}
		else if(filterPars[i].name == "nameLog")
		{
			strcpy(nameLog,filterPars[i].ValStr.c_str());
		}
		else if(filterPars[i].name == "verbosity")
		{
			verbosity = filterPars[i].ValInt;
		}
		else if (filterPars[i].name == "clobber")
		{
			strcpy(clobberStr, filterPars[i].ValStr.c_str());
			if(strcmp(clobberStr,"yes")==0){
			  clobber=1;
			}else{
			  clobber=0;
			}
		}
		
		// Check if parameter value is in allowed range
		if (filterPars[i].type == "int" &&
				(filterPars[i].ValInt < filterPars[i].minValInt ||
						filterPars[i].ValInt > filterPars[i].maxValInt))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(filterPars[i].minValInt) + "," + boost::lexical_cast<std::string>(filterPars[i].maxValInt) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		else if (filterPars[i].type == "double" &&
				(filterPars[i].ValReal < filterPars[i].minValReal ||
						filterPars[i].ValReal > filterPars[i].maxValReal))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(filterPars[i].minValReal) + "," + boost::lexical_cast<std::string>(filterPars[i].maxValReal) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
	}// end loop for parameters

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
* readLib: This function reads the pulse templates library input FITS file
*
* If calibration mode and the column MATCHEDF column does not exist => inDataIteratorLib0
* If production mode and the column MATCHEDF column does not exist => Error
* If production mode and the column MATCHEDF column does exist => inDataIteratorLib1
******************************************************************************/
int readLib ()
{	
	string message = "";
	int status = EPOK, extver=0;
	
	// Open pulse templates library file (LIBRARY extension)
	if (fits_open_file(&inLibObject, inLibName,READWRITE,&status))
	{
  	      message = "Input library file (" + string(inLibName) + ") cannot be open ";
          writeLog(fileRef,"Error", verbosity,message);
	      EP_PRINT_ERROR(message,status); return(EPFAIL);  
	}
	strcpy(extname,"LIBRARY");
	if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_get_num_rows(inLibObject,&nummodels, &status))
	{
	    message = "Cannot access HDU " + string(extname) + " in " + string(inLibName) + " file (cannot get number of rows)";
	    writeLog(fileRef,"Error", verbosity,message);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Get structure of pulse templates library input FITS file columns
	strcpy(straux,"Energy");
	if (fits_get_colnum(inLibObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get number for column " + string(straux) + " in " + string(inLibName) + " file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(straux,"EstEnergy");
	if (fits_get_colnum(inLibObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get number for column " + string(straux) + " in " + string(inLibName) + " file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(straux,"Pulse");
	if (fits_get_colnum(inLibObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get number for column " + string(straux) + " in " + string(inLibName) + " file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(straux,"MATCHEDF");
	if (fits_get_colnum(inLibObject,0,straux,&colnum,&status))
	{
	    if (mode == 0)
	    {
	    	matchedf_exist = 0;
	    	status = 0;
	    }
	    else if (mode == 1)
	    {
	    	message = "Cannot get number for column " + string(straux) + " in " + string(inLibName) + " file";
	    	EP_PRINT_ERROR(message,status); return(EPFAIL);
	    }
	}

	// When matchedf_exist = 0
	extern int inDataIteratorLib0(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);
	// When matchedf_exist = 1
	extern int inDataIteratorLib1(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

	// Create structure to run Iteration
	long rows_per_loop = nummodels; 		// 0: Use default: Optimum number of rows
	long offset=0; 							// 0: Process all the rows

	if (matchedf_exist == 1)
	{
		iteratorCol colsLib1 [4]; 			// Structure of Iteration
		int n_colsLib1 = 4; 				// Number of columns:  Energy + EstEnergy + Pulse + MATCHEDF

		// Read Energy Column
		strcpy(straux,"Energy");
		status = fits_iter_set_by_name(&colsLib1[0], inLibObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate by name column " + string(straux) + " when MATCHEDF column is already present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Read EstEnergy Column
		strcpy(straux,"EstEnergy");
		status = fits_iter_set_by_name(&colsLib1[1], inLibObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate by name column " + string(straux) + " when MATCHEDF column is already present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Read Pulse Column
		strcpy(straux,"Pulse");
		status = fits_iter_set_by_name(&colsLib1[2], inLibObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate by name column " + string(straux) + " when MATCHEDF column is already present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Read MATCHEDF
		strcpy(straux,"MATCHEDF");
		status = fits_iter_set_by_name(&colsLib1[3], inLibObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate by name column " + string(straux) + " when MATCHEDF column is already present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Called iteration function: inDataIteratorLib1
		if (fits_iterate_data(n_colsLib1, colsLib1, offset, rows_per_loop, inDataIteratorLib1,0L,&status))
		{
		    message = "Cannot iterate data in columns when MATCHEDF column is already present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}
	else if (matchedf_exist == 0)
	{
		iteratorCol colsLib0 [3]; 				// Structure of Iteration DAL
		int n_colsLib0 = 3; 					// Number of columns:  Energy + EstEnergy + PULSE

		// Read Energy Column
		strcpy(straux,"Energy");
		status = fits_iter_set_by_name(&colsLib0[0], inLibObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate by name column " + string(straux) + " when MATCHEDF column is not present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(straux,"EstEnergy");
		status = fits_iter_set_by_name(&colsLib0[1], inLibObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate by name column " + string(straux) + " when MATCHEDF column is not present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(straux,"Pulse");
		status = fits_iter_set_by_name(&colsLib0[2], inLibObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate by name column " + string(straux) + " when MATCHEDF column is not present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Called iteration function: inDataIteratorLib0
		if (fits_iterate_data(n_colsLib0, colsLib0, offset, rows_per_loop, inDataIteratorLib0,0L,&status))
		{
		    message = "Cannot iterate data in columns when MATCHEDF column is not present in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return (EPOK);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* inDataIteratorLib0: This function runs in calibration mode when the MATCHEDF column does not exist in the library
*                     FITS file. It takes the optimum number of rows to read the pulse templates library and works iteratively.
*                     It reads ENERGY, ESTENERGY and PULSE columns from the pulse templates library and stores them in
*                     'energylibrary', 'estenergylibrary' and 'templateslibrary'
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each row of the pulse templates library
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorLib0 (long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	int status = EPOK, extver=0;
	string message = "";
	  
	// Declare variables
	double *energy, *energyin;				// Vector of ENERGY column
	double *estenergy, *estenergyin;		// Vector of ESTENERGY column
	double *pulse, *pulsein;				// Vector of PULSE column

	// Allocate input GSL vectors
	gsl_vector *energygsl_i = gsl_vector_alloc(nrows);
	gsl_vector *estenergygsl_i = gsl_vector_alloc(nrows);
	gsl_matrix *pulsegsl = gsl_matrix_alloc(nrows, eventsz);
	energylibrary = gsl_vector_alloc(nummodels);					// Energy
	estenergylibrary = gsl_vector_alloc(nummodels);					// EstEnergy
	templateslibrary = gsl_matrix_alloc(nummodels,eventsz);			// Pulse

	gsl_vector *rowaux = gsl_vector_alloc(eventsz);

	strcpy(extname,"LIBRARY");
	if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU  " + string(extname) + " in library file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"EVENTSZ");
	if (fits_read_key(inLibObject,TLONG,keyname, &eventsz,comment,&status))
	{
	    message = "Cannot read keyword " + string(keyname) + " in library file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	// Read iterator
	energyin    = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	energy = &energyin[1];
	if (toGslVector(((void **)&energy), &energygsl_i, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert ENERGY column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	estenergyin = (double *) fits_iter_get_array(&cols[1]);
	estenergy = &estenergyin[1];
	if (toGslVector(((void **)&estenergy), &estenergygsl_i, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert ESTENERGY column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	pulsein = (double *) fits_iter_get_array(&cols[2]);
	pulse = &pulsein[1];
	if (toGslMatrix(((void **)&pulse), &pulsegsl, eventsz, nrows, (int)TDOUBLE,0))
	{
	    message = "Cannot convert PULSE column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Processing each row of the pulse templates library
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set (energylibrary, ntotalrows, gsl_vector_get(energygsl_i, i));
		gsl_vector_set (estenergylibrary, ntotalrows, gsl_vector_get(estenergygsl_i, i));

		gsl_matrix_get_row(rowaux,pulsegsl,i);
		gsl_matrix_set_row(templateslibrary,ntotalrows,rowaux);

		ntotalrows ++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free(energygsl_i);
	gsl_vector_free(estenergygsl_i);
	gsl_matrix_free(pulsegsl);

	return (EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* inDataIteratorLib1: This function runs in production mode when the MATCHEDF column exists in the library FITS file.
*                     It takes the optimum number of rows to read the pulse templates library and works iteratively.
*                     It reads ENERGY, ESTENERGY, PULSE and MATCHEDF columns from the pulse templates library and stores them in
*                     'energylibrary', 'estenergylibrary', 'templateslibrary' and 'matchedfilters'
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each row of the pulse templates library
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorLib1 (long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	int status = EPOK, extver=0;
	string message = "";
  
	// Declare variables
	double *energy, *energyin;							// Vector of ENERGY column
	double *estenergy, *estenergyin;					// Vector of EstENERGY column
	double *pulse, *pulsein;							// Vector of PULSE column
	double *matchedfiltermodel, *matchedfiltermodelin;  // Vector of MATCHEDF column

	// Allocate input GSL vectors
	gsl_vector *energygsl_i = gsl_vector_alloc(nrows);
	gsl_vector *estenergygsl_i = gsl_vector_alloc(nrows);
	gsl_matrix *pulsegsl = gsl_matrix_alloc(nrows, eventsz);
	gsl_matrix *matchedfiltergsl = gsl_matrix_alloc(nummodels,eventsz);
	energylibrary = gsl_vector_alloc(nummodels);			// Energy
	estenergylibrary = gsl_vector_alloc(nummodels);			// EstEnergy
	templateslibrary = gsl_matrix_alloc(nummodels,eventsz);	// Pulse
	matchedfilters = gsl_matrix_alloc(nummodels,eventsz);
	gsl_vector *rowaux = gsl_vector_alloc(eventsz);			// Auxiliary local variable
	matchedfilter = gsl_vector_alloc(eventsz);

	// Read iterator
	energyin = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	energy = &energyin[1];
	if (toGslVector(((void **)&energy), &energygsl_i, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert ENERGY column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	estenergyin = (double *) fits_iter_get_array(&cols[1]);
	estenergy = &estenergyin[1];
	if (toGslVector(((void **)&estenergy), &estenergygsl_i, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert ESTENERGY column toGslVector";	  
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	pulsein = (double *) fits_iter_get_array(&cols[2]);
	pulse = &pulsein[1];
	if (toGslMatrix(((void **)&pulse), &pulsegsl, eventsz, nrows, TDOUBLE,0))
	{
	    message = "Cannot convert PULSE column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	matchedfiltermodelin = (double *) fits_iter_get_array(&cols[3]);
	matchedfiltermodel = &matchedfiltermodelin[1];
	if (toGslMatrix(((void **)&matchedfiltermodel), &matchedfiltergsl, eventsz, nrows, TDOUBLE,0))
	{
	    message = "Cannot convert MATCHEDF column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Processing each row of pulses
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set (energylibrary, ntotalrows, gsl_vector_get(energygsl_i, i));
		gsl_vector_set (estenergylibrary, ntotalrows, gsl_vector_get(estenergygsl_i, i));

		gsl_matrix_get_row(rowaux,pulsegsl,i);
		gsl_matrix_set_row(templateslibrary,ntotalrows,rowaux);

		gsl_matrix_get_row(rowaux,matchedfiltergsl,i);
		gsl_matrix_set_row(matchedfilters,ntotalrows,rowaux);

		ntotalrows ++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free(energygsl_i);
	gsl_vector_free(estenergygsl_i);
	gsl_matrix_free(pulsegsl);
	gsl_matrix_free(matchedfiltergsl);
	gsl_vector_free(rowaux);

	return (EPOK);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* createFilterFile function: In calibration mode this function creates the structure of the output _flt FITS file
*                            In production mode the PROCESS keyword in _trg FITS file is updated
***************************************************************************/
int createFilterFile ()
{
	int status = EPOK, extver=0;
	string message = "";
	char *tt[1];
	char *tf[1];
	char *tu[1];
	
	if (mode == 0)	//Calibration mode
	{
		// Create output FITS files: If it does not exist yet
		//Create Filter File
		if ( fileExists(string(fltName)) && clobber==1)
		{
		    if (remove(fltName)){
			message = "Filter file " + string(fltName) + " already exists: & cannot be deleted ("+string(strerror(errno))+")";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    }
		}
		else if(fileExists(string(fltName)) && clobber==0)
		{
		    message = "Filter file " + string(fltName) + " already exists: must not be overwritten";
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
		
		message = "Create Filter FITS File: " + string(fltName);
		writeLog(fileRef,"Log", verbosity,message);
		
		// Create extensions: FILTER
		strcpy(extname,"FILTER");
		if (fits_open_file(&fltObject,fltName,READWRITE,&status))
		{
		    message = "Cannot open created output file " + string(fltName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL); 
		}
		if (fits_create_tbl(fltObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
		    message = "Cannot crate table " + string(extname) + " in file " + string(fltName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL); 
		}
		if (fits_movnam_hdu(fltObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) + " in file " + string(fltName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL); 
		}

		// Create keywords
		strcpy(keyname,"CREATOR");
		strcpy(keyvalstr,create);
		if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		      message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
		      EP_PRINT_ERROR(message,status); return(EPFAIL);  
		}
	}

	// Set PROCESS keyword
	char str_verb[125];			sprintf(str_verb,"%d",verbosity);
	char str_Hp[125];			sprintf(str_Hp,"%d",Hp);
	char str_Mp[125];			sprintf(str_Mp,"%d",Mp);
	char str_Ms[125];			sprintf(str_Ms,"%d",Ms);
	char str_Lp[125];			sprintf(str_Lp,"%d",Lp);
	char str_Ls[125];			sprintf(str_Ls,"%d",Ls);

	string processoutFLT (string("FILTER") 	+ ' ' +
	string(trgName)  + ' ' + string(noisespecName) + ' ' + string(inLibName) + ' ' + string(fltName) + ' ' +
	string(str_Hp)   + ' ' + string(str_Mp)        + ' ' + string(str_Ms)  + ' ' +
	string(str_Lp)   + ' ' + string(str_Ls)        + ' ' +
	string(nameLog)  + ' ' + string(str_verb)      + ' ' + string(clobberStr) + ' ' +
	string("(")		 +      (string) create 	   +   	   string(")"));

	if (mode == 0)
	{
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
 		strcpy(keyname,"FTYPE");
		strcpy(keyvalstr,"filter");	
		if (fits_write_key(fltObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
 		strcpy(keyname,"MODE");
		if (fits_write_key(fltObject,TINT,keyname,&mode,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(fltName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}
	else if (mode == 1)
	{
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) + " in file " + string(trgName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL); 
		}

		strcpy(keyname,"PROC2");
		strcpy(keyvalstr,processoutFLT.c_str());
		if (fits_write_key_longwarn(trgObject,&status))
		{
		    message = "Cannot update keyword " + string(keyname) + " in " + string(trgName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);  
		}
		if (fits_update_key_longstr(trgObject,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(trgName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return EPOK;
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* inDataIteratorTPS: This function reads F and CSD columns from the _noisespec input FITS file (NOISEALL extension)
*                    and stores them in 'freqgsl' and 'csdgsl'. It takes the optimum number of rows and works iteratively.
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Storing each row of the current noise spectral density
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorTPS (long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct)
{
	int status = EPOK;
	string message = "";

	/////////////////////// Declare variables
	double *f, *fin;				// Vector of FREQ column
	double *csd, *csdin;			// Vector of CSD column

	// Allocate input GSL vectors:
	gsl_vector *freqgsl_i = gsl_vector_alloc(nrows);	// GSL of FREQ column
	gsl_vector *csdgsl_i = gsl_vector_alloc(nrows);		// GSL of CSD column

	// Read iterator
	// NOTE: fits_iter_get_array: in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	fin = (double *) fits_iter_get_array(&cols[0]);
	f = &fin[1];
	if (toGslVector(((void **)&f),&freqgsl_i,nrows,0,TDOUBLE))
	{
	    message = "Cannot convert FREQ column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	csdin = (double *) fits_iter_get_array(&cols[1]);
	csd = &csdin[1];
	if (toGslVector(((void **)&csd),&csdgsl_i,nrows,0,TDOUBLE))
	{
	    message = "Cannot convert CSD column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	//  Storing each row of current noise spectral density
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set(freqgsl,ntotalrows,gsl_vector_get(freqgsl_i,i));
		gsl_vector_set(csdgsl,ntotalrows,gsl_vector_get(csdgsl_i,i));
		ntotalrows ++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free (freqgsl_i);
	gsl_vector_free (csdgsl_i);

	return(EPOK);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* find_estenergy: This function uses the ENERGY of the pulses (it is always used in calibration mode, when
*                 the pulses are monochromatic) in order to choose the proper row of the pulse templates library
*                 to write the calculated matched filter.
*                 If ENERGY does not match to any energy of the pulse models library, rowFound will be -1.
*
* Parameters:
* - energyKeyword: ENERGY keyword of the pulses whose place (row) in the pulse templates library is looking for
* - energygsl: Vector where the values of Energy of each pulse template are stored
* - rowFound: Row of the pulse templates library where the matched filter will be written
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
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ************************************************************
* inDataIteratorTrg: This function runs when the operation 'mode' is 1 (production).
*                    It takes the optimum number of rows to read the _trg input FITS file and works iteratively
*                    It reads QUALITY, GRADE and ESTENERGY columns from the _trg input FITS file and stores it in
*                    'qualitygsl', 'gradegsl' and 'energygsl'
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each row of pulses
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorTrg (long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	int status = EPOK;
	string message= "";
  
	// Declare variables
	short *quality, *qualityin;		// Vector of QUALITY column
	short *grade, *gradein;			// Vector of GRADE column
	double *energy, *energyin;		// Vector of ESTENRGY column

	// Allocate input GSL vectors:
	gsl_vector *qualitygsl_i = gsl_vector_alloc(nrows);		// GSL of Quality column
	gsl_vector *gradegsl_i = gsl_vector_alloc(nrows);		// GSL of Grade column
	gsl_vector *energygsl_i = gsl_vector_alloc(nrows);		// GSL of EstEnrgy column

	// Read iterator
	energyin = (double *) fits_iter_get_array(&cols[0]);
	energy = &energyin[1];
	if (toGslVector(((void **)&energy),&energygsl_i,nrows,0,TDOUBLE))
	{
	    message = "Cannot convert ESTENERGY column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	qualityin = (short *) fits_iter_get_array(&cols[1]);
	quality = &qualityin[1];
	if (toGslVector(((void **)&quality),&qualitygsl_i,nrows,0,TSHORT))
	{
	    message = "Cannot convert QUALITY column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gradein = (short *) fits_iter_get_array(&cols[2]);
	grade = &gradein[1];
	if (toGslVector(((void **)&grade),&gradegsl_i,nrows,0,TSHORT))
	{
	    message = "Cannot convert GRADE column to GSL vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Processing each found pulse
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set (qualitygsl, ntotalrows, gsl_vector_get(qualitygsl_i, i));
		gsl_vector_set (gradegsl, ntotalrows, gsl_vector_get(gradegsl_i, i));
		gsl_vector_set (energygsl, ntotalrows, gsl_vector_get(energygsl_i, i));

		ntotalrows ++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free (qualitygsl_i);
	gsl_vector_free (gradegsl_i);
	gsl_vector_free (energygsl_i);

	return(EPOK);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
* calculus_macthedFilter: This function calculates the matched filter by using the average pulse read from the pulse
*                         templates library input FITS file
*
* - Divide the matched filter by the energy(eV)
* - Filter template with equal areas below and above the x-axis (in order to be baseline insensitive when calculating
*   the pulses energy): Call area0
* - Free allocate of GSL vectors
****************************************/
int calculus_matchedFilter(gsl_vector **matchedfiltergsl, long mf_size)
{
	int status = EPOK;
	string message = "";
	
	gsl_vector *div = gsl_vector_alloc (mf_size);

	// Divide the matched filter by the energy(eV)
	//gsl_vector_set_all(div,energy*1.602176462e-19);
	gsl_vector_set_all(div,energy);
	gsl_vector_div(*matchedfiltergsl,div);

	// Filter template with equal areas below and above the x-axis
	// (in order to be baseline insensitive when calculating the pulses energy)
	if (area0 (matchedfiltergsl))
	{
	    message = "Cannot run area0 routine";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	gsl_vector_free(div);

    return(EPOK);
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
* writeLib: This function is used when the ENERGY of the pulses is not in the library FITS file (LIBRARY extension).
*           Then, the corresponding matched filter will be written in another extension of the library FITS file
*           (ANX-LIB). (Used only in calibration mode)
*
* Parameters:
*      - energyK: Energy keyword of the pulses
*      - matchedfiltermatrix: Matched filter calculated and to be written
****************************************/
int writeLib (double energyK, gsl_matrix *matchedfiltermatrix)
{
	int status = EPOK, extver=0;
	string message = "";
	long eventcntLib;
	char *tt[1];
	char *tf[1];
	char *tu[1];

	// Library FITS file is already opened
	// If ANX-LIB extension is found, then it is opened. If it is not found, it is created.
	strcpy(extname,"ANX-LIB");
	if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
	{
		// ANX-LIB extension not found => Create it
		status = 0;
		if (fits_create_tbl(inLibObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
		{
		    message = "Cannot crate table " + string(extname) + " in library file ";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);  
		}		  
		if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) + " in library file ";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		eventcntLib = 0;
	}
	else
	{
		if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) + " in library file ";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(keyname,"EVENTCNT");
		if (fits_read_key(inLibObject,TLONG,keyname, &eventcntLib,comment,&status))
		{
		    message = "Cannot read key " + string(keyname) + " in library file ";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}
	strcpy(keyname,"EVENTCNT");
	keyvalint = eventcntLib+1;
	if (fits_update_key(inLibObject,TLONG,keyname,&keyvalint,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in library file";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);	  
	}
	
	gsl_vector *energyout = gsl_vector_alloc(1);
	gsl_vector *estenergyout = gsl_vector_alloc(1);
	gsl_matrix *igsl_matrixout = gsl_matrix_alloc(1,eventsz);
	gsl_vector_set(energyout, 0, energyK);
	gsl_vector_set(estenergyout, 0, -1.0);
	gsl_matrix_set_all(igsl_matrixout,-1.0);

	//  Creating ENERGY column
	obj.inObject = inLibObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"ANX-LIB");
	obj.iniRow = eventcntLib+1;
	obj.endRow = eventcntLib+1;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	strcpy(obj.nameCol,"ENERGY");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"eV");
	if (writeFitsSimple(obj, energyout))
	{
	    message = "Cannot create " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Creating ESTENERGY column
	strcpy(obj.nameCol,"ESTENERGY");
	strcpy(obj.unit," ");
	if (writeFitsSimple(obj, estenergyout))
	{
	    message = "Cannot create " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Creating PULSE column
	strcpy(obj.nameCol,"PULSE");
	strcpy(obj.unit,"Amps");
	if (writeFitsComplex(obj, igsl_matrixout))
	{
	    message = "Cannot create " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Creating MATCHEDF column
	strcpy(obj.nameCol,"MatchedF");
	strcpy(obj.unit," ");
	if (writeFitsComplex(obj, matchedfilters))
	{
	    message = "Cannot create " + string(obj.nameCol) + " column";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	return (EPOK);
}
/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 13 ************************************************************
* calculus_optimalFilter: This function calculates the optimal filter by using the matched filter and the noise spectra.
*
*                     MatchedFilter(f)
* OptimalFilter(f) = ------------------
*                           N^2(f)
*
* - FFT calculus of the matched filter
* 	- Declare variables
* 	- FFT calculus
* 	- Generation of the f's (positive and negative) (different depending on mf_size being odd or even)
* 	- Magnitude and argument for positive and negative f's
* 	- Free allocate of GSL vectors
* - Order the f's and magnitudes according [-fmax,...,0,...,fmax]
* - Calculus of the optimal filter
* 	- Also order the f's and csd's of N(f) according [-fmax,...,0,...,fmax]
* 	- To divide MatchedFilter(f)/N^2(f) => Interpolate the one what has a frequency step bigger (matched filter or N(f)):
* 		- Af{MatchedFilter(f)} > Af{N(f)} => Interpolation of MatchedFilter(f)
* 		- Af{MatchedFilter(f)} < Af{N(f)} => Interpolation of N(f)
* 	- Calculation: MatchedFilter(f)/N^2(f)
* - Calculus of the normalization factor
* - It is necessary to go backwards in order to have the same number of points than MatchedFilter(f):
* 	- Af{MatchedFilter(f)} > Af{N(f)} => Interpolation of OptimalFilter(f)
* 	- Af{MatchedFilter(f)} < Af{N(f)} => Already same number of points MatchedFilter(f)
* - Disorder the f's and magnitudes according [0,...,fmax,-fmax,...] in order to store the OptimalFilter(f) in the output file
* - Inverse FFT (to get the expression of the optimal filter in the time domain)
*	- Complex OptimalFilter(f) => Taking into account the magnitude (MatchedFilter(f)/N^2(f)) and the phase
*	  (given by MatchedFilter(f))
*	- OptimalFilter(t) (It is necessary to delete the last argfilterFFT->size/2 elements)
*
* Parameters:
* - matchedfiltergsl: Input vector which is the matched filter
* - mf_size: Matched filter size
* - optimalfiltergsl: Input/OUTPUT vector which is the optimal filter expressed in the time domain (optimalfilter(t))
* - f_optimal_filtergsl_aux: Input/OUTPUT vector which contains the disordered optimal filter f's according [0,...,fmax,-fmax,...] (frequency domain)
* - optimal_filterFFTgsl_aux: Input/OUTPUT vector which contains the disordered optimal filter magnitudes according [0,...,fmax,-fmax,...] (frequency domain)
****************************************/
int calculus_optimalFilter(gsl_vector *matchedfiltergsl, long mf_size, gsl_vector **optimal_filtergsl, gsl_vector **f_optimal_filtergsl_aux, gsl_vector **optimal_filterFFTgsl_aux)
{
	// FFT calculus of the filter template (MATCHED FILTER->matchedfiltergsl)
	// Declare variables
	double SelectedTimeDuration = mf_size/samprate;
	int status = EPOK;
	string message = "";

	//Complex FFT values for positive and negative f's
	gsl_vector_complex *filterFFTcomp = gsl_vector_complex_alloc(mf_size);
	gsl_vector *freqfiltergsl = gsl_vector_alloc(mf_size);		// f's (positive and negative)
	gsl_vector *filterFFT = gsl_vector_alloc(mf_size);			// Magnitude for positive and negative f's
	gsl_vector *argfilterFFT = gsl_vector_alloc(mf_size);		// Argument for positive and negative f's
	gsl_vector *freqfiltergsl_aux = gsl_vector_alloc(mf_size);	// Ordered f's according [-fmax,...,0,...,fmax]
	gsl_vector *filterFFT_aux = gsl_vector_alloc(mf_size);		// Ordered magnitude according [-fmax,...,0,...,fmax]
	gsl_vector_view temp;

	// FFT calculus
	if (FFT(matchedfiltergsl,filterFFTcomp,SelectedTimeDuration))
	{
	    message = "Cannot run FFT routine to calculate filter template";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

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
			gsl_vector_set(freqfiltergsl,i,i/SelectedTimeDuration);
		}
		gsl_vector_set(freqfiltergsl,mf_size/2,(mf_size/2)/SelectedTimeDuration);
		for (int i=1; i<(mf_size/2); i++)
		{
			gsl_vector_set(freqfiltergsl,i+mf_size/2,(-1.0)*(i+mf_size/2-i*2)/SelectedTimeDuration);
		}
	}
	else	//Odd
	{
		for (int i=0; i< (mf_size)/2; i++)
		{
			gsl_vector_set(freqfiltergsl,i,i/SelectedTimeDuration);
		}
		gsl_vector_set(freqfiltergsl,mf_size/2,(mf_size/2)/SelectedTimeDuration);
		gsl_vector_set(freqfiltergsl,mf_size/2+1,(-1.0)*(mf_size/2)/SelectedTimeDuration);
		for (int i=2; i<=(mf_size/2); i++)
		{
			gsl_vector_set(freqfiltergsl,i+mf_size/2,(-1.0)*(1+mf_size/2-i)/SelectedTimeDuration);
		}
	}

	// Magnitude and argument for positive and negative f's
	gsl_vector_complex_absRIFCA(filterFFT,filterFFTcomp);		// Magnitude
	gsl_vector_complex_argIFCA(argfilterFFT,filterFFTcomp);	    // Argument

	// Free allocate of GSL vectors
	gsl_vector_complex_free(filterFFTcomp);

	// Order the f's and magnitudes according [-fmax,...,0,...,fmax]
	if (reorderFFT(freqfiltergsl,&freqfiltergsl_aux))
	{
	    message = "Cannot run routine reorderFFT for freqfilter";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (reorderFFT(filterFFT,&filterFFT_aux))
	{
	    message = "Cannot run routine reorderFFT for filterFFT";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Optimal filter (Matched_filter(f)/N^2(f))
	// To divide MatchedFilter(f)/N^2(f) => Also order the f's and csd's of N(f) according [-fmax,...,0,...,fmax]
	gsl_vector *freqgsl_aux = gsl_vector_alloc(freqgsl->size);
	gsl_vector *csdgsl_aux = gsl_vector_alloc(freqgsl->size);
	if (reorderFFT(freqgsl,&freqgsl_aux))
	{
	    message = "Cannot run routine reorderFFT for freq";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (reorderFFT(csdgsl,&csdgsl_aux))
	{
	    message = "Cannot run routine reorderFFT for cds";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// To divide MatchedFilter(f)/N^2(f) => Interpolate the one what has a frequency step bigger (matched filter or N(f))
	gsl_vector *f1;
	gsl_vector *FFT1;
	gsl_vector *f2;
	gsl_vector *FFT2;
	gsl_vector *FFT1_2;
	gsl_vector *FFT2_2;
	long n0start = 0;
	long n0end = 0;
	if ((gsl_vector_get(freqfiltergsl,1)-gsl_vector_get(freqfiltergsl,0)) > (gsl_vector_get(freqgsl,1)-gsl_vector_get(freqgsl,0)))
	// Af{MatchedFilter(f)} > Af{N(f)} => Interpolation of MatchedFilter(f)
	{
		gsl_vector *freqfiltergsl_interp = gsl_vector_alloc(freqgsl->size);
	   	gsl_vector *filterFFT_interp = gsl_vector_alloc(freqgsl->size);

	   	if (interpolate (freqfiltergsl_aux, filterFFT_aux, freqgsl->size, gsl_vector_get(freqgsl,1)-gsl_vector_get(freqgsl,0), &freqfiltergsl_interp, &filterFFT_interp, &n0start, &n0end))
	   	{
		    message = "Cannot run routine interpolate for matched filter interpolation";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		f1 = gsl_vector_alloc(freqgsl->size-n0start-n0end);		// Matched filter f
	   	FFT1 = gsl_vector_alloc(freqgsl->size-n0start-n0end);	// Matched filter magnitude
	   	f2 = gsl_vector_alloc(freqgsl->size-n0start-n0end);		// N(f) f
	   	FFT2 = gsl_vector_alloc(freqgsl->size-n0start-n0end);	// N(f) magnitude
	   	temp = gsl_vector_subvector(freqfiltergsl_interp,n0start,freqgsl->size-n0start-n0end);
	   	gsl_vector_memcpy(f1,&temp.vector);
	   	temp = gsl_vector_subvector(filterFFT_interp,n0start,freqgsl->size-n0start-n0end);
	   	gsl_vector_memcpy(FFT1,&temp.vector);
	   	temp = gsl_vector_subvector(freqgsl_aux,n0start,freqgsl->size-n0start-n0end);
		gsl_vector_memcpy(f2,&temp.vector);
	   	temp = gsl_vector_subvector(csdgsl_aux,n0start,freqgsl->size-n0start-n0end);
	   	gsl_vector_memcpy(FFT2,&temp.vector);

	   	gsl_vector_free(freqfiltergsl_interp);
	   	gsl_vector_free(filterFFT_interp);
	}
	else
	// Af{MatchedFilter(f)} < Af{N(f)} => Interpolation of N(f)
	{
	   	gsl_vector *freqgsl_interp = gsl_vector_alloc(freqfiltergsl->size);
		gsl_vector *csdgsl_interp = gsl_vector_alloc(freqfiltergsl->size);

	   	if (interpolate (freqgsl_aux, csdgsl_aux, freqfiltergsl->size, gsl_vector_get(freqfiltergsl,1)-gsl_vector_get(freqfiltergsl,0), &freqgsl_interp, &csdgsl_interp, &n0start, &n0end))
	   	{
		    message = "Cannot run routine interpolate for N(f) filter interpolation";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

	   	f1 = gsl_vector_alloc(freqfiltergsl->size-n0start-n0end);	// Matched filter f
	   	FFT1 = gsl_vector_alloc(freqfiltergsl->size-n0start-n0end);	// Matched filter magnitude
	   	f2 = gsl_vector_alloc(freqfiltergsl->size-n0start-n0end);	// N(f) f
	   	FFT2 = gsl_vector_alloc(freqfiltergsl->size-n0start-n0end);	// N(f) magnitude
	   	temp = gsl_vector_subvector(freqfiltergsl_aux,n0start,freqfiltergsl->size-n0start-n0end);
	   	gsl_vector_memcpy(f1,&temp.vector);
	   	temp = gsl_vector_subvector(filterFFT_aux,n0start,freqfiltergsl->size-n0start-n0end);
	   	gsl_vector_memcpy(FFT1,&temp.vector);
	   	temp = gsl_vector_subvector(freqgsl_interp,n0start,freqfiltergsl->size-n0start-n0end);
	   	gsl_vector_memcpy(f2,&temp.vector);
	   	temp = gsl_vector_subvector(csdgsl_interp,n0start,freqfiltergsl->size-n0start-n0end);
	   	gsl_vector_memcpy(FFT2,&temp.vector);

	   	gsl_vector_free(freqgsl_interp);
	   	gsl_vector_free(csdgsl_interp);
	}
	gsl_vector_free(freqgsl_aux);
	gsl_vector_free(csdgsl_aux);
	gsl_vector_free(freqfiltergsl_aux);
	gsl_vector_free(filterFFT_aux);
	gsl_vector_free(f2);

	FFT1_2= gsl_vector_alloc(FFT1->size);
	FFT2_2= gsl_vector_alloc(FFT2->size);
	gsl_vector_memcpy(FFT1_2,FFT1);

	// Calculation: MatchedFilter(f)/N^2(f)
	gsl_vector_mul(FFT2,FFT2);	// N^2(f)

	gsl_vector_div(FFT1,FFT2);	// Matched filter FFT / N^2(f) => FF1=Optimal filter(f)

	gsl_vector_mul(FFT1_2,FFT1_2);
	gsl_vector_memcpy(FFT2_2,FFT2);

	// Calculus of the normalization factor
	normalizationFactor = 0;
	for (int k=0; k<FFT1->size; k++)
	{
		normalizationFactor = normalizationFactor + gsl_vector_get(FFT1_2,k)/gsl_vector_get(FFT2_2,k);
	}
	gsl_vector_free(FFT1_2);
	gsl_vector_free(FFT2_2);
	gsl_vector_free(FFT2);

   	// It is necessary to go backwards in order to have the same number of points than MatchedFilter(f)
	gsl_vector *f_optimal_filtergslPREV;	// Ordered Optimal filter(f) f
	gsl_vector *optimal_filterFFTgslPREV;	// Ordered Optimal filter(f) magnitude
	if ((gsl_vector_get(freqfiltergsl,1)-gsl_vector_get(freqfiltergsl,0)) > (gsl_vector_get(freqgsl,1)-gsl_vector_get(freqgsl,0)))
	// Af{MatchedFilter(f)} > Af{N(f)} => Interpolation of OptimalFilter(f)
	{
		f_optimal_filtergslPREV = gsl_vector_alloc(freqfiltergsl->size);	// Ordered Optimal filter(f) f
		optimal_filterFFTgslPREV = gsl_vector_alloc(freqfiltergsl->size);
		if (interpolate (f1, FFT1, freqfiltergsl->size, gsl_vector_get(freqfiltergsl,1)-gsl_vector_get(freqfiltergsl,0), &f_optimal_filtergslPREV, &optimal_filterFFTgslPREV, &n0start, &n0end))
		{
		    message = "Cannot run routine interpolate for optimal filter interpolation";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		gsl_vector *v1_aux = gsl_vector_alloc(f_optimal_filtergslPREV->size);
		gsl_vector_memcpy(v1_aux,f_optimal_filtergslPREV);
		gsl_vector_free(f_optimal_filtergslPREV);
		f_optimal_filtergslPREV = gsl_vector_alloc((optimal_filterFFTgslPREV->size)-n0start-n0end);
		temp = gsl_vector_subvector(v1_aux,n0start,v1_aux->size-n0start-n0end);
		gsl_vector_memcpy(f_optimal_filtergslPREV,&temp.vector);
		gsl_vector_memcpy(v1_aux,optimal_filterFFTgslPREV);
		gsl_vector_free(optimal_filterFFTgslPREV);
		optimal_filterFFTgslPREV = gsl_vector_alloc(f_optimal_filtergslPREV->size);
		temp = gsl_vector_subvector(v1_aux,n0start,v1_aux->size-n0start-n0end);
		gsl_vector_memcpy(optimal_filterFFTgslPREV,&temp.vector);
		gsl_vector_free(v1_aux);
	}
	else
	// Af{MatchedFilter(f)} < Af{N(f)} => Already the same number of points than MatchedFilter(f)
	{
		f_optimal_filtergslPREV = gsl_vector_alloc(f1->size);	// Ordered Optimal filter(f) f
		optimal_filterFFTgslPREV = gsl_vector_alloc(f1->size);
		gsl_vector_memcpy(f_optimal_filtergslPREV,f1);
	   	gsl_vector_memcpy(optimal_filterFFTgslPREV,FFT1);
	}
	gsl_vector_free(f1);
	gsl_vector_free(FFT1);

	// Disorder the f's and magnitudes according [0,...,fmax,-fmax,...] in order to store the OptimalFilter(f) in the output file
	*f_optimal_filtergsl_aux = gsl_vector_alloc(f_optimal_filtergslPREV->size);		// Disordered f's according [0,...,fmax,-fmax,...]
	*optimal_filterFFTgsl_aux = gsl_vector_alloc(f_optimal_filtergslPREV->size);	// Disordered magnitude according [0,...,fmax,-fmax,...]
	if (disorderFFT(f_optimal_filtergslPREV, f_optimal_filtergsl_aux))
	{
	    message = "Cannot run routine disorderFFT for filter";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (disorderFFT(optimal_filterFFTgslPREV, optimal_filterFFTgsl_aux))
	{
	    message = "Cannot run routine disorderFFT for filterFFT";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Inverse FFT (to get the expression of the optimal filter in time domain)
	// Complex OptimalFilter(f) => Taking into account magnitude (MatchedFilter(f)/N^2(f)) and phase (given by MatchedFilter(f))
	gsl_vector_complex *optimal_filterFFTgslcomp = gsl_vector_complex_alloc(f_optimal_filtergslPREV->size);
	// OptimalFilter(t) (It is necessary to delete the last argfilterFFT->size/2 elements)
	*optimal_filtergsl = gsl_vector_alloc(f_optimal_filtergslPREV->size);
	gsl_vector *argfilterFFT_aux = gsl_vector_alloc(f_optimal_filtergslPREV->size);
	for (int i=0;i<f_optimal_filtergslPREV->size;i++)
	{
		for (int j=0;j<argfilterFFT->size;j++)
		{
			if (fabs(gsl_vector_get(*f_optimal_filtergsl_aux,i)-gsl_vector_get(freqfiltergsl,j)) < 1e-5)
			{
				gsl_vector_set(argfilterFFT_aux,i,gsl_vector_get(argfilterFFT,j));
				break;
			}
		}
	}
	for (int i=0;i<f_optimal_filtergslPREV->size;i++)
	{
		gsl_vector_complex_set(optimal_filterFFTgslcomp,i,gsl_complex_polar(gsl_vector_get(*optimal_filterFFTgsl_aux,i),gsl_vector_get(argfilterFFT_aux,i)));
	}
	if (FFTinverse(optimal_filterFFTgslcomp,*optimal_filtergsl,(optimal_filterFFTgslcomp->size)/samprate))
	{
	    message = "Cannot run routine FFTinverse to get optimal filter in time domain";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	gsl_vector_complex_free(optimal_filterFFTgslcomp);
	gsl_vector_free(optimal_filterFFTgslPREV);
	gsl_vector_free(f_optimal_filtergslPREV);

	gsl_vector_free(filterFFT);
	gsl_vector_free(argfilterFFT);
	gsl_vector_free(argfilterFFT_aux);
	gsl_vector_free(freqfiltergsl);

    return(EPOK);
}
/*xxxx end of SECTION 13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 14 ************************************************************
* reorderFFT: This function changes the order of the elements of an input vector organized as
*             [-fmax,...,0,...,fmax] according to [0,...,fmax,-fmax,...]
*
*****************************************************************************/
int reorderFFT(gsl_vector *invector, gsl_vector **outvector)
{
	int status = EPOK;
	string message = "";

	if ((invector->size)%2 == 0)	//EVEN
	{
		for (int i=0; i<(invector->size/2-1);i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i+(invector->size/2)+1));
		}
		gsl_vector_set(*outvector,invector->size/2-1,gsl_vector_get(invector,0));
		for (int i=invector->size/2; i<invector->size;i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i-(invector->size/2)+1));
		}
	}
	else
	{
		for (int i=0; i<(invector->size/2);i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i+(invector->size/2)+1));
		}
		for (int i=invector->size/2; i<invector->size;i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i-(invector->size/2)));
		}
	}

	return EPOK;
}
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 15 ************************************************************
* interpolate: This function interpolates an input vector, creating an output vector with the size and
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
int interpolate (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out, long *numzerosstart, long *numzerosend)
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
	xi = -1.0*(size/2-1)*step;
	for (int i=0; i<size; i++)
	{
	    if (((xi < 0) && (fabs(xi) < fabs(x[0]))) || ((xi >= 0) && (fabs(xi) < fabs(x[N-1]))))
	    {
    		yi = gsl_spline_eval (spline, xi, acc);

    		gsl_vector_set(*x_out,i,xi);
    		gsl_vector_set(*y_out,i,yi);
	    }
	    xi = xi+step;
	}
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

	// Free memory
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);

	return EPOK;
}
/*xxxx end of SECTION 15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 16 ************************************************************
* disorderFFT: This function changes the order of the elements of an input vector organized as
*             [0,...,fmax,-fmax,...] according to [-fmax,...,0,...,fmax]
*
*****************************************************************************/
int disorderFFT(gsl_vector *invector, gsl_vector **outvector)
{
	int status = EPOK;

	if ((invector->size)%2 == 0)	//EVEN
	{
		for (int i=0; i<=invector->size/2;i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i+(invector->size/2)-1));
		}
		for (int i=(invector->size/2)+1; i<invector->size;i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i-(invector->size/2)-1));
		}
	}
	else	//ODD
	{
		for (int i=0; i<=invector->size/2;i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i+(invector->size/2)));
		}
		for (int i=(invector->size/2)+1; i<invector->size;i++)
		{
			gsl_vector_set(*outvector,i,gsl_vector_get(invector,i-(invector->size/2)-1));
		}
	}
	gsl_vector_set(*outvector,0,fabs(gsl_vector_get(*outvector,0)));

	return EPOK;
}
/*xxxx end of SECTION 16 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 17 ************************************************************
* find_matchedfilter: This function uses the pulse height (or energy) in order to choose the proper matched filter
*                     of the pulse templates library
*
* It finds the two closest energies in the pulse templates library and interpolates
*  - If pulse height (energy) is lower than the lowest pulse height (energy) in the pulse model library => The matched
* 	 filter corresponding to the lowest pulse height (energy) in the pulse model library is chosen
*  - If pulse height (energy) is higher than the highest pulse height (energy) in the pulse model library => The matched
*    filter corresponding to the highest pulse height (energy) in the pulse model library is chosen
*
* Parameters:
* - ph: Pulse height (or energy) of the pulse whose matched filter is looking for
* - energiesvalues: Vector where the values of EstEnergy of each pulse model are stored
* - matchedfilters: Matrix where all the matched filters of the pulse model library are stored
* - matchedfilterFound: Found matched filter of the pulse whose pulse height is 'ph'
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
/*xxxx end of SECTION 17 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
