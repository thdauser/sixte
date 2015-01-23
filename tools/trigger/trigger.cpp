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
*                      TRIGGER
*
*  File:      trigger.cpp
*  Version:   16.0.0
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
*  version 1.0.0: 	29/02/08	First version
*  version 1.0.1:	12/03/08	Adding keyword: "id_dmx"
* 								Deleting variable: "infile".
* 								Changing name of the extension "RAW_ADC_DATA" to "EUR-TRG"
* 								Adding input parameter "deltaV".
*  version 1.0.2:	01/04/08	Adding function "truncate"	   cout << m << endl;
* 								Adding function "getQuality"
* 								Adding function "setQuality"
* 								Adding input parameter "numBitsQuality"
* 	version 1.1.0	22/04/08	Moving functions "getQuality" and "set Quality" to miscellaneous module
* 								included in Utils package.
* 	version 1.1.3	29/04/08	Changing documentation of the module
*  version 1.3.0	26/05/08	Adding input parameter "ql"
* 								Adding keyword in the output FITS file "ql"
* 								Adding output column "sigma"
* 	version 2.0.0	10/09/08	Adding input parametrer "nameLog"
* 								Adding input parametrer "verbosity" and processing it using RIL libraries.
* 								Adding new function "readInputKeywords"
* 								Deleting in/output keyword "channelcount"
* 								Including error codes.
*  version 2.1.0 	29/09/08	Included "CREATE" and "PROCESS" keyword in the output FITS file.
* 								Included error processing of "INPUT_PARAMETER_NON_DEFINE"
*  version 2.2.0	24/09/08	Choose of input column Name
*  								If the module fail. It will be warning in the output.
* 								Included parameter time to measure.
*								Included new input parameter: writePulse
*								Included new input parameter: selectRow
*								Deleted FindPulses and FindMean function and added to Utils module
* 								Included function "writePulses"
*  version 2.3.0	19/11/08	Obtained tau value in rise slope and fall slope for each pulse
* 								Include obtainTau function
* version 3.0.0		14/01/09 	Changed units of W and wb input parameters from bin to seconds.
* 								Included delDuplicated function.status = DALattributePutReal(trgExten, "B_CF", b_cF, unit, comment, status);
* 								Included oldTstart global variable.
* 								Resolved bug in findMean
* 								deleted the input parameter "sizeArrayPulse"
* 								Deleted input parameter "columnNameI"
* 								Changed the output column name of Current from "V" to "I0"
* 								Included new input keywords
* version 3.1.0		22/01/09	New function "bins2Seconds"
* 								Resolved bug of baseline Extension write.
* 								Changed default value of input parameters: "verbosity -> 1" and "dt_before and dt_after" -> 5e-05
* 								Modificated error code of input parameters.
* version 3.2.0		02/02/09	Changed method to obtain taurise and taufall.
* 								Changed name of output keyword "CREATE" to "CREATOR"
* version 3.2.1		04/02/09	Resolved bug in findPulse.
* 								Checked free memory allocate.
* version 4.0.0		20/02/09	Included, modificated and removed input keywords
* version 4.1.0		24/02/09	Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
* 								Included warning: "wb vector size larger than Row size of the input FITS file. Resize wb vector size to Row size"
* 								Solved bug 4: stops on zero value for the keyword TIMEZERO, error -62436.
* 								Included FTYPE keywords in output FITS file.
* version 4.2.0 	12/03/09	Solved correctly right of tstart values.
* 								The task write truncated pulses and set "1" value to quality column
* 								Change error -62803 by warning: The sizePulse value is lower than the size of the pulse.
* version 5.0.0					Change algorith to obtain tau's.
* 								Convertion of decrease pulses to increase pulses (changed sign).
* 								Included new columns in EUR-TRG extension: "Baseline and Sigma" of each pulses.
* 								Used the keywords libraries
* version 6.0.0     14/04/09    The input parameter selectRow is not used
*                               Delete "PILInit(argc, argv);" and "PILClose(PIL_OK);"
* 								Using IVCAL
*                               Documentation updated
*                               Units in EUR-TRG and BASELINE extensions
* version 7.0.0		04/08/09	Changed inDataIterator, procSegment, findMean2 and findPulses2
* version 7.0.1		18/09/09	Error in the previous tags version
* version 8.0.0		05/08/10    A lot of modifications to change the method to find pulses (tstart)
*								Deleted functions "findMean2" & "findPulses2"
* 								New functions added: "lpf_boxcar", "derMTH", "findTstart" & "findDerPoslength" and moved to Utils Library
*								Deleted input parameters: "ColumnNameI", "W", "wb", "sizePulse", "dt_before", "dt_after", "k", "k2" & "n2"
*								New input parameters added: "tauFALL" & "ntaus" also updated .par file
*								New output keywords added: "TAUFALL" & "NTAUS"
*								Output FITS file column "difTstrt" added. Output FITS file columns "TIME", "EndPulse", "Baseline" & "Sigma" deleted
*								Restructuring of "inDataIterator" & "writePulses" functions
*								Delete references to _dmx files
*								Delete BASELINE extension
* 								Documentation updated
* version 8.1.0		28/09/10	Resolved bug 43
* 								Tend=Tstart+sizePulze_b-1 instead Tend=tstart+sizePulse_b
* 								pulses with the same sign of the derivative (+1 or -1) from "Tstart" to "Tstart+sizePulze_b-1" are truncated
* version 9.0.0		08/11/10	Included Warning --> "lpf_boxcar: tauFALL too small => Cut-off frequency too high => Equivalent to not filter.")
* 								Resolved bug 47
* 								Output FITS file column "IO_NotFiltered" added. Output FITS file column "IO" renamed --> "I0_Filtered".
* 								Documentation updated
* version 10.0.0	16/11/10	Modifications to assig*
*
* n positive polarity to the pulses
* 					29/11/10	Included Warning --> "lpf_boxcar: tauFALL too high => Cut-off frequency too low"). Program finishes.
*					14/12/10	'Log' argument instead 'OK' in some "writeLog" functions
*					??/??/??	findMean not used (use ASQUID)
*					18/02/11	I0_Filtered deleted in the output FITS file (not used in PULSESHAPE)
*					        	I0_NotFiltered renamed as I0
*					23/02/11    If tendi>tstartnexti => Not tendi=tstartnexti. tendi=tstarti+ntaus*tauFALL unless tendi>end of the event
*					28/02/11    PROCESS also includes TRIGGER version (CREATOR)
* version 11.0.0    25/03/11    Quick Look mode no longer used
*					04/04/11    ENERGY=0 => No pulses => It is not going to look for pulses
*					09/05/11    First derivative and threshold used to look for pulses
*					            (findThreshold, kappaClipping, findTtstart1_1)
*					17/05/11    seconds2Bins instead bins2Seconds
*					18/05/11    medianKappaClipping
*					19/05/11    vectorFIL not included in procRow
*					            obtainTau uses vectorNOTFIL instead vectorFIL
*					24/05/11    findMeanSigma used instead findMean (from Utils)status = DALattributePutReal(trgExten, "B_CF", b_cF, unit, comment, status);
*					            New input parameters: samplesUp and nSgms
 *                              Deleted input parameter n (and n_b)
*					25/05/11	obtainTau uses as a pulse tend, the minimum between tend_(i) and tstart_(i+1)
*								obtainTau uses tstart and tend without taking into account the safetyMargin of the tstart
*                   09/06/11    New keyword: PLSPOLAR
*                               Assign positive polarity to the pulses not only using ASQUID but also PLSPOLAR
*                   10/06/11    obtainTau: Fixed to not past end of vector
*                   17/06/11    First derivative is not normalized any more
*                     /  /      Alternative method for filtering and deriving (Savitsky-Golay)
*                   23/06/11    New pulses models library input FITS file (readLib, inDataIteratorLib)
*                               New function derMTHSimple used from pulseProcess (in order to work with the pulse model)
*                   04/08/11    Different input/output parameters in procRow
*                   			Piled-up primary pulses and secondary pulses looked for (used findTstartPrimary from pulseProcess)
*                   18/10/11    New output keyword "chngplrt"
*                   24/10/11    Different input/output parameters in findTstart (tmaxDER)
*						        Pulse model anticipated pulse and pulse maximum & next sample differ less than 5%
*						        adjusted first derivative is fixed as 0 in tstartDERgsl+tmaxDERgsl+"20"
*							    Secondary pulse larger than 1/100 times the preceding pulse
*							    New parameters to handle with code hardpoints:
*               					stopCriteriaMKC
*               					kappaMKC
*               					limitCriteriaMAX
*						            nsAftrtstart
*               					levelPrvPulse
*               					primaryThresholdCriteria
*                   07/11/11    If trgName does not finish as '.fits' and the file trgName+'.fits' already exists =>
*								Data are appended to trgName file => Not allowed in createTriggerFilestatus = DALattributePutReal(trgExten, "B_CF", b_cF, unit, comment, status);
*					15/02/12    Keyword "chngplrt" was written before changing of value
*					07/03/12    New input parameter "scaleFactor"
*					12/03/12    Two different operation modes:
*					                Calibration (to calculate the calibration factors to convert pseudoenergies into energies at the end of the XRAY chain)
*					                Normal (calibration factor are supposed to be known)
*					26/02/12    Deleted variables related to Savitsky-Golay method (methodUsed, nL, nR)
*					xx/04/12    Threshold to look for primary pulses changed depending on the original pulse being single or secondary
*					            (new variable SorSeorPrgsl)
*					xx/04/12    Process to find pulses divided into functions:
*					               findSEPrPulses
*					14/12/12    Changes in 'procRow' and 'findSePrPulses' to iteratively look for pulses
*                   22/01/13    find_model function renamed as 'find_modelOLD'
*                               New 'find_model' function in order to interpolate between pulse shapes
*                               New 'interpolate_model' function
*                               readLib and inDataIteratorLib functions renamed as 'readLibOLD' and 'inDataIteratorLibOLD'
*                               readLibNew and inDataIteratorLibNew functions renamed as 'readLib' and 'inDataIteratorLib'
*                   30/01/13    getEnergy function renamed as 'getPulseHeight'
*                               Comments related to energy changed to comments related to pulse height
*                   05/02/13    New 'Energy' column in the output FITS file where the estimated energy is written
*                   06/02/13    New output keyword "mode"
*                   12/02/13    Solved some particular cases (related to limitCriteriaMAX): pulse_max_index and pulse_max in findSePrPulses
* version 12.0.0    14/02/13    Improved the modularity to find pulses and moved the functions used to find pulses to Utils
*                   15/02/13    Deleted in/out parameter tend in findPulses
*                   02/04/13    Changed some errors codes in readLib (new codes: 62503-62506)
*                   09/04/13    ENERGY column in library FITS file is renamed as ESTENERGY
*                   26/04/13    tAftrtstart = 0
*                   02/05/13    'Atstart' to calculate the tstart precisely (pulse broadening due to the filtering)
*                               It is supposed that PRETRIGGS of the input files to build the library is 1000
* version 13.0.0	14/05/13	PRETRIGS column is also read from the library => Delete previous hardpoint (supposing PRETRIGS=1000)
*                               (in the calibration mode the tstarts are not going to be corrected, it does not matter, so the library is not read)
*                   15/05/13    'b_cF' and 'c_cF' input parameters are written in the output FITS file
*                   17/05/13	In 'readLib, 'colsLib [4]' instead of 'colsLib [5]'
*                   24/05/13    'Atstart = gsl_vector_alloc(nummodels);' only if mode=1 (normal operation)
*                   21/06/13	New method to get the tstart as precisely as possible (iteratively filter with different scaleFactor)
*                   01/10/13	Variables related to PRETRIGS have been deleted or modified because PRETRIGS info is
*                               already not used in XRAYCHAIN
*                   09/10/13    'gsl_vector_set(tstartgsl,i,+gsl_vector_get(tstartgsl,i)+tstart0/j);' instead of 'gsl_vector_set(tstartgsl,i,+gsl_vector_get(tstartgsl,i)+tstart0/(tstartjgsl->size));'
*                   06/11/13    tAftrtstart has changed from 'constant' to 'input parameter'
*                   05/12/13    Process to recalculate precisely the tstarts:
*                                - 'if ((fabs(tstartj-tstart0) <= 2) && ((tstart0LIMITED != 1) || (tstartjLIMITED != 1)) && (tstart0-tstartj != 0))'
*                                  instead 'if (fabs(tstartj-tstart0) <= 0.2*safetyMarginTstart)'
*                                - If the iteration runs 10 times (10 different LPFs) and tstart0 and tstartj are not enough close => The tstarts of
*                                  the last 5 iterations are averaged (not the tstarts of all the iterations)
*                                - If a pulse is piled-up in the tail of a previous one => The process to recalculate precisely the
*                                  tstart is not run
*                   18/12/13    New input parameter 'getTaus'
* version 14.0.0    24/03/14	Added a new output extension EUR-TEST to IFCA tests
*                               (the low-pass filtered and its derivative of the first event are stored)
*                   23/05/14    New output keyword SCLFCTR (with the scaleFactor input parameter)
* 		    		03/07/2014  Removed PIL,RIL,Common dependencies
* version 14.0.1    07/07/2014  Assignation of default parameters re-done. 
* version 14.1.0    07/07/2014  Adapted for new parameter in interativePars function
* version 14.2.0    08/07/2014  Solved bug preventing the task from reading the full command line
* version 14.3.0    11/07/2014  Solved bug preventing a segmentation fault
* version 14.3.1    11/07/14    Comments of the 'initModule' function modified
*                   15/07/14    'initModule' modified in order to accept negative int or double (if parameters are read from command line)
*                   31/07/14    Comments referring to PIL and RIL deleted
*                   05/08/14    Free memory revised
* version 15.0.0    05/09/14    ADC2016 instead RAW_ADC_DATA
*                               PXL02016 instead I0 column
*                               Differences when reading or writing some keywords
*                   10/09/14    50 null initial samples are added to the pulse templates read from the library ('models')
*                   18/09/14    RECORDS instead ADC2016
*                               ADC instead PXL02016 column
*                   29/09/14    Improved the calculus of tend (after obtaining the tstart in a precise way)
*                   03/10/14    Three operation modes: Strict calibration -> Non piled-up monochromatic pulses
*                                                      Calibration -> Piled-up monochromatic pulses
*                                                      Production mode -> Piled-up non monochromatic pulses
*                               Strict calibration mode => Pulse templates library is built
*                               New parameters for 'writePulses'
*                               New functions 'createLibrary' and 'writeLibrary'
*                   15/10/14    New column TAILBFR, just in case there is a tail before the first pulse of the row
*                   22/10/14    Modified the calculus of tstart in a precise way
* version 16.0.0      Nov/14    Migrated to CFITSIO (removal of ISDC DAL)
*                     Dec/14	Errors processing changed
*                   			Deleted some unnecessary input parameters
*                   			Deleted some unnecessary output keywords
*                   			Coming back again, two operation mode:
*                   				Calibration -> Monochromatic pulses (piled up is allowed)
*                   			    Production -> Piled up non monochromatic pulses
*                   			New functions 'calculateTemplate', 'createHisto', 'align', 'shiftm' and 'shift_m' in order
*                   			to build the pulse templates library by averaging some of the found pulses. A pulseheight
*                   			histogram is used to identify the piled up pulses which are not going to be used to average
*                   			Deleted the calculus of the tstart precisely
*                   12/01/15    'PROCESS' output keyword renamed as 'PROC0'
*                   13/01/15    Added a new input parameter 'crtlib' (create library or not):
*                   				New 'crtlib' plays the same role as the old 'mode' (mode=0 => crtLib=1, mode=1 => crtLib=0)
*                   				'mode' refers to calibration or production and it is not used by TRIGGER but in the following tasks
*                   				"creationlib run mode" instead of "calibration mode" and "notcreationlib run mode" instead of "production mode"
*                   20/01/15    Bug corrected: EVENTSZ keyword was being written as TLONG using int sizePulse_b value (corrected to TINT)
*                               Clobber parameter added. Some renameing of variables/FITS extensions (EUR-TRG -> TRIGGER, EUR-LIB --> LIBRARY, EUR-TEST->TEST)
*                               Quality column was read as int,TDOUBLE in inDataIteratorOutTrg => short, TSHORT
*                   23/01/15    'inDataIterator': Just in case the last record has been filled in with 0's => Re-allocate 'invector'
*
*********************************************************************************************************************/

/******************************************************************************
DESCRIPTION:

The TRIGGER task has two possible run modes, in order to create or not a pulse templates library (i.e., if it is necessary to look for
all the pulses or not): creationlib or notcreationlib run modes.

In the notcreationlib run mode, the goal of the TRIGGER task is try to find all the pulses in the input FITS file. Each record is
low-pass filtered and derived. Then, the process of looking for pulses is scanning the first derivative of the filtered
record and if there is a number of consecutive samples over a pre-established threshold, a pulse has been found. To look
for secondary pulses (due to the piled-up), the adjusted derivative of the record is used. The adjusted derivative is
built by using the difference between the derivative of the filtered found pulses and the corresponding scaled and shifted
derivative of the filtered pulse template.

In the creationlib run mode, the input FITS file only contains monochromatic pulses and, the goal is to build
the pulse templates library by averaging some of the found pulses in order to get the pulse template. The secondary pulses
are not going to be searched for, so the pulse templates library is not going to be used (to make the adjusted derivative)
but to be build. Once the pulses have been found (and written in the output FITS file), a pulseheights histogram is used
to discard the piled up pulses. Pulses have to be aligned before being added and averaged.

The task provides some information in the output FITS file, among others the starting and the ending time of each found
pulse (the task finds the tstart and calculate the tend by using ntaus and tauFALL). Optionally, the task could estimate
the rise and fall times of the pulses. Moreover, the task determines if the initial or last pulse in a record1 is a
truncated pulse and identifies the saturated pulses.

The user must supply the following input parameters:

- inFile: Name of the input FITS file
- outFile: Name of the output FITS file
- inLibFile: Name of the pulses templates library FITS file
- mode: Operation mode, calibration mode (0) or production mode (1)
- crtLib: Run mode, create the pulse templates library (1) or not (0)
- LrsT: Running sum length in seconds (only in notcreationlib mode)
- LbT: Baseline averaging length in seconds (only in notcreationlib mode)
- tauFALL: Fall time of the pulses in seconds
- scaleFactor: Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length
- samplesUp: Consecutive samples over the threshold to locate a pulse
- nSgms: Number of Sigmas to establish the threshold
- ntaus: Number of tauFALLs after starting the pulse to calculate its tend
- writePulses: Write pulses (yes) or not (no) in the output FITS file
- getTaus: Calculate (1) or not (0) the approximate rise and fall times of the pulses
- namelog: Output log file name
- verbosity: Verbosity level of the output log file

The RECORDS extension of the input FITS file must contain 2 columns [TIME, ADC]. The TIME column is a scalar column with
the starting time of each record. The ADC column is a vectorial column with the values of the output current of the TES.

The results are written into an output FITS file (_trg.fits). In the output FITS file there is one row per found pulse and
it contains 8 or 9 columns (depending on if the pulses are written or not) with pulse data as the starting/end time,
pulseheight and quality among other parameters. If the rise and fall times of the pulses are not calculated, the
corresponding columns in the output FITS file will be filled with 0's.

MAP OF SECTIONS IN THIS FILE:

 - 1.- INCLUDE's
 - 2.- MAIN
 - 3.- initModule
 - 4.- readLib
 - 5.- inDataIteratorLib
 - 6.- createLibrary
 - 7.- createTriggerFile
 - 8.- inDataIterator
 - 9.- procRow
 - 10.- obtainTau
 - 11.- writePulses
 - 12.- calculateTemplate
 - 13.- inDataIteratorOutTrg
 - 14.- createHisto
 - 15.- align
 - 16.- shiftm
 - 17.- shift_m
 - 18.- writeLibrary

********************************************************************************/

/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include <trigger.h>


/***** SECTION 2 ************************************
* MAIN function: This function is the main function of the TRIGGER task
*
* - Read input parameters (call initModule)
* - Open input FITS file
* - Read input keywords and check their values
* - If notcreationlib run mode => Read the pulse templates library input FITS file
* 	- Open pulses templates library file (EUR-LIB extension)
* 	- Read the pulse templates library input FITS file and store the whole library in "library" and "models" (call readLib)
* - elseIf creationlib run mode => Build the pulse templates library FITS file
*   (Secondary pulses are not searched for => Pulse templates library not used)
* 	- Call createLibrary
* - If MONOEN = 0 (there are no pulses):
*  	- Create output FITS file
* - elseIf MONOEN > 0 (there are pulses):
*	- Initialize variables and transform from seconds to samples
*	- Get structure of input FITS file columns
* 	- Create output FITS file (call createTriggerFile)
*	- Create structure to run Iteration: inDataIterator
*		- Read columns (TIME and ADC)
*	- If notcreationlib run mode => Pulse templates: Low-pass filtering and first derivative
*	- Called iteration function: inDataIterator
*   - Write EVENTCNT output keywords (because 'totalpulses' was not ok yet when the rest of the keywords were written)
* - If creationlib run mode => Calculate the pulse template by averaging some found pulses
* 	- Call calculateTemplate
* 	- Call writeLibrary
* - Close FITS files
* - Free memory
* - Finalize the task
****************************************************/
int main(int argc, char **argv)
{
	create = "trigger v.16.0.0";			//Set "CREATOR" keyword of output FITS file
	time_t t_start = time(NULL);
	string message="";
	int status=EPOK, extver=0;
	double cutFreq = 0.;
	int boxLength = 0;
	
	sprintf(temporalFileName,"TRIGGERauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if (temporalFile == NULL)
	{
		message = "Cannot open auxiliary file TRIGGERauxfile.txt";
		writeLog(fileRef,"Error",verbosity,message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	sprintf(temporalFileName2,"ppr");
	strcat(temporalFileName2,".txt");
	temporalFile2 = fopen (temporalFileName2,"w");
	if (temporalFile2 == NULL)
	{
		message = "Cannot open auxiliary file ppr.txt";
		writeLog(fileRef,"Error",verbosity,message);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Read input parameters
	if (initModule(argc, argv))
	{
	    message = "Cannot run initModule routine to get input parameters";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	writeLog(fileRef,"Log", verbosity, "Into Trigger task");
  
	// Open input FITS file
	if (fits_open_file(&inObject, inName,0,&status))
	{
	    message = "Cannot open file " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(extname,"RECORDS");
	if (fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname) + " in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}

	// Read	input keywords and check their values
	strcpy(keyname,"TRIGGSZ");
	if (fits_read_key(inObject,TLONG,keyname, &eventsz,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	if (eventsz <= 0)
	{
		message = "Legal values for TRIGGSZ (RECORDS) are integer numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	strcpy(keyname,"DELTAT");
	if (fits_read_key(inObject,TDOUBLE,keyname, &samprate,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	if (samprate <= 0)
	{
		message = "Legal values for DELTAT (RECORDS) are real numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	samprate = 1/samprate;
	
	if (fits_get_num_rows(inObject,&eventcnt, &status))
	{
	    message = "Cannot get number of rows in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	ivcal=1.0;
	asquid = 1.0;
	strcpy(keyname,"MONOEN");
	if (fits_read_key(inObject,TDOUBLE,keyname, &energy,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	if (energy < 0)
	{
		message = "Legal values for MONOEN (RECORDS) are non negative real numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	energy = energy*1e3;
	plspolar = 1.0;

	// If CREATIONLIB run mode (crtLib = 1) => Secondary pulses are not searched for => Pulse templates library not used => Library has to be created
	// If CREATIONLIB mode (crtLib = 0) => Read the pulse templates library input FITS file
	if (crtLib == 0)
	{
	    // Open pulses templates library file (LIBRARY extension)
	    if (fits_open_file(&inLibObject, inLibName,0,&status))
	    {
	    	message = "Cannot open file " +  string(inLibName);
	    	EP_EXIT_ERROR(message,status);
	    }
	    strcpy(extname,"LIBRARY");
	    if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
	    {
	    	message = "Cannot move to HDU " + string(extname) + " in " + string(inLibName);
	    	EP_EXIT_ERROR(message,status);
	    }

	    // Read the pulse templates library input FITS file and store the whole library in "library" and "models"
	    if (readLib())
	    {
	    	message = "Cannot run routine readLib to read pulses library";
	    	EP_EXIT_ERROR(message,EPFAIL);
	    }
	    ntotalrows = 1;
	}
	else
	{
		if (createLibrary())
		{
			message = "Cannot run routine createLibrary to create pulses library";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
		if (fits_open_file(&inLibObject,inLibName,READWRITE,&status))
		{
			message = "Cannot open file " +  string(inLibName);
		    EP_EXIT_ERROR(message,status);
		}
	}

	if (energy == 0)	// There are no pulses
	{
		// Create output FITS file: Trigger file (*_trg.fits)
		if (createTriggerFile())
		{
		    message = "Cannot run routine createTriggerFile to create output file";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
		if (fits_open_file(&trgObject,trgName,READWRITE,&status))
		{
		    message = "Cannot open file " +  string(trgName);
		    EP_EXIT_ERROR(message,status);
		}

		// Write output keywords (their values have been previously checked)
		evtcnt = 0;
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(keyname,"EVENTCNT");
		if (evtcnt < 0)
		{
			message = "Legal values for EVENTCNT (TRIGGER) are integer numbers greater than or equal to 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if (fits_write_key(trgObject,TLONG,keyname,&evtcnt,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}
	}
	else				// There are pulses
	{
		// Initialize variables and transform from seconds to samples
		sizePulse = ntaus * tauFALL;
		safetyMarginTstart = safetyMarginTstart*samprate;
		sizePulse_b = (int)(sizePulse * samprate);
		Lrs = (int)(LrsT*samprate);
		Lb = (int)(LbT*samprate);

		if (sizePulse_b > eventsz)
		{
			writeLog(fileRef, "Warning", verbosity, "Standard vector size larger than Record size of the input FITS file. Resize sizePulse to Record size");
			sizePulse_b = eventsz;
		}
		if (sizePulse_b == 0)
		{
			message = "Calculated pulse size is 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}

		// Get structure of input FITS file columns
		strcpy(extname,"RECORDS");
		if (fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(inName);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"Time");
		if (fits_get_colnum(inObject,0,straux,&colnum,&status))
		{
		    message = "Cannot get column number for " + string(straux) +" in " + string(inName);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"ADC");
		if (fits_get_colnum(inObject,0,straux,&colnum,&status))
		{
		    message = "Cannot get column number for " + string(straux) +" in " + string(inName);
		    EP_EXIT_ERROR(message,status);
		}

		message = "Open InFits: " +  string(inName) + " " ;
		writeLog(fileRef,"Log",verbosity,message);

		// Create output FITS file: Trigger file (*_trg.fits)
		if (createTriggerFile())
		{
		    message = "Cannot create file " +  string(trgName);
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		if (fits_open_file(&obj.inObject,trgName,1,&status))
		{
		    message = "Cannot open file " +  string(trgName);
		    EP_EXIT_ERROR(message,status);
		}

		extern int inDataIterator(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

		// Create structure to run Iteration
		iteratorCol cols [2]; 				// Structure of Iteration
		int n_cols = 2; 					// Number of columns:  TIME + ADC
		long rows_per_loop = 0; 			// 0: Use default: Optimum number of rows
		long offset=0; 						// 0: Process all the rows

		strcpy(extname,"RECORDS");
		if (fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(inName);
		    EP_EXIT_ERROR(message,status);
		}
		
		// Read TIME column
		strcpy(straux,"TIME");
		status = fits_iter_set_by_name(&cols[0], inObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(inName);
		    EP_EXIT_ERROR(message,status);
		}

		// Read ADC column
		strcpy(straux,"ADC");
		status = fits_iter_set_by_name(&cols[1], inObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(inName);
		    EP_EXIT_ERROR(message,status);
		}

		if (crtLib == 0)	// NOTCREATIONLIB mode
		{
			// Check Boxlength
			cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
			boxLength = (int) ((1/cutFreq) * samprate);
			if (boxLength <= 1)
			{
			      message = "lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
			      writeLog(fileRef,"Warning", verbosity,message);
			}
			
			// Once the pulse templates read from the library are low-pass filtered and derived, they are stored in "models"
			for (int i=0; i<models->size1; i++)
			{
				gsl_matrix_get_row(model,models,i);

				// PULSE TEMPLATE: Low-pass filtering
				status = lpf_boxcar(&model,model->size,tauFALL*scaleFactor,samprate);
				if (status == 1)
				{
				    message = "Cannot run routine lpf_boxcar for low-pass filtering";
				    EP_EXIT_ERROR(message,EPFAIL);
				}
				if (status == 3)
				{
					status = EPOK;
				}
				if (status == 4)
				{
					message = "lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low";
					writeLog(fileRef,"Error", verbosity,message);
					EP_EXIT_ERROR(message,EPFAIL);
				}

				// PULSE TEMPLATE: Derivative after filtering
				modelSGN = gsl_vector_alloc(model->size);
				if (derMTHSimple (&model,&modelSGN,model->size))
				{
				    message = "Cannot run routine derMTHSimple to calculate derivative";
				    EP_EXIT_ERROR(message,EPFAIL);
				}

				gsl_matrix_set_row(models,i,model);
			}
			gsl_vector_free(model);
		}

		// Check Boxlength
		cutFreq = 2 * (1/(2*pi*tauFALL*scaleFactor));
		boxLength = (int) ((1/cutFreq) * samprate);
		
		if (boxLength <= 1)
		{
		      message = "lpf_boxcar: tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.";
		      writeLog(fileRef,"Warning", verbosity,message);
		}

		// Called iteration function
		if (fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIterator,0L,&status))
		{
		    message = "Cannot iterate data using InDataTterator";
		    EP_EXIT_ERROR(message,status);
		}

		// Write output keywords (their values have been previously checked)
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}
		
		ttpls1 = totalpulses-1;
		strcpy(keyname,"EVENTCNT");
		if (ttpls1 < 0)
		{
			message = "Legal values for EVENTCNT (TRIGGER) are integer numbers greater than or equal to 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if(fits_write_key(trgObject,TINT,keyname,&ttpls1,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}

		// EVENTCNT = totalpulses-1 because totalpulses has been initialized to 1 in trigger.h

		if ((energy > 0) && (totalpulses == 0))
		{
			writeLog(fileRef, "Warning", verbosity, "Some pulses have not been found");
		}
	}

	if (crtLib == 1)	// CREATIONLIB run mode => Calculate the pulse template by averaging some found pulses
	{
		gsl_vector *pulsetemplate = gsl_vector_alloc(sizePulse_b);
		double pulseheighttemplate = 0;

		if (calculateTemplate (totalpulses,&pulsetemplate,&pulseheighttemplate))
		{
		    message = "Cannot run routine calculateTemplate in creationlib run mode";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
		
		if (writeLibrary(pulseheighttemplate, &pulsetemplate))
		{
		    message = "Cannot run routine writeLibrary in crationlib run mode";
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		gsl_vector_free(pulsetemplate);
	}

	// Close output FITS files
	if (fits_close_file(inObject,&status))
	{
	    message = "Cannot close file " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	if (fits_close_file(trgObject,&status))
	{
	    message = "Cannot close file " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	if (fits_close_file(inLibObject,&status))
	{
	    message = "Cannot close file " + string(inLibName);
	    EP_EXIT_ERROR(message,status);
	}
	
	if (fclose(temporalFile))
	{
	    message = "Cannot close file " + string(temporalFileName);
	    EP_EXIT_ERROR(message,EPFAIL);
	}
	if (fclose(temporalFile2))
	{
	    message = "Cannot close file " + string(temporalFileName2);
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Free memory
	if (energy != 0)
	{
		//delete [] straux;
		delete [] obj.nameTable;
		delete [] obj.nameCol;
		delete [] obj.unit;
	}

	if (crtLib == 0) // Notcreationlib mode
	{
		//gsl_vector_free(modelSGN);

		gsl_matrix_free(library);
		gsl_matrix_free(models);
	}

	// Finalize the task
	time_t t_end = time(NULL);
	sprintf(straux,"%f",(double) (t_end - t_start));
	message = "Time:" + string(straux);
	writeLog(fileRef,"Log", verbosity,message);
	
	writeLog(fileRef,"OK", verbosity,"Trigger Module OK");
	
	if(fclose(fileRef))
	{
	    message = "Cannot close log file";
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
int initModule(int argc, char **argv)
{
	// Define TRIGGER input parameters and assign values to variables
	// Parameter definition and assignation of default values
	const int npars = 17, npars1 = 18;
	inparam triggerPars[npars];
	int optidx =0, par=0; 
	string message="";
	string task="trigger";

	triggerPars[0].name = "inFile";
    triggerPars[0].description = "Input file name";
    triggerPars[0].defValStr = "a.fits";
    triggerPars[0].type =  "char";
	triggerPars[0].ValStr = triggerPars[0].defValStr;

	triggerPars[1].name ="outFile"; 
    triggerPars[1].description = "Output file name";
    triggerPars[1].defValStr = "a_trg.fits";
    triggerPars[1].type = "char";
	triggerPars[1].ValStr = triggerPars[1].defValStr;

	triggerPars[2].name = "inLibFile"; 
    triggerPars[2].description = "Pulse models library file name";
    triggerPars[2].defValStr = "library.fits";
    triggerPars[2].type =  "char";
	triggerPars[2].ValStr = triggerPars[2].defValStr;

	triggerPars[3].name = "mode"; 
    triggerPars[3].description = "Calibration mode (0) or production mode (1)";
    triggerPars[3].defValInt = 1;
    triggerPars[3].type = "int";
    triggerPars[3].minValInt = 0;
    triggerPars[3].maxValInt = 1;
	triggerPars[3].ValInt = triggerPars[3].defValInt;

	triggerPars[4].name = "crtLib";
	triggerPars[4].description = "Create pulse templates library (1) or not (0)";
	triggerPars[4].defValInt = 0;
	triggerPars[4].type = "int";
	triggerPars[4].minValInt = 0;
	triggerPars[4].maxValInt = 1;
	triggerPars[4].ValInt = triggerPars[4].defValInt;

	triggerPars[5].name = "LrsT";
    triggerPars[5].description = "Running sum length (in the RS filter case) (seconds) (only in notcreationmode mode) [>0]";
    triggerPars[5].defValReal = 30.E-6;
    triggerPars[5].type = "double";
    triggerPars[5].minValReal = 1.E-50;
    triggerPars[5].maxValReal = +1.E+50;
	triggerPars[5].ValReal = triggerPars[5].defValReal;

	triggerPars[6].name = "LbT";
    triggerPars[6].description = "Baseline averaging length (in the RS filter case) (seconds) (only in notcreation mode) [>0]";
    triggerPars[6].defValReal = 1.E-3;
    triggerPars[6].type = "double";
    triggerPars[6].minValReal = +1.E-50;
    triggerPars[6].maxValReal = +1.E+50;
	triggerPars[6].ValReal = triggerPars[6].defValReal;

	triggerPars[7].name = "tauFALL";
    triggerPars[7].description = "Fall time of the pulses (seconds) [>0]";
    triggerPars[7].defValReal = 3.5E-4;
    triggerPars[7].type = "double";
    triggerPars[7].minValReal = +1.E-10;
    triggerPars[7].maxValReal = +1.E+10;
	triggerPars[7].ValReal = triggerPars[7].defValReal;

	triggerPars[8].name = "scaleFactor";
    triggerPars[8].description = "Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length [>0]";
    triggerPars[8].defValReal = 0.1;
    triggerPars[8].type = "double";
    triggerPars[8].minValReal = +1.E-10;
    triggerPars[8].maxValReal = +1.E+10;
	triggerPars[8].ValReal = triggerPars[8].defValReal;

	triggerPars[9].name = "samplesUp";
    triggerPars[9].description = "Consecutive samples over the threshold to locate a pulse (bins) [>0]";
    triggerPars[9].defValInt = 20;
    triggerPars[9].type = "int";
    triggerPars[9].minValInt = 1;
    triggerPars[9].maxValInt = int(1E20);
	triggerPars[9].ValInt = triggerPars[9].defValInt;

	triggerPars[10].name = "nSgms";
    triggerPars[10].description = "Number of Sigmas to establish the threshold [>0]";
    triggerPars[10].defValInt = 3;
    triggerPars[10].type = "int";
    triggerPars[10].minValInt = 1;
    triggerPars[10].maxValInt = 100;
	triggerPars[10].ValInt = triggerPars[10].defValInt;

	triggerPars[11].name = "ntaus";
    triggerPars[11].description = "Number of tauFALLs after starting the pulse to calculate its tend [>0]";
    triggerPars[11].defValInt = 15;
    triggerPars[11].type = "int";
    triggerPars[11].minValInt = 1;
    triggerPars[11].maxValInt = 1000;
	triggerPars[11].ValInt = triggerPars[11].defValInt;

	triggerPars[12].name = "writePulse";
    triggerPars[12].description = "Write pulses in the output FITS file? (1:TRUE, 0:FALSE)";
    triggerPars[12].defValInt = 1;
    triggerPars[12].type = "int";
    triggerPars[12].minValInt = 0;
    triggerPars[12].maxValInt = 1;
	triggerPars[12].ValInt = triggerPars[12].defValInt;

	triggerPars[13].name = "getTaus";
    triggerPars[13].description = "Calculate the approximate rise and fall times of the pulses? (1:TRUE, 0:FALSE)";
    triggerPars[13].defValInt = 0;
    triggerPars[13].type = "int";
    triggerPars[13].minValInt = 0;
    triggerPars[13].maxValInt = 1;
	triggerPars[13].ValInt = triggerPars[13].defValInt;

	triggerPars[14].name = "nameLog";
    triggerPars[14].description = "Output log file name";
    triggerPars[14].defValStr = "trg_log.txt";
    triggerPars[14].type = "char";
	triggerPars[14].ValStr = triggerPars[14].defValStr;

	triggerPars[15].name = "verbosity";
    triggerPars[15].description = "Verbosity level of the output log file (in [0,3])";
    triggerPars[15].defValInt = 3;
    triggerPars[15].type = "int";
    triggerPars[15].minValInt = 0;
    triggerPars[15].maxValInt = 3;
	triggerPars[15].ValInt = triggerPars[15].defValInt;
	
	triggerPars[16].name = "clobber";
    triggerPars[16].description = "Re-write output files if clobber=yes";
    triggerPars[16].defValStr = "no";
    triggerPars[16].type = "char";
	triggerPars[16].ValStr = triggerPars[16].defValStr;
	
	// Define structure for command line options
	static struct option long_options[npars1];
	for (optidx = 0; optidx<npars; optidx++)
	{
		long_options[optidx].name= triggerPars[optidx].name.c_str();
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
	    			if (long_options[optidx].name == triggerPars[i].name.c_str())
	    			{
	    				if (triggerPars[i].type == "char") //save char value for par
	    				{
	    					triggerPars[i].ValStr = optarg;
	    				}
	    				else // check if numeric value
	    				{
	    					if ((!isdigit(optarg[0]) && (optarg[0] != '-')) ||
	    							(!isdigit(optarg[0]) && (optarg[0] == '-') && (!isdigit(optarg[1]))))
	    					{
	    						message = "Invalid value for input argument " + string(triggerPars[i].name);
	    						EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	    					}
	    					if (triggerPars[i].type == "int")
	    					{
	    						triggerPars[i].ValInt = atoi(optarg);
	    					}
	    					else
	    					{
	    						triggerPars[i].ValReal= atof(optarg);
	    					}
	    				}
	    				break;
	    			} // endif
	    		} // endfor
	    		break;
	    	default:
	    		message = "Invalid parameter name ";
	    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}//switch
	}//while
	// If command line is empty: ask for params interactively
	if (commandLine == 0)
	{
		if (interactivePars(triggerPars,npars,task))
		{
		    message = "Error reading parameters interactively";
		    EP_PRINT_ERROR(message,EPFAIL); 
		}
	}

	// Save parameter values into meaningful variables
	for (int i=0;i<npars; i++)
	{
		if (triggerPars[i].name == "inFile")
		{
			strcpy(inName, triggerPars[i].ValStr.c_str());
		}
		else if (triggerPars[i].name == "inLibFile")
		{
			strcpy(inLibName, triggerPars[i].ValStr.c_str());
		}
		else if (triggerPars[i].name == "outFile")
		{
			strcpy(trgName, triggerPars[i].ValStr.c_str());
		}
		else if (triggerPars[i].name == "mode")
		{
			mode = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "crtLib")
		{
			crtLib = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "LrsT")
		{
			LrsT = triggerPars[i].ValReal;
		}
		else if (triggerPars[i].name == "LbT")
		{
			LbT = triggerPars[i].ValReal;
		}
		else if (triggerPars[i].name == "tauFALL")
		{
			tauFALL = triggerPars[i].ValReal;
		}
		else if (triggerPars[i].name == "scaleFactor")
		{
			scaleFactor = triggerPars[i].ValReal;
		}
		else if (triggerPars[i].name == "samplesUp")
		{
			samplesUp = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "nSgms")
		{
			nSgms = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "ntaus")
		{
			ntaus = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "writePulse")
		{
			writePulse = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "getTaus")
		{
			getTaus = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "nameLog")
		{
			strcpy(nameLog,triggerPars[i].ValStr.c_str());
		}
		else if (triggerPars[i].name == "verbosity")
		{
			verbosity = triggerPars[i].ValInt;
		}
		else if (triggerPars[i].name == "clobber")
		{
			strcpy(clobberStr, triggerPars[i].ValStr.c_str());
			if(strcmp(clobberStr,"yes")==0){
			  clobber=1;
			}else{
			  clobber=0;
			}
		}

		// Check if parameter value is in allowed range
		if( triggerPars[i].type == "int" &&
				(triggerPars[i].ValInt < triggerPars[i].minValInt ||
						triggerPars[i].ValInt > triggerPars[i].maxValInt))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(triggerPars[i].minValInt) + "," + boost::lexical_cast<std::string>(triggerPars[i].maxValInt) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		else if ( triggerPars[i].type == "double" &&
				(triggerPars[i].ValReal < triggerPars[i].minValReal ||
						triggerPars[i].ValReal > triggerPars[i].maxValReal))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(triggerPars[i].minValReal) + "," + boost::lexical_cast<std::string>(triggerPars[i].maxValReal) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
	} // loop for parameters

	// Open the log file of the task
	fileRef = fopen(nameLog,"w+");	// Remove file if it already exists and open a new file to save log messages
	if (fileRef == NULL)
	{
	    message = "Cannot open Log file " + string(nameLog);
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	return (EPOK);
}
/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
* readLib: This function loads all the pulse templates from the pulses templates library input FITS file
*
******************************************************************************/
int readLib()
{
	string message="";
	int extver=0, status = EPOK;

	strcpy(extname,"LIBRARY");
	if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_get_num_rows(inLibObject,&nummodels, &status))
	{
	    message = "Cannot get number of rows in HDU " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Get structure of pulse templates library input FITS file columns
	strcpy(straux,"ENERGY");
	if (fits_get_colnum(inLibObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get number for column " + string(straux);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(straux,"ESTENERGY");
	if (fits_get_colnum(inLibObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get number for column " + string(straux);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(straux,"PULSE");
	if (fits_get_colnum(inLibObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get number for column " + string(straux);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	extern int inDataIteratorLib(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

	// Create structure to run Iteration
	iteratorCol colsLib [3]; 			// Structure of Iteration
	int n_cols = 3; 					// Number of columns:  Energy + EstEnergy + PULSE
	long rows_per_loop = nummodels; 	// 0: Use default: Optimum number of rows
	long offset=0; 						// 0: Process all the rows

	// Read Energy Column
	strcpy(straux,"Energy");
	status = fits_iter_set_by_name(&colsLib[0], inLibObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate by name column " + string(straux);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Read EstEnergy Column
	strcpy(straux,"EstEnergy");
	status = fits_iter_set_by_name(&colsLib[1], inLibObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate by name column " + string(straux);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Read PULSE
	strcpy(straux,"PULSE");
	status = fits_iter_set_by_name(&colsLib[2], inLibObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate by name column " + string(straux);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Called iteration function
	if (fits_iterate_data(n_cols, colsLib, offset, rows_per_loop, inDataIteratorLib,0L,&status))
	{
	    message = "Cannot iterate data in columns";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	return (EPOK);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* inDataIteratorLib function: This function takes the optimum number of rows to read the pulse templates library input FITS file
*                             and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each row of the pulse templates library
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorLib(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
  	int status = EPOK;
	int extver=0;
	string message = "";

	// Declare variables
	double *energy, *energyin;			// Vector of ENERGY column
	double *estenergy, *estenergyin;	// Vector of ESTENERGY column
	double *pulsemodel, *pulsemodelin;  // Vector of PULSE column
	long eventsz;

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
	if (eventsz <= 0)
	{
		message = "Legal values for EVENTSZ (LIBRARY) are integer numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Allocate input GSL vectors
	gsl_vector *energygsl = gsl_vector_alloc(nummodels);
	gsl_vector *estenergygsl = gsl_vector_alloc(nummodels);
	gsl_matrix *pulsemodelgsl = gsl_matrix_alloc(nummodels,eventsz);
	library = gsl_matrix_alloc(nummodels,2);		// Energy + EstEnergy
	models = gsl_matrix_alloc(nummodels,eventsz);	// All the pulse templates
	gsl_vector *rowaux = gsl_vector_alloc(eventsz);	// Auxiliary variable
	model = gsl_vector_alloc(eventsz);

	// Read iterator
	energyin = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element
	//  of the output array is the null pixel value!
	energy = &energyin[1];
	if (toGslVector(((void **)&energy), &energygsl, nrows, 0, TDOUBLE))
	{
	    message = "Cannot run routine toGslVector for energy vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	estenergyin = (double *) fits_iter_get_array(&cols[1]);
	estenergy = &estenergyin[1];
	if (toGslVector(((void **)&estenergy), &estenergygsl, nrows, 0, TDOUBLE))
	{
	    message = "Cannot run routine toGslVector for estenergy vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	pulsemodelin = (double *) fits_iter_get_array(&cols[2]);
	pulsemodel = &pulsemodelin[1];
	if (toGslMatrix(((void **)&pulsemodel), &pulsemodelgsl, eventsz, nrows, (int)TDOUBLE,0))
	{
	    message = "Cannot run routine toGslVector for pulsemodel vector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Processing each row of the pulse templates library
	for (int i=0; i< nrows; i++)
	{
		gsl_matrix_set(library,ntotalrows-1,0,gsl_vector_get(energygsl,i));
		gsl_matrix_set(library,ntotalrows-1,1,gsl_vector_get(estenergygsl,i));

		gsl_matrix_get_row(rowaux,pulsemodelgsl,i);
		gsl_matrix_set_row(models,ntotalrows-1,rowaux);

		ntotalrows ++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free(energygsl);
	gsl_vector_free(estenergygsl);
	gsl_matrix_free(pulsemodelgsl);
	gsl_vector_free(rowaux);

	return (EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* createLibrary: This function builds the pulse templates library.
*                If the pulse templates library does not exist, it is created
*                If the pulse templates library exists, new information will be appended to it.
*
* - If the pulse templates library exists => Call readLib
* - If the pulse templates does not exist => It is created
****************************************************************************/
int createLibrary()
{
	int extver=0, status=EPOK;
	string message = "";

	// Create pulse templates library FITS file: If it does not exist yet
	// Create pulse templates library file or open it
	fits_open_file(&inLibObject, inLibName,READWRITE,&status);
	if (status == EPOK)
	{
		append = true;
		strcpy(extname,"LIBRARY");
		if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to HDU " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(keyname,"EVENTCNT");
		if (fits_read_key(inLibObject,TLONG,keyname, &eventcntLib,comment,&status))
		{
		    message = "Cannot read keyword " + string(keyname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		if (eventcntLib <= 0)
		{
			message = "Legal values for read EVENTCNT (LIBRARY) are integer numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		eventcntLib1 = eventcntLib + 1;
		strcpy(keyname,"EVENTCNT");
		if (eventcntLib1 <= eventcntLib)
		{
			message = "Legal value for written EVENTCNT (LIBRARY) is read EVENTCNT (LIBRARY) plus 1";
			writeLog(fileRef, "Error", verbosity, message);
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_update_key(inLibObject,TLONG,keyname, &eventcntLib1,comment,&status))
		{
		    message = "Cannot update keyword " + string(keyname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		if (readLib())
		{
		    message = "Cannot run routine readLib to read pulse templates library";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		ntotalrows = 1;
	}
	else
	{
		append = false;
		status = EPOK;
		if (fits_create_file(&inLibObject, inLibName, &status))
		{
		    message = "Cannot create library " + string(inLibName);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		// Create extension LIBRARY
		strcpy(extname,"LIBRARY");
		if (fits_create_tbl(inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status))
		{
		    message = "Cannot create table " + string(extname);
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		evtcnt = 1;
		strcpy(keyname,"EVENTCNT");
		if (evtcnt <= 0)
		{
		    message = "Legal values for EVENTCNT (LIBRARY) are integer numbers greater than 0";
		    writeLog(fileRef, "Error", verbosity, message);
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_write_key(inLibObject,TLONG,keyname,&evtcnt,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"CREATOR");
		strcpy(keyvalstr,create);
		if (fits_write_key(inLibObject,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in library";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	char str_energy[125];       sprintf(str_energy,"%f",energy);

	string process (string("TRIGGER") 	+ ' ' +
	string(inName) 		+ ' ' + string(inLibName) 	  + ' ' +
	string(str_energy)      + ' ' +
	string("(")				+      (string) create 		  +   	  string(")"));

	strcpy(keyname,"PROC0");
	strcpy(keyvalstr,process.c_str());
	if (fits_write_key_longwarn(inLibObject,&status))
	{
	    message = "Cannot write keyword long warn in library";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_write_key_longstr(inLibObject,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in library";
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	} 

	// Close file if newly created
	if (!append)
	{
		if (fits_close_file(inLibObject, &status))
		{
			message = "Cannot close library file " + string(inLibName);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
	}

	return (EPOK);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* createTriggerFile function: This function creates the structure of the _trg output FITS file
*
* - Create _trg output FITS file (if it does not exist yet)
* - Create extension TRIGGER
* - Create keywords
***************************************************************************/
int createTriggerFile()
{
	string message = "";
	int extver=0, status = EPOK;
  
	// Create output FITS file: If it does not exist yet
	// If trgName does not finish as '.fits' and the file trgName+'.fits' already exists =>
	// => Data are appended to trgName file => Must not be allowed
	// '.fits' => 5 characters
	if (strlen(trgName)<6)
	{
		// trgName has 5 or less characters => Does not contain '.fits' =>Append '.fits' to trgName
		char trgNameaux[255];
		sprintf(trgNameaux,trgName);
		strcat(trgNameaux,".fits");
		strcpy(trgName,trgNameaux);
	}
	else if (strlen(trgName)>=6)
	{
		// Check if trgName has '.fits' and if not, append it
		if (strncmp(strndup(trgName+strlen(trgName)-5, 5),".fits",5) != 0)
		{
			// trgName does not finish as '.fits' => Append '.fits' to trgName
			char trgNameaux[255];
			sprintf(trgNameaux,trgName);
			strcat(trgNameaux,".fits");
			strcpy(trgName,trgNameaux);
		}
	}

	// Create _trg file (if file already exists => check clobber)
	
	if (fileExists(string(trgName)) && clobber==1)
	{
	      if (remove(trgName)){
		  message = "Output trigger file already exists & cannot be deleted ("+string(strerror(errno))+")";
		  EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	      }
	}
	else if(fileExists(string(trgName)) && clobber==0)
	{
	      message = "Output trigger file already exists: must not be overwritten";
	      EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if(!fileExists(string(trgName)))
	{
	      if(fits_create_file(&trgObject, trgName, &status))
	      {
		  message = "Cannot create output trigger file " + string(trgName);
		  EP_PRINT_ERROR(message,status); return(EPFAIL);
	      }
	}
	message = "Create Trigger Fits File: " + string(trgName);
	writeLog(fileRef,"Log", verbosity,message);

	// Create extension TRIGGER
	char *tt[1];
	char *tf[1];
	char *tu[1];
	strcpy(extname,"TRIGGER");
	if (fits_open_file(&trgObject,trgName,READWRITE,&status))
	{
	    message = "Cannot open output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_create_tbl(trgObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))
	{
	    message = "Cannot create table " + string(extname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	// Provisional (in order to write useful info to test)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Create extension EUR-TEST and its columns
	strcpy(extname,"EUR-TEST");
	if (fits_create_tbl(trgObject,BINARY_TBL,0,0,tt,tf,tu,extname,&status))
	{
	    message = "Cannot create table " + string(extname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	// Create keywords
	strcpy(extname,"TRIGGER");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
		message = "Cannot move to HDU " + string(extname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	
	strcpy(keyname,"MODE");
	if ((mode != 0) && (mode != 1))
	{
		message = "Legal values for MODE (TRIGGER) are 0 or 1";
		writeLog(fileRef, "Error", verbosity, message);
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (fits_write_key(trgObject,TINT,keyname,&mode,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"EVENTSZ");
	if (sizePulse_b <= 0)
	{
		message = "Legal values for EVENTSZ (TRIGGER) are integer numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (fits_write_key(trgObject,TINT,keyname,&sizePulse_b,comment,&status))
	{
  	    message = "Cannot write keyword " + string(keyname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"ENERGY");
	if (energy < 0)
	{
		message = "Legal values for ENERGY (TRIGGER) are non negative real numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (fits_write_key(trgObject,TDOUBLE,keyname,&energy,comment,&status))
	{
  	    message = "Cannot write keyword " + string(keyname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"SAMPRATE");
	if (samprate <= 0)
	{
		message = "Legal values for SAMPRATE (TRIGGER) are non negative real numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	if (fits_write_key(trgObject,TDOUBLE,keyname,&samprate,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,create);
	if (fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write keyword " + string(keyname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Set PROCESS keyword
	char str_verb[125];			sprintf(str_verb,"%d",verbosity);
	char str_ntaus[125];		sprintf(str_ntaus,"%d",ntaus);
	char str_tauFall[125];		sprintf(str_tauFall,"%e",tauFALL);
	char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",scaleFactor);
	char str_samplesUp[125];	sprintf(str_samplesUp,"%d",samplesUp);
	char str_nSgms[125];	    sprintf(str_nSgms,"%f",nSgms);
	char str_wp[125];			if (writePulse) sprintf(str_wp,"y"); else  sprintf(str_wp,"n");
	char str_getTaus[125];		if (getTaus) sprintf(str_getTaus,"y"); else  sprintf(str_getTaus,"n");
	char str_mode[125];			sprintf(str_mode,"%d",mode);
	char str_crtLib[125];		sprintf(str_crtLib,"%d",crtLib);
	char str_LrsT[125];			sprintf(str_LrsT,"%e",LrsT);
	char str_LbT[125];			sprintf(str_LbT,"%e",LbT);

	string process (string("TRIGGER") 	+ ' ' +
	string(inName) 			+ ' ' + string(inLibName) 	    + ' ' + string(trgName) 	+ ' ' +
	string(str_mode) 		+ ' ' + string(str_crtLib) 		+ ' ' +
	string(str_LrsT)   	    + ' ' + string(str_LbT)         + ' ' +
	string(str_tauFall)     + ' ' + string(str_scaleFactor) + ' ' +
	string(str_samplesUp)   + ' ' + string(str_nSgms)       + ' ' +
	string(str_ntaus) 		+ ' ' +
	string(str_wp)		    + ' ' + string(str_getTaus) + ' ' +
	string(nameLog)			+ ' ' + string(str_verb)	    + ' ' + string(clobberStr) + ' ' +
	string("(")				+      (string) create 		    +   	  string(")"));

	strcpy(keyname,"PROC0");
	strcpy(keyvalstr,process.c_str());
	if (fits_write_key_longwarn(trgObject,&status))
	{
		message = "Cannot write long warn in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_write_key_longstr(trgObject,keyname,keyvalstr,comment,&status))
	{
  	    message = "Cannot write keyword " + string(keyname) + " in output trigger file " + string(trgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);  
	}

	return EPOK;
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* inDataIterator: This function takes the optimum number of rows to read the input FITS file
*                 and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each record
*   - Just in case the last record has been filled in with 0's => Re-allocate invector
*   - Assign positive polarity to the pulses (by using asquid and plspolar)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*   - Low-pass filtering and first derivative
*   - Process each record (call procRow)
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIterator(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	char val[256];
	char val_aux[256];

	string message = "";
	int status = EPOK;
	int extver=0;

	// Declare variables
	double *time, *timein;	// Vector of TIME column
	double *v, *vin;		// Vector of ADC column
	nPulsesRow = 0;			// Number of pulses in each record initialized to 0

	// Allocate input GSL vectors
	gsl_vector *timegsl = gsl_vector_alloc(nrows); 				 // Input TIME column
	gsl_matrix *vgsl = gsl_matrix_alloc(nrows, eventsz); 		 // Input ADC column (matrix)
	gsl_vector *invector = gsl_vector_alloc(eventsz); 			 // Each record
	gsl_vector *invectorNOTFILTERED; 							 // Record without having been filtered
	gsl_vector *invectorFILTERED;	 	                         // Filtered (LPF) record
	gsl_vector *invectorDERIVATIVE;                              // Derivative of invectorFILTERED
	gsl_vector *SGN;
	gsl_vector *derSGN;				                             // Sign of the invectorDERIVATIVE

	// Read iterator
	timein = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	time = &timein[1];
	if (toGslVector(((void **)&time), &timegsl, nrows, 0, TDOUBLE))
	{
	    message = " Cannot run routine toGslVector for time array";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	vin = (double *) fits_iter_get_array(&cols[1]);
	v = &vin[1];
	if (toGslMatrix(((void **)&v), &vgsl, eventsz, nrows, (int)TDOUBLE, 0))
	{
	    message = " Cannot run routine toGslMatrix for v array";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	
	// Processing each record
	for (int i=0; i< nrows; i++)
	{
		sprintf(straux,"%d",ntotalrows);
		message = "-------------> Record: " + string(straux);
		sprintf(straux,"%d",eventcnt);		
		message += " of " + string(straux) + " <------------------ ";
		writeLog(fileRef,"Log", verbosity,message);
		sprintf(val,"-------------> Record: %d of %d <------------------ ",ntotalrows,eventcnt);
		strcat(val,"\n");
		fputs(val,temporalFile);

		double time0 = gsl_vector_get(timegsl, i);	// Initial time of the record
		gsl_matrix_get_row(invector, vgsl, i);

		// Just in case the last record has been filled in with 0's => Re-allocate invector
		if ((ntotalrows == eventcnt) && (gsl_vector_ispos(invector) != 1))
		{
			// Know the new dimension of the last record (elements different from 0)
			long eventszLastRecord;
			eventszLastRecord = gsl_vector_min_index(invector);

			// Auxiliary vector
			gsl_vector *vector_aux = gsl_vector_alloc(eventszLastRecord);
			gsl_vector_view temp;
			temp = gsl_vector_subvector(invector,0,eventszLastRecord);
			gsl_vector_memcpy(vector_aux,&temp.vector);

			// Free invector, allocate with the new dimension and free the auxiliary vector
			gsl_vector_free(invector);
			invector = gsl_vector_alloc(eventszLastRecord);
			gsl_vector_memcpy(invector,vector_aux);
			gsl_vector_free(vector_aux);

			invectorNOTFILTERED = gsl_vector_alloc(eventszLastRecord); 	// Record without having been filtered
			invectorFILTERED = gsl_vector_alloc(eventszLastRecord);	 	// Filtered (LPF) record
			invectorDERIVATIVE = gsl_vector_alloc(eventszLastRecord);  	// Derivative of invectorFILTERED
			SGN = gsl_vector_alloc(eventszLastRecord);
			derSGN = gsl_vector_alloc(eventszLastRecord);				// Sign of the invectorDERIVATIVE
		}
		else
		{
			invectorNOTFILTERED = gsl_vector_alloc(eventsz); 	// Record without having been filtered
			invectorFILTERED = gsl_vector_alloc(eventsz);	 	// Filtered (LPF) record
			invectorDERIVATIVE = gsl_vector_alloc(eventsz);  	// Derivative of invectorFILTERED
			SGN = gsl_vector_alloc(eventsz);
			derSGN = gsl_vector_alloc(eventsz);			    	// Sign of the invectorDERIVATIVE
		}

		gsl_vector_scale(invector,ivcal);			// IVCAL to change arbitrary units of voltage to non-arbitrary
		                                            // units of current (Amps)

		// Assign positive polarity to the pulses
		if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))
		{
			gsl_vector_scale(invector,-1);
			chngplrt = 1;
		}
		strcpy(extname,"TRIGGER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to HDU " + string(extname);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}
		strcpy(keyname,"CHNGPLRT");
		if ((chngplrt != 0) && (chngplrt != 1))
		{
			message = "Legal values for CHNGPLRT (TRIGGER) are 0 or 1";
			writeLog(fileRef, "Error", verbosity, message);
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if (fits_update_key(trgObject,TINT,keyname,&chngplrt,comment,&status))
		{
		    message = "Cannot update keyword " + string(keyname) + " in output trigger file ";
		    EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		// Low-pass filtering
		gsl_vector_memcpy(invectorNOTFILTERED,invector);
		status = lpf_boxcar(&invector,invector->size,scaleFactor*tauFALL,samprate);
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
	  	  writeLog(fileRef,"Error", verbosity,message);
	  	  EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		gsl_vector_memcpy(invectorFILTERED,invector);

		// Derivative after filtering
		if (derMTHSimple (&invector,&SGN,invector->size))
		{
		    message = "Cannot run routine derMTHSimple for derivative after filtering";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		gsl_vector_memcpy(invectorDERIVATIVE,invector);
		gsl_vector_memcpy(derSGN,SGN);

		// Provisional (EUR-TEST extension)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if (indice+1 == 1)
		{
			//  Creating NTFltRow Column
			obj.inObject = trgObject;
			obj.nameTable = new char [255];
			strcpy(obj.nameTable,"EUR-TEST");
			obj.iniRow = 1;
			obj.endRow = eventsz;
			obj.iniCol = 0;
			obj.nameCol = new char [255];
			strcpy(obj.nameCol,"NTFltRow");
			obj.type = TDOUBLE;
			obj.unit = new char [255];
			strcpy(obj.unit," ");
			
			if (writeFitsSimple(obj, invectorNOTFILTERED))
			{
				message = "Cannot run routine writeFitsSimple for invectorNOTFILTERED";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			//  Creating FltRow Column
			strcpy(obj.nameCol,"FltRow");
			if (writeFitsSimple(obj, invectorFILTERED))
			{
				message = "Cannot run routine writeFitsSimple for invectorFILTERED";
				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}

			//  Creating DerRow Column
			strcpy(obj.nameCol,"DerRow");
			if (writeFitsSimple(obj, invectorDERIVATIVE))
			{
				message = "Cannot run routine writeFitsSimple for invectorDERIVATIVE";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		// Process each record
		initialtime = time0;
		if (procRow(invectorNOTFILTERED, invectorDERIVATIVE, &nPulsesRow))
		{
		    message = "Cannot run routine procRow for record processing";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		ntotalrows ++;

		indice++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free(timegsl);
	gsl_matrix_free(vgsl);
	gsl_vector_free(invector);
	gsl_vector_free(invectorNOTFILTERED);
	gsl_vector_free(invectorFILTERED);
	gsl_vector_free(invectorDERIVATIVE);
	gsl_vector_free(SGN);
	gsl_vector_free(derSGN);

	return (EPOK);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* procRow function:  This function processes each record of the input FITS file (ADC column)
*
* - Initialize variables
* - Declare and allocate GSL vectors
* - Find pulses of the record (call findPulses)
* - Calculate the tend of each found pulse
* - If 'getTaus' is true => Obtain the approximate rise and fall times of each pulse
* - Write each pulse into output FITS file (call writePulses)
* - Free allocate of GSL vectors
*
* Parameters:
* - vectorNOTFIL: ADC not filtered
* - vectorDER: Derivative of filtered ADC
* - npulses: Number of found pulses in the record
****************************************************************************/
int procRow(gsl_vector *vectorNOTFIL, gsl_vector *vectorDER, int *npulses)
{
	char val[256];
	char val_aux[256];

	int status = EPOK;
	string message = "";

	// Initialize variables
	nPulsesRow = 0;

	// Allocate GSL vectors
	// To look for pulses
	gsl_vector *tstartgsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *tendgsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *qualitygsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *energygsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector_set_zero(qualitygsl);
	gsl_vector_set_zero(energygsl);								// In order to choose the proper pulse model to calculate
	                                                            // the adjusted derivative and to fill in the Energy column
	                                                            // in the output FITS file

	// To apply obtainTau
	gsl_vector *tstartWOutSftMrggsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *tendWOutSftMrggsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *tauRisegsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *tauFallgsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector_set_zero(tauRisegsl);
	gsl_vector_set_zero(tauFallgsl);

	// Find pulses of the record
	if (crtLib == 0)
	{
		if (findPulses (vectorNOTFIL, vectorDER, &tstartgsl, &qualitygsl, &energygsl,
				npulses,
				1,
				tauFALL, scaleFactor, sizePulse_b, samprate,
				samplesUp, nSgms,
				Lb, Lrs,
				library, models,
				safetyMarginTstart,
				stopCriteriaMKC,
				kappaMKC,
				levelPrvPulse,
				temporalFile,
				indice+1))
		{
			message = "Cannot run routine findPulses";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}
	else if (crtLib == 1)
	{
		if (findPulses (vectorNOTFIL, vectorDER, &tstartgsl, &qualitygsl, &energygsl,
				npulses,
				0,
				tauFALL, scaleFactor, sizePulse_b, samprate,
				samplesUp, nSgms,
				Lb, Lrs,
				library, models,
				safetyMarginTstart,
				stopCriteriaMKC,
				kappaMKC,
				levelPrvPulse,
				temporalFile,
				indice+1))
		{
			message = "Cannot run routine findPulses";
			EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	// Provisional!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	sprintf(val,"%d %d",indice+1,nPulsesRow);
	strcat(val,"\n");
	fputs(val,temporalFile2);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// Calculate the tend of each found pulse
	gsl_vector *tendAUXgsl;
	for (int i=0;i<nPulsesRow;i++)
	{
		if (i==0)	tendAUXgsl = gsl_vector_alloc(nPulsesRow);

		gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i)+sizePulse_b-1); 	//tend_i = tstart_i + (ntaus*tauFALL*samprate)-1

		if ((gsl_vector_get(qualitygsl,i) != 1) || (gsl_vector_get(qualitygsl,i) != 3))
		// If it is already truncated at the beginning, it is not taken into account to classify it again as truncated (at the end)
		{
			if (gsl_vector_get(tendgsl,i) >= vectorDER->size)	// Truncated pulses at the end of the row
			{
				gsl_vector_set(tendgsl,i,(vectorDER->size)-1);
				gsl_vector_set (qualitygsl,i,gsl_vector_get(qualitygsl,i)+1);
			}
		}
		sprintf(val,"Energy: ");
		sprintf(val_aux,"%e",gsl_vector_get(energygsl,i));
		strcat(val,val_aux);
		strcat(val,"\n");
		fputs(val, temporalFile);
		sprintf(val,"Quality: ");
		sprintf(val_aux,"%e",gsl_vector_get(qualitygsl,i));
		strcat(val,val_aux);
		strcat(val,"\n");
		fputs(val, temporalFile);

		if (i != nPulsesRow-1)	// Not last pulse in the event
		{
			if (gsl_vector_get(tendgsl,i) >= gsl_vector_get(tstartgsl,i+1))
			{
				gsl_vector_set(tendAUXgsl,i,gsl_vector_get(tstartgsl,i+1)-1);
			}
			else
			{
				gsl_vector_set(tendAUXgsl,i,gsl_vector_get(tendgsl,i));
			}
		}
		else	// Last pulse in the event
		{
			gsl_vector_set(tendAUXgsl,i,gsl_vector_get(tendgsl,i));
		}
	}

	// The safetyMarginTstart is taken into account
	for (int i=0;i<nPulsesRow;i++)
	{
		gsl_vector_set(tstartgsl,i,gsl_vector_get(tstartgsl,i)+safetyMarginTstart);
	}
	for (int i=0;i<nPulsesRow;i++)
	{
		sprintf(val,"Pulsesfound(");
		sprintf(val_aux,"%d",i);
		strcat(val,val_aux);
		sprintf(val_aux,"): ");
		strcat(val,val_aux);
		sprintf(val_aux,"%.12e",gsl_vector_get(tstartgsl,i));
		strcat(val,val_aux);
		strcat(val,"\n");
		fputs(val, temporalFile);

		if (gsl_vector_get(tstartgsl,i)-safetyMarginTstart < 0)
		{
			gsl_vector_set(tstartgsl,i,0.0);
			if (gsl_vector_get(qualitygsl,i) == 0)
			{
				gsl_vector_set(qualitygsl,i,1.0);
			}
		}
		else
		{
			gsl_vector_set(tstartgsl,i,gsl_vector_get(tstartgsl,i)-safetyMarginTstart);
		}

	}
	for (int i=0;i<nPulsesRow;i++)
	{
		gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i)+sizePulse_b-1); 	//tend_i = tstart_i + (ntaus*tauFALL*samprate)-1

		if (gsl_vector_get(tendgsl,i) >= vectorDER->size)	// Truncated pulses at the end of the record
		{
			gsl_vector_set(tendgsl,i,(vectorDER->size)-1);
			gsl_vector_set (qualitygsl,i,gsl_vector_get(qualitygsl,i)+1);
		}
		if ((nPulsesRow !=1) && (i != nPulsesRow-1)) 		// More than one pulse in the record and not the last one
		{
			if (gsl_vector_get(tendgsl,i) > gsl_vector_get(tstartgsl,i+1))
			{
				gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i+1));
			}
		}
	}

	if (getTaus == 1)
	{
		// Obtain the approximate rise and fall times of each pulse
		// The tstarts and tends which are going to be used do not take into account the safety margin before the tstart
		gsl_vector *safetyMarginVector = gsl_vector_alloc(tstartgsl->size);
		gsl_vector_set_all(safetyMarginVector,1.0);
		gsl_vector_scale(safetyMarginVector,safetyMarginTstart);
		gsl_vector_memcpy(tstartWOutSftMrggsl,tstartgsl);
		gsl_vector_memcpy(tendWOutSftMrggsl,tendgsl);
		gsl_vector_add(tstartWOutSftMrggsl,safetyMarginVector);
		gsl_vector_add(tendWOutSftMrggsl,safetyMarginVector);
		if (obtainTau (vectorNOTFIL, tstartWOutSftMrggsl, tendWOutSftMrggsl, nPulsesRow, &tauRisegsl, &tauFallgsl))
		{
		    message = "Cannot run routine obtainTau to calculate pulses slopes";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
	}

	// Write pulses and other info in output FITS file
	if (nPulsesRow != 0)	pulsesgsl = gsl_matrix_alloc(nPulsesRow,sizePulse_b);
	if (writePulses (vectorNOTFIL, vectorDER, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, energygsl, &pulsesgsl))
	{
	    message = "Cannot run routine writePulses to write pulses in output file";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Free allocate of GSL vectors
	gsl_vector_free(tstartgsl);
	gsl_vector_free(tendgsl);
	gsl_vector_free(qualitygsl);
	gsl_vector_free(energygsl);

	gsl_vector_free(tstartWOutSftMrggsl);
	gsl_vector_free(tendWOutSftMrggsl);
	gsl_vector_free(tauRisegsl);
	gsl_vector_free(tauFallgsl);

	if (nPulsesRow != 0)	gsl_matrix_free (pulsesgsl);

	return EPOK;
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


// TO BE REVISED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/***** SECTION 10 ************************************************************
* obtainTau function: This function obtains the approximate value of the rise and fall time of each pulse
*
* The method to obtain each tau value is to transform the input array values.
* The result of the transformation is a straight line, so it is easier to make an analysis and to obtain its tau values.
*
* Steps:
*       - Declare variables
* 		- Initial values
* 		- Obtained pulse and processed it
* 		- Obtained tau rise
* 		- Obtained tau fall
*
* Parameters:
* 		- invector: Array of current values of an event
* 		- tstartgsl: Start indexes of the pulses
* 		- tendgsl: End indexes of the pulses
* 		- nPulses: Number of pulses in the event
* 		- taurisegsl: Tau value of rise slope of pulse (I/O parameter)
* 		- taufallgsl: Tau value of fall slope of pulse (I/O parameter)
 ******************************************************************************/
int obtainTau (gsl_vector *invector, gsl_vector *tstartgsl, gsl_vector *tendgsl, int nPulses, gsl_vector **taurisegsl, gsl_vector **taufallgsl)
{
	int status = EPOK;

	// Declare variables
	int tstart_index, tend_index;			// tstart and tend of each pulse
	gsl_vector *pulse;						// Each pulse
	int max_index;							// Bin where the maximum of the input value is found
											// It is the limit between rise and fall slopes
	gsl_vector_view temp;

	// Obtained pulse
	for (int i=0; i<nPulses; i++){ 		// Processing each pulse
		tstart_index = (int) gsl_vector_get(tstartgsl, i);

		// tend will be the minimum between tend_(i) and tstart_(i+1)
		if (i != nPulses-1)		// No-last pulse of the event
		{
			if (gsl_vector_get(tstartgsl,i+1) < gsl_vector_get(tendgsl,i))
			{
				tend_index = (int) gsl_vector_get(tstartgsl, i+1);
			}
			else
			{
				tend_index = (int) gsl_vector_get(tendgsl, i);
			}
		}
		else					// Last pulse of the event
		{
			tend_index = (int) gsl_vector_get(tendgsl, i);
			if (tend_index >= invector->size)	tend_index = invector->size-1;
		}

		temp = gsl_vector_subvector(invector,tstart_index,tend_index-tstart_index+1);
		pulse = gsl_vector_alloc(tend_index-tstart_index+1);
		gsl_vector_memcpy(pulse,&temp.vector);	//Each pulse

		tend_index = tend_index - tstart_index;
		tstart_index = 0;
		max_index = gsl_vector_max_index (pulse);

		gsl_vector *y = gsl_vector_alloc(pulse->size);
		gsl_vector *t = gsl_vector_alloc(pulse->size);
		double st=0,sy=0,st2=0,syt=0,n=0;

		///////////////////////// Obtained t and y vectors
		for (int j=0; j<pulse->size; j++)
		{
			gsl_vector_set (y,j,log (gsl_vector_get(pulse,j)));
			gsl_vector_set (t,j,j/1e6);
		}

		// Tau Rise
		n = max_index - tstart_index;
		st=0,sy=0,st2=0,syt=0;

		for (int j = tstart_index; j<max_index; j++)
		{
			if (gsl_vector_get (pulse,j) > 0)
			{
				st = st + gsl_vector_get (t,j);
				sy = sy + gsl_vector_get (y,j);
				st2 = st2 +  pow (gsl_vector_get (t,j),2);
				syt = syt + (gsl_vector_get (y,j) * gsl_vector_get (t,j));
			}
			else n--;
		}

		gsl_vector_set (*taurisegsl,i,-((pow (st,2)/n) - st2)  / (syt - (sy*st/n)));

		// Tau Fall
		n = tend_index - max_index;
		//n = tend_index - max_index - 300;
		//n = tend_index - max_index-200;
		//n = floor((tend_index - max_index)*0.1);

		st=0,sy=0,st2=0,syt=0;

		for (int j = max_index; j<tend_index; j++)
		//for (int j = max_index; j<tend_index-ceil(0.9*(tend_index-max_index)); j++)
		//for (int j = max_index+100; j<tend_index-200; j++)
		//for (int j = max_index+100; j<tend_index-100; j++)
		{
			if (gsl_vector_get (pulse,j) > 0)
			{
				st = st + gsl_vector_get (t,j);
				sy = sy + gsl_vector_get (y,j);
				st2 = st2 +  pow (gsl_vector_get (t,j),2);
				syt = syt + (gsl_vector_get (y,j) * gsl_vector_get (t,j));
			}
			else n--;
		}

		gsl_vector_set (*taufallgsl,i,((pow (st,2)/n) - st2)  / (syt - (sy*st/n)));

		/*tauRISE=-1.0;

		if (gsl_isnan(gsl_vector_get(*taufallgsl,i)) != 1)
		{
			status = obtainTauRise (max_index*1e-6, gsl_vector_max(pulse), gsl_vector_get(*taufallgsl,i), &tauRISE, status);
		}
		gsl_vector_set(*taurisegsl,i,tauRISE);*/

		gsl_vector_free(pulse);
		gsl_vector_free(y);
		gsl_vector_free(t);
	}

	return (EPOK);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


/***** SECTION 11 ************************************************************
* writePulses function: It writes each pulse data into the TRIGGER extension of the output FITS file
*
******************************************************************************/
int writePulses(gsl_vector *invectorNOTFIL, gsl_vector *invectorDER, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, gsl_vector *energy, gsl_matrix **pulses)
{
	int status = EPOK;
	char val[256];
	string message = "";

	// Declare variables
	int t0;   					// First value of index of pulse
	gsl_matrix *vgslout2;
	gsl_vector *difTstart; 		// Difference between the tstart of two consecutive pulses
	                       	   	// For the first pulse of the record, difTtsart will be -1
	gsl_vector *tailBeforegsl;
	bool tailBefore = true;
	int min = 5;

	gsl_vector_view temp;

	if (nPulsesRow!=0)
	{
		vgslout2 = gsl_matrix_alloc(nPulsesRow,sizePulse_b);
		difTstart = gsl_vector_alloc(nPulsesRow);
		gsl_vector_set_all(difTstart,-1);
		tailBeforegsl = gsl_vector_alloc(nPulsesRow);
		gsl_vector_set_all(tailBeforegsl,0);

		// Determine if there is the tail of a previous pulse at the beginning of the record
		if (gsl_vector_get(tstart,0) < min)		min = gsl_vector_get(tstart,0)-1;
		if (gsl_vector_get(tstart,0) != 0)
		{
			for (int j=0;j<min;j++)
			{
				if (gsl_vector_get(invectorDER,gsl_vector_get(tstart,0)-j) > 0)
				{
					tailBefore = false;
					break;
				}
			}
		}

		// Converting bins to time
		for (int i=0; i<nPulsesRow; i++)
		{
			t0 = gsl_vector_get (tstart,i);
			gsl_vector_set(tstart,i,initialtime + (gsl_vector_get (tstart,	i) * (1/samprate)));
			gsl_vector_set(tend,i,initialtime + (gsl_vector_get (tend,	i) * (1/samprate)));

			if (i != 0)		gsl_vector_set(difTstart,i,gsl_vector_get(tstart,i)-gsl_vector_get(tstart,i-1));
			if (i == 0)
			{
				if (tailBefore == true)		gsl_vector_set(tailBeforegsl,i,1);
			}

			if (invectorNOTFIL->size - t0 > sizePulse_b)	//The invectorNOTFIL has more bins than sizePulse
			{
				temp = gsl_vector_subvector(invectorNOTFIL,t0+1,sizePulse_b);
				gsl_matrix_set_row(vgslout2, i, &temp.vector);
			}
			else 									// The invectorNOTFIL has less bins than sizePulse (truncated)
			{
				for (int j=0; j<(invectorNOTFIL->size) - t0; j++)
				{
					if (t0 == -1) t0 = 0;
					gsl_matrix_set (vgslout2,i,j,gsl_vector_get(invectorNOTFIL,j+t0));
				}

				for (int j=(invectorNOTFIL->size)-t0; j< sizePulse_b; j++) {gsl_matrix_set (vgslout2,i,j,0.0);}
			}
		}

		// Creating Tstart Column
		obj.inObject = trgObject;
		obj.nameTable = new char [255];
		strcpy(obj.nameTable,"TRIGGER");
		obj.iniRow = totalpulses;
		obj.endRow = totalpulses+nPulsesRow-1;
		obj.iniCol = 0;
		obj.nameCol = new char [255];
		strcpy(obj.nameCol,"Tstart");
		obj.type = TDOUBLE;
		obj.unit = new char [255];
		strcpy(obj.unit,"seconds");
		
		temp = gsl_vector_subvector(tstart,0,nPulsesRow);
		if (writeFitsSimple(obj, &temp.vector))
		{
		    message = "Cannot run routine writeFitsSimple for column Tstart";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating I0 Column
		if (writePulse)		// If the input parameter writePulses = yes then it puts in output FITS file
		{
			strcpy(obj.nameCol,"I0");
			obj.type = TDOUBLE;
			strcpy(obj.unit,"Amps");
			if (writeFitsComplex(obj, vgslout2))
			{
				message = "Cannot run routine writeFitsComplex for column IO";
				EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
			}
		}
		gsl_matrix_memcpy(*pulses,vgslout2);

		// Creating Tend Column
		strcpy(obj.nameCol,"Tend");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"seconds");
		temp = gsl_vector_subvector(tend,0,nPulsesRow);
		if (writeFitsSimple(obj, &temp.vector))
		{
		    message = "Cannot run routine writeFitsSimple for column Tend";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating tauRise Column
		strcpy(obj.nameCol,"tauRise");
		temp = gsl_vector_subvector(taurise,0,nPulsesRow);
		if (writeFitsSimple(obj, &temp.vector))
		{
		    message = "Cannot run routine writeFitsSimple for column tauRise";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating tauFall Column
		strcpy(obj.nameCol,"tauFall");
		temp = gsl_vector_subvector(taufall,0,nPulsesRow);
		if (writeFitsSimple(obj, &temp.vector))
		{
			message = "Cannot run routine writeFitsSimple for column tauFall";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		// Creating difTstrt Column
		strcpy(obj.nameCol,"difTstrt");
		if (writeFitsSimple(obj, difTstart))
		{
		    message = "Cannot run routine writeFitsSimple for column diffTstrt";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating EstEnrgy Column
		strcpy(obj.nameCol,"EstEnrgy");
		strcpy(obj.unit,"a.u.");
		temp = gsl_vector_subvector(energy,0,nPulsesRow);
		if (writeFitsSimple(obj, &temp.vector))
		{
		    message = "Cannot run routine writeFitsSimple for column EstEnrgy";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}
		
		// Creating Quality Column
		strcpy(obj.nameCol,"Quality");
		obj.type = TSHORT;
		strcpy(obj.unit,"bits");
		temp = gsl_vector_subvector(quality,0,nPulsesRow);
		if (writeFitsSimple(obj, &temp.vector))
		{
		    message = "Cannot run routine writeFitsSimple for column Quality";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		// Creating TailBfr Column
		strcpy(obj.nameCol,"TailBfr");
		obj.type = TSHORT;
		strcpy(obj.unit," ");
		if (writeFitsSimple(obj, tailBeforegsl))
		{
		    message = "Cannot run routine writeFitsSimple for column TailBfr";
		    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
		}

		totalpulses = totalpulses + nPulsesRow;
		sprintf(val,"totalpulses: %d",totalpulses-1);
		strcat(val,"\n");
		fputs(val,temporalFile);

		// Free allocate GSL vectors
		gsl_matrix_free (vgslout2);
		gsl_vector_free (difTstart);
		gsl_vector_free (tailBeforegsl);
	}

	return (EPOK);
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 12 ************************************************************
* calculateTemplate function: In creationlib run mode (monochromatic pulses) this function averages some of the found pulses
*                             in order to get the pulse template and build the pulse templates library FITS file.
*                             It works when the pulses have already been found (and written in the output FITS file)
*                             because it needs to do a pulseheights histogram to discard the piled up pulses.
*                             Only the non piled up pulses are going to be read from the _trg output FITS file (in order
*                             to not handle a long matrix of pulses, found pulses x sizePulse_b).
*                             Pulses have to be aligned before being added and averaged.
*
* - Declare (and initialize) variables
* - Get structure of output FITS file columns
* - Create structure to run Iteration: inDataIteratorOutTrg
* - Call iteration function: inDataIteratorOutTrg
* - Create pulseheights histogram
* - Check if the pulse is piled up or not
* 	-If non piled up pulses => Read, align and average them
* - Free allocate of GSL vectors
******************************************************************************/
int calculateTemplate (int totalPulses, gsl_vector **pulseaverage, double *pulseaverageHeight)
{
	char val[256];

	int status = EPOK;
	string message = "";
	int extver=0;

	// Declare (and initialize) variables
	tstartout = gsl_vector_alloc(totalPulses-1);	// Tstart column from the output trg FITS file
	pulseheight = gsl_vector_alloc(totalPulses-1);	// EstEnrgy column from the output trg FITS file
	qualityout = gsl_vector_alloc(totalPulses-1);	// Quality column from the output trg FITS file

	gsl_vector *nonpileup = gsl_vector_alloc(totalPulses-1);	// Piled up pulse => Not taking into account to calculate the template
	long nonpileupPulses = totalPulses-1;						// A priori, all the found pulses are considered as non piled up
	gsl_vector_set_all(nonpileup,1);

	int nBins = floor(sqrt(totalPulses-1));			// Number of bins of the pulseheights histogram
													// Square-root choice (used by Excel and many others)
	gsl_vector *xhisto = gsl_vector_alloc(nBins);	// X-axis of the pulseheights histogram
	gsl_vector *yhisto = gsl_vector_alloc(nBins);	// Y-axis of the pulseheights histogram
	int index_maximumpulseheight;					// Index where the maximum of the pulseheights histogram is
	double maximumpulseheight;						// Maximum of the pulseheights histogram

	bool firstnonpileupPulse = true;
	gsl_matrix *pulse_matrix = gsl_matrix_alloc(1,sizePulse_b);
	gsl_vector *pulse = gsl_vector_alloc(sizePulse_b);
	gsl_vector *pulseshifted = gsl_vector_alloc(sizePulse_b);

	int Shift= 0;
	double tstartnext;

	// Only if one of the found pulses is validated as non piled up pulse by using EstEnrgy
	// (pulseheights histogram) and Tstart, and Quality => The pulse, I0, is going
	// to be read from the trg output FITS file (in order to not handle a long matrix of
	// pulses, found pulses x sizePulse_b)
	strcpy(obj.nameTable,"TRIGGER");
	strcpy(obj.nameCol,"I0");
	obj.inObject=trgObject;
	obj.type=TDOUBLE;

	// Get structure of output FITS file columns
	if (fits_movnam_hdu(obj.inObject, ANY_HDU,obj.nameTable, extver, &status))
	{
	    message = "Cannot move to HDU " + string(obj.nameTable);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(straux,"Tstart");
	if (fits_get_colnum(obj.inObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column " + string(straux);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(straux,"EstEnrgy");
	if (fits_get_colnum(obj.inObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column " + string(straux);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(straux,"Quality");
	if (fits_get_colnum(obj.inObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column " + string(straux);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	extern int inDataIteratorOutTrg(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

	// Create structure to run Iteration
	iteratorCol cols [3]; 					// Structure of Iteration
	int n_cols = 3; 						// Number of columns:  Tstart + EstEnrgy + Quality
	long rows_per_loop = 0; 				// 0: Use default: Optimum number of rows
	long offset=0; 							// 0: Process all the rows

	strcpy(straux,"Tstart");
	status = fits_iter_set_by_name(&cols[0], obj.inObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate by name for column " + string(straux);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(straux,"EstEnrgy");
	status = fits_iter_set_by_name(&cols[1], obj.inObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate by name for column " + string(straux);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}
	strcpy(straux,"Quality");
	status = fits_iter_set_by_name(&cols[2], obj.inObject, straux, TSHORT, InputCol);
	if (status)
	{
	    message = "Cannot iterate by name for column " + string(straux);
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	// Call inDataIteratorOutTrg function
	ntotalrows = 1;	// Because it has already been used in other inDataIterators functions
	if (fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIteratorOutTrg,0L,&status))
	{
	    message = "Cannot iterate data for columns by inDataIteratorOutTrg";
	    EP_PRINT_ERROR(message,status);return(EPFAIL);
	}

	gsl_vector_scale(tstartout,samprate); 	//tstarts not in sec but in samples

	// Create pulseheights histogram
	if (createHisto(pulseheight, nBins, &xhisto, &yhisto))
	{
	    message = "Cannot run createHisto routine";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	index_maximumpulseheight = gsl_vector_max_index(yhisto);
	maximumpulseheight = gsl_vector_get(xhisto,index_maximumpulseheight);

	for (int i=0;i<totalPulses-1;i++)
	{
		if (i == totalPulses-2)		tstartnext = gsl_vector_get(tstartout,i)+2*sizePulse_b;
		else 						tstartnext = gsl_vector_get(tstartout,i+1);

		// Check if the pulse is piled up or not
		if ((gsl_vector_get(pulseheight,i) < maximumpulseheight-0.1*maximumpulseheight) || (gsl_vector_get(pulseheight,i) > maximumpulseheight+0.1*maximumpulseheight)
			|| (tstartnext-gsl_vector_get(tstartout,i) <= sizePulse_b) || (gsl_vector_get(qualityout,i) >= 1))
		{
 			gsl_vector_set(nonpileup,i,0);
			nonpileupPulses --;
		}
		else
		{
			// Non piled up pulses => Read, align and average them
			if (firstnonpileupPulse == true)
			{
				obj.iniRow = i+1;
				obj.endRow = i+1;
				if (readFitsComplex (obj,&pulse_matrix))
				{
				    message = "Cannot run writeFitsComplex for pulse " + boost::lexical_cast<std::string>(i) + " when 1st pulse is not piled-up";
				    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_matrix_get_row(pulse,pulse_matrix,0);
				gsl_vector_memcpy(*pulseaverage,pulse);
				*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
			}
			else
			{
				obj.iniRow = i+1;
				obj.endRow = i+1;
				if (readFitsComplex (obj,&pulse_matrix))
				{
				    message = "Cannot run writeFitsComplex for pulse " + boost::lexical_cast<std::string>(i) + " when 1st pulse is piled-up";
				    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
				}
				gsl_matrix_get_row(pulse,pulse_matrix,0);

				if (align(pulseaverage,&pulse))
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
	*pulseaverageHeight = *pulseaverageHeight/nonpileupPulses;

	// Free allocate of GSL vectors
	gsl_vector_free(tstartout);
	gsl_vector_free(pulseheight);
	gsl_vector_free(qualityout);

	return (EPOK);
}
/*xxxx end of SECTION 12 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 13 ************************************************************
* inDataIteratorOutTrg function: This function takes the optimum number of rows to read the output FITS file
*                                and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each found pulse
* - Free allocate of GSL vectors
******************************************************************************/
int inDataIteratorOutTrg(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	// Declare variables
	int status = EPOK;
	double *tstart, *tstartin;		// Vector of TSTART column
	double *estenrgy, *estenrgyin;	// Vector of ENERGY column
	//static double *quality;		// Vector of QUALITY column
	short *quality, *qualityin;		// Vector of QUALITY column
	string message = "";

	// Allocate input GSL vectors
	gsl_vector *tstartoutgsl = gsl_vector_alloc(totalpulses-1);
	gsl_vector *estenrgygsl = gsl_vector_alloc(totalpulses-1);
	gsl_vector *qualityoutgsl = gsl_vector_alloc(totalpulses-1);

	// Read iterator
	tstartin = (double *) fits_iter_get_array(&cols[0]);
	tstart = &tstartin[1];
	if (toGslVector(((void **)&tstart), &tstartoutgsl, nrows, 0, TDOUBLE))
	{
	    message = "Cannot run toGslVector routine for vector tstart";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	estenrgyin = (double *) fits_iter_get_array(&cols[1]);
	estenrgy = &estenrgyin[1];
	if (toGslVector(((void **)&estenrgy), &estenrgygsl, nrows, 0, TDOUBLE))
	{
	    message = "Cannot run toGslVector routine for vector estenrgy";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}
	//quality = (double *) fits_iter_get_array(&cols[2]);
	qualityin = (short *) fits_iter_get_array(&cols[2]);
	quality = &qualityin[1];
	if (toGslVector(((void **)&quality), &qualityoutgsl, nrows, 0, TSHORT))
	{
	    message = "Cannot run toGslVector routine for vector quality";
	    EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
	}

	// Processing each found pulse
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set(tstartout,ntotalrows-1,gsl_vector_get(tstartoutgsl,i));
		gsl_vector_set(pulseheight,ntotalrows-1,gsl_vector_get(estenrgygsl,i));
		gsl_vector_set(qualityout,ntotalrows-1,gsl_vector_get(qualityoutgsl,i));

		ntotalrows ++;
	}

	// Free allocate of GSL vectors
	gsl_vector_free(tstartoutgsl);
	gsl_vector_free(estenrgygsl);
	gsl_vector_free(qualityoutgsl);

	return (EPOK);
}
/*xxxx end of SECTION 13 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 14 ************************************************************
* createHisto function: This function builds the histogram of the pulseheights
*
* Parameters:
* - invector: Pulseheights vector
* - nbins: Number of bins of the histogram (it is used an optimal value found in the literature)
* - xhistogsl: X-axis of the histogram
* - yhistogsl: Y-axis of the histogram
******************************************************************************/
int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhistogsl, gsl_vector **yhistogsl)
{
	char val[256];

    int status = EPOK;
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
/*xxxx end of SECTION 14 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 15 ************************************************************
* align function: This function aligns the input vectors by using the shift property of the Fourier transform
*
* - Declare variables
* - FFT of vector1
* - FFT of vector2
* - Phases of the FFT_vector1 and FFT_vector2, *size/(2*pi)
* - Shift between the input vectors
* - shiftdouble into shiftint
* - Move forward or delay vector2 depending on positive or negative shift
******************************************************************************/
int align(gsl_vector **vector1, gsl_vector **vector2)
{
	char val[256];

	int status = EPOK;
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
/*xxxx end of SECTION 15 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 16 ************************************************************
* shiftm function: This function returns (in vectorout) the vectorin delayed m samples
*
******************************************************************************/
int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	char val[256];

	int status = EPOK;

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
/*xxxx end of SECTION 16 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 17 ************************************************************
* shift_m function: This function returns (in vectorout) the vectorin moved forward m samples
*
******************************************************************************/
int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	char val[256];

	int status = EPOK;

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
/*xxxx end of SECTION 17 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 18 ************************************************************
* writeLibrary function: This function writes the pulse template into the LIBRARY extension of the output FITS file
*
******************************************************************************/
int writeLibrary(double estenergy, gsl_vector **pulsetemplate)
{
	int status = EPOK,extver=0;
	string message = "";

	// Declare variables
	gsl_vector *energyoutgsl = gsl_vector_alloc(1);
	gsl_vector *estenergyoutgsl = gsl_vector_alloc(1);
	gsl_matrix *pulsetemplates_matrix = gsl_matrix_alloc(1,sizePulse_b);

    if (append == true)
    {
    	gsl_vector *energycolumn = gsl_vector_alloc(eventcntLib+1);
    	gsl_vector *estenergycolumn = gsl_vector_alloc(eventcntLib+1);

    	gsl_vector *modelsrow = gsl_vector_alloc(sizePulse_b);
    	gsl_matrix *modelsaux = gsl_matrix_alloc(eventcntLib+1,sizePulse_b);
    	gsl_matrix *modelsaux1 = gsl_matrix_alloc(eventcntLib+1,sizePulse_b);

    	for (int i=0;i<eventcntLib;i++)
    	{
    		gsl_vector_set(energycolumn,i,gsl_matrix_get(library,i,0));
    		gsl_vector_set(estenergycolumn,i,gsl_matrix_get(library,i,1));
    		gsl_matrix_get_row(modelsrow,models,i);
    		gsl_matrix_set_row(modelsaux,i,modelsrow);
    	}
    	gsl_vector_set(energycolumn,eventcntLib,energy);
    	gsl_vector_set(estenergycolumn,eventcntLib,estenergy);
    	gsl_matrix_set_row(modelsaux,eventcntLib,*pulsetemplate);

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
    	}
    	gsl_vector_memcpy(energycolumn,energycolumnaux);
    	gsl_vector_memcpy(estenergycolumn,estenergycolumnaux);
    	gsl_matrix_memcpy(modelsaux,modelsaux1);

    	gsl_permutation_free(perm);
    	gsl_vector_free(energycolumnaux);
    	gsl_vector_free(estenergycolumnaux);
    	gsl_matrix_free(modelsaux1);

    	obj.inObject = inLibObject;
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

    	    strcpy(obj.nameCol,"ESTENERGY");
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
    	}
    	gsl_vector_free(modelsrow);

    	gsl_vector_free (energycolumn);
    	gsl_vector_free (estenergycolumn);
    	gsl_matrix_free(modelsaux);
    }
    else
    {
    	strcpy(extname,"LIBRARY");
    	if (fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status))
    	{
    		message = "Cannot move to HDU  " + string(extname) + " in library";
    		EP_PRINT_ERROR(message,status);return(EPFAIL);
    	}
    	if (sizePulse_b <= 0)
    	{
    		message = "Legal values for EVENTSZ (TRIGGER) are integer numbers greater than 0";
    		writeLog(fileRef, "Error", verbosity, message);
    		EP_PRINT_ERROR(message,EPFAIL);return(EPFAIL);
    	}
    	if  (fits_write_key(inLibObject,TINT,"EVENTSZ",&sizePulse_b,comment,&status))
    	{
    		message = "Cannot write keyword EVENTSZ in library";
    		EP_PRINT_ERROR(message,status);return(EPFAIL);
    	}
    	gsl_vector_set (energyoutgsl,0,energy);
    	gsl_vector_set (estenergyoutgsl,0,estenergy);

    	gsl_matrix_set_row(pulsetemplates_matrix,0,*pulsetemplate);

    	// Creating ENERGY Column
    	obj.inObject = inLibObject;
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
	
    	// Creating ESTENERGY Column
    	strcpy(obj.nameCol,"ESTENERGY");
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
    }

    // Free allocate of GSL vectors
    gsl_vector_free(energyoutgsl);
    gsl_vector_free(estenergyoutgsl);

    gsl_matrix_free (pulsetemplates_matrix);

    return (EPOK);
}
/*xxxx end of SECTION 18 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
