/*******************************************************************************************************************
* (C) This product has been produced by IFCA with funding from the Spanish
* Ministry of Science and Innovation (MICINN) under project
* ESP2006-13608-C02-01, as an in-kind contribution to the EURECA
* project.  It remains the property of CSIC (Spanish Council for
* Scientific Research) and it cannot be modified, adapted or used for
* purposes other than the original ones except by its owner.
*
*                  INSTITUTO DE FISICA DE CANTABRIA
*
*                      TRIGGERING PULSES
*
*  File:      trigger.cpp
*  Version:   14.0.1
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
* version 10.0.0	16/11/10	Modifications to assign positive polarity to the pulses
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
*                   Nov/14      Migrated to CFITSIO (removal of ISDC DAL)
*********************************************************************************************************************/

/******************************************************************************
 DESCRIPTION:

The goal of the task TRIGGER is to find the Tstart of all the pulses in the input FITS file (I0 column) and
calculate the Tend (by using ntaus and tauFALL), get the practical tauRise and tauFall of each pulse and
write an output FITS file with pulses data. Additionally, the task determines if the initial or last pulse
into an event is a truncated pulse and identifies the saturated pulses.

The user must supply the following items:
- The name of the input FITS file 
- The name of the pulses models library FITS file
- The operation mode
- The value of some variables to find Tstart and Tend and to analyze the pulses.

If the function scans the first derivative of the filtered event and it finds a number of consecutive bins
(window) over a calculated threshold, a SINGLE pulse is found.

To look for secondary pulses the adjusted first derivative of the whole event is used. The adjusted first derivative
is built by using the differences between each pulse first derivative and the scaled and shifted pulse model
first derivative. At the part where the difference between the first derivatives of the pulse and the model
could be higher because the slope is steeper, the corresponding part of the adjusted first derivative is
fixed as 0. Therefore, if a pulse is piled up in the rising part of another pulse they could not be detected as two
different pulses.

Meanwhile, there is to look for primary pulses which are piled up in the rising part of another one.

If calibration mode (non piled-up monochromatic pulses), the task does not build the adjusted first derivative
and consequently it is not going to look for secondary pulses (no problems to select the pulse model from the
pulse model library). At the end of the XRAY chain, b and c are going to be calculated.

The input FITS file must contain 2 columns [TIME, I0]. The first one is a scalar with the time of the start
of data intake. The second one is an array with the values of the output current.

Input parameters related to the operation mode:

mode: To decide among calibration mode (0) or normal mode (1)
b_cF: Linear calibration factor (only in normal mode)
c_cF: Quadratic calibration factor (only in normal mode)

Input parameters to apply the running sum filter (only in normal mode):

LrsT: Running sum length (seconds)
LbT: Baseline averaging length (seconds)

Input parameters for analyzing the pulse:

tauFALL: Fall time of the pulses (seconds)
scaleFactor: Scale factor to apply to the fall time of the pulses in order to calculate the box-car length
             (For the moment, in order to low-pass filter the event, tauFALL/100 has been shown as optimum) => scaleFactor = 0.01
samplesUp: Consecutive points over the threshold to locate a pulse (bins)
nSgms: Number of Sigmas to establish the threshold
tAftrtstart: Time after the alignment between a pulse and its model to be fixed as 0 in the adjusted derivative

ntaus: Number of tauFALLs after starting the pulse to its end (integer)

The results are written into an output FITS file.
In the output FITS file (type _trg.fits) there is one row per pulse and it contains 7 or 8 columns
(depending on if the pulses are written or not) with pulse data as the starting/end time of the data and
the pulse quality among other parameters (the QUALITY column will be fill in by the PULSESHAPE task).
If the rise and fall times of the pulses are not calculated, the corresponding columns in the output FITS
files will be filled with 0's.
The estimated energy of each found pulse is also written into a column (EstEnrgy).

MAP OF SECTIONS IN THIS FILE:
 - 1.- INCLUDE's
 - 2.- MAIN
 - 3.- initModule
 - 4.- seconds2Bins
 - 5.- createTriggerFile
 - 6.- inDataIterator
 - 7.- procRow
 - 8.- obtainTau
 - 9.- writePulses
 - 10. readLib
 - 11. inDataIteratorLib

 UTILITIES CALLED:
 keywords: getKg, getKi, readKeywords, goodValuesKeywords, writeKeywords
           getKeywordsValueString, getKeywordsValueReal, getKeywordsValueInt, setKeywordsValueInt
 inoutUtils: toGslvector, toGslMatrix, writeFitsSimple, writeFitsComplex, writeLog
 pulseProcess: lpf_boxcar, derMTHSimple, findMeanSigma, medianKappaClipping, gsl_vector_Sumsubvector, getB,
               getPulseHeight, RS_filter, find_model, interpolate_model, findTstart, findPulses,
               findSePrPulses, findTstartPrimary

 *******************************************************************************/

/***** SECTION 1 ************************************
 *       INCLUDE's
 ****************************************************/

#include <trigger.h>

/***** SECTION 2 ************************************
 *       MAIN function: This function is the main function of the TRIGGER task
 *
 * 	- initModule: Read input parameters
 *	- If normal operation mode => Read the pulse models library input FITS file
 *	  (If calibration mode => No secondary pulses searched for => No pulse models library used)
 *	- Open input FITS file
 *	- Read and check input keywords
 *  - Error handling depending on the parameters valid range
 *  - If ENERGY = 0 (there are no pulses):
 *  	* Create output FITS file
 *		* Write output keywords (only EVENTCNT because 'totalpulses' was not ok yet when the rest of the keywords have been written)
 *  - If ENERGY > 0 (there are pulses):
 *      * Initialize input parameters: Convert sizePulse from seconds to bins
 * 		* Get structure of input FITS file columns
 * 		* Create output FITS file
 *		* Iterator DAL
 * 			* Create structure to run Iteration DAL: inDataIterator
 *			* Read Columns
 *				Read TIME
 *				Read V0
 *		    * If normal operation mode => Pulse model: Low-pass filtering and first derivative
 *	 		* Called iteration function
 *      * Write output keywords (only EVENTCNT because 'totalpulses' was not ok yet when the rest of the keywords have been written)
 *  - Close input FITS file
 *	- Close output FITS file
 *  - Free memory
 ****************************************************/
int main(int argc, char **argv) {

        create = "trigger v.15.0.0";			//Set "CREATE" keyword of output FITS file
	time_t t_start = time(NULL);
	char str_stat[8];
	char val[256];
	string message="";
	int status=EPOK, extver=0;
	
	sprintf(temporalFileName,"TRIGGERauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if(temporalFile == NULL)
	{
		writeLog(fileRef,"Error",verbosity,"Cannot open auxiliary file TRIGGERauxfile.txt");
		EPexit(1009);
	}
	sprintf(temporalFileName2,"ppr");
	strcat(temporalFileName2,".txt");
	temporalFile2 = fopen (temporalFileName2,"w");
	if(temporalFile2 == NULL)
	{
		writeLog(fileRef,"Error",verbosity,"Cannot open auxiliary file ppr.txt");
		EPexit(1009);
	}

	// Read input parameters
	status = initModule(argc, argv);
	if(status)EPexit(status);

	writeLog(fileRef,"Log", verbosity, "Into Trigger task");
  
	////////////////////////  	Open input FITS file
	if(fits_open_file(&inObject, inName,0,&status)) printerror(status);
	strcpy(extname,"RECORDS");
	if(fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status)) printerror(status);	
	
	///////////////////////	Input keywords
	strcpy(keyname,"TRIGGSZ");
	if(fits_read_key(inObject,TLONG,keyname, &eventsz,comment,&status)) printerror(status);
	strcpy(keyname,"DELTAT");
	if(fits_read_key(inObject,TDOUBLE,keyname, &samprate,comment,&status)) printerror(status);
	samprate = 1/samprate;
	if(fits_get_num_rows(inObject,&eventcnt, &status)) printerror(status);
	ivcal=1.0;
	asquid = 1.0;
	strcpy(keyname,"MONOEN");
	if(fits_read_key(inObject,TDOUBLE,keyname, &energy,comment,&status)) printerror(status);
	energy = energy*1e3;
	plspolar = 1.0;

	// If CALIBRATION mode (mode = 0) => No secondary pulses searched for => No pulse models library used => Library has to be created
	// If PRODUCTION mode (mode = 1) => Library is read
	if (mode == 1)
	{
	    // Open pulses models library file (EUR-LIB extension)
	    if(fits_open_file(&inLibObject, inLibName,0,&status)) printerror(status);
	    strcpy(extname,"EUR-LIB");
	    if(fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status)) printerror(status);

	    // Read the pulse models library input FITS file and store the whole library in "library" and "models"
	    status = readLib();
	    if(status)EPexit(status);
	    ntotalrows = 1;
	}
	else
	{
		status = createLibrary();
		if(status)EPexit(status);
		if(fits_open_file(&inLibObject,inLibName,READWRITE,&status))printerror(status);
	}

	if (energy == 0)	// There are no pulses
	{
		/////////////////////// Create output FITS file: Trigger file (*_trg.fits)
		status = createTriggerFile();
		if(status)EPexit(status);
		if(fits_open_file(&trgObject,trgName,READWRITE,&status))printerror(status);

		/////////////////////// Write output keywords
		evtcnt = 0;
		strcpy(extname,"EUR-TRG");
		if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(trgObject,TLONG,keyname,&evtcnt,comment,&status)) printerror(status);
	}
	else				// There are pulses
	{
		//////////////////////// Initialize input parameters (and transform from seconds to samples)
		sizePulse = ntaus * tauFALL;
		safetyMarginTstart = safetyMarginTstart*samprate;
		nsAftrtstart = tAftrtstart*samprate;
		sizePulse_b = (int)(sizePulse * samprate);
		Lrs = (int)(LrsT*samprate);
		Lb = (int)(LbT*samprate);

		if (sizePulse_b > eventsz){
			writeLog(fileRef, "Warning", verbosity, "Standard vector size larger than Row size of the input FITS file. Resize sizePulse to Row size");
			sizePulse_b = eventsz;
		}
		if (sizePulse_b == 0){
			writeLog(fileRef, "Error", verbosity, "Calculated pulse size is 0");
			EPexit(2021);
		}

		/////////////////////// Get structure of input FITS file columns

		strcpy(straux,"Time");
		strcpy(extname,"RECORDS");
		if(fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status)) printerror(status);
		if(fits_get_colnum(inObject,0,straux,&colnum,&status)) printerror(status);
		strcpy(straux,"ADC");
		if(fits_get_colnum(inObject,0,straux,&colnum,&status)) printerror(status);

		sprintf(str_stat,"%d",status);
		message = "Open InFits: " +  string(inName) + " " + string(str_stat);
		writeLog(fileRef,"Log",verbosity,message);
		
		/////////////////////// Create output FITS file: Trigger file (*_trg.fits)
		status = createTriggerFile();
		if(status)EPexit(status);
		if(fits_open_file(&obj.inObject,trgName,1,&status))printerror(status);

		///////////////////////// Iteration DAL
		extern int inDataIterator(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

			// Create structure to run Iteration DAL
		iteratorCol cols [2]; 					// Structure of Iteration DAL
		int n_cols = 2; 					// Number of columns:  TIME + ADC
		long rows_per_loop = 0; 				// 0: Use default: Optimum number of rows
		long offset=0; 						// 0: Process all the rows

		strcpy(extname,"RECORDS");
		if(fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status)) printerror(status);
		
		// Read Time column
		strcpy(straux,"TIME");
		if(fits_iter_set_by_name(&cols[0], inObject, straux, TDOUBLE, InputCol)) printerror(status);

		// Read ADC column
		strcpy(straux,"ADC");
		if(fits_iter_set_by_name(&cols[1], inObject, straux, TDOUBLE, InputCol))printerror(status);

		if (mode == 1)	// Production mode
		{
			// Once the pulse models read from the library are low-pass filtered and derived, they are stored in "models"
			for (int i=0; i<models->size1; i++)
			{
				gsl_matrix_get_row(model,models,i);

				// PULSE MODEL: Low-pass filtering
				status = lpf_boxcar(&model,model->size,tauFALL*scaleFactor,samprate);
				if (status == 2023)
				{
					writeLog(fileRef,"Warning", verbosity,"lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.");
					status = EPOK;
				}
				if (status == 2024)
				{
					writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
					return(status);
				}
				// PULSE MODEL: Derivative after filtering
				modelSGN = gsl_vector_alloc(model->size);
				if(derMTHSimple (&model,&modelSGN,model->size)) return (status);

				gsl_matrix_set_row(models,i,model);
			}
			gsl_vector_free(model);
		}

			// Called iteration function
		if(fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIterator,0L,&status))printerror(status);

		/////////////////////// Write output keywords
		strcpy(extname,"EUR-TRG");
		if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
		ttpls1 = totalpulses-1;
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(trgObject,TINT,keyname,&ttpls1,comment,&status)) printerror(status);

		// EVENTCNT = totalpulses-1 because totalpulses has been initialized to 1 in trigger.h

		if ((energy > 0) && (totalpulses == 0))
		{
			writeLog(fileRef, "Warning", verbosity, "Some pulses have not been found");
		}
	}

	if (mode == 0)	// Calibration mode
	{
		gsl_vector *pulsetemplate = gsl_vector_alloc(sizePulse_b);
		double pulseheighttemplate = 0;

		status = calculateTemplate (totalpulses,&pulsetemplate,&pulseheighttemplate);
		if(status)EPexit(status);
		
		status = writeLibrary(pulseheighttemplate, &pulsetemplate);
		if(status)EPexit(status);

		gsl_vector_free(pulsetemplate);
	}

	/////////////////////// Close output FITS files
	// Close pulses models library file (EUR-LIB extension)
	if(fits_close_file(inObject,&status) || fits_close_file(trgObject,&status) || 
	  fits_close_file(inLibObject,&status))  printerror(status);
	
	if(fclose(temporalFile)||fclose(temporalFile2))EPexit(1010);

	/////////////////////// Free memory
	if (energy != 0)
	{
		//delete [] straux;
		delete [] obj.nameTable;
		delete [] obj.nameCol;
		delete [] obj.unit;
	}

	if (mode == 1) //Production mode
	{
		//gsl_vector_free(modelSGN);

		gsl_matrix_free(library);
		gsl_matrix_free(models);
	}

	/////////////////////// Finalize the task
	time_t t_end = time(NULL);
	sprintf(straux,"%f",(double) (t_end - t_start));
	message = "Time:" + string(straux);
	writeLog(fileRef,"Log", verbosity,message);
	
	if (status == EPOK)
	{
		writeLog(fileRef,"OK", verbosity,"Trigger Module OK");
	}
	else
	{
		writeLog(fileRef,"Error", verbosity,"Trigger Module FAILS");
		EPexit(1008);
	}

	if(fclose(fileRef))EPexit(1005);

	return status;
}
/*xxxx end of SECTION 2 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 3 ************************************************************
 * InitModule function: This function reads and processes input parameters
 *
 * - It reads all input parameter of the task
 * 		- Name of the FITS files (input and output) included the input FITS file containing the pulse models
 *      - Parameter to decide which operation mode is going to be used (mode)
 *        (and other parameters to put it into practice, b_cF and c_cF)
 *      - Parameter to apply the running sum filter (only in normal mode) (LrsT, LbT)
 * 		- Parameters to analyze the pulses (tauFALL, scaleFactor, samplesUp, nSgms, ntaus,tAftrtstart)
 * 		- Parameters to create output FITS file (writePulses, getTaus)
 *      - General parameters (nameLog, verbosity)
 * - It initializes log file of the task
 ****************************************************************************/
int initModule(int argc, char **argv)
{
	int status = EPOK;
	//
	// Define Trigger Input parameters and assign values to variables
	//
	// Parameter Definition
	
	// SECTION 1/2 TO BE MODIFIED IF PARAMETERS NEED TO BE MODIFIED
        const int npars = 19, npars1 = 20;
	inparam triggerPars[npars];
	int optidx =0, par=0, fst=0, ipar; 
	int fstLib=0;
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
	triggerPars[3].ValInt = triggerPars[3].defValInt;

	triggerPars[4].name = "b_cF"; 
    triggerPars[4].description = "Linear calibration factor (only in normal mode)";
    triggerPars[4].defValReal = 1.;
    triggerPars[4].type = "double";
	triggerPars[4].ValReal = triggerPars[4].defValReal;

	triggerPars[5].name = "c_cF"; 
    triggerPars[5].description = "Quadratic calibration factor (only in normal mode)";
    triggerPars[5].defValReal = 0;
    triggerPars[5].type = "double";
	triggerPars[5].ValReal = triggerPars[5].defValReal;
	
	triggerPars[6].name = "LrsT"; 
    triggerPars[6].description = "Running sum length (in the RS filter case) (seconds) (only in normal mode)";
    triggerPars[6].defValReal = 30.E-6;
    triggerPars[6].type = "double";
	triggerPars[6].ValReal = triggerPars[6].defValReal;

	triggerPars[7].name = "LbT"; 
    triggerPars[7].description = "Baseline averaging length (in the RS filter case) (seconds) (only in normal mode)";
    triggerPars[7].defValReal = 1E-3;
    triggerPars[7].type = "double";
	triggerPars[7].ValReal = triggerPars[7].defValReal;

	triggerPars[8].name = "tauFALL"; 
    triggerPars[8].description = "Fall time of the pulses (seconds)";
    triggerPars[8].defValReal = 3.5E-4;
    triggerPars[8].type = "double";
	triggerPars[8].ValReal = triggerPars[8].defValReal;

	triggerPars[9].name = "scaleFactor"; 
    triggerPars[9].description = "Scale factor to apply to the fall time of the pulses in order to calculate the box-car length";
    triggerPars[9].defValReal = 0.1;
    triggerPars[9].type = "double";
	triggerPars[9].ValReal = triggerPars[9].defValReal;

	triggerPars[10].name = "samplesUp"; 
    triggerPars[10].description = "Consecutive points over the threshold to locate a pulse (bins)";
    triggerPars[10].defValInt = 20;
    triggerPars[10].type = "int";
	triggerPars[10].ValInt = triggerPars[10].defValInt;

	triggerPars[11].name = "nSgms"; 
    triggerPars[11].description = "Number of Sigmas to establish the threshold";
    triggerPars[11].defValInt = 3;
    triggerPars[11].type = "int";
	triggerPars[11].ValInt = triggerPars[11].defValInt;

	triggerPars[12].name = "ntaus"; 
    triggerPars[12].description = "Number of tauFALLs after starting the pulse to ending it";
    triggerPars[12].defValInt = 15;
    triggerPars[12].type = "int";
	triggerPars[12].ValInt = triggerPars[12].defValInt;
	
	triggerPars[13].name = "tAftrtstart"; 
    triggerPars[13].description = "Time after the alignment between a pulse and its model to be fixed as 0 in the adjusted derivative";
    triggerPars[13].defValReal = 0.;
    triggerPars[13].type, "double";
	triggerPars[13].ValReal = triggerPars[13].defValReal;

    triggerPars[14].name = "numBitsQuality";
    triggerPars[14].description = "Number of bits using for Quality";
    triggerPars[14].defValInt = 16;
    triggerPars[14].type = "int";
	triggerPars[14].ValInt = triggerPars[14].defValInt;

	triggerPars[15].name = "writePulse"; 
    triggerPars[15].description = "Write pulses in the output FITS file? (1:TRUE, 0:FALSE)";
    triggerPars[15].defValInt = 1;
    triggerPars[15].type = "int";
	triggerPars[15].ValInt = triggerPars[15].defValInt;

	triggerPars[16].name = "getTaus"; 
    triggerPars[16].description = "Calculate the approximate rise and fall times of the pulses? (1:TRUE, 0:FALSE)";
    triggerPars[16].defValInt = 0;
    triggerPars[16].type = "int";
	triggerPars[16].ValInt = triggerPars[16].defValInt;

	triggerPars[17].name = "nameLog"; 
    triggerPars[17].description = "Output log file Name";
    triggerPars[17].defValStr = "trg_log.txt";
    triggerPars[17].type = "char";
	triggerPars[17].ValStr = triggerPars[17].defValStr;

	triggerPars[18].name = "verbosity"; 
    triggerPars[18].description = "Verbosity Level of the output log file";
    triggerPars[18].defValInt = 3;
    triggerPars[18].type = "int";
	triggerPars[18].ValInt = triggerPars[18].defValInt;
	// END OF SECTION 1/2 TO BE MODIFIED IF PARAMETERS NEED TO BE MODIFIED
	
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
		for(int i=0;i<npars; i++)
		{
			if(long_options[optidx].name == triggerPars[i].name.c_str())
			{
				if(triggerPars[i].type == "char")
				{
				  triggerPars[i].ValStr = optarg;
				}
				else if(triggerPars[i].type == "int")
				{
				  if ((!isdigit(optarg[0]) && (optarg[0] != '-'))) return(1001);
				  if ((!isdigit(optarg[0]) && (optarg[0] == '-')))
				  {
				    if (!isdigit(optarg[1])) return(1001);
				  }
				  triggerPars[i].ValInt = atoi(optarg);
				}
				else
				{
				  if ((!isdigit(optarg[0]) && (optarg[0] != '-'))) return(1001);
				  if ((!isdigit(optarg[0]) && (optarg[0] == '-')))
				  {
				    if (!isdigit(optarg[1]))	return(1001);
				  }
				  triggerPars[i].ValReal= atof(optarg);
				}
				break;
			} // endif
		} // endfor
		break;
	    default:
	    	cout << "Error: invalid parameter name "<< long_options[optidx].name << endl;
	    	EPexit(1002);
	  }//switch
	}//while
	// If command line is empty: ask for params interactively
	if(commandLine == 0)
	{
		status = interactivePars(triggerPars,npars,task);
		if(status)EPexit(status);
	}

	/*
	//Print parameters read (for debugging purposes)
	cout << "Parameters read from command line:" << endl;
	for(int i=0;i<npars; i++)
	{
	    if(triggerPars[i].type == "char")
	    {
			cout << triggerPars[i].description << "=" << triggerPars[i].ValStr << endl;
	    }
	    else if(triggerPars[i].type == "int")
	    {
			cout << triggerPars[i].description << "=" << triggerPars[i].ValInt << endl;
	    }
	    else
	    {
			cout << triggerPars[i].description << "=" << triggerPars[i].ValReal << endl;
	    }
	}
	*/

	// SAVE parameter values into meaningful variables
	
	// SECTION 2/2 TO BE MODIFIED IF PARAMETERS NEED TO BE MODIFIED
	for(int i=0;i<npars; i++){
	  if(triggerPars[i].name == "inFile") {
	    strcpy(inName, triggerPars[i].ValStr.c_str()); 					
	  }else if(triggerPars[i].name == "inLibFile"){
	    strcpy(inLibName, triggerPars[i].ValStr.c_str()); 
	    fst = access(inLibName, R_OK);
	  }else if(triggerPars[i].name == "outFile"){
	    strcpy(trgName, triggerPars[i].ValStr.c_str()); 
	  }else if(triggerPars[i].name == "mode"){
	    mode = triggerPars[i].ValInt;
	    if ((mode != 0) && (mode != 1))			return (2009);
	  }else if(triggerPars[i].name == "b_cF"){
	    b_cF = triggerPars[i].ValReal; 
	    if (b_cF < 0)					return (2010);
	  }else if(triggerPars[i].name == "c_cF"){
	    c_cF = triggerPars[i].ValReal; 
	    if (c_cF < 0)					return (2011);
	  }else if(triggerPars[i].name == "LrsT"){
	    LrsT = triggerPars[i].ValReal; 
	    if (LrsT <= 0)					return (2012);
	  }else if(triggerPars[i].name == "LbT"){
	    LbT = triggerPars[i].ValReal; 
	    if (LbT <= 0)					return (2013);
	  }else if(triggerPars[i].name == "tauFALL"){
	    tauFALL = triggerPars[i].ValReal; 
	    if (tauFALL <= 0)					return (2003);
	  }else if(triggerPars[i].name == "scaleFactor"){
	    scaleFactor = triggerPars[i].ValReal; 
	    if (scaleFactor <= 0)				return (2008);
	  }else if(triggerPars[i].name == "samplesUp"){
	    samplesUp = triggerPars[i].ValInt; 
	    if (samplesUp <= 0)					return (2006);
	  }else if(triggerPars[i].name == "nSgms"){
	    nSgms = triggerPars[i].ValInt; 
	    if (nSgms <= 0)					return (2005);
	  }else if(triggerPars[i].name == "ntaus"){
	    ntaus = triggerPars[i].ValInt; 
	    if (ntaus <= 0)					return (2004);
	  }else if(triggerPars[i].name == "tAftrtstart"){
	    tAftrtstart = triggerPars[i].ValReal; 
	    if (tAftrtstart < 0)				return (2007);
	  }else if(triggerPars[i].name == "numBitsQuality"){
	    numBitsQual = triggerPars[i].ValInt; 
	    if (numBitsQual< 4 || numBitsQual > 32)		return (2002);
	  }else if(triggerPars[i].name == "writePulse"){
	    writePulse = triggerPars[i].ValInt; 
	  }else if(triggerPars[i].name == "getTaus"){
	    getTaus = triggerPars[i].ValInt;
	  }else if(triggerPars[i].name == "nameLog"){
	    strcpy(nameLog,triggerPars[i].ValStr.c_str()); 
	  }else if(triggerPars[i].name == "verbosity"){
	    verbosity = triggerPars[i].ValInt; 
	    if (verbosity<0 || verbosity>3)			return (2001);
	  }
	}
	// END OF SECTION 2/2 TO BE MODIFIED IF PARAMETERS NEED TO BE MODIFIED

	fileRef = fopen(nameLog,"w+");	// Remove file if it already exists and open a new file to save log messages
	if (fileRef == NULL) return(1003);
	
	return (status);
}


/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
 * createTriggerFile function: This function creates the structure of the output FITS file TRIGGER
 *
 * - Create output FITS file: If it does not exist already
 * - Create extension EUR-TRG
 * - Create keywords
 *
 ***************************************************************************/
int createTriggerFile()
{
	string message;
	int extver=0, status = EPOK;
  
	/////////////////////// Create output FITS files: If it does not exist already
	//If trgName does not finish as '.fits' and the file trgName+'.fits' already exists =>
	//=> Data are appended to trgName file => Must not be allowed
	//'.fits' => 5 characters
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
		//Check if trgName has '.fits' and if not, append it
		if (strncmp(strndup(trgName+strlen(trgName)-5, 5),".fits",5) != 0)
		{
			// trgName does not finish as '.fits' => Append '.fits' to trgName
			char trgNameaux[255];
			sprintf(trgNameaux,trgName);
			strcat(trgNameaux,".fits");
			strcpy(trgName,trgNameaux);
		}
	}

	//Create Trigger File
	// if file already exists => error 
	fits_open_file(&trgObject, trgName,0,&status);
	if (status == EPOK){
	  return(1001);
	}else {
	  status = EPOK;
	}

	if(fits_create_file(&trgObject, trgName, &status)) printerror(status);

	sprintf(straux, "%d", status);
	message = "Create Trigger Fits File: " + string(trgName) + string(straux);
	writeLog(fileRef,"Log", verbosity,message);

	/////////////////////// Create extension EUR-TRG
	char *tt[1];
	char *tf[1];
	char *tu[1];
	strcpy(extname,"EUR-TRG");
	if(fits_open_file(&trgObject,trgName,READWRITE,&status))printerror(status);
	if(fits_create_tbl(trgObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))printerror(status);
	
	/////////////////////// Create extension EUR-TEST and its columns
	strcpy(extname,"EUR-TEST");
	if(fits_create_tbl(trgObject,BINARY_TBL,0,0,tt,tf,tu,extname,&status))printerror(status);

	/////////////////////// Create keywords
	strcpy(extname,"EUR-TRG");
	if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
	
	strcpy(keyname,"EVENTSZ");
	if(fits_write_key(trgObject,TLONG,keyname,&sizePulse_b,comment,&status)) printerror(status);
	strcpy(keyname,"ENERGY");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&energy,comment,&status)) printerror(status);
	strcpy(keyname,"SAMPRATE");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&samprate,comment,&status)) printerror(status);
	
	strcpy(keyname,"NBQUAL");
	if(fits_write_key(trgObject,TINT,keyname,&numBitsQual,comment,&status)) printerror(status);
	strcpy(keyname,"TRG_ID");
	if(fits_write_key(trgObject,TINT,keyname,&trg_id,comment,&status)) printerror(status);
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,create);
	if(fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);

	// Set process keyword
	char str_verb[125];		sprintf(str_verb,"%d",verbosity);
	char str_numB[125];		sprintf(str_numB,"%d",numBitsQual);
	char str_ntaus[125];		sprintf(str_ntaus,"%d",ntaus);
	char str_tauFall[125];		sprintf(str_tauFall,"%e",tauFALL);
	char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",scaleFactor);
	char str_samplesUp[125];	sprintf(str_samplesUp,"%d",samplesUp);
	char str_nSgms[125];	    	sprintf(str_nSgms,"%f",nSgms);
	char str_wp[125];		if (writePulse) sprintf(str_wp,"y"); else  sprintf(str_wp,"n");
	char str_getTaus[125];		if (getTaus) sprintf(str_getTaus,"y"); else  sprintf(str_getTaus,"n");
	char str_mode[125];		sprintf(str_mode,"%d",mode);
	char str_b_cF[125];		sprintf(str_b_cF,"%f",b_cF);
	char str_c_cF[125];		sprintf(str_c_cF,"%f",c_cF);
	char str_LrsT[125];		sprintf(str_LrsT,"%e",LrsT);
	char str_LbT[125];		sprintf(str_LbT,"%e",LbT);
	char str_tAftrtstart[125];	sprintf(str_tAftrtstart,"%e",tAftrtstart);

	string process (string("TRIGGER") 	+ ' ' +
	string(inName) 			+ ' ' + string(inLibName) 	    + ' ' + string(trgName) 	+ ' ' +
	string(str_mode) 		+ ' ' + string(str_b_cF)   	    + ' ' + string(str_c_cF)    + ' ' +
	string(str_LrsT)   	    + ' ' + string(str_LbT)         + ' ' +
	string(str_tauFall)     + ' ' + string(str_scaleFactor) + ' ' +
	string(str_samplesUp)   + ' ' + string(str_nSgms)       + ' ' +
	string(str_ntaus) 		+ ' ' + string(str_tAftrtstart)	+ ' ' +
	string(str_numB)		+ ' ' + string(str_wp)		    + ' ' + string(str_getTaus) + ' ' +
	string(nameLog)			+ ' ' + string(str_verb)	    + ' ' +
	string("(")				+      (string) create 		    +   	  string(")"));

	strcpy(keyname,"PROCESS");
	strcpy(keyvalstr,process.c_str());
	if(fits_write_key_longwarn(trgObject,&status))printerror(status);
	if(fits_write_key_longstr(trgObject,keyname,keyvalstr,comment,&status)) printerror(status);
	
	strcpy(keyname,"FTYPE");
	strcpy(keyvalstr,"pulse");
	if(fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);
	
	strcpy(keyname,"NTAUS");
	if(fits_write_key(trgObject,TINT,keyname,&ntaus,comment,&status)) printerror(status);
	strcpy(keyname,"TAUFALL");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&tauFALL,comment,&status)) printerror(status);
	strcpy(keyname,"SCLFCTR");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&scaleFactor,comment,&status)) printerror(status);
	strcpy(keyname,"MODE");
	if(fits_write_key(trgObject,TINT,keyname,&mode,comment,&status)) printerror(status);	
	strcpy(keyname,"B_CF");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&b_cF,comment,&status)) printerror(status);
	strcpy(keyname,"C_CF");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&c_cF,comment,&status)) printerror(status);

	return status;
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
 * Iteration: This function takes the optimum number of rows to read and iterate them.
 *
 * - Declare variables
 * - Allocate input GSL vectors
 * - Read iterator
 * - Processing each row of pulses
 *   * Assign positive polarity to the pulses (by using ASQUID and PLSPOLAR)
 *   * Low-pass filtering and first derivative
 *   * Rows (events) are processed in only one step (without cutting into several segments)
 * - Free allocate of GSL vectors
 ****************************************************************************/
int inDataIterator(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	string message;
	int status = EPOK;
	int extver=0;

	/////////////////////// Declare variables
	static double *time;	// Vector of TIME column
	static double *v;	// Vector of V column
	nPulsesRow = 0;		// Number of pulses in each row

	/////////////////////// Allocate input GSL vectors
	gsl_vector *timegsl = gsl_vector_alloc(nrows); 				 // Input TIME column
	gsl_matrix *vgsl = gsl_matrix_alloc(nrows, eventsz); 		 // Input I0 column (matrix)
	gsl_vector *invector = gsl_vector_alloc(eventsz); 			 // Each row of vgsl
	gsl_vector *invectorNOTFILTERED = gsl_vector_alloc(eventsz); // vgsl hasn't been filtered
	gsl_vector *invectorFILTERED = gsl_vector_alloc(eventsz);	 // vgsl has been filtered (LPF)
	gsl_vector *invectorDERIVATIVE = gsl_vector_alloc(eventsz);  // Derivative of invectorFILTERED
	gsl_vector *SGN = gsl_vector_alloc(eventsz);
	gsl_vector *derSGN = gsl_vector_alloc(eventsz);				 // Sign of the invectorDERIVATIVE

	/////////////////////// Read iterator
	time = (double *) fits_iter_get_array(&cols[0]);
	status = toGslVector(((void **)&time), &timegsl, nrows, 0, TDOUBLE);
	if(status)EPexit(status);
	v = (double *) fits_iter_get_array(&cols[1]);
	status = toGslMatrix(((void **)&v), &vgsl, eventsz, nrows, (int)TDOUBLE, 0);
	if(status)EPexit(status);

	char val[256];
	char val_aux[256];

	///////////////////////  Processing each row of pulses
	for (int i=0; i< nrows; i++)
	{
		sprintf(straux,"%d",ntotalrows);
		message = "-------------> Row: " + string(straux);
		sprintf(straux,"%d",eventcnt);		
		message += " of " + string(straux) + " <------------------ ";
		writeLog(fileRef,"Log", verbosity,message);
		sprintf(val,"-------------> Row: %d of %d <------------------ ",ntotalrows,eventcnt);
		strcat(val,"\n");
		fputs(val,temporalFile);


		double time0 = gsl_vector_get(timegsl, i);	// Time in the begin of the row
		gsl_matrix_get_row(invector, vgsl, i);
		gsl_vector_scale(invector,ivcal);		//IVCAL to change arbitrary units of voltage to non-arbitrary units of current (Amps)

		// Assign positive polarity to the pulses
		if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))
		{
			gsl_vector_scale(invector,-1);
			chngplrt = 1;
		}
		strcpy(extname,"EUR-TRG");
		if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
		strcpy(keyname,"CHNGPLRT");
		if(fits_update_key(trgObject,TINT,keyname,&chngplrt,comment,&status))printerror(status);
		
		gsl_vector_memcpy(invectorNOTFILTERED,invector);

		// Low-pass filtering
		status = lpf_boxcar(&invector,invector->size,scaleFactor*tauFALL,samprate);
		if (status == 2023)
		{
	   	   writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactot too small => Cut-off frequency too high => Equivalent to not filter.");
	   	   status = EPOK;
		}
		if (status == 2024)
		{
	  	  writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
	  	  return(status);
		}
		gsl_vector_memcpy(invectorFILTERED,invector);

		// Derivative after filtering
		if(derMTHSimple (&invector,&SGN,invector->size))return (status);
		gsl_vector_memcpy(invectorDERIVATIVE,invector);
		gsl_vector_memcpy(derSGN,SGN);

		if (indice+1 == 1)
		{
			//  Creating NTFltRow Column
			if (status == EPOK) {
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
			}
			status = writeFitsSimple(obj, invectorNOTFILTERED);
			if(status)EPexit(status);

			//  Creating FltRow Column
			strcpy(obj.nameCol,"FltRow");
			status = writeFitsSimple(obj, invectorFILTERED);
			if(status)EPexit(status);

			//  Creating DerRow Column
			strcpy(obj.nameCol,"DerRow");
			status = writeFitsSimple(obj, invectorDERIVATIVE);
			if(status)EPexit(status);
		}

		// Rows (events) are processed in only one step (without cutting into several segments)
		initialtime = time0;

		status = procRow(invectorNOTFILTERED, invectorDERIVATIVE, &nPulsesRow);
		if(status)EPexit(status);
		ntotalrows ++;

		indice++;
	}

	/////////////////////// Free allocate of GSL vectors
	gsl_vector_free(timegsl);
	gsl_matrix_free(vgsl);
	gsl_vector_free(invector);
	gsl_vector_free(invectorNOTFILTERED);
	gsl_vector_free(invectorFILTERED);
	gsl_vector_free(invectorDERIVATIVE);
	gsl_vector_free(SGN);
	gsl_vector_free(derSGN);

	return (status);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
 * procRow function:  This function processes each row (event) of the input FITS file (I0 column)
 *
 * - Initialize variables
 * - Declare variables
 * - Allocate GSL vectors
 * - Find pulses of the row
 * - Calculate the tend of each found pulse
 * - Obtain in a precise way the tstart of the founded pulses
 *      - 'for' loop for each pulse
 *          - Extract the current pulse from the event
 *	 		- Low-pass filtering by using sF = 0.1scaleFactor
 * 			- Calculate tstart0 by calculating the maximum of the derivative of the filtered data
 * 		- 'for' loop for different cut frequencies of the low-pass filter (sFj=0.2scaleFactor, 0.3scaleFactor...)
 * 		 	- Low-pass filtering by using sFj
 * 		 	- Calculate tstartj by calculating the maximum of the derivative of the filtered data
 * 		 	- Compare tstart0 and tstartj: if tstart0 and tstartj are similar => Stop ('break') and tstart=(tstart0+tstartj)/2
 *                                         else => Another iteration sFj=sFj+sF*j
 * - If 'getTaus' is true => Obtain Tau of rise slope and fall slope of each pulse
 * - Write each pulse into output FITS file (call writePulses)
 * - Free allocate of GSL vectors
 *
 * Function parameters:
 * - vectorNOTFIL: I0 not filtered
 * - vectorDER: Derivative of vectorFIL
 * - npulses: Number of pulses found in the row (event)
 * - status: Auxiliary variable which has the function status: Error, Warning or OK
 ****************************************************************************/
int procRow(gsl_vector *vectorNOTFIL, gsl_vector *vectorDER, int *npulses)
{
	char val[256];
	char val_aux[256];
	int status = EPOK;

	/////////////////////// Initialize variables
	nPulsesRow = 0;

	/////////////////////// Allocate GSL vectors
	// To look for pulses
	gsl_vector *tstartgsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *tendgsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *qualitygsl = gsl_vector_alloc(vectorDER->size);
	gsl_vector *energygsl = gsl_vector_alloc(vectorDER->size);
	//gsl_vector *tendDERgsl = gsl_vector_alloc(vectorDER->size);
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

	/////////////////////// Find pulses of the row
	if(findPulses (vectorNOTFIL, vectorDER, &tstartgsl, &qualitygsl, &energygsl,
			npulses,
			mode,
			tauFALL, scaleFactor, sizePulse_b, samprate,
			samplesUp, nSgms,
			Lb, Lrs, b_cF, c_cF,
			library, models,
			safetyMarginTstart,
			stopCriteriaMKC,
			kappaMKC,
			limitCriteriaMAX,
			nsAftrtstart,
			levelPrvPulse,
			primaryThresholdCriteria,
			temporalFile,
			indice+1))EPexit(status);

	sprintf(val,"%d %d",indice+1,nPulsesRow);
	strcat(val,"\n");
	fputs(val,temporalFile2);

	// Calculate the tend of each found pulse
	gsl_vector *tendAUXgsl;

	for (int i=0;i<nPulsesRow;i++)
	{
		/*sprintf(val,"Pulsesfound( ");
		sprintf(val_aux,"%d",i);
		strcat(val,val_aux);
		sprintf(val_aux,"): ");
		strcat(val,val_aux);
		sprintf(val_aux,"%e",gsl_vector_get(tstartgsl,i));
		strcat(val,val_aux);
		strcat(val,"\n");
		fputs(val, temporalFile);*/

		if (i==0)	tendAUXgsl = gsl_vector_alloc(nPulsesRow);

		gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i)+sizePulse_b-1); 	//tend_i = tstart_i + (ntaus*tauFALL*samprate)-1

		if ((gsl_vector_get(qualitygsl,i) != 1) || (gsl_vector_get(qualitygsl,i) != 3))
		// If it is already truncated at the beginning, it is not taken into account to classify it again as truncated (at the end)
		{
			if (gsl_vector_get(tendgsl,i) >= vectorDER->size)	// Truncated pulses at the end of the row
			{
				gsl_vector_set(tendgsl,i,(vectorDER->size)-1);
				gsl_vector_set (qualitygsl,i,gsl_vector_get(qualitygsl,i)+1);
				//status = writeLog(fileRef,"Warning", verbosity,status,"Truncated at the end");
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

		if (gsl_vector_get(tendgsl,i) >= vectorDER->size)	// Truncated pulses at the end of the row
		{
			gsl_vector_set(tendgsl,i,(vectorDER->size)-1);
			gsl_vector_set (qualitygsl,i,gsl_vector_get(qualitygsl,i)+1);
		}
		if ((nPulsesRow !=1) && (i != nPulsesRow-1)) // More than one pulse in the event and not the last one
		{
			if (gsl_vector_get(tendgsl,i) > gsl_vector_get(tstartgsl,i+1))
			{
				gsl_vector_set(tendgsl,i,gsl_vector_get(tstartgsl,i+1));
			}
		}
	}

	if (getTaus == 1)
	{
		/////////////////////// Obtain Tau of rise slope and fall slope of each pulse
		// The tstarts and tends which are going to be used do not take into account the safety margin before the tstart
		gsl_vector *safetyMarginVector = gsl_vector_alloc(tstartgsl->size);
		gsl_vector_set_all(safetyMarginVector,1.0);
		gsl_vector_scale(safetyMarginVector,safetyMarginTstart);
		gsl_vector_memcpy(tstartWOutSftMrggsl,tstartgsl);
		gsl_vector_memcpy(tendWOutSftMrggsl,tendgsl);
		gsl_vector_add(tstartWOutSftMrggsl,safetyMarginVector);
		gsl_vector_add(tendWOutSftMrggsl,safetyMarginVector);
		status = obtainTau (vectorNOTFIL, tstartWOutSftMrggsl, tendWOutSftMrggsl, nPulsesRow, &tauRisegsl, &tauFallgsl);
		if(status)EPexit(status);
	}

	/////////////////////// Write pulses and measurements in output FITS files
	if (nPulsesRow != 0)	pulsesgsl = gsl_matrix_alloc(nPulsesRow,sizePulse_b);
	status = writePulses (vectorNOTFIL, vectorDER, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, energygsl, &pulsesgsl);
	if(status)EPexit(status);

	/////////////////////// Free allocate of GSL vectors
	gsl_vector_free(tstartgsl);
	gsl_vector_free(tendgsl);
	gsl_vector_free(qualitygsl);
	gsl_vector_free(energygsl);

	gsl_vector_free(tstartWOutSftMrggsl);
	gsl_vector_free(tendWOutSftMrggsl);
	gsl_vector_free(tauRisegsl);
	gsl_vector_free(tauFallgsl);

	if (nPulsesRow != 0)	gsl_matrix_free (pulsesgsl);

	return status;
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
 * ObtainTau function: This function obtains the tau value of rise and fall slopes of each pulse.
 *
 * The method to obtain each tau value is to transform the input array values.
 * The result of the transformation is a straight line, so it is easier to make an analysis and to obtain its tau values.
 *
 * Arguments:
 * 		- invector: Array of current values of an event
 * 		- tstartgsl: Start indexes of the pulses
 * 		- tendgsl: End indexes of the pulses
 * 		- nPulses: Number of pulses in the event
 * 		- taurisegsl: Tau value of rise slope of pulse (I/O parameter)
 * 		- taufallgsl: Tau value of fall slope of pulse (I/O parameter)
 * 		- status: Auxiliary variable which has the function status: Error, Warning or OK.
 *
 * Steps:
 *      - Declare variables
 * 		- Initial values
 * 		- Obtained pulse and processed it
 * 		- Obtained tau rise
 * 		- Obtained tau fall
 ******************************************************************************/
int obtainTau (gsl_vector *invector, gsl_vector *tstartgsl, gsl_vector *tendgsl, int nPulses, gsl_vector **taurisegsl, gsl_vector **taufallgsl)
{
	int status = EPOK;

	/////////////////////// Declare variables
	int tstart_index, tend_index;			// tstart and tend of each pulse
	gsl_vector *pulse;						// Each pulse
	int max_index;							// Bin where the maximum of the input value is found
											// It is the limit between rise and fall slopes
	gsl_vector_view temp;

	/////////////////////// Obtained pulse
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

		/////////////////////// Tau Rise
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

		/////////////////////// Tau Fall
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

	return (status);
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* writePulses function: It writes each pulse data into the EUR-TRG extension of the output FITS file
*
* - Creating Tstart Column
* - Creating I0 Column (non-filtered pulses)
* - Creating Tend Column
* - Creating tauRise Column
* - Creating tauFall Column
* - Creating difTstrt Column
* - Creating EstEnrgy Column
* - Creating Quality Column
*
******************************************************************************/
int writePulses(gsl_vector *invectorNOTFIL, gsl_vector *invectorDER, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, gsl_vector *energy, gsl_matrix **pulses)
{
	int status = EPOK;
	char val[256];

	/////////////////////// TRG Extension
	int t0;   			// First value of index of pulse
	gsl_matrix *vgslout2;
	gsl_vector *difTstart; // Difference between the tstart of two consecutive pulses
	                       // First pulse row difTtsart will be -1
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

		// Determine if there is the tail of a previous pulse at the beginning of the row
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
			else 									//The invectorNOTFIL has less bins than sizePulse (truncated)
			{
				for (int j=0; j<(invectorNOTFIL->size) - t0; j++)
				{
					if (t0 == -1) t0 = 0;
					gsl_matrix_set (vgslout2,i,j,gsl_vector_get(invectorNOTFIL,j+t0));
				}

				for (int j=(invectorNOTFIL->size)-t0; j< sizePulse_b; j++) {gsl_matrix_set (vgslout2,i,j,0.0);}
			}
		}

		//  Creating Tstart Column
		if (status == EPOK) {
			obj.inObject = trgObject;
			obj.nameTable = new char [255];
			strcpy(obj.nameTable,"EUR-TRG");
			obj.iniRow = totalpulses;
			obj.endRow = totalpulses+nPulsesRow-1;
			obj.iniCol = 0;
			obj.nameCol = new char [255];
			strcpy(obj.nameCol,"Tstart");
			obj.type = TDOUBLE;
			obj.unit = new char [255];
			strcpy(obj.unit,"seconds");
		}
		temp = gsl_vector_subvector(tstart,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector);
		if(status)EPexit(status);

		// Creating I0 Column
		if (writePulse){	// If the input parameter writePulses = yes then it puts in output FITS file
			strcpy(obj.nameCol,"I0");
			obj.type = TDOUBLE;
			strcpy(obj.unit,"Amps");
			if(writeFitsComplex(obj, vgslout2))EPexit(status);
		}
		//gsl_matrix_get_row(*pulsetemplate,vgslout2,0);
		gsl_matrix_memcpy(*pulses,vgslout2);

		// Creating Tend Column
		strcpy(obj.nameCol,"Tend");
		obj.type = TDOUBLE;
		strcpy(obj.unit,"seconds");
		temp = gsl_vector_subvector(tend,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector);
		if(status)EPexit(status);

		// Creating tauRise Column
		strcpy(obj.nameCol,"tauRise");
		temp = gsl_vector_subvector(taurise,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector);
		if(status)EPexit(status);

		// Creating tauFall Column
		strcpy(obj.nameCol,"tauFall");
		temp = gsl_vector_subvector(taufall,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector);
		if(status)EPexit(status);

		// Creating difTstrt Column
		strcpy(obj.nameCol,"difTstrt");
		status = writeFitsSimple(obj, difTstart);
		if(status)EPexit(status);

		// Creating EstEnrgy Column
		strcpy(obj.nameCol,"EstEnrgy");
		strcpy(obj.unit,"eV");
		temp = gsl_vector_subvector(energy,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector);
		if(status)EPexit(status);

		// Creating Quality Column
		strcpy(obj.nameCol,"Quality");
		obj.type = TSHORT;
		strcpy(obj.unit,"bits");
		temp = gsl_vector_subvector(quality,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector);
		if(status)EPexit(status);

		// Creating TailBfr Column
		strcpy(obj.nameCol,"TailBfr");
		obj.type = TSHORT;
		strcpy(obj.unit," ");
		status = writeFitsSimple(obj, tailBeforegsl);
		if(status)EPexit(status);

		totalpulses = totalpulses + nPulsesRow;
		sprintf(val,"totalpulses: %d",totalpulses-1);
		strcat(val,"\n");
		fputs(val,temporalFile);

		/////////////////////// Free allocate GSL vectors
		gsl_matrix_free (vgslout2);
		gsl_vector_free (difTstart);
		gsl_vector_free (tailBeforegsl);
	}

	return (status);
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/




/***** SECTION 10 ************************************************************
* readLib: This function loads all the pulse models from the pulses models library input FITS file.
******************************************************************************/
int readLib ()
{
	string message="";
	int extver=0, status = EPOK;
	
	strcpy(extname,"EUR-LIB");
	if(fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status)) printerror(status);
	if(fits_get_num_rows(inLibObject,&nummodels, &status)) printerror(status);
	
	/////////////////////// Get structure of pulse model library input FITS file columns
	strcpy(straux,"ENERGY");
	if(fits_get_colnum(inLibObject,0,straux,&colnum,&status)) printerror(status);
	strcpy(straux,"ESTENERGY");
	if(fits_get_colnum(inLibObject,0,straux,&colnum,&status)) printerror(status);
	strcpy(straux,"PULSE");
	if(fits_get_colnum(inLibObject,0,straux,&colnum,&status)) printerror(status);		

	extern int inDataIteratorLib(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

		// Create structure to run Iteration DAL
	iteratorCol colsLib [3]; 				// Structure of Iteration DAL
	int n_cols = 3; 					// Number of columns:  Energy + EstEnergy + PULSE
	long rows_per_loop = nummodels; 			// 0: Use default: Optimum number of rows
	long offset=0; 						// 0: Process all the rows

	// Read Energy Column
	strcpy(straux,"Energy");
	if(fits_iter_set_by_name(&colsLib[0], inLibObject, straux, TDOUBLE, InputCol))printerror(status);
		// Read Energy Column
	strcpy(straux,"EstEnergy");
	if(fits_iter_set_by_name(&colsLib[1], inLibObject, straux, TDOUBLE, InputCol))printerror(status);
		// Read PULSE
	strcpy(straux,"PULSE");
	if(fits_iter_set_by_name(&colsLib[2], inLibObject, straux, TDOUBLE, InputCol))printerror(status);

		// Called iteration function
	if(fits_iterate_data(n_cols, colsLib, offset, rows_per_loop, inDataIteratorLib,0L,&status)) printerror(status);

	return (status);
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
 * Iteration: This function takes the optimum number of rows to read and iterate them.
 *
 * - Declare variables
 * - Allocate input GSL vectors
 * - Read iterator
 * - Processing each row of pulses
 * - Free allocate of GSL vectors
 ****************************************************************************/
int inDataIteratorLib(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
  	int status = EPOK;
	int extver=0;

	/////////////////////// Declare variables
	static double *energy;		// Vector of ENERGY column
	static double *estenergy;	// Vector of ESTENERGY column
	static double *pulsemodel;  	// Vector of PULSE column
	long eventsz;

	strcpy(extname,"EUR-LIB");
	if(fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status)) printerror(status);
	strcpy(keyname,"EVENTSZ");
	if(fits_read_key(inLibObject,TLONG,keyname, &eventsz,comment,&status)) printerror(status);
	
	/////////////////////// Allocate input GSL vectors
	gsl_vector *energygsl = gsl_vector_alloc(nummodels);
	gsl_vector *estenergygsl = gsl_vector_alloc(nummodels);
	gsl_matrix *pulsemodelgsl = gsl_matrix_alloc(nummodels,eventsz);
	library = gsl_matrix_alloc(nummodels,2);		// Energy + EstEnergy
	models = gsl_matrix_alloc(nummodels,eventsz);
	gsl_vector *rowaux = gsl_vector_alloc(eventsz);	// Auxiliary local variable
	model = gsl_vector_alloc(eventsz);

	/////////////////////// Read iterator
	energy = (double *) fits_iter_get_array(&cols[0]);
	status = toGslVector(((void **)&energy), &energygsl, nrows, 0, TDOUBLE);
	if(status)EPexit(status);
	estenergy = (double *) fits_iter_get_array(&cols[1]);
	status = toGslVector(((void **)&estenergy), &estenergygsl, nrows, 0, TDOUBLE);
	if(status)EPexit(status);
	pulsemodel = (double *) fits_iter_get_array(&cols[2]);
	status = toGslMatrix(((void **)&pulsemodel), &pulsemodelgsl, eventsz, nrows, (int)TDOUBLE,0);
	if(status)EPexit(status);

	///////////////////////  Processing each row of pulses
	for (int i=0; i< nrows; i++)
	{
		gsl_matrix_set(library,ntotalrows-1,0,gsl_vector_get(energygsl,i));
		gsl_matrix_set(library,ntotalrows-1,1,gsl_vector_get(estenergygsl,i));

		gsl_matrix_get_row(rowaux,pulsemodelgsl,i);
		gsl_matrix_set_row(models,ntotalrows-1,rowaux);

		ntotalrows ++;
	}

	/////////////////////// Free allocate of GSL vectors
	gsl_vector_free(energygsl);
	gsl_vector_free(estenergygsl);
	gsl_matrix_free(pulsemodelgsl);
	gsl_vector_free(rowaux);

	return (status);
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
 * createLibrary: This function
 *
 * - xxx
 * - xxx
 ****************************************************************************/
int createLibrary()
{
	int extver=0, status=EPOK;

	/////////////////////// Create output FITS files: If it does not exist already
	//Create Library file or open it
	fits_open_file(&inLibObject, inLibName,0,&status);

	if (status == EPOK)
	{
		append = true;

		strcpy(extname,"EUR-LIB");
		if(fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status)) printerror(status);
		strcpy(keyname,"EVENTCNT");
		if(fits_read_key(inLibObject,TLONG,keyname, &eventcntLib,comment,&status)) printerror(status);
		eventcntLib1 = eventcntLib + 1;
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(inLibObject,TLONG,keyname, &eventcntLib1,comment,&status))printerror(status);
		
		status = readLib();
		if(status)EPexit(status);
		ntotalrows = 1;
	}
	else
	{
		append = false;
		status = EPOK;
		if(fits_create_file(&inLibObject, inLibName, &status))printerror(status);

		/////////////////////// Create extension EUR-LIB
		strcpy(extname,"EUR-LIB");
		if(fits_create_tbl(inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status)) printerror(status);
		evtcnt = 1;
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(inLibObject,TLONG,keyname,&evtcnt,comment,&status)) printerror(status);
		/////////////////////// Create keywords
		lib_id = 1;
		strcpy(keyname,"LIB_ID");
		if(fits_write_key(inLibObject,TLONG,keyname,&lib_id,comment,&status)) printerror(status);
		strcpy(keyname,"CREATOR");
		strcpy(keyvalstr,create);
		if(fits_write_key(inLibObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);
		strcpy(keyname,"FTYPE");
		strcpy(keyvalstr,"lib");
		if(fits_write_key(inLibObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);
	}

	char str_energy[125];       sprintf(str_energy,"%f",energy);

	string process (string("TRIGGER") 	+ ' ' +
	string(inName) 		+ ' ' + string(inLibName) 	  + ' ' +
	string(str_energy)      + ' ' +
	string("(")				+      (string) create 		  +   	  string(")"));

	strcpy(keyname,"PROCESS");
	strcpy(keyvalstr,process.c_str());
	if(fits_write_key_longwarn(inLibObject,&status))printerror(status);
	if(fits_write_key_longstr(inLibObject,keyname,keyvalstr,comment,&status)) printerror(status);
	//close file if newly created
	if(!append){
	  if (fits_close_file(inLibObject, &status))printerror(status);
	}
	return (status);
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* writeLibrary function: It writes each pulse model into the EUR-LIB extension of the output FITS file
*
******************************************************************************/
int writeLibrary(double estenergy, gsl_vector **pulsetemplate)
{
	int status = EPOK;
	/////////////////////// LIB Extension
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
    	strcpy(obj.nameTable,"EUR-LIB");
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
    	    status = writeFitsSimple(obj, energyoutgsl);
	    if(status)EPexit(status);

    	    strcpy(obj.nameCol,"ESTENERGY");
    	    strcpy(obj.unit," ");
    	    gsl_vector_set (estenergyoutgsl,0,gsl_vector_get(estenergycolumn,i));
    	    status = writeFitsSimple(obj, estenergyoutgsl);
	    if(status)EPexit(status);

    	    strcpy(obj.nameCol,"PULSE");
    	    strcpy(obj.unit,"Amps");
    	    gsl_matrix_get_row(modelsrow,modelsaux,i);
    	    gsl_matrix_set_row(pulsetemplates_matrix,0,modelsrow);
    	    status = writeFitsComplex(obj, pulsetemplates_matrix);
	    if(status)EPexit(status);
    	}
    	gsl_vector_free(modelsrow);

    	gsl_vector_free (energycolumn);
    	gsl_vector_free (estenergycolumn);
    	gsl_matrix_free(modelsaux);
    }
    else
    {
	if(fits_write_key(inLibObject,TINT,"EVENTSZ",&sizePulse_b,comment,&status)) printerror(status);

    	gsl_vector_set (energyoutgsl,0,energy);
    	gsl_vector_set (estenergyoutgsl,0,estenergy);

    	gsl_matrix_set_row(pulsetemplates_matrix,0,*pulsetemplate);

    	// Creating ENERGY Column
   	obj.inObject = inLibObject;
   	obj.nameTable = new char [255];
   	strcpy(obj.nameTable,"EUR-LIB");
   	obj.iniRow = 1;
   	obj.endRow = 1;
   	obj.iniCol = 0;
   	obj.nameCol = new char [255];
   	strcpy(obj.nameCol,"ENERGY");
   	obj.type = TDOUBLE;
   	obj.unit = new char [255];
   	strcpy(obj.unit,"eV");
   	status = writeFitsSimple(obj, energyoutgsl);
	if(status)EPexit(status);

   	// Creating ESTENERGY Column
    	strcpy(obj.nameCol,"ESTENERGY");
    	strcpy(obj.unit,"u1");
    	status = writeFitsSimple(obj, estenergyoutgsl);
	if(status)EPexit(status);

   	// Creating PULSE Column
   	strcpy(obj.nameCol,"PULSE");
   	strcpy(obj.unit,"Amps");
   	status = writeFitsComplex(obj, pulsetemplates_matrix);
	if(status)EPexit(status);
    }

    // Free memory
    gsl_vector_free(energyoutgsl);
    gsl_vector_free(estenergyoutgsl);

    gsl_matrix_free (pulsetemplates_matrix);

    return (status);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhistogsl, gsl_vector **yhistogsl)
{
    int status = EPOK;

    /////////////////////// Declare variables
    int size = invector->size;
    double invectormax= 0;         // Maximum of invector
    double invectormin=1e10;  // Minimum of invector
    double binSize;
    int ind = 0;                       // Index of the bin which contains each invector element
    gsl_vector *invectoraux = gsl_vector_alloc(size);
    gsl_vector *invectoraux2;
    gsl_vector_view temp;					// In order to handle with gsl_vector_view (subvectors)

    char val[256];

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

    ///Obtain invector_max
    for (int i=0; i<size; i++)
    {
    	if (invectormax < gsl_vector_get (invectoraux2,i))	invectormax = gsl_vector_get (invectoraux2,i);
    }
    // Obtain invector_min
    for (int i=0; i<size; i++)
    {
    	if (invectormin > gsl_vector_get (invectoraux2,i))	invectormin = gsl_vector_get (invectoraux2,i);
    }

    if (invectormax == invectormin) EPexit(1013);

    if ((invectormax != 0) && (invectormin != 1e10))
    {

    	// Obtain binSize
        binSize = (invectormax - invectormin) / nbins;        // Size of each bin of the histogram

        /////////////////////// Create histogram x-axis
        for (int i=0; i<nbins; i++)	{gsl_vector_set (*xhistogsl,i,binSize*i+invectormin+binSize/2);}

        /////////////////////// Create histogram y-axis
        gsl_vector_set_zero (*yhistogsl);
        for (int i=0; i<size; i++)
        {
        	ind = (int) ((gsl_vector_get(invectoraux2,i) - invectormin) / binSize);
            if (ind == nbins) ind--;
            gsl_vector_set (*yhistogsl,ind, gsl_vector_get(*yhistogsl,ind) + 1);
        }
    }

    gsl_vector_free(invectoraux);
    gsl_vector_free(invectoraux2);

    return status;
}


int inDataIteratorOutTrg(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	/////////////////////// Declare variables
	int status = EPOK;
	static double *tstart;		// Vector of TSTART column
	static double *estenrgy;	// Vector of ENERGY column
	static double *quality;		// Vector of QUALITY column

	/////////////////////// Allocate input GSL vectors
	gsl_vector *tstartoutgsl = gsl_vector_alloc(totalpulses-1);
	gsl_vector *estenrgygsl = gsl_vector_alloc(totalpulses-1);
	gsl_vector *qualityoutgsl = gsl_vector_alloc(totalpulses-1);

	/////////////////////// Read iterator
	tstart = (double *) fits_iter_get_array(&cols[0]);
	status = toGslVector(((void **)&tstart), &tstartoutgsl, nrows, 0, TDOUBLE);
	if(status)EPexit(status);
	estenrgy = (double *) fits_iter_get_array(&cols[1]);
	status = toGslVector(((void **)&estenrgy), &estenrgygsl, nrows, 0, TDOUBLE);
	if(status)EPexit(status);
	quality = (double *) fits_iter_get_array(&cols[2]);
	status = toGslVector(((void **)&quality), &qualityoutgsl, nrows, 0, TDOUBLE);
	if(status)EPexit(status);

	///////////////////////  Processing each row of pulses
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set(tstartout,ntotalrows-1,gsl_vector_get(tstartoutgsl,i));
		gsl_vector_set(pulseheight,ntotalrows-1,gsl_vector_get(estenrgygsl,i));
		gsl_vector_set(qualityout,ntotalrows-1,gsl_vector_get(qualityoutgsl,i));

		ntotalrows ++;
	}

	/////////////////////// Free allocate of GSL vectors
	gsl_vector_free(tstartoutgsl);
	gsl_vector_free(estenrgygsl);
	gsl_vector_free(qualityoutgsl);

	return (status);
}

int align(gsl_vector **vector1, gsl_vector ** vector2)
{
	int status = EPOK;
	const double pi = 4.0 * atan(1.0);

	char val[256];

	int size = (*vector1)->size;

	double SelectedTimeDuration = size/samprate;
	gsl_vector_complex *vector1fft = gsl_vector_complex_alloc(size);
	gsl_vector_complex *vector2fft = gsl_vector_complex_alloc(size);
	double vector1fft_ph;
	double vector2fft_ph;

	double shiftdouble;
	int shiftint;

	gsl_vector *vector2shifted = gsl_vector_alloc(size);

	status = FFT(*vector1,vector1fft,SelectedTimeDuration);
	if(status)EPexit(status);
	status = FFT(*vector2,vector2fft,SelectedTimeDuration);
	if(status)EPexit(status);

	vector1fft_ph= gsl_complex_arg(gsl_vector_complex_get(vector1fft,1))*size/(2*pi);
	vector2fft_ph= gsl_complex_arg(gsl_vector_complex_get(vector2fft,1))*size/(2*pi);
	shiftdouble = vector1fft_ph-vector2fft_ph;

	if ((shiftdouble > -1) && (shiftdouble < 1)) shiftint = 0;
	else if (shiftdouble > 1)	shiftint = floor(shiftdouble);
	else if (shiftdouble < -1)	shiftint = ceil(shiftdouble);


	if (shiftint > 0)
	{
		status = shift_m(*vector2,vector2shifted,shiftint);
		if(status)EPexit(status);
		gsl_vector_memcpy(*vector2,vector2shifted);
	}
	else if (shiftint < 0)
	{
		status = shiftm(*vector2,vector2shifted,fabs(shiftint));
		if(status)EPexit(status);
		gsl_vector_memcpy(*vector2,vector2shifted);
	}

	/*sprintf(val,"ph1: %.12e ph2: %.12e shift: %.12e shiftint: %d",vector1fft_ph,vector2fft_ph,shiftdouble,shiftint);
	strcat(val,"\n");
	fputs(val,temporalFile);*/

	return (status);
}


int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	int status = EPOK;
	char val[256];
	int size = vectorin->size;

	for (int i=0;i<size-m;i++)
	{
		gsl_vector_set(vectorout,i+m,gsl_vector_get(vectorin,i));
	}
	for (int i=0;i<m;i++)
	{
		gsl_vector_set(vectorout,i,gsl_vector_get(vectorin,0));
	}

	return (status);
}

int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m)
{
	int status = EPOK;
	char val[256];
	int size = vectorin->size;

	for (int i=m;i<size;i++)
	{
		gsl_vector_set(vectorout,i-m,gsl_vector_get(vectorin,i));
	}
	for (int i=size-m;i<size;i++)
	{
		gsl_vector_set(vectorout,i,gsl_vector_get(vectorin,size-1));
	}

	return (status);
}

int calculateTemplate (int totalPulses, gsl_vector **pulseaverage, double *pulseaverageHeight)
{
	int status = EPOK;

	char val[256];
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

	//gsl_vector *pulseaverage = gsl_vector_alloc(sizePulse_b);	// Template (calculated by averaging pulses)
	bool firstnonpileupPulse = true;
	gsl_matrix *pulse_matrix = gsl_matrix_alloc(1,sizePulse_b);
	gsl_vector *pulse = gsl_vector_alloc(sizePulse_b);
	gsl_vector *pulseshifted = gsl_vector_alloc(sizePulse_b);
	gsl_vector *pulseoriginal = gsl_vector_alloc(sizePulse_b);	// Only to info purposes (TO BE DELETED)

	int Shift= 0;
	double tstartnext;

	// When one of the found pulses is validated as non piled up pulse by using EstEnrgy
	strcpy(obj.nameTable,"EUR-TRG");    // (pulseheights histogram) and Tstart, and Quality => The pulse, I0, is going
	strcpy(obj.nameCol,"I0");           // to be read from the trg output FITS file (in order to not handle a long matrix of
	                                    // pulses, found pulses x sizePulse_b)

	/////////////////////// Get structure of output FITS file columns
	//long numBins = 0L;	dal_dataType type;	int varLength=0; long sizeT, sizeV;
	if(fits_movnam_hdu(obj.inObject, ANY_HDU,obj.nameTable, extver, &status)) printerror(status);
	strcpy(straux,"Tstart");
	if(fits_get_colnum(obj.inObject,0,straux,&colnum,&status)) printerror(status);
	//status = DALtableGetColStruct (trgExten,0,&numBins,straux,&type, &sizeT,&varLength,status);	if (status == -1219){status = -62511; EPexit(status);}
	strcpy(straux,"EstEnrgy");
	if(fits_get_colnum(obj.inObject,0,straux,&colnum,&status)) printerror(status);
	//status = DALtableGetColStruct (trgExten,0,&numBins,straux,&type, &sizeT,&varLength,status);	if (status == -1219){status = -62512; EPexit(status);}
	strcpy(straux,"Quality");
	if(fits_get_colnum(obj.inObject,0,straux,&colnum,&status)) printerror(status);
	//status = DALtableGetColStruct (trgExten,0,&numBins,straux,&type, &sizeT,&varLength,status);	if (status == -1219){status = -62513; EPexit(status);}

	///////////////////////// Iteration DAL
	extern int inDataIteratorOutTrg(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

		// Create structure to run Iteration DAL
	iteratorCol cols [3]; 					// Structure of Iteration DAL
	int n_cols = 3; 						// Number of columns:  Tstart + EstEnrgy + Quality
	long rows_per_loop = 0; 				// 0: Use default: Optimum number of rows
	long offset=0; 							// 0: Process all the rows

	strcpy(straux,"Tstart");
	//status = DALiterSetByName(&cols[0], trgExten, straux, DAL_DOUBLE, InputCol,	status);
	if(fits_iter_set_by_name(&cols[0], obj.inObject, straux, TDOUBLE, InputCol)) printerror(status);
	strcpy(straux,"EstEnrgy");
	//status = DALiterSetByName(&cols[1], trgExten, straux, DAL_DOUBLE, InputCol,	status);
	if(fits_iter_set_by_name(&cols[1], obj.inObject, straux, TDOUBLE, InputCol)) printerror(status);
	strcpy(straux,"quality");
	//status = DALiterSetByName(&cols[2], trgExten, straux, DAL_DOUBLE, InputCol,	status);
	if(fits_iter_set_by_name(&cols[2], obj.inObject, straux, TDOUBLE, InputCol)) printerror(status);


	ntotalrows = 1;	// Because it has already been used in other inDataIterators functions
	//status = DALiterateData(n_cols, cols, offset, rows_per_loop, inDataIteratorOutTrg, 0L, status);
	if(fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIterator,0L,&status))printerror(status);

	gsl_vector_scale(tstartout,samprate); 	//tstarts not in sec but in samples

	status = createHisto(pulseheight, nBins, &xhisto, &yhisto);
	if(status)EPexit(status);
	index_maximumpulseheight = gsl_vector_max_index(yhisto);
	maximumpulseheight = gsl_vector_get(xhisto,index_maximumpulseheight);

	for (int i=0;i<totalPulses-1;i++)
	{
		if (i == totalPulses-2)		tstartnext = gsl_vector_get(tstartout,i)+2*sizePulse_b;
		else 						tstartnext = gsl_vector_get(tstartout,i+1);

		if ((gsl_vector_get(pulseheight,i) < maximumpulseheight-0.1*maximumpulseheight) || (gsl_vector_get(pulseheight,i) > maximumpulseheight+0.1*maximumpulseheight)
			|| (tstartnext-gsl_vector_get(tstartout,i) <= sizePulse_b) || (gsl_vector_get(qualityout,i) != 0))
		{
			gsl_vector_set(nonpileup,i,0);
			nonpileupPulses --;
		}
		else
		{
			if (firstnonpileupPulse == true)
			{
				obj.iniRow = i+1;
				obj.endRow = i+1;
				status=readFitsComplex (obj,&pulse_matrix);
				if(status)EPexit(status);
				/*for (int j=0;j<pulse_matrix->size2;j++)
				{
					sprintf(val,"%d %.12e",j,gsl_matrix_get(pulse_matrix,0,j));
					strcat(val,"\n");
					fputs(val,temporalFile);
				}*/
				gsl_matrix_get_row(pulse,pulse_matrix,0);
				gsl_vector_memcpy(*pulseaverage,pulse);
				*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
			}
			else
			{
				obj.iniRow = i+1;
				obj.endRow = i+1;
				status=readFitsComplex (obj,&pulse_matrix);
				if(status)EPexit(status);
				gsl_matrix_get_row(pulse,pulse_matrix,0);

				/*sprintf(val,"%d %d %d",i,gsl_vector_max_index(pulseaverage),gsl_vector_max_index(pulse));
				strcat(val,"\n");
				fputs(val,temporalFile);*/

				gsl_vector_memcpy(pulseoriginal,pulse);
				status=align(pulseaverage,&pulse);
				if(status)EPexit(status);
				gsl_vector_add(*pulseaverage,pulse);
				*pulseaverageHeight = *pulseaverageHeight + gsl_vector_get(pulseheight,i);
				/*for (int k=0;k<(*pulseaverage)->size;k++)
				{
					sprintf(val,"%.12e %.12e %.12e",gsl_vector_get(*pulseaverage,k),gsl_vector_get(pulseoriginal,k),gsl_vector_get(pulse,k));
					strcat(val,"\n");
					fputs(val,temporalFile);
				}*/
				/*status = shiftm(pulse,pulseshifted,i,status);
				status = align(&pulse,&pulseshifted,&Shift,status);
				sprintf(val,"%d",(-1)*i);
				strcat(val,"\n");
				fputs(val,temporalFile);
				status = shift_m(pulse,pulseshifted,i,status);
				status = align(&pulse,&pulseshifted,&Shift,status);*/
			}
			if (firstnonpileupPulse == true)	firstnonpileupPulse = false;
		}
	}

	gsl_vector_scale(*pulseaverage,1.0/(nonpileupPulses));
	*pulseaverageHeight = *pulseaverageHeight/nonpileupPulses;

	gsl_vector_free(tstartout);
	gsl_vector_free(pulseheight);
	gsl_vector_free(qualityout);

	return (status);
}
