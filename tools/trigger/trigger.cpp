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
*                                - If a pulse if pile-up in the tail of a previous one => The process to recalculate precisely the
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
	extern int inDataIterator(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

	
	sprintf(temporalFileName,"TRIGGERauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if(temporalFile == NULL)
	{
		writeLog(fileRef,"Error",verbosity,"Cannot open auxiliary file TRIGGERauxfile.txt");
		EPexit(-62900);
	}
	sprintf(temporalFileName2,"ppr");
	strcat(temporalFileName2,".txt");
	temporalFile2 = fopen (temporalFileName2,"w");
	if(temporalFile2 == NULL)
	{
		writeLog(fileRef,"Error",verbosity,"Cannot open auxiliary file ppr.txt");
		EPexit(-62900);
	}

	// Read input parameters
	status = initModule(argc, argv, status);
	if (!status == EPOK)  EPexit(status);

	writeLog(fileRef,"Log", verbosity, "Into Trigger task");
  
	////////////////////////  	Open input FITS file
	//MC status = DALobjectOpen(inName, &inObject, status); 				if (status == -2004){EPexit(-62102);}
	//status = DALobjectFindElement(inObject, "RAW_ADC_DATA", &inExten, status);	if (status == -2104){EPexit(-62201);}
	//MC status = DALobjectFindElement(inObject, "RECORDS", &inExten, status);	if (status == -2104){EPexit(-62201);}
	if(fits_open_file(&inObject, inName,0,&status)) printerror(status);
	strcpy(extname,"RECORDS");
	if(fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status)) printerror(status);	
	
	///////////////////////	Input keywords
	/*getKg (&kg, &kgNumber, status);
	if (!status == EPOK) EPexit(status);
	getKi (&ki, &kiNumber, status);
	if (!status == EPOK) EPexit(status);
	status = readKeywords(kg, kgNumber, inName, inExten, fileRef,  status);
	if (!status == EPOK) EPexit(status);
	status = readKeywords(ki, kiNumber, inName, inExten, fileRef,  status);
	if (!status == EPOK) EPexit(status);
	status = goodValuesKeywords (kg, kgNumber, inName, inExten, fileRef,  status);
	if (!status == EPOK) EPexit(status);
	
	string datatype;
	status = getKeywordsValueString (kg, kgNumber, "DATATYPE", &datatype, fileRef, status);
	if ((datatype != "XRAY") & (datatype != "TESNOISE")){
	  message = "Illegal value of DATATYPE keyword in " + string(inName) + " file and " + inExten->name + " extension";
	  writeLog(fileRef,"Error",verbosity,message);
	  writeLog(fileRef,"Error", verbosity,"Legal values 'XRAY' or 'TESNOISE'");
	  status = -72002;
	  EPexit(status);
	}
	string bias;
	double oscfrekw;
	status = getKeywordsValueString (kg, kgNumber, "BIAS", &bias, fileRef, status);
	if (bias == "A"){
	  status = getKeywordsValueReal (kg, kgNumber, "OSCFREKW", &oscfrekw, fileRef, status);
	  if (oscfrekw == 0){
	    message = "Illegal value of OSCFREKW keyword in " + string(inName) + " file and " + inExten->name + "extension";
	    writeLog(fileRef,"Error",verbosity,message);
	    writeLog(fileRef,"Error", verbosity,"Legal values greater than 0");
	    status = -72002;
	    EPexit(status);
	  }
	}*/

	//MC status = DALattributeGetInt(inExten, "TRIGGSZ", &eventsz, unit, comment, status);
	strcpy(keyname,"TRIGGSZ");
	if(fits_read_key(inObject,TLONG,keyname, &eventsz,comment,&status)) printerror(status);
	//MC status = DALattributeGetReal(inExten, "DELTAT", &samprate, unit, comment, status);
	strcpy(keyname,"DELTAT");
	if(fits_read_key(inObject,TDOUBLE,keyname, &samprate,comment,&status)) printerror(status);
	samprate = 1/samprate;
	//MC status = DALtableGetNumRows(inExten,&eventcnt,status);
	if(fits_get_num_rows(inObject,&eventcnt, &status)) printerror(status);
	ivcal=1.0;
	asquid = 1.0;
	//MC status = DALattributeGetReal(inExten, "MONOEN", &energy, unit, comment, status);
	strcpy(keyname,"MONOEN");
	if(fits_read_key(inObject,TDOUBLE,keyname, &energy,comment,&status)) printerror(status);
	
	energy = energy*1e3;
	plspolar = 1.0;
	/*status = getKeywordsValueInt (kg, kgNumber, "EVENTSZ", &eventsz, fileRef, status);
	status = getKeywordsValueReal (kg, kgNumber, "SAMPRATE", &samprate, fileRef, status);
	status = getKeywordsValueInt (kg, kgNumber, "EVENTCNT", &eventcnt, fileRef, status);
	status = getKeywordsValueReal (kg, kgNumber, "IVCAL", &ivcal, fileRef, status);
	status = getKeywordsValueReal (kg, kgNumber, "ASQUID", &asquid, fileRef, status);
	status = getKeywordsValueReal (kg, kgNumber, "ENERGY", &energy, fileRef, status);
	status = getKeywordsValueReal (kg, kgNumber, "PLSPOLAR", &plspolar, fileRef, status);*/

	// If STRICT calibration mode => No secondary pulses searched for => No pulse models library used
	if ((mode == 1) || (mode == 2))
	{
		// Open pulses models library file (EUR-LIB extension)
		//MC status = DALobjectOpen(inLibName, &inLibObject, status);
		//MC if (status == -2004)
		//MC {
		//MC 	status = 0;
		//MC 	message = "Input library FITS file either not found or file name too long: " + string(inLibName);
		//MC 	writeLog(fileRef,"Error", verbosity,message);
		//MC 	EPexit(-2004);
		//MC }
		if(fits_open_file(&inLibObject, inLibName,0,&status)) printerror(status);

		//MC status = DALobjectFindElement(inLibObject, "EUR-LIB", &inLibExten, status);
		//MC if (status == -2104)
		//MC {
		//MC 	status = 0;
		//MC 	message = "EUR-LIB extension not found in " + string(inLibName) + " file";
		//MC 	writeLog(fileRef,"Error", verbosity,message);
		//MC 	EPexit(-2104);
		//MC }
		strcpy(extname,"EUR-LIB");
		if(fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status)) printerror(status);

		// Read the pulse models library input FITS file and store the whole library in "library" and "models"

		status = readLib(status);

		ntotalrows = 1;
		/*gsl_matrix *modelsAux = gsl_matrix_alloc(models->size1,models->size2+50);
	    for (int i=0;i<nummodels;i++)
	    {
	    	for (int j=0;j<50;j++)
	    	{
	    		gsl_matrix_set(modelsAux,i,j,0.0);
	    	}
	    	for (int j=50;j<modelsAux->size2;j++)
	    	{
	    		gsl_matrix_set(modelsAux,i,j,gsl_matrix_get(models,i,j-50));
	    	}
	    }
	    gsl_matrix_free(models);
	    gsl_matrix *models = gsl_matrix_alloc(modelsAux->size1,modelsAux->size2);
	    gsl_matrix_memcpy(models,modelsAux);
	    gsl_matrix_free(modelsAux);*/
	}
	else
	{
		status = createLibrary(status);
		if (!status == EPOK) EPexit (status);
		if(fits_open_file(&inLibObject,inLibName,READWRITE,&status))printerror(status);
	}

	if (energy == 0)	// There are no pulses
	{
		/////////////////////// Create output FITS file: Trigger file (*_trg.fits)
		status = createTriggerFile(status);
		if (!status == EPOK) EPexit (status);
		//MC 
		if(fits_open_file(&trgObject,trgName,READWRITE,&status))printerror(status);

		/////////////////////// Write output keywords
		//MC status = DALattributePutInt(trgExten, "EVENTCNT", 0, unit, comment,	status);
		evtcnt = 0;
		strcpy(extname,"EUR-TRG");
		if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(trgObject,TLONG,keyname,&evtcnt,comment,&status)) printerror(status);

	}
	else				// There are pulses
	{
		//////////////////////// Initialize input parameters
		sizePulse = ntaus * tauFALL;
		safetyMarginTstart = safetyMarginTstart*samprate;
		nsAftrtstart = tAftrtstart*samprate;

		// Transform from seconds to bins
		status = seconds2Bins (status);

		Lrs = (int)(LrsT*samprate);
		Lb = (int)(LbT*samprate);

		if (sizePulse_b > eventsz){
			writeLog(fileRef, "Warning", verbosity, "Standard vector size larger than Row size of the input FITS file. Resize sizePulse to Row size");
			sizePulse_b = eventsz;
		}
		if (sizePulse_b == 0){
			writeLog(fileRef, "Error", verbosity, "Calculated pulse size is 0");
			status = -62801;
			EPexit(status);
		}
		//pulsegsl = gsl_vector_alloc(sizePulse_b);

		/////////////////////// Get structure of input FITS file columns
		//MC long numBins = 0L;	dal_dataType type;	int varLength=0; long sizeT, sizeV;

		strcpy(straux,"Time");
		//MC status = DALtableGetColStruct (inExten,0,&numBins,straux,&type, &sizeT,&varLength,status);	if (status == -1219){status = -62501; EPexit(status);}
		strcpy(extname,"RECORDS");
		if(fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status)) printerror(status);
		if(fits_get_colnum(inObject,0,straux,&colnum,&status)) printerror(status);
		//strcpy(straux,"I0");
		strcpy(straux,"ADC");
		//MC status = DALtableGetColStruct (inExten,0,&numBins,straux,&type, &sizeV,&varLength,status);	if (status == -1219){status = -62502; EPexit(status);}
		if(fits_get_colnum(inObject,0,straux,&colnum,&status)) printerror(status);
		
		sprintf(str_stat,"%d",status);
		message = "Open InFits: " +  string(inName) + " " + string(str_stat);
		writeLog(fileRef,"Log",verbosity,message);


		/////////////////////// Create output FITS file: Trigger file (*_trg.fits)
		status = createTriggerFile(status);

		if (status != EPOK) EPexit (status);
		if(fits_open_file(&obj.inObject,trgName,1,&status))printerror(status);

		///////////////////////// Iteration DAL
		//MC extern int inDataIterator(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

			// Create structure to run Iteration DAL
		iteratorCol cols [2]; 					// Structure of Iteration DAL
		int n_cols = 2; 						// Number of columns:  TIME + I
		long rows_per_loop = 0; 				// 0: Use default: Optimum number of rows
		long offset=0; 							// 0: Process all the rows
		
		//MC
		strcpy(extname,"RECORDS");
		if(fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status)) printerror(status);
		
			// Read Time Column
		strcpy(straux,"TIME");
		//MC status = DALiterSetByName(&cols[0], inExten, straux, DAL_DOUBLE, InputCol,	status);
		if(fits_iter_set_by_name(&cols[0], inObject, straux, TDOUBLE, InputCol)) printerror(status);
		

			// Read I
		//strcpy(straux,"I0");
		strcpy(straux,"ADC");
		//MC status = DALiterSetByName(&cols[1], inExten, straux,DAL_DOUBLE, InputCol, status);
		if(fits_iter_set_by_name(&cols[1], inObject, straux, TDOUBLE, InputCol))printerror(status);
		

		if ((mode == 1) || (mode == 2))
		{
			// Once the pulse models read from the library are low-pass filtered and derived, they are stored in "models"
			for (int i=0; i<models->size1; i++)
			{
				gsl_matrix_get_row(model,models,i);

				// PULSE MODEL: Low-pass filtering
				status = lpf_boxcar(&model,model->size,tauFALL*scaleFactor,samprate,status);
				if (status == 1)
				{
					writeLog(fileRef,"Warning", verbosity,"lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.");
					status = 0;
				}
				if (status == 2)
				{
					writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
					return(status);
				}
				if (!status == EPOK) return (status);

				// PULSE MODEL: Derivative after filtering
				modelSGN = gsl_vector_alloc(model->size);
				status = derMTHSimple (&model,&modelSGN,model->size,status);
				if (!status == EPOK) return (status);

				gsl_matrix_set_row(models,i,model);
			}
			gsl_vector_free(model);
		}

			// Called iteration function
		//MC status = DALiterateData(n_cols, cols, offset, rows_per_loop, inDataIterator, 0L, status);
		if(fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIterator,0L,&status))printerror(status);

		fclose(temporalFile);

		/////////////////////// Write output keywords
		//MC status = DALattributePutInt(trgExten, "EVENTCNT", totalpulses-1, unit, comment,	status);
		strcpy(extname,"EUR-TRG");
		if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
		ttpls1 = totalpulses-1;
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(trgObject,TINT,keyname,&ttpls1,comment,&status)) printerror(status);


		if ((energy > 0) && (totalpulses == 0))
		{
			writeLog(fileRef, "Warning", verbosity, "Some pulses have been not found");
		}
	}

	/////////////////////// Close input FITS files
	//MC status = DALobjectClose(inObject, DAL_SAVE, status);

	/////////////////////// Close output FITS files
	//MC status = DALobjectClose(trgObject, DAL_SAVE, status);

	// Close pulses models library file (EUR-LIB extension)
	//MC status = DALobjectClose(inLibObject, DAL_SAVE, status);
	if(fits_close_file(inObject,&status) || fits_close_file(trgObject,&status) || 
	  fits_close_file(inLibObject,&status))  printerror(status);
	

	//gsl_vector_free(pulsegsl);

	if (mode != 0)
	{
		//gsl_vector_free(modelSGN);

		gsl_matrix_free(library);
		gsl_matrix_free(models);
	}
	fclose(temporalFile2);

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
		EPexit(-8000);
	}
	
	if (energy != 0)
	{
		/////////////////////// Free memory
		//delete [] straux;
		delete [] obj.nameTable;
		delete [] obj.nameCol;
		delete [] obj.unit;
	}

	status=fclose(fileRef);
	if(!status == EPOK) EPexit(-4003);
	EPexit(status);

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
int initModule(int argc, char **argv, int status)
{
	if (!status == EPOK)	return status;
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
    triggerPars[3].description = "Strict calibration mode (0), calibration (1) or production mode (2)";
    triggerPars[3].defValInt = 2;
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
							if ((!isdigit(optarg[0]) && (optarg[0] != '-'))) return(-3001);
							if ((!isdigit(optarg[0]) && (optarg[0] == '-')))
							{
								if (!isdigit(optarg[1])) return(-3001);
							}

							triggerPars[i].ValInt = atoi(optarg);
						}
						else
						{
							if ((!isdigit(optarg[0]) && (optarg[0] != '-')))	return(-3001);
							if ((!isdigit(optarg[0]) && (optarg[0] == '-')))
							{
								if (!isdigit(optarg[1]))	return(-3001);
							}

							triggerPars[i].ValReal= atof(optarg);
						}
						break;
					} // endif
				} // endfor
				break;
		    default:
		    	cout << "Error: invalid parameter name "<< long_options[optidx].name << endl;
		    	EPexit(-3007);
		}//switch
	}//while
	// If command line is empty: ask for params interactively
	if(commandLine == 0)
	{
		status=interactivePars(triggerPars,npars,task);
		if(status) EPexit(status);
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
	    fst = access(inName, R_OK);
	    if (fst < 0) 						return(-62102); 						
	  }else if(triggerPars[i].name == "inLibFile"){
	    strcpy(inLibName, triggerPars[i].ValStr.c_str()); 
	    fst = access(inLibName, R_OK);
	    fstLib=fst;
	    //if (fst < 0) 						return(-62103);
	  }else if(triggerPars[i].name == "outFile"){
	    strcpy(trgName, triggerPars[i].ValStr.c_str()); 
	  }else if(triggerPars[i].name == "mode"){
	    mode = triggerPars[i].ValInt;
	    if ((mode != 0) && (mode != 1) && (mode != 2))	return (-62620);
	    if ((mode == 1) || (mode == 2))
	    {
	    	if (fstLib < 0) 						return(-62103);
	    }
	  }else if(triggerPars[i].name == "b_cF"){
	    b_cF = triggerPars[i].ValReal; 
	    if (b_cF < 0)						return (-62621);
	  }else if(triggerPars[i].name == "c_cF"){
	    c_cF = triggerPars[i].ValReal; 
	    if (c_cF < 0)						return (-62622);
	  }else if(triggerPars[i].name == "LrsT"){
	    LrsT = triggerPars[i].ValReal; 
	    if (LrsT <= 0)						return (-62623);
	  }else if(triggerPars[i].name == "LbT"){
	    LbT = triggerPars[i].ValReal; 
	    if (LbT <= 0)						return (-62624);
	  }else if(triggerPars[i].name == "tauFALL"){
	    tauFALL = triggerPars[i].ValReal; 
	    if (tauFALL <= 0)						return (-62612);
	  }else if(triggerPars[i].name == "scaleFactor"){
	    scaleFactor = triggerPars[i].ValReal; 
	    if (scaleFactor <= 0)					return (-62619);
	  }else if(triggerPars[i].name == "samplesUp"){
	    samplesUp = triggerPars[i].ValInt; 
	    if (samplesUp <= 0)						return (-62615);
	  }else if(triggerPars[i].name == "nSgms"){
	    nSgms = triggerPars[i].ValInt; 
	    if (nSgms <= 0)						return (-62614);
	  }else if(triggerPars[i].name == "ntaus"){
	    ntaus = triggerPars[i].ValInt; 
	    if (ntaus <= 0)						return (-62613);
	  }else if(triggerPars[i].name == "tAftrtstart"){
	    tAftrtstart = triggerPars[i].ValReal; 
	    if (tAftrtstart < 0)					return (-62616);
	  }else if(triggerPars[i].name == "numBitsQuality"){
	    numBitsQual = triggerPars[i].ValInt; 
	    if (numBitsQual< 4 || numBitsQual > 32)			return (-62609);
	  }else if(triggerPars[i].name == "writePulse"){
	    writePulse = triggerPars[i].ValInt; 
	  }else if(triggerPars[i].name == "getTaus"){
	    getTaus = triggerPars[i].ValInt;
	  }else if(triggerPars[i].name == "nameLog"){
	    strcpy(nameLog,triggerPars[i].ValStr.c_str()); 
	  }else if(triggerPars[i].name == "verbosity"){
	    verbosity = triggerPars[i].ValInt; 
	    if (verbosity<0 || verbosity>3)				return (-62606);
	  }
	}
	// END OF SECTION 2/2 TO BE MODIFIED IF PARAMETERS NEED TO BE MODIFIED

	fileRef = fopen(nameLog,"w+");	// Remove file if it already exists and open a new file to save log messages
	if (fileRef == NULL) return(-4001);
	
	return (status);
}


/*xxxx end of SECTION 3 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 4 ************************************************************
 * seconds2Bins function: This function transforms the units of the input parameters from seconds to bins.
 *
 ****************************************************************************/
int seconds2Bins (int status)
{
	if (status != EPOK) return status;

	sizePulse_b = (int)(sizePulse * samprate);

	return status;
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
 * createTriggerFile function: This function creates the structure of the output FITS file TRIGGER
 *
 * - Create output FITS file: If it does not exist already
 * - Create extension EUR-TRG
 * - Create keywords
 *
 ***************************************************************************/
int createTriggerFile(int status)
{
	string message;
	int extver=0;
	
	//MC if (!status == 0)		return status;

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
	if (status == EPOK) {// if file already exists => error 
		//MC status = DALobjectOpen(trgName, &trgObject, status);
		fits_open_file(&trgObject, trgName,0,&status);
		if (status == EPOK) {status = -62101; return status;}
		else {status = EPOK;}
	}

	//MC status = DALobjectCreate(trgName, DAL_DISK, NULL, &trgObject, status);
	if(fits_create_file(&trgObject, trgName, &status)) printerror(status);

	sprintf(straux, "%d", status);
	message = "Create Trigger Fits File: " + string(trgName) + string(straux);
	writeLog(fileRef,"Log", verbosity,message);

	/////////////////////// Create extension EUR-TRG
	//MC status = DALtableCreate(trgObject, DAL_DISK, 0, "EUR-TRG", &trgExten, status);
	char *tt[1];
	char *tf[1];
	char *tu[1];
	strcpy(extname,"EUR-TRG");
	if(fits_open_file(&trgObject,trgName,READWRITE,&status))printerror(status);
	if(fits_create_tbl(trgObject, BINARY_TBL,0,0,tt,tf,tu,extname,&status))printerror(status);

	/////////////////////// Create extension EUR-TEST and its columns
	//MC status = DALtableCreate(trgObject, DAL_DISK, 0, "EUR-TEST", &trgExten1, status);
	strcpy(extname,"EUR-TEST");
	if(fits_create_tbl(trgObject,BINARY_TBL,0,0,tt,tf,tu,extname,&status))printerror(status);

	/////////////////////// Create keywords
	strcpy(extname,"EUR-TRG");
	if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
	
	//MC status = DALattributePutInt(trgExten, "EVENTSZ", sizePulse_b, unit, comment, status);
	strcpy(keyname,"EVENTSZ");
	if(fits_write_key(trgObject,TLONG,keyname,&sizePulse_b,comment,&status)) printerror(status);
	//MC status = DALattributePutReal(trgExten, "ENERGY", energy, unit, comment, status);
	strcpy(keyname,"ENERGY");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&energy,comment,&status)) printerror(status);
	/*status = writeKeywords (kg, kgNumber, trgName, trgExten, fileRef,  status);
	status = writeKeywords (ki, kiNumber, trgName, trgExten, fileRef,  status);*/
	//MC status = DALattributePutReal(trgExten, "SAMPRATE", samprate, unit, comment, status);
	strcpy(keyname,"SAMPRATE");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&samprate,comment,&status)) printerror(status);

	//MC status = DALattributePutInt(trgExten, "NBQUAL", numBitsQual, unit, comment, status);
	strcpy(keyname,"NBQUAL");
	if(fits_write_key(trgObject,TINT,keyname,&numBitsQual,comment,&status)) printerror(status);
	//MC status = DALattributePutInt(trgExten, "TRG_ID", 1, unit, comment, status);
	strcpy(keyname,"TRG_ID");
	if(fits_write_key(trgObject,TINT,keyname,&trg_id,comment,&status)) printerror(status);
	//MC status = DALattributePutChar(trgExten, "CREATOR", create, unit, comment, status);
	strcpy(keyname,"CREATOR");
	strcpy(keyvalstr,create);
	if(fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);

	// Set process keyword
	char str_verb[125];			sprintf(str_verb,"%d",verbosity);
	char str_numB[125];			sprintf(str_numB,"%d",numBitsQual);
	char str_ntaus[125];		sprintf(str_ntaus,"%d",ntaus);
	char str_tauFall[125];		sprintf(str_tauFall,"%e",tauFALL);
	char str_scaleFactor[125];	sprintf(str_scaleFactor,"%f",scaleFactor);
	char str_samplesUp[125];	sprintf(str_samplesUp,"%d",samplesUp);
	char str_nSgms[125];	    sprintf(str_nSgms,"%f",nSgms);
	char str_wp[125];			if (writePulse) sprintf(str_wp,"y"); else  sprintf(str_wp,"n");
	char str_getTaus[125];		if (getTaus) sprintf(str_getTaus,"y"); else  sprintf(str_getTaus,"n");
	char str_mode[125];			sprintf(str_mode,"%d",mode);
	char str_b_cF[125];			sprintf(str_b_cF,"%f",b_cF);
	char str_c_cF[125];			sprintf(str_c_cF,"%f",c_cF);
	char str_LrsT[125];			sprintf(str_LrsT,"%e",LrsT);
	char str_LbT[125];			sprintf(str_LbT,"%e",LbT);
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

	//MC status = DALattributePutChar(trgExten, "PROCESS",process.c_str(), unit, comment, status);
	strcpy(keyname,"PROCESS");
	strcpy(keyvalstr,process.c_str());
	if(fits_write_key_longwarn(trgObject,&status))printerror(status);
	if(fits_write_key_longstr(trgObject,keyname,keyvalstr,comment,&status)) printerror(status);
	//MC status = DALattributePutChar(trgExten, "FTYPE", "pulse", unit, comment, status);
	strcpy(keyname,"FTYPE");
	strcpy(keyvalstr,"pulse");
	if(fits_write_key(trgObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);

	//MC status = DALattributePutInt(trgExten, "NTAUS", ntaus, unit, comment, status);
	strcpy(keyname,"NTAUS");
	if(fits_write_key(trgObject,TINT,keyname,&ntaus,comment,&status)) printerror(status);
	//MC status = DALattributePutReal(trgExten, "TAUFALL", tauFALL, unit, comment, status);
	strcpy(keyname,"TAUFALL");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&tauFALL,comment,&status)) printerror(status);
	//MC status = DALattributePutReal(trgExten, "SCLFCTR", scaleFactor, unit, comment, status);
	strcpy(keyname,"SCLFCTR");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&scaleFactor,comment,&status)) printerror(status);
	
	//MC if ((mode == 0) || (mode == 1))
	//MC {
	//MC 	status = DALattributePutInt(trgExten, "MODE", 0, unit, comment, status);
	//MC }
	//MC else
	//MC {
	//MC	status = DALattributePutInt(trgExten, "MODE", 1, unit, comment, status);
	//MC }
	modeval = 1;
	if ((mode == 0) || (mode == 1)) modeval = 0;
	strcpy(keyname,"MODE");
	if(fits_write_key(trgObject,TINT,keyname,&modeval,comment,&status)) printerror(status);	
	
	//MC status = DALattributePutReal(trgExten, "B_CF", b_cF, unit, comment, status);
	strcpy(keyname,"B_CF");
	if(fits_write_key(trgObject,TDOUBLE,keyname,&b_cF,comment,&status)) printerror(status);
	//MC status = DALattributePutReal(trgExten, "C_CF", c_cF, unit, comment, status);
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
int inDataIterator(long totalrows, long offset, long firstrow, long nrows,	int ncols, iteratorCol *cols, void *user_strct)
{
	string message;
	//MC 
	int status=EPOK;
	int extver=0;
	
	//MC if (!status == 0)	return status;

	/////////////////////// Declare variables
	static double *time;	// Vector of TIME column
	static double *v;		// Vector of V column
	nPulsesRow = 0;			// Number of pulses in each row

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
	//MC status = DALiterGetArray(&cols[0], ((void **)&time), status);
	time = (double *) fits_iter_get_array(&cols[0]);
	//MC status = toGslVector(((void **)&time), &timegsl, nrows, 0, DAL_DOUBLE, status);
	status = toGslVector(((void **)&time), &timegsl, nrows, 0, TDOUBLE, status);
	//MC status = DALiterGetArray(&cols[1], ((void **)&v), status);
	v = (double *) fits_iter_get_array(&cols[1]);
	//MC status = toGslMatrix(((void **)&v), &vgsl, eventsz, nrows, (int)DAL_DOUBLE,	0, status);
	status = toGslMatrix(((void **)&v), &vgsl, eventsz, nrows, (int)TDOUBLE, 0, status);

	char val[256];
	char val_aux[256];
	//tstartevents = gsl_vector_alloc(eventcnt);
	//event = gsl_vector_alloc(eventsz);

	
	///////////////////////  Processing each row of pulses
	for (int i=0; i< nrows; i++) {

		sprintf(straux,"%d",ntotalrows);
		message = "-------------> Row: " + string(straux);
		sprintf(straux,"%d",eventcnt);		
		message += " of " + string(straux) + " <------------------ ";
		writeLog(fileRef,"Log", verbosity,message);
		sprintf(val,"-----------> Row: ");
		sprintf(val_aux,"%d",ntotalrows);
		strcat(val,val_aux);
		sprintf(val_aux," of ");
		strcat(val,val_aux);
		sprintf(val_aux,"%d",eventcnt);
		strcat(val,val_aux);
		strcat(val,"\n");
		fputs(val,temporalFile);

		double time0 = gsl_vector_get(timegsl, i);	// Time in the begin of the row
		gsl_matrix_get_row(invector, vgsl, i);
		gsl_vector_scale(invector,ivcal);		//IVCAL to change arbitrary units of voltage to non-arbitrary units of current (Amps)

		if (i == 0)		timefirstevent = time0;
		//gsl_vector_set(tstartevents,ntotalrows-1,gsl_vector_get(timegsl,i));

		// Assign positive polarity to the pulses
		if (((asquid>0) && (plspolar<0)) || ((asquid<0) && (plspolar>0)))
		{
			gsl_vector_scale(invector,-1);
			chngplrt = 1;
		}
		//MC status = DALattributePutInt(trgExten, "CHNGPLRT", chngplrt, unit, comment, status);
		strcpy(extname,"EUR-TRG");
		if(fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status)) printerror(status);
		strcpy(keyname,"CHNGPLRT");
		if(fits_update_key(trgObject,TINT,keyname,&chngplrt,comment,&status))printerror(status);

		gsl_vector_memcpy(invectorNOTFILTERED,invector);

		//gsl_vector_memcpy(event, invectorNOTFILTERED);

		// Low-pass filtering
		status = lpf_boxcar(&invector,invector->size,scaleFactor*tauFALL,samprate,status);
		if (status == 1)
		{
	   	   writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.");
	   	   status = 0;
		}
		if (status == 2)
		{
	  	  writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
	  	  return(status);
		}
		if (!status == EPOK) return (status);
		gsl_vector_memcpy(invectorFILTERED,invector);

		/*if (indice == 26)
		{
			writeLog(fileRef,"Warning", verbosity,status,"Va a escribir el fichero");
			FILE * temporalFile;
			char val[256];
			char val_aux[256];
			char temporalFileName[255];
			sprintf(temporalFileName,trgName);
			strcat(temporalFileName,".txt");
			//temporalFile = fopen ("temporal.txt","w");
			temporalFile = fopen (temporalFileName,"w");
			//cout<<"temporalFileName: "<<temporalFileName<<endl;
			for (int i=0;i<eventsz;i++)
			{
				sprintf(val,"%e",gsl_vector_get(invectorNOTFILTERED,i));
				strcat(val," ");
				sprintf(val_aux,"%e",gsl_vector_get(invectorFILTERED,i));
				strcat(val,val_aux);
				//strcat(val," ");
				//sprintf(val_aux,"%e",thresholdmediankappa);
				//strcat(val,val_aux);
				strcat(val,"\n");
				fputs (val,temporalFile);
			}
			fclose(temporalFile);
		}*/

		// Derivative after filtering
		status = derMTHSimple (&invector,&SGN,invector->size,status);
		if (!status == EPOK) return (status);
		gsl_vector_memcpy(invectorDERIVATIVE,invector);
		gsl_vector_memcpy(derSGN,SGN);

		
		if (indice+1 == 86)
		{
			//  Creating NTFltRow Column
			if (status == EPOK) {
				//MC obj.inObject = trgExten1;
				obj.inObject = trgObject;
				obj.nameTable = new char [255];
				strcpy(obj.nameTable,"EUR-TEST");
				obj.iniRow = 1;
				obj.endRow = eventsz;
				obj.iniCol = 0;
				obj.nameCol = new char [255];
				strcpy(obj.nameCol,"NTFltRow");
				//MC obj.type = DAL_DOUBLE;
				obj.type = TDOUBLE;
				obj.unit = new char [255];
				strcpy(obj.unit," ");
			}
			status = writeFitsSimple(obj, invectorNOTFILTERED, status);

			//  Creating FltRow Column
			if (status == EPOK)
			{
				strcpy(obj.nameCol,"FltRow");
			}
			status = writeFitsSimple(obj, invectorFILTERED, status);

			//  Creating DerRow Column
			if (status == EPOK) {
				strcpy(obj.nameCol,"DerRow");
			}
			status = writeFitsSimple(obj, invectorDERIVATIVE, status);
		}
		
		/*sprintf(val,"vectorinNOTFIL vectorinFIL vectorinDER");
		strcat(val,"\n");
		fputs (val,temporalFile);
		for (int i=0; i<invectorFILTERED->size;i++)
		{
			sprintf(val,"%e %e %e",gsl_vector_get(invectorNOTFILTERED,i),gsl_vector_get(invectorFILTERED,i),gsl_vector_get(invectorDERIVATIVE,i));
			strcat(val,"\n");
			fputs (val,temporalFile);
		}*/

		// Rows (events) are processed in only one step (without cutting into several segments)
		initialtime = time0;

		status = procRow(invectorNOTFILTERED, invectorDERIVATIVE, &nPulsesRow, status);

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
int procRow(gsl_vector *vectorNOTFIL, gsl_vector *vectorDER, int *npulses, int status)
{
	char val[256];
	char val_aux[256];

	if (!status == 0)	return status;

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
	if (mode == 2)
	{
		status = findPulses (vectorNOTFIL, vectorDER, &tstartgsl, &qualitygsl, &energygsl,
				npulses,
				1,
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
				status,
				indice+1);
	}
	else
	{
		status = findPulses (vectorNOTFIL, vectorDER, &tstartgsl, &qualitygsl, &energygsl,
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
				status,
				indice+1);
	}
	sprintf(val,"%d %d",indice+1,nPulsesRow);
	strcat(val,"\n");
	fputs(val,temporalFile2);
	/*sprintf(val,"vectorDER");
		strcat(val,"\n");
		fputs(val,temporalFile);
	if (indice+1 == 207)
	{
		for (int j=0;j<vectorDER->size;j++)
		{
			sprintf(val,"%.12f",gsl_vector_get(vectorDER,j));
			strcat(val,"\n");
			fputs(val,temporalFile);
		}
	}*/

	// Calculate the tend of each found pulse
	/*sprintf(val,"nPulsesRowTOT: ");
	sprintf(val_aux,"%d",nPulsesRow);
	strcat(val,val_aux);
	strcat(val,"\n");
	fputs(val, temporalFile);*/

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

	/*for (int i=0;i<nPulsesRow;i++)
		{
			sprintf(val,"Pulsesfound( ");
			sprintf(val_aux,"%d",i);
			strcat(val,val_aux);
			sprintf(val_aux,"): ");
			strcat(val,val_aux);
			sprintf(val_aux,"%e",gsl_vector_get(tstartgsl,i));
			strcat(val,val_aux);
			strcat(val,"\n");
			fputs(val, temporalFile);
		}*/

	/*/////////////////////// Obtain in a precise way the tstart of the founded pulses
		// Declare variables
	gsl_vector *pulsevector;		// Current pulse
	gsl_vector *pulsevectorFIL;		// Low-pass filtered current pulse
	long pulsevectorFILmax_index;	// Index where the maximum of the low-pass filtered current pulse is
	gsl_vector *SGN1;
	long sizePulse_bi;				// Current pulse size
	                                // (lower than sizePulse_b if the following pulse is closer than sizePulse_b)
	long tstart0;					// From the beginning of 'pulsevector' => Add to tstartgsl (which takes into account the safety margin)
	double tstartj;                 // From the beginning of 'pulsevector' => Add to tstartgsl (which takes into account the safety margin)
	double sF = scaleFactor/10;		// The first scale factor to use is 10 times lower the input parameter 'scaleFactor'
									// HARPOINT!!!
	double sFj;						// 0.1*scaleFactor, 0.2*scaleFactor, 0.3*scaleFactor...
	bool iterate = true;
	gsl_vector *tstartjgsl=gsl_vector_alloc(10);	// It stores the tstart's calculated by using different cut
													// frequencies (different sFj) of the low-pass filter

	int tstart0LIMITED = 0;			// Not limited
	int tstartjLIMITED = 0;

		// Auxiliary variables
	gsl_vector *pulsevector_aux;
	double t0;	//tstartgsl_i (only to shorten the code)
	gsl_vector_view temp;	// In order to handle with gsl_vector_view (subvectors)

	for (int i=0; i<nPulsesRow; i++)
	// For each pulse...
	{
		iterate = true;
		sizePulse_bi = gsl_vector_get(tendAUXgsl,i)-gsl_vector_get(tstartgsl,i)+1;
		pulsevector = gsl_vector_alloc(sizePulse_bi);
		pulsevector_aux = gsl_vector_alloc(sizePulse_bi);
		pulsevectorFIL = gsl_vector_alloc(sizePulse_bi);
		SGN1 = gsl_vector_alloc(sizePulse_bi);

		t0 = gsl_vector_get (tstartgsl,i);
		sprintf(val,"tstart: %.12e",t0);
		strcat(val,"\n");
		fputs(val,temporalFile);

		// Extract the current pulse from the event
		if (vectorNOTFIL->size - t0 > sizePulse_bi)
		// Not truncated pulse
		{
			temp = gsl_vector_subvector(vectorNOTFIL,t0+1,sizePulse_bi);

			gsl_vector_memcpy(pulsevector,&temp.vector);
		}
		else
		// Truncated pulse
		{
			for (int j=0; j<(vectorNOTFIL->size) - t0; j++)
			{
				if (t0 == -1) t0 = 0;
				gsl_vector_set(pulsevector,j,gsl_vector_get(vectorNOTFIL,j+t0));
			}
			for (int j=(vectorNOTFIL->size)-t0; j< sizePulse_bi; j++) // Fill in with 0's
			{
				gsl_vector_set(pulsevector,j,0.0);
			}
		}

		sprintf(val,"gsl_vector_max_index(pulsevector): %d",gsl_vector_max_index(pulsevector));
		strcat(val,"\n");
		fputs(val,temporalFile);
		if (gsl_vector_max_index(pulsevector)<10)					// HARDPOINT!!!!!
		//if (gsl_vector_max_index(pulsevector)<3)					// HARDPOINT!!!!!
		{
			iterate = false;

		}

		// Low-pass filtering by using sF  = scaleFactor/10
		gsl_vector_memcpy(pulsevector_aux,pulsevector);
		status = lpf_boxcar(&pulsevector_aux,pulsevector_aux->size,sF*tauFALL,samprate,status);
		if (status == 1)
		{
			writeLog(fileRef,"Warning", verbosity,"lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.");
			//sprintf(val,"1st sF lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.");
			//strcat(val,"\n");
			//fputs(val,temporalFile);
			status = 0;
		}
		if (status == 2)
		{
			writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
			//sprintf(val,"1st sF lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
			//strcat(val,"\n");
			//fputs(val,temporalFile);
			iterate = false;
			status = 0;
		}

		if (iterate == true)
		{
			gsl_vector_memcpy(pulsevectorFIL,pulsevector_aux);
			pulsevectorFILmax_index = gsl_vector_max_index(pulsevectorFIL);

			// Calculate tstart0 by calculating the maximum of the derivative of the filtered data
			status = derMTHSimple (&pulsevector_aux,&SGN1,pulsevector_aux->size,status);
			if (!status == EPOK) return (status);
			tstart0 = gsl_vector_max_index(pulsevector_aux);
			// The first and simplest correction of the tstart
			if (tstart0 > pulsevectorFILmax_index)
			{
				tstart0 = pulsevectorFILmax_index;
				tstart0LIMITED = 1;
			}
			gsl_vector_set(tstartjgsl,0,tstart0);

			sFj = sF+sF;
			for (int j=1;j<10;j++)
			// For different cut frequencies of the low-pass filter...
			// sFj = 0.2scaleFactor, 0.3scaleFactor,...
			// sFj goes increasing => Cut frequency decreases => boxLength increases
			{
				tstartjLIMITED = 0;

				// Low-pass filtering by using sF
				gsl_vector_memcpy(pulsevector_aux,pulsevector);
				status = lpf_boxcar(&pulsevector_aux,pulsevector_aux->size,sFj*tauFALL,samprate,status);
				if (status == 1)
				{
					writeLog(fileRef,"Warning", verbosity,"lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.");
					sprintf(val,"lpf_boxcar(Model): tauFALL*scaleFactor too small => Cut-off frequency too high => Equivalent to not filter.");
					strcat(val,"\n");
					fputs(val,temporalFile);
					status = 0;
				}
				if (status == 2)
				// At the end => break, because boxLength would go on increasing
				{
					writeLog(fileRef,"Warning", verbosity,"lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
					sprintf(val,"lpf_boxcar: tauFALL*scaleFactor too high => Cut-off frequency too low");
					strcat(val,"\n");
					fputs(val,temporalFile);

					if (j == 1)
					{
						// tstart0/2 (averaging) + tstartgsl
						gsl_vector_set(tstartgsl,i,floor((tstart0+gsl_vector_get(tstartgsl,i)+gsl_vector_get(tstartgsl,i))/2));
						status = 0;
						break;
					}
					else
					{
						tstart0 = gsl_vector_get(tstartjgsl,0);
						for (int k=1;k<j;k++)
						{
							tstart0 = tstart0+gsl_vector_get(tstartjgsl,k);
						}
						// (tstart0+tstartjgsl_1+tstartjgsl_2 +...+tstartjgsl_j)/j + tstartgsl
						//gsl_vector_set(tstartgsl,i,+gsl_vector_get(tstartgsl,i)+tstart0/(tstartjgsl->size));
						gsl_vector_set(tstartgsl,i,gsl_vector_get(tstartgsl,i)+tstart0/j);
						status = 0;
						break;
					}
				}
				gsl_vector_memcpy(pulsevectorFIL,pulsevector_aux);
				pulsevectorFILmax_index = gsl_vector_max_index(pulsevectorFIL);

				// Calculate tstartj by calculating the maximum of the derivative of the filtered data
				status = derMTHSimple (&pulsevector_aux,&SGN1,pulsevector_aux->size,status);
				tstartj = gsl_vector_max_index(pulsevector_aux);
				if (tstartj > pulsevectorFILmax_index)
				{
					tstartj = pulsevectorFILmax_index;
					tstartjLIMITED = 1;
				}
				gsl_vector_set(tstartjgsl,j,tstartj);

				// Compare tstart0 and tstartj
				// if tstart0 and tstartj are similar => Stop ('break') and tstart=(tstart0+tstartj)/2
				// else => Another iteration sFj=sFj+sF*j
				if ((fabs(tstartj-tstart0) <= 2) && ((tstart0LIMITED != 1) || (tstartjLIMITED != 1)) && (tstart0-tstartj != 0))
				//if ((fabs(tstartj-tstart0) <= 0.1*safetyMarginTstart) && ((tstart0LIMITED != 1) || (tstartjLIMITED != 1)) && (tstart0-tstartj != 0))
				//if ((fabs(tstartj-tstart0) <= 0.1*safetyMarginTstart) && ((tstart0LIMITED != 1) || (tstartjLIMITED != 1)))
				//if (fabs(tstartj-tstart0) <= 0.1*safetyMarginTstart)
				//if (fabs(tstartj-tstart0) <= 0.2*safetyMarginTstart)
				{
					sprintf(val,"Entra1");
					strcat(val,"\n");
					fputs(val,temporalFile);
					if ((i != nPulsesRow-1) && (floor((tstart0+tstartj)/2)+gsl_vector_get(tstartgsl,i) < gsl_vector_get(tstartgsl,i+1))
					|| ((i != 0) && (floor((tstart0+tstartj)/2)+gsl_vector_get(tstartgsl,i) > gsl_vector_get(tstartgsl,i-1))) || ((i==0) && (i==nPulsesRow-1)))
					{
						sprintf(val,"Entra1_1");
						strcat(val,"\n");
						fputs(val,temporalFile);
						gsl_vector_set(tstartgsl,i,floor((tstart0+tstartj)/2)+gsl_vector_get(tstartgsl,i));

						break;
					}
				}
				else
				{
					sprintf(val,"Entra2");
					strcat(val,"\n");
					fputs(val,temporalFile);
					tstart0 = tstartj;
					tstart0LIMITED = tstartjLIMITED;
					sFj = sFj + sF;
					if (j == 10-1)
					{
						//tstart0 = gsl_vector_get(tstartjgsl,0);
						tstart0 = gsl_vector_get(tstartjgsl,5);
						//for (int k=1;k<tstartjgsl->size;k++)
						for (int k=6;k<tstartjgsl->size;k++)
						{
							tstart0 = tstart0+gsl_vector_get(tstartjgsl,k);
						}
						//gsl_vector_set(tstartgsl,i,gsl_vector_get(tstartgsl,i)+tstart0/(tstartjgsl->size));
						if ((i != nPulsesRow-1) && (gsl_vector_get(tstartgsl,i)+tstart0/5 < gsl_vector_get(tstartgsl,i+1))
						|| ((i != 0) && (gsl_vector_get(tstartgsl,i)+tstart0/5 > gsl_vector_get(tstartgsl,i-1))) || ((i==0) && (i==nPulsesRow-1)))
						{
							sprintf(val,"Entra2_1");
							strcat(val,"\n");
							fputs(val,temporalFile);
							if (i!=nPulsesRow-1)
							{
								sprintf(val,"Siguiente: %.12e",gsl_vector_get(tstartgsl,i+1));
								strcat(val,"\n");
								fputs(val,temporalFile);
							}
							if (i!= 0)
							{
								sprintf(val,"Anterior: %.12e",gsl_vector_get(tstartgsl,i-1));
								strcat(val,"\n");
								fputs(val,temporalFile);
							}
							gsl_vector_set(tstartgsl,i,gsl_vector_get(tstartgsl,i)+tstart0/5);
							sprintf(val,"El del pulso i: %.12e",gsl_vector_get(tstartgsl,i));
							strcat(val,"\n");
							fputs(val,temporalFile);
						}
					}
				}
			}

			gsl_vector_free(pulsevector);
			gsl_vector_free(pulsevector_aux);
			gsl_vector_free(pulsevectorFIL);
			gsl_vector_free(SGN1);
		}
		else
		{

			gsl_vector_set(tstartgsl,i,gsl_vector_get(tstartgsl,i)+safetyMarginTstart);
		}
	}
	gsl_vector_free(tstartjgsl);*/

	for (int i=0;i<nPulsesRow;i++)
	{
		gsl_vector_set(tstartgsl,i,gsl_vector_get(tstartgsl,i)+safetyMarginTstart);
	}
	//if (nPulsesRow!=0) gsl_vector_free(tendAUXgsl);

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
	//...............................................

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
		status = obtainTau (vectorNOTFIL, tstartWOutSftMrggsl, tendWOutSftMrggsl, nPulsesRow, &tauRisegsl, &tauFallgsl, status);
	}

	/////////////////////// Write pulses and measurements in output FITS files
	if (nPulsesRow != 0)	pulsesgsl = gsl_matrix_alloc(nPulsesRow,sizePulse_b);

	status = writePulses (vectorNOTFIL, vectorDER, tstartgsl, tendgsl, qualitygsl, tauRisegsl, tauFallgsl, energygsl, &pulsesgsl, status);


	/*if ((mode == 0) && (indice == 0))
	{
		status = writeLibrary(gsl_vector_get(energygsl,0), &pulsegsl, status);
	}*/
	if (mode == 0)
	{
		if (indice == 0)
		{
			pulseheighttemplate = 0.0;
			pulseaverage = gsl_vector_alloc(sizePulse_b);
			gsl_vector_set_all(pulseaverage,0.0);
			row_aux = gsl_vector_alloc(sizePulse_b);
		}

		for (int i=0;i<nPulsesRow;i++)
		{
			pulseheighttemplate = pulseheighttemplate + gsl_vector_get(energygsl,i);
			gsl_matrix_get_row(row_aux,pulsesgsl,i);
			gsl_vector_add(pulseaverage,row_aux);
		}

		if (indice == eventcnt-1)
		{
			gsl_vector_scale(pulseaverage,1.0/((double) totalpulses));
			pulseheighttemplate = pulseheighttemplate/totalpulses;

			//status = getB(pulseaverage, tstartgsl, nPulsesRow, &Lbgsl, sizePulse_b, &Bgsl, temporalFile,status);
			//status = getPulseHeight(pulseaverage, gsl_vector_get(tstartgsl,0), 0, 1, Lrs, gsl_vector_get(Lbgsl,0), gsl_vector_get(Bgsl,0), sizePulse_b, &pulseheighttemplate, temporalFile, status);
			status = writeLibrary(pulseheighttemplate, &pulseaverage, status);
			gsl_vector_free (pulseaverage);
			gsl_vector_free (row_aux);
		}

	}

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
int obtainTau (gsl_vector *invector, gsl_vector *tstartgsl, gsl_vector *tendgsl, int nPulses, gsl_vector **taurisegsl, gsl_vector **taufallgsl, int status)
{
	if (!status == EPOK || nPulses == 0)	return status;

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
int writePulses(gsl_vector *invectorNOTFIL, gsl_vector *invectorDER, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, gsl_vector *energy, gsl_matrix **pulses, int status)
{

	if (!status == EPOK)	return status;
	char val[256];

	/////////////////////// TRG Extension
	int t0;   			// First value of index of pulse
	gsl_matrix *vgslout2;
	gsl_vector *difTstart; // Difference between the tstart of two consecutive pulses
	                       // First pulse row difTtsart will be -1
	gsl_vector *tailBeforegsl;
	bool tailBefore = true;
	int min = 5;
	double tmp0;

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

		/*if (indice == 40)	//Uno menos del evento que se quiere manejar
		{
			char inNameaux[255];
			char * pch;
			sprintf(inNameaux,inName);
			if (strlen(inNameaux)>=6)
			{
				//Check if inName has '.fits' and if yes, delete it
				if (strncmp(strndup(inNameaux+strlen(inNameaux)-5, 5),".fits",5) == 0)
				{
					// inName finishes as '.fits' => Delete '.fits' to inName => inNameaux
					pch = strtok (inNameaux,".");
					sprintf(inNameaux,pch);
				}
			}
			FILE * generalFile;
			FILE * eventFile;
			FILE * tstartFile;
			FILE * pulseFile;
			char val[256];
			char filegeneral[255];
			char fileevent[255];
			char filetstart[255];
			char filepulse[255];

			sprintf(filegeneral,inNameaux);
			sprintf(val,"_event%d_general.txt",indice+1);
			strcat(filegeneral,val);
			generalFile = fopen (filegeneral,"w");
			//sprintf(val,"%.12e",timefirstevent);
			//strcat(val,"\n");
			//fputs (val,generalFile);
			sprintf(val,"%.14e",gsl_vector_get(tstartevents,indice));
			strcat(val,"\n");
			fputs (val,generalFile);
			sprintf(val,"%d",nPulsesRow);
			strcat(val,"\n");
			fputs (val,generalFile);
			sprintf(val,"%e",samprate);
			strcat(val,"\n");
			fputs (val,generalFile);
			if (chngplrt == 1) biasvolt = biasvolt*(-1.0);
			sprintf(val,"%.14e",biasvolt);
			strcat(val,"\n");
			fputs (val,generalFile);

			sprintf(fileevent,inNameaux);
			sprintf(val,"_event%d.txt",indice+1);
			strcat(fileevent,val);
			eventFile = fopen (fileevent,"w");
			for (int i=0; i<eventsz; i++)
			{
				sprintf(val,"%.14e",gsl_vector_get(event,i));
				strcat(val,"\n");
				fputs (val,eventFile);
			}

			sprintf(filetstart,inNameaux);
			sprintf(val,"_event%d_tstartpulses.txt",indice+1);
			strcat(filetstart,val);
			tstartFile = fopen (filetstart,"w");

			for (int i=0; i<nPulsesRow; i++)
			{
				sprintf(val,"%.14e",gsl_vector_get(tstart,i));
				strcat(val,"\n");
				fputs (val,tstartFile);

				sprintf(filepulse,inNameaux);
				sprintf(val,"_event%d_pulse%d.txt",indice+1,i+1);
				strcat(filepulse,val);
				pulseFile = fopen (filepulse,"w");
				for (int j=0; j<vgslout2->size2; j++)
				{
					sprintf(val,"%.14e",gsl_matrix_get(vgslout2,i,j));
					strcat(val,"\n");
					fputs (val,pulseFile);
				}
				fclose(pulseFile);
			}
			fclose(tstartFile);
		}*/

		//  Creating Tstart Column
		if (status == EPOK) {
			//MC obj.inObject = trgExten;
			obj.inObject = trgObject;
			obj.nameTable = new char [255];
			strcpy(obj.nameTable,"EUR-TRG");
			obj.iniRow = totalpulses;
			obj.endRow = totalpulses+nPulsesRow-1;
			obj.iniCol = 0;
			obj.nameCol = new char [255];
			strcpy(obj.nameCol,"Tstart");
			//MC obj.type = DAL_DOUBLE;
			obj.type = TDOUBLE;
			obj.unit = new char [255];
			strcpy(obj.unit,"seconds");
		}
		temp = gsl_vector_subvector(tstart,0,nPulsesRow);

		status = writeFitsSimple(obj, &temp.vector, status);

		// Creating I0 Column
		if (writePulse){	// If the input parameter writePulses = yes then it puts in output FITS file
			if (status == EPOK) {
				strcpy(obj.nameCol,"I0");
				//MC obj.type = DAL_DOUBLE;
				obj.type = TDOUBLE;
				strcpy(obj.unit,"Amps");
			}
			status = writeFitsComplex(obj, vgslout2, status);

		}

		//gsl_matrix_get_row(*pulsetemplate,vgslout2,0);
		gsl_matrix_memcpy(*pulses,vgslout2);

		// Creating Tend Column
		if (status == EPOK) {
			strcpy(obj.nameCol,"Tend");
			//MC obj.type = DAL_DOUBLE;
			obj.type = TDOUBLE;
			strcpy(obj.unit,"seconds");
		}
		temp = gsl_vector_subvector(tend,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector, status);


		// Creating tauRise Column
		if (status == EPOK) {
			strcpy(obj.nameCol,"tauRise");
		}
		temp = gsl_vector_subvector(taurise,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector, status);


		// Creating tauFall Column
		if (status == EPOK) {
			strcpy(obj.nameCol,"tauFall");
		}
		temp = gsl_vector_subvector(taufall,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector, status);

		// Creating difTstrt Column
		if (status == EPOK) {
			strcpy(obj.nameCol,"difTstrt");
		}
		status = writeFitsSimple(obj, difTstart, status);

		// Creating EstEnrgy Column
		if (status == EPOK) {
			strcpy(obj.nameCol,"EstEnrgy");
			strcpy(obj.unit,"eV");
		}
		temp = gsl_vector_subvector(energy,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector, status);

		// Creating Quality Column
		if (status == EPOK) {
			strcpy(obj.nameCol,"Quality");
			//MC obj.type = DAL_SHORT;
			obj.type = TSHORT;
			strcpy(obj.unit,"bits");
		}
		temp = gsl_vector_subvector(quality,0,nPulsesRow);
		status = writeFitsSimple(obj, &temp.vector, status);

		// Creating TailBfr Column
		if (status == EPOK) {
			strcpy(obj.nameCol,"TailBfr");
			//MC obj.type = DAL_SHORT;
			obj.type = TSHORT;
			strcpy(obj.unit," ");
		}
		status = writeFitsSimple(obj, tailBeforegsl, status);

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


/*int obtainTauRise (double Tmax, double Imax, double tauFall, double *tauRise, int status)
{
	if (!status == 0)	return status;

	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	//double x_lo = 1/tauFall, x_hi = Imax/1e-6;
	double x_lo = 1e5, x_hi = 2e5;
	cout<<"Imax: "<<Imax<<endl;
	cout<<"x_lo: "<<x_lo<<" x_hi: "<<x_hi<<endl;
	cout<<"x_loOK: "<<1/tauFall<<" x_hiOK: "<<Imax/1e-6<<endl;
	double r = (x_hi-x_lo)/2.0;
	gsl_function F;
	struct my_f_params params = { 1.0/tauFall, Tmax };

    F.function = &my_f;
	F.params = &params;

	T = gsl_root_fsolver_bisection;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);

	writeLog(fileRef,"Warning", verbosity,"Using %s method", gsl_root_fsolver_name (s));

	writeLog(fileRef,"Warning", verbosity,"%5s [%9s, %9s] %9s","iter", "lower", "upper", "root");

	do
	{
		iter++;
	    status = gsl_root_fsolver_iterate (s);
	    r = gsl_root_fsolver_root (s);
	    x_lo = gsl_root_fsolver_x_lower (s);
	    x_hi = gsl_root_fsolver_x_upper (s);
	    status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);

	    if (status == GSL_SUCCESS)	writeLog(fileRef,"Warning", verbosity,"Converged: ");

	    //status = writeLog(fileRef,"Warning", verbosity,status,"%5d [%.7f, %.7f] %.7f",iter, x_lo, x_hi, r);
	    writeLog(fileRef,"Warning", verbosity,"%d [%f, %f] %f",iter, x_lo, x_hi, r);

	}while (status == GSL_CONTINUE && iter < max_iter);

	*tauRise = 1/r;

	gsl_root_fsolver_free (s);*/

	/*int iter = 0, max_iter = 100;
	const gsl_root_fdfsolver_type *T;
	gsl_root_fdfsolver *s;
	double alpha = 1.0/tauFall;
	//double x0, x = 0.005e-3, r_expected = Tmax-(1/x)*log((alpha+x)/alpha);
	double x0, x = alpha;
	//double x0, x = 145000;
	gsl_function_fdf FDF;
	struct my_f_params params = {alpha, Tmax};

	FDF.f = &my_f;
	FDF.df = &my_f_deriv;
	FDF.fdf = &my_f_fdf;
	FDF.params = &params;

	T = gsl_root_fdfsolver_newton;
	s = gsl_root_fdfsolver_alloc (T);
	gsl_root_fdfsolver_set (s, &FDF, x);

	//status = writeLog(fileRef,"Warning", verbosity,status,"Using %s method", gsl_root_fsolver_name (s));

	writeLog(fileRef,"Warning", verbosity,"%5s %9s","iter", "root");

	do
	{
		iter++;
	    status = gsl_root_fdfsolver_iterate (s);
	    x0 = x;
	    x = gsl_root_fdfsolver_root (s);
	    status = gsl_root_test_delta (x, x0, 0, 1e-3);

	    if (status == GSL_SUCCESS)	status = writeLog(fileRef,"Warning", verbosity,status,"Converged: ");

	    status = writeLog(fileRef,"Warning", verbosity,status,"%d %f",iter, x);
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	*tauRise = 1/x;

	gsl_root_fdfsolver_free (s);*/

	/*return(status);
}*/


/***** SECTION 10 ************************************************************
* readLib: This function loads all the pulse models from the pulses models library input FITS file.
******************************************************************************/
int readLib (int status)
{
	string message="";
	//MC 
	int extver=0;

	//MC status = DALtableGetNumRows(inLibExten, &nummodels,status);
	strcpy(extname,"EUR-LIB");
	if(fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status)) printerror(status);
	if(fits_get_num_rows(inLibObject,&nummodels, &status)) printerror(status);

	/////////////////////// Get structure of pulse model library input FITS file columns
	//MC long numBins = 0L;	dal_dataType type;	int varLength=0; long sizeT, sizeV;
	strcpy(straux,"ENERGY");
	//MC status = DALtableGetColStruct (inLibExten,0,&numBins,straux,&type, &sizeT,&varLength,status);	if (status == -1219){status = -62507; EPexit(status);}
	if(fits_get_colnum(inLibObject,0,straux,&colnum,&status)) printerror(status);

	strcpy(straux,"ESTENERGY");
	//MC status = DALtableGetColStruct (inLibExten,0,&numBins,straux,&type, &sizeT,&varLength,status);	if (status == -1219){status = -62503; EPexit(status);}
	if(fits_get_colnum(inLibObject,0,straux,&colnum,&status)) printerror(status);
	strcpy(straux,"PULSE");
	//MC status = DALtableGetColStruct (inLibExten,0,&numBins,straux,&type, &sizeV,&varLength,status);	if (status == -1219){status = -62506; EPexit(status);}
	if(fits_get_colnum(inLibObject,0,straux,&colnum,&status)) printerror(status);		

	extern int inDataIteratorLib(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

		// Create structure to run Iteration DAL
	iteratorCol colsLib [3]; 				// Structure of Iteration DAL
	int n_cols = 3; 						// Number of columns:  Energy + EstEnergy + PULSE
	long rows_per_loop = nummodels; 		// 0: Use default: Optimum number of rows
	long offset=0; 							// 0: Process all the rows

	// Read Energy Column
	strcpy(straux,"Energy");
	//MC status = DALiterSetByName(&colsLib[0], inLibExten, straux, DAL_DOUBLE, InputCol, status);
	if(fits_iter_set_by_name(&colsLib[0], inLibObject, straux, TDOUBLE, InputCol))printerror(status);
	
		// Read Est Energy Column
	strcpy(straux,"EstEnergy");
	//MC status = DALiterSetByName(&colsLib[1], inLibExten, straux, DAL_DOUBLE, InputCol, status);
	if(fits_iter_set_by_name(&colsLib[1], inLibObject, straux, TDOUBLE, InputCol))printerror(status);
	
		// Read PULSE
	strcpy(straux,"PULSE");
	//MC status = DALiterSetByName(&colsLib[2], inLibObject, straux, DAL_DOUBLE, InputCol, status);
	if(fits_iter_set_by_name(&colsLib[2], inLibObject, straux, TDOUBLE, InputCol))printerror(status);

		// Called iteration function
	//MC status = DALiterateData(n_cols, colsLib, offset, rows_per_loop, inDataIteratorLib, 0L, status);
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
	//MC if (!status == 0)	return status;
	int status=EPOK;
	int extver=0;

	/////////////////////// Declare variables
	static double *energy;		// Vector of ENERGY column
	static double *estenergy;		// Vector of ESTENERGY column
	static double *pulsemodel;  // Vector of PULSE column
	long eventsz;

	//MC status = DALattributeGetInt(inLibExten, "EVENTSZ", &eventsz, unit, comment, status);
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
	//MC status = DALiterGetArray(&cols[0], ((void **)&energy), status);
	energy = (double *) fits_iter_get_array(&cols[0]);
	//MC status = toGslVector(((void **)&energy), &energygsl, nrows, 0, DAL_DOUBLE, status);
	status = toGslVector(((void **)&energy), &energygsl, nrows, 0, TDOUBLE, status);
	//MC status = DALiterGetArray(&cols[1], ((void **)&estenergy), status);
	estenergy = (double *) fits_iter_get_array(&cols[1]);
	//MC status = toGslVector(((void **)&estenergy), &estenergygsl, nrows, 0, DAL_DOUBLE, status);
	status = toGslVector(((void **)&estenergy), &estenergygsl, nrows, 0, TDOUBLE, status);
	//MC status = DALiterGetArray(&cols[2], ((void **)&pulsemodel), status);
	pulsemodel = (double *) fits_iter_get_array(&cols[2]);
	//MC status = toGslMatrix(((void **)&pulsemodel), &pulsemodelgsl, eventsz, nrows, (int)DAL_DOUBLE,	0, status);
	status = toGslMatrix(((void **)&pulsemodel), &pulsemodelgsl, eventsz, nrows, (int)TDOUBLE,0,status);

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
int createLibrary(int status)
{
	//MC 
	int extver=0;

	//MC if (!status == 0)	return status;

	/////////////////////// Create output FITS files: If it does not exist already
	//Create Library file or open it
	
	//MC status = DALobjectOpen(inLibName, &inLibObject, status);
	fits_open_file(&inLibObject, inLibName,0,&status);

	if (status == EPOK)
	{
		append = true;
		//MC status = DALobjectFindElement(inLibObject, "EUR-LIB", &inLibExten, status);
		//if (status == -2104){EPexit(-62202);}
		strcpy(extname,"EUR-LIB");
		if(fits_movnam_hdu(inLibObject, ANY_HDU,extname, extver, &status)) printerror(status);
					
		//MC status = DALattributeGetInt(inLibExten, "EVENTCNT", &eventcntLib, unit, comment, status);
		strcpy(keyname,"EVENTCNT");
		if(fits_read_key(inLibObject,TLONG,keyname, &eventcntLib,comment,&status)) printerror(status);
		//MC status = DALattributePutInt(inLibExten, "EVENTCNT", eventcntLib+1, unit, comment, status);
		eventcntLib1 = eventcntLib + 1;
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(inLibObject,TLONG,keyname, &eventcntLib1,comment,&status))printerror(status);

		status = readLib(status);
		ntotalrows = 1;
	}
	else
	{
		append = false;
		status = EPOK;
		//MC status = DALobjectCreate(inLibName, DAL_DISK, NULL, &inLibObject, status);
		if(fits_create_file(&inLibObject, inLibName, &status))printerror(status);

		/////////////////////// Create extension EUR-LIB
		//MC status = DALtableCreate(inLibObject, DAL_DISK, 0, "EUR-LIB", &inLibExten, status);
		strcpy(extname,"EUR-LIB");
		if(fits_create_tbl(inLibObject,BINARY_TBL,0,0,ttype,tform,tunit,extname,&status)) printerror(status);
		//MC status = DALattributePutInt(inLibExten, "EVENTCNT", 1, unit, comment, status);
		evtcnt = 1;
		strcpy(keyname,"EVENTCNT");
		if(fits_write_key(inLibObject,TLONG,keyname,&evtcnt,comment,&status)) printerror(status);

		/////////////////////// Create keywords
		//MC status = DALattributePutInt(inLibExten, "LIB_ID", 1, unit, comment, status);
		lib_id = 1;
		strcpy(keyname,"LIB_ID");
		if(fits_write_key(inLibObject,TLONG,keyname,&lib_id,comment,&status)) printerror(status);
		//MC status = DALattributePutChar(inLibExten, "CREATOR", create, unit, comment, status);
		strcpy(keyname,"CREATOR");
		strcpy(keyvalstr,create);
		if(fits_write_key(inLibObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);
		//MC status = DALattributePutChar(inLibExten, "FTYPE", "lib", unit, comment, status);
		strcpy(keyname,"FTYPE");
		strcpy(keyvalstr,"lib");
		if(fits_write_key(inLibObject,TSTRING,keyname,keyvalstr,comment,&status)) printerror(status);
	}

	char str_energy[125];       sprintf(str_energy,"%f",energy);

	string process (string("TRIGGER") 	+ ' ' +
	string(inName) 		+ ' ' + string(inLibName) 	  + ' ' +
	string(str_energy)      + ' ' +
	string("(")				+      (string) create 		  +   	  string(")"));

	//MC status = DALattributePutChar(inLibExten, "PROCESS",process.c_str(), unit, comment, status);
	strcpy(keyname,"PROCESS");
	strcpy(keyvalstr,process.c_str());
	if(fits_write_key_longwarn(inLibObject,&status))printerror(status);
	if(fits_write_key_longstr(inLibObject,keyname,keyvalstr,comment,&status)) printerror(status);
	
	//MC close file if newly created
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
int writeLibrary(double estenergy, gsl_vector **pulsetemplate, int status)
{
	//MC if (!status == ISDC_OK)	return status;

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

    	//MC obj.inObject = inLibExten;
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
    	    //MC obj.type = DAL_DOUBLE;
	    obj.type = TDOUBLE;
    	    strcpy(obj.unit,"eV");
    	    gsl_vector_set (energyoutgsl,0,gsl_vector_get(energycolumn,i));
    	    status = writeFitsSimple(obj, energyoutgsl, status);

    	    strcpy(obj.nameCol,"ESTENERGY");
    	    strcpy(obj.unit," ");
    	    gsl_vector_set (estenergyoutgsl,0,gsl_vector_get(estenergycolumn,i));
    	    status = writeFitsSimple(obj, estenergyoutgsl, status);

    	    strcpy(obj.nameCol,"PULSE");
    	    strcpy(obj.unit,"Amps");
    	    gsl_matrix_get_row(modelsrow,modelsaux,i);
    	    gsl_matrix_set_row(pulsetemplates_matrix,0,modelsrow);
    	    status = writeFitsComplex(obj, pulsetemplates_matrix, status);
    	}
    	gsl_vector_free(modelsrow);

    	gsl_vector_free (energycolumn);
    	gsl_vector_free (estenergycolumn);
    	gsl_matrix_free(modelsaux);
    }
    else
    {
    	//MC status = DALattributePutInt(inLibExten, "EVENTSZ", sizePulse_b, unit, comment, status);
	if(fits_write_key(inLibObject,TINT,"EVENTSZ",&sizePulse_b,comment,&status)) printerror(status);

    	gsl_vector_set (energyoutgsl,0,energy);
    	gsl_vector_set (estenergyoutgsl,0,estenergy);

    	gsl_matrix_set_row(pulsetemplates_matrix,0,*pulsetemplate);

    	// Creating ENERGY Column
   		//MC obj.inObject = inLibExten;
		obj.inObject = inLibObject;
   		obj.nameTable = new char [255];
   		strcpy(obj.nameTable,"EUR-LIB");
   		obj.iniRow = 1;
   		obj.endRow = 1;
   		obj.iniCol = 0;
   		obj.nameCol = new char [255];
   		strcpy(obj.nameCol,"ENERGY");
   		//MC obj.type = DAL_DOUBLE;
		obj.type = TDOUBLE;
   		obj.unit = new char [255];
   		strcpy(obj.unit,"eV");
   		status = writeFitsSimple(obj, energyoutgsl, status);

   		// Creating ESTENERGY Column
		strcpy(obj.nameCol,"ESTENERGY");
		strcpy(obj.unit,"u1");
		status = writeFitsSimple(obj, estenergyoutgsl, status);

   		// Creating PULSE Column
   		strcpy(obj.nameCol,"PULSE");
   		strcpy(obj.unit,"Amps");
   		status = writeFitsComplex(obj, pulsetemplates_matrix, status);
    }

		// Free memory
    gsl_vector_free(energyoutgsl);
    gsl_vector_free(estenergyoutgsl);

    gsl_matrix_free (pulsetemplates_matrix);

    return (status);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
