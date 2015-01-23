/************************************************************************************************

   This file is part of SIXTE/SIRENA software

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
*                      PSEUDOENERGY
*
*  File:      pseudoenergy.cpp
*  Version:   15.0.0
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
* version beta: 	20/09/06	First version
* version 1.0.0 	04/10/06	Error Processing
* version 1.1.0		18/10/06	Energyresol only calculate pulse parameters. NOT calculate fit parameters.
*  							   	Add Base Line value to output fits
* 							   	Change in Base Line function.
* 								Now, it's based in mean of values of pretrigs
*							   	Adding new keyword in output fits file: ERS_ID
* 					 		   	Modify shiftright function
* version 1.1.1 	01/12/06 	Delete shiftleft function
* version 1.1.2 	05/12/06 	Delete shiftright function
* version 1.2.0		12/12/06 	The holzgauss module is added to energy resolution module
* version 1.3.0   	19/02/07 	The holzgauss module is deleted to energy resolution module
* version 1.3.1   	06/03/07 	Change the name of Group Extension of input Fits file for "EUR-DMX"
* 							   	Change the name of Group Extension of output Fits file
* 							  	for "EUR-FILTER" and "EUR-ENERGY"
* version 1.3.2 	20/04/07 	Delete environment value "FITS_DIR", it's incompatible with the browse of param_gui.
* 							   	"FIT_DIR" imposes a input directory which contains the input fits file (That is incompatible)
* 							   	with the browse).
*	  						   	Change input Name of the fits file. Now, the Name include the .fits extension. Now, the file
* 							   	Name is compatible with the browse of param_gui.
* 							   	Change in the Call to "readFitsComplex", "readFitsSimple", "writeFitsComplex" and
*                			   	"writeFitsSimple". Now is used the struct "OIData". See documentation of Utils module. (inoutUtils).
* version 1.3.3		24/04/20 	Include new input parameter in the module. Adding in .par file "output File". Now,
* 							   	the user can choose the name of the output .fits file.
* 							   	-jjs-
* 							    Change maximum number of data channels in an input file (4)
* 							    Parametrize the names for the input ADC channels, and
*          						the names for the result data column names.
* 								If input file is NOT module	RILlogMessage(fileRef,Log_3,);ted, do not demodulate, but do all else, eg. scaling the data
* 								ADC channel names are now parametrized
* version 1.3.4 	07/05/07	If exist more than one columns in the input fits file --> There is a ERROR (new error)
* 								Add error: Status
* 								Otherside is to answer about the name of the columns which you want to use.
* version 1.4.0		16/05/07 	Modification of the module using the new module inoutFits included in Utils. Now, Reads and Writes
* 								of Fits use block of rows.
* 								Including comments using RIL logs.
* version 1.4.1		12/06/07	Modify error processing
* version 2.0.0   	12/03/08	Optimum read and write to fits file
* 								Energyresol reads from trg fits file and flt fits file
* 								Adding the function "init module"
*								Adding two inputs parameters: "Central Energy" and "house keeping"
* 								Deleting the keyword "yscale"
*								Deleting the column "Baseline"
* version 3.0.0		07/04/08	Filtering of pulse using the column "Quality"
* 								Adding the output column "quality"
* version 3.1.0		26/05/08	Adding input/output keyword "ql"
* 								Adding Error Processing: "The input FITS file hasn't run in PSH module."
* 								New function to filter pulses not valid of TRG fits file
* 								Adding function "polyFit"
* 								New Software documentation.
* 								Adding function "energyCalibration"
* 								Create a HISTORY in the pulse shape file.
* 								Deleting input/output keyword "ql";
* 								Adding input parameter "ql"
* version 4.0.0		10/09/08	New documentation in the source code
* 								Deleting PolyFit function. It has been include in "miscellaneous" module.
* 								Adding input parameter "nameLog"
* 								Adding input parameter "verbosity" and processing it using RIL libraries.
* 								Adding new function "readInputKeywords"
* 								Deleting in/output keyword "channelcount"
* 								Including error codes.
* version 4.1.0 	29/09/08	Included "CREATE" and "PROCESS" keyword in the output FITS file.
* 								Include error processing of "INPUT_PARAMETER_NON_DEFINE"
* version 4.2.0 	24/10/08	If the module fail. It will be warning in the output.
*  								Included parameter time to measure.
* version 5.0.0		14/09/08 	Changed the output column name of Current from "V" to "I0"
* 								Included new input keywords.
* version 5.1.0		22/01/09	Modified error code of input parameters
* version 5.2.0		02/02/09	Remove remove Hkfile input parameter.
* 								Remove Central Energy input parameter.
* 								Remove functionality of this task: convert from pseudoenergies to eV.
* 								now, it's included in holzgauss task.
* 								Removed energyResol and energyCalibration functions.
* version 5.2.1		04/02/09	Checked free memory allocate.
* version 6.0.0		20/02/09	Included, modified and removed input keywords
* version 6.1.0		24/02/09	Xraychain can run input FITS files of type (XRAY, TESNOISE and IV)
* 								Solved bug 4: stops on zero value for the keyword TIMEZERO, error -62436.
* version 6.2.0		12/03/09	Changed the documentation.
* version 7.0.0                 Used the keywords libraries
* version 7.0.1     16/04/09    Delete "PILInit(argc, argv);" and "PILClose(PIL_OK);"
* 								Documentation updated
* version 7.2.0     21/05/09    Cmono3_10cps_pix2016_tf150_wn_156khz_200s_trgNEW2015.fitsorrected energy calculation in findEnergy
* version 8.0.0		27/05/09    Pseudoenergy_INT column added in _psh FITS file calculated by integrating the Norris function
* version 9.0.0     02/06/09	Pseudoenergy_INT column only added in _psh FITS file if a Norris fit has been made in PULSESHAPE
*                               If Norris fit columns not in PSH file, not read here
* version 9.1.0     07/07/09    Bug resolved: Gamma function argument could not be greater than 171.
* version 9.1.1     15/07/09    Changed order in delete's
* version 10.0.0    14/09/09    Pseudoenergy_INT is now well-written
* 								findEnergy changed (gsl_vector_set (multgsl,sgSzL - i-1,mult); instead of gsl_vector_set (multgsl,sgSzL - i,mult);)
* version 11.0.0	xx/xx/09	If error file exists, its contests will be deleted
*                   30/11/09    The PROCESS keyword written in the modified PSH file keeps the PROCESS keyword read from the
*                               PSH input FITS file and the PROCESS generated in ENERGYRESOL and the ENERGYRESOL version are added
*                   30/09/10    Restructuring of the columns of the input FITS files from TRIGGER:
*                               	Old: TIME-I0-Tstart-Tend-tauRise-tauFall-Baseline-Sigma-EndPulse-Quality
*                                   New: Tstart-I0-Tend-tauRise-tauFall-difTstrt-Quality
*                               Restructuring the columns of the input FITS files from PULSESHAPE:
*                               	Old: TIME-Tstart-Tend-Baseline-Sigma-MaxTime-MaxCurrent-TRise-Nrise-TFall-Nfall-...
*                                        ...-MaxTimeFit-MaxCurrentFit-TRiseFit-NriseFit-TFallFit-NfallFit-ChisqFit-...
*                                        ...-BaselineFit-Quality
*                                   New: Tstart-Tend-MaxTime-MaxCurrent-TRise-Nrise-TFall-Nfall-Quality
*                               Energy is not calculated by integrating the Norris function because the Norris fitting
*                               is not made in PULSESHAPE (delete inDataIteratorPsh)
*                   01/10/10    'eventsz_flt' is not necessary (EVENTSZ from TRIGGER is used)
*                   08/11/10	As well as calculating the pseudoenergy of pulses, now the pseudoenergy of intervals without pulses is calculated
*                   			The pseudoenergy of each pulse is estimated using non-filtered pulses instead of filtered pulses
*                   			New input parameters "Tinitial" & "Tfinal"
*                   			Now is possible to select the filter template (not necessarily the calculated one in the previous task)
*                   			New input parameter "inFile" --> inputfile to XRAY is necessary
*                   			"difTstrt" column in EUR-TRG extension of TRG file is necessary
*                   			Added new functions: "inDataIterator", "writeEnergy" & "writeKeywords"
*                   			The pseudoenergy of each pulse or interval without pulses is written in a new extension in PULSESHAPE FITS (EUR-ENR)
*                   			--> Pseudoenergy" column
*                   			Three columns more are added: "Tstart", "Tend" & "value"
*                   			Documentation updated
* version 12.0.0	21/12/10	"EVENTCNT" Keyword is added in the header of "EUR-ENR" Extension
* 								'Log' argument instead 'OK' in some "writeLog" functions
* 								Documentation updated
*                   23/02/11    I0 instead I0_NotFiltered
*                   28/02/11    PROCESS also includes ENERGYRESOL version (CREATOR)
*                   01/03/11    CREATOR not overwritten in _psh (ERS version info in HISTORY)
*                   08/03/11	'row_input' allocated in 'main' instead in 'inDataIterator'
*                   21/03/11    "columnName" input parameter not used
* version 13.0.0    25/03/11    Quick Look mode no longer used
*                   05/04/11    No pulses in _trg input FITS file => It only reads the input FITS file and calculates the pseudoenergy
*                   24/08/11    'gsl_vector_view temp;' added to handle with subvectors
*                   18/10/11    New possible way of calculating the convolution in findEnergyNEW (new 'calculateConvolution' function)
*                               New input keyword "chngplrt" from EUR_TRG
*                               Depending on chngplrt pulse-free intervals are multiplied or not by -1.0
*                               (according to the pulse polarity modified in TRIGGER)
*                   04/11/11    Added 'if (gsl_vector_get(Tendrow,0)-gsl_vector_get(Tstartrow,0) > 0)' to write
*                               the corresponding row in the _psh output FITS file
*                   			It is not necessary that a pulse-free interval was longer than the filter
*                               template to calculate the pseudoenergy
*                               ('if(size_trg>fltSz)' and 'if(eventsz_in>fltSz)' deleted)
*                               Pulse-free intervals must have longer than 0 to calculate the pseudoenergy
*                               ('if(size_trg!=0)' changed to 'if(size_trg>0)')
*                               Pseudoenergy of non-valid pulses ('quality' different from 0) is calculated
*                               Pulses pseudoenergy is not obliged to be positive
*                               (deleted 'if (gsl_vector_get(energygsl,0) < 0){gsl_vector_set(valuegsl,0,-1);}')
*                   09/11/11    New column 'Quality' in the EUR-ENR extension:
*                               	Pulse-free intervals will always have quality=0
*	                                Pulses will have the same quality as in EUR-PSH extension
*	                18/11/11    Change the point where an event of the input FITS file is read
*	                            (when it is the first pulse in an event, the last pulse-free interval of the
*	                            previous event has to be analyzed => New event read after that)
*	                19/12/11    Changes in calculateConvol:
*	                                 Pulse-free intervals => Middle-point of the convolution => Pseudoenergy
*	                                 Pulses => Maximum of the convolution => Pseudoenergy
*	                26/01/12    Choose the domain to filter (time or frequency) => New functionality in 'findEnergy'
*	                20/02/12    Running sum filter
*	                13/03/14    VALUE column renamed as PLSORNOT (PuLSe OR NOT pulse) column
*	                17/04/14    Read the column OF_L from the _trg FITS file
*	                14/05/14    Read the column NRMFCTR from the _trg FITS file
*	                05/07/14    Deleted input parameter 'numshift' (or the whole convolution is calculated or the convolution is not used)
*	                            FREQ and OPTIMALFF columns from _flt FITS file (not used)
*	                            Pseudoenergies of pulse-free intervals are not calculated
*	                            Deleted 'inDataIterator_NoPulsesTrg' => If there are no pulses, the routine finishes
*	                            Deleted output PLSORNOT column because there is no necessary distinguish between pulses
*	                            and pulse-free intervals
*	                            ANNALS added to the _psh FITS file
* 			        07/07/14    PIL/RIL/Common dependencies removed. new functions used to read input parameters
* version 14.0.0    10/07/14    EUR-ENR no longer used => Pseudenergy column added to the EUR-PSH extension
*                   11/07/14    Solved bug preventing the task from reading the full command line
* version 14.0.1    11/07/14    Comments of the 'initModule' function modified
*                   15/07/14    'initModule' modified in order to accept negative int or double (if parameters are read from command line)
*                   31/07/14    Comments referring to PIL and RIL deleted
*                   16/09/14    ADC2016 instead RAW_ADC_DATA
*                               PXL02016 instead I0 column
*                               Differences when reading or writing some keywords
*                   18/09/14    RECORDS instead ADC2016
*                               ADC instead PXL02016 column
* version 15.0.0 	  Dec/14    DAL->CFITSIO migration
*			                    Deleted some unnecessary functions
*			                    Deleted some unnecessary input parameters
*			                    Deleted some unnecessary keywords
*		09/01/15    'short' instead of 'int' when QUALITY and TAILBFR columns are read (in 'inDataIteratorTrg')
*		12/01/15    'PROCESS' input keyword from _psh input FITS file is not read
*		            'PROCESS' output keyword in _psh output FITS file divided in two ones ('PROC0' generated by PULSEGRADE and 'PROC1'
*		            generated by PSEUDOENERGY)
*		            'ANNALS' input keyword from _psh input FITS file is not read
*		            'ANNALS' output keyword in _psh output FITS file renamed as 'MOD0'
* 		20/01/15    Renaming extensions: EUR-TRG->TRIGGER, EUR-PSH->PSGRADE
* 		            Renaming variables (from energyresol->pseudoenergy change)
*                           Add clobbering
*  		22/01/15    Renaming of some input parameters: 
* 				    TorF ->filterDomain ;  Hp_OForRS -> filterHp ; Mp_OForRS -> filterMp ; Ms_OForRS -> filterMs
* 				    Lp_OForRS -> filterLp ; Ls_OForRS -> filterLs
*                   Corrected bug in dimensioning a matrix using eventcnt_in instead of eventcnt_flt
*******************************************************************************************/

/******************************************************************************
DESCRIPTION:

The goal of the PSEUDOENERGY task is to estimate the pseusoenergy of each found pulse.

Depending on the grade of the pulses the pseudoenergy could be calculated by using the running sum filter or by
filtering in time domain (by calculating the appropriate integral) or in frequency domain (by multiplying the
Fourier transform of the pulse and the Fourier transform of the optimal filter).

In time domain, E = int(0,inf){p(t)f(t)df} when we assume that the pulse and the template are well aligned (p(t) -> Pulse
and f(t) -> Optimal filter).

In the frequency domain, the product of the Fourier transforms of the pulse and the optimal filter is employed.
The optimal filter expression in time domain is read from _trg(.fits) and converted into frequency domain.

	E = int(-inf,inf){|P(f)|^2df} = int(-inf,inf){|P(f)||F(f)|df}
	P(f) -> FFT of the pulse
	F(f) -> FFT of the optimal filter (matched to the pulse => P(f)~F(f)))
	|P(f)|^2 ~ |P(f)||F(f)| => E = int(-inf,inf){|P(f)||F(f)|df}

The calculated pseudoenergy is written in the _psg input/output FITS file (adding a new column Pseudoenergy)

The user must supply the following input parameters:

- inFile: Name of the input FITS file to the chain
- trgFile: Name of the _trg input FITS file which contains the found pulses
- fltFile: Name of the _flt input FITS file which contains the filter template (only in calibration mode)
- psgFile: Name of the _psg input/output FITS file which contains the characteristics of the found pulses
- filterDomain: Optimal filtering will be done in time (T) or in frequency (F) domain
- filterHp: Optimal filter (OP) or running sum filter (RS) will be applied for High-Res Primary pulses
- filterMp: Optimal filter (OP) or running sum filter (RS) will be applied for Medium-Res Primary pulses
- filterMs: Optimal filter (OP) or running sum filter (RS) will be applied for Medium-Res Secondary pulses
- filterLp: Optimal filter (OP) or running sum filter (RS) will be applied for Low-Res Primary pulses
- filterLs: Optimal filter (OP) or running sum filter (RS) will be applied for Low-Res Secondary pulses
- namelog: Output log file name
- verbosity: Verbosity level of the output log file

MAP OF SECTIONS IN THIS FILE:

 - 1. INCLUDE's
 - 2. MAIN
 - 3. initModule
 - 4. inDataIteratorTrgFLT
 - 5. inDataIteratorTrg
 - 6. inDataIteratorFlt
 - 7. inDataIterator
 - 8. findEnergy
 - 9. RS_filter
 - 10. writeEnergy
 - 11. writeKeywords

*******************************************************************************/

/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include "pseudoenergy.h"


/***** SECTION 2 ************************************
* MAIN function: This function is the main function of the PSEUDOENERGY task
*
* - Read input parameters (call initModule)
* - Open input FITS fileS
* 		- Input FITS file to the chain
* 		- Input _trg FITS file
* 		- Input/output _psg FITS file
* - Read input keywords and check their values
* - If there are no pulses in _trg input FITS file (empty TRIGGER extension) => The task finishes
* - If calibration mode => Open _flt input FITS file
* - Convert Seconds into Samples
* - Get structure of input FITS fileS columns
* - If calibration mode => Iterator for _flt FITS file
* 	- Create structure to run Iteration: inDataIteratorFlt
*	- Read columns of Filter Template file (OPTIMALF column)
* 	- Called iteration function: inDataIteratorFlt
* - If production mode => Iterator for FILTER in TRIGGER
* 	- Create structure to run Iteration: inDataIteratorTrgFLT
*	- Read columns OPTIMALF and NRMFCTR
* 	- Called iteration function: inDataIteratorTrgFLT
* - Iterator for TRIGGER in TRIGGER
* 	- Create structure to run Iteration: inDataIteratorTrg
*	- Read columns I0, Tstart, Tend, difTstrt, Quality and Grade
*	- Called iteration function: inDataIteratorTrg
* 		- Call findEnergy
*       - Call writeEnergy
* - Create MOD0 in PULSESHAPE file
* - Call writeKeywords
* - Close input FITS files
* - Close output FITS files
* - Free allocate of memory
* - Finalize the task
****************************************************/
int main (int argc, char **argv)
{
	string message = "";
	create = "pseudoenergy v.15.0.0";			//Set "CREATOR" keyword of output FITS file
	time_t t_start = time(NULL);
	int status = EPOK, extver = 0;
	
	sprintf(temporalFileName,"PSEUDOENERGYauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if (temporalFile == NULL)
	{
	    message = "Cannot open auxiliary file " + string(temporalFileName);
	    writeLog(fileRef,"Error",verbosity,message);
	    EP_EXIT_ERROR(message,EPFAIL);	    
	}

	// Read input parameters
	if (initModule(argc, argv))
	{
		message = "Error in initModule";
		EP_EXIT_ERROR(message,EPFAIL);
	}

	writeLog(fileRef,"Log", verbosity,"Into Pseudoenergy Task");

	//  Open input FITS file to the chain
	if (fits_open_file(&inObject, inName,READWRITE,&status))
	{
	    message = "Cannot open file " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(extname,"RECORDS");
	if (fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	message = "Open inputFits: " + string(inName);
	writeLog(fileRef,"Log", verbosity,message);

	// Open _trg input FITS file
	if (fits_open_file(&trgObject, trgName,READWRITE,&status))
	{
	    message = "Cannot open file " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}	
	strcpy(extname,"TRIGGER");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	message = "Open TriggerFits: " + string(trgName);
	writeLog(fileRef,"Log", verbosity,message);

	// Open  input/output _psg FITS file
	if (fits_open_file(&psgObject, psgName,READWRITE,&status))
	{
	    message = "Cannot open file " + string(psgName);
	    EP_EXIT_ERROR(message,status);
	}	
	strcpy(extname,"PSGRADE");
	if (fits_movnam_hdu(psgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(psgName);
	    EP_EXIT_ERROR(message,status);
	}
	message = "Open PulseFits: " + string(psgName);
	writeLog(fileRef,"Log", verbosity,message);
	
	// Read input keywords and check their values
	strcpy(extname,"TRIGGER");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
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
	
	if (eventcnt == 0)	// There are no pulses in _trg input FITS file (empty TRIGGER extension)
	{
		time_t t_end = time(NULL);
		message = "There are no pulses in the input FITS file " + string(trgName) +" There is no output FITS file ";
		writeLog(fileRef,"Warning", verbosity,message);
		writeLog(fileRef,"OK", verbosity,"Pseudoenergy Module OK");
			
		sprintf(straux,"%f",(double) (t_end - t_start));
		message = "Time: " + string(straux);
		writeLog(fileRef,"Log", verbosity,message);
		
		if(fclose(fileRef))
		{
		    message = "Cannot close log file";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
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
	strcpy(keyname,"CHNGPLRT");
	if (fits_read_key(trgObject,TLONG,keyname, &chngplrt,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	if ((chngplrt != 0) && (chngplrt !=1))
	{
		message = "Legal values for CHNGPLRT (TRIGGER) are 0 or 1";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	
	strcpy(extname,"RECORDS");
	if (fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(keyname,"TRIGGSZ");
	if (fits_read_key(inObject,TLONG,keyname, &eventsz_in,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(inName);
	    EP_EXIT_ERROR(message,status);
	}
	if (eventsz_in <= 0)
	{
		message = "Legal values for TRIGGSZ (RECORDS) are integer numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	if (fits_get_num_rows(inObject,&eventcnt_in, &status))
	{
	    message = "Cannot access HDU " + string(extname) + " in " + string(inName) + " file (cannot get number of rows)";
	    EP_EXIT_ERROR(message,status); 
	}
	
	row_input = gsl_vector_alloc(eventsz_in);

	if (mode == 0)		// Calibration mode
	{
		// Open FILTER FITS
		if (fits_open_file(&fltObject, fltName,READWRITE,&status))
		{
		    message = "Cannot open file " + string(fltName);
		    EP_EXIT_ERROR(message,status);
		}	
		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(fltObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to " + string(extname) + " in " + string(fltName);
		    EP_EXIT_ERROR(message,status);
		}
		message = "Open FilterFits: " + string(fltName);
		writeLog(fileRef,"Log", verbosity,message);

		strcpy(keyname,"EVENTCNT");
		if (fits_read_key(fltObject,TLONG,keyname, &eventcnt_flt,comment,&status))
		{
		    message = "Cannot read key " + string(keyname) + " in " + string(fltName);
		    EP_EXIT_ERROR(message,status);
		}
		if (eventcnt_flt < 0)
		{
			message = "Legal values for EVENTCNT (FILTER in _flt) are non negative integer numbers";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		strcpy(keyname,"NRMFCTR");
		if (fits_read_key(fltObject,TDOUBLE,keyname, &normalizationFactor,comment,&status))
		{
		    message = "Cannot read key " + string(keyname) + " in " + string(fltName);
		    EP_EXIT_ERROR(message,status);
		}
		if (normalizationFactor <= 0)
		{
			message = "Legal values for NRMFCTR (FILTER in _flt) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}
	else if (mode == 1)		// Production mode
	{
		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to " + string(extname) + " in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}	
	}

	// Convert Seconds into Samples
	Lrs = LrsT*samprate;
	Lb = LbT*samprate;

	// Get structure of input FITS fileS columns
	strcpy(extname,"RECORDS");
	if (fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(inName);
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

	strcpy(extname,"TRIGGER");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"I0");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"TSTART");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"TEND");
	if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"difTstrt");
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

	if (mode == 0)
	{
		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(fltObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to " + string(extname) + " in " + string(fltName);
			EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"OPTIMALF");
		if (fits_get_colnum(fltObject,0,straux,&colnum,&status))
		{
			message = "Cannot get column number for " + string(straux) +" in " + string(fltName);
			EP_EXIT_ERROR(message,status);
		}
	}
	else if (mode == 1)
	{
		strcpy(extname,"FILTER");
		if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
		{
			message = "Cannot move to " + string(extname) + "in " + string(trgName);
			EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"OPTIMALF");
		if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
		{
			message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
			EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"NRMFCTR");
		if (fits_get_colnum(trgObject,0,straux,&colnum,&status))
		{
			message = "Cannot get column number for " + string(straux) +" in " + string(trgName);
			EP_EXIT_ERROR(message,status);
		}
	}

	int n_cols;
	long rows_per_loop;
	long offset;
	if (mode == 0)		// Calibration mode
	{
		// Iteration for FILTER
		filtergsl=gsl_vector_alloc(eventcnt_flt);
		i_filter=0;
		extern int inDataIteratorFlt(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols,void *user_strct );

		// Create structure to run Iteration: inDataIteratorFlt
		iteratorCol cols_flt [1];				// Structure of Iteration
		n_cols = 1; 							// Number of columns:  OPTIMALF
		rows_per_loop = 0;						// 0: Use default: Optimum number of rows
		offset = 0;								// 0: Process all the rows

		// Read OPTIMALF Column
		strcpy(straux,"OPTIMALF");
		status = fits_iter_set_by_name(&cols_flt[0], fltObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(fltName);
		    EP_EXIT_ERROR(message,status);
		}

		// Called iteration function: inDataIteratorFlt

		if (fits_iterate_data(n_cols, cols_flt, offset, rows_per_loop, inDataIteratorFlt,0L,&status))
		{
		    message = "Cannot iterate data using InDataIteratorFlt";
		    EP_EXIT_ERROR(message,status);
		}
		if (fits_close_file(fltObject,&status))
		{
		    message = "Cannot close file " + string(fltName);
		    EP_EXIT_ERROR(message,status);
		}
		
	}
	else if (mode == 1)		// Production mode
	{
		// Iteration for optimal filter of the triggered pulses
		extern int inDataIteratorTrgFLT(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols_trgFLT,void *user_strct );

		optfiltergslMATRIX = gsl_matrix_alloc(eventcnt,eventsz);		// GSL of input OPTIMALF column
		nrmfctrgsl = gsl_vector_alloc(eventcnt);

		// Create structure to run Iteration: inDataIteratorTrgFLT
		iteratorCol cols_trgFLT [2];				// Structure of Iteration
		n_cols = 2; 								// Number of columns: OPTIMALF + NRMFCTR
		rows_per_loop = 0;							// 0: Use default: Optimum number of rows
		offset  = 0;								// 0: Process all the rows

			// Read OPTIMALF
		strcpy(straux,"OPTIMALF");
		status = fits_iter_set_by_name(&cols_trgFLT[0], trgObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}

		strcpy(straux,"NRMFCTR");
		status = fits_iter_set_by_name(&cols_trgFLT[1], trgObject, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
		    EP_EXIT_ERROR(message,status);
		}

		// Called iteration function: InDataIteratorTrgFLT
		if (fits_iterate_data(n_cols, cols_trgFLT, offset, rows_per_loop, inDataIteratorTrgFLT,0L,&status))
		{
		    message = "Cannot iterate data using inDataIteratorTrgFLT";
		    EP_EXIT_ERROR(message,status);
		}
	}

	// Iteration for triggered pulses
	extern int inDataIteratorTrg(long totalrows,long offset,long firstrow,long nrows,int ncols,iteratorCol *cols_trg,void *user_strct );

	// Create structure to run Iteration: inDataIteratorTrg
	iteratorCol cols_trg [6];					// Structure of Iteration
	n_cols = 6; 								// Number of columns: I0, Tstart, Tend, difTstrt, Quality and Grade
	rows_per_loop = 0;							// 0: Use default: Optimum number of rows
	offset  = 0;								// 0: Process all the rows

	strcpy(extname,"TRIGGER");
	if (fits_movnam_hdu(trgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	
	// Read I0
	strcpy(straux,"I0");
	status = fits_iter_set_by_name(&cols_trg[0], trgObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}

	// Read Tstart, Tend, difTstrt, Quality and Grade columns
	strcpy(straux,"Tstart");
	status = fits_iter_set_by_name(&cols_trg[1], trgObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"Tend");
	status = fits_iter_set_by_name(&cols_trg[2], trgObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"difTstrt");
	status = fits_iter_set_by_name(&cols_trg[3], trgObject, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"Quality");
	status = fits_iter_set_by_name(&cols_trg[4], trgObject, straux, TSHORT, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"Grade");
	status = fits_iter_set_by_name(&cols_trg[5], trgObject, straux, TSHORT, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(trgName);
	    EP_EXIT_ERROR(message,status);
	}
	
	// Called iteration function: inDataIteratorTrg
	if (fits_iterate_data(n_cols, cols_trg, offset, rows_per_loop, inDataIteratorTrg,0L,&status))
	{
	    message = "Cannot iterate data using inDataIteratorTrg";
	    EP_EXIT_ERROR(message,status);
	}
	
	// Free allocate of GSL vectors
	gsl_vector_free (filtergsl);

	// Create MOD0 in _psg FITS file
	string mod0 (string("File MODIFIED by") + ' ' +	(string) create);
	strcpy(extname,"PSGRADE");
	if (fits_movnam_hdu(psgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(psgName);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(keyname,"MOD0");
	strcpy(keyvalstr,mod0.c_str());
	if (fits_write_key(psgObject,TSTRING,keyname, keyvalstr,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(psgName);
	    EP_EXIT_ERROR(message,status);
	}

	// Write keywords in _psg FITS file
	if (writeKeywords ())
	{
	    message = "Cannot run routine writeKeywords";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Close input FITS files
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

	// Close output FITS files
	if (fits_close_file(psgObject,&status))
	{
	    message = "Cannot close file " + string(psgName);
	    EP_EXIT_ERROR(message,status);
	}

	gsl_vector_free (row_input);

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
	
	writeLog(fileRef,"OK", verbosity,"Pseudoenergy Module OK");

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
int initModule(int argc, char **argv)
{
	int status = EPOK;
	string message = "";
	
	// Define PSEUDOENERGY input parameters and assign values to variables
	// Parameter definition and assignation of default values
	const int npars = 15, npars1 = 16;
	inparam pseudoPars[npars];
	int optidx =0, par=0, fst=0, ipar; 
	string task="pseudoenergy ";

	pseudoPars[0].name = "inFile";
    pseudoPars[0].description = "Input fits filename";
    pseudoPars[0].defValStr = "a.fits";
    pseudoPars[0].type =  "char";
	pseudoPars[0].ValStr = pseudoPars[0].defValStr;
	
	pseudoPars[1].name = "trgFile";
    pseudoPars[1].description = "TRIGGER fits filename";
    pseudoPars[1].defValStr = "a_trg.fits";
    pseudoPars[1].type =  "char";
	pseudoPars[1].ValStr = pseudoPars[1].defValStr;
	
	pseudoPars[2].name = "fltFile";
    pseudoPars[2].description = "FILTER fits filename (only in calibration mode)";
    pseudoPars[2].defValStr = "a_flt.fits";
    pseudoPars[2].type =  "char";
	pseudoPars[2].ValStr = pseudoPars[2].defValStr;
	
	pseudoPars[3].name = "psgFile";
    pseudoPars[3].description = "PULSESHAPE fits filename";
    pseudoPars[3].defValStr = "a_psg.fits";
    pseudoPars[3].type =  "char";
	pseudoPars[3].ValStr = pseudoPars[3].defValStr;
	
	pseudoPars[4].name = "filterDomain";
    pseudoPars[4].description = "(T/F):Optimal filtering will be done in time (T) or in frequency (F)";
    pseudoPars[4].defValStr = "T";
    pseudoPars[4].type =  "char";
	pseudoPars[4].ValStr = pseudoPars[4].defValStr;
	
	pseudoPars[5].name = "filterHp";
    pseudoPars[5].description = "(OP/RS): Filter applied to Hp events is OPTIMAL FILTER (OP) or RUNNING SUM (RS)";
    pseudoPars[5].defValStr = "OP";
    pseudoPars[5].type =  "char";
	pseudoPars[5].ValStr = pseudoPars[5].defValStr;

	pseudoPars[6].name = "filterMp";
    pseudoPars[6].description = "(OP/RS): Filter applied to Mp events is OPTIMAL FILTER (OP) or RUNNING SUM (RS)";
    pseudoPars[6].defValStr = "OP";
    pseudoPars[6].type =  "char";
	pseudoPars[6].ValStr = pseudoPars[6].defValStr;	
	
	pseudoPars[7].name = "filterMs";
    pseudoPars[7].description = "(OP/RS): Filter applied to Ms events is OPTIMAL FILTER (OP) or RUNNING SUM (RS)";
    pseudoPars[7].defValStr = "OP";
    pseudoPars[7].type =  "char";
	pseudoPars[7].ValStr = pseudoPars[7].defValStr;
	
	pseudoPars[8].name = "filterLp";
    pseudoPars[8].description = "(OP/RS): Filter applied to Lp events is OPTIMAL FILTER (OP) or RUNNING SUM (RS)";
    pseudoPars[8].defValStr = "OP";
    pseudoPars[8].type =  "char";
	pseudoPars[8].ValStr = pseudoPars[8].defValStr;
	
	pseudoPars[9].name = "filterLs";
    pseudoPars[9].description = "(OP/RS): Filter applied to Ls events is OPTIMAL FILTER (OP) or RUNNING SUM (RS)";
    pseudoPars[9].defValStr = "OP";
    pseudoPars[9].type =  "char";
	pseudoPars[9].ValStr = pseudoPars[9].defValStr;
	
	pseudoPars[10].name = "LrsT";
    pseudoPars[10].description = "Running sum length (in the RS filter case) (seconds)";
    pseudoPars[10].defValReal = 30.E-6;
    pseudoPars[10].type = "double";
    pseudoPars[10].minValReal = 0.;
    pseudoPars[10].maxValReal = 10.;
	pseudoPars[10].ValReal = pseudoPars[10].defValReal;
	
	pseudoPars[11].name = "LbT";
    pseudoPars[11].description = "Baseline averaging length (in the RS filter case) (seconds)";
    pseudoPars[11].defValReal = 1.E-3;
    pseudoPars[11].type = "double";
    pseudoPars[11].minValReal = 0.;
    pseudoPars[11].maxValReal = 10.;
	pseudoPars[11].ValReal = pseudoPars[11].defValReal;

	pseudoPars[12].name = "nameLog";
    pseudoPars[12].description = "Output log file name";
    pseudoPars[12].defValStr = "ers_log.txt";
    pseudoPars[12].type = "char";
	pseudoPars[12].ValStr = pseudoPars[12].defValStr;
	
	pseudoPars[13].name = "verbosity";
    pseudoPars[13].description = "Verbosity level of the output log file (in [0,3])";
    pseudoPars[13].defValInt = 3;
    pseudoPars[13].type = "int";
    pseudoPars[13].minValInt = 0;
    pseudoPars[13].maxValInt = 3;
	pseudoPars[13].ValInt = pseudoPars[13].defValInt;
	
	pseudoPars[14].name = "clobber";
    pseudoPars[14].description = "Re-write output files if clobber=yes";
    pseudoPars[14].defValStr = "no";
    pseudoPars[14].type = "char";
	pseudoPars[14].ValStr = pseudoPars[14].defValStr;

	// Define structure for command line options
	static struct option long_options[npars1];
	for (optidx = 0; optidx<npars; optidx++)
	{
		long_options[optidx].name= pseudoPars[optidx].name.c_str();
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
					if(long_options[optidx].name == pseudoPars[i].name.c_str())
					{
						if(pseudoPars[i].type == "char") //save char value for par
						{
						    pseudoPars[i].ValStr = optarg;
						}
						else // check if numeric value
						{
						    if ( (!isdigit(optarg[0]) && (optarg[0] != '-')) ||
						    		(!isdigit(optarg[0]) && (optarg[0] == '-') && (!isdigit(optarg[1]))))
						    {
						    	message = "Invalid value for input argument " + string(pseudoPars[i].name);
						    	EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
						    }
						    if (pseudoPars[i].type == "int")
						    {
						    	pseudoPars[i].ValInt = atoi(optarg);
						    }
						    else
						    {
						    	pseudoPars[i].ValReal= atof(optarg);
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
		if (interactivePars(pseudoPars,npars,task))
		{
		    message = "Error reading parameters interactively";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL); 
		}
	}

	// Save parameter values into meaningful variables
	for(int i=0;i<npars; i++)
	{
		if(pseudoPars[i].name == "inFile")
		{
			strcpy(inName, pseudoPars[i].ValStr.c_str());
		}
		else if(pseudoPars[i].name == "trgFile")
		{
			strcpy(trgName, pseudoPars[i].ValStr.c_str());
		}
		else if(pseudoPars[i].name == "fltFile")
		{
			strcpy(fltName, pseudoPars[i].ValStr.c_str());
		}
		else if(pseudoPars[i].name == "psgFile")
		{
			strcpy(psgName, pseudoPars[i].ValStr.c_str());
		}
		else if(pseudoPars[i].name == "filterDomain")
		{
			if(strcmp(pseudoPars[i].ValStr.c_str(),"T")==0){
			    TorF=0;
			}else if(strcmp(pseudoPars[i].ValStr.c_str(),"F")==0){
			    TorF=1;
			}else{
			    message = "Parameter name " + string(pseudoPars[i].name) + " out of range: [T/F]";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
		else if(pseudoPars[i].name == "filterHp")
		{
			if(strcmp(pseudoPars[i].ValStr.c_str(),"OP")==0){
			    Hp_OForRS = 1;
			}else if(strcmp(pseudoPars[i].ValStr.c_str(),"RS")==0){
			    Hp_OForRS = 0;
			}else{	
			    message = "Parameter name " + string(pseudoPars[i].name) + " out of range: [OP/RS]";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
		else if(pseudoPars[i].name == "filterMp")
		{
			if(strcmp(pseudoPars[i].ValStr.c_str(),"OP")==0){
			    Mp_OForRS = 1;
			}else if(strcmp(pseudoPars[i].ValStr.c_str(),"RS")==0){
			    Mp_OForRS = 0;
			}else{	
			    message = "Parameter name " + string(pseudoPars[i].name) + " out of range: [OP/RS]";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
		else if(pseudoPars[i].name == "filterMs")
		{
			if(strcmp(pseudoPars[i].ValStr.c_str(),"OP")==0){
			    Ms_OForRS = 1;
			}else if(strcmp(pseudoPars[i].ValStr.c_str(),"RS")==0){
			    Ms_OForRS = 0;
			}else{	
			    message = "Parameter name " + string(pseudoPars[i].name) + " out of range: [OP/RS]";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
		else if(pseudoPars[i].name == "filterLp")
		{
			if(strcmp(pseudoPars[i].ValStr.c_str(),"OP")==0){
			    Lp_OForRS = 1;
			}else if(strcmp(pseudoPars[i].ValStr.c_str(),"RS")==0){
			    Lp_OForRS = 0;
			}else{	
			    message = "Parameter name " + string(pseudoPars[i].name) + " out of range: [OP/RS]";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
		else if(pseudoPars[i].name == "filterLs")
		{
			if(strcmp(pseudoPars[i].ValStr.c_str(),"OP")==0){
			    Ls_OForRS = 1;
			}else if(strcmp(pseudoPars[i].ValStr.c_str(),"RS")==0){
			    Ls_OForRS = 0;
			}else{	
			    message = "Parameter name " + string(pseudoPars[i].name) + " out of range: [OP/RS]";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
		else if(pseudoPars[i].name == "LrsT")
		{
			LrsT = pseudoPars[i].ValReal;
		}
		else if(pseudoPars[i].name == "LbT")
		{
			LbT = pseudoPars[i].ValReal;
		}
		else if(pseudoPars[i].name == "nameLog")
		{
			strcpy(nameLog,pseudoPars[i].ValStr.c_str());
		}
		else if(pseudoPars[i].name == "verbosity")
		{
			verbosity = pseudoPars[i].ValInt;
		}
		else if (pseudoPars[i].name == "clobber")
		{
			strcpy(clobberStr, pseudoPars[i].ValStr.c_str());
			if(strcmp(clobberStr,"yes")==0){
			  clobber=1;
			}else{
			  clobber=0;
			}
		}
		// Check if parameter value is in allowed range
		if (pseudoPars[i].type == "int" &&
				(pseudoPars[i].ValInt < pseudoPars[i].minValInt ||
						pseudoPars[i].ValInt > pseudoPars[i].maxValInt))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(pseudoPars[i].minValInt) + "," + boost::lexical_cast<std::string>(pseudoPars[i].maxValInt) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		else if ( pseudoPars[i].type == "double" &&
				(pseudoPars[i].ValReal < pseudoPars[i].minValReal ||
						pseudoPars[i].ValReal > pseudoPars[i].maxValReal))
		{
			message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" +
					boost::lexical_cast<std::string>(pseudoPars[i].minValReal) + "," + boost::lexical_cast<std::string>(pseudoPars[i].maxValReal) + "]";
			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
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
* inDataIteratorTrgFLT function: In production mode, this function takes the optimum number of rows to read
*                                the FILTER extension in the _trg input FITS file and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each row of the optimal filters of the found pulses
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorTrgFLT (long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct)
{
	char val[256];

	int status = EPOK;
	string message = "";
	
	// Declare variables
	double *optfilter, *optfilterin;			// Vector of OPTIMALF column
	double *nrmfctr, *nrmfctrin;				// Vector of NRMFCTR column

	// Allocate GSL vectors and matrix
	gsl_matrix *optfiltergsl_i = gsl_matrix_alloc(nrows,eventsz);
	gsl_vector *row_aux = gsl_vector_alloc(eventsz);
	gsl_vector *nrmfctrgsl_i = gsl_vector_alloc(nrows);				// GSL of input NRMFCTR column

	// Read iterator
	optfilterin = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	optfilter = &optfilterin[1];
	if (toGslMatrix(((void **)&optfilter),&optfiltergsl_i,eventsz,nrows,TDOUBLE,0))
	{
	    message = "Cannot convert OPTIMALF column toGslMatrix";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	nrmfctrin = (double *) fits_iter_get_array(&cols[1]);
	nrmfctr = &nrmfctrin[1];
	if (toGslVector(((void **)&nrmfctr),&nrmfctrgsl_i,nrows,0,TDOUBLE))
	{
	    message = "Cannot convert NRMFCTR column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	// Processing each row of the optimal filters of the found pulses
	for (int i=0; i< nrows; i++)
	{
		gsl_matrix_get_row(row_aux,optfiltergsl_i,i);
		gsl_matrix_set_row(optfiltergslMATRIX,ind,row_aux);
		gsl_vector_set(nrmfctrgsl,ind,gsl_vector_get(nrmfctrgsl_i,i));

		ind ++;
	}

	//Free allocate of GSL vectors
	gsl_matrix_free(optfiltergsl_i);
	gsl_vector_free(row_aux);
	gsl_vector_free(nrmfctrgsl_i);

	return(EPOK);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* inDataIteratorTrg function: This function takes the optimum number of rows to read the input _trg FITS file
*                             and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each found pulse
*	- Initialization of variables
*		- If valid pulse (quality==0)
*			- If calibration mode: 'filter2use' is what has been read from the _flt FITS file
*			- If production mode: 'filter2use' is a row of the matrix where all the optimal filters of the pulses are stored
* 	- If first pulse in the record
*	 	- If no first pulse in the _trg FITS file (state A)
* 		  (The last pulse-free segment in the last record with pulses)
* 			- Calculate RS filter parameters (just in case)
* 		- Iterator for the input FITS file to the chain
* 			- Create structure to run Iteration: inDataIterator
*			- Read columns: Time and ADC
* 			- Called iteration function: inDataIterator
* 		- If there are whole rows without pulses (in xray.fits) before the first pulse (in xray_trg.fits) (state B)
* 			- Calculate RS filter parameters (just in case)
* 			- Iterator for input FITS file
* 				- Create structure to run Iteration DAL: inDataIterator
*				- Read columns: Time and ADC
* 				- Called iteration function: inDataIterator
* 		- Pulse-free segment before the current pulse (state C)
* 			- Calculate RS filter parameters (just in case)
* 		- Current pulse (state D)
* 			- Call RS_filter or findEnergy
* 			- Call writeEnergy
*	- elseif no first pulse in the record
* 		- Previous segment without pulses (state E)
* 			- Calculate RS filter parameters (just in case)
* 		- Current pulse (state F)
* 	 		- Call RS_filter or findEnergy
* 			- Call writeEnergy
* 	- If last pulse in _trg FITS file
* 		- The last pulse-free segment in the current row (state F)
* 			- Calculate RS filter parameters (just in case)
* 		- If there are whole rows without pulses (in xray.fits) after the last pulse
* 			- Calculate RS filter parameters (just in case)
* 			- Iterator for input FITS file
* 				- Create structure to run Iteration: inDataIterator
*				- Read columns: Time and ADC
* 				- Called iteration function: inDataIterator
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorTrg (long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct)
{
	char val2[256];

	int status = EPOK, extver=0;
	string message = "";

	// Declare variables
	double *i0, *i0in;				// Vector of I0 column
	double *tstart, *tstartin;		// Vector of Tstart column
	double *tend, *tendin;			// Vector of Tend column
	short *quality, *qualityin;		// Vector of QUALITY column
	double *dif, *difin;			// Vector of difTstrt column
	short *grade, *gradein;			// Vector of GRADE column

	gsl_vector_view temp;			// In order to handle with gsl_vector_view (subvectors)

	// Allocate input GSL vectors
	pulses = gsl_matrix_alloc(nrows,eventsz);			// GSL of input PULSE column
	qualitygsl = gsl_vector_alloc(nrows);				// GSL of input QUALITY column
	tstartgsl = gsl_vector_alloc(nrows);				// GSL of input Tstart column
	tendgsl = gsl_vector_alloc(nrows);					// GSL of input Tend column
   	difgsl = gsl_vector_alloc(nrows);					// GSL of input difTstrt column
   	gradegsl = gsl_vector_alloc(nrows);					// GSL of input GRADE column

   	energygsl = gsl_vector_alloc(1);					// Pseudoenergy of each pulse

   	gsl_vector *ofMATRIX_row = gsl_vector_alloc(eventsz);

   	int pulseSize;
   	gsl_vector *pulse_aux = gsl_vector_alloc(eventsz);	// GSL of each row

	// Read iterator
   	i0in = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	i0 = &i0in[1];
	if (toGslMatrix(((void **)&i0),&pulses,eventsz,nrows,TDOUBLE,0))
	{
	    message = "Cannot convert I0 column toGslMatrix";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	tstartin = (double *) fits_iter_get_array(&cols[1]);
	tstart = &tstartin[1];
   	if (toGslVector(((void **)&tstart),&tstartgsl,nrows,0,TDOUBLE))
   	{
	    message = "Cannot convert Tstart column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
 	tendin = (double *) fits_iter_get_array(&cols[2]);
	tend = &tendin[1];
	if (toGslVector(((void **)&tend),&tendgsl,nrows,0,TDOUBLE))
	{
	    message = "Cannot convert Tend column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	difin = (double *) fits_iter_get_array(&cols[3]);
	dif = &difin[1];
	if (toGslVector(((void **)&dif),&difgsl,nrows,0,TDOUBLE))
	{
	    message = "Cannot convert difTstrt column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
   	qualityin = (short *) fits_iter_get_array(&cols[4]);
	quality = &qualityin[1];
	if (toGslVector(((void **)&quality),&qualitygsl,nrows,0,TSHORT))
	{
	    message = "Cannot convert QUALITY column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gradein = (short *) fits_iter_get_array(&cols[5]);
	grade   = &gradein[1];
	if (toGslVector(((void **)&grade),&gradegsl,nrows,0,TSHORT))
	{
	    message = "Cannot convert GRADE column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	sprintf(val2,"indataiteratortrg antes de empezar el for %d",nrows);
	strcat(val2,"\n");
	fputs(val2,temporalFile);
		
	// Processing each found pulse
	for (int i=0; i<nrows; i++)
	{
		sprintf(val2,"-------------> Pulse: %d of %d <------------------ ", ntotalrows, eventcnt);
		writeLog(fileRef,"Log", verbosity,string(val2));
		strcat(val2,"\n");
		fputs(val2,temporalFile);

		pulseSize = int((gsl_vector_get(tendgsl,i)-gsl_vector_get(tstartgsl,i))*samprate)-1;
		if (pulseSize == 0)
		{
		    pulseSize = 2;
		}
		pulse = gsl_vector_alloc(pulseSize);

		// Initialize variables
		if (gsl_vector_get(qualitygsl,i) == 0)
		{
		    if (mode == 0)
		    {
		    	filter2use = gsl_vector_alloc(filtergsl->size);
		    	gsl_vector_memcpy(filter2use,filtergsl);
		    }
		    else if (mode == 1)
		    {
		    	filter2use = gsl_vector_alloc(eventsz);
		    	gsl_matrix_get_row(filter2use,optfiltergslMATRIX,totalpulses);
		    	normalizationFactor = gsl_vector_get(nrmfctrgsl,totalpulses);
		    }
		}

		if (gsl_vector_get(difgsl,i) == -1) //First pulse in the record (xray_t.fits)
		{
		    if (tend_p!=-10)  //No first pulse in xray_t.fits
		    {
		    	// The last pulse-free segment in the last event with pulses (state A)
		    	size_trg=int((TIME_p+(eventsz_in/samprate)-tend_p)*samprate);

		    	if(size_trg>0)
		    	{
		    		// No pulse
		    		//cout<<"State A"<<endl;
		    		/*sprintf(val2,"State A");
		    		strcat(val2,"\n");
		    		fputs(val2,temporalFile);*/
		    		gsl_vector *input=gsl_vector_alloc(size_trg);
		    		temp = gsl_vector_subvector(row_input,(tend_p-TIME_p)*samprate,size_trg);
		    		gsl_vector_memcpy(input,&temp.vector);

		    		if (chngplrt == 1)	gsl_vector_scale(input,-1.0);

		    		if (Lb >= input->size)
		    		{
		    			// Sum all the elements of input
		    			if (gsl_vector_Sumsubvector(input,0,input->size,&B))
		    			{
		    				message = "Cannot run routine gsl_vector_Sumsubvector for State A & Lb>=input->size";
		    				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    			}
		    		}
		    		else
		    		{
		    			// Sum the last Lb samples of input
		    			if (gsl_vector_Sumsubvector(input,input->size-Lb,Lb,&B))
		    			{
		    				message = "Cannot run routine gsl_vector_Sumsubvector for State A & Lb<input->size";
		    				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    			}
		    		}

		    		gsl_vector_free(input);
		    	}
		    }

		    iterator_trg++;
		    iterator_in = -1;

		    // Iteration
		    extern int inDataIterator(long totalrows_b,long offset_b,long firstrow_b,long nrows_b,int ncols_b,iteratorCol *cols_b,void *user_strct_b );

		    // Create structure to run Iteration
		    iteratorCol cols_in [2];					// Structure of Iteration
		    int n_cols_b = 2; 							// Number of columns: Time and I0
		    long rows_per_loop_b = 1;					// 1: The iterator will always process one row of the table at a time
		    long offset_b = 0;							// 0: Process all the rows

		    strcpy(extname,"RECORDS");
		    if (fits_movnam_hdu(inObject, ANY_HDU,extname, extver, &status))
		    {
		    	message = "Cannot move to HDU " + string(extname) + " (first pulse in record)";
		    	EP_PRINT_ERROR(message,status); return(EPFAIL);
		    }
		    // Read TIME Column
		    strcpy(straux,"TIME");
		    status = fits_iter_set_by_name(&cols_in[0], inObject, straux, TDOUBLE, InputCol);
		    if	(status)
		    {
		    	message = "Cannot iterate by name column " + string(straux) + " in State A";
		    	EP_PRINT_ERROR(message,status); return(EPFAIL);
		    }
		    // Read ADC Column
		    strcpy(straux,"ADC");
		    status = fits_iter_set_by_name(&cols_in[1], inObject, straux, TDOUBLE, InputCol);
		    if (status)
		    {
		    	message = "Cannot iterate by name column " + string(straux) + " in State A";
		    	EP_PRINT_ERROR(message,status); return(EPFAIL);
		    }

		    // Called iteration function
		    if (fits_iterate_data(n_cols_b, cols_in, offset_b, rows_per_loop_b, inDataIterator,0L,&status))
		    {
		    	message = "Cannot iterate data in columns (first pulse in record)";
		    	EP_PRINT_ERROR(message,status); return(EPFAIL);
		    }

		    while (gsl_vector_get(tstartgsl,i)>(time_input+(eventsz_in/samprate)))
		    {
		    	// There are whole rows without pulses (in xray.fits) before the first pulse (in xray_t.fits) (state B)

		    	// No pulse
		    	//cout<<"State B"<<endl;
		    	/*sprintf(val2,"State B %d",i);
		    	strcat(val2,"\n");
		    	fputs(val2,temporalFile);*/

		    	gsl_vector *input = gsl_vector_alloc(eventsz_in);
		    	gsl_vector_memcpy(input,row_input);

		    	if (chngplrt == 1)	gsl_vector_scale(input,-1.0);

		    	if (Lb >= input->size)
		    	{
		    		// Sum all the elements of input
		    		if (gsl_vector_Sumsubvector(input,0,input->size,&B))
		    		{
		    			message = " Cannot run gsl_vector_Sumsubvector for state B & Lb>=input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}
		    	else
		    	{
		    		// Sum the last Lb samples of input
		    		if (gsl_vector_Sumsubvector(input,input->size-Lb,Lb,&B))
		    		{
		    			message = " Cannot run gsl_vector_Sumsubvector for state B & Lb<input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}

		    	gsl_vector_free(input);

		    	iterator_trg++;
		    	iterator_in = -1;

		    	// Iteration
		    	extern int inDataIterator(long totalrows_b,long offset_b,long firstrow_b,long nrows_b,int ncols_b,iteratorCol *cols_b,void *user_strct_b );

				// Create structure to run Iteration
		    	iteratorCol cols_in [2];					// Structure of Iteration
		    	int n_cols_b = 2; 							// Number of columns: Time and ADC
		    	long rows_per_loop_b = 1;					// 1: The iterator will always process one row of the table at a time
		    	long offset_b = 0;							// 0: Process all the rows

				// Read TIME Column
		    	strcpy(straux,"TIME");
		    	status = fits_iter_set_by_name(&cols_in[0], inObject, straux, TDOUBLE, InputCol);
		    	if (status)
		    	{
		    		message = "Cannot iterate by name column " + string(straux) + " in State B";
		    		EP_PRINT_ERROR(message,status); return(EPFAIL);
		    	}
				// Read ADC Column
		    	strcpy(straux,"ADC");
		    	status = fits_iter_set_by_name(&cols_in[1], inObject, straux, TDOUBLE, InputCol);
		    	if (status)
		    	{
		    		message = "Cannot iterate by name column " + string(straux) + " in State B";
		    		EP_PRINT_ERROR(message,status); return(EPFAIL);
		    	}

				// Called iteration function
		    	if (fits_iterate_data(n_cols_b, cols_in, offset_b, rows_per_loop_b, inDataIterator,0L,&status))
		    	{
		    		message = "Cannot iterate data in columns in State B";
		    		EP_PRINT_ERROR(message,status); return(EPFAIL);
		    	}
		    }

		    // Pulse-free segment before the current pulse (state C)
		    size_trg=int((gsl_vector_get(tstartgsl,i)-time_input)*samprate);

		    if(size_trg>0)
		    {
		    	//No pulse
		    	//cout<<"State C"<<endl;
		    	/*sprintf(val2,"State C");
		    	strcat(val2,"\n");
		    	fputs(val2,temporalFile);*/
		    	gsl_vector *input=gsl_vector_alloc(size_trg);
		    	temp = gsl_vector_subvector(row_input,0,size_trg);
		    	gsl_vector_memcpy(input,&temp.vector);

		    	if (chngplrt == 1)	gsl_vector_scale(input,-1.0);

		    	if (Lb >= input->size)
		    	{
		    		// Sum all the elements of input
		    		if (gsl_vector_Sumsubvector(input,0,input->size,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector in State C when Lb>=input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}
		    	else
		    	{
		    		// Sum the last Lb samples of input
		    		if (gsl_vector_Sumsubvector(input,input->size-Lb,Lb,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector in State C when Lb<input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}

		    	gsl_vector_free(input);
		    }

		    // Current pulse (state D)
		    gsl_vector_set(energygsl,0,0);

		    if (gsl_vector_get(qualitygsl,i) == 0)
		    {
			//if((gsl_vector_get(tstartgsl,i)>Tinitial || Tinitial==-1) && (gsl_vector_get(tendgsl,i)<Tfinal || Tfinal==-1))
			//{
			    gsl_matrix_get_row(pulse_aux,pulses,i);
			    //cout<<"State D"<<endl;
			    /*sprintf(val2,"State D");
			    strcat(val2,"\n");
			    fputs(val2,temporalFile);*/
			    temp = gsl_vector_subvector(pulse_aux,0,pulseSize);
			    gsl_vector_memcpy(pulse,&temp.vector);

			    //Pulse
			    if ((gsl_vector_get(gradegsl,i) == 1) && (Hp_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state D for 1 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 1) && (Hp_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state D for 1 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 21) && (Mp_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state D for 21 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 21) && (Mp_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state D for 21 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 22) && (Ms_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state D for 22 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 22) && (Ms_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state D for 22 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 31) && (Lp_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state D for 31 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 31) && (Lp_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state D for 31 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 32) && (Ls_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state D for 32 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 32) && (Ls_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy(pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state D for 32 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    gsl_vector_set(energygsl,0,pseudoenergy);

			    gsl_vector_free(pulse);
			//}
		    }

		    // Call writeEnergy
		    if (writeEnergy())
		    {
		    	message = "Cannot run routine writeEnergy for state D";
		    	EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    }
		    totalpulses++;
		    totalsegments++;
		    ntotalrows++;
		}

		if(gsl_vector_get(difgsl,i) != -1)	//No first pulse in the record (xray_t.fits)
		{
		    // Previous segment without pulses (state E)
		    size_trg=int((gsl_vector_get(tstartgsl,i)-tend_p)*samprate);

		    if(size_trg>0)
		    {
		    	//No pulse
		    	//cout<<"State E"<<endl;
		    	/*sprintf(val2,"State E");
		    	strcat(val2,"\n");
		    	fputs(val2,temporalFile);*/
		    	gsl_vector *input=gsl_vector_alloc(size_trg);
		    	temp = gsl_vector_subvector(row_input,(tend_p-time_input)*samprate,size_trg);
		    	gsl_vector_memcpy(input,&temp.vector);

		    	if (chngplrt == 1)	gsl_vector_scale(input,-1.0);

		    	if (Lb >= input->size)
		    	{
		    		// Sum all the elements of input
		    		if (gsl_vector_Sumsubvector(input,0,input->size,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector for State E when Lb>=input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}
		    	else
		    	{
		    		// Sum the last Lb samples of input
		    		if (gsl_vector_Sumsubvector(input,input->size-Lb,Lb,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector for State E Lb<input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}

		    	gsl_vector_free(input);
		    }

		    // Current pulse (state F)
		    gsl_vector_set(energygsl,0,0);

		    //if((gsl_vector_get(tstartgsl,i)>Tinitial || Tinitial==-1) && (gsl_vector_get(tendgsl,i)<Tfinal || Tfinal==-1))
		    //{
			//cout<<"State F"<<endl;
			/*sprintf(val2,"State F");
			strcat(val2,"\n");
			fputs(val2,temporalFile);*/
			gsl_matrix_get_row(pulse_aux,pulses,i);
			temp = gsl_vector_subvector(pulse_aux,0,pulseSize);
			gsl_vector_memcpy(pulse,&temp.vector);

			if (gsl_vector_get(qualitygsl,i) == 0)
			{
			    if ((gsl_vector_get(gradegsl,i) == 1) && (Hp_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state F for 1 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 1) && (Hp_OForRS == 1))
			    {
			    	// Call findEnergy
				if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
				{
					message = "Cannot run routine findEnergy in state F for 1 & 1 conditions";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 21) && (Mp_OForRS == 0))
			    {
				if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
				{
					message = "Cannot run routine RS_filter in state F for 21 & 0 conditions";
					EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
				}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 21) && (Mp_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state F for 21 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 22) && (Ms_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state F for 22 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 22) && (Ms_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state F for 22 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 31) && (Lp_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state F for 31 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 31) && (Lp_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state F for 31 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 32) && (Ls_OForRS == 0))
			    {
			    	if (RS_filter (pulse, Lrs, Lb, B, &pseudoenergy))
			    	{
			    		message = "Cannot run routine RS_filter in state F for 32 & 0 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    else if ((gsl_vector_get(gradegsl,i) == 32) && (Ls_OForRS == 1))
			    {
			    	// Call findEnergy
			    	if (findEnergy (pulse,filter2use,TorF,normalizationFactor,&pseudoenergy))
			    	{
			    		message = "Cannot run routine findEnergy in state F for 32 & 1 conditions";
			    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    	}
			    }
			    gsl_vector_set(energygsl,0,pseudoenergy);
			}
			gsl_vector_free(pulse);
		    //}

		    //Call writeEnergy
		    if (writeEnergy())
		    {
		    	message = "Cannot run routine writeEnergy for state F";
		    	EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    }
		    totalpulses++;
		    totalsegments++;
		    ntotalrows++;
		}

		tend_p = gsl_vector_get(tendgsl,i);
		TIME_p = time_input;

		if ((ntotalrows-1)==eventcnt)	// Last pulse in _trg.fits FITS file
		{
		    // The last pulse-free segment in the current row (state G)
		    size_trg=int((time_input+(eventsz_in/samprate)-gsl_vector_get(tendgsl,i))*samprate);
		    if (size_trg>0)
		    {
		    	//No pulse
		    	//cout<<"State G"<<endl;
		    	/*sprintf(val2,"State G");
		    	strcat(val2,"\n");
		    	fputs(val2,temporalFile);*/
		    	gsl_vector *input=gsl_vector_alloc(size_trg);
		    	temp = gsl_vector_subvector(row_input,(gsl_vector_get(tendgsl,i)-time_input)*samprate,size_trg);
		    	gsl_vector_memcpy(input,&temp.vector);

		    	if (chngplrt == 1)	gsl_vector_scale(input,-1.0);

		    	if (Lb >= input->size)
		    	{
		    		// Sum all the elements of input
		    		if (gsl_vector_Sumsubvector(input,0,input->size,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector for State G and Lb >= input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}
		    	else
		    	{
		    		// Sum the last Lb samples of input
		    		if (gsl_vector_Sumsubvector(input,input->size-Lb,Lb,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector for State G and Lb < input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}
		    	gsl_vector_free(input);
		    }

		    while((iterator_trg)!=eventcnt_in-1)
		    {
		    	// There are whole rows without pulses (in xray.fits) after the last pulse (in xray_t.fits) (state H)
		    	//No pulse
		    	//cout<<"State H"<<endl;
		    	/*sprintf(val2,"State H");
		    	strcat(val2,"\n");
		    	fputs(val2,temporalFile);*/
		    	gsl_vector *input = gsl_vector_alloc(eventsz_in);
		    	gsl_vector_memcpy(input,row_input);

		    	if (chngplrt == 1)	gsl_vector_scale(input,-1.0);

		    	if (Lb >= input->size)
		    	{
		    		// Sum all the elements of input
		    		if (gsl_vector_Sumsubvector(input,0,input->size,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector for State H and Lb >= input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}
		    	else
		    	{
		    		// Sum the last Lb samples of input
		    		if (gsl_vector_Sumsubvector(input,input->size-Lb,Lb,&B))
		    		{
		    			message = "Cannot run routine gsl_vector_Sumsubvector for State H and Lb < input->size";
		    			EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		    		}
		    	}

		    	gsl_vector_free(input);

		    	iterator_trg++;
		    	iterator_in = -1;

		    	// Iteration
		    	extern int inDataIterator(long totalrows_b,long offset_b,long firstrow_b,long nrows_b,int ncols_b,iteratorCol *cols_b,void *user_strct_b );

			// Create structure to run Iteration
		    	iteratorCol cols_in [2];					// Structure of Iteration
		    	int n_cols_b = 2; 							// Number of columns: TIME and ADC
		    	long rows_per_loop_b = 1;					// 1: The iterator will always process one row of the table at a time
		    	long offset_b = 0;							// 0: Process all the rows

			// Read TIME Column
		    	strcpy(straux,"TIME");
		    	status = fits_iter_set_by_name(&cols_in[0], inObject, straux, TDOUBLE, InputCol);
		    	if (status)
		    	{
		    		message = "Cannot iterate by name column " + string(straux) + " in State H";
		    		EP_PRINT_ERROR(message,status); return(EPFAIL);
		    	}
			// Read V Column
		    	strcpy(straux,"ADC");
		    	status = fits_iter_set_by_name(&cols_in[1], inObject, straux, TDOUBLE, InputCol);
		    	if (status)
		    	{
		    		message = "Cannot iterate by name column " + string(straux) + " in State H";
		    		EP_PRINT_ERROR(message,status); return(EPFAIL);
		    	}

			// Called iteration function
		    	if (fits_iterate_data(n_cols_b, cols_in,offset_b, rows_per_loop_b, inDataIterator,0L,&status))
		    	{
		    		message = "Cannot iterate data in columns (first pulse in record)";
		    		EP_PRINT_ERROR(message,status); return(EPFAIL);
		    	}
		    }
		}
	}
	

	// Free allocate of GSL vectors
	gsl_vector_free (pulse_aux);
	gsl_matrix_free (pulses);
	gsl_vector_free (qualitygsl);
	gsl_vector_free (tstartgsl);
	gsl_vector_free (tendgsl);
	gsl_vector_free (difgsl);
	gsl_vector_free (gradegsl);
	gsl_vector_free (energygsl);

	return(EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* inDataIteratorFlt function: In calibration mode, this function takes the optimum number of rows to read
*                             the _flt input FITS file and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each row of the _flt.fits file
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorFlt (long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct)
{
	int status = EPOK;
	string message = "";

	// Declare variables
	double *filter, *filterin;			// Vector of filter column

	// Allocate input GSL vectors
	gsl_vector *filtergsl_i=gsl_vector_alloc(nrows);

	// Read iterator
	filterin = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	filter = &filterin[1];

	if (toGslVector(((void **)&filter),&filtergsl_i,nrows,0,TDOUBLE))
	{
	    message = "Cannot convert Filter column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Processing each row of the _flt.fits file
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set(filtergsl,i_filter,gsl_vector_get(filtergsl_i,i));

		i_filter ++;
	}

	//Free allocate of GSL vectors
	gsl_vector_free(filtergsl_i);

	return(EPOK);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* inDataIterator function: This function takes the optimum number of rows to read the input FITS file to the chain
*                 and works iteratively
*
* - Declare variables
* - Allocate input GSL vectors
* - Read iterator
* - Processing each record of pulses
* - Free allocate of GSL vectors
*****************************************/
int inDataIterator (long totalrows_b, long offset_b, long firstrow_b, long nrows,int ncols_b, iteratorCol *cols_b, void *user_strct_b)
{
	int status = EPOK;
	string message = "";

	// Declare variables
	double *time, *timein;
	double *v, *vin;

	double time_input_aux;
	gsl_vector *row_input_aux = gsl_vector_alloc(eventsz_in);

	// Allocate input GSL vectors
	timegsl = gsl_vector_alloc(1);					// GSL of input TIME column
	eventsgsl = gsl_matrix_alloc(1,eventsz_in);		// GSL of input I0 column

	// Read iterator
	timein = (double *) fits_iter_get_array(&cols_b[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	time = &timein[1];
	if (toGslVector(((void **)&time),&timegsl,1,0,TDOUBLE))
	{
	    message = "Cannot convert TIME column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	vin = (double *) fits_iter_get_array(&cols_b[1]);
	v = &vin[1];
	if (toGslMatrix(((void **)&v),&eventsgsl,eventsz_in,1,TDOUBLE,0))
	{
	    message = "Cannot convert I0 column toGslMatrix";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Processing each record of pulses
	for (int i=0; i< nrows; i++)
	{
		time_input_aux=gsl_vector_get(timegsl,i);
		gsl_matrix_get_row(row_input_aux,eventsgsl,0);
		iterator_in++;
	}

	if(iterator_in==iterator_trg)
	{
		time_input=time_input_aux;
		gsl_vector_memcpy(row_input,row_input_aux);
	}

	// Free allocate of GSL vectors
	gsl_vector_free(timegsl);
	gsl_vector_free(row_input_aux);
	gsl_matrix_free(eventsgsl);

	return(EPOK);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 8 ************************************************************
* findEnergy function: This function obtains the pseudoenergy of each pulse in time or frequency domain
*
* Time domain: E = int(0,inf){p(t)f(t)df} when we assume that the pulse and the template are well aligned
* 			   p(t) -> Pulse
*              f(t) -> Optimal filter
* Frequency domain: E = int(-inf,inf){|P(f)|^2df} = int(-inf,inf){|P(f)||F(f)|df}
*                   P(f) -> FFT of the pulse
*                   F(f) -> FFT of the optimal filter (matched to the pulse => P(f)~F(f)))
*                   |P(f)|^2 ~ |P(f)||F(f)| => E = int(-inf,inf){|P(f)||F(f)|df}
*
* - If TIME DOMAIN filtering (domain = 0):
*  - Declare variables
*  - Calculate the pseudoenergy by calculating the appropriate integral
*  - Apply the normalization factor
* - If FREQUENCY DOMAIN filtering (domain = 1):
*   - Declare variables
*   - FFT of pulse
*   - FFT of filter
*   - Multiply FFT's
*   - Abs
*   - Pseudoenergy value: Sum of the bins of the FFTabs (integral)
*   - Apply the normalization factor
*   - Free allocate of GSL vectors
*
* Parameters:
* - vector: Pulse
* - filter: Optimal filter
* - domain: Time domain or frequency domain to filter
* - pseudoenergy: Calculated pseudoenenergy
****************************************************************************/
int findEnergy (gsl_vector *vector, gsl_vector *filter, int domain, double nrmfctr, double *pseudoenergy)
{
	char val1[256];

	int status = EPOK;
	string message = "";

	if (domain == 0)	// Time domain filtering
	{
		int low_size = min(vector->size,filter->size);
		double sum = 0;

		for (int i=0; i<low_size; i++)
		{
			sum = sum + gsl_vector_get(vector,i)*gsl_vector_get(filter,i);
		}

		*pseudoenergy = sum;
		*pseudoenergy = *pseudoenergy/nrmfctr;
	}
	else if (domain == 1)	// Frequency domain filtering (multiply vectorFFT and filterFFT)
    {
    	// Declare variables
    	gsl_vector_complex *vectorFFT = gsl_vector_complex_alloc(vector->size);
    	gsl_vector_complex *filterFFT = gsl_vector_complex_alloc(filter->size);

    	long low_size, large_size;
    	if (vector->size < filter->size)
    	{
    		low_size = vector->size;
    		large_size = filter->size;
    	}
    	else if (vector->size >= filter->size)
    	{
    		low_size = filter->size;
    		large_size = vector->size;
    	}
    	double SelectedTimeDuration = large_size/samprate;	// In order to have the minimum frequency step (1/SelectedTimeDuration)

    	gsl_vector_complex *prodFFT = gsl_vector_complex_alloc(large_size);
    	gsl_vector *prodFFT_R = gsl_vector_alloc(large_size);

    	// FFT calculus (vectorFFT)
    	if (FFT(vector,vectorFFT,SelectedTimeDuration))
    	{
    		message = "Cannot run routine FFT";
    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
    	}
    	gsl_vector_complex_scaleIFCA(vectorFFT,gsl_complex_rect(sqrt(vector->size),0.0));
    	gsl_vector_complex_scaleIFCA(vectorFFT,gsl_complex_rect(1.0/sqrt(large_size),0.0)); //Factor 1/sqrt(N), in the FFT expression

    	// FFT calculus (filterFFT)
    	if (FFT(filter,filterFFT,SelectedTimeDuration))
    	{
    		message = "Cannot run routine FFT when domain=1)";
    		EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
    	}
    	gsl_vector_complex_scaleIFCA(filterFFT,gsl_complex_rect(sqrt(filter->size),0.0));
    	gsl_vector_complex_scaleIFCA(filterFFT,gsl_complex_rect(1.0/sqrt(large_size),0.0)); //Factor 1/sqrt(N), in the FFT expression

    	// Frequency domain filtering => Multiply vectorFFT and filterFFT
    	for (int i=0;i<prodFFT->size;i++)
    	{
    		gsl_vector_complex_set(prodFFT,i,gsl_complex_rect(0.0,0.0));
    	}
    	gsl_vector *vectorFFT_R = gsl_vector_alloc(large_size);
    	gsl_vector *filterFFT_R = gsl_vector_alloc(large_size);
    	gsl_vector_complex_absRIFCA(vectorFFT_R,vectorFFT);
    	gsl_vector_complex_absRIFCA(filterFFT_R,filterFFT);

    	for (int i=0;i<low_size;i++)
    	{
    		gsl_vector_set(prodFFT_R,i,gsl_vector_get(vectorFFT_R,i)*gsl_vector_get(filterFFT_R,i));
    	}
    	*pseudoenergy = 0.0;
    	for (int i=0;i<prodFFT_R->size/2;i++)
    	{
    		*pseudoenergy = *pseudoenergy + gsl_vector_get(prodFFT_R,i)*(1/SelectedTimeDuration);
    	}
    	*pseudoenergy = *pseudoenergy/nrmfctr;

    	// Free allocate of GSL vectors
    	gsl_vector_complex_free (vectorFFT);
    	gsl_vector_complex_free (filterFFT);
    	gsl_vector_complex_free (prodFFT);
    	gsl_vector_free (prodFFT_R);
    	gsl_vector_free (vectorFFT_R);
    	gsl_vector_free (filterFFT_R);
    }

    return EPOK;
}
/*xxxx end of SECTION 8 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 9 ************************************************************
* RS_filter function: This function uses the running sum filter to find the pseudoenergy.
*            It always works in time domain.
****************************************/
int RS_filter (gsl_vector *vector, long Lrs, long Lb, double B, double *pseudoenergy)
{
	char val2[256];

	int status = EPOK;
	string message = "";

	// Declare variables
	double Rs;
	double Rs_max = -1e10;
	double Bp;

	for (int i=0;i<vector->size;i++)
	{
		if (i+Lrs > vector->size) break;

		if (gsl_vector_Sumsubvector(vector,i,Lrs,&Rs))
		{
		    message = "Cannot run routine gsl_vector_Sumsubvector for iteration i=" + boost::lexical_cast<std::string>(i);
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}

		if (Rs > Rs_max)	  Rs_max = Rs;
	}

	Bp = B*Lrs/Lb;
	*pseudoenergy = (Rs_max-Bp)/Lrs;

	return EPOK;
}
/*xxxx end of SECTION 9 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 10 ***************************************
* writeEnergy function: This function adds information (Pseudoenergy column) to the "PSGRADE" extension of _psg input/output FITS file
***********************************************/
int writeEnergy ()
{
	char val2[256];

	string message = "";
	int status = EPOK;

	// Add Pseudoenergy Column to PSG FITS
	IOData obj;
	
	obj.inObject = psgObject;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"PSGRADE");
	obj.iniRow = totalpulses+1;
	obj.endRow = totalpulses+1;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	strcpy(obj.nameCol,"Pseudoenergy");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit," ");
	
	if (writeFitsSimple (obj,energygsl))
	{
	    message = "Cannot run routine writeFitsSimple for " + string(obj.nameCol);
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

 	// Free allocate of memory
 	delete [] obj.nameTable;
 	delete [] obj.nameCol;
 	delete [] obj.unit;

 	return EPOK;
}
/*xxxx end of SECTION 10 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 11 ************************************************************
* writeKeywords function: This function writes keywords in "PSGRADE" extension of _psg input/output FITS file
****************************************/
int writeKeywords ()
{
	int status = EPOK, extver=0;
	string message = "";
	
	strcpy(extname,"PSGRADE");
	if (fits_movnam_hdu(psgObject, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to HDU " + string(extname);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

 	// Set PROCESS keyword
	char str_verb[125];			sprintf(str_verb,"%d",verbosity);
	char str_TorF[125];			sprintf(str_TorF,"%d",TorF);
	char str_Hp_OForRS[125];	sprintf(str_Hp_OForRS,"%d",Hp_OForRS);
	char str_Mp_OForRS[125];	sprintf(str_Mp_OForRS,"%d",Mp_OForRS);
	char str_Ms_OForRS[125];	sprintf(str_Ms_OForRS,"%d",Ms_OForRS);
	char str_Lp_OForRS[125];	sprintf(str_Lp_OForRS,"%d",Lp_OForRS);
	char str_Ls_OForRS[125];	sprintf(str_Ls_OForRS,"%d",Ls_OForRS);
	char str_Lrs[125];			sprintf(str_Lrs,"%f",Lrs);
	char str_Lb[125];			sprintf(str_Lb,"%f",Lb);

	string processout (string("PSEUDOENERGY") 	+ ' ' +
	string(trgName) 		+ ' ' + string(fltName)      + ' ' + string(psgName)      + ' ' + string(inName) + ' ' +
	string(str_TorF) 	    + ' ' +
	string(str_Hp_OForRS) 	+ ' ' + string(str_Mp_OForRS)+ ' ' + string(str_Ms_OForRS)+ ' ' +
	string(str_Lp_OForRS)   + ' ' + string(str_Ls_OForRS)+ ' ' +
	string(str_Lrs)         + ' ' + string(str_Lb)       + ' ' +
	string(nameLog)		    + ' ' + string(str_verb)     + ' ' + string(clobberStr) + ' ' +
	string("(")				+      (string) create       + string(")"));

	strcpy(keyname,"PROC1");
	strcpy(keyvalstr,processout.c_str());
	if (fits_write_key_longwarn(psgObject,&status))
	{
	    message = "Cannot update key " + string(keyname) + " in " + string(psgName);
	    EP_PRINT_ERROR(message,status); return(EPFAIL);
	}
	if (fits_update_key_longstr(psgObject,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(psgName);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	return EPOK;
}
/*xxxx end of SECTION 11 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
