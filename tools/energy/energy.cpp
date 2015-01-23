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
*                 ENERGY
*
*  File:      energy.cpp
*  Version:   12.0.0
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
* version 1.0.0		19/10/06    First version
* version 1.0.1		29/01/06	Adding new function cutFunction
* version 1.0.2   	16/02/06	Change value of binsize for scale of Histogram to 200
* version 1.0.3   	08/03/07	Change the name of Group Extention of input Fits file
* 							    for "EUR-DMX"
* 							    Change the name of Group Extention of output Fits file
* 							    for "EUR-HOLZGAUSS"
* version 1.0.4 	20/04/07 	Delete environment value "FITS_DIR", it's incompatible with the browse of param_gui.
* 							   	"FIT_DIR" imposes a input directory wich contains the input fits file (That is incompatible)
* 							   	with the browse).
*  							   	Change input Name of the fits file. Now, the Name include the .fits extension. Now, the file
* 							   	Name is compatible with the browse of param_gui.
* 							   	Change in the Call to "readFitsComplex", "readFitsSimple", "writeFitsComplex" and
* 								"writeFitsSimple". Now is used the struct "OIData". See documentation of Utils module. (inoutUtils).
* version 1.0.5		24/04/07 	Include new input parameter in the module. Adding in .par file "output File". Now,
* 								the user can choose the name of the output .fits file.
* 								-jjs-
* 								Change maximum number of data channels in an input file (4)
*  								Parameterise the names for the input ADC channels, and
*          					 	the names for the result data column names.
*  								If inputfile is NOT modulated, do not demodulate, but do all else, eg. scaling the data
* 								ADC channel names are now parameterised
* version 1.1.0		23/05/07 	Modification of the module using the new module inoutFits included in Utils. Now, Reads and Writes
* 								of Fits use block of rows.
* 								Including coments using RIL logs.
* version 2.0.0		12/02/08	Holzgauss reads from ers fits file 2.0.0
* version 2.1.0		29/05/08	Adding input/output parameters "ql" and "nBins".
* 								Optimezed read and write of the FITS files.
* 								Adding output keyword: The result of the fits "Scale Factor, Offset of Energy and fwhm"
* 								Adding output keyword: "ENERGYRESOL"
* 								Adding column in the hzg FITS file "Histo" and "HistoFit"
* 								Change fuction fitGauss by getParams
* 								Adding fuction fit_poly
* version 3.0.0		10/09/08	Adding input parametrer "nameLog"
* 								Adding input parametrer "verbosity" and processing it using RIL libraries.
* 								Adding new function "readInputKeywords"
* 								Deleting in/output keyword "channelcount"
* 								Including error codes.
* version 3.1.0 	29/09/08	Included "CREATE" and "PROCESS" keyword in the output FITS file
* 								Include error processing of "INPUT_PARAMETER_NON_DEFINE"
* version 3.2.0		24/09/08	If the module fail. It will be warning in the output.
* 								Included parameter time to measure.
* version 3.3.0		19/11/08	Resolved bug.
* version 4.0.0		14/01/09	Included new input parameter: deltaBase
* 								Included new input keywords
* version 4.1.0		22/01/09	Modificated error code of input parameters
* version 5.1.0		02/02/09	Include new functinality: Convert the energies of each pulses in eV.
* 								Using linear or non lineal method to obtain the coeficient.
* 								Its use energyCalibration function
* 								Implement gaussian fit for PULSESRC=MONOCHROMATIC
* 								Remove the function "cutHisto"
* 								Included new input keywords to MONOCHROMATIC and PSEUDO modes
* 								Included new input keywords "ec". Energy Calibration.
* 								Included new Warning: Non Valid Gauss Fit and new codes of quality columns.
* 								Included new Warning: Non Valid Holzer Fit and new codes of quality columns.
* 								Changed name of output keyword "CREATE" to "CREATOR"
* version 5.1.1		04/02/09	Resolved bug of gauss and holzer fit (nan values).
* 								Used sigma of baseline to obtain chisq of gauss and holzer fit.
* 								Resolved bug in createHisto.
* 								Checked free memory allocate.
* version 6.0.0		20/02/09	Included, modificated and removed input keywords
* version 6.1.0		24/02/09	Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
* 								Solved bug 4: stops on zero value for the keyword TIMEZERO, error -62436.
* 								Included FTYPE keywords in output FITS file.
* 								Rename output keywords: NFACTOR, EFWHM and OFFSET.
* version 6.2.0		12/03/09	Resolved bug in EnergyCalibration function.
* version 7.0.0                 Used the keywords libraries
* version 8.0.0     17/04/09    The input parameter fitabstol is used
* 							    Delete "PILInit(argc, argv);" and "PILClose(PIL_OK);"
*                               Documentation updated
*                               Added Emin and Emax and change deltaBase to deltaKbeta (input parameters)
* 								Added two new functions: checkEnergyRange and checkNullBars
* 								Changed gsl_vector_set xhistogsl in createHisto (+binSize/2 added)
* version 9.0.0		03/06/09    Added Tmin and Tmax (input parameters)
* 								Added analyzeDrift, fitPseudo, fcnChisq, gaussShort and trapezeRule functions
*                               Modified getParams function
* 								polyFit not used in getParams
* version 10.0.0    23/06/09	Added pCentdown and pCentup input parameters (instead pCent)
*                               FEK treatment: fitPseudo and fcnChisq modified
*                                              lorenzian, humlicek_v12, voigt and holzer added
*                               Deleted non-used functions
* 								Changed order in delete's
* version 11.0.0	xx/xx/09	If error file exists, its contests will be deleted
*                   30/11/09	UEFWHM really written in the HZG output FITS file
*								The PROCESS keyword written in the modified PSH file keeps the PROCESS keyword read from the
*                               PSH input FITS file and the PROCESS generated in HOLZGAUSS and the HOLZGAUSS version are added
*                               66801 applied also when the histogram peak is just in the beginning.
*                   29/09/10    Restructuring the columns of the input FITS files from PULSESHAPE:
*                               	Old: TIME-Tstart-Tend-Baseline-Sigma-MaxTime-MaxCurrent-TRise-Nrise-TFall-Nfall-...
*                                        ...-MaxTimeFit-MaxCurrentFit-TRiseFit-NriseFit-TFallFit-NfallFit-ChisqFit-...
*                                        ...-BaselineFit-Quality
*                                   New: Tstart-Tend-MaxTime-MaxCurrent-Quality (EUR-PSH extension)
*                                   	 Tstart-Tend-Pseudoenergy-Value (EUR-ENR extension)
*                   25/02/11	Analyzing pulses-free intervals (not only pulses)
*                               Functions: "checkNullBars", "getParams", "fitPseudo" changed
*                               Not necessary to read "Tend", "MaxTime", "MaxCurrent", "Quality" columns from xray_psh.fits
*								The energy of interval without pulses (not only pulses) is written in EUR-ENR extension in PULSESHAPE FITS --> "Energy" column
*								Renamed EUR-HZG extension (now EUR-HZGPULSES) in xray_hzg.fits
*								New extension (EUR-HZGNONPULSES) in xray_hzg.fits
*					28/02/11    PROCESS also includes HOLZGAUSS version (CREATOR)
*					11/03/11    "valuegsl2" added
*					25/03/11    Quick Look mode no longer used
*					25/03/11    Handling an empty EUR-PSH extension (no pulses)
*					24/08/11    'gsl_vector_view temp;' added to handle with subvectors
*					09/11/11	_psh input FITS file has a new column QUALITY in EUR-ENR => It is read
*					            Instead of checking if energygsl is different from 0, quality=0 is checked:
*    					            Pulse-free intervals will always have quality=0 => Always taken into account
*	                                Pulses whose quality is different from 0 are rejected
*	                14/11/11    New input parameters: 'b_cF' and 'c_cF' (calibration factors)
*	                            If calibration factor are provided, they are not going to be calculated in the task
*	                09/07/14    PIL/RIL/Common dependencies removed. New functions used to read input parameters
*	                10/07/14    Pulse-free intervals no longer used
*	                24/07/14	Histogram without null bars used in the calculus BUT histogram with null bars
*	                            in the output FITS file
*	                31/07/14    Code modified in order to work with non monochromatic input FITS files
*	                01/08/14    c (input parameter) is allowed to be < 0
* version 12.0.0	  Sep/14    Analytical calculus of b and c and sigma (no fitting)
*                     Dec/14    DAL->CFITSIO migration + error routines update
*                   09/01/15	'short' instead of 'int' when QUALITY and TAILBFR columns are read (in 'inDataIteratorA' and 'inDataIteratorB')
*                   12/01/15    'PROCESS' input keywords from _psh input FITS files are not read
*					            'PROCESS' output keywords in _psh output FITS files divided in three ones ('PROC0' generated by PULSEGRADE, 'PROC1'
*					            generated by PSEUDOENERGY and 'PROC2' generated by ENERGY)
*					            'ANNALS' input keywords from _psh input FITS files are not read
*					            'ANNALS' output keywords in _psh output FITS files divided in two ones ('MOD0' generated by PSEUDOENERGY and 'MOD1'
*            					generated by ENERGY)
* 		            21/01/15	Added clobbering
*           			        Renaming of extensions: EUR-PSH -> PSGRADE
*                               Renaming of parameters variables (from holzgauss to energy)
* 	                22/01/15    Renaming of input parameter HorMorL --> gradesForER
* **********************************************************************************/

/******************************************************************************
DESCRIPTION:

The ENERGY task has two operation modes, calibration and production.

In calibration mode, the goal is calculate the linear and quadratic calibration factors to transform pseudoenergies
(a.u.), e, into energies (eV), E, by solving according to Cramer a system of analytical equations. There must be two
monochromatic input FITS files. Then, pseudoenergies are transformed into energies. Moreover, energy resolution
information, given by the standard deviation or the full width at half maximum, is calculated.

	E = be +ce^2                                |
	Monochromatic pulses                        | => SUM(i=1,nx)[bxi+cxi²]=nx·E0x
	2 input FITS files:                         |    SUM(j=1,ny)[byj+cyj²]=ny·E0y
  	  inFileA: nx values of pseudoenergies xi	|
           	   Energy = E0x                     |    SUM(i=1,nx)[(bxi+cxi²-E0x)²]=(nx-1)·sigmax²
  	  inFileB: ny values of pseudoenergies yj   |    SUM(j=1,ny)[(byj+cyj²-E0y)²]=(ny-1)·sigmay²
           	   Energy = E0y                     |
	E[(X-media)²]=sigma²                        |

In production mode, there is only one input FITS file (non monochromatic) and the linear and quadratic calibration factors
are provided by the user as input parameters. Then, pseudoenergies are transformed into energies.

The user must supply the following input parameters:

 - optmode: Calibration (0) or production (1) mode
 - inFileA: Name of the _psg FITS file which contains the info of the pulses pseudoenergies
 - inFileB: Name of the _psg FITS file which contains the info of the pulses pseudoenergies (only in calibration mode)
 - b_cF: Linear calibration factor to transform from pseudoenergies into energies (only in production mode)
 - c_cF: Quadratic calibration factor to transform from pseudoenergies into energies (only in production mode)
 - gradesForER: Calculate energy resolution for None(N) or Hp(H) or Mp+Ms(M) or Lp+Ls(3) or Hp+Mp+Ms+Lp+Ls(A) events
 - namelog: Output log file name
 - verbosity: Verbosity level of the output log file

 MAP OF SECTIONS IN THIS FILE:

 - 1. INCLUDE's
 - 2. MAIN
 - 3. initModule
 - 4. inDataIteratorA
 - 5. inDataIteratorB
 - 6. calculus_bc
 - 7. calculus_sigma

*******************************************************************************/


/***** SECTION 1 ************************************
*       INCLUDE's
****************************************************/
#include "energy.h"


/***** SECTION 2 ************************************
 MAIN function: This function is the main function of the ENERGY task

- Read input parameters (call initModule)
- Open input FITS file (inFileA)
- Read and check input keywords (inFileA)
- If there are no pulses in input FITS file (inFileA) => The task finishes
- Get structure of input FITS file columns (inFileA)
- Allocate input GSL vectors (inFileA)
- Create structure to run Iteration: inDataIteratorA
- Read columns: Tstart, Pseudoenergy Quality and Grade (inFileA)
- Called iteration function: inDataIteratorA
- If calibration mode:
	- Open input FITS file (inFileB)
	- Read and check input keywords (inFileB)
	- If there are no pulses in input FITS file (inFileB) => The task finishes
	- Get structure of input FITS file columns (inFileB)
	- Allocate input GSL vectors (inFileB)
	- Create structure to run Iteration: inDataIteratorB
	- Read columns: Tstart, Pseudoenergy Quality and Grade (inFileB)
	- Called iteration function: inDataIteratorB
- Calculus of b_cF and c_cF
- Converting input pseudoenergy (arbitrary units) into output energy (eVolts) (inFileA)
- Call calculus_sigma (inFileA)
- Write ENERGY column in inFileA
- Write keywords in inFileA
- Close input FITS file (inFileA)
- If calibration mode:
	- Converting input pseudoenergy (arbitrary units) into output energy (eVolts) (inFileB)
	- Call calculus_sigma (inFileB)
	- Write ENERGY column in inFileB
	- Write keywords in inFileB
	- Close input FITS file (inFileB)
- Free allocate of memory
- Finalize the task
****************************************************/
int main(int argc, char **argv)
{
	char val[256];

	create = "energy v.12.0.0";			//Set "CREATOR" keyword of output FITS file
	time_t t_start = time(NULL);
	string message = "";
	int status = EPOK, extver=0;
	double keyval;

	sprintf(temporalFileName,"ENERGYauxfile");
	strcat(temporalFileName,".txt");
	temporalFile = fopen (temporalFileName,"w");
	if (temporalFile == NULL)
	{
	    message = "Cannot open auxiliary file ENERGYauxfile.txt";
	    writeLog(fileRef,"Error",verbosity,message);
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Read input parameters
	if (initModule(argc, argv))
	{
	    message = "Error in initModule";
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	writeLog(fileRef,"Log", verbosity,"Into Energy Module");
	  
	if ((optmode == 1) && (HorMorL != 0))
	{
		writeLog(fileRef,"Log",verbosity,"Production mode => Be careful! Resolution is going to be calculated");
	}
	else if ((optmode == 0) && (HorMorL == 0))
	{
		writeLog(fileRef,"Log",verbosity,"Calibration mode => Be careful! Resolution is NOT going to be calculated");
	}

	//	Open A input FITS file
	if (fits_open_file(&psgObjectA, psgNameA,READWRITE,&status))
	{
	    message = "Cannot open file " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(extname,"PSGRADE");
	if (fits_movnam_hdu(psgObjectA, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	message = "Open InFits1: " + string(psgNameA);
	writeLog(fileRef,"Log",verbosity,message);

	// Read	input keywords and check their values
	strcpy(keyname,"EVENTCNT");
	if (fits_read_key(psgObjectA,TLONG,keyname, &eventcntA,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	eventcntA_OK = eventcntA;
	if (eventcntA < 0)
	{
		message = "Legal values for EVENTCNT (PSGRADE inFileA) are non negative integer numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	if (eventcntA == 0)	// There are no pulses in the PSGRADE extension
	{
	    time_t t_end = time(NULL);
	    message = "There are no pulses in the input psg FITS file1 " + string(psgNameA);
	    writeLog(fileRef,"Warning", verbosity,message);
	    writeLog(fileRef,"OK", verbosity,"Energy Module OK");
			
	    sprintf(straux,"%f",(double) (t_end - t_start));
	    message = "Time: " + string(straux);
	    writeLog(fileRef,"Log", verbosity,message);
	    
	    if (fclose(fileRef))
	    {
	    	message = "Cannot close log file";
	    	EP_EXIT_ERROR(message,EPFAIL);
	    }
	}
	strcpy(keyname,"ENERGY");
	if (fits_read_key(psgObjectA,TDOUBLE,keyname, &energyA,comment,&status))
	{
	    message = "Cannot read key " + string(keyname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);	  
	}
	if (energyA <= 0) // After checking if EVENTCNT = 0 (no pulses) => ENERGY>0
	{
		message = "Legal values for ENERGY (PSGRADE inFileA) are real numbers greater than 0";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}

	// Get structure of input FITS file columns
	strcpy(straux,"TSTART");
	if (fits_get_colnum(psgObjectA,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"PSEUDOENERGY");
	if (fits_get_colnum(psgObjectA,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"QUALITY");
	if (fits_get_colnum(psgObjectA,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"GRADE");
	if (fits_get_colnum(psgObjectA,0,straux,&colnum,&status))
	{
	    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}

	// Allocate input GSL vectors
	tstartgslA = gsl_vector_alloc (eventcntA);			// GSL of input TSTART column
	energygslA = gsl_vector_alloc (eventcntA); 			// GSL of input ENERGY column
	pseudoenergygslA = gsl_vector_alloc (eventcntA);	// GSL of input PSEUDOENERGY column
	pseudoenergygslA_mod = gsl_vector_alloc (eventcntA);
	qualitygslA = gsl_vector_alloc (eventcntA); 		// GSL of input QUALITY column
	gradegslA = gsl_vector_alloc (eventcntA);			// GSL of input GRADE column

	extern int inDataIteratorA(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

	// Create structure to run Iteration
	iteratorCol cols [4];			// Structure of Iteration
	int n_cols = 4; 				// Number of columns:  Tstart + Pseudoenergy + Quality + Grade
	long rows_per_loop = 0; 	 	// 0: Use default: Optimum number of rows
	long offset = 0;				// 0: Process all the rows

	// Read columns: Tstart, Pseudoenergy Quality and Grade
	strcpy(straux,"Tstart");
	status = fits_iter_set_by_name(&cols[0], psgObjectA, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"Pseudoenergy");
	status = fits_iter_set_by_name(&cols[1], psgObjectA, straux, TDOUBLE, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"Quality");
	status = fits_iter_set_by_name(&cols[2], psgObjectA, straux, TSHORT, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	strcpy(straux,"Grade");
	status = fits_iter_set_by_name(&cols[3], psgObjectA, straux, TSHORT, InputCol);
	if (status)
	{
	    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}

	// Called iteration function: inDataIteratorA
	if (fits_iterate_data(n_cols, cols, offset, rows_per_loop, inDataIteratorA,0L,&status))
	{
	    message = "Cannot iterate data using InDataIteratorA";
	    EP_EXIT_ERROR(message,status);
	}

	if (optmode == 0)
	{
		if (fits_open_file(&psgObjectB, psgNameB,READWRITE,&status))
		{
		    message = "Cannot open file " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(extname,"PSGRADE");
		if (fits_movnam_hdu(psgObjectB, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to " + string(extname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		message = "Open InFits2: " + string(psgNameB);
		writeLog(fileRef,"Log",verbosity,message);

		strcpy(keyname,"EVENTCNT");
		if (fits_read_key(psgObjectB,TLONG,keyname, &eventcntB,comment,&status))
		{
		    message = "Cannot read key " + string(keyname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		eventcntB_OK = eventcntB;
		if (eventcntB < 0)
		{
			message = "Legal values for EVENTCNT (PSGRADE inFileB) are non negative integer numbers";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		if (eventcntB == 0)	// There are no pulses in the PSGRADE extension
		{
		    time_t t_end = time(NULL);
		    message = "There are no pulses in the input psg FITS file2 " + string(psgNameB);
		    writeLog(fileRef,"Warning", verbosity,message);
		    writeLog(fileRef,"OK", verbosity,"Energy Module OK");

		    sprintf(straux,"%f",(double) (t_end - t_start));
		    message = "Time: " + string(straux);
		    writeLog(fileRef,"Log", verbosity,message);

		    if (fclose(fileRef))
		    {
		    	message = "Cannot close log file";
		    	EP_EXIT_ERROR(message,EPFAIL);
		    }
		}

		strcpy(keyname,"ENERGY");
		if (fits_read_key(psgObjectB,TDOUBLE,keyname, &energyB,comment,&status))
		{
		    message = "Cannot read key " + string(keyname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);	  
		}
		if (energyB <= 0) // After checking if EVENTCNT = 0 (no pulses) => ENERGY>0
		{
			message = "Legal values for ENERGY (PSGRADE inFileB) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		
		drift = 0;

		strcpy(straux,"TSTART");
		if (fits_get_colnum(psgObjectB,0,straux,&colnum,&status))
		{
		    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"PSEUDOENERGY");
		if (fits_get_colnum(psgObjectB,0,straux,&colnum,&status))
		{
		    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"QUALITY");
		if (fits_get_colnum(psgObjectB,0,straux,&colnum,&status))
		{
		    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"GRADE");
		if (fits_get_colnum(psgObjectB,0,straux,&colnum,&status))
		{
		    message = "Cannot get column number for " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}

		tstartgslB = gsl_vector_alloc (eventcntB);
		energygslB = gsl_vector_alloc (eventcntB);
		pseudoenergygslB = gsl_vector_alloc (eventcntB);
		pseudoenergygslB_mod = gsl_vector_alloc (eventcntB);
		qualitygslB = gsl_vector_alloc (eventcntB);
		gradegslB = gsl_vector_alloc (eventcntB);

		extern int inDataIteratorB(long totalrows, long offset, long firstrow,long nrows, int ncols, iteratorCol *cols, void *user_strct);

		strcpy(straux,"Tstart");
		status = fits_iter_set_by_name(&cols[0], psgObjectB, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"Pseudoenergy");
		status = fits_iter_set_by_name(&cols[1], psgObjectB, straux, TDOUBLE, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"Quality");
		status = fits_iter_set_by_name(&cols[2], psgObjectB, straux, TSHORT, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		strcpy(straux,"Grade");
		status = fits_iter_set_by_name(&cols[3], psgObjectB, straux, TSHORT, InputCol);
		if (status)
		{
		    message = "Cannot iterate in column " + string(straux) +" in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}

		// Called iteration function: inDataIteratorB
		if (fits_iterate_data(n_cols,cols,offset,rows_per_loop,inDataIteratorB,0L,&status))
		{
		    message = "Cannot iterate data through inDataIteratorB ";
		    EP_EXIT_ERROR(message,status);	
		}

		// Calculus of b_cF and c_cF
		if (calculus_bc (eventcntA, eventcntA_OK, pseudoenergygslA_mod, energyA, eventcntB, eventcntB_OK, pseudoenergygslB_mod, energyB, &b_cF, &c_cF))
		{
		    message = "Cannot run routines calculus_bc in mode = 0";
		    EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	// Converting input pseudoenergy (arbitrary units) into output energy (eVolts) (inFileA)
	eventcntA_OK = eventcntA;
	for (int i=0;i<eventcntA;i++)
	{
		gsl_vector_set(energygslA,i,b_cF*gsl_vector_get(pseudoenergygslA,i)+c_cF*pow(gsl_vector_get(pseudoenergygslA,i),2));
		if (gsl_vector_get(qualitygslA,i) != 0)
		{
			gsl_vector_set(energygslA,i,0);
			eventcntA_OK = eventcntA_OK-1;
		}
	}

	if ((HorMorL == 1) || (HorMorL == 6))
	{
		energygslA1 = gsl_vector_alloc(eventcntA);
		gsl_vector_memcpy(energygslA1,energygslA);
	}
	if ((HorMorL == 2) || (HorMorL == 6))
	{
		energygslA2 = gsl_vector_alloc(eventcntA);
		gsl_vector_memcpy(energygslA2,energygslA);
	}
	if ((HorMorL == 3) || (HorMorL == 6))
	{
		energygslA3 = gsl_vector_alloc(eventcntA);
		gsl_vector_memcpy(energygslA3,energygslA);
	}

	eventcntA1_OK = eventcntA_OK;
	eventcntA2_OK = eventcntA_OK;
	eventcntA3_OK = eventcntA_OK;
	for (int i=0;i<eventcntA;i++)
	{
		if ((gsl_vector_get(qualitygslA,i) == 0) && ((gsl_vector_get(gradegslA,i) != 1)	&& ((HorMorL == 1) || (HorMorL == 6))))
		{
			gsl_vector_set(energygslA1,i,0);
			eventcntA1_OK = eventcntA1_OK-1;
		}
		if ((gsl_vector_get(qualitygslA,i) == 0) && (((gsl_vector_get(gradegslA,i) != 21) && (gsl_vector_get(gradegslA,i) != 22)) && ((HorMorL == 2) || (HorMorL == 6))))
		{
			gsl_vector_set(energygslA2,i,0);
			eventcntA2_OK = eventcntA2_OK-1;
		}
		if ((gsl_vector_get(qualitygslA,i) == 0) && (((gsl_vector_get(gradegslA,i) != 31) && (gsl_vector_get(gradegslA,i) != 32)) && ((HorMorL == 3) || (HorMorL == 6))))
		{
			gsl_vector_set(energygslA3,i,0);
			eventcntA3_OK = eventcntA3_OK-1;
		}
	}

	// If Production mode => Only with monochromatic pulses it makes sense calculating the energy resolution
	if ((HorMorL == 1) || (HorMorL == 6))
	{
		if (calculus_sigma(eventcntA, eventcntA1_OK, energygslA1, energyA, b_cF, c_cF, &sigmaA1))
		{
		    message = "Cannot run routine calculus_sigma for (HorMorL == 1) || (HorMorL == 6) in file A";
		    EP_EXIT_ERROR(message, EPFAIL);
		}
	}
	if ((HorMorL == 2) || (HorMorL == 6))
	{
		if (calculus_sigma(eventcntA, eventcntA2_OK, energygslA2, energyA, b_cF, c_cF, &sigmaA2))
		{
		    message = "Cannot run routine calculus_sigma for (HorMorL == 2) || (HorMorL == 6) in file A";
		    EP_EXIT_ERROR(message, EPFAIL);
		}
	}
	if ((HorMorL == 3) || (HorMorL == 6))
	{
		if (calculus_sigma(eventcntA, eventcntA3_OK, energygslA3, energyA, b_cF, c_cF, &sigmaA3))
		{
		    message = "Cannot run routine calculus_sigma for (HorMorL == 3) || (HorMorL == 6) in file A";
		    EP_EXIT_ERROR(message, EPFAIL);
		}
	}

	//  Creating Energy Column in the psg file
	obj.inObject = psgObjectA;
	obj.nameTable = new char [255];
	strcpy(obj.nameTable,"PSGRADE");
	obj.iniRow = 1;
	obj.endRow = eventcntA;
	obj.iniCol = 0;
	obj.nameCol = new char [255];
	strcpy(obj.nameCol,"Energy");
	obj.type = TDOUBLE;
	obj.unit = new char [255];
	strcpy(obj.unit,"eV");
	strcpy(extname,"PSGRADE");
	if (fits_movnam_hdu(psgObjectA, ANY_HDU,extname, extver, &status))
	{
	    message = "Cannot move to " + string(extname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	if (writeFitsSimple (obj,energygslA))
	{
	    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,EPFAIL);
	}

	// Set process keyword
	char str_verb[125];		sprintf(str_verb,"%d",verbosity);
	char str_b_cF[125];		sprintf(str_b_cF,"%f",b_cF);
	char str_c_cF[125];		sprintf(str_c_cF,"%f",c_cF);
	char str_optmode[125];		sprintf(str_optmode,"%d",optmode);

	string processoutenr (string("ENERGY") 	        	+ ' ' +
	string(str_optmode)		+ ' ' + string(psgNameA)	+ ' ' + string(psgNameB)	+ ' ' +
	string(str_b_cF)	    + ' ' + string(str_c_cF)	+ ' ' +
	string(nameLog)         + ' ' + string(str_verb)	+ ' ' + string(clobberStr) + ' ' +
	string("(")		        +      (string) create   	+   	string(")"));

	strcpy(keyname,"PROC2");
	strcpy(keyvalstr,processoutenr.c_str());
	if (fits_write_key_longwarn(psgObjectA,&status))
	{
	    message = "Cannot update keyword " + string(keyname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);
	}
	if (fits_update_key_longstr(psgObjectA,keyname,keyvalstr,comment,&status))
	{
		message = "Cannot write keyword " + string(keyname) + " in " + string(psgNameA);
		EP_PRINT_ERROR(message,status); return(EPFAIL);
	}

	// Energy Calibration in case of a linear relation between pseudoenergies and energies (E= b*e)
	strcpy(keyname,"EC_LIN");	// b
	if (fits_update_key(psgObjectA,TDOUBLE,keyname, &b_cF,comment,&status))
	{
	    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);	  
	}
	if (b_cF < 0)
	{
		message = "Legal values for EC_LIN (output PSGRADE) are non negative real numbers";
		writeLog(fileRef, "Error", verbosity, message);
		EP_EXIT_ERROR(message,EPFAIL);
	}
	
	// Energy Calibration in case of a quadratic relation between pseudoenergies and energies (E= b*e+c*e²)
	strcpy(keyname,"EC_NLIN");	// c
	if (fits_update_key(psgObjectA,TDOUBLE,keyname, &c_cF,comment,&status))
	{
	    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);	  
	}

	if ((HorMorL == 1) || (HorMorL == 6))
	{
		strcpy(keyname,"SIGMA_H");
		if (fits_update_key(psgObjectA,TDOUBLE,keyname, &sigmaA1,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA) + " for (HorMorL == 1) || (HorMorL == 6)";
		    EP_EXIT_ERROR(message,status);	  
		}
		if (sigmaA1 <= 0)
		{
			message = "Legal values for SIGMA_H (output PSGRADE inFileA) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		strcpy(keyname,"FWHM_H");
		keyval = 2.35*sigmaA1;
		if (fits_update_key(psgObjectA,TDOUBLE,keyname, &keyval,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA) + " for (HorMorL == 1) || (HorMorL == 6)";
		    EP_EXIT_ERROR(message,status);	  
		}
		if (keyval <= 0)
		{
			message = "Legal values for FWHM_H (output PSGRADE inFileA) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}
	if ((HorMorL == 2) || (HorMorL == 6))
	{
		strcpy(keyname,"SIGMA_M");
		if (fits_update_key(psgObjectA,TDOUBLE,keyname, &sigmaA2,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA) + " for (HorMorL == 2) || (HorMorL == 6)";
		    EP_EXIT_ERROR(message,status);	  
		}
		if (sigmaA2 <= 0)
		{
			message = "Legal values for SIGMA_M (output PSGRADE inFileA) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		strcpy(keyname,"FWHM_M");
		keyval = 2.35*sigmaA2;
		if (fits_update_key(psgObjectA,TDOUBLE,keyname, &keyval,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA) + " for (HorMorL == 2) || (HorMorL == 6)";
		    EP_EXIT_ERROR(message,status);	  
		}
		if (keyval <= 0)
		{
			message = "Legal values for FWHM_M (output PSGRADE inFileA) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}
	if ((HorMorL == 3) || (HorMorL == 6))
	{
		strcpy(keyname,"SIGMA_L");
		if (fits_update_key(psgObjectA,TDOUBLE,keyname, &sigmaA3,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA) + " for (HorMorL == 3) || (HorMorL == 6)";
		    EP_EXIT_ERROR(message,status);	  
		}
		if (sigmaA3 <= 0)
		{
			message = "Legal values for SIGMA_L (output PSGRADE inFileA) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		strcpy(keyname,"FWHM_L");
		keyval = 2.35*sigmaA3;
		if (fits_update_key(psgObjectA,TDOUBLE,keyname, &keyval,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameA) + " for (HorMorL == 3) || (HorMorL == 6)";
		    EP_EXIT_ERROR(message,status);	  
		}
		if (keyval <= 0)
		{
			message = "Legal values for FWHM_L (output PSGRADE inFileA) are real numbers greater than 0";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
	}

	string mod1 (string("File MODIFIED by") + ' ' +	(string) create);
	strcpy(keyname,"MOD1");
	strcpy(keyvalstr,mod1.c_str());
	if (fits_write_key(psgObjectA,TSTRING,keyname,keyvalstr,comment,&status))
	{
	    message = "Cannot write key " + string(keyname) + " in " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);	  
	}

	// Close input FITS files
	if (fits_close_file(psgObjectA,&status))
	{
	    message = "Cannot close file " + string(psgNameA);
	    EP_EXIT_ERROR(message,status);	  
	}

	if (optmode == 0)
	{
		// Converting input pseudoenergy (arbitrary units) into output energy (eVolts)
		eventcntB_OK = eventcntB;
		for (int i=0;i<eventcntB;i++)
		{
			gsl_vector_set(energygslB,i,b_cF*gsl_vector_get(pseudoenergygslB,i)+c_cF*pow(gsl_vector_get(pseudoenergygslB,i),2));
			if (gsl_vector_get(qualitygslB,i) != 0)
			{
				gsl_vector_set(energygslB,i,0);
				eventcntB_OK = eventcntB_OK-1;
			}
		}

		if ((HorMorL == 1) || (HorMorL == 6))
		{
			energygslB1 = gsl_vector_alloc(eventcntB);
			gsl_vector_memcpy(energygslB1,energygslB);
		}
		if ((HorMorL == 2) || (HorMorL == 6))
		{
			energygslB2 = gsl_vector_alloc(eventcntB);
			gsl_vector_memcpy(energygslB2,energygslB);
		}
		if ((HorMorL == 3) || (HorMorL == 6))
		{
			energygslB3 = gsl_vector_alloc(eventcntB);
			gsl_vector_memcpy(energygslB3,energygslB);
		}

		eventcntB1_OK = eventcntB_OK;
		eventcntB2_OK = eventcntB_OK;
		eventcntB3_OK = eventcntB_OK;
		for (int i=0;i<eventcntB;i++)
		{
			if ((gsl_vector_get(qualitygslB,i) == 0) && ((gsl_vector_get(gradegslB,i) != 1)	&& ((HorMorL == 1) || (HorMorL == 6))))
			{
				gsl_vector_set(energygslB1,i,0);
				eventcntB1_OK = eventcntB1_OK-1;
			}
			if ((gsl_vector_get(qualitygslB,i) == 0) && (((gsl_vector_get(gradegslB,i) != 21) && (gsl_vector_get(gradegslB,i) != 22)) && ((HorMorL == 2) || (HorMorL == 6))))
			{
				gsl_vector_set(energygslB2,i,0);
				eventcntB2_OK = eventcntB2_OK-1;
			}
			if ((gsl_vector_get(qualitygslB,i) == 0) && (((gsl_vector_get(gradegslB,i) != 31) && (gsl_vector_get(gradegslB,i) != 32)) && ((HorMorL == 3) || (HorMorL == 6))))
			{
				gsl_vector_set(energygslB3,i,0);
				eventcntB3_OK = eventcntB3_OK-1;
			}
		}

		// If Production mode => Only with monochromatic pulses it makes sense calculating the energy resolution
		if ((HorMorL == 1) || (HorMorL == 6))
		{
			if (calculus_sigma(eventcntB, eventcntB1_OK, energygslB1, energyB, b_cF, c_cF, &sigmaB1))
			{
			    message = "Cannot run routine calculus_sigma for (HorMorL == 1) || (HorMorL == 6) in file B";
			    EP_EXIT_ERROR(message, EPFAIL);
			}
		}
		if ((HorMorL == 2) || (HorMorL == 6))
		{
			if (calculus_sigma(eventcntB, eventcntB2_OK, energygslB2, energyB, b_cF, c_cF, &sigmaB2))
			{
			    message = "Cannot run routine calculus_sigma for (HorMorL == 2) || (HorMorL == 6) in file B";
			    EP_EXIT_ERROR(message, EPFAIL);
			}
		}
		if ((HorMorL == 3) || (HorMorL == 6))
		{
			if (calculus_sigma(eventcntB, eventcntB3_OK, energygslB3, energyB, b_cF, c_cF, &sigmaB3))
			{
			    message = "Cannot run routine calculus_sigma for (HorMorL == 3) || (HorMorL == 6) in file B";
			    EP_EXIT_ERROR(message, EPFAIL);
			}
		}

		//  Creating Energy Column in the psg file
		obj.inObject = psgObjectB;
		obj.nameTable = new char [255];
		strcpy(obj.nameTable,"PSGRADE");
		obj.iniRow = 1;
		obj.endRow = eventcntB;
		obj.iniCol = 0;
		obj.nameCol = new char [255];
		strcpy(obj.nameCol,"Energy");
		obj.type = TDOUBLE;
		obj.unit = new char [255];
		strcpy(obj.unit,"eV");
		strcpy(extname,"PSGRADE");
		if (fits_movnam_hdu(psgObjectB, ANY_HDU,extname, extver, &status))
		{
		    message = "Cannot move to " + string(extname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		if (writeFitsSimple (obj,energygslB))
		{
		    message = "Cannot run routine writeFitsSimple to write " + string(obj.nameCol);
		    EP_EXIT_ERROR(message,EPFAIL);
		}

		strcpy(keyname,"PROC2");
		strcpy(keyvalstr,processoutenr.c_str());
		if (fits_write_key_longwarn(psgObjectB,&status))
		{
		    message = "Cannot write keyword " + string(keyname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);
		}
		if (fits_update_key_longstr(psgObjectB,keyname,keyvalstr,comment,&status))
		{
			message = "Cannot write keyword " + string(keyname) + " in " + string(psgNameB);
			EP_PRINT_ERROR(message,status); return(EPFAIL);
		}

		strcpy(keyname,"EC_LIN");
		if (fits_update_key(psgObjectB,TDOUBLE,keyname, &b_cF,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);	  
		}
		if (b_cF < 0)
		{
			message = "Legal values for EC_LIN (output PSGRADE) are non negative real numbers";
			writeLog(fileRef, "Error", verbosity, message);
			EP_EXIT_ERROR(message,EPFAIL);
		}
		strcpy(keyname,"EC_NLIN");
		if (fits_update_key(psgObjectB,TDOUBLE,keyname, &c_cF,comment,&status))
		{
		    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);	  
		}
		if ((HorMorL == 1) || (HorMorL == 6))
		{
			strcpy(keyname,"SIGMA_H");
			if (fits_update_key(psgObjectB,TDOUBLE,keyname, &sigmaB1,comment,&status))
			{
			    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB) + " for (HorMorL == 1) || (HorMorL == 6)";
			    EP_EXIT_ERROR(message,status);	  
			}
			if (sigmaB1 <= 0)
			{
				message = "Legal values for SIGMA_H (output PSGRADE inFileB) are real numbers greater than 0";
				writeLog(fileRef, "Error", verbosity, message);
				EP_EXIT_ERROR(message,EPFAIL);
			}
			strcpy(keyname,"FWHM_H");
			keyval = 2.35*sigmaB1;
			if (fits_update_key(psgObjectB,TDOUBLE,keyname, &keyval,comment,&status))
			{
			    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB) + " for (HorMorL == 1) || (HorMorL == 6)";
			    EP_EXIT_ERROR(message,status);	  
			}
			if (keyval <= 0)
			{
				message = "Legal values for FWHM_H (output PSGRADE inFileB) are real numbers greater than 0";
				writeLog(fileRef, "Error", verbosity, message);
				EP_EXIT_ERROR(message,EPFAIL);
			}
		}
		if ((HorMorL == 2) || (HorMorL == 6))
		{
			strcpy(keyname,"SIGMA_M");
			if (fits_update_key(psgObjectB,TDOUBLE,keyname, &sigmaB2,comment,&status))
			{
			    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB) + " for (HorMorL == 1) || (HorMorL == 6)";
			    EP_EXIT_ERROR(message,status);	  
			}
			if (sigmaB2 <= 0)
			{
				message = "Legal values for SIGMA_M (output PSGRADE inFileB) are real numbers greater than 0";
				writeLog(fileRef, "Error", verbosity, message);
				EP_EXIT_ERROR(message,EPFAIL);
			}
			strcpy(keyname,"FWHM_M");
			keyval = 2.35*sigmaB2;
			if (fits_update_key(psgObjectB,TDOUBLE,keyname, &keyval,comment,&status))
			{
			    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB) + " for (HorMorL == 1) || (HorMorL == 6)";
			    EP_EXIT_ERROR(message,status);	  
			}
			if (keyval <= 0)
			{
				message = "Legal values for FWHM_M (output PSGRADE inFileB) are real numbers greater than 0";
				writeLog(fileRef, "Error", verbosity, message);
				EP_EXIT_ERROR(message,EPFAIL);
			}
		}
		if ((HorMorL == 3) || (HorMorL == 6))
		{
			strcpy(keyname,"SIGMA_L");
			if (fits_update_key(psgObjectB,TDOUBLE,keyname, &sigmaB3,comment,&status))
			{
			    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB) + " for (HorMorL == 1) || (HorMorL == 6)";
			    EP_EXIT_ERROR(message,status);	  
			}
			if (sigmaB3 <= 0)
			{
				message = "Legal values for SIGMA_L (output PSGRADE inFileB) are real numbers greater than 0";
				writeLog(fileRef, "Error", verbosity, message);
				EP_EXIT_ERROR(message,EPFAIL);
			}
			strcpy(keyname,"FWHM_L");
			keyval = 2.35*sigmaB3;
			if (fits_update_key(psgObjectB,TDOUBLE,keyname, &keyval,comment,&status))
			{
			    message = "Cannot update key " + string(keyname) + " in " + string(psgNameB) + " for (HorMorL == 1) || (HorMorL == 6)";
			    EP_EXIT_ERROR(message,status);	  
			}
			if (keyval <= 0)
			{
				message = "Legal values for FWHM_L (output PSGRADE inFileB) are real numbers greater than 0";
				writeLog(fileRef, "Error", verbosity, message);
				EP_EXIT_ERROR(message,EPFAIL);
			}
		}

		string mod1 (string("File MODIFIED by") + ' ' +	(string) create);
		strcpy(keyname,"MOD1");
		strcpy(keyvalstr,mod1.c_str());
		if (fits_write_key(psgObjectB,TSTRING,keyname,keyvalstr,comment,&status))
		{
		    message = "Cannot write key " + string(keyname) + " in " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);	  
		}

		if (fits_close_file(psgObjectB,&status))
		{
		    message = "Cannot close file " + string(psgNameB);
		    EP_EXIT_ERROR(message,status);	  
		}
	}

	// Free allocate of memory
	delete [] obj.nameTable;
	delete [] obj.nameCol;
	delete [] obj.unit;

	gsl_vector_free(tstartgslA);
	gsl_vector_free(energygslA);
	gsl_vector_free(pseudoenergygslA);
	gsl_vector_free(pseudoenergygslA_mod);
	gsl_vector_free(qualitygslA);
	gsl_vector_free(gradegslA);
	if (optmode == 0)
	{
		gsl_vector_free(tstartgslB);
		gsl_vector_free(energygslB);
		gsl_vector_free(pseudoenergygslB);
		gsl_vector_free(pseudoenergygslB_mod);
		gsl_vector_free(qualitygslB);
		gsl_vector_free(gradegslB);
	}

	// Finalize the task
	time_t t_end = time(NULL);
	sprintf(straux,"%f",(double) (t_end - t_start));
	message = "Time:" + string(straux);
	writeLog(fileRef,"Log", verbosity,message);

	writeLog(fileRef,"OK", verbosity,"Energy Module OK");
	
	if (fclose(temporalFile))
	{
	    message = "Cannot close temporalFile " + string(temporalFileName);
	    EP_EXIT_ERROR(message,EPFAIL);
	}

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
  
	// Define ENERGY input parameters and assign values to variables
	// Parameter definition and assignation of default values
	const int npars = 9, npars1 = 10;
	inparam energyPars[npars];
	int optidx =0, par=0, fst=0, ipar;
	string message, task="energy";

	energyPars[0].name = "optmode";
    energyPars[0].description = "Operation mode: Calibration(0) or Production(1)";
    energyPars[0].defValInt = 0;
    energyPars[0].type =  "int";
    energyPars[0].minValInt = 0;
    energyPars[0].maxValInt = 1;
	energyPars[0].ValStr = energyPars[0].defValInt;

	energyPars[1].name = "inFileA";
    energyPars[1].description = "Enter psg file name";
    energyPars[1].defValStr = "a_psg.fits";
    energyPars[1].type =  "char";
	energyPars[1].ValStr = energyPars[1].defValStr;

	energyPars[2].name = "inFileB";
    energyPars[2].description = "Enter psg file name";
    energyPars[2].defValStr = "b_psg.fits";
    energyPars[2].type =  "char";
	energyPars[2].ValStr = energyPars[2].defValStr;

	energyPars[3].name = "b_cF";
    energyPars[3].description = "Linear calibration factor";
    energyPars[3].defValReal = -1.;
    energyPars[3].type =  "double";
    energyPars[3].minValReal = -1;
    energyPars[3].maxValReal = 1.E50;
	energyPars[3].ValReal = energyPars[3].defValReal;

	energyPars[4].name = "c_cF";
    energyPars[4].description = "Quadratic calibration factor";
    energyPars[4].defValReal = -1.;
    energyPars[4].type = "double";
    energyPars[4].minValReal = -1;
    energyPars[4].maxValReal = 1.E50;
	energyPars[4].ValReal = energyPars[4].defValReal;

	energyPars[5].name = "gradesForER";
    energyPars[5].description = "Event grades to be used for Energy Resolution calculation: Hp (H), Mp+Ms (M), Lp+Ls (L), all (A) or none (N)";
    energyPars[5].defValStr = "H";
    energyPars[5].type = "char";
	energyPars[5].ValStr = energyPars[5].defValStr;

	energyPars[6].name = "nameLog";
    energyPars[6].description = "Output log file name";
    energyPars[6].defValStr = "enr_log.txt";
    energyPars[6].type = "char";
	energyPars[6].ValStr = energyPars[6].defValStr;

	energyPars[7].name = "verbosity";
    energyPars[7].description = "Verbosity level of the output log file (in [0,3])";
    energyPars[7].defValInt = 3;
    energyPars[7].type = "int";
    energyPars[7].minValInt = 0;
    energyPars[7].maxValInt = 3;
	energyPars[7].ValInt = energyPars[7].defValInt;
	
	energyPars[8].name = "clobber";
    energyPars[8].description = "Re-write output files if clobber=yes";
    energyPars[8].defValStr = "no";
    energyPars[8].type = "char";
	energyPars[8].ValStr = energyPars[8].defValStr;

	// Define structure for command line options
	static struct option long_options[npars1];
	for (optidx = 0; optidx<npars; optidx++)
	{
		long_options[optidx].name= energyPars[optidx].name.c_str();
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
			    	if(long_options[optidx].name == energyPars[i].name.c_str())
			    	{
			    		if(energyPars[i].type == "char") //save char value for par
			    		{
			    			energyPars[i].ValStr = optarg;
			    		}
			    		else // check if numeric value
			    		{
			    			if ( (!isdigit(optarg[0]) && (optarg[0] != '-')) ||
			    					(!isdigit(optarg[0]) && (optarg[0] == '-') && (!isdigit(optarg[1]))))
			    			{
			    				message = "Invalid value for input argument " + string(energyPars[i].name);
			    				EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			    			}
			    			if (energyPars[i].type == "int")
			    			{
			    				energyPars[i].ValInt = atoi(optarg);
			    			}
			    			else
			    			{
			    				energyPars[i].ValReal= atof(optarg);
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
		if (interactivePars(energyPars,npars,task))
		{
		    message = "Error reading parameters interactively";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL); 
		}
	}

	// Save parameter values into meaningful variables
	for(int i=0;i<npars; i++)
	{
		if (energyPars[i].name == "optmode")
		{
			optmode = energyPars[i].ValInt;
		}
		else if (energyPars[i].name == "inFileA")
		{
			strcpy(psgNameA, energyPars[i].ValStr.c_str());
		}
		else if (energyPars[i].name == "inFileB" &&  (optmode == 0))
		{
			strcpy(psgNameB, energyPars[i].ValStr.c_str());
		}
		else if (energyPars[i].name == "b_cF")
		{
			b_cF = energyPars[i].ValReal;
		}
		else if (energyPars[i].name == "c_cF")
		{
			c_cF = energyPars[i].ValReal;
		}
		else if (energyPars[i].name == "gradesForER")
		{
			if(strcmp(energyPars[i].ValStr.c_str(),"H")==0){
			    HorMorL = 1;
			}else if(strcmp(energyPars[i].ValStr.c_str(),"M")==0){
			    HorMorL = 2;
			}else if(strcmp(energyPars[i].ValStr.c_str(),"L")==0){
			    HorMorL = 3;
			}else if(strcmp(energyPars[i].ValStr.c_str(),"A")==0){
			    HorMorL = 6;
			}else if(strcmp(energyPars[i].ValStr.c_str(),"N")==0){
			    HorMorL = 0;
			}else{
			    message = "Parameter name " + string(energyPars[i].name) + " out of range: [H|M|L|A|N]";
			    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
			}
		}
		else if (energyPars[i].name == "nameLog")
		{
			strcpy(nameLog,energyPars[i].ValStr.c_str());
		}
		else if (energyPars[i].name == "verbosity")
		{
			verbosity = energyPars[i].ValInt;
		}
		else if (energyPars[i].name == "clobber")
		{
			strcpy(clobberStr, energyPars[i].ValStr.c_str());
			if(strcmp(clobberStr,"yes")==0){
			  clobber=1;
			}else{
			  clobber=0;
			}
		}
		
		// Check if parameter value is in allowed range
		if (energyPars[i].type == "int" &&
				(energyPars[i].ValInt < energyPars[i].minValInt ||
						energyPars[i].ValInt > energyPars[i].maxValInt))
		{
		    message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" + 
		    		boost::lexical_cast<std::string>(energyPars[i].minValInt) + "," + boost::lexical_cast<std::string>(energyPars[i].maxValInt) + "]";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);

		}
		else if (energyPars[i].type == "double" &&
				(energyPars[i].ValReal < energyPars[i].minValReal ||
						energyPars[i].ValReal > energyPars[i].maxValReal))
		{
		    message = "Parameter name " + string(long_options[optidx].name) + " out of range: [" + 
		    		boost::lexical_cast<std::string>(energyPars[i].minValReal) + "," + boost::lexical_cast<std::string>(energyPars[i].maxValReal) + "]";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		} 
		if ((HorMorL>3) && (HorMorL!=6))
		{
		    message = "Parameter name HorMorL out of range: [0,1,2,3,6]";
		    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if ((optmode == 1) && (b_cF < 0) && (b_cF != -1))
		{
		      message = "Parameter name b_cF cannot be a negative number /= -1";
		      EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
		}
		if ((optmode == 1) && (b_cF == -1.) && (c_cF == -1.))
		{
		    message = "Parameters b_cF & c_cF cannot be both negative number in production mode";
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
* inDataIteratorA function: This function takes the optimum number of rows to read the found pulses in inFileA
*                           and works iteratively
*
* - Declare variables
* - Read iterator
* - Processing each row of found pulses
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorA(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	int status = EPOK;
	string message = "";

	// Declare variables
	double *tstart, *tstartin;				// Vector of Tstart column
	double *pseudoenergy, *pseudoenergyin; 	// Vector of Pseudoenergy column
	short *quality, *qualityin; 			// Vector of Quality column
	short *grade, *gradein; 				// Vector of Grade column
	gsl_vector *tstart_block;
	gsl_vector *pseudoenergy_block;
	gsl_vector *quality_block;
	gsl_vector *grade_block;

	tstart_block = gsl_vector_alloc(nrows);
	pseudoenergy_block = gsl_vector_alloc(nrows);
	quality_block = gsl_vector_alloc(nrows);
	grade_block = gsl_vector_alloc(nrows);

	// Read iterator
	tstartin = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	tstart = &tstartin[1];
	if (toGslVector(((void **)&tstart), &tstart_block, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert Tstart column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	pseudoenergyin = (double *) fits_iter_get_array(&cols[1]);
	pseudoenergy = &pseudoenergyin[1];
	if (toGslVector(((void **)&pseudoenergy), &pseudoenergy_block, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert Pseudoenergy column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	qualityin = (short *) fits_iter_get_array(&cols[2]);
	quality = &qualityin[1];
	if (toGslVector(((void **)&quality), &quality_block, nrows, 0, TSHORT))
	{
	    message = "Cannot convert Quality column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gradein = (short *) fits_iter_get_array(&cols[3]);
	grade = &gradein[1];
	if (toGslVector(((void **)&grade), &grade_block, nrows, 0, TSHORT))
	{
	    message = "Cannot convert Grade column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}

	// Processing each row of found pulses
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set(tstartgslA,i+drift,gsl_vector_get(tstart_block,i));
		gsl_vector_set(pseudoenergygslA,i+drift,gsl_vector_get(pseudoenergy_block,i));
		gsl_vector_set(pseudoenergygslA_mod,i+drift,gsl_vector_get(pseudoenergy_block,i));
		gsl_vector_set(qualitygslA,i+drift,gsl_vector_get(quality_block,i));
		gsl_vector_set(gradegslA,i+drift,gsl_vector_get(grade_block,i));

		// Because the calculus of b and c only needs info of good pulses (quality=0 and grade=1)
		if ((gsl_vector_get(qualitygslA,i+drift) != 0) || (gsl_vector_get(gradegslA,i+drift) != 1))
		{
			gsl_vector_set(pseudoenergygslA_mod,i+drift,0);
			eventcntA_OK = eventcntA_OK-1;
		}
	}

	drift =drift + nrows;

	// Free allocate of GSL vectors
	gsl_vector_free(tstart_block);
	gsl_vector_free(pseudoenergy_block);
	gsl_vector_free(quality_block);
	gsl_vector_free(grade_block);

	return (EPOK);
}
/*xxxx end of SECTION 4 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 5 ************************************************************
* inDataIteratorB function: This function takes the optimum number of rows to read the found pulses in inFileB
*                           and works iteratively
*
* - Declare variables
* - Read iterator
* - Processing each row of found pulses
* - Free allocate of GSL vectors
****************************************************************************/
int inDataIteratorB(long totalrows, long offset, long firstrow, long nrows, int ncols, iteratorCol *cols, void *user_strct)
{
	int status = EPOK;
	string message = "";
	
	// Declare variables
	double *tstart, *tstartin;				// Vector of Tstart column
	double *pseudoenergy, *pseudoenergyin; 	// Vector of Pseudoenergy column
	short *quality, *qualityin; 			// Vector of Quality column
	short *grade, *gradein; 				// Vector of Grade column
	gsl_vector *tstart_block;
	gsl_vector *pseudoenergy_block;
	gsl_vector *quality_block;
	gsl_vector *grade_block;

	tstart_block = gsl_vector_alloc(nrows);
	pseudoenergy_block = gsl_vector_alloc(nrows);
	quality_block = gsl_vector_alloc(nrows);
	grade_block = gsl_vector_alloc(nrows);

	// Read iterator
	tstartin = (double *) fits_iter_get_array(&cols[0]);
	// NOTE: fits_iter_get_array because in this fits function the 1st element 
	//  of the output array is the null pixel value! 
	tstart = &tstartin[1];
	if (toGslVector(((void **)&tstart), &tstart_block, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert Tstart column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	pseudoenergyin = (double *) fits_iter_get_array(&cols[1]);
	pseudoenergy = &pseudoenergyin[1];
	if (toGslVector(((void **)&pseudoenergy), &pseudoenergy_block, nrows, 0, TDOUBLE))
	{
	    message = "Cannot convert Pseudoenergy column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	qualityin = (short *) fits_iter_get_array(&cols[2]);
	quality = &qualityin[1];
	if (toGslVector(((void **)&quality), &quality_block, nrows, 0, TSHORT))
	{
	    message = "Cannot convert Quality column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	gradein = (short *) fits_iter_get_array(&cols[3]);
	grade = &gradein[1];
	if (toGslVector(((void **)&grade), &grade_block, nrows, 0, TSHORT))
	{
	    message = "Cannot convert Grade column toGslVector";
	    EP_PRINT_ERROR(message,EPFAIL); return(EPFAIL);
	}
	
	// Processing each row of found pulses
	for (int i=0; i< nrows; i++)
	{
		gsl_vector_set(tstartgslB,i+drift,gsl_vector_get(tstart_block,i));
		gsl_vector_set(pseudoenergygslB,i+drift,gsl_vector_get(pseudoenergy_block,i));
		gsl_vector_set(pseudoenergygslB_mod,i+drift,gsl_vector_get(pseudoenergy_block,i));
		gsl_vector_set(qualitygslB,i+drift,gsl_vector_get(quality_block,i));
		gsl_vector_set(gradegslB,i+drift,gsl_vector_get(grade_block,i));

		if ((gsl_vector_get(qualitygslB,i+drift) != 0) || (gsl_vector_get(gradegslB,i+drift) != 1))
		{
			gsl_vector_set(pseudoenergygslB_mod,i+drift,0);
			eventcntB_OK = eventcntB_OK-1;
		}
	}

	drift =drift + nrows;

	// Free allocate of GSL vectors
	gsl_vector_free(tstart_block);
	gsl_vector_free(pseudoenergy_block);
	gsl_vector_free(quality_block);
	gsl_vector_free(grade_block);

	return (EPOK);
}
/*xxxx end of SECTION 5 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 6 ************************************************************
* calculus_bc function: This function calculates the linear and quadratic calibration factors in order to transfom from
*                       pseudoenergies (a.u.), e, into energies (eV), E.
*
* E = be +ce^2                                |
* In calibration mode => Monochromatic pulses | => SUM(i=1,nx)[bxi+cxi²]=nx·E0x
* 2 input FITS files:                         |    SUM(j=1,ny)[byj+cyj²]=ny·E0y
*   inFileA: nx values of pseudoenergies xi	  |
*            Energy = E0x                     |
*   inFileB: ny values of pseudoenergies yj   |
*            Energy = E0y                     |               |
*
* Solving by Cramer, b and c are obtained
****************************************************************************/
int calculus_bc (long nx, long nx_OK, gsl_vector *xi, double E0x, long ny, long ny_OK, gsl_vector *yj, double E0y, double *b_cF, double *c_cF)
{
	int status = EPOK;

	// Declare variables
	double x_= 0.0;
	double x2_= 0.0;
	double y_= 0.0;
	double y2_ = 0.0;

	for (int i=0; i < nx; i++)
	{
		if (gsl_vector_get(xi,i) != 0.0)
		{
			x_ = x_ + gsl_vector_get(xi,i);
			x2_ = x2_ + pow(gsl_vector_get(xi,i),2.0);
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
		}
	}
	y_ = y_/ny_OK;
	y2_ = y2_/ny_OK;

	*b_cF = (E0x*y2_ - E0y*x2_)/(x_*y2_ - y_*x2_);
	*c_cF = (E0y*x_ - E0x*y_)/(x_*y2_ - y_*x2_);

	return (EPOK);
}
/*xxxx end of SECTION 6 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/


/***** SECTION 7 ************************************************************
* calculus_sigma function: This function calculates the standard deviation of the calculated energies with respect to
*                          the known energy (E0, input keyword) (when b and c have been already calculated)
*
* E[(X-media)²]=sigma²                        |
* In calibration mode => Monochromatic pulses | => SUM(i=1,nx)[(bxi+cxi²-E0x)²]=(nx-1)·sigmax²
* 2 input FITS files:                         |    SUM(j=1,ny)[(byj+cyj²-E0y)²]=(ny-1)·sigmay²
*   inFileA: nx values of pseudoenergies xi	  |
*            Energy = E0x                     |
*   inFileB: ny values of pseudoenergies yj   |
*            Energy = E0y                     |
****************************************************************************/
int calculus_sigma (long n, long n_OK, gsl_vector *energygsl, double E0, double b, double c, double *sigma)
{
	int status = EPOK;

	*sigma = 0;

	for (int i=0; i < n; i++)
	{
		if (gsl_vector_get(energygsl,i) != 0)
		{
			*sigma = *sigma + pow(gsl_vector_get(energygsl,i)-E0,2.0);
		}
	}

	if (n_OK == 1) 	*sigma = *sigma/n_OK;
	else			*sigma = *sigma/(n_OK-1);

	*sigma = sqrt(*sigma);

	return (EPOK);
}
/*xxxx end of SECTION 7 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

