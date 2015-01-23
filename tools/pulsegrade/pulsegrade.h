/***********************************************************************

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

/***********************************************************************
*                      	PULSEGRADE
*
*  File:      pulsegrade.h
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
* 02/10/06   	First version
* 24/04/07		Include MAX_CHANNEL_COUNT: a value for the maximum number of data channels in an input file
*
* 07/04/08		Pulse Shape reads from trg fits file.
*				Restructuration of .h file
*				Adding output keyword "psh_id"
*				Adding the output column "quality"
* 26/06/08		Adding new package "gsl_fit.h"
*				Adding new input keywork "NBQUAL"
*				Adding findMultipulse function
*				Adding transformDatas function
*				Adding fitNorris function
*				Adding input/output keyword "ql"
*				Adding output Keyword in PSH FITS "NP-FIT"
*				Adding baseline and fit_baseline columns in output file
*				Adding getPulseParams function
*          		Adding input parameter "wRise" and "wFall"
*				Adding new columns "tstart" and "tend" in the output FITS file.
*				Deleting column "chisq " from output FITS file
* 10/09/08		Changing documentacion of file
* 				Adding variable called "version". It contains the number version of this module.
* 				Adding input parametrer "nameLog"
* 				Adding input parametrer "verbosity" and processing it using RIL libraries.
* 				Deleting in/output keyword "channelcount"
* 29/09/08		Included "CREATE" and "PROCESS" keyword in the output FITS file.
* 14/01/09		Changed the unit of wRise and wFall input parameter from bins to time
* 				Included new input keywords
* 19/01/09		Changed fixed value of wFall to input parameter wFall.
* 02/02/09		Included new error process:
* 				-63804 PSH_NON_FOUND_PULSES: Non found pulses. Check the input paramenters of trigger task.
* 				Include non valid parameter
* 20/02/09 		Included, modificated and removed input keywords
* 24/02/09		Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
* 12/03/09		Included new function called findSaturatedPulses.
* --/--/--		Used the keywords libraries
* 15/04/09      Deleted some non-used libraries
* 30/11/09  	Adding "processin" and "history" variables
* 10/09/10		Saturated pulses are now found. New function added findSaturatedPulses
*				Deleted input parameter fitabstol. Added input parameter CurrentLevel
*				Not necessary to cut most part of the pulse (tstart-Tend): pulse_cut
*				Deleted extension BASELINE in input FITS file
*				Input column difTstrt added. Input columns TIME, EndPulse, Baseline, Sigma deleted
*				Deleted functions fit_norris and norris
* 28/09/10		"fit_linear2" function has been modified and moved to Utils Library with name "fit_linear"
* 08/11/10		Now "proccesin" variable has more memory space
* 16/12/10		Deleted input parameters "wRise" & "wfall"
* 				Deleted output columns: "Trise", "Nrise", "Tfall" & "Nfall"
* 				Not necessary to read "tauRise", "tauFall" and "difTstrt" columns from xray_t.fits
* 				Deleted functions: "getPulseParams" & "transformDatas"
* 18/02/11      findMultipulses not used (inDataIteratorTrg and pulseShape modified):
*     				Necessary variables updated
*                   Non-used variables deleted
* 22/03/11      Deleted some non necessary libraries
*               "columnName" input parameter not used
* 25/03/11      Deleted "ql" variable
* 29/03/11      Updated .h
* 04/04/11      New variable "noPulsesTrg"
* 22/12/11		New function "eventGrade"
*               New input parameters: "tauFALL", "inTaus" and "outTaus"
* 06/02/13      New output keyword "MODE"
*               Saturated pulses are no longer searched for (already done in TRIGGER):
*               	'CurrentLevel' variable deleted
*               	'findSaturatedPulses' function deleted
* 08/02/13      Deleted 'goodPulse' variable
* 15/05/13		New output keywords "B_CF" and "C_CF"
* 27/11/13      ANNALS keyword substitutes HISTORY 'pseudokeyword'
* 07/07/14      Adapted to PIL/RIL/Common removed dependencies & removed some redundant includes
* 24/09/14      inTaus and outTaus input parameters have been changed from 'int' to 'double'
* 15/10/14      TailBfr column is read from _trg and if there is a tail before the pulse, its grade is fixed at 32
*   Dec/14      DAL -> CFITSIO migration
*    			'pulsegrade' instead of 'pulseshape'
*    			Deleted some unnecessary output keywords
* 12/01/15      Deleted 'processin'
* 20/01/15      Replacing psh related variables by psg names
***********************************************************************/

#ifndef PULSEGRADE_H_
#define PULSEGRADE_H_

// Utils module

	#include <inoututils.h>

// CFITSIO helpers

	int  colnum=0, felem=0;
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];

// INPUT FILE

	fitsfile *trgObject = NULL;	// Object which contains information of the input FITS file
	char trgName[255];			// Name of the input FITS file

	FILE *fileRef;				// Pointer for file which contains errors and warnings

	FILE * temporalFile;
	char temporalFileName[255];

// INPUT KEYWORDS

	long eventsz;
	long eventcnt;
	double energy;

//INPUT VECTORS

	gsl_vector *pulseF;	// *F: To store a pulse of a nrows-inDataIteratorTrg-set
	double tstartF;		//     and use it in the next nrows-inDataIteratorTrg-set
	double tendF;
	int qualityF;
	int tailbeforeF;
	gsl_vector *pulseG;	// *G: To store a pulse of a nrows-inDataIteratorTrg-set
	double tstartG;		//     and use it in the next nrows-inDataIteratorTrg-set
	double tendG;
	int qualityG;
	int tailbeforeG;

//INPUT PARAMETERS

	double tauFALL;		// Fall time of the pulses
	double inTaus;		// Times the fall time of the pulses to establish the inner interval for event grading
	double outTaus;		// Times the fall time of the pulses to establish the outer interval for event grading

	char nameLog[255];	// Output log file name
	int verbosity;		// Verbosity level of the output log file
	char clobberStr[4];	// Clobber=yes then overwritte output files	
	int clobber=0;
	
//AUXILIARY VARIABLES

	//To avoid the deprecate conversions
	char *straux = new char[255];

	// Parameters to use in the iteration
	int iteration = 0;		// Number of iteration
	int totalpulses = 1;	// Total number of pulses processed
	int ntotalrows = 1;
	int indexToWrite = 1;	// Index to write in the output FITS file (it must be different from totalpulses
							// because different nrows-inDataIteratorTrg-sets have to be used)

// OUTPUT FILE

	fitsfile *psgObject = NULL;		// Object which contains information of the output FITS file
	char psgName[255];				// Name of the output FITS file

	char *unit=NULL, *comment=NULL;

// OUTPUT KEYWORDS

	const char *create;				   		// Name and version of the module: name-0.0.0

// FUNCTIONS

	int initModule(int argc, char **argv);

	int createPulseShapeFile();

	int inDataIteratorTrg (long totalrows, long offset, long firstrow, long nrows,int ncols, iteratorCol *cols, void *user_strct);

	int pulseShape (gsl_vector *pulse, double tstart, double tend, double qual, double tail, double tstartprev, double tendprev, double tstartnext);

	int eventGrade (double Tstart, double Tstartprev, double Tstartnext, gsl_vector **Gradinggsl);

	using namespace std;

#endif /*PULSEGRADE_H_*/
