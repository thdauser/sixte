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
*  File:      trigger.h
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene Gonzélez Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
*  29/02/08    First version
*  12/03/08    Deleting variable "infile"
* 								Adding input parameter "deltaV"
*  01/04/08	   Adding function "truncate"
* 			   Adding function "findUpperK"
* 			   Adding function "findLowK"
*  			   Adding function "getQuality"
* 			   Adding function "setQuality"
* 			   Adding input parameter "numBitsQuality"
*  22/04/08    Moving functions "getQuality" and "set Quality" to miscellaneous module included in Utils package.
*  29/04/08	   Changing documentation of the module
*  26/06/08	   Adding input parameter "ql"
* 		       Adding output column "sigma"
*  10/09/08	   Adding input parametrer "nameLog"
* 			   Adding input parametrer "verbosity" and processing it using RIL libraries.
* 			   Adding new function "readInputKeywords"
* 			   Deleting in/output keyword "channelcount"
*  29/09/08	   Included "CREATE" and "PROCESS" keyword in the output FITS file
*			   Included new input parameter: writePulse
*			   Included new input parameter: selectRow
*  19/11/08    Included obtainTau function.
*  14/01/09    Included delDuplicated function.
* 			   Included oldTstart global variable.
*			   Included new input keywords
*  02/02/09    New function "Bins2Seconds".
* 			   Resolved bug of baseline Extension write.
*  20/02/09    Included, modificated and removed input keywords
*  24/02/09    Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
*  12/03/09    Changed documentation
*			   Included new columns in EUR-TRG extension: "Baseline and Sigma" of each pulses.
* 			   Used the keywords libraries.
*  14/04/09    The input parameter selectRow is not used
*              Added IVCAL
*  15/07/09	   Deleted non-used variables
*  05/08/10    Deleted functions "findMean2" & "findPulses2"
* 			   Deleted input parameters: "ColumnNameI", "W", "wb", "sizePulse", "dt_before", "dt_after", "k", "k2" & "n2"
*			   New input parameters added: "tauFALL" & "ntaus"
*			   New parameters in "writePulses" function
*			   New parameters in "procSegment" function and renamed ("procRow")
*  08/11/10	   A new parameter in "writePulses" function
*  			   A new parameter in "procRow" function
*  18/02/11	   I0_Filtered deleted in the output FITS file (not used in PULSESHAPE)(writePulses)
*  14/03/11    Deleted non used variables
*  22/03/11    Deleted some non necessary libraries
*  23/03/11	   Added new variable: "pi"
*  25/03/11    Deleted "ql" variable
*  29/03/11    Updated .h
*  04/04/11    "energy" added
*  17/05/11    seconds2Bins instead bins2Seconds
*  19/05/11    vectorFIL not included in procRow
*  25/05/11    New input parameters: samplesUp and nSgms
*              Deleted input parameter n (and n_b)
*              New parameter: safetyMarginTstart
*  09/06/11    "plspolar" added
*  23/06/11    New pulses models library input FITS file ("readLib" function and some parameters added)
*              sign not included in procRow
*  18/10/11    New output keyword "chngplrt"
*  27/10/11    New parameters to handle with code hardpoints:
*               stopCriteriaMKC
*               kappaMKC
*               limitCriteriaMAX
*               nsAftrtstart
*               levelPrvPulse
*               primaryThresholdCriteria
*  07/03/12    New input parameter "scaleFactor"
*  12/03/12    New input parameters "mode", "b_cF" and "c_cF"
*  26/02/12    Deleted variables related to Savitsky-Golay method (methodUsed, nL, nR)
*  14/12/12    Changes in 'procRow' and 'findSePrPulses' to iteratively look for pulses
*  22/01/13    find_model function renamed as 'find_modelOLD'
*              New 'find_model' function in order to interpolate between pulse shapes
*              New 'interpolate_model' function
*  30/01/13    getEnergy function renamed as 'getPulseHeight'
*  05/02/13    New 'Energy' column in the output FITS file where the estimated energy is written
*  14/02/13    Improved the modularity to find pulses and moved the functions used to find pulses to Utils
*              Deleted a no necessary library: <gsl/gsl_sort_vector.h>
*  26/04/13    tAftrtstart = 0
*  02/05/13    'Atstart' to calculate the tstart precisely (pulse broadening due to the filtering)
*  03/07/13	   tAftrtstart != 0
*  01/10/13	   Variables related to PRETRIGS have been deleted or modified because PRETRIGS info is
*              already not used in XRAYCHAIN
*  06/11/13    tAftrtstart has changed from 'constant' to 'input parameter'
*  18/12/13	   New input parameter 'getTaus'
*  24/03/14	   Added a new output extension EUR-TEST to IFCA tests
*  07/07/14    Adapted to PIL/RIL/Common removed dependencies & removed some includes (set in Utils)
*  03/10/14    New parameters for 'writePulses'
*              New functions 'createLibrary' and 'writeLibrary'
*  15/10/14    New column TAILBFR, just in case there is a tail before the first pulse of the row
*  Dec/14      Deleted some unnecessary input parameters
*  	           Deleted some unnecessary output keywords
*  	           Deleted some constants and variables not used
*  	           New functions 'calculateTemplate', 'createHisto', 'align', 'shiftm' and 'shift_m'
* ********************************************************************************************/

#ifndef TRIGGER_H_
#define TRIGGER_H_

// Utils module

	#include <inoututils.h>
	#include <unistd.h>
	#include <pulseprocess.h>

//MC FOR DAL -> CFITSIO migration

	int  colnum=0, felem=0;
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];
	int evtcnt=0, ttpls1=0, modeval=0,eventcntLib1=0, lib_id=0;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Provisional => To be deleted in future
	int indice = 0;		// Provisional to export the event data to a file!!!!!!!!!!!!!!!!!!!!!!!!!
                    	// Functions in pulseprocess (Utils) make use of it to write info to an auxiliary file
	FILE * temporalFile;
	char temporalFileName[255];
	FILE * temporalFile2;
	char temporalFileName2[255];
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Constants

	double safetyMarginTstart = 50e-6;
	double stopCriteriaMKC = 1.0;  			// Used in medianKappaClipping
                               	   	   		// Given in %
	double kappaMKC = 3.0;					// Used in medianKappaClipping
	double levelPrvPulse = 100.0;  		    // Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse

// INPUT FILES & OUTPUT FILES

	fitsfile *inObject=NULL; 		// Object which contains information of the input FITS file
	char inName[255];				// Name of the input FITS file

	fitsfile *inLibObject=NULL; 	// Object which contains information of the pulses templates library FITS file
	char inLibName[255];			// Name of the pulses templates library FITS file

	FILE *fileRef;					// Pointer for file which contains errors and warnings

	fitsfile *trgObject=NULL;	    // Object which contains information of the output FITS file
	char trgName[255];				// Name of the output FITS file

	char *unit=NULL, *comment=NULL;

// INPUT KEYWORDS

	long eventsz;		// TRIGGSZ
	double samprate;	// Related to DELTAT
	long eventcnt;		// Number of rows
	double energy;		// Related to MONOEN

	double ivcal;		// Just in case it would be necessary
	double asquid;		// Just in case it would be necessary
	double plspolar;	// Just in case it would be necessary

// INPUT PARAMETERS (and others closely related to them)
	
	int mode;				// Operation mode

	double LrsT;			// Running sum length (in the RS filter case): T -> Time => Seconds
	double LbT;				// Baseline averaging length (in the RS filter case): T -> Time => Seconds
	double Lrs;				// LrsT in samples
	double Lb;				// LbT in samples
	gsl_vector *Bgsl;
	gsl_vector *Lbgsl;

	double tauFALL;			// Fall time of the pulses
	double scaleFactor; 	// Scale factor to apply to the fall time of the pulses in order to calculate the LPF box-car length
	int samplesUp;			// Consecutive samples over the threshold to locate a pulse
	double nSgms;			// Number of Sigmas to establish the threshold

	int ntaus;				// Number of tauFALLs after starting the pulse to calculate its tend
	double sizePulse; 		// Size of pulse (seconds) (ntaus*tauFALL)
	int sizePulse_b;		// sizePulse in bins

	int writePulse = 1;		// Write pulses in the output FITS file. Default value = true
	int getTaus = 0;		// Calculate the approximate rise and fall times of the pulses. Default value = false
	char nameLog[255];		// Output log file name
	int verbosity;			// Verbosity level of the output log file

// OUTPUT VECTORS (in calibration mode to write the pulse templates library)

	gsl_vector *energygsl = gsl_vector_alloc(1);
	gsl_vector *estenergygsl = gsl_vector_alloc(1);
	gsl_matrix *pulsesgsl;

// AUXILIARY VARIABLES

	//To avoid the deprecate conversions
	char *straux = new char[255];

	bool append;				// Pulse templates library FITS file new (append=false) or not (append=true)
	long eventcntLib;			// Rows in the EUR-LIB extension of the pulse templates library output FITS file

	int totalpulses = 1; 		// Total number of founded pulses. Be careful because it is initialized to 1
	int ntotalrows = 1;			// Used in inDataIterator functions
	double initialtime = 0;

	int nPulsesRow = 0;			// Number of pulses of each record (record is equivalent to row)

	long nummodels;				// Number of the pulse templates (number of rows) included in the pulses templates library file
	gsl_vector *model;			// Pulse template which is going to be used as model
	                            // (selected row of the PULSE column of the EUR-LIB extension from the pulses templates library file)
								// It will be overwritten with the its first derivative
	gsl_vector *modelSGN;       // Sign of the first derivative of the pulse template
	gsl_matrix *library;		// Energy - EstEnergy
	gsl_matrix *models;			// Matrix where all the pulse templates of the library are going to be stored

	// Parameters used to write output FITS files
	IOData obj;

	gsl_vector *tstartout;		// Tstart column from the output trg FITS file
	gsl_vector *pulseheight;	// EstEnrgy column from the output trg FITS file
	gsl_vector *qualityout;		// Quality column from the output trg FITS file

// OUTPUT KEYWORDS

	const char *create;		// Name and version of the module: name-0.0.0
	long chngplrt = 0; 		// 1 => Polarity changed (pulses multiplied by -1)
	                        // 0 => Polarity not changed

//FUNCTIONS

	int initModule(int argc, char **argv);
	
	int readLib ();
	int createLibrary();

	int createTriggerFile();

	int procRow(gsl_vector *vectorNOTFIL, gsl_vector *vectorDER, int *npulses);
	int obtainTau (gsl_vector *invector, gsl_vector *tstartgsl, gsl_vector *tendgsl, int nPulses, gsl_vector **taurisegsl, gsl_vector **taufallgsl);
	int writePulses(gsl_vector *invectorNOTFIL, gsl_vector *invectorDER, gsl_vector *tstart, gsl_vector *tend, gsl_vector *quality, gsl_vector *taurise, gsl_vector *taufall, gsl_vector *energy, gsl_matrix **pulses);

	int calculateTemplate (int totalPulses, gsl_vector **pulseaverage, double *pulseaverageHeight);
	int createHisto (gsl_vector *invector, int nbins, gsl_vector **xhisto, gsl_vector **yhisto);
	int align(gsl_vector **vector1, gsl_vector ** vector2);
	int shiftm(gsl_vector *vectorin, gsl_vector *vectorout, int m);
	int shift_m(gsl_vector *vectorin, gsl_vector *vectorout, int m);

	int writeLibrary(double estenergy, gsl_vector **pulsetemplate);

	using namespace std;

#endif /*TRIGGER_H_*/
