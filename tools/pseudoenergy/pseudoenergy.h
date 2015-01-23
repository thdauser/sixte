/*******************************************************************************************
* This file is part of SIXTE/SIRENA software

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
*
*				  
*                      ENERGY RESOLUTION APLICATION			          					      
*												     					
*                                                                  						            
*  File:      pseudoenergy.h
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*                                                                   					
*  Revision History:                                              						                                                                  					   
*  20/09/06   	First version                        					
*  04/10/06   	Error Processing                   					     
*  10/10/06   	Energy resol only calculate pulse parameters. NOT calculate fit parameters 	
*  12/12/06 	Add to energy resolution module holzgauss module
* 			   	Delete shiftright and shiftleft functions			
*  24/04/07 	Include MAX_CHANNEL_COUNT :a value for the maximum number of data channels in an input file				
*  12/03/08		Adding the function "init module"
*				Adding two inputs parameters: "Central Energy" and "house keeping"
* 				Deleting the keyword "yscale"
*				Deleting the column "Baseline"
*  07/04/08		Adding the output keyword "ers_id"
* 				Adding the output column "quality"
*  26/06/08		Adding input/output keyword "ql"
* 				Adding function "fit_poly"
* 				Adding function "energyCalibration"
*  10/09/08		New documentation in the source code
* 				Deleting poly_fit function. It has been include in "miscellaneous" module.
* 				Adding input parametrer "nameLog"
* 				Adding input parametrer "verbosity" and processing it using RIL libraries.
* 				Deleting in/output keyword "channelcount"
*  29/09/08		Included "CREATE" and "PROCESS" keyword in the output FITS file
* 				Included new input keywords.
*  02/02/09		Removed energyResol and Energycalibration functions.
*  20/02/09		Included, modificated and removed input keywords
*  24/02/09		Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
*  12/03/09		Changed documentation
*  --/--/--		Used the keywords libraries
*  16/04/09     Deleted some non-used libraries
*  27/05/09     New library added (gsl_sf_gamma) to integrate the Norris function
*  02/06/09	    Pseudoenergy_INT column only added in _psh FITS file if a Norris fit has been made in PULSESHAPE
*  30/11/09 	Adding "processin" and "history" variables
*  30/09/10     Restructuring the columns of the input FITS files from TRIGGER & Restructuring the columns of the
*               input FITS files from PULSESHAPE: delete "t0gsl" and "endpulsegsl"
*               Energy is not calculated by integrating the Norris function because the Norris fitting
*               is not made in PULSESHAPE (delete inDataIteratorPsh)
*  08/11/10		Now "proccesin" variable has more memory space
*  				New input parameters "Tinitial", "Tfinal" & "inFile"
*				Added new functions: "writeEnergy" & "writeKeywords"
*  21/03/11     "columnName" input parameter not used
*  25/03/11     Deleted "ql" variable
*  29/03/11     Updated .h
*  05/04/11     New variable "noPulsesTrg"
*  18/10/11     New function "calculateConvolution"
*               New input keyword "chngplrt" from EUR_TRG
*  09/11/11     New output vector 'qualityoutgsl' to write a Quality column also in the EUR-ENR extension
*  26/01/12     Choose the domain to filter (time or frequency):
*                 New variable "TorF"
*                 New input parameter "domain" in "findEnergy" function
*                 New parameters in 'findEnergy' (new variable 'pseudoenergy') and new functionality
*  20/02/12    Running sum filter
*  17/04/14    Read the column OF_L from the _trg FITS file
*  14/05/14    Read the column NRMFCTR from the _trg FITS file
*  05/07/14    Deleted input parameter 'numshift' (or the whole convolution is calculated or the convolution is not used)
*  07/07/14    Adapted to removal of PIL/RIL/Common dependencies. Removed redundant includes
*  10/07/14    EUR-ENR no longer used => Pseudenergy column added to the EUR-PSH extension => Some variables deleted
*     12/14    DAL -> CFITSIO migration; Restructuring error functionality
*              Deleted some unnecessary functions
*              Deleted some unnecessary input parameters
*              Deleted some unnecessesary keywords
*  12/01/15    Deleted 'processin'
*  21/01/15    Variable renaming for parameters
*  22/01/15    Renaming of some input parameters:
* 				    TorF ->filterDomain ;  Hp_OForRS -> filterHp ; Mp_OForRS -> filterMp ; Ms_OForRS -> filterMs
* 				    Lp_OForRS -> filterLp ; Ls_OForRS -> filterLs
* ******************************************************************************************/

#ifndef PSEUDOENERGY_H_
#define PSEUDOENERGY_H_

// Utils module

	#include <inoututils.h>
	#include <pulseprocess.h>

// General

	char val[256];

// CFITSIO helpers

	int  colnum=0, felem=0, keyvalint=0;
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];

// INPUT FILES

	fitsfile *inObject=NULL; 	// Object which contains information of the input FITS file to the chain
	char inName[255];			// Name of the input FITS file

	fitsfile *trgObject = NULL;	// Object which contains information of the _trg input FITS file
	char trgName[255];			// Name of the _trg input FITS file
	
	fitsfile *fltObject = NULL;	// Object which contains information of the _flt input FITS file
	char fltName[255];			// Name of the _flt input FITS file
	
	FILE * fileRef;				// Pointer for file which contains errors and warnings

	FILE * temporalFile;
	char temporalFileName[255];

// INPUT KEYWORDS

	//EUR-TRG extension (from the _trg input FITS file)
	long chngplrt;
	long mode;
	long eventcnt;
	long eventsz;
	double samprate;
	
	//RAW_ADC_DATA extension (from the input FITS file to the chain)
	long eventsz_in;
	long eventcnt_in;

	//EUR-FLT extension (from the _flt input FITS file)
	long fltid;
	long eventcnt_flt;
	
//INPUT VECTORS
	
	gsl_matrix *pulses;				// Gsl vector that contains values of pulses column
	gsl_vector *tstartgsl;			// Gsl vector that contains values of tstart column
	gsl_vector *tendgsl;			// Gsl vector that contains values of tend column
	gsl_vector *qualitygsl;			// Gsl vector that contains values of quality column
	gsl_vector *difgsl;				// Gsl vector that contains values of difTstrt column
	gsl_vector *gradegsl;			// Gsl vector that contains values of grade column

	gsl_vector *filtergsl; 			// Gsl vector that contains values of filter column
	gsl_matrix *optfiltergslMATRIX;	// Gsl vector that contains values of optimalf column
	long ind = 0;
	gsl_vector *filter2use;
	gsl_vector *filter2use1;
	gsl_vector *nrmfctrgsl;			// Gsl vector that contains values of NRMFCTR column
	double normalizationFactor;

	gsl_vector *timegsl;			// Gsl vector that contains values of TIME column.
	gsl_matrix *eventsgsl;			// Gsl vector that contains values of I0 column.

//INPUT PARAMETERS

	char nameLog[255];		// Output log file name
	int verbosity;			// Verbosity level of the output log file
	char clobberStr[4];		// Clobber=yes then overwritte output files	
	int clobber=0;
	int TorF;			// Time domain(0) or frequency domain(1)
	char filterDomain[2];		// T -> time domain; F-> freq. domain	


	int Hp_OForRS;			// Optimal filter (1) or running sum filter (0) for High Primary pulses
	int Mp_OForRS;			// Optimal filter (1) or running sum filter (0) for Med Primary pulses
	int Ms_OForRS;			// Optimal filter (1) or running sum filter (0) for Med Secondary pulses
	int Lp_OForRS;			// Optimal filter (1) or running sum filter (0) for Low Primary pulses
	int Ls_OForRS;			// Optimal filter (1) or running sum filter (0) for Low Secondary pulses
	char filterHp[3];		// Optimal filter (OP) or running sum filter (RS) for High Primary pulses
	char filterMp[3];		// Optimal filter (OP) or running sum filter (RS) for Med Primary pulses
	char filterMs[3];		// Optimal filter (OP) or running sum filter (RS) for Med Secondary pulses
	char filterLp[3];		// Optimal filter (OP) or running sum filter (RS) for Low Primary pulses
	char filterLs[3];		// Optimal filter (OP) or running sum filter (RS) for Low Secondary pulses
	
	
	double LrsT;			// Running sum length (in the RS filter case): T -> Time => Seconds
	double LbT;				// Baseline averaging length (in the RS filter case): T -> Time => Seconds
	long Lrs;				// LrsT in samples
	long Lb;				// LbT in samples
	double B; 				// Summation of Lb samples into a pulse-free interval

//AUXILIARY VARIABLES
	
	//To avoid the deprecate conversions
	char *straux = new char[255];
		
	gsl_vector *pulse;
	gsl_vector *row_input;
	
	int totalsegments = 1; 		// Total number of segments processed
	int totalpulses = 0; 		// Total number of pulses processed
	int ntotalrows = 1;			// Total number of row processed
	int i_filter=0;				// i_filter processed
	int fltSz;					// filter size
	
	int iterator_trg = -1;		// iterator in inDataIterator_trg
	int iterator_in;			// iterator in inDataIterator
	double tend_p = -10;		// previous tend (i-1)
	double TIME_p;				// previous TIME (i-1)
	double time_input;			// TIME column (xray.fits)
	int size_trg;				// size of each pulse-free interval in inDataIterator_trg

	double pseudoenergy;

// IN/OUTPUT FILE
	
	fitsfile *psgObject = NULL;		// Object which contains information of _psg input/output FITS file
	char psgName[255];				// Name of the _psg input/output FITS file
	
	char *unit=NULL, *comment=NULL;

// OUTPUT KEYWORDS
	
	const char *create;					// Name and version of the module: name-0.0.0
	//char *processin = new char[1023];	// Shell command (including input parameters) used to create FITS file
	
// OUTPUT VECTORS
	
	gsl_vector *energygsl;				// Pseudoenergy of each pulse
	
// FUNCTIONS
	
	int initModule(int argc, char **argv);
	
	int findEnergy (gsl_vector *vector, gsl_vector *filter, int domain, double nrmfctr, double *pseudoenergy);
	int RS_filter (gsl_vector *vector, long Lrs, long Lb, double B, double *pseudoenergy);

	int writeEnergy();
	int writeKeywords();
	
	using namespace std;

#endif /*PSEUDOENERGY_H_*/
