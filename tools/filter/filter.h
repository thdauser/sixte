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

/************************************************************************
*                 				FILTER
*
*  File:      filter.h
*  Developer: Beatriz Cobo Martín
* 			  cobo@ifca.unican.es
*             IFCA
*             Irene González Pérez
*             José Ramón Rodón Ortiz
*
*  Revision History:
*
* 12/03/08	First version
* 07/04/08	Adding the output keyword "flt_id"
* 			Adding the output column "quality"
* 26/06/08	Adding input/output keyword "ql"
* 			Adding input keyword "psh_id"
* 			Adding baseline and normalize function
* 			Adding input parameters "dt_before", "cut_filter" and "dt_after"
* 			New input FITS file "Pulse Shape"
* 10/09/08	Ad//+ding new function "readInputKeywords"
* 			Deleting in/output keyword "channelcount"
* 29/09/08	Included "CREATE" and "PROCESS" keyword in the output FITS file
* 			Include input parameters: "nameLog" and "verbosity"
* 14/01/09	Included new input keywords
* 20/02/09	Included, modificated and removed input keywords
* 24/02/09	Xray chain cans run input FITS files of type (XRAY, TESNOISE and IV)
* 12/03/09	Included new column called "filterError" in the output FITS file.
* 			Created obtainFilterError function.
* --/--/--  Used the keywords libraries
* 15/04/09  Deleted some non-used libraries
* 24/11/09  endPulse is not a hardpoint but a keyword read from _trg input FITS file (EVENTSZ)
* 25/11/09  EVENTSZ from _trg not included in the khaux
* 30/11/09 	Adding "processin" and "history" variables
* 24/09/10	Added input parameter "tmin"
*			Deleted input parameters "dt_after" & "dt_before"
*			Normalize the filter template ("normalize" function added)
*			Deleted "cut_filter" function
*			Deleted "writeFitsSimple2" function
* 30/09/10	khauxNumber=3 instead 2. Included keyword "EVENTCNT" in getKhaux
* 08/11/10	Now "proccesin" variable has more memory space
* 01/03/11  Deleted "processin" and "history"
* 17/03/11  _psh input FITS file is no longer used: QUALITY will be read from _trg
*              	inDataIteratorPsh deleted and inDataIteratorTrg modified
*               getKh deleted => getKhaux renamed as getKh
*               Some variable deleted: pshObject, pshExten...
* 25/03/11  Deleted some non necessary libraries
*           Deleted "ql" variable
* 29/03/11  Updated .h
* 05/04/11  New variable "noPulsesTrg"
* 31/01/12  New _noisespec input FITS file (to calculate the optimal filter):
*              "noisespecObject"
*              "noisespecExten"
*              "noisespecName"
*              "eventcnt_noisespec"
*              "freqgsl" and "csdgsl"
*           New libraries ("vector" and "complex")
* 05/03/12  difTstrt is not going to be used => Delete (difgsl)
*           tmin is not going to be used
*           Grade column of _trg is read (Hp, Mp, Ms, Lp and Ls new input parameters)
* 12/02/13  New output keyword "MODE"
* 25/03/13  New input parameters to take into account an energy window: Emin, Emax
*           Necessary changes to divide the thread between CALIBRATION and NORMAL mode
*           NORMAL mode: Read EstEnrgy column from _trg input FITS file
*           Added as input file the pulses library
*           Added pulseProcess library
* xx/03/13  Cut the ending part of the sum of all the pulses to avoid an undesirable minimum: 'filterlength'
* 10/05/13  New function "find_matchedfilter"
* 14/05/13  Delete some functions definitions in order to 'clean' the file
* 15/05/13  New output keywords "B_CF" and "C_CF"
* 07/01/14  New output column MF_L (matched filter length) to know where the artifact at the end of the
*           matched filter starts (mode = 0)
* 20/02/14  New output column OF_L (optimal filter length) to know where the artifact at the end of the
*           optimal filter starts
* 14/05/14  New output column NRMFCTR (normalization factor)
* 22/05/14  NOT USED -> create_filterNEW1, cretae_filterNEW2: Instead of aligning the pulses to be
*           averaged (in order to build the matched filter) by using the maximum of the pulse
*           (no filtered), it is going to low-pass filter the pulses, locate the maximum, align and
*           average the aligned NOT filtered pulses.
* 01/07/14  Deleted the 'normalize' function
*           New 'area0' function
* 07/07/14  Adapted for PIL/RIL/Common removed dependencies and command line options reading
*   Dec/14  DAL->CFITSIO migration
*   		Deleted inDataIteratorTrg0
*   		Deleted some input/output keywords
*   		Deleted some unnecessary input parameters (Emin, Emax)
* 12/01/15  Deleted 'annalsin'
********************************************************************************************/

#ifndef FILTER_H_
#define FILTER_H_

// Utils module

	#include <inoututils.h>
	#include <pulseprocess.h>

// CFITSIO helpers

	int  colnum=0, felem=0, keyvalint=0;
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];

	int status = EPOK;

// INPUT FILE

	fitsfile *trgObject	  = NULL; 		// Object which contains information of the _trg input FITS file
	char trgName[255];					// Name of the _trg FITS file

	fitsfile *noisespecObject = NULL; 	// Object which contains information of _noisespec input FITS file
	char noisespecName[255];			// Name of the _noisespec FITS file

	fitsfile *inLibObject     = NULL; 	// Object that contains information of pulses templates library input FITS file
	char inLibName[255];				// Name of the pulses templates library FITS file

	FILE *fileRef;						// Pointer for file which contains errors and warnings

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Provisional => To be deleted in future
	char val[256];
	char val_aux[256];
	FILE * temporalFile;
	char temporalFileName[255];
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// INPUT KEYWORDS

	//EUR-TRG extension (from the _trg input FITS file)
	long eventcnt;
	long eventsz;
	double samprate;
	long mode;

	//TPSrepr extension (from the _noisespec input FITS file)
	long eventcnt_noisespec;

//INPUT VECTORS

	//From the _trg input FITS file
	gsl_matrix *pulses;				// Input Pulses
	gsl_vector *qualitygsl;			// GSL vector that contains values of QUALITY column
	gsl_vector *gradegsl;			// GSL vector that contains values of GRADE column
	gsl_vector *energygsl;			// GSL vector that contains values of ENERGY column

	//From the _noisespec input FITS file
	gsl_vector *freqgsl;			// Frequency
	gsl_vector *csdgsl;				// Current noise spectral density

//INPUT PARAMETERS

	char nameLog[255];				// Output log file name
	int verbosity;					// Verbosity level of the output log file

	int Hp, Mp, Ms, Lp, Ls;			// To consider a pulse as a valid one to calculate the optimal filter

//AUXILIARY VARIABLES

	//To avoid the deprecate conversions
	char *straux = new char[255];

	gsl_vector *pulse;				// Input pulses

	int ntotalrows = 0;

	// Parameters used to read the library FITS file
	long nummodels;					// Number of pulse modes (rows number) included in the pulses templates library file
	gsl_vector *matchedfilter;		// Matched filter (selected row of the MATCHEDF column of the EUR-LIB extension from the library file)
									// It will be overwritten with the its first derivative
	gsl_vector *energylibrary;		// ENERGY column of the EUR-LIB extension from the library file
	gsl_vector *estenergylibrary;	// ESTENERGY column of the EUR-LIB extension from the library file
	gsl_matrix *templateslibrary;	// PULSE column of the EUR-LIB extension from the library file
	gsl_matrix *matchedfilters;		// Matrix where all the matched filters of the library are going to be stored
	int matchedf_exist = 1;			// 1: MATCHEDF column in library FITS file already exists
	                            	// 0: MATCHEDF column in library FITS file does not exist yet

	int iter;								// Number of pulses whose QUALITY is different from 0
	int iter1, iter2, iter3, iter4, iter5;	// Number of pulses whose GRADE is 1 if Hp=1 (etc...)

	// Parameters used to call the function 'find_energy'
	double energy;
	long energyInLibrary_row;

	// To express the matched filter (previous to the optimal filter)
	gsl_vector *filtergsl;			// Matched filter values (time domain)
	gsl_vector *freqfiltergsl;		// Matched filter frequency (frequency domain)

	// Parameters used to write data (all set to 0) from non valid pulses in EUR-FLT in _trg FITS file
	long nValidPulses = 0;			// Valid pulses number
	long sizeoptimalfiltergsl;		// Size of the optimal filter
	long sizeoptimalfilterFFTgsl;	// Size of the optimal filter expressed in frequency domain
	int alreadyValidPulse = 0;		// If 1 means than there has already been some valid pulse
	long nonValidPulsesBeginning = 0;	// If !=0 means there are some non valid pulses at the beginning

	// Parameters used to write library FITS file
	IOData obj;

// OUTPUT FILE

	fitsfile *fltObject 	 = NULL;    // Object which contains information of the _flt output FITS file
	char fltName[255];					// Name of the _flt FITS file

	char *unit=NULL, *comment= NULL;

// OUTPUT KEYWORDS

	const char *create;					// Name and version of the module: name-0.0.0

// OUTPUT VECTORS

	gsl_vector *optimalfiltergsl;		// Optimal filter expressed in the time domain (optimalfilter(t))
	gsl_vector *optimal_filterFFTgsl;	// Disordered optimal filter when f's are according to [0,...fmax,-fmax,...]
	                                    // (frequency domain)
	gsl_vector *f_optimal_filtergsl;	// Disordered optimal filter when magnitudes are according to [0,...fmax,-fmax,...]
	                                    // (frequency domain)
	gsl_vector *foptimalfiltergslaux;	// Auxiliary disordered optimal filter when f's are according to [0,...fmax,-fmax,...]
	                                    // (frequency domain)
	gsl_vector *optimalfilterFFTgslaux;	// Auxiliary disordered optimal filter when magnitudes are according [0,...fmax,-fmax,...]
	                                    // (frequency domain)

	double normalizationFactor;
	gsl_vector *nrmfctrgsl = gsl_vector_alloc(1);	// Normalization factor

// FUNCIONS

	int initModule(int argc, char **argv);

	int readLib ();

	int createFilterFile();

	int find_energy(double energyKeyword, gsl_vector *energygsl, long *rowFound);

	int calculus_matchedFilter(gsl_vector **matchedfiltergsl, long mf_size);

	int writeLib (double energyK, gsl_matrix *matchedfiltermatrix);

	int calculus_optimalFilter(gsl_vector *matchedfiltergsl, long mf_size, gsl_vector **optimal_filtergsl, gsl_vector **f_optimal_filtergsl_aux, gsl_vector **optimal_filterFFTgsl_aux);
	int reorderFFT(gsl_vector *invector, gsl_vector **outvector);
	int interpolate (gsl_vector *x_in, gsl_vector *y_in, long size, double step, gsl_vector **x_out, gsl_vector **y_out, long *numzerosstart, long *numzerosend);
	int disorderFFT(gsl_vector *invector, gsl_vector **outvector);

	int find_matchedfilter(double ph, gsl_vector *energiesvalues, gsl_matrix *matchedfilters, gsl_vector **matchedfilterFound, FILE * temporalFile);

	using namespace std;

#endif /*FILTER_H_*/
