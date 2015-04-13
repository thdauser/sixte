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

   Copyright 2014:  GENNOISESPEC has been developed by the INSTITUTO DE FISICA DE
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01 and
   ESP2013-48637-C2-1-P.

/************************************************************************************************
*                      GENNOISESPEC
*
*  File:      gennoisespec.h
*  Developer: Beatriz Cobo Martí­n
* 			  cobo@ifca.unican.es
*             IFCA
*
*  Revision History:
*
* 02/06/09    First version
* 02/17/09    Changing input (and output) keywords
* 03/02/09    Deleting AmplificationFactor
* 03/09/09    Used the keywords libraries
* 07/01/09    Changed functions names to avoid problems with new GSl functions
* 11/13/09	  Added findIntervalN function
* 20/10/10    include's updated
*             File FITS resulting from TRIGGER is not used (not read)
*             Input parameters updated
*             getKhaux function deleted
*             Parameters of the finInterval and findIntervalN functions updated
* 18/11/10    include's updated
* 10/01/11    Included "miscellaneous.h" => Deleted functions to handle vectors, complex vectors and matrix
*             Added input keyword ASQUID
* 25/03/11    Deleted "ql" variable
* 29/03/11    Updated .h
* 05/04/11    "energy" added
* 26/05/11    New input parameters: samplesUp and nSgms
*             Deleted input parameter n (and n_b)
*             New parameter: safetyMarginTstart
* 25/08/11    "plspolar" added
*             New pulses models library input FITS file ("readLib" function and some parameters added)
* 13/02/13    New input parameter scaleFactor
*             New parameters to handle with code hardpoints:
*               stopCriteriaMKC
*               kappaMKC
*               limitCriteriaMAX
*               nsAftrtstart
*               levelPrvPulse
*               primaryThresholdCriteria
* 14/02/13    Functions 'readLib' and 'inDataIteratorLib' updated=> New variables 'library' and 'models'
*             New input parameters: b_cF, c_cF, LrsT, LbT
* 26/11/13    tAftrtstart has changed from 'constant' to 'input parameter'
* 13/01/15    Migrated to CFITSIO (removal of ISDC DAL)
*             Run dependency on datatype files (xray, iv or tesnoise) deleted
* 20/01/15    Parameter renaming for task renaming
* 27/01/15    The auxiliary file GENNOISESPECauxfile.txt not created
********************************************************************************************/

#ifndef GENNOISE_H_
#define GENNOISE_H_

// Utils module

	#include "genutils.h"
	#include "pulseprocess.h"
	#include "inoututils.h"

// CFITSIO helpers

	int  colnum=0, felem=0;
	char extname[20];
	char keyname[10];
	char keyvalstr[1000];
	char *tform[1];
	char *ttype[1];
	char *tunit[1];
	int evtcnt=0, keyvalint=0;

// Constants

	//double safetyMarginTstart = 50e-6;
	double stopCriteriaMKC = 1.0;  			// Used in medianKappaClipping
											// Given in %
	double kappaMKC = 3.0;					// Used in medianKappaClipping
	double levelPrvPulse = 100.0;  		    // Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse

// INPUT FILES

	fitsfile *infileObject = NULL;		// Object which contains information of the input FITS file
	char infileName[255];				// Name of the input FITS file

	FILE *fileRef;						// Pointer for file which contains errors and warnings

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Provisional => To be deleted in future
	FILE * temporalFile;
	/*char temporalFileName[255];*/
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// INPUT KEYWORDS

	long eventcnt;		// Number of rows
	long eventsz;		// TRIGGSZ
	double samprate;	// Related to DELTAT
	double energy;		// Related to MONOEN

	double ivcal;		// Just in case it would be necessary
	double asquid;		// Just in case it would be necessary
	double plspolar;	// Just in case it would be necessary

// INPUT VECTORS

	gsl_vector *timegsl;	// GENNOISESPEC has to look at RECORDS
	gsl_vector *ioutgsl;

//INPUT PARAMETERS

	double intervalMin;		// Minimum length of a pulse-free interval (seconds)
	int intervalMinBins;	// Minimum length of a pulse-free interval (bins)
	int intervalMinSamples;	// Minimum length of a pulse-free interval (samples)
	int ntausPF;			// Number of tauFALLs after ending the pulse to start the pulse-free interval
	int nintervals; 		// Number of pulse-free intervals to use
	double tauFALL;			// Fall time of the pulses (seconds)
	int tauFALL_b;			// Fall time of the pulses (samples)
	double scaleFactor; 	// Scale factor to apply to the fall time of the pulses in order to calculate the box-car length
	int samplesUp;			// Consecutive samples over the threshold to locate a pulse
	double nSgms;			// Number of Sigmas to establish the threshold
	int pulse_length;		// Pulse length (samples)
	double LrsT;			// Running sum length (in the RS filter case): T -> Time => Seconds
	double LbT;				// Baseline averaging length (in the RS filter case): T -> Time => Seconds
	double Lrs;				// LrsT in samples
	double Lb;				// LbT in samples
	double baseline;

	char nameLog[255];		// Output log file name
	int verbosity;			// Verbosity level of the output log file
	char clobberStr[4];		// Clobber=yes then overwritte output files	
	int clobber=0;
	
//AUXILIARY VARIABLES

	//To avoid the deprecate conversions
	char *straux = new char[255];

	// Parameters used to inDataIterator
	int ntotalrows = 1;				// Total number of rows processed

	//Used in relation with findInterval
	int nIntervals;					// Number of free-pulse intervals in an event
	gsl_vector *startIntervalgsl;	// In samples

	// Parameters used to write output FITS files
	IOData obj;

	int NumMeanSamples = 0;
	gsl_vector *EventSamplesFFTMean;
	gsl_matrix *EventSamplesFFT;

	gsl_matrix *library;		// Not used. Necessary only in order to be used as input parameter in findPulses
	gsl_matrix *models;			// Not used. Necessary only in order to be used as input parameter in findPulses

// OUTPUT FILE

	fitsfile *gnoiseObject = NULL;		// Object which contains information of output FITS file
	char gnoiseName[255];

	char *unit=NULL, *comment=NULL;

// OUTPUT KEYWORDS

	const char *creator;			// Name and version of the module: name v.0.0.0

// OUTPUT VECTORS

	gsl_vector *freqgsl;
	gsl_vector *csdgsl;			//Amount of current per unit (density) of frequency (spectral), as a function of the frequency
	gsl_vector *sigmacsdgsl;

// FUNCTIONS

	int initModule (int argc, char **argv);

	int findInterval(gsl_vector *invector, gsl_vector *startpulse, int npin, int pulse_length, int nPF, double tau, int interval, int *ni, gsl_vector **startinterval);
	int findIntervalN(gsl_vector *invector, int interval, int *ni, gsl_vector **startinterval);

	int createTPSreprFile ();
	int writeTPSreprExten ();

	int find_baseline(gsl_vector *invector, double kappa, double stopCriteria, int boxLPF, double *baseline, FILE *temporalFile);

	using namespace std;

#endif /*GENNOISESPEC_H_*/

