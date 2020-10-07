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

   Copyright 2014:  GENNOISESPEC has been developed by the INSTITUTO DE FISICA DE 
   CANTABRIA (CSIC-UC) with funding from the Spanish Ministry of Science and 
   Innovation (MICINN) under project  ESP2006-13608-C02-01, and Spanish 
   Ministry of Economy (MINECO) under projects AYA2012-39767-C02-01, 
   ESP2013-48637-C2-1-P and ESP2014-53672-C3-1-P.

/***********************************************************************
*                      GENNOISESPEC
*
*  File:       gennoisespec.h
*  Developers: Beatriz Cobo
* 	       cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef GENNOISE_H
#define GENNOISE_H 1

// Utils module

	#include "genutils.h"
	#include "pulseprocess.h"
	#include "inoututils.h"

	#define TOOLSUB gennoisespec_main
        #include "headas_main.c"
	
	#include <time.h>
        
struct Parameters {
    
	//Input FITS file
	char inFile[MAXFILENAME];
        
    //Output FITS file
	char outFile[MAXFILENAME];

	//Base-2 minimum length of a pulse-free interval (samples)
	int intervalMinSamples;
        
    //Number of pulse lengths after the end of the pulse to start the pulse-free interval searching
	int nplPF;
        
        //Number of pulse-free intervals to use for the noise average
	int nintervals;
        
    //Scale Factor
	double scaleFactor;
	
	//Number of samples for threshold trespassing
	int samplesUp;
        	
	//Number of standard deviations in the kappa-clipping process for threshold estimation
	double nSgms;
        
    //Pulse length in samples
	int pulse_length;
        
    //Running sum (RS) length for the RS-filtering for raw energy estimation, in seconds
	double LrsT;
	
	//Baseline averaging length for the RS-filtering for raw energy estimation, in seconds
	double LbT;
        
    //Calculate and write the weight matrixes if weightMS=yes
	char weightMS;
        
    //Transform to resistance space (I2R, I2RALL, I2RNOL, I2RFITTED) or not (I)
	char I2R[10];
        
    //Boolean to choose whether to erase an already existing event list
	char clobber;
        
    //Size of noise matrix if only one to be created
	int matrixSize;
        
    //Remove some noise intervals before calculating the noise spectrum if rmNoiseIntervals=yes
	char rmNoiseIntervals;
       
	// END GENNOISESPEC input parameters
};
        	
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

	double stopCriteriaMKC = 1.0;  		// Used in medianKappaClipping
						// Given in %
	double kappaMKC = 3.0;			// Used in medianKappaClipping
	double levelPrvPulse = 100.0;  		// Secondary pulses must be 1/levelPrvPulse times larger than the preceding pulse

// INPUT FILES

	fitsfile *infileObject = NULL;		// Object which contains information of the input FITS file
	//char infileName[255];			// Name of the input FITS file

	FILE *fileRef;				// Pointer for file which contains errors and warnings

// INPUT KEYWORDS	

	long eventcnt;		        // Number of rows
	long eventsz;	         	// TRIGGSZ
	double samprate = -999.0;       // Related to DELTAT or TCLOCK+DEC_FAC or NROW+P_ROW
	double energy;		        // Related to MONOEN

	double aducnv;
	double ivcal;		// Just in case it would be necessary
	double asquid;		// Just in case it would be necessary
	double plspolar;	// Just in case it would be necessary
	
	double adu_cnv = -999.0;
    double i_bias = -999.0;
    double adu_bias = -999.0;
    int adu_cnv_exists = 0;
	
	double Imin = -999.0;
    double Imax = -999.0;
    double R0 = -999.0;
    double V0 = -999.0;
    double Ibias = -999.0;
    double RPARA = -999.0;
    double TTR = -999.0;
    double LFILTER = -999.0;
    int R0_error = 0;
    int V0_error = 0;
    int Ibias_error = 0;
    int RPARA_error = 0;
    int TTR_error = 0;
    int LFILTER_error = 0;

// INPUT VECTORS

	gsl_vector *timegsl;	// GENNOISESPEC has to look at RECORDS
	gsl_vector *ioutgsl;
    
//INPUT PARAMETERS
    
    struct Parameters par;
    
	int intervalMinBins;		// Minimum length of a pulse-free interval (bins)
	double Lrs;			// LrsT in samples
	double Lb;			// LbT in samples
	
	int weightMS=0;
        
        //char I2R[10];		        // 

	//char clobberStr[4];		// Clobber=yes then overwritte output files	
	//int clobber=0;
	
//AUXILIARY VARIABLES

    string message = " ";
    
	//To avoid the deprecate conversions
	char *straux = new char[255];

	// Parameters used to inDataIterator
	int ntotalrows = 1;		// Total number of rows processed

	//Used in relation with findInterval
	int nIntervals;			// Number of free-pulse intervals in an event
	gsl_vector *startIntervalgsl;	// In samples

	// Parameters used to write output FITS files
	IOData obj;

	int NumMeanSamples = 0;
	gsl_vector *EventSamplesFFTMean;
	//gsl_matrix *EventSamplesFFT;

	gsl_matrix *library;		// Not used. Necessary only in order to be used as input parameter in findPulses
	gsl_matrix *models;		// Not used. Necessary only in order to be used as input parameter in findPulses
	
	int indexBaseline = 0;
	gsl_vector *baseline;
	gsl_vector *sigma;
	
	gsl_matrix *noiseIntervals;
    gsl_vector *weightpoints;
    gsl_matrix *weightMatrixes;
    
    int tessimOrxifusim = -999;     // 0: tessim, 1: xifusim
    
    double deltat;
    double tclock;
    double dec_fac;
    
// OUTPUT FILE

	fitsfile *gnoiseObject = NULL;	// Object which contains information of output FITS file
	char gnoiseName[MAXFILENAME];

	char *unit=NULL, *comment=NULL;

// OUTPUT KEYWORDS

	const char *creator;		// Name and version of the module: name v.0.0.0

// OUTPUT VECTORS

	gsl_vector *freqgsl;
	gsl_vector *csdgsl;		//Amount of current per unit (density) of frequency (spectral), as a function of the frequency
	gsl_vector *sigmacsdgsl;

// FUNCTIONS

	int findInterval(int tail_duration, gsl_vector *invector, gsl_vector *startpulse, int npin, int pulse_length, int nPF, int interval, int *ni, gsl_vector **startinterval);
	int findIntervalN(gsl_vector *invector, int interval, int *ni, gsl_vector **startinterval);

	int createTPSreprFile ();
	int writeTPSreprExten ();

	int findPulsesNoise
	(
		gsl_vector *vectorin,
		gsl_vector *vectorinDER,
		gsl_vector **tstart,
		gsl_vector **quality,
		gsl_vector **energy,

		int *nPulses,
		double *threshold,

		double scalefactor,
		int sizepulsebins,
		double samplingRate,

		int samplesup,
		double nsgms,

		double lb,
		double lrs,

		double stopcriteriamkc,
		double kappamkc);
	
	int findTstartNoise (int maxPulsesPerRecord, gsl_vector *der, double adaptativethreshold, int nSamplesUp,
		int *numberPulses, gsl_vector **tstartgsl, gsl_vector **flagTruncated, gsl_vector **maxDERgsl);
	
	int weightMatrixNoise (gsl_matrix *intervalMatrix, gsl_matrix **weight);
        
        int medianKappaClipping_noiseSigma (gsl_vector *invector, double kappa, double stopCriteria, double nSigmas, double *mean, double *sigma);
        
        int getpar_noiseSpec(struct Parameters* const par);

        void MyAssert(int expr, char* msg);
        
        using namespace std;

#endif /*GENNOISESPEC_H_*/

