#ifndef INITSIRENA_H
#define INITSIRENA_H 1

#include <stdio.h>
#include <stdlib.h>

#include <fitsio.h>

#include "sixt.h"
#include "advdet.h"
#include "testriggerfile.h"
#include "testrigger.h"

#include "integraSIRENA.h"

struct Parameters {
	// File containing the optimal filter
	char OptimalFilterFile[MAXFILENAME];

	// File to reconstruct
	char RecordFile[MAXFILENAME];

	// Ouput event list
	char TesEventFile[MAXFILENAME];

	// File containing the pulse template
	char PulseTemplateFile[MAXFILENAME];

	// Filter length not padded with 0s (only necessary when reconstructing with 0-padding)
	int OFLengthNotPadded;

	// Threshold level
	double Threshold;

	// Calibration factor
	double Calfac;

	// Default size of the event list
	int EventListSize;

	// Minimal distance before using OFs after a mireconstruction
	int NormalExclusion;

	// Minimal distance before reconstructing any event after a mireconstruction
	int DerivateExclusion;

	// Saturation level of the ADC curves
	double SaturationValue;

	// Boolean to choose whether to erase an already existing event list
	char clobber;

	// Boolean to choose to save the run parameters in the output file
	char history;

	// File containing the library
	char LibraryFile[MAXFILENAME];

	// Scale Factor for initial filtering
	double scaleFactor;

	// Number of samples for threshold trespassing
	int samplesUp;

    // Number of samples below the threshold to look for other pulse
	int samplesDown;

	// Number of standard deviations in the kappa-clipping process for threshold estimation
	double nSgms;

    // Detect secondary pulses or not
    int detectSP;

	// Calibration run (0) or energy reconstruction run (1)?
	int opmode;

    // DetectionMode: Adjusted Derivative(AD) or Single Threshold Crossing(STC)
	char detectionMode[4];

	// Monochromatic energy for library creation
	double monoenergy;

	// Boolean to choose whether to add the PRECALWN HDU in the library file
	char hduPRECALWN;

	// Boolean to choose whether to add the PRCLOFWM HDU in the library file
	char hduPRCLOFWM;

	// Length of the longest fixed filter for library creation
	int largeFilter;

	// Running sum length for the RS raw energy estimation, in seconds (only in CALIBRATION)
	double LrsT;

	// Baseline averaging length for the RS raw energy estimation, in seconds (only in CALIBRATION)
	double LbT;

	// Noise filename
	char NoiseFile[MAXFILENAME];

	// Filtering Domain: Time(T) or Frequency(F)
	char FilterDomain[2];

	// Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline) of F0B0 (deleting always the baseline)
	char FilterMethod[5];

	// Energy Method: OPTFILT, WEIGHT, WEIGHTN, I2R, I2RFITTED or PCA
	char EnergyMethod[10];

    // Energy of the filters of the library to be used to calculate energy (only for OPTFILT, I2R and I2RFITTED)
    double filtEev;

    // Constant to apply the I2RFITTED conversion
    double Ifit;

	// Noise to use with Optimal Filtering: NSD (Noise Spectral Density) or WEIGHTM (weight matrix)
	char OFNoise[8];

	// LagsOrNot: LAGS == 1 or NOLAGS == 0
	int LagsOrNot;

    // Number of lags (odd number)
	int nLags;

    // Using 3 lags to analytically calculate a parabola (3) or using 5 lags to fit (5)
	int Fitting35;

	// OFIter: Iterate == 1 or NOTIterate == 0
	int OFIter;

	// Boolean to choose whether to use a library with optimal filters or calculate the optimal filter to each pulse
	char OFLib;

	// Optimal Filter by using the Matched Filter (MF) or the DAB as matched filter (MF, DAB)
    char OFInterp[4];

	// Optimal Filter length Strategy: FREE, BASE2, BYGRADE or FIXED
	char OFStrategy[8];

	// Optimal Filter length (taken into account if OFStrategy=FIXED)
	int OFLength;

    // 0-padding filter if 0 (from pulseLength to OFLength filter filled in with 0's) or filter with a filter+preBuffer if different from 0
	char preBuffer;

	// Write intermediate files (Yes:1, No:0)
	int intermediate;

	// File with the output detections
	char detectFile[256];

	// File with the output filter (only in calibration)
	char filterFile[256];

    // Additional error (in samples) added to the detected time"  (Logically, it changes the reconstructed energies)
	int errorT;

    // Sum0Filt: 0-padding: Subtract the sum of the filter (1) or not (0)
	int Sum0Filt;

	// Tstart of the pulses (to be used instead of calculating them if tstartPulse1 =! 0)
	//int tstartPulse1;
    char tstartPulse1[MAXFILENAME]; // Integer number: Sample where the first pulse starts
                                    // or
                                    // nameFile: File where is the tstart (in seconds) of every pulse
	int tstartPulse2;
	int tstartPulse3;

	// Energies for PCA
	double energyPCA1;
	double energyPCA2;

	// XML file with instrument definition
	char XMLFile[MAXFILENAME];
};

void MyAssert(int expr, char* msg);

int checkXmls(struct Parameters* const par);

char* subString (const char* input, int offset, int len, char* dest);

int getSamplingrate_trigreclength (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength, int* numfits);
int getSamplingrate_trigreclength_Filei (char* inputFile, struct Parameters par, double* samplingrate, int* trigreclength);

int fillReconstructInitSIRENAGrading (struct Parameters par, AdvDet *det, ReconstructInitSIRENA** reconstruct_init_sirena);

int callSIRENA_Filei(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFile * outfile);
int callSIRENA(char* inputFile, SixtStdKeywords* keywords, ReconstructInitSIRENA* reconstruct_init_sirena,struct Parameters par, double sampling_rate, int *trig_reclength, PulsesCollection* pulsesAll, TesEventFile * outfile);

#endif /* INITSIRENA_H */
