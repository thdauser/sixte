/*
   This file is part of SIXTE.

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


   Copyright 2015 Philippe Peille, IRAP
*/

#ifndef TESRECONSTRUCTION_H
#define TESRECONSTRUCTION_H 1

#include "optimalfilters.h"
#include "testriggerfile.h"
#include "teseventlist.h"
#include "tesproftemplates.h"
#include "testrigger.h"
#include "gti.h"
#include "integraSIRENA.h"

#define TOOLSUB tesreconstruction_main
#include "headas_main.c"

struct Parameters {
	//File containing the optimal filter
	char OptimalFilterFile[MAXFILENAME];

	//File to reconstruct
	char RecordFile[MAXFILENAME];

	//Ouput event list
	char TesEventFile[MAXFILENAME];

	//File containing the pulse template
	char PulseTemplateFile[MAXFILENAME];

	//Pulse Length
	int PulseLength;

	//Threshold level
	double Threshold;

	//Calibration factor
	double Calfac;

	//Default size of the event list
	int EventListSize;

	//Minimal distance before using OFs after a mireconstruction
	int NormalExclusion;

	//Minimal distance before reconstructing any event after a mireconstruction
	int DerivateExclusion;

	//Saturation level of the ADC curves
	double SaturationValue;

	//Boolean to choose whether to erase an already existing event list
	char clobber;

	//Boolean to choose to save the run parameters in the output file
	char history;

	//Reconstruction Method (PP or SIRENA)
	char Rcmethod[7];
	
	//
	// SIRENA parameters
	//
	//File containing the library
	char LibraryFile[MAXFILENAME];
	
	//Scale Factor for initial filtering
	double scaleFactor;
	
	//Fall time of the pulses
	double tauFall;
	
	//Number of samples for threshold trespassing
	double samplesUp;
	
	//Number of standard deviations in the kappa-clipping process for threshold estimation
	double nSgms;
	
	//Run for Library creation?: Y(1), N(0)
	int crtLib;
	
	//Last energy to be included in the library file?: Y(1), N(0)
	int lastELibrary;

	//Calibration run (0) or energy reconstruction run (1)?
	int mode;

	/** Monochromatic energy for library creation **/
	double monoenergy;
	
	/** Running sum length for the RS raw energy estimation, in seconds (only in crtLib=0) **/
	double LrsT;
	
	/** Baseline averaging length for the RS raw energy estimation, in seconds (only in crtLib=0) **/
	double LbT;
	
	/** Baseline (in ADC units) **/
	double baseline;

	//Noise filename
	char NoiseFile[MAXFILENAME];
	
	//Pixel Type: SPA, LPA1, LPA2 or LPA3
	char PixelType[5];

	//Filtering Domain: Time(T) or Frequency(F)
	char FilterDomain[2];

	//Filtering Method: F0 (deleting the zero frequency bin) or F0 (deleting the baseline) **/
	char FilterMethod[3];
	
	//Energy Method: NOLAGS, LAGS or WEIGHT **/
	char EnergyMethod[7];

	//Linear (1) or Quadratic calibration of the energies (2)
	int calibLQ;
	
	//Linear calibration factor
	double b_cF;
	
	//Quadratic calibration factor
	double c_cF;
	
	//Write intermediate files (Yes:1, No:0)
	int intermediate;
	
	// File with the output detections 
	char detectFile[256];
	
	// File with the output filter (only in calibration)
	char filterFile[256];
	
	// Second calibration file to calculate calibration factors b,c (if mode=0 & crtLib=0)
	char RecordFileCalib2[256];
	
	/** Monochromatic energy of the second calibration file to calculate calibration factors b,c (if mode=0 & crtLib=0)**/
	double monoenergy2;
	
	// END SIRENA PARAMETERS
};

int getpar(struct Parameters* const par);


#endif /* TESRECONSTRUCTION_H */
